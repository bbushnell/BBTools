package test.jgi;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import jgi.BBMerge;
import parse.Parse;
import stream.Read;

/** Measures the real retained-quality helper with a fixed ring of owned read pairs.
 * Reference/metadata restoration, ring lookup and an output checksum are timed too.
 * The helper allocates consensus and sliced outputs; this is not allocation-free
 * timing or a whole-pipeline benchmark. Inputs and phase oracles are generated before
 * timing, checked in full outside the clocks, and pinned to the baseline pilot.
 * @author Shinobu
 */
public final class BBMergeRetainedQualityTiming{

	/** mode=pilot|measure, phase=baseline|fixed, shape=full|short, outdir, measure-only counts. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "shape", "outdir", "counts").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		final String mode=options.get("mode"), phase=options.get("phase"), shape=options.get("shape");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require("full".equals(shape) || "short".equals(shape), "Set shape=full|short.");
		final boolean pilot=mode.equals("pilot"), baseline=phase.equals("baseline");
		require(!pilot || baseline, "Only the baseline selects matched call counts.");
		require(pilot ? !options.containsKey("counts") : options.containsKey("counts"), "Only measure requires counts=BASELINE_PILOT/counts.tsv.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use the developer launcher with assertions enabled.");
		final Path out=path(options.get("outdir")), saved=pilot ? null : path(options.get("counts"));
		final byte[][] reference=saved==null ? null : readReference(saved);
		final boolean oldQuality=BBMerge.changeQuality, oldReadQuality=Read.CHANGE_QUALITY;
		String summary;
		try{
			BBMerge.changeQuality=Read.CHANGE_QUALITY=false;
			final BBMergeRetainedQualityTiming test=new BBMergeRetainedQualityTiming(out, shape, baseline);
			final long calls=pilot ? INITIAL_CALLS : test.savedCalls(reference[0]);
			final byte[] inputs=test.inputs().getBytes(StandardCharsets.US_ASCII);
			final byte[] oracles=test.oracles().getBytes(StandardCharsets.US_ASCII);
			if(reference!=null){
				require(Arrays.equals(inputs, reference[1]) && Arrays.equals(oracles, reference[2]), "Input ring or coordinate oracles differ from the baseline pilot.");
			}
			Files.write(out.resolve("ring.tsv"), inputs, StandardOpenOption.CREATE_NEW);
			Files.write(out.resolve("oracles.tsv"), oracles, StandardOpenOption.CREATE_NEW);
			final String validation=test.validateRing();
			if(reference!=null && (baseline || shape.equals("full"))){
				require(Arrays.equals(validation.getBytes(StandardCharsets.US_ASCII), reference[3]), "Baseline/full-length output bytes differ from the saved pilot.");
			}
			write(out.resolve("validation.tsv"), validation);
			test.sample("warmup", 0, INITIAL_CALLS);
			Result selected=null;
			long candidate=calls;
			for(int attempt=1; attempt<=(pilot ? MAX_ATTEMPTS : 1); attempt++){
				selected=test.sample(mode, attempt, candidate);
				if(!pilot || selected.cpu>=MIN_PILOT_CPU){break;}
				if(attempt<MAX_ATTEMPTS){candidate=nextCalls(candidate, selected.cpu);}
			}
			require(selected!=null, "Every invocation must retain its measured result.");
			final String after=test.validateRing();
			require(validation.equals(after), "Full ring outputs changed after timing.");
			write(out.resolve("validation-after.tsv"), after);
			if(reference!=null){
				final byte[][] afterReference=readReference(saved);
				for(int i=0; i<reference.length; i++){
					require(Arrays.equals(reference[i], afterReference[i]), "Saved baseline input "+i+" changed during measurement.");
				}
			}
			final boolean enough=selected.cpu>=(pilot ? MIN_PILOT_CPU : MIN_CPU);
			if(pilot && enough){
				write(out.resolve("counts.tsv"), COUNTS_HEADER+"\n"+RECIPE+"\t"+shape+"\tbaseline\tpilot\t64\t"+selected.calls+
					"\t150\t"+test.bLength+"\t150\t"+selected.cpu+"\n");
			}
			summary="BBMERGE_RETAINED_QUALITY_TIMING_"+(enough ? "PASS" : "SHORT")+" mode="+mode+" phase="+phase+" shape="+shape+
				" ring=64 validationCalls=128 calls="+selected.calls+" r1Length=150 r2Length="+test.bLength+" insert=150"+
				" wallNanos="+selected.wall+" cpuNanos="+selected.cpu+" checksum="+selected.checksum+" originalArraysUnchanged=true\n";
			if(!enough){
				write(out.resolve("summary.txt"), summary);
				throw new AssertionError("Compute region too short; recalibrate matched BASELINE counts, never increase only fixed.\n"+summary);
			}
		}finally{
			BBMerge.changeQuality=oldQuality;
			Read.CHANGE_QUALITY=oldReadQuality;
		}
		require(BBMerge.changeQuality==oldQuality && Read.CHANGE_QUALITY==oldReadQuality, "Owned quality settings were not restored.");
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	private BBMergeRetainedQualityTiming(final Path out_, final String shape_, final boolean baseline_) throws IOException{
		out=out_; shape=shape_; baseline=baseline_;
		bLength=shape.equals("full") ? LENGTH : LENGTH+PREFIX.length();
		clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Current-thread CPU time is required; wall time is not a substitute.");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		require(Read.MAX_CALLED_QUALITY>=39, "The frozen input qualities must not be capped.");
		Files.createDirectory(out);
		for(int i=0; i<ring.length; i++){ring[i]=new Pair(i, bLength-LENGTH);}
		write(out.resolve("measurements.tsv"), "stage\tattempt\tphase\tshape\tcalls\tchecksum\texpectedChecksum\twallNanos\tcpuNanos\toriginalArraysUnchanged\n");
	}

	/** Validates all output bytes and source-array preservation outside the timed loop. */
	private String validateRing(){
		final StringBuilder table=new StringBuilder("ring\treturn\tA\tB\tAquality\tBquality\tAerrors\tBerrors\tAchrom\tBchrom\tchecksum\n");
		cycleChecksum=0;
		for(int i=0; i<ring.length; i++){
			final Pair p=ring[i]; p.restore();
			final int result=BBMerge.errorCorrectWithInsert(p.a, p.b, LENGTH);
			final byte[] wanted=baseline ? p.sourceQuality : p.fixedQuality;
			require(result==0 && p.a.errors==0 && p.b.errors==0, "Matching fragment calls must reset the sentinel error counters: slot "+i);
			require(Arrays.equals(p.a.bases, p.fragment) && Arrays.equals(p.b.bases, p.fragment), "Fragment bases or lengths differ: slot "+i);
			require(Arrays.equals(p.a.quality, p.aqSnapshot) && Arrays.equals(p.b.quality, wanted), "Retained qualities differ from the phase coordinate oracle: slot "+i);
			require(p.originalsUnchanged() && p.metadataMatches(), "Original arrays or retained metadata changed unexpectedly: slot "+i);
			final long expected=p.expectedChecksum(baseline), actual=p.checksum(result);
			require(actual==expected, "Checksum does not consume the expected output: slot "+i);
			cycleChecksum=Math.addExact(cycleChecksum, expected);
			table.append(i).append('\t').append(result).append('\t').append(new String(p.a.bases, StandardCharsets.US_ASCII)).append('\t')
				.append(new String(p.b.bases, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(p.a.quality)).append('\t')
				.append(Arrays.toString(p.b.quality)).append('\t').append(p.a.errors).append('\t').append(p.b.errors).append('\t')
				.append(p.a.chrom).append('\t').append(p.b.chrom).append('\t').append(actual).append('\n');
		}
		return table.toString();
	}

	/** Times restoration, ring lookup, the actual allocating helper and a small output sink. */
	private Result sample(final String stage, final int attempt, final long calls) throws IOException{
		require(calls>0 && calls<=MAX_CALLS && calls%RING==0, "The checksum model requires bounded whole ring cycles: "+calls);
		long checksum=0;
		final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
		for(long i=0; i<calls; i++){
			final Pair p=ring[(int)i&(RING-1)];
			p.restore();
			final int returned=BBMerge.errorCorrectWithInsert(p.a, p.b, LENGTH);
			checksum+=p.checksum(returned);
		}
		final long wall=System.nanoTime()-wallStart, cpu=clock.getCurrentThreadCpuTime()-cpuStart;
		final long expected=Math.multiplyExact(cycleChecksum, calls/RING);
		boolean unchanged=true;
		for(final Pair p : ring){unchanged&=p.originalsUnchanged();}
		Files.write(out.resolve("measurements.tsv"), (stage+"\t"+attempt+"\t"+(baseline ? "baseline" : "fixed")+"\t"+shape+"\t"+calls+
			"\t"+checksum+"\t"+expected+"\t"+wall+"\t"+cpu+"\t"+unchanged+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
		require(checksum==expected && unchanged && cpu>0 && wall>0, "Timed output checksum, original arrays or CPU clock failed.");
		return new Result(calls, checksum, wall, cpu);
	}

	private String inputs(){
		final StringBuilder s=new StringBuilder("ring\tA\tB\tAquality\tBquality\tAchrom\tBchrom\tAflags\tBflags\tAstart\tAstop\tBstart\tBstop\tAmapScore\tBmapScore\tAinsert\tBinsert\tAcopies\tBcopies\tAerrors\tBerrors\n");
		for(int i=0; i<ring.length; i++){
			final Pair p=ring[i];
			s.append(i).append('\t').append(new String(p.ab, StandardCharsets.US_ASCII)).append('\t').append(new String(p.bb, StandardCharsets.US_ASCII))
				.append('\t').append(Arrays.toString(p.aq)).append('\t').append(Arrays.toString(p.bq)).append('\t')
				.append(p.aChrom).append('\t').append(p.bChrom).append('\t').append(p.aFlags).append('\t').append(p.bFlags).append('\t')
				.append(p.aStart).append('\t').append(p.aStop).append('\t').append(p.bStart).append('\t').append(p.bStop).append('\t')
				.append(p.aMapScore).append('\t').append(p.bMapScore).append('\t').append(p.aInsert).append('\t').append(p.bInsert).append('\t')
				.append(p.aCopies).append('\t').append(p.bCopies).append("\t7\t9\n");
		}
		return s.toString();
	}

	private String oracles(){
		final StringBuilder s=new StringBuilder("ring\tfragment\tAquality\tsourceBquality\tfixedBquality\n");
		for(int i=0; i<ring.length; i++){
			final Pair p=ring[i];
			s.append(i).append('\t').append(new String(p.fragment, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(p.aqSnapshot))
				.append('\t').append(Arrays.toString(p.sourceQuality)).append('\t').append(Arrays.toString(p.fixedQuality)).append('\n');
		}
		return s.toString();
	}

	private long savedCalls(final byte[] bytes){
		final String[] rows=new String(bytes, StandardCharsets.US_ASCII).split("\n");
		require(rows.length==2 && rows[0].equals(COUNTS_HEADER), "Saved baseline counts have the wrong schema.");
		final String[] f=rows[1].split("\t", -1);
		require(f.length==10 && f[0].equals(RECIPE) && f[1].equals(shape) && f[2].equals("baseline") && f[3].equals("pilot") &&
			f[4].equals("64") && f[6].equals("150") && f[7].equals(Integer.toString(bLength)) && f[8].equals("150"), "Saved counts do not describe this frozen workload.");
		final long calls=Long.parseLong(f[5]), cpu=Long.parseLong(f[9]);
		require(calls>0 && calls<=MAX_CALLS && calls%RING==0 && cpu>=MIN_PILOT_CPU, "Saved pilot must meet the CPU minimum and whole-ring count bounds.");
		return calls;
	}

	private static byte[][] readReference(final Path counts) throws IOException{
		final Path parent=counts.getParent();
		return new byte[][]{Files.readAllBytes(counts), Files.readAllBytes(parent.resolve("ring.tsv")),
			Files.readAllBytes(parent.resolve("oracles.tsv")), Files.readAllBytes(parent.resolve("validation.tsv"))};
	}
	private static long nextCalls(final long calls, final long cpu){
		require(cpu>0, "Calibration requires measured thread CPU time.");
		final long scaled=(long)Math.min(calls*16.0, Math.ceil(calls*(TARGET_CPU/(double)cpu)));
		final long next=((Math.max(calls+RING, scaled)+RING-1)/RING)*RING;
		require(next<=MAX_CALLS, "Pilot call cap exceeded; review sizing instead of silently expanding work.");
		return next;
	}
	private static Path path(final String s){require(s!=null && !s.isEmpty(), "Missing required path."); return Paths.get(s).toAbsolutePath().normalize();}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static void write(final Path p, final String s) throws IOException{Files.write(p, s.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}

	private static final class Pair{
		Pair(final int index, final int prefix){
			require(index>=0 && index<RING && (prefix==0 || prefix==PREFIX.length()), "Only the declared ring slots and two input geometries are supported.");
			fragment=new byte[LENGTH];
			for(int i=0; i<LENGTH; i++){fragment[i]=(byte)FRAGMENT.charAt((i+index)%FRAGMENT.length());}
			ab=fragment.clone(); bb=new byte[LENGTH+prefix];
			for(int i=0; i<prefix; i++){bb[i]=(byte)PREFIX.charAt(i);}
			System.arraycopy(fragment, 0, bb, prefix, LENGTH);
			aq=new byte[LENGTH]; Arrays.fill(aq, (byte)30);
			bq=new byte[LENGTH+prefix];
			for(int i=0; i<bq.length; i++){bq[i]=(byte)(20+(index+i)%20);}
			sourceQuality=new byte[LENGTH]; fixedQuality=new byte[LENGTH];
			for(int i=0; i<LENGTH; i++){
				sourceQuality[i]=(byte)(20+(index+i)%20);
				fixedQuality[i]=(byte)(20+(index+prefix+i)%20);
			}
			abSnapshot=ab.clone(); bbSnapshot=bb.clone(); aqSnapshot=aq.clone(); bqSnapshot=bq.clone();
			id="retained-quality-timing-"+index;
			a=new Read(ab, aq, id, index, false); b=new Read(bb, bq, id, index, false);
			a.mate=b; b.mate=a; b.setPairnum(1);
			// Match the public overlap caller's validated input state outside the clocks.
			a.validate(true); b.validate(true);
			aChrom=a.chrom; bChrom=b.chrom; aFlags=a.flags; bFlags=b.flags;
			aStart=a.start; aStop=a.stop; bStart=b.start; bStop=b.stop; aMapScore=a.mapScore; bMapScore=b.mapScore;
			aInsert=a.insert(); bInsert=b.insert(); aCopies=a.copies; bCopies=b.copies;
			numericID=index; position=(index*7)%LENGTH;
		}
		/** Restore references and the source metadata that the helper may consume/mutate. */
		void restore(){
			a.bases=ab; b.bases=bb; a.quality=aq; b.quality=bq;
			a.errors=7; b.errors=9; a.chrom=aChrom; b.chrom=bChrom; a.flags=aFlags; b.flags=bFlags;
		}
		boolean originalsUnchanged(){return Arrays.equals(ab, abSnapshot) && Arrays.equals(bb, bbSnapshot) && Arrays.equals(aq, aqSnapshot) && Arrays.equals(bq, bqSnapshot);}
		boolean metadataMatches(){
			return a.mate==b && b.mate==a && a.id.equals(id) && b.id.equals(id) && a.numericID==numericID && b.numericID==numericID &&
				a.flags==aFlags && b.flags==bFlags && a.chrom==0 && b.chrom==bChrom &&
				a.start==aStart && a.stop==aStop && b.start==bStart && b.stop==bStop &&
				a.mapScore==aMapScore && b.mapScore==bMapScore && a.insert()==aInsert && b.insert()==bInsert &&
				a.copies==aCopies && b.copies==bCopies;
		}
		long checksum(final int returned){return returned+17L*a.bases.length+19L*b.bases.length+31L*a.bases[position]+131L*b.bases[position]+257L*a.quality[position]+65537L*b.quality[position]+a.errors+b.errors;}
		long expectedChecksum(final boolean baseline){return 36L*LENGTH+162L*fragment[position]+257L*30+65537L*(baseline ? sourceQuality[position] : fixedQuality[position]);}
		final Read a, b;
		final byte[] fragment, ab, bb, aq, bq, abSnapshot, bbSnapshot, aqSnapshot, bqSnapshot, sourceQuality, fixedQuality;
		final String id;
		final int numericID, position, aChrom, bChrom, aFlags, bFlags, aStart, aStop, bStart, bStop;
		final int aMapScore, bMapScore, aInsert, bInsert, aCopies, bCopies;
	}
	private static final class Result{
		Result(final long calls_, final long checksum_, final long wall_, final long cpu_){calls=calls_; checksum=checksum_; wall=wall_; cpu=cpu_;}
		final long calls, checksum, wall, cpu;
	}

	private final Path out;
	private final String shape;
	private final boolean baseline;
	private final int bLength;
	private final ThreadMXBean clock;
	private final Pair[] ring=new Pair[RING];
	private long cycleChecksum;
	private static final int RING=64, LENGTH=150, MAX_ATTEMPTS=6;
	private static final long INITIAL_CALLS=32768, MAX_CALLS=536870912;
	private static final long MIN_CPU=10000000000L, MIN_PILOT_CPU=15000000000L, TARGET_CPU=18000000000L;
	private static final String RECIPE="retained-quality-ring-v1";
	private static final String COUNTS_HEADER="recipe\tshape\tphase\tmode\tring\tcalls\tr1Length\tr2Length\tinsert\tpilotCpuNanos";
	private static final String FRAGMENT="ACGTTGCATGTCAGTACGATCGTACCTGAGCTAGTCGATC";
	private static final String PREFIX="TGTGCTCTTCCGATCT";
}
