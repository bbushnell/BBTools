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

import jgi.BBMergeOverlapper;
import parse.Parse;
import stream.Read;

/** One-shape timing of the real expected-error helper, never a merger/model benchmark.
 * A baseline pilot saves calls; later baseline/fixed measurements must read that file.
 * All validation, allocation, input checks and evidence writes are outside timing.
 * @author Sucrose
 */
public final class BBMergeExpectedMismatchesTiming{

	/** Uses mode=pilot|measure, phase=baseline|fixed, shape, new outdir and measure-only counts. */
	public static void main(final String[] args) throws Exception{
		String mode=null, phase=null, name=null, directory=null, counts=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("shape")){name=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("counts")){counts=pair[1];}
			else{throw new IllegalArgumentException("Unknown helper timing option: "+arg);}
		}
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		final boolean pilot=mode.equals("pilot"), baseline=phase.equals("baseline");
		require(!pilot || baseline, "Only baseline pilots may select call counts.");
		require(pilot ? counts==null : counts!=null && !counts.isEmpty(),
			"Pilot takes no counts; measure requires counts=<saved baseline pilot counts.tsv>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergeexpectedmismatchestiming.sh with assertions enabled.");
		final BBMergeExpectedMismatchesTiming test=new BBMergeExpectedMismatchesTiming(
			Paths.get(directory).toAbsolutePath().normalize(), shape(name), baseline);
		final Path saved=counts==null ? null : Paths.get(counts).toAbsolutePath().normalize();
		final byte[] savedBytes=saved==null ? null : Files.readAllBytes(saved);
		final Path savedRing=saved==null ? null : saved.getParent().resolve("ring.tsv");
		final byte[] savedRingBytes=savedRing==null ? null : Files.readAllBytes(savedRing);
		final long calls=pilot ? INITIAL_CALLS : test.savedCalls(savedBytes);
		test.validateRing();
		if(savedRing!=null){require(Arrays.equals(savedRingBytes, Files.readAllBytes(test.out.resolve("ring.tsv"))), "Input ring differs from the saved baseline pilot.");}
		test.sample("warmup", 0, INITIAL_CALLS);
		Result selected=null;
		long candidate=calls;
		for(int attempt=1; attempt<=(pilot ? MAX_ATTEMPTS : 1); attempt++){
			selected=test.sample(mode, attempt, candidate);
			if(!pilot || selected.cpuNanos>=MIN_CPU_NANOS){break;}
			if(attempt<MAX_ATTEMPTS){candidate=nextCalls(candidate, selected.cpuNanos);}
		}
		require(selected!=null, "A timing invocation must retain its measured result.");
		if(saved!=null){
			require(Arrays.equals(savedBytes, Files.readAllBytes(saved)) && Arrays.equals(savedRingBytes, Files.readAllBytes(savedRing)),
				"Saved baseline counts or input ring changed during measurement.");
		}
		final boolean longEnough=selected.cpuNanos>=MIN_CPU_NANOS;
		if(pilot && longEnough){test.writeCounts(selected);}
		final String summary="BBMERGE_EXPECTED_MISMATCHES_TIMING_"+(longEnough ? "PASS" : "SHORT")+
			" mode="+mode+" phase="+phase+" shape="+test.spec.name+" ring="+RING_SIZE+
			" validationCalls="+RING_SIZE+" calls="+selected.calls+" positionsPerCall="+test.positions+
			" evaluatedPositions="+Math.multiplyExact(selected.calls, (long)test.positions)+
			" positionCountEvidence=literal-geometry"+
			" sourcePositionsPerCall="+test.spec.sourceCount+" correctedPositionsPerCall="+test.spec.correctedCount+
			" correctedVsSourcePositionRatio="+(test.spec.correctedCount/(double)test.spec.sourceCount)+
			" wallNanos="+selected.wallNanos+" cpuNanos="+selected.cpuNanos+
			" checksum="+selected.checksum+" inputsUnchanged=true\n";
		Files.write(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		System.out.print(summary);
		if(!longEnough){throw new AssertionError("Helper compute below ten CPU seconds; increase baseline counts in a new matched pilot, not fixed alone.");}
	}

	/** Builds one deterministic ring and evidence table without calling the timed helper. */
	private BBMergeExpectedMismatchesTiming(final Path out_, final Shape spec_, final boolean baseline_) throws IOException{
		out=out_; spec=spec_; baseline=baseline_; positions=baseline ? spec.sourceCount : spec.correctedCount;
		clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Current-thread CPU accounting is required for the ten-second compute minimum.");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		Files.createDirectory(out);
		Files.write(out.resolve("measurements.tsv"), ("stage\tattempt\tphase\tshape\tcalls\tpositionsPerCall"+
			"\tevaluatedPositions\tpositionCountEvidence\tchecksum\texpectedChecksum\twallNanos\tcpuNanos"+
			"\tatLeastTenCpuSeconds\tinputsUnchanged\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		Files.write(out.resolve("runtime.txt"), ("java.version="+System.getProperty("java.version")+
			"\njava.vendor="+System.getProperty("java.vendor")+"\njava.vm.name="+System.getProperty("java.vm.name")+
			"\nrecipe="+RECIPE+"\npositionCountEvidence=literal-geometry; helper is not instrumented\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		for(int r=0; r<RING_SIZE; r++){
			final byte[] ab=new byte[spec.alen], bb=new byte[spec.blen];
			final byte[] aq=new byte[spec.alen], bq=new byte[spec.blen];
			Arrays.fill(ab, (byte)'A'); Arrays.fill(bb, (byte)'A');
			for(int p=0; p<aq.length; p++){aq[p]=QUALITIES[(p+r+(p*r)%7)&3];}
			for(int p=0; p<bq.length; p++){bq[p]=QUALITIES[(3*p+5*r+(p*r)%11)&3];}
			ring[r]=new Pair(new Read(ab, aq, "timing-A", r, false), new Read(bb, bq, "timing-B-rc", r, false));
		}
	}

	/** Checks each ring result against independently indexed decimal arithmetic, before timing. */
	private void validateRing() throws IOException{
		final StringBuilder inputs=new StringBuilder("ring\taBases\tbBases\taQualities\tbQualities\n");
		final StringBuilder validation=new StringBuilder("ring\tsourceAStart\tsourceBStart\tsourcePositions"+
			"\tcorrectedAStart\tcorrectedBStart\tcorrectedPositions\tsourceDecimalOracle\tcorrectedDecimalOracle"+
			"\tobservedFloat\tobservedRawBitsHex\toracleMatches\tinputsUnchanged\n");
		int failures=0, different=0;
		for(int r=0; r<RING_SIZE; r++){
			final Pair pair=ring[r];
			final double source=oracle(pair, spec.sourceA, spec.sourceB, spec.sourceCount);
			final double corrected=oracle(pair, spec.correctedA, spec.correctedB, spec.correctedCount);
			if(Math.abs(source-corrected)>ORACLE_TOLERANCE*2){different++;}
			final float observed=BBMergeOverlapper.expectedMismatches(pair.a, pair.b, spec.insert);
			final boolean match=Float.isFinite(observed) && Math.abs((double)observed-(baseline ? source : corrected))<=ORACLE_TOLERANCE;
			final boolean unchanged=pair.unchanged();
			if(!match || !unchanged){failures++;}
			ringBits[r]=Float.floatToRawIntBits(observed);
			inputs.append(r).append('\t').append(new String(pair.ab, StandardCharsets.US_ASCII)).append('\t')
				.append(new String(pair.bb, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(pair.aq))
				.append('\t').append(Arrays.toString(pair.bq)).append('\n');
			validation.append(r).append('\t').append(spec.sourceA).append('\t').append(spec.sourceB).append('\t').append(spec.sourceCount)
				.append('\t').append(spec.correctedA).append('\t').append(spec.correctedB).append('\t').append(spec.correctedCount)
				.append('\t').append(source).append('\t').append(corrected).append('\t').append(Float.toString(observed))
				.append('\t').append(Integer.toHexString(ringBits[r])).append('\t').append(match).append('\t').append(unchanged).append('\n');
		}
		Files.write(out.resolve("ring.tsv"), inputs.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		Files.write(out.resolve("validation.tsv"), validation.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		require(failures==0, "Ring oracle/input failures="+failures+"; inspect validation.tsv before timing.");
		require(spec.changes ? different>0 : different==0, "Literal ring must distinguish affected shapes and preserve controls: "+different);
	}

	/** Decimal oracle uses explicit shape ranges, never the helper's start/stop expressions. */
	private static double oracle(final Pair pair, final int aStart, final int bStart, final int count){
		require(aStart>=0 && bStart>=0 && count>0 && aStart+count<=pair.aq.length && bStart+count<=pair.bq.length,
			"Literal oracle indices must remain inside both quality arrays.");
		double sum=0;
		for(int offset=0; offset<count; offset++){
			sum+=1-correctness(pair.aq[aStart+offset])*correctness(pair.bq[bStart+offset]);
		}
		return sum;
	}

	/** These are decimal table literals at Q10,Q20,Q30,Q40, not a phred formula. */
	private static double correctness(final byte quality){
		switch(quality){
			case 10:return .9;
			case 20:return .99;
			case 30:return .999;
			case 40:return .9999;
			default:throw new IllegalArgumentException("Unsupported ring quality: "+quality);
		}
	}

	/** Times helper calls plus ring lookup/raw-bit checksum overhead; checks happen after both clocks stop. */
	private Result sample(final String stage, final int attempt, final long calls) throws IOException{
		require(calls>0 && calls<=MAX_CALLS && calls%RING_SIZE==0, "Calls must be bounded complete ring cycles: "+calls);
		long checksum=0;
		final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
		for(long i=0; i<calls; i++){
			final Pair pair=ring[(int)i&(RING_SIZE-1)];
			checksum+=Float.floatToRawIntBits(BBMergeOverlapper.expectedMismatches(pair.a, pair.b, spec.insert));
		}
		final long wallNanos=System.nanoTime()-wallStart, cpuNanos=clock.getCurrentThreadCpuTime()-cpuStart;
		final Result result=new Result(calls, checksum, wallNanos, cpuNanos);
		long cycleChecksum=0;
		for(final int bits : ringBits){cycleChecksum+=bits;}
		final long expected=Math.multiplyExact(cycleChecksum, calls/RING_SIZE);
		boolean unchanged=true;
		for(final Pair pair : ring){unchanged&=pair.unchanged();}
		Files.write(out.resolve("measurements.tsv"), (stage+"\t"+attempt+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+spec.name+"\t"+calls+"\t"+positions+"\t"+Math.multiplyExact(calls, (long)positions)+
			"\tliteral-geometry\t"+checksum+"\t"+expected+"\t"+wallNanos+"\t"+cpuNanos+
			"\t"+(cpuNanos>=MIN_CPU_NANOS)+"\t"+unchanged+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
		require(checksum==expected && unchanged, "Timed checksum/input mismatch; inspect measurements.tsv.");
		require(wallNanos>0 && cpuNanos>0, "Both measured clocks must advance; wall="+wallNanos+", cpu="+cpuNanos);
		return result;
	}

	/** Reads only a successful baseline pilot's matching recipe/shape/count row. */
	private long savedCalls(final byte[] bytes){
		final String[] lines=new String(bytes, StandardCharsets.US_ASCII).split("\n");
		require(lines.length==2 && lines[0].equals(COUNTS_HEADER), "Expected the two-line saved baseline counts.tsv.");
		final String[] fields=lines[1].split("\t", -1);
		require(fields.length==9 && fields[0].equals(RECIPE) && fields[1].equals(spec.name) &&
			fields[2].equals("baseline") && fields[3].equals("pilot") && Integer.parseInt(fields[4])==RING_SIZE &&
			Integer.parseInt(fields[6])==spec.sourceCount && Integer.parseInt(fields[7])==spec.correctedCount,
			"Saved baseline recipe, phase, ring or geometry differs from this shape.");
		final long calls=Long.parseLong(fields[5]);
		require(calls>0 && calls<=MAX_CALLS && calls%RING_SIZE==0 && Long.parseLong(fields[8])>=MIN_CPU_NANOS,
			"Saved pilot must contain bounded full-ring calls and at least ten CPU seconds.");
		return calls;
	}

	/** Publishes selected counts only after a sufficiently long, checked baseline sample. */
	private void writeCounts(final Result result) throws IOException{
		require(baseline && result.cpuNanos>=MIN_CPU_NANOS, "Only a validated baseline pilot can select matched call counts.");
		Files.write(out.resolve("counts.tsv"), (COUNTS_HEADER+"\n"+RECIPE+"\t"+spec.name+"\tbaseline\tpilot\t"+
			RING_SIZE+"\t"+result.calls+"\t"+spec.sourceCount+"\t"+spec.correctedCount+"\t"+result.cpuNanos+"\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Grows pilot work toward eleven CPU seconds, at most sixteenfold, in full ring cycles. */
	private static long nextCalls(final long calls, final long cpuNanos){
		require(cpuNanos>0, "Cannot size a pilot without measured CPU time.");
		final double scaled=Math.ceil(calls*(TARGET_CPU_NANOS/(double)cpuNanos));
		final long next=(long)Math.min(calls*16.0, scaled);
		final long aligned=((Math.max(calls+RING_SIZE, next)+RING_SIZE-1)/RING_SIZE)*RING_SIZE;
		require(aligned<=MAX_CALLS, "Pilot would exceed the call cap; review sizing instead of launching larger work.");
		return aligned;
	}

	/** Selects one of four explicit source/corrected coordinate shapes. */
	private static Shape shape(final String name){
		for(final Shape candidate : SHAPES){if(candidate.name.equals(name)){return candidate;}}
		throw new IllegalArgumentException("Set shape=equal-long|unequal-long|short-a|short-b.");
	}

	/** Rejects invalid setup or results loudly; never reports a partial pilot as passing. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Snapshots only the base/quality arrays; the hot loop creates no Read or array objects. */
	private static final class Pair{
		Pair(final Read a_, final Read b_){
			a=a_; b=b_; ab=a.bases.clone(); bb=b.bases.clone(); aq=a.quality.clone(); bq=b.quality.clone();
			abRef=a.bases; bbRef=b.bases; aqRef=a.quality; bqRef=b.quality;
		}
		boolean unchanged(){
			return a.bases==abRef && b.bases==bbRef && a.quality==aqRef && b.quality==bqRef &&
				Arrays.equals(ab, a.bases) && Arrays.equals(bb, b.bases) &&
				Arrays.equals(aq, a.quality) && Arrays.equals(bq, b.quality);
		}
		final Read a, b;
		final byte[] ab, bb, aq, bq, abRef, bbRef, aqRef, bqRef;
	}

	/** Positive finite helper values keep the bounded sum of raw float bits inside a long. */
	private static final class Result{
		Result(final long calls_, final long checksum_, final long wallNanos_, final long cpuNanos_){
			calls=calls_; checksum=checksum_; wallNanos=wallNanos_; cpuNanos=cpuNanos_;
		}
		final long calls, checksum, wallNanos, cpuNanos;
	}

	/** Starts/counts are literal paired-index ranges, not copies of the helper's coordinate formula. */
	private static final class Shape{
		Shape(final String name_, final int alen_, final int blen_, final int insert_,
				final int sa, final int sb, final int sn, final int ca, final int cb, final int cn, final boolean changes_){
			name=name_; alen=alen_; blen=blen_; insert=insert_;
			sourceA=sa; sourceB=sb; sourceCount=sn; correctedA=ca; correctedB=cb; correctedCount=cn; changes=changes_;
		}
		final String name;
		final int alen, blen, insert, sourceA, sourceB, sourceCount, correctedA, correctedB, correctedCount;
		final boolean changes;
	}

	private final Path out;
	private final Shape spec;
	private final boolean baseline;
	private final int positions;
	private final ThreadMXBean clock;
	private final Pair[] ring=new Pair[RING_SIZE];
	private final int[] ringBits=new int[RING_SIZE];
	private static final int RING_SIZE=64, MAX_ATTEMPTS=6;
	private static final long INITIAL_CALLS=65536, MAX_CALLS=536870912;
	private static final long MIN_CPU_NANOS=10000000000L, TARGET_CPU_NANOS=11000000000L;
	/** At most150 terms of at most.19: float summation bound plus table/product rounding is below.0005. */
	private static final double ORACLE_TOLERANCE=0.0005;
	private static final byte[] QUALITIES={10, 20, 30, 40};
	private static final String RECIPE="expected-mismatches-ring-v1";
	private static final String COUNTS_HEADER="recipe\tshape\tphase\tmode\tringSize\tcalls\tsourcePositionsPerCall\tcorrectedPositionsPerCall\tpilotCpuNanos";
	private static final Shape[] SHAPES={
		new Shape("equal-long", 150, 150, 200, 50, 0, 100, 50, 0, 100, false),
		new Shape("unequal-long", 150, 250, 300, 50, 0, 100, 50, 0, 100, false),
		new Shape("short-a", 150, 250, 150, 0, 0, 150, 0, 100, 150, true),
		new Shape("short-b", 250, 150, 150, 0, 100, 50, 0, 0, 150, true)
	};
}
