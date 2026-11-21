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
import java.util.Locale;

import jgi.BBMergeOverlapper;
import parse.Parse;
import stream.Read;

/** Matched kernel timing for the real public probability helper.
 * A baseline pilot selects calls; both revisions reuse that immutable count and input
 * ring. Normal150/Q30 cases retain baseline score bits; subnormal1024/Q4 cases use
 * endpoint witnesses before repair and grouped ratio oracles after repair. Setup,
 * allocation, validation, snapshots and output are outside the measured region.
 * Measured work includes ring lookup and a raw-bit checksum, not model/merger work.
 * Uses the retained expected-mismatch fixture's CPU-time and saved-count conventions.
 * @author Shinobu
 */
public final class BBMergeProbabilityTiming{

	/** Requires mode=pilot|measure, phase=baseline|fixed, shape, new outdir and measure-only counts. */
	public static void main(final String[] args) throws Exception{
		String mode=null, phase=null, shape=null, directory=null, counts=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("mode")){require(mode==null, "Duplicate mode"); mode=pair[1];}
			else if(pair[0].equals("phase")){require(phase==null, "Duplicate phase"); phase=pair[1];}
			else if(pair[0].equals("shape")){require(shape==null, "Duplicate shape"); shape=pair[1];}
			else if(pair[0].equals("outdir")){require(directory==null, "Duplicate outdir"); directory=pair[1];}
			else if(pair[0].equals("counts")){require(counts==null, "Duplicate counts"); counts=pair[1];}
			else{throw new IllegalArgumentException("Unknown probability timing option: "+arg);}
		}
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require("normal".equals(shape) || "subnormal".equals(shape), "Set shape=normal|subnormal.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		final boolean pilot=mode.equals("pilot"), baseline=phase.equals("baseline");
		require(!pilot || baseline, "Only the original implementation selects matched call counts.");
		require(pilot ? counts==null : counts!=null && !counts.isEmpty(), "Only measure requires counts=BASELINE_PILOT/counts.tsv.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergeprobabilitytiming.sh with assertions enabled.");
		final BBMergeProbabilityTiming test=new BBMergeProbabilityTiming(
			Paths.get(directory).toAbsolutePath().normalize(), shape, baseline);
		final Path saved=counts==null ? null : Paths.get(counts).toAbsolutePath().normalize();
		final byte[] savedCounts=saved==null ? null : Files.readAllBytes(saved);
		final byte[] savedRing=saved==null ? null : Files.readAllBytes(saved.getParent().resolve("ring.tsv"));
		final byte[] savedValidation=saved==null ? null : Files.readAllBytes(saved.getParent().resolve("validation.tsv"));
		final long calls=pilot ? INITIAL_CALLS : test.savedCalls(savedCounts);
		test.validateRing(savedRing, savedValidation);
		test.sample("warmup", 0, INITIAL_CALLS);
		Result selected=null;
		long candidate=calls;
		for(int attempt=1; attempt<=(pilot ? MAX_ATTEMPTS : 1); attempt++){
			selected=test.sample(mode, attempt, candidate);
			if(!pilot || selected.cpuNanos>=MIN_PILOT_NANOS){break;}
			if(attempt<MAX_ATTEMPTS){candidate=nextCalls(candidate, selected.cpuNanos);}
		}
		require(selected!=null, "Every timing invocation must retain its measured result.");
		if(saved!=null){
			require(Arrays.equals(savedCounts, Files.readAllBytes(saved)) &&
				Arrays.equals(savedRing, Files.readAllBytes(saved.getParent().resolve("ring.tsv"))) &&
				Arrays.equals(savedValidation, Files.readAllBytes(saved.getParent().resolve("validation.tsv"))),
				"Saved baseline counts, inputs or score bits changed during measurement.");
		}
		final boolean enough=selected.cpuNanos>=(pilot ? MIN_PILOT_NANOS : MIN_CPU_NANOS);
		if(pilot && enough){test.writeCounts(selected);}
		final String summary="BBMERGE_PROBABILITY_TIMING_"+(enough ? "PASS" : "SHORT")+
			" mode="+mode+" phase="+phase+" shape="+shape+" ring="+RING_SIZE+
			" validationCalls="+RING_SIZE+" calls="+selected.calls+" coveredPositionsPerCall="+test.length+
			" wallNanos="+selected.wallNanos+" cpuNanos="+selected.cpuNanos+
			" checksum="+selected.checksum+" inputsUnchanged=true\n";
		write(test.out.resolve("summary.txt"), summary); System.out.print(summary);
		if(!enough){throw new AssertionError("Compute region too short; select new matched BASELINE counts, never increase fixed alone.");}
	}

	/** Creates a deterministic rotating ring; no helper call occurs during construction. */
	private BBMergeProbabilityTiming(final Path out_, final String shape_, final boolean baseline_) throws IOException{
		out=out_; shape=shape_; baseline=baseline_; normal=shape.equals("normal");
		length=normal ? 150 : 1024; quality=normal ? 30 : 4;
		clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Thread CPU time is required for the compute-duration minimum.");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		Files.createDirectory(out);
		write(out.resolve("measurements.tsv"), "stage\tattempt\tphase\tshape\tcalls\tcoveredPositionsPerCall\tchecksum\texpectedChecksum\twallNanos\tcpuNanos\tatLeastTenCpuSeconds\tinputsUnchanged\n");
		write(out.resolve("runtime.txt"), "java.version="+System.getProperty("java.version")+
			"\njava.vendor="+System.getProperty("java.vendor")+"\njava.vm.name="+System.getProperty("java.vm.name")+
			"\nrecipe="+RECIPE+"\npositions=covered input pairs, not instrumented visits; fixed fallback may rescan\n");
		for(int r=0; r<RING_SIZE; r++){
			final byte[] ab=new byte[length], bb=new byte[length], aq=new byte[length], bq=new byte[length];
			Arrays.fill(ab, (byte)'A'); Arrays.fill(bb, (byte)'A');
			Arrays.fill(aq, (byte)quality); Arrays.fill(bq, (byte)quality);
			final int mismatches=normal ? r%3 : 1+(r/2)%3;
			final boolean tail=(r&1)!=0;
			for(int j=0; j<mismatches; j++){
				final int p=normal ? (r*37+11+j*47)%length : tail ? length-mismatches+j : j;
				bb[p]='C';
			}
			ring[r]=new Pair(new Read(ab, aq, "timing-A", r, false), new Read(bb, bq, "timing-B-aligned", r, false), mismatches, tail);
		}
	}

	/** Checks grouped oracles and original endpoint witnesses before warming or timing. */
	private void validateRing(final byte[] savedRing, final byte[] savedValidation) throws IOException{
		final StringBuilder inputs=new StringBuilder("ring\tmismatches\tplacement\tA\tB-aligned\tA-quality\tB-quality\n");
		final StringBuilder values=new StringBuilder(VALIDATION_HEADER+"\n");
		final String[] reference=savedValidation==null ? null : new String(savedValidation, StandardCharsets.US_ASCII).split("\n");
		if(reference!=null){require(reference.length==RING_SIZE+1 && reference[0].equals(VALIDATION_HEADER), "Saved baseline score table has wrong schema/row count.");}
		final float correct=normal ? 0.9990f : 0.6019f;
		final float pairCorrect=correct*correct, match=pairCorrect+(1-pairCorrect)*0.25f, error=1-match;
		require(match>error && error>0, "The selected Q30/Q4 grouped oracle requires M>E>0.");
		for(int r=0; r<RING_SIZE; r++){
			final Pair pair=ring[r];
			final String placement=normal ? "rotating" : pair.tail ? "last" : "first";
			final double oracle=Math.pow((double)error/match, pair.mismatches*0.5);
			final float observed=BBMergeOverlapper.probability(pair.a, pair.b, length);
			final boolean unchanged=pair.unchanged();
			require(Float.isFinite(observed) && observed>=0 && observed<=1 && unchanged, "Invalid score or mutated input in ring slot "+r);
			if(baseline && !normal){
				require(observed==(pair.tail ? 0 : 1), "Baseline endpoint witness changed before timing at slot "+r+": "+observed);
			}else{
				// Normal150's two accumulated float products have a conservative relative
				// rounding bound below2e-5. Fixed fallback has only1–3 adverse factors.
				final double tolerance=Math.max(8.0*Math.ulp((float)oracle), (normal ? 0.00002 : 0.000002)*oracle);
				require(observed>0 && Math.abs(observed-oracle)<=tolerance, "Grouped probability oracle mismatch at slot "+r+": "+observed+" vs "+oracle);
			}
			ringBits[r]=Float.floatToRawIntBits(observed);
			if(reference!=null){
				final String[] fields=reference[r+1].split("\t", -1);
				require(fields.length==7 && fields[0].equals(Integer.toString(r)) && fields[1].equals(Integer.toString(pair.mismatches)) &&
					fields[2].equals(placement) && fields[5].matches("[0-9a-f]{8}") && fields[6].equals("true"), "Malformed baseline score row "+r);
				if(normal || baseline){require(ringBits[r]==Integer.parseUnsignedInt(fields[5], 16), "Baseline/ordinary score bits changed at slot "+r);}
			}
			inputs.append(r).append('\t').append(pair.mismatches).append('\t').append(placement)
				.append('\t').append(new String(pair.ab, StandardCharsets.US_ASCII)).append('\t').append(new String(pair.bb, StandardCharsets.US_ASCII))
				.append('\t').append(Arrays.toString(pair.aq)).append('\t').append(Arrays.toString(pair.bq)).append('\n');
			values.append(r).append('\t').append(pair.mismatches).append('\t').append(placement).append('\t').append(oracle)
				.append('\t').append(observed).append('\t').append(hex(ringBits[r])).append('\t').append(unchanged).append('\n');
		}
		final byte[] inputBytes=inputs.toString().getBytes(StandardCharsets.US_ASCII);
		if(savedRing!=null){require(Arrays.equals(savedRing, inputBytes), "Input ring differs from the saved baseline pilot.");}
		Files.write(out.resolve("ring.tsv"), inputBytes, StandardOpenOption.CREATE_NEW);
		write(out.resolve("validation.tsv"), values.toString());
	}

	/** Measures only the public helper, rotating lookup and checksum; no loop allocations. */
	private Result sample(final String stage, final int attempt, final long calls) throws IOException{
		require(calls>0 && calls<=MAX_CALLS && calls%RING_SIZE==0, "Calls must be bounded complete ring cycles: "+calls);
		long checksum=0;
		final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
		for(long i=0; i<calls; i++){
			final Pair pair=ring[(int)i&(RING_SIZE-1)];
			checksum+=Float.floatToRawIntBits(BBMergeOverlapper.probability(pair.a, pair.b, length));
		}
		final long wallNanos=System.nanoTime()-wallStart, cpuNanos=clock.getCurrentThreadCpuTime()-cpuStart;
		final Result result=new Result(calls, checksum, wallNanos, cpuNanos);
		long cycle=0;
		for(final int bits : ringBits){cycle+=bits;}
		final long expected=Math.multiplyExact(cycle, calls/RING_SIZE);
		boolean unchanged=true;
		for(final Pair pair : ring){unchanged&=pair.unchanged();}
		Files.write(out.resolve("measurements.tsv"), (stage+"\t"+attempt+"\t"+(baseline ? "baseline" : "fixed")+"\t"+shape+
			"\t"+calls+"\t"+length+"\t"+checksum+"\t"+expected+"\t"+wallNanos+"\t"+cpuNanos+
			"\t"+(cpuNanos>=MIN_CPU_NANOS)+"\t"+unchanged+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
		require(checksum==expected && unchanged && cpuNanos>0 && wallNanos>0, "Timed checksum/input/clock failure; inspect measurements.tsv.");
		return result;
	}

	/** Reads only a sufficiently long baseline pilot with the same recipe and shape. */
	private long savedCalls(final byte[] bytes){
		final String[] lines=new String(bytes, StandardCharsets.US_ASCII).split("\n");
		require(lines.length==2 && lines[0].equals(COUNTS_HEADER), "Expected the two-line baseline counts.tsv.");
		final String[] f=lines[1].split("\t", -1);
		require(f.length==9 && f[0].equals(RECIPE) && f[1].equals(shape) && f[2].equals("baseline") &&
			f[3].equals("pilot") && Integer.parseInt(f[4])==RING_SIZE && Integer.parseInt(f[6])==length &&
			Integer.parseInt(f[7])==quality, "Baseline recipe/shape/phase/ring differs.");
		final long calls=Long.parseLong(f[5]);
		require(calls>0 && calls<=MAX_CALLS && calls%RING_SIZE==0 && Long.parseLong(f[8])>=MIN_PILOT_NANOS,
			"Baseline pilot must have bounded calls and at least15 CPU seconds.");
		return calls;
	}

	private void writeCounts(final Result r) throws IOException{
		require(baseline && r.cpuNanos>=MIN_PILOT_NANOS, "Only a checked baseline pilot selects matched counts.");
		write(out.resolve("counts.tsv"), COUNTS_HEADER+"\n"+RECIPE+"\t"+shape+"\tbaseline\tpilot\t"+RING_SIZE+
			"\t"+r.calls+"\t"+length+"\t"+quality+"\t"+r.cpuNanos+"\n");
	}

	/** Pilot growth is bounded sixteenfold and targets eighteen thread-CPU seconds. */
	private static long nextCalls(final long calls, final long cpuNanos){
		require(cpuNanos>0, "Pilot sizing requires measured thread CPU time.");
		final long scaled=(long)Math.min(calls*16.0, Math.ceil(calls*(TARGET_CPU_NANOS/(double)cpuNanos)));
		final long next=((Math.max(calls+RING_SIZE, scaled)+RING_SIZE-1)/RING_SIZE)*RING_SIZE;
		require(next<=MAX_CALLS, "Pilot call cap exceeded; review sizing rather than silently expanding work.");
		return next;
	}

	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static String hex(final int bits){return String.format(Locale.ROOT, "%08x", bits);}
	private static void write(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}

	private static final class Pair{
		Pair(final Read a_, final Read b_, final int mismatches_, final boolean tail_){
			a=a_; b=b_; mismatches=mismatches_; tail=tail_;
			ab=a.bases.clone(); bb=b.bases.clone(); aq=a.quality.clone(); bq=b.quality.clone();
			abRef=a.bases; bbRef=b.bases; aqRef=a.quality; bqRef=b.quality;
		}
		boolean unchanged(){return a.bases==abRef && b.bases==bbRef && a.quality==aqRef && b.quality==bqRef &&
			Arrays.equals(ab, a.bases) && Arrays.equals(bb, b.bases) && Arrays.equals(aq, a.quality) && Arrays.equals(bq, b.quality);}
		final Read a, b;
		final int mismatches;
		final boolean tail;
		final byte[] ab, bb, aq, bq, abRef, bbRef, aqRef, bqRef;
	}
	private static final class Result{
		Result(final long calls_, final long checksum_, final long wall_, final long cpu_){calls=calls_; checksum=checksum_; wallNanos=wall_; cpuNanos=cpu_;}
		final long calls, checksum, wallNanos, cpuNanos;
	}

	private final Path out;
	private final String shape;
	private final boolean normal, baseline;
	private final int length, quality;
	private final ThreadMXBean clock;
	private final Pair[] ring=new Pair[RING_SIZE];
	private final int[] ringBits=new int[RING_SIZE];
	private static final int RING_SIZE=64, MAX_ATTEMPTS=6;
	private static final long INITIAL_CALLS=32768, MAX_CALLS=536870912;
	private static final long MIN_CPU_NANOS=10000000000L, MIN_PILOT_NANOS=15000000000L, TARGET_CPU_NANOS=18000000000L;
	private static final String RECIPE="probability-ring-v1";
	private static final String COUNTS_HEADER="recipe\tshape\tphase\tmode\tring\tcalls\tlength\tquality\tpilotCpuNanos";
	private static final String VALIDATION_HEADER="ring\tmismatches\tplacement\toracle\tobserved\tbits\tinputsUnchanged";
}
