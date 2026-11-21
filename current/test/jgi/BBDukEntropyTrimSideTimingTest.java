package test.jgi;

import java.io.IOException;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import jgi.BBDuk;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import tracker.EntropyTracker;

/** Synthetic timing of the high-entropy left-only trim control.
 * A normal parsed owner and one NEW worker reuse an unmapped ACGT/Q40 read.
 * Every complete k1/window4 window has maximal entropy, so cutoff0.1 removes
 * nothing. The proposed side-selection and left-cap branches must precede the
 * zero-count return for this workload to measure their cost after repair.
 * Each real helper call retains its null-BitSet allocation and diagnostic print.
 * Warmup and measurement share an exact primitive bound-handle loop. Full bases,
 * qualities, array identity, metadata, tracker positions and zero direct-call
 * statistics are checked outside the clocks. No fields are seeded and no worker,
 * processing loop, streams or reports run; normal constructor allocators may run.
 * Bash owns calibration and three fresh matched JVM samples per source phase.
 * This measures helper/window/no-trim overhead including allocation and printing,
 * not retained-base copying, mapped reads, pipeline behavior or whole-tool speed.
 * @author Sucrose
 */
public final class BBDukEntropyTrimSideTimingTest{

	/** phase=baseline|fixed mode=pilot|measure warmup=32 calls=POS length=MULTIPLE4 outdir=NEW. */
	public static void main(final String[] args) throws Throwable{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option: "+arg);
		}
		require(options.size()==6 && options.keySet().containsAll(Arrays.asList("phase", "mode", "warmup", "calls", "length", "outdir")), "Use exactly phase, mode, warmup, calls, length and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukentropytrimsidetimingtest.sh with assertions enabled.");
		final String phase=options.get("phase"), mode=options.get("mode");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		final long warmup=positive(options.get("warmup")), calls=positive(options.get("calls"));
		require(warmup==32, "Use the agreed small warmup=32; do not reuse the old100000-call short-kmer warmup.");
		final int length=Math.toIntExact(positive(options.get("length")));
		require(length>=4 && (length&3)==0, "Complete ACGT repetitions require length>=4 and a multiple of4; length="+length);
		final long total=Math.addExact(warmup, calls), windows=(long)length-4+1;
		Math.multiplyExact(total, (long)length); Math.multiplyExact(total, windows);
		final long diagnosticBytes=Math.multiplyExact(total, (long)bytes(diagnostic(length)).length);
		require(diagnosticBytes<=MAX_DIAGNOSTIC_BYTES, "Bound all helper diagnostics including warmup to1MiB; use longer reads and recalibrate rather than suppress output. bytes="+diagnosticBytes);
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported() && clock.isThreadCpuTimeEnabled() && clock.getCurrentThreadCpuTime()>=0, "Require supported, already enabled main-thread CPU time; never substitute process accounting.");
		final BBDukEntropyTrimSideTimingTest test=new BBDukEntropyTrimSideTimingTest(path(options.get("outdir")), length);
		test.run(clock, phase, mode, warmup, calls, total, diagnosticBytes);
	}

	private BBDukEntropyTrimSideTimingTest(final Path out_, final int length_) throws Exception{
		out=out_; length=length_; Files.createDirectory(out);
		bases=new byte[length]; qualities=new byte[length]; final byte[] asciiQuality=new byte[length];
		for(int i=0; i<length; i++){bases[i]=DNA[i&3];}
		Arrays.fill(qualities, (byte)40); Arrays.fill(asciiQuality, (byte)'I');
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("constructor-input.fq"), "@constructor_input\n"+new String(bases, StandardCharsets.US_ASCII)+"\n+\n"+new String(asciiQuality, StandardCharsets.US_ASCII)+"\n");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("literal=ACGA"); args.add("in="+out.resolve("constructor-input.fq"));
		write(out.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		owner=new BBDuk(args.toArray(new String[args.size()]));
		require(field(BBDuk.class, "entropyTrim").getInt(owner)==2 && field(BBDuk.class, "entropyCutoff").getFloat(owner)==0.1f && field(BBDuk.class, "calcEntropy").getBoolean(owner) && field(BBDuk.class, "entropyHighpass").getBoolean(owner) && !field(BBDuk.class, "entropyMask").getBoolean(owner) && !field(BBDuk.class, "entropyMark").getBoolean(owner), "Normal public arguments must select LEFT/high-pass0.1, with no entropy masking or marking.");
		final Class<?> type=Class.forName("jgi.BBDuk$ProcessThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(BBDuk.class, ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, boolean.class); constructor.setAccessible(true);
		worker=(Thread)constructor.newInstance(owner, null, null, null, null, Boolean.FALSE);
		tracker=(EntropyTracker)field(type, "eTrackerT").get(worker);
		require(tracker!=null && tracker.k()==1 && tracker.windowBases()==4 && tracker.cutoff()==0.1f && field(EntropyTracker.class, "highPass").getBoolean(tracker), "Use the normal worker's actual k1/window4/high-pass0.1 tracker; no substitute or field seed.");
		final Method method=type.getDeclaredMethod("trimLowEntropy", Read.class, java.util.BitSet.class, EntropyTracker.class); method.setAccessible(true);
		bound=MethodHandles.lookup().unreflect(method).bindTo(worker);
		require(bound.type().equals(MethodType.methodType(int.class, Read.class, java.util.BitSet.class, EntropyTracker.class)), "Require the exact primitive real helper entry, bound once; no per-call reflection arrays or boxing.");
		workerReads=field(type, "readsEFilteredT"); workerBases=field(type, "basesEFilteredT"); ownerReads=field(BBDuk.class, "readsEFiltered"); ownerBases=field(BBDuk.class, "basesEFiltered");
		read=new Read(bases, qualities, NAME, 0);
		write(out.resolve("reflection.tsv"), "operation\tmember\nCONSTRUCT_ONLY\tBBDuk.ProcessThread(null streams,localArrays=false)\nREAD_ONLY\tBBDuk.entropyTrim,entropyCutoff,calcEntropy,entropyHighpass,entropyMask,entropyMark,readsEFiltered,basesEFiltered\nREAD_ONLY\tProcessThread.eTrackerT,readsEFilteredT,basesEFilteredT;EntropyTracker.highPass\nUNREFLECT_BIND_ONCE\tProcessThread.trimLowEntropy(Read,BitSet,EntropyTracker):int\nFIELD_WRITES\tNONE\n");
		write(out.resolve("controls.tsv"), "mode\tLEFT\npattern\tACGT\nquality\t40\nk\t1\nwindow\t4\ncutoff\t0.1\nhighpass\ttrue\nexpected_removed_per_call\t0\nbitset_argument\tnull\nhelper_allocation\tINCLUDED\nhelper_diagnostics\tRETAINED\nworker_starts\t0\n");
		write(out.resolve("diagnostic.txt"), diagnostic(length));
		write(out.resolve("observations.tsv"), "stage\tcompleted_calls\tread_length\tarrays_identical\tfull_payload_checked\tworker_state\tworker_reads\tworker_bases\towner_reads\towner_bases\ttracker_left\ttracker_right\ttracker_ns\n");
	}

	/** Clock only the same helper loop used for warmup; validate whole state outside it. */
	private void run(final ThreadMXBean clock, final String phase, final String mode, final long warmup, final long calls, final long total, final long diagnosticBytes) throws Throwable{
		final long main=Thread.currentThread().getId(), windows=(long)length-4+1;
		write(out.resolve("config.tsv"), "version\t1\nphase\t"+phase+"\nmode\t"+mode+"\nwarmup\t"+warmup+"\ncalls\t"+calls+"\nlength\t"+length+"\nfull_windows_per_call\t"+windows+"\nmax_helper_diagnostic_bytes\t"+MAX_DIAGNOSTIC_BYTES+"\nexpected_helper_diagnostic_bytes\t"+diagnosticBytes+"\n");
		write(out.resolve("clock.tsv"), "metric\tcurrent_thread_cpu_nanoseconds\nmain_thread_id\t"+main+"\nsupported\ttrue\nenabled\ttrue\n");
		observe("before_warmup", 0);
		final long warmupChecksum=loop(warmup);
		require(warmupChecksum==0, "All high-entropy warmup calls must remove zero bases; checksum="+warmupChecksum);
		observe("after_warmup", warmup);
		final long cpu0=clock.getCurrentThreadCpuTime(), wall0=System.nanoTime();
		final long checksum=loop(calls);
		final long wall=System.nanoTime()-wall0, cpu1=clock.getCurrentThreadCpuTime(), cpu=cpu1-cpu0;
		require(Thread.currentThread().getId()==main && cpu0>=0 && cpu1>cpu0 && wall>0, "Measure on the original main thread with positive CPU/wall deltas.");
		require(checksum==0, "All high-entropy measured calls must remove zero bases; checksum="+checksum);
		observe("after_measure", total);
		write(out.resolve("totals.tsv"), "warmup_checksum\t"+warmupChecksum+"\nmeasured_checksum\t"+checksum+"\ntotal_helper_calls\t"+total+"\nwarmup_bases\t"+Math.multiplyExact(warmup, (long)length)+"\nmeasured_bases\t"+Math.multiplyExact(calls, (long)length)+"\ntotal_bases\t"+Math.multiplyExact(total, (long)length)+"\nwarmup_windows\t"+Math.multiplyExact(warmup, windows)+"\nmeasured_windows\t"+Math.multiplyExact(calls, windows)+"\ntotal_windows\t"+Math.multiplyExact(total, windows)+"\nexpected_helper_diagnostic_lines\t"+total+"\nexpected_helper_diagnostic_bytes\t"+diagnosticBytes+"\n");
		require("pilot".equals(mode) || cpu>=10000000000L, "Measured loop below ten main-thread CPU seconds; recalibrate the common count before matched source cohorts.");
		final String summary="BBDUK_ENTROPY_TRIM_SIDE_TIMING_PASS version=1 phase="+phase+" mode="+mode+" warmup="+warmup+" calls="+calls+" length="+length+" checksum="+checksum+" cpuNanos="+cpu+" wallNanos="+wall+" owners=1 workers=1 workerStarts=0 failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** Helper allocation and printing stay inside each real call; no per-call checks or reset. */
	private long loop(final long calls) throws Throwable{
		assert(calls>0) : "Warmup/measurement require positive calls for tracker position checks; calls="+calls;
		long checksum=0;
		for(long i=0; i<calls; i++){checksum+=(int)bound.invokeExact(read, (java.util.BitSet)null, tracker);}
		return checksum;
	}
	private void observe(final String stage, final long completed) throws Exception{
		require(worker.getState()==Thread.State.NEW && !worker.isAlive(), "Construct the worker only; no start/run/processing or streams at "+stage);
		require(workerReads.getLong(worker)==0 && workerBases.getLong(worker)==0 && ownerReads.getLong(owner)==0 && ownerBases.getLong(owner)==0, "Only the surrounding pipeline aggregates entropy statistics; direct calls leave all four zero at "+stage);
		require(read.bases==bases && read.quality==qualities && read.length()==length, "The no-trim control must preserve both original arrays and length at "+stage);
		for(int i=0; i<length; i++){if(bases[i]!=DNA[i&3] || qualities[i]!=40){require(false, "Check every ACGT/Q40 position outside timing; stage="+stage+" index="+i);}}
		require(NAME.equals(read.id) && read.numericID==0 && read.mate==null && read.obj==null && read.match==null && read.samline==null && read.sites==null && !read.discarded(), "Keep this fresh unmapped read's identity/null metadata/non-discarded state at "+stage);
		final int left=completed==0 ? -4 : length-4, right=completed==0 ? -1 : length-1;
		require(tracker.leftPos()==left && tracker.rightPos()==right && tracker.ns()==0, "Each real call clears and scans every defined base; stage="+stage+" left="+tracker.leftPos()+" right="+tracker.rightPos());
		append(out.resolve("observations.tsv"), stage+"\t"+completed+"\t"+length+"\ttrue\ttrue\t"+worker.getState()+"\t0\t0\t0\t0\t"+left+"\t"+right+"\t0\n");
	}
	private static String diagnostic(final int length){return "Trimming "+NAME+", len "+length+", window 4"+System.lineSeparator();}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field result=type.getDeclaredField(name); result.setAccessible(true); return result;}
	private static long positive(final String text){require(text!=null && text.matches("[0-9]+"), "Use a positive decimal integer: "+text); final long value=Long.parseLong(text); require(value>0, "Counts/length must be positive."); return value;}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0, "Use a nonempty path without equals,tabs,newlines for the legacy constructor."); return Paths.get(text).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.UTF_8);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private final Path out;
	private final int length;
	private final BBDuk owner;
	private final Thread worker;
	private final EntropyTracker tracker;
	private final MethodHandle bound;
	private final Field workerReads, workerBases, ownerReads, ownerBases;
	private final Read read;
	private final byte[] bases, qualities;
	private static final String NAME="measured_entropy_left";
	private static final byte[] DNA={'A', 'C', 'G', 'T'};
	private static final long MAX_DIAGNOSTIC_BYTES=1048576L;
	private static final String[] SETTINGS={"k=4", "mink=-1", "edist=0", "hdist=0", "ktrim=f", "maskmiddle=f", "rcomp=f", "findbestmatch=f", "rename=f", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "array=t", "prealloc=f", "ow=f", "app=f", "entropy=0.1", "entropyk=1", "entropywindow=4", "entropymask=f", "entropymark=f", "entropytrim=l"};
}
