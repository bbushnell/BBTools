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
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLongArray;

import jgi.BBDuk;
import kmer.AbstractKmerTable;
import parse.Parse;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Synthetic timing of legacy BBDuk's histogram-enabled covered-base bookkeeping.
 * A reused 150-base all-A read has 120 matching length31 windows. Threshold75 is
 * crossed at match45, leaving76 baseline scaffold updates versus one after repair;
 * both versions return full coverage150 and increment histogram bin150 once per call.
 * An independent no-histogram worker preserves early coverage75/read1/base150.
 * Two ordinary owners select histogram state through outduk arguments. Processing
 * workers remain NEW; only joined constructor table allocators may start and finish.
 * A bound primitive handle invokes the real helper on the main thread. Setup,
 * validation, snapshots and output stay outside the identical warmup/measurement loop.
 * This measures a cache-hot synthetic bookkeeping case, not whole-tool or biological
 * performance. Bash owns explicit calibration, matched counts and fresh JVM samples.
 * @author Sucrose
 */
public final class BBDukCoveredBasesTimingTest{

	/** phase=baseline|fixed mode=pilot|measure warmup=N calls=N outdir=NEW. */
	public static void main(final String[] args) throws Throwable{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && !pair[1].isEmpty(), "Use documented nonempty flag=value arguments: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		require(options.size()==5 && options.keySet().containsAll(Arrays.asList("phase", "mode", "warmup", "calls", "outdir")), "Use exactly phase, mode, warmup, calls and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukcoveredbasestimingtest.sh with assertions enabled.");
		final String phase=options.get("phase"), mode=options.get("mode");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		final long warmup=positive(options.get("warmup")), calls=positive(options.get("calls"));
		require(warmup<=MAX_TOTAL_CALLS && calls<=MAX_TOTAL_CALLS-warmup, "Warmup plus measured calls must be at most "+MAX_TOTAL_CALLS+" for bounded counter/checksum arithmetic.");
		final long total=Math.addExact(warmup, calls);
		Math.multiplyExact(Math.multiplyExact(total, 76L), 150L);// Largest expected counter: baseline full-read base credit.
		Math.multiplyExact(calls, 150L); Math.multiplyExact(calls, 120L);
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported() && clock.isThreadCpuTimeEnabled() && clock.getCurrentThreadCpuTime()>=0, "Supported and already enabled main-thread CPU timing is required; do not substitute another metric.");
		final BBDukCoveredBasesTimingTest test=new BBDukCoveredBasesTimingTest(path(options.get("outdir")));
		test.run(clock, phase, mode, warmup, calls, total);
	}

	private BBDukCoveredBasesTimingTest(final Path out_) throws Exception{
		out=out_; Files.createDirectory(out);
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		final Class<?> type=Class.forName("jgi.BBDuk$ProcessThread");
		constructor=type.getDeclaredConstructor(BBDuk.class, ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class,
			ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, boolean.class); constructor.setAccessible(true);
		method=type.getDeclaredMethod("countCoveredBases", Read.class, AbstractKmerTable[].class, int.class); method.setAccessible(true);
		readCounts=field(type, "scaffoldReadCountsT"); baseCounts=field(type, "scaffoldBaseCountsT"); histogram=field(type, "hitCountsT"); ownerHistogram=field(BBDuk.class, "hitCounts");
		on=newOwner(true); off=newOwner(false);
		worker=newWorker(on); covered=bind(worker);
		expectedBases=bytes(repeat('A', 150)); expectedQuality=new byte[150]; Arrays.fill(expectedQuality, (byte)40);
		read=new Read(expectedBases.clone(), expectedQuality.clone(), "measured", 0);
		control=new Read(expectedBases.clone(), expectedQuality.clone(), "no_histogram", 1);
		write(out.resolve("inputs.tsv"), "id\tnumeric_id\tbases\tqualities\nmeasured\t0\t"+repeat('A', 150)+"\t"+repeat('I', 150)+"\nno_histogram\t1\t"+repeat('A', 150)+"\t"+repeat('I', 150)+"\n");
		write(out.resolve("oracles.tsv"), "mode\tcoverage\treads_per_call\tbases_per_call\thistogram\nmeasured_baseline\t150\t76\t11400\tbin150=1,others0\nmeasured_fixed\t150\t1\t150\tbin150=1,others0\nno_histogram_control\t75\t1\t150\tnull\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nREAD_CONTENTS\tBBDuk.keySets,BBDuk.scaffoldNames\nWRITE_NONFINAL\tBBDuk.storedKmers,BBDuk.scaffoldReadCounts,BBDuk.scaffoldBaseCounts\nCONSTRUCT_ONLY\tBBDuk.ProcessThread\nUNREFLECT_BIND_ONCE\tProcessThread.countCoveredBases\nOBSERVE\tProcessThread.scaffoldReadCountsT,scaffoldBaseCountsT,hitCountsT,BBDuk.hitCounts\n");
		write(out.resolve("observations.tsv"), "stage\tscaffold_reads\tscaffold_bases\thistogram\tworker_state\n");
	}

	/** Histogram allocation uses ordinary outduk selection; only three named nonfinal fields are seeded. */
	private Owner newOwner(final boolean enabled) throws Exception{
		final Path directory=out.resolve(enabled ? "histogram_on" : "histogram_off"); Files.createDirectory(directory);
		write(directory.resolve("constructor-input.fq"), "@constructor_input\nA\n+\nI\n");
		final Path report=directory.resolve("unused-duk.txt");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("literal="+repeat('A', 31)); args.add("in="+directory.resolve("constructor-input.fq")); if(enabled){args.add("outduk="+report);}
		write(directory.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		final BBDuk tool=new BBDuk(args.toArray(new String[args.size()]));
		final AbstractKmerTable[] sets=(AbstractKmerTable[])field(BBDuk.class, "keySets").get(tool);
		@SuppressWarnings("unchecked")
		final ArrayList<String> names=(ArrayList<String>)field(BBDuk.class, "scaffoldNames").get(tool);
		require(sets.length==7 && names.size()==1 && "".equals(names.get(0)), "Retain the normal reserved scaffold ID0 and seven empty tables.");
		long size=0; for(final AbstractKmerTable set : sets){size+=set.size();}
		require(size==0, "No reference loader may populate this direct-helper setup.");
		names.add("target"); sets[(int)(KEY%7)].set(KEY, 1);
		require(sets[(int)(KEY%7)].getValue(KEY)==1, "Forward unmasked A^31 uses key2^62 in shard4, scaffold1.");
		final Owner owner=new Owner(tool, sets, enabled, report);
		setNonFinal(tool, field(BBDuk.class, "storedKmers"), Long.valueOf(1));
		setNonFinal(tool, field(BBDuk.class, "scaffoldReadCounts"), owner.reads);
		setNonFinal(tool, field(BBDuk.class, "scaffoldBaseCounts"), owner.bases);
		ownerUnchanged(owner);
		write(directory.resolve("seed.tsv"), "key\t"+KEY+"\nshard\t"+(KEY%7)+"\nscaffold_id\t1\nstored_kmers\t1\nscaffold_name_0\t"+names.get(0)+"\nscaffold_name_1\t"+names.get(1)+"\nscaffold_counter_length\t2\nhistogram_enabled\t"+enabled+"\nhistogram_length\t"+(enabled ? 1001 : 0)+"\nfinal_or_histogram_assignments\t0\n");
		return owner;
	}

	/** Warmup and measured calls share the identical loop; all assertions/output surround it. */
	private void run(final ThreadMXBean clock, final String phase, final String mode, final long warmup, final long calls, final long total) throws Throwable{
		final boolean baseline=phase.equals("baseline");
		final long main=Thread.currentThread().getId(), expectedWarmup=Math.multiplyExact(warmup, 150L), expectedChecksum=Math.multiplyExact(calls, 150L);
		final long warmupReads=expectedReads(warmup, baseline), totalReads=expectedReads(total, baseline);
		final long warmupBases=Math.multiplyExact(warmupReads, 150L), totalBases=Math.multiplyExact(totalReads, 150L);
		write(out.resolve("config.tsv"), "version\t1\nphase\t"+phase+"\nmode\t"+mode+"\nwarmup\t"+warmup+"\ncalls\t"+calls+"\nread_length\t150\nk\t31\nwindows_per_call\t120\ncovered_bases_per_call\t150\nthreshold\t75\nmax_total_calls\t"+MAX_TOTAL_CALLS+"\n");
		write(out.resolve("clock.tsv"), "metric\tcurrent_thread_cpu_nanoseconds\nmain_thread_id\t"+main+"\nsupported\ttrue\nenabled\ttrue\n");
		noHistogramControl(); observe("before_warmup", on, worker, 0, 0, 0); payloads();
		final long warmupChecksum=loop(warmup);
		require(warmupChecksum==expectedWarmup, "Warmup must return full literal coverage150 per call; checksum="+warmupChecksum);
		observe("after_warmup", on, worker, warmupReads, warmupBases, warmup); payloads();
		final long cpu0=clock.getCurrentThreadCpuTime(), wall0=System.nanoTime();
		final long checksum=loop(calls);
		final long wall=System.nanoTime()-wall0, cpu1=clock.getCurrentThreadCpuTime(), cpu=cpu1-cpu0;
		require(Thread.currentThread().getId()==main && cpu0>=0 && cpu1>cpu0 && wall>0, "Complete the loop on the original main thread with positive CPU and wall-clock deltas.");
		require(checksum==expectedChecksum, "Every measured call must return full coverage150; checksum="+checksum);
		observe("after_measure", on, worker, totalReads, totalBases, total); payloads(); ownerUnchanged(off);
		write(out.resolve("totals.tsv"), "warmup_checksum\t"+warmupChecksum+"\nmeasured_checksum\t"+checksum+"\ntotal_calls\t"+total+"\nexpected_scaffold_reads\t"+totalReads+"\nexpected_scaffold_bases\t"+totalBases+"\nexpected_histogram_bin150\t"+total+"\n");
		require("pilot".equals(mode) || cpu>=10000000000L, "Measured loop below ten main-thread CPU seconds; enlarge the common count and redo BOTH phases.");
		final String summary="BBDUK_COVERED_BASES_TIMING_PASS version=1 phase="+phase+" mode="+mode+" warmup="+warmup+" calls="+calls+
			" windows="+Math.multiplyExact(calls, 120L)+" checksum="+checksum+" cpuNanos="+cpu+" wallNanos="+wall+" owners=2 workers=2 processingStarts=0 failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** No boxing, reflection arrays, per-call allocations, snapshots, output or counter resets. */
	private long loop(final long calls) throws Throwable{
		assert(calls>0 && calls<=MAX_TOTAL_CALLS) : "Bounded positive calls protect long coverage/statistics arithmetic; calls="+calls;
		long checksum=0;
		for(long i=0; i<calls; i++){checksum+=(int)covered.invokeExact(read, on.sets, 75);}
		return checksum;
	}
	private void noHistogramControl() throws Throwable{
		final Thread target=newWorker(off); final MethodHandle handle=bind(target);
		final int result=(int)handle.invokeExact(control, off.sets, 75);
		require(result==75, "No-histogram control intentionally returns early coverage75, not full coverage150.");
		observe("no_histogram_control", off, target, 1, 150, 0); payloads();
		write(out.resolve("controls.tsv"), "id\tcoverage\tscaffold_reads\tscaffold_bases\thistogram\tworker_state\nno_histogram\t75\t[0, 1]\t[0, 150]\tnull\tNEW\n");
	}
	private void observe(final String stage, final Owner owner, final Thread target, final long expectedReads, final long expectedBases, final long completedCalls) throws Exception{
		final long[] reads=(long[])readCounts.get(target), bases=(long[])baseCounts.get(target), counts=(long[])histogram.get(target);
		require(target.getState()==Thread.State.NEW && !target.isAlive() && reads!=null && reads.length==2 && bases!=null && bases.length==2, "Observe only a never-started worker with two-scaffold local statistics.");
		require(Arrays.equals(reads, new long[]{0, expectedReads}) && Arrays.equals(bases, new long[]{0, expectedBases}), "Literal phase-specific scaffold counter mismatch at "+stage);
		if(owner.enabled){
			final long[] expected=new long[1001]; expected[150]=completedCalls;// Independent oracle; never write the actual histogram.
			require(Arrays.equals(counts, expected), "All1001 histogram bins must count full coverage150 once per call at "+stage);
		}else{require(counts==null, "No-histogram control must retain null worker histogram.");}
		ownerUnchanged(owner);
		append(out.resolve("observations.tsv"), stage+"\t"+Arrays.toString(reads)+"\t"+Arrays.toString(bases)+"\t"+Arrays.toString(counts)+"\t"+target.getState()+"\n");
	}
	private void ownerUnchanged(final Owner owner) throws IllegalAccessException{
		require(owner.reads.get(0)==0 && owner.reads.get(1)==0 && owner.bases.get(0)==0 && owner.bases.get(1)==0, "No pipeline aggregation may update owner atomic statistics.");
		final long[] counts=(long[])ownerHistogram.get(owner.tool);
		require(owner.enabled ? counts!=null && counts.length==1001 && allZero(counts) : counts==null, "Ordinary outduk selection must leave the owner histogram allocated/zero or null.");
		require(!Files.exists(owner.report), "Never create or run a DUK report writer in this direct-helper timing fixture.");
	}
	private Thread newWorker(final Owner owner) throws Exception{final Thread target=(Thread)constructor.newInstance(owner.tool, null, null, null, null, Boolean.TRUE); require(target.getState()==Thread.State.NEW && !target.isAlive(), "Construct only a NEW processing worker; never start it."); return target;}
	private MethodHandle bind(final Thread target) throws IllegalAccessException{final MethodHandle result=MethodHandles.lookup().unreflect(method).bindTo(target); require(result.type().equals(MethodType.methodType(int.class, Read.class, AbstractKmerTable[].class, int.class)), "Bind the exact primitive helper signature for allocation-free invokeExact."); return result;}
	private void payloads(){require("measured".equals(read.id) && read.numericID==0 && read.mate==null && Arrays.equals(read.bases, expectedBases) && Arrays.equals(read.quality, expectedQuality), "Preserve measured literal identity/bases/Q40."); require("no_histogram".equals(control.id) && control.numericID==1 && control.mate==null && Arrays.equals(control.bases, expectedBases) && Arrays.equals(control.quality, expectedQuality), "Preserve independent no-histogram literal identity/bases/Q40.");}
	private void setNonFinal(final BBDuk tool, final Field field, final Object value) throws IllegalAccessException{require(!Modifier.isFinal(field.getModifiers()) && Arrays.asList("storedKmers", "scaffoldReadCounts", "scaffoldBaseCounts").contains(field.getName()), "Only three named NONFINAL owner fields may be seeded: "+field.getName()); field.set(tool, value);}
	private static long expectedReads(final long calls, final boolean baseline){return Math.multiplyExact(calls, baseline ? 76L : 1L);}
	private static boolean allZero(final long[] values){for(final long value : values){if(value!=0){return false;}}return true;}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field field=type.getDeclaredField(name); field.setAccessible(true); return field;}
	private static String repeat(final char base, final int length){final char[] bases=new char[length]; Arrays.fill(bases, base); return new String(bases);}
	private static long positive(final String text){require(text!=null && text.matches("[0-9]+"), "Expected positive call count: "+text); final long value=Long.parseLong(text); require(value>0, "Call counts must be positive."); return value;}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0 && text.indexOf('\n')<0 && text.indexOf('\t')<0, "Use a nonempty path without equals, tabs or newlines for the legacy constructor."); return Paths.get(text).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final class Owner{
		Owner(final BBDuk tool_, final AbstractKmerTable[] sets_, final boolean enabled_, final Path report_){tool=tool_; sets=sets_; enabled=enabled_; report=report_;}
		final BBDuk tool; final AbstractKmerTable[] sets; final boolean enabled; final Path report;
		final AtomicLongArray reads=new AtomicLongArray(2), bases=new AtomicLongArray(2);
	}
	private final Path out;
	private final Constructor<?> constructor;
	private final Method method;
	private final MethodHandle covered;
	private final Thread worker;
	private final Owner on, off;
	private final Read read, control;
	private final byte[] expectedBases, expectedQuality;
	private final Field readCounts, baseCounts, histogram, ownerHistogram;
	private static final long KEY=1L<<62, MAX_TOTAL_CALLS=1000000000L;
	private static final String[] SETTINGS={"k=31", "maskmiddle=f", "rcomp=f", "qskip=1", "speed=0", "qhdist=0", "hdist=0", "edist=0", "useshortkmers=f", "mcf=0.5",
		"findbestmatch=f", "rename=f", "skipr1=f", "skipr2=f", "restrictleft=0", "restrictright=0", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "prealloc=f", "ow=f", "app=f"};
}
