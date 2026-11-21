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
import structures.IntList;

/** Synthetic timing of legacy BBDuk's below-threshold best-match cleanup.
 * A bound handle invokes the actual helper synchronously on one never-started worker.
 * A reused 150-base read has 120 length31 windows and one reference hit; every call
 * returns -1 on both sources. Setup, the independent qualifying control, snapshots
 * and file output are outside the timed loop. No processing pipeline starts.
 * Normal constructor allocators may start and join; no final field is assigned.
 * Bash owns calibration, identical host/runtime settings and three fresh samples
 * per source. Loop CPU uses the main thread's supported, already enabled CPU clock;
 * pilot samples may be short, measured samples require ten CPU seconds.
 * This measures one synthetic branch, not biological accuracy or whole-tool speed.
 * @author Sucrose
 */
public final class BBDukBestMatchTimingTest{

	/** phase=baseline|fixed mode=pilot|measure warmup=N calls=N outdir=NEW. */
	public static void main(final String[] args) throws Throwable{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && !pair[1].isEmpty(), "Use nonempty documented flag=value arguments: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		require(options.size()==5 && options.keySet().containsAll(Arrays.asList("phase", "mode", "warmup", "calls", "outdir")), "Use exactly phase, mode, warmup, calls and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukbestmatchtimingtest.sh with assertions enabled.");
		final String phase=options.get("phase"), mode=options.get("mode");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		final long warmup=positive(options.get("warmup")), calls=positive(options.get("calls"));
		require(warmup<=MAX_TOTAL_CALLS && calls<=MAX_TOTAL_CALLS-warmup, "Warmup plus measured calls must be at most "+MAX_TOTAL_CALLS+"; baseline scratch increments once per call and must not approach int overflow.");
		final Path out=path(options.get("outdir"));
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported() && clock.isThreadCpuTimeEnabled() && clock.getCurrentThreadCpuTime()>=0, "Supported and already enabled main-thread CPU timing is required; do not substitute another metric.");
		final BBDukBestMatchTimingTest test=new BBDukBestMatchTimingTest(out);
		test.run(clock, phase, mode, warmup, calls);
	}

	private BBDukBestMatchTimingTest(final Path out_) throws Exception{
		out=out_; Files.createDirectory(out);
		write(out.resolve("constructor-input.fq"), "@constructor_input\nA\n+\nI\n");
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("literal="+repeat('A', 31)); args.add("in="+out.resolve("constructor-input.fq"));
		write(out.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		owner=new BBDuk(args.toArray(new String[args.size()]));
		final Class<?> type=Class.forName("jgi.BBDuk$ProcessThread");
		constructor=type.getDeclaredConstructor(BBDuk.class, ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class,
			ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, boolean.class); constructor.setAccessible(true);
		method=type.getDeclaredMethod("findBestMatch", Read.class, AbstractKmerTable[].class, int.class); method.setAccessible(true);
		countArray=field(type, "countArray"); idList=field(type, "idList"); countList=field(type, "countList");
		readCounts=field(type, "scaffoldReadCountsT"); baseCounts=field(type, "scaffoldBaseCountsT");
		sets=(AbstractKmerTable[])field(BBDuk.class, "keySets").get(owner);
		@SuppressWarnings("unchecked")
		final ArrayList<String> names=(ArrayList<String>)field(BBDuk.class, "scaffoldNames").get(owner);
		require(sets.length==7 && names.size()==1 && "".equals(names.get(0)), "Retain the ordinary constructor's reserved empty scaffold ID0 and seven tables.");
		long size=0; for(final AbstractKmerTable set : sets){size+=set.size();}
		require(size==0, "No reference loader may populate this direct-helper setup.");
		names.add("target"); sets[(int)(KEY%7)].set(KEY, 1);
		require(sets[(int)(KEY%7)].getValue(KEY)==1, "A^31 with forward unmasked DNA lookup must map key2^62 to scaffold1.");
		setNonFinal(field(BBDuk.class, "storedKmers"), Long.valueOf(1));
		setNonFinal(field(BBDuk.class, "scaffoldReadCounts"), outerReads);
		setNonFinal(field(BBDuk.class, "scaffoldBaseCounts"), outerBases);
		worker=newWorker(); find=bind(worker);
		final String measured=repeat('A', 31)+repeat('C', 119), control=repeat('A', 32)+repeat('C', 118);
		read=newRead(measured, "measured", 0); qualifying=newRead(control, "qualifying", 1);
		expectedBases=bytes(measured); expectedQuality=new byte[150]; Arrays.fill(expectedQuality, (byte)40);
		controlBases=bytes(control); controlQuality=expectedQuality.clone();
		write(out.resolve("inputs.tsv"), "id\tnumeric_id\tbases\tqualities\nmeasured\t0\t"+measured+"\t"+repeat('I', 150)+"\nqualifying\t1\t"+control+"\t"+repeat('I', 150)+"\n");
		write(out.resolve("seed.tsv"), "key\t"+KEY+"\nshard\t"+(KEY%7)+"\nscaffold_id\t1\nstored_kmers\t1\nscaffold_name_0\t"+names.get(0)+"\nscaffold_name_1\t"+names.get(1)+"\nscaffold_counter_length\t2\nfinal_field_assignments\t0\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nREAD_CONTENTS\tBBDuk.keySets,BBDuk.scaffoldNames\nWRITE_NONFINAL\tBBDuk.storedKmers,BBDuk.scaffoldReadCounts,BBDuk.scaffoldBaseCounts\nCONSTRUCT_ONLY\tBBDuk.ProcessThread\nUNREFLECT_BIND_ONCE\tProcessThread.findBestMatch\nOBSERVE\tProcessThread.countArray,idList,countList,scaffoldReadCountsT,scaffoldBaseCountsT\n");
		write(out.resolve("observations.tsv"), "stage\tcounts\tids\thits\tscaffold_reads\tscaffold_bases\tworker_state\n");
	}

	/** Only loop() is measured; validation observes scratch without imposing a cleanup strategy. */
	private void run(final ThreadMXBean clock, final String phase, final String mode, final long warmup, final long calls) throws Throwable{
		final long main=Thread.currentThread().getId();
		write(out.resolve("config.tsv"), "version\t1\nphase\t"+phase+"\nmode\t"+mode+"\nwarmup\t"+warmup+"\ncalls\t"+calls+"\nread_length\t150\nk\t31\nwindows_per_call\t120\nhits_per_call\t1\nthreshold\t1\nmax_total_calls\t"+MAX_TOTAL_CALLS+"\n");
		write(out.resolve("clock.tsv"), "metric\tcurrent_thread_cpu_nanoseconds\nmain_thread_id\t"+main+"\nsupported\ttrue\nenabled\ttrue\n");
		qualifyingControl(); observe("before_warmup", worker); payloads();
		final long warmupChecksum=loop(warmup);
		require(warmupChecksum==-warmup, "Every warmup return must be -1 for the literal single-hit read; checksum="+warmupChecksum);
		observe("after_warmup", worker); payloads();
		final long cpu0=clock.getCurrentThreadCpuTime(), wall0=System.nanoTime();
		final long checksum=loop(calls);
		final long wall=System.nanoTime()-wall0, cpu1=clock.getCurrentThreadCpuTime(), cpu=cpu1-cpu0;
		require(Thread.currentThread().getId()==main && cpu0>=0 && cpu1>cpu0 && wall>0, "Loop must complete on the original main thread with positive CPU and wall-clock deltas.");
		require(checksum==-calls, "Every measured return must be -1 for the literal single-hit read; checksum="+checksum);
		observe("after_measure", worker); payloads();
		require(outerReads.get(0)==0 && outerReads.get(1)==0 && outerBases.get(0)==0 && outerBases.get(1)==0, "Direct local-array calls must not aggregate owner counts or run processing workers.");
		require("pilot".equals(mode) || cpu>=10000000000L, "Measured loop is below ten main-thread CPU seconds; enlarge the frozen count and redo BOTH phases.");
		final String summary="BBDUK_BEST_MATCH_TIMING_PASS version=1 phase="+phase+" mode="+mode+" warmup="+warmup+" calls="+calls+
			" windows="+Math.multiplyExact(calls, 120L)+" checksum="+checksum+" cpuNanos="+cpu+" wallNanos="+wall+" workers=2 processingStarts=0 failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** Same hot loop for warmup and measurement: no reflection arrays, boxing, snapshots or output. */
	private long loop(final long calls) throws Throwable{
		assert(calls>0 && calls<=MAX_TOTAL_CALLS) : "Positive bounded call counts protect the baseline int counter and long checksum; calls="+calls;
		long checksum=0;
		for(long i=0; i<calls; i++){checksum+=(int)find.invokeExact(read, sets, 1);}
		return checksum;
	}

	/** A different NEW worker owns the qualifying call, never the stale baseline timing scratch. */
	private void qualifyingControl() throws Throwable{
		final Thread control=newWorker(); final MethodHandle handle=bind(control);
		final int result=(int)handle.invokeExact(qualifying, sets, 1);
		final long[] reads=(long[])readCounts.get(control), bases=(long[])baseCounts.get(control);
		require(result==1 && Arrays.equals(reads, new long[]{0, 1}) && Arrays.equals(bases, new long[]{0, 150}), "Independent two-hit qualifying control must select scaffold1/read1/base150.");
		observe("qualifying_control", control);
		write(out.resolve("controls.tsv"), "id\treturn\tscaffold_reads\tscaffold_bases\tworker_state\nqualifying\t1\t[0, 1]\t[0, 150]\tNEW\n");
	}
	private Thread newWorker() throws Exception{
		final Thread result=(Thread)constructor.newInstance(owner, null, null, null, null, Boolean.TRUE);
		require(result.getState()==Thread.State.NEW && !result.isAlive(), "Processing worker is construction-only and must remain NEW.");
		return result;
	}
	private MethodHandle bind(final Thread target) throws IllegalAccessException{
		final MethodHandle result=MethodHandles.lookup().unreflect(method).bindTo(target);
		require(result.type().equals(MethodType.methodType(int.class, Read.class, AbstractKmerTable[].class, int.class)), "Bound helper must retain the exact primitive-return signature for allocation-free invokeExact.");
		return result;
	}
	private void observe(final String stage, final Thread target) throws Exception{
		final int[] counts=(int[])countArray.get(target);
		final long[] reads=(long[])readCounts.get(target), bases=(long[])baseCounts.get(target);
		require(target.getState()==Thread.State.NEW && !target.isAlive() && counts.length==2 && reads.length==2 && bases.length==2, "Observe only a NEW worker with two-scaffold local scratch.");
		require(target!=worker || (reads[0]==0 && reads[1]==0 && bases[0]==0 && bases[1]==0), "Below-threshold timing calls must not update scaffold statistics.");
		append(out.resolve("observations.tsv"), stage+"\t"+Arrays.toString(counts)+"\t"+list((IntList)idList.get(target))+"\t"+list((IntList)countList.get(target))+
			"\t"+Arrays.toString(reads)+"\t"+Arrays.toString(bases)+"\t"+target.getState()+"\n");
	}
	private void payloads(){
		require("measured".equals(read.id) && read.numericID==0 && read.mate==null && Arrays.equals(read.bases, expectedBases) && Arrays.equals(read.quality, expectedQuality), "Timing helper must preserve measured read identity, bases and Q40 qualities.");
		require("qualifying".equals(qualifying.id) && qualifying.numericID==1 && qualifying.mate==null && Arrays.equals(qualifying.bases, controlBases) && Arrays.equals(qualifying.quality, controlQuality), "Timing helper must preserve the independent qualifying payload.");
	}
	private void setNonFinal(final Field field, final Object value) throws IllegalAccessException{
		require(!Modifier.isFinal(field.getModifiers()) && Arrays.asList("storedKmers", "scaffoldReadCounts", "scaffoldBaseCounts").contains(field.getName()), "Only three named NONFINAL owner fields may be seeded: "+field.getName());
		field.set(owner, value);
	}
	private static Read newRead(final String bases, final String id, final long numeric){final byte[] q=new byte[bases.length()]; Arrays.fill(q, (byte)40); return new Read(bytes(bases), q, id, numeric);}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field field=type.getDeclaredField(name); field.setAccessible(true); return field;}
	private static String list(final IntList values){final StringBuilder out=new StringBuilder("["); for(int i=0; i<values.size; i++){if(i>0){out.append(',');}out.append(values.get(i));} return out.append(']').toString();}
	private static String repeat(final char value, final int count){final char[] text=new char[count]; Arrays.fill(text, value); return new String(text);}
	private static long positive(final String text){require(text!=null && text.matches("[0-9]+"), "Expected positive call count: "+text); final long value=Long.parseLong(text); require(value>0, "Call counts must be positive."); return value;}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0 && text.indexOf('\n')<0 && text.indexOf('\t')<0, "Set a nonempty path without equals, tabs or newlines."); return Paths.get(text).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private final Path out;
	private final BBDuk owner;
	private final Constructor<?> constructor;
	private final Method method;
	private final MethodHandle find;
	private final Thread worker;
	private final Read read, qualifying;
	private final byte[] expectedBases, expectedQuality, controlBases, controlQuality;
	private final AbstractKmerTable[] sets;
	private final Field countArray, idList, countList, readCounts, baseCounts;
	private final AtomicLongArray outerReads=new AtomicLongArray(2), outerBases=new AtomicLongArray(2);
	private static final long KEY=1L<<62, MAX_TOTAL_CALLS=1000000000L;
	private static final String[] SETTINGS={"k=31", "maskmiddle=f", "rcomp=f", "qskip=1", "speed=0", "qhdist=0", "hdist=0", "edist=0", "useshortkmers=f",
		"findbestmatch=t", "rename=f", "skipr1=f", "skipr2=f", "restrictleft=0", "restrictright=0", "maxbadkmers=1", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "prealloc=f"};
}
