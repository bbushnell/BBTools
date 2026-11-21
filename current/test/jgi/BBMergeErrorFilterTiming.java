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
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jgi.BBMerge;
import jgi.BBMergeOverlapper;
import parse.Parse;
import stream.Read;
import stream.Streamer;
import stream.Writer;

/** One non-neural worker timing, using the accepted error-filter fixture's inputs.
 * Selected retained correctness classes must precede the shared timing build, then
 * current protected dependencies on the classpath. Actual CodeSource locations and
 * loaded byte pins are checked before/after; the runner verifies accepted inventories. No production class
 * is installed or rewritten. The strict case changes from an erroneous acceptance
 * to rejection; report it separately from the unchanged perfect-overlap control.
 * The loop contains one typed worker call and return-checksum addition. Reflection,
 * allocation, input checks, provenance and file I/O remain outside both clocks.
 * Three repetitions per case/phase are orchestrated by the external runner; a
 * successful baseline pilot can serve as repetition one. This is not CLI throughput,
 * NN timing, default accuracy, or a custom/training runtime test.
 * @author Sucrose
 */
public final class BBMergeErrorFilterTiming{

	/** mode=pilot|measure, phase=baseline|fixed, case=perfect|strict, root, classes,
	 * timingclasses, new outdir; measure requires calls=&lt;saved baseline count&gt;. */
	public static void main(final String[] args) throws Throwable{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "case", "root", "classes", "timingclasses", "outdir", "calls").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate timing option: "+pair[0]);
		}
		final String mode=options.get("mode"), phase=options.get("phase"), name=options.get("case");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require("perfect".equals(name) || "strict".equals(name), "Only the accepted perfect/strict non-NN cases are timed.");
		final boolean pilot=mode.equals("pilot");
		require(!pilot || phase.equals("baseline"), "Only baseline pilots may select shared call counts.");
		require(pilot ? !options.containsKey("calls") : options.containsKey("calls"), "Pilot selects calls; measure requires calls=<saved baseline count>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergeerrorfiltertiming.sh with assertions enabled.");
		final BBMergeErrorFilterTiming test=new BBMergeErrorFilterTiming(options, phase, name);
		long calls=pilot ? INITIAL_CALLS : Long.parseLong(options.get("calls"));
		test.sample("warmup", 0, WARMUP_CALLS);
		Result result=null;
		for(int attempt=1; attempt<=(pilot ? MAX_ATTEMPTS : 1); attempt++){
			result=test.sample(mode, attempt, calls);
			if(!pilot || result.cpu>=TARGET_CPU_NANOS){break;}
			if(attempt<MAX_ATTEMPTS){
				final long next=(long)Math.ceil(calls*(TARGET_CPU_NANOS/(double)result.cpu)*1.05);
				calls=Math.min(MAX_CALLS, Math.max(calls+1, next));
			}
		}
		require(result!=null, "A timing invocation must retain its measured result.");
		test.verifyProvenance("after");
		final boolean longEnough=result.cpu>=MIN_CPU_NANOS && (!pilot || result.cpu>=TARGET_CPU_NANOS);
		write(test.out.resolve("metrics.tsv"), "case\tphase\tmode\tcalls\texpected_return\tchecksum\tcpu_ns\twall_ns\n"+
			name+"\t"+phase+"\t"+mode+"\t"+result.calls+"\t"+test.expectedReturn+"\t"+result.checksum+"\t"+result.cpu+"\t"+result.wall+"\n");
		final String summary="BBMERGE_ERROR_FILTER_TIMING_"+(longEnough ? "PASS" : "SHORT")+" mode="+mode+" phase="+phase+" case="+name+
			" calls="+result.calls+" expectedReturn="+test.expectedReturn+" checksum="+result.checksum+" cpuNanos="+result.cpu+" wallNanos="+result.wall+
			" inputsUnchanged="+test.unchanged()+" kernel=one-worker-call-plus-return-checksum\n";
		write(test.out.resolve("summary.txt"), summary); System.out.print(summary);
		require(longEnough, "Measured worker region is too short; adjust a new baseline pilot, never fixed-only counts.");
	}

	/** Resolves fixture setup and typed worker invocation before the measured region. */
	@SuppressWarnings("unchecked")
	private BBMergeErrorFilterTiming(final Map<String, String> options, final String phase_, final String name_) throws Throwable{
		phase=phase_; name=name_; expectedReturn=phase.equals("fixed") && name.equals("strict") ? BBMerge.RET_AMBIG : 100;
		root=existing(options.get("root")); build=existing(options.get("classes")); timingBuild=existing(options.get("timingclasses"));
		out=path(options.get("outdir"));
		require(Files.isDirectory(root) && Files.isDirectory(build) && Files.isDirectory(timingBuild) && !build.equals(timingBuild), "Use distinct existing selected/shared build directories and one runtime root.");
		Files.createDirectory(out);
		verifyProvenance("before");
		clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Current-thread CPU time is required, not process wall time.");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		final Class<?> fixture=Class.forName("test.jgi.BBMergeErrorFilterTest"), specClass=Class.forName("test.jgi.BBMergeErrorFilterTest$Spec");
		final Object spec=invoke(fixture, null, "spec", new Class<?>[]{String.class}, name);
		a=(Read)invoke(specClass, spec, "first", new Class<?>[0]);
		final Read rawB=(Read)invoke(specClass, spec, "second", new Class<?>[0]);
		b=(Read)invoke(specClass, spec, "alignedSecond", new Class<?>[0]);
		a.mate=b; b.mate=a; b.setPairnum(1);
		invoke(fixture, null, "writeFastq", new Class<?>[]{Path.class, Read.class}, out.resolve("A.fq"), a);
		invoke(fixture, null, "writeFastq", new Class<?>[]{Path.class, Read.class}, out.resolve("B.fq"), rawB);
		ab=a.bases.clone(); aq=a.quality.clone(); bb=b.bases.clone(); bq=b.quality.clone();
		write(out.resolve("inputs.tsv"), "read\tid\tbases\tqualities\n"+record("A", a)+record("alignedB", b));
		final List<String> arguments=(List<String>)invoke(fixture, null, "arguments", new Class<?>[]{Path.class, Path.class, specClass}, root, out, spec);
		final BBMerge merge=new BBMerge(arguments.toArray(new String[0]));
		require(Runtime.getRuntime().availableProcessors()==1 && !shared.Shared.SIMD, "Use one active processor and SIMD disabled, as in the accepted fixture.");
		require(!bool(merge, "MAKE_VECTOR") && !bool(merge, "TAG_CUSTOM") && !bool(merge, "parseCustom") && !bool(merge, "ONLY_OUTPUT_FALSE_POSITIVE"), "Timing is only the ordinary non-NN path.");
		require(bool(merge, "useEfilter") && number(merge, "efilterRatio")==.06f && number(merge, "efilterOffset")==.05f && number(merge, "pfilterRatio")==0 && bool(merge, "useQuality"), "Actual strict quality-filter gates differ from the accepted fixture.");
		require(!bool(merge, "overlapUsingQuality") && bool(merge, "overlapWithoutQuality") && bool(merge, "useRatioMode") && !bool(merge, "useFlatMode") && !bool(merge, "useMEEfilter") && !bool(merge, "useEntropy"), "An alternate candidate/prefilter path invalidates this bounded timing.");
		final Class<?> workerClass=Class.forName("jgi.BBMerge$MateThread");
		final Constructor<?> constructor=workerClass.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, Writer.class, boolean.class, boolean.class);
		constructor.setAccessible(true);
		final Object worker=constructor.newInstance(merge, null, null, null, null, true, false);
		require(((Thread)worker).getState()==Thread.State.NEW, "Timing must not start a worker thread or stream.");
		final Method method=workerClass.getDeclaredMethod("mateByOverlap", Read.class, Read.class); method.setAccessible(true);
		call=MethodHandles.lookup().unreflect(method).bindTo(worker).asType(MethodType.methodType(int.class, Read.class, Read.class));
		final int returned=(int)call.invokeExact(a, b);
		final int[] vector=(int[])value(worker, "rvector");
		final float expected=BBMergeOverlapper.expectedMismatches(a, b, 100);
		require(vector[2]==(name.equals("strict") ? 1 : 0) && vector[4]==0 && Math.abs(expected-(name.equals("strict") ? 1.98109 : 1.99))<.0005,
			"Actual mismatch/candidate quality prerequisite differs from the accepted truth case.");
		require(returned==expectedReturn && unchanged(), "Preflight worker result or input preservation differs from phase oracle.");
		write(out.resolve("preflight.tsv"), "case\tphase\texpectedReturn\treturned\tbad\thelperExpected\n"+name+"\t"+phase+"\t"+expectedReturn+"\t"+returned+"\t"+vector[2]+"\t"+expected+"\n");
		write(out.resolve("measurements.tsv"), "stage\tattempt\tphase\tcase\tcalls\texpectedReturn\tchecksum\texpectedChecksum\twallNanos\tcpuNanos\tatLeastTenCpuSeconds\tinputsUnchanged\n");
		write(out.resolve("runtime.txt"), "java.version="+System.getProperty("java.version")+"\njava.vendor="+System.getProperty("java.vendor")+
			"\njava.vm.name="+System.getProperty("java.vm.name")+"\nactiveProcessors="+Runtime.getRuntime().availableProcessors()+"\nsimd="+shared.Shared.SIMD+
			"\nrecipe="+RECIPE+"\nkernel=one-worker-call-plus-return-checksum\n");
	}

	/** Exact current-thread region; no reflective invocation, parsing or output in loop. */
	private Result sample(final String stage, final int attempt, final long calls) throws Throwable{
		require(calls>0 && calls<=MAX_CALLS, "Calls must remain positive and bounded: "+calls);
		long checksum=0;
		final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
		for(long i=0; i<calls; i++){checksum+=(int)call.invokeExact(a, b);}
		final long wall=System.nanoTime()-wallStart, cpu=clock.getCurrentThreadCpuTime()-cpuStart;
		final long expectedChecksum=Math.multiplyExact(calls, (long)expectedReturn);
		final boolean unchanged=unchanged();
		append(out.resolve("measurements.tsv"), stage+"\t"+attempt+"\t"+phase+"\t"+name+"\t"+calls+"\t"+expectedReturn+"\t"+checksum+"\t"+expectedChecksum+"\t"+wall+"\t"+cpu+"\t"+(cpu>=MIN_CPU_NANOS)+"\t"+unchanged+"\n");
		require(cpu>0 && wall>0 && checksum==expectedChecksum && unchanged, "Worker checksum/input/clock failure; inspect retained measurements.tsv.");
		return new Result(calls, checksum, wall, cpu);
	}

	/** Proves loader order from actual classes; the runner pins accepted build inventories. */
	private void verifyProvenance(final String tag) throws Exception{
		final StringBuilder text=new StringBuilder("class\tcodeSource\tsha80\n");
		for(final String name : new String[]{"jgi.BBMerge", "jgi.BBMerge$MateThread", "test.jgi.BBMergeErrorFilterTest", "test.jgi.BBMergeErrorFilterTest$Spec",
			"test.jgi.BBMergeErrorFilterTiming", "test.jgi.BBMergeErrorFilterTiming$Result", "stream.Read", "align2.QualityTools"}){
			final Class<?> type=Class.forName(name, false, BBMergeErrorFilterTiming.class.getClassLoader());
			final Path actual=Paths.get(type.getProtectionDomain().getCodeSource().getLocation().toURI()).toRealPath();
			final Path wanted=name.startsWith("test.jgi.BBMergeErrorFilterTiming") ? timingBuild : name.equals("stream.Read") || name.equals("align2.QualityTools") ? root.resolve("current").toRealPath() : build;
			require(actual.equals(wanted), "Wrong CodeSource for "+name+": "+actual+" instead of "+wanted+"; selected build must come first.");
			text.append(name).append('\t').append(actual).append('\t').append(pin(actual.resolve(name.replace('.', '/')+".class"))).append('\n');
		}
		require(pin(root.resolve("networks/bbmerge.bbnet.gz")).equals(MODEL_PIN), "Runtime model bytes differ from the accepted correctness reference.");
		write(out.resolve(tag+".provenance.tsv"), text.toString());
		if(tag.equals("after")){require(read(out.resolve("before.provenance.tsv")).equals(text.toString()), "Loaded class bytes/locations changed during timing.");}
	}

	/** Input metadata may change inside the real worker; its bases and qualities may not. */
	private boolean unchanged(){return Arrays.equals(ab, a.bases) && Arrays.equals(aq, a.quality) && Arrays.equals(bb, b.bases) && Arrays.equals(bq, b.quality);}
	private static String record(final String name, final Read r){return name+"\t"+r.id+"\t"+new String(r.bases, StandardCharsets.US_ASCII)+"\t"+Arrays.toString(r.quality)+"\n";}
	private static Object invoke(final Class<?> type, final Object owner, final String name, final Class<?>[] types, final Object... args) throws Exception{final Method m=type.getDeclaredMethod(name, types); m.setAccessible(true); return m.invoke(owner, args);}
	private static Object value(final Object owner, final String name) throws Exception{final Field f=owner.getClass().getDeclaredField(name); f.setAccessible(true); return f.get(owner);}
	private static boolean bool(final Object owner, final String name) throws Exception{return (Boolean)value(owner, name);}
	private static float number(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).floatValue();}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required path option."); return Paths.get(value).toAbsolutePath().normalize();}
	private static Path existing(final String value) throws IOException{return path(value).toRealPath();}
	private static void write(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);}
	private static String pin(final Path path) throws Exception{
		final byte[] bytes=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(path)); final StringBuilder text=new StringBuilder();
		for(int i=bytes.length-10; i<bytes.length; i++){text.append(Character.forDigit((bytes[i]>>>4)&15, 16)).append(Character.forDigit(bytes[i]&15, 16));}
		return text.toString();
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static final class Result{
		Result(final long calls_, final long checksum_, final long wall_, final long cpu_){calls=calls_; checksum=checksum_; wall=wall_; cpu=cpu_;}
		final long calls, checksum, wall, cpu;
	}

	private final String phase, name;
	private final int expectedReturn;
	private final Path root, build, timingBuild, out;
	private final ThreadMXBean clock;
	private final MethodHandle call;
	private final Read a, b;
	private final byte[] ab, aq, bb, bq;
	private static final String RECIPE="bbmerge-worker-errorfilter-v1";
	private static final String MODEL_PIN="24998a8aa257552d4a22";
	private static final long INITIAL_CALLS=4096, WARMUP_CALLS=4096, MAX_CALLS=1L<<29;
	private static final int MAX_ATTEMPTS=6;
	private static final long MIN_CPU_NANOS=10000000000L, TARGET_CPU_NANOS=12000000000L;
}
