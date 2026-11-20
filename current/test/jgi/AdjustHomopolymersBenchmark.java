package test.jgi;

import java.io.File;
import java.io.InputStream;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.regex.Pattern;

import jgi.AdjustHomopolymers;
import parse.Parse;
import shared.Shared;
import stream.Read;

/**
 * Matched direct-helper timing, using retained baseline/fixed production classes.
 * Each fresh JVM runs one variant/profile. Pilot is calibration only and measures
 * at least ten seconds of current-thread CPU; sample uses explicit fixed counts.
 * Literal full output/quality checks and immutable-input checks occur outside timing.
 * The timed region includes input reset, bound invokeExact and cheap result consumption,
 * not constructors, stream processing, full hashing or output validation.
 * @author Sucrose
 */
public final class AdjustHomopolymersBenchmark{

	/**
	 * Requires mode=pilot|sample, variant=baseline|fixed|candidate|reference, profile=singleton|mixed,
	 * targetclasses=retained-class-directory and a new outdir. Samples also require common
	 * iterations and warmup counts; pilot warmup defaults to50000.
	 * Candidate requires candidatepin=20-lowercase-hex and uses fixed-correct output oracles.
	 * Reference requires referencepin=20-lowercase-hex and uses fixed-correct output oracles.
	 */
	public static void main(final String[] args) throws Throwable{
		String mode=null, variant=null, profile=null, targetName=null, directory=null, candidatePin=null, referencePin=null;
		long iterations=-1, warmup=-1;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("variant")){variant=pair[1];}
			else if(pair[0].equals("candidatepin")){candidatePin=pair[1];}
			else if(pair[0].equals("referencepin")){referencePin=pair[1];}
			else if(pair[0].equals("profile")){profile=pair[1];}
			else if(pair[0].equals("targetclasses")){targetName=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("iterations")){iterations=Long.parseLong(pair[1]);}
			else if(pair[0].equals("warmup")){warmup=Long.parseLong(pair[1]);}
			else{throw new IllegalArgumentException("Unknown timing argument: "+arg);}
		}
		require("pilot".equals(mode) || "sample".equals(mode), "Set mode=pilot|sample.");
		require("baseline".equals(variant) || "fixed".equals(variant) || "candidate".equals(variant) || "reference".equals(variant), "Set variant=baseline|fixed|candidate|reference.");
		require("singleton".equals(profile) || "mixed".equals(profile), "Set profile=singleton|mixed.");
		require(targetName!=null && directory!=null, "Set targetclasses and a new outdir.");
		final boolean pilot=mode.equals("pilot"), baseline=variant.equals("baseline"), candidate=variant.equals("candidate"), reference=variant.equals("reference");
		require(candidate ? candidatePin!=null && candidatePin.matches("[0-9a-f]{20}") : candidatePin==null,
			"candidatepin=20 lowercase hex digits is required only for variant=candidate.");
		require(reference ? referencePin!=null && referencePin.matches("[0-9a-f]{20}") : referencePin==null,
			"referencepin=20 lowercase hex digits is required only for variant=reference.");
		if(pilot){require(iterations==-1, "Pilot chooses its count from measured CPU, not iterations=."); if(warmup==-1){warmup=DEFAULT_WARMUP;}}
		else{require(iterations>0 && warmup>0, "Samples require explicit common iterations and warmup.");}
		require(warmup>0, "Warmup must be positive.");
		final Path target=Paths.get(targetName).toRealPath();
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		require(safePath(target) && safePath(root), "Use plain evidence paths without equals, tabs or newlines.");
		require(Files.isDirectory(target), "Target must be the retained compiled package directory: "+target);
		final Path loaded=codeSource(AdjustHomopolymers.class), harness=codeSource(AdjustHomopolymersBenchmark.class);
		require(loaded.equals(target), "Production class was shadowed: requested="+target+" loaded="+loaded);
		final String[] cp=System.getProperty("java.class.path").split(Pattern.quote(File.pathSeparator), -1);
		require(cp.length==2 && Paths.get(cp[0]).toRealPath().equals(target) && Paths.get(cp[1]).toRealPath().equals(harness),
			"Use the dev launcher: classpath must be requested target first, identical shared harness second.");
		final Path productionClass=target.resolve("jgi/AdjustHomopolymers.class");
		final String productionPin=sha80(productionClass), expectedPin=baseline ? BASELINE_CLASS : (candidate ? candidatePin : (reference ? referencePin : FIXED_CLASS));
		require(productionPin.equals(expectedPin), "Wrong retained production bytes: expected="+expectedPin+" actual="+productionPin);
		require(System.getProperty("java.version").equals("21.0.8"), "Matched protocol requires recorded JDK21.0.8.");
		require(Runtime.getRuntime().availableProcessors()==1, "Matched protocol requires one visible CPU.");
		require(AdjustHomopolymersBenchmark.class.desiredAssertionStatus() && AdjustHomopolymers.class.desiredAssertionStatus(), "Enable assertions for harness and target.");
		final List<String> vmArgs=ManagementFactory.getRuntimeMXBean().getInputArguments();
		require(lastOption(vmArgs, "-Xmx").equalsIgnoreCase("-Xmx256m") && lastOption(vmArgs, "-Xms").equalsIgnoreCase("-Xms256m"),
			"Matched protocol requires -Xmx256m and -Xms256m without later overrides.");
		final ThreadMXBean cpu=ManagementFactory.getThreadMXBean();
		require(cpu.isCurrentThreadCpuTimeSupported(), "Current-thread CPU accounting is required; wall time is not a substitute.");
		if(!cpu.isThreadCpuTimeEnabled()){cpu.setThreadCpuTimeEnabled(true);}
		require(cpu.isThreadCpuTimeEnabled() && cpu.getCurrentThreadCpuTime()>=0, "Current-thread CPU accounting is unavailable.");
		Files.createDirectory(root);
		final Data data=new Data(profile, baseline);
		final Path input=root.resolve("constructor-input.fq"), output=root.resolve("unused-output.fq");
		write(input, INPUT); write(root.resolve("constructor-input-before.fq"), INPUT);
		final String[] constructorArgs={"in="+input, "out="+output, "rate=0.5", "int=f", "t=1", "simd=f", "qin=33", "qout=33", "ow=t", "app=f"};
		Files.write(root.resolve("constructor.argv"), Arrays.asList(constructorArgs), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final long constructorStartCpu=cpu.getCurrentThreadCpuTime(), constructorStartWall=System.nanoTime();
		final AdjustHomopolymers app=new AdjustHomopolymers(constructorArgs);
		final long constructorElapsed=System.nanoTime()-constructorStartWall, constructorCpu=cpu.getCurrentThreadCpuTime()-constructorStartCpu;
		require(Shared.threads()==1 && !Shared.SIMD, "Constructor must establish t=1/simd=f.");
		final Method method=AdjustHomopolymers.class.getDeclaredMethod("processRead", Read.class);
		method.setAccessible(true);
		final MethodHandle edit=MethodHandles.lookup().unreflect(method).bindTo(app);
		require(edit.type().equals(MethodType.methodType(void.class, Read.class)), "Expected bound (Read)void editing handle: "+edit.type());
		final AdjustHomopolymersBenchmark bench=new AdjustHomopolymersBenchmark(data, edit);
		bench.validCount(warmup);
		if(!pilot){bench.validCount(iterations);}
		final Path harnessClass=harness.resolve("test/jgi/AdjustHomopolymersBenchmark.class");
		require(codeSource(Data.class).equals(harness) && codeSource(Measurement.class).equals(harness), "Nested harness classes must come from the same shared directory.");
		write(root.resolve("configuration.tsv"), "setting\tvalue\nmode\t"+mode+"\nvariant\t"+variant+"\nprofile\t"+profile+
			"\nrequestedTarget\t"+target+"\nloadedProduction\t"+loaded+"\nproductionSha80\t"+productionPin+
			"\nexpectedProductionSha80\t"+expectedPin+
			"\nharnessSource\t"+harness+"\nharnessSha80\t"+sha80(harnessClass)+
			"\nharnessDataSha80\t"+sha80(harness.resolve("test/jgi/AdjustHomopolymersBenchmark$Data.class"))+
			"\nharnessMeasurementSha80\t"+sha80(harness.resolve("test/jgi/AdjustHomopolymersBenchmark$Measurement.class"))+
			"\nclasspath\t"+System.getProperty("java.class.path")+
			"\njavaVersion\t"+System.getProperty("java.version")+"\njavaVendor\t"+System.getProperty("java.vendor")+
			"\nvmArguments\t"+String.join(" ", vmArgs)+"\nvisibleCPUs\t"+Runtime.getRuntime().availableProcessors()+
			"\nmaxHeapBytes\t"+Runtime.getRuntime().maxMemory()+"\nassertions\ttrue\nworkerThreads\t"+Shared.threads()+"\nsimd\t"+Shared.SIMD+
			"\nconstructor_cpu_ns\t"+constructorCpu+"\nconstructor_elapsed_ns\t"+constructorElapsed+"\ninputBasesPerCall\t"+data.input.length+
			"\noutputBasesPerCall\t"+data.expected.length+"\nwarmupCalls\t"+warmup+"\nrequestedCalls\t"+iterations+"\n");
		Files.write(root.resolve("input-bases.bin"), data.input, StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("input-qualities.bin"), data.inputQual, StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("expected-bases.bin"), data.expected, StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("expected-qualities.bin"), data.expectedQual, StandardOpenOption.CREATE_NEW);
		write(root.resolve("oracles.tsv"), "stage\tinputUnchanged\tfullOutputMatches\n");
		bench.calls(1);
		bench.verify(root, "initial");
		final Measurement warmupResult=bench.sample(cpu, warmup);
		BLACKHOLE=warmupResult.checksum;
		require(BLACKHOLE==Math.multiplyExact(bench.unit, warmup), "Warmup checksum differs from literal output oracle.");
		bench.verify(root, "warmup");
		System.err.println("BENCH_BEGIN mode="+mode+" variant="+variant+" profile="+profile+" warmupCalls="+warmup);
		final Measurement result=pilot ? bench.pilot(cpu) : bench.sample(cpu, iterations);
		BLACKHOLE=result.checksum;
		bench.verify(root, "measured");
		require(result.checksum==Math.multiplyExact(bench.unit, result.calls), "Timed checksum differs from the phase-specific literal oracle.");
		require(sha80(productionClass).equals(expectedPin), "Retained target class changed during the JVM run.");
		require(!Files.exists(output) && Arrays.equals(Files.readAllBytes(input), bytes(INPUT)), "Benchmark must preserve preflight input and leave data output absent.");
		Files.write(root.resolve("constructor-input-after.fq"), Files.readAllBytes(input), StandardOpenOption.CREATE_NEW);
		require(result.cpuNanos>0 && result.elapsedNanos>0, "Measured CPU and elapsed clocks must advance.");
		final long proposed=roundUp((long)Math.ceil(result.calls*(double)TARGET_CPU/result.cpuNanos), PILOT_BATCH);
		final long proposedWarmup=roundUp(Math.max(DEFAULT_WARMUP, (proposed+9)/10), PILOT_BATCH);
		final boolean qualified=result.cpuNanos>=MIN_CPU;
		write(root.resolve("results.tsv"), "mode\tvariant\tprofile\titerations\twarmup\tcpu_ns\telapsed_ns\twarmup_cpu_ns\twarmup_elapsed_ns\tinput_bases_per_call\toutput_bases_per_call\tchecksum\tcandidate_iterations\tcandidate_warmup\n"+
			mode+"\t"+variant+"\t"+profile+"\t"+result.calls+"\t"+warmup+"\t"+result.cpuNanos+"\t"+result.elapsedNanos+
			"\t"+warmupResult.cpuNanos+"\t"+warmupResult.elapsedNanos+"\t"+data.input.length+"\t"+data.expected.length+"\t"+result.checksum+
			"\t"+(pilot ? Long.toString(proposed) : "NA")+"\t"+(pilot ? Long.toString(proposedWarmup) : "NA")+"\n");
		require(qualified, "Timed current-thread CPU must reach10seconds; increase the common count before a matched sample set.");
		System.out.println("ADJUST_HOMOPOLYMERS_BENCH_PASS mode="+mode+" variant="+variant+" profile="+profile+
			" calls="+result.calls+" cpuNanos="+result.cpuNanos+" elapsedNanos="+result.elapsedNanos+" checksum="+result.checksum+
			" suggestedIterations="+proposed+" suggestedWarmup="+proposedWarmup);
	}

	private AdjustHomopolymersBenchmark(final Data data_, final MethodHandle edit_){
		data=data_; edit=edit_;
		read=new Read(data.input, data.inputQual, "timing-read", 7, false);
		unit=consume(data.expected, data.expectedQual);
	}

	/** No per-call varargs, reflection arrays, input clone, constructor or validation work. */
	private long calls(final long count) throws Throwable{
		long sink=0;
		for(long i=0; i<count; i++){
			read.bases=data.input; read.quality=data.inputQual;
			edit.invokeExact(read);
			sink+=consume(read.bases, read.quality);
		}
		return sink;
	}
	private Measurement pilot(final ThreadMXBean cpu) throws Throwable{
		long done=0, sink=0, now;
		final long startCpu=cpu.getCurrentThreadCpuTime(), startWall=System.nanoTime();
		do{
			validCount(Math.addExact(done, PILOT_BATCH));
			sink+=calls(PILOT_BATCH); done=Math.addExact(done, PILOT_BATCH);
			now=cpu.getCurrentThreadCpuTime();
			if(now<startCpu){throw new IllegalStateException("Current-thread CPU clock became unavailable.");}
		}while(now-startCpu<MIN_CPU);
		final long elapsed=System.nanoTime()-startWall;
		return new Measurement(done, sink, now-startCpu, elapsed);
	}
	private Measurement sample(final ThreadMXBean cpu, final long count) throws Throwable{
		final long startCpu=cpu.getCurrentThreadCpuTime(), startWall=System.nanoTime();
		final long sink=calls(count);
		final long elapsed=System.nanoTime()-startWall, used=cpu.getCurrentThreadCpuTime()-startCpu;
		return new Measurement(count, sink, used, elapsed);
	}
	private void validCount(final long count){
		require(count>0 && count<=Long.MAX_VALUE/Math.max(unit, Math.max(data.input.length, data.expected.length)),
			"Call count must keep checksum and work counters representable: "+count);
	}
	private void verify(final Path root, final String stage) throws Exception{
		final boolean immutable=Arrays.equals(data.input, data.original) && Arrays.equals(data.inputQual, data.originalQual);
		final boolean correct=Arrays.equals(read.bases, data.expected) && Arrays.equals(read.quality, data.expectedQual);
		Files.write(root.resolve("oracles.tsv"), bytes(stage+"\t"+immutable+"\t"+correct+"\n"), StandardOpenOption.APPEND);
		require(immutable && correct, stage+": full literal output/quality oracle or immutable input failed.");
	}
	private static long consume(final byte[] bases, final byte[] quality){
		return bases.length+(bases[0]&255)+257L*(bases[bases.length-1]&255)+65537L*(quality[0]&255)+1048583L*(quality[quality.length-1]&255);
	}

	/** Independent literal repeated patterns, including distinct old/fixed mixed outputs. */
	private static final class Data{
		Data(final String profile, final boolean baseline){
			if(profile.equals("singleton")){
				input=repeat(bytes("ACGT"), 640); inputQual=repeat(new byte[]{30, 30, 30, 30}, 640);
				expected=repeat(bytes("ACGT"), 640); expectedQual=repeat(new byte[]{30, 30, 30, 30}, 640);
			}else{
				input=repeat(bytes("AACCN"), 512); inputQual=repeat(new byte[]{30, 30, 30, 30, 0}, 512);
				expected=repeat(bytes(baseline ? "AACACNC" : "AAACCCN"), 512);
				expectedQual=repeat(baseline ? new byte[]{30, 30, 30, 30, 30, 0, 30} : new byte[]{30, 30, 30, 30, 30, 30, 0}, 512);
			}
			original=input.clone(); originalQual=inputQual.clone();
			require(input.length==2560 && inputQual.length==2560 && expected.length==expectedQual.length,
				"Literal timing profiles require2560 input bases and one quality per expected base.");
		}
		final byte[] input, inputQual, expected, expectedQual, original, originalQual;
	}
	private static byte[] repeat(final byte[] pattern, final int count){
		final byte[] out=new byte[Math.multiplyExact(pattern.length, count)];
		for(int i=0; i<count; i++){System.arraycopy(pattern, 0, out, i*pattern.length, pattern.length);}
		return out;
	}
	private static final class Measurement{
		Measurement(final long calls_, final long checksum_, final long cpuNanos_, final long elapsedNanos_){
			calls=calls_; checksum=checksum_; cpuNanos=cpuNanos_; elapsedNanos=elapsedNanos_;
		}
		final long calls, checksum, cpuNanos, elapsedNanos;
	}
	private static Path codeSource(final Class<?> type) throws Exception{
		return Paths.get(type.getProtectionDomain().getCodeSource().getLocation().toURI()).toRealPath();
	}
	private static String lastOption(final List<String> args, final String prefix){
		String found="";
		for(final String arg : args){if(arg.toLowerCase(Locale.ROOT).startsWith(prefix.toLowerCase(Locale.ROOT))){found=arg;}}
		return found;
	}
	private static String sha80(final Path path) throws Exception{
		final MessageDigest hash=MessageDigest.getInstance("SHA-256");
		final byte[] buffer=new byte[8192];
		try(final InputStream in=Files.newInputStream(path)){
			for(int n=in.read(buffer); n>=0; n=in.read(buffer)){hash.update(buffer, 0, n);}
		}
		final byte[] digest=hash.digest();
		final StringBuilder out=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){out.append(HEX[(digest[i]&255)>>>4]).append(HEX[digest[i]&15]);}
		return out.toString();
	}
	private static long roundUp(final long value, final long quantum){
		require(value>0 && value<=Long.MAX_VALUE-(quantum-1), "Iteration recommendation overflow.");
		return ((value+quantum-1)/quantum)*quantum;
	}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, bytes(text), StandardOpenOption.CREATE_NEW);
	}
	private static boolean safePath(final Path path){
		final String text=path.toString();
		return text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0 && text.indexOf('\r')<0;
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	private static final String BASELINE_CLASS="54641e08c28d1efd38f7", FIXED_CLASS="400bc0eee927e2855855";
	private static final String INPUT="@timing-constructor\nG\n+\nI\n";
	private static final long MIN_CPU=10000000000L, TARGET_CPU=12000000000L, PILOT_BATCH=256, DEFAULT_WARMUP=50000;
	private static final char[] HEX="0123456789abcdef".toCharArray();
	private static volatile long BLACKHOLE;
	private final Data data;
	private final MethodHandle edit;
	private final Read read;
	private final long unit;
}
