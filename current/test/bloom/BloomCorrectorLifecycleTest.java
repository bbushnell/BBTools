package test.bloom;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import bloom.BloomFilterCorrectorWrapper;
import shared.Shared;
import shared.Timer;
import stream.Read;

/** Linux-only owned-stream cleanup checks for the correction wrapper (BLM214).
 * Each case has a fresh child JVM and a finite parent deadline. Failure cases replace
 * an empty output directory after wrapper construction, before process opens streams.
 * This exercises real output construction failures independently of path preflight.
 * No stream factory is replaced; the observer never cleans up product resources.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorLifecycleTest{

	/** Runs three isolated cases; child mode is private to this fixture's launcher. */
	public static void main(final String[] args) throws Exception{
		if(args.length==2 && args[0].startsWith("child=") && args[1].startsWith("outdir=")){
			child(args[0].substring(6), Paths.get(args[1].substring(7)));
			return;
		}
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path dir=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(dir), "Result directory must be new.");
		check(Files.isDirectory(Paths.get("/proc/self/fd")), "This fixture requires Linux /proc/self/fd.");
		Files.createDirectories(dir);
		int failures=0;
		for(String name : new String[]{"success", "primary-construction", "rejected-construction"}){
			failures+=runChild(dir, name);
		}
		if(failures>0){throw new AssertionError("Correction wrapper lifecycle: "+failures+" failures of 3 cases.");}
		System.out.println("BLOOM_CORRECTOR_LIFECYCLE_TEST_PASS cases=3");
	}

	/** Captures a child and terminates only that owned process if its deadline expires.
	 * Natural exit, returned observations and all checks must pass; halting is not cleanup. */
	private static int runChild(final Path dir, final String name) throws Exception{
		final Path caseDir=Files.createDirectory(dir.resolve(name));
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(BloomCorrectorLifecycleTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx256m");
		if(!System.getProperty("java.specification.version").startsWith("1.")){
			command.add("--add-modules=jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(BloomCorrectorLifecycleTest.class.getName());
		command.add("child="+name);
		command.add("outdir="+caseDir);
		Files.write(caseDir.resolve("command.txt"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).redirectInput(Paths.get("/dev/null").toFile())
				.redirectOutput(caseDir.resolve("stdout.txt").toFile()).redirectError(caseDir.resolve("stderr.txt").toFile()).start();
		boolean completed=false;
		try{
			completed=process.waitFor(15, TimeUnit.SECONDS);
		}finally{
			if(process.isAlive()){
				process.destroy();
				if(!process.waitFor(2, TimeUnit.SECONDS)){
					process.destroyForcibly();
					check(process.waitFor(5, TimeUnit.SECONDS), "Owned child did not terminate after forced destruction.");
				}
			}
		}
		final int exit=process.exitValue();
		final Properties result=new Properties();
		final Path resultFile=caseDir.resolve("result.properties");
		if(Files.exists(resultFile)){
			try(InputStream input=Files.newInputStream(resultFile)){result.load(input);}
		}
		final boolean pass=completed && exit==0 && "true".equals(result.getProperty("pass"));
		final String summary="case="+name+" completed="+completed+" exit="+exit+" postReturn="+result+" pass="+pass;
		Files.write(caseDir.resolve("outcome.txt"), Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return pass ? 0 : 1;
	}

	/** Builds the real table, changes only a fixture-owned directory, invokes process,
	 * and records observations before natural JVM shutdown. Reflection accesses the
	 * package-private lifecycle boundary without exposing a new production API. */
	private static void child(final String name, final Path dir) throws Exception{
		final boolean success=name.equals("success");
		check(success || name.equals("primary-construction") || name.equals("rejected-construction"), "Unknown child case.");
		final Path input=dir.resolve("input.fq"), reference=dir.resolve("reference.fa");
		final Path primaryParent=Files.createDirectory(dir.resolve("primary"));
		final Path rejectedParent=Files.createDirectory(dir.resolve("rejected"));
		final Path primary=primaryParent.resolve("output.fq"), rejected=rejectedParent.resolve("output.fq");
		final Path brokenParent=(name.equals("primary-construction") ? primaryParent : rejectedParent);
		final Path brokenOutput=(name.equals("primary-construction") ? primary : rejected);
		final StringBuilder text=new StringBuilder();
		//One-read batches and two buffers force unfinished consumption to remain visible.
		for(int i=0; i<128; i++){text.append("@read").append(i).append("\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");}
		final byte[] expected=text.toString().getBytes(StandardCharsets.US_ASCII);
		Files.write(input, expected);
		Files.write(reference, REFERENCE);
		final Path probe=dir.resolve("descriptor-probe");
		try(OutputStream held=Files.newOutputStream(probe)){
			held.flush();
			check(openDescriptors(probe)>0, "Descriptor observer missed its positive control.");
		}
		check(openDescriptors(probe)==0, "Descriptor observer failed its closed control.");
		final String[] flags={"in="+input, "ref="+reference, "out="+primary, "outbad="+rejected,
				"ow=t", "k=15", "ksmall=15", "bits=4", "hashes=1", "cells=10007", "seed=1", "minprob=0",
				"t=1", "readbufferlen=1", "readbuffers=2", "interleaved=f", "pigz=f", "unpigz=f", "qin=33", "qout=33",
				"ecc=f", "merge=f", "ecco=f", "markerrors=f", "mincount=0", "tossjunk=f", "loglog=f"};
		Files.write(dir.resolve("wrapper-args.txt"), Arrays.asList(flags), StandardCharsets.UTF_8);
		final BloomFilterCorrectorWrapper wrapper=new BloomFilterCorrectorWrapper(flags);
		check(Shared.threads()==1, "One worker is required to change validation from false to true during process.");
		check(liveStreamThreads(false)==0 && openDescriptors(input)==0 && openDescriptors(reference)==0,
				"Table construction must finish with no stream workers or fixture input descriptors.");
		check(!Files.exists(primary) && !Files.exists(rejected), "Wrapper construction must not open output writers.");
		if(!success){
			//This is an empty directory created above, never an existing user artifact.
			Files.delete(brokenParent);
			Files.write(brokenParent, new byte[]{1});
		}
		Read.VALIDATE_IN_CONSTRUCTOR=false;
		final Method process=BloomFilterCorrectorWrapper.class.getDeclaredMethod("process", Timer.class);
		process.setAccessible(true);
		Throwable failure=null;
		try{process.invoke(wrapper, new Timer());}
		catch(InvocationTargetException invocation){failure=invocation.getCause(); failure.printStackTrace(System.err);}
		int live=liveStreamThreads(false);
		for(int i=0; i<20 && live>0; i++){
			Thread.sleep(25);
			live=liveStreamThreads(false);
		}
		liveStreamThreads(true);
		final int inputOpen=openDescriptors(input), referenceOpen=openDescriptors(reference);
		final int primaryOpen=(Files.exists(primary) ? openDescriptors(primary) : 0);
		final int rejectedOpen=(Files.exists(rejected) ? openDescriptors(rejected) : 0);
		final boolean validationRestored=!Read.VALIDATE_IN_CONSTRUCTOR;
		final boolean originalFailure=success ? failure==null : retainsPathFailure(failure, brokenOutput);
		final boolean originalBytesIntact=Arrays.equals(expected, Files.readAllBytes(input)) &&
				Arrays.equals(REFERENCE, Files.readAllBytes(reference));
		final boolean outputMatches=!success || (Files.exists(primary) && Arrays.equals(expected, Files.readAllBytes(primary)) &&
				Files.exists(rejected) && Files.size(rejected)==0);
		final boolean pass=originalFailure && originalBytesIntact && outputMatches && validationRestored &&
				inputOpen==0 && referenceOpen==0 && primaryOpen==0 && rejectedOpen==0 && live==0;
		final Properties result=new Properties();
		result.setProperty("exception", failure==null ? "none" : failure.getClass().getName());
		result.setProperty("originalFailure", ""+originalFailure);
		result.setProperty("validationRestored", ""+validationRestored);
		result.setProperty("inputOpen", ""+inputOpen);
		result.setProperty("referenceOpen", ""+referenceOpen);
		result.setProperty("primaryOpen", ""+primaryOpen);
		result.setProperty("rejectedOpen", ""+rejectedOpen);
		result.setProperty("liveStreamThreads", ""+live);
		result.setProperty("originalBytesIntact", ""+originalBytesIntact);
		result.setProperty("outputMatches", ""+outputMatches);
		result.setProperty("pass", ""+pass);
		try(OutputStream out=Files.newOutputStream(dir.resolve("result.properties"))){result.store(out, "Post-return observations, before JVM shutdown");}
		//No flag reset, stream close or System.exit here: product cleanup must do the work.
	}

	/** Requires the real failed output open in the causal chain, not an unrelated exception. */
	private static boolean retainsPathFailure(final Throwable failure, final Path output){
		for(Throwable cause=failure; cause!=null; cause=cause.getCause()){
			if(cause instanceof FileNotFoundException && cause.getMessage()!=null &&
					cause.getMessage().contains(output.toString()) && cause.getMessage().contains("Not a directory")){return true;}
		}
		return false;
	}

	/** Counts stream workers by concrete class/active stack, including wrapper workers and BF2. */
	private static int liveStreamThreads(final boolean report){
		int count=0;
		for(Map.Entry<Thread, StackTraceElement[]> entry : Thread.getAllStackTraces().entrySet()){
			final Thread thread=entry.getKey();
			if(thread==Thread.currentThread() || !thread.isAlive()){continue;}
			boolean owned=isStreamClass(thread.getClass().getName());
			for(StackTraceElement frame : entry.getValue()){owned|=isStreamClass(frame.getClassName());}
			if(owned){
				count++;
				if(report){
					System.err.println("LIVE_STREAM_THREAD name="+thread.getName()+" class="+thread.getClass().getName()+" state="+thread.getState());
					for(StackTraceElement frame : entry.getValue()){System.err.println("  at "+frame);}
				}
			}
		}
		return count;
	}

	/** Class families used by the real wrapper's input/output/processing workers. */
	private static boolean isStreamClass(final String name){
		return name.startsWith("stream.") || name.startsWith("fileIO.ByteFile") || name.startsWith("bloom.BloomFilterCorrectorWrapper$");
	}

	/** Counts open descriptors targeting the exact existing fixture path. */
	private static int openDescriptors(final Path path) throws IOException{
		final Path wanted=path.toRealPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){
			for(Path descriptor : descriptors){
				try{if(Files.readSymbolicLink(descriptor).equals(wanted)){count++;}}
				catch(NoSuchFileException vanished){/* Descriptor closed during enumeration. */}
			}
		}
		return count;
	}

	/** Keeps the independent fixture checks active in either assertion mode. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final byte[] REFERENCE=">reference\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n".getBytes(StandardCharsets.US_ASCII);
}
