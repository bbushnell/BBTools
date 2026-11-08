package test.bloom;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
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

import bloom.ErrorCorrect;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import fileIO.ReadWrite;
import shared.Shared;
import stream.FASTQ;

/** Linux-only lifecycle checks for the legacy file-owning ErrorCorrect.detect overload.
 * Each case runs in an owned, bounded child JVM so a leaked reader cannot hang the runner.
 * Retains commands, stdout, stderr and post-return descriptor/thread observations.
 * No production algorithm or stream factory is replaced.
 * @author Brian Bushnell, Sucrose */
public final class ErrorCorrectLifecycleTest{

	/** Runs three isolated cases; child mode is used only by this runner. */
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
		check(Files.isDirectory(Paths.get("/proc/self/fd")) && Files.exists(Paths.get("/dev/full")),
				"This finite test requires Linux /proc/self/fd and /dev/full.");
		Files.createDirectories(dir);
		int failures=0;
		for(String name : new String[]{"success", "full", "construction"}){
			failures+=runChild(dir, name);
		}
		if(failures>0){throw new AssertionError("Legacy stream lifecycle: "+failures+" failures of 3 cases.");}
		System.out.println("ERROR_CORRECT_LIFECYCLE_TEST_PASS cases=3");
	}

	/** Runs only this fixture's child; a timeout is failed evidence, never accepted cleanup. */
	private static int runChild(final Path dir, final String name) throws Exception{
		final Path caseDir=dir.resolve(name);
		Files.createDirectory(caseDir);
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(ErrorCorrectLifecycleTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx256m");
		if(!System.getProperty("java.specification.version").startsWith("1.")){
			command.add("--add-modules=jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(ErrorCorrectLifecycleTest.class.getName());
		command.add("child="+name);
		command.add("outdir="+caseDir);
		Files.write(caseDir.resolve("command.txt"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).redirectInput(Paths.get("/dev/null").toFile()).redirectOutput(caseDir.resolve("stdout.txt").toFile())
				.redirectError(caseDir.resolve("stderr.txt").toFile()).start();
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
		final String stderr=new String(Files.readAllBytes(caseDir.resolve("stderr.txt")), StandardCharsets.UTF_8);
		//A writer may halt the JVM under -ea; natural nonzero exit plus the real device
		//error is fail-loud evidence. OS process exit closes descriptors; no in-process
		//cleanup claim is made for that route. An owned timeout never qualifies.
		final boolean haltedFull=name.equals("full") && completed && exit!=0 && !Files.exists(resultFile) &&
				stderr.contains("No space left on device");
		final boolean pass=completed && ((exit==0 && "true".equals(result.getProperty("pass"))) || haltedFull);
		final String summary="case="+name+" completed="+completed+" exit="+exit+" haltedFull="+haltedFull+
				" postReturn="+result+" pass="+pass;
		Files.write(caseDir.resolve("outcome.txt"), Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return pass ? 0 : 1;
	}

	/** Calls the real file-owning API, then observes its resources without closing them for it. */
	private static void child(final String name, final Path dir) throws Exception{
		check(name.equals("success") || name.equals("full") || name.equals("construction"), "Unknown child case.");
		Shared.setThreads(1);
		Shared.setBuffers(2);
		Shared.setBufferLen(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
		FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT=33;
		ErrorCorrect.OUTPUT_INFO=ErrorCorrect.verbose=false;
		final Path input=dir.resolve("input.fq");
		final StringBuilder text=new StringBuilder();
		final int records=(name.equals("construction") ? 128 : 1);
		for(int i=0; i<records; i++){text.append("@read").append(i).append("\nCCCCCC\n+\n??????\n");}
		final byte[] expected=text.toString().getBytes(StandardCharsets.US_ASCII);
		Files.write(input, expected);
		final Path output;
		if(name.equals("full")){output=Paths.get("/dev/full");}
		else if(name.equals("construction")){
			final Path blocker=dir.resolve("regular-file");
			Files.write(blocker, new byte[]{1});
			output=blocker.resolve("output.fq");
		}else{output=dir.resolve("output.fq");}
		//Validate the descriptor observer with a deliberately open file in this child.
		final Path probe=dir.resolve("descriptor-probe");
		try(OutputStream held=Files.newOutputStream(probe)){
			held.flush();
			check(openDescriptors(probe)>0, "Descriptor observer missed its positive control.");
		}
		check(openDescriptors(probe)==0, "Descriptor observer failed its closed control.");
		KCountArray7MTA.setSeed(20261105);
		final KCountArray table=KCountArray.makeNew(2000003, 4, 3);
		table.increment(341, 5);//CCCCC
		table.increment(682, 5);//GGGGG: support both legacy query conventions.
		check(table.read(341)==5 && table.read(682)==5, "Trusted literal windows must be present.");
		Throwable failure=null;
		try{ErrorCorrect.detect(input.toString(), null, table, 5, 2, -1, output.toString(), true, false, true);}
		catch(Throwable caught){failure=caught; caught.printStackTrace(System.err);}
		int live=liveStreamThreads(false);
		for(int i=0; i<20 && live>0; i++){
			Thread.sleep(25);
			live=liveStreamThreads(false);
		}
		liveStreamThreads(true);
		final int inputOpen=openDescriptors(input);
		final int outputOpen=(Files.exists(output) ? openDescriptors(output) : 0);
		final boolean errorExpected=!name.equals("success");
		final boolean outputMatches=errorExpected || (Files.exists(output) && Arrays.equals(expected, Files.readAllBytes(output)));
		final boolean pass=(errorExpected ? failure!=null : failure==null) && inputOpen==0 && outputOpen==0 && live==0 && outputMatches;
		final Properties result=new Properties();
		result.setProperty("exception", failure==null ? "none" : failure.getClass().getName());
		result.setProperty("inputOpen", ""+inputOpen);
		result.setProperty("outputOpen", ""+outputOpen);
		result.setProperty("liveStreamThreads", ""+live);
		result.setProperty("outputMatches", ""+outputMatches);
		result.setProperty("pass", ""+pass);
		try(OutputStream out=Files.newOutputStream(dir.resolve("result.properties"))){result.store(out, "Post-return observations, before JVM shutdown");}
		//Return naturally. Leaked non-daemon readers remain visible to the parent's deadline.
	}

	/** Counts stream workers using their concrete class or active Java stack, excluding this thread. */
	private static int liveStreamThreads(final boolean report){
		int count=0;
		for(Map.Entry<Thread, StackTraceElement[]> entry : Thread.getAllStackTraces().entrySet()){
			final Thread thread=entry.getKey();
			if(thread==Thread.currentThread() || !thread.isAlive()){continue;}
			boolean stream=thread.getClass().getName().startsWith("stream.");
			for(StackTraceElement frame : entry.getValue()){stream|=frame.getClassName().startsWith("stream.");}
			if(stream){
				count++;
				if(report){
					System.err.println("LIVE_STREAM_THREAD name="+thread.getName()+" class="+thread.getClass().getName()+" state="+thread.getState());
					for(StackTraceElement frame : entry.getValue()){System.err.println("  at "+frame);}
				}
			}
		}
		return count;
	}

	/** Counts only descriptors targeting the exact fixture path. */
	private static int openDescriptors(final Path path) throws IOException{
		final Path wanted=path.toRealPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){
			for(Path descriptor : descriptors){
				try{if(Files.readSymbolicLink(descriptor).equals(wanted)){count++;}}
				catch(NoSuchFileException vanished){/* Transient descriptor closed during enumeration. */}
			}
		}
		return count;
	}

	/** Keeps fixture checks active in either assertion mode. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
