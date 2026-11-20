package test.jgi;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import jgi.AddAdapters;
import parse.Parse;

/**
 * Finite AddAdapters generation cleanup baseline/fixed check for owner-thread
 * assertion after input and output streams start. The absent SAM case reaches
 * a null-base read through a public SAM representation; the control proves the
 * same child harness can generate one ordinary one-base FASTQ record.
 * @author Jean
 */
public final class AddAdaptersGenerationCleanupTest{

	/** Requires mode=parent plus phase/tree/outdir, or mode=child plus kind/input/output/work. */
	public static void main(final String[] args) throws Throwable{
		String mode=null, phase=null, sourceTree=null, directory=null, kind=null, input=null, output=null, work=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("tree")){sourceTree=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("kind")){kind=pair[1];}
			else if(pair[0].equals("input")){input=pair[1];}
			else if(pair[0].equals("output")){output=pair[1];}
			else if(pair[0].equals("work")){work=pair[1];}
			else{throw new IllegalArgumentException("Unknown AddAdapters cleanup fixture argument: "+arg);}
		}
		if("child".equals(mode)){
			child(kind, input, output, work);
			return;
		}
		require("parent".equals(mode), "Set mode=parent|child.");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(sourceTree!=null && directory!=null, "Set tree and a new outdir.");
		final Path tree=Paths.get(sourceTree).toAbsolutePath().normalize();
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		require(safePath(tree) && safePath(root), "Paths must survive argument parsing and fit one evidence row.");
		require(Files.isDirectory(tree.resolve("current")), "Missing BBTools current classpath: "+tree);
		require(Files.isExecutable(Paths.get("/usr/bin/timeout")), "Owned child supervision requires GNU timeout.");
		Files.createDirectory(root);
		final AddAdaptersGenerationCleanupTest test=new AddAdaptersGenerationCleanupTest(tree, root, phase.equals("baseline"));
		test.scenario("control-fastq", false);
		test.scenario("absent-sam", true);
		test.check(test.cases==2 && test.children==2, "Frozen generation-cleanup matrix requires2 cases/2 children: "+test.cases+"/"+test.children);
		if(test.failures>0){throw new AssertionError("AddAdapters generation cleanup failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADD_ADAPTERS_GENERATION_CLEANUP_PASS phase="+phase+" cases="+test.cases+" children="+test.children+" checks="+test.checks);
	}

	private AddAdaptersGenerationCleanupTest(final Path tree_, final Path root_, final boolean baseline_){
		tree=tree_; root=root_; baseline=baseline_;
	}

	/** Runs one supervised child that invokes AddAdapters.main in-process. */
	private void scenario(final String name, final boolean absentSam) throws Exception{
		final Path work=Files.createDirectory(root.resolve(name));
		final Path input=work.resolve(absentSam ? "absent.sam" : "control.fq");
		final Path output=work.resolve("generated.fq");
		write(input, absentSam ? ABSENT_SAM : CONTROL_FASTQ);
		write(work.resolve("input-before."+(absentSam ? "sam" : "fq")), absentSam ? ABSENT_SAM : CONTROL_FASTQ);
		final int exit=run(work, absentSam ? "absent" : "control", input, output);
		final String stdout=readIfExists(work.resolve("stdout"));
		final String stderr=readIfExists(work.resolve("stderr"));
		final String primary=readIfExists(work.resolve("primary-stack.txt"));
		final String snapshot=readIfExists(work.resolve("thread-snapshot.tsv"));
		final boolean gnuTimeout=Files.exists(work.resolve("gnu-timeout.txt"));
		final boolean outerTimeout=Files.exists(work.resolve("outer-timeout.txt"));
		final long outputSize=size(output);
		final boolean outputExists=Files.exists(output);
		final boolean outputRegular=Files.isRegularFile(output);
		final boolean noFastqData=!outputExists || (outputRegular && outputSize==0);
		final boolean inputUnchanged=saveBytes(work, input, "input-after."+(absentSam ? "sam" : "fq"), absentSam ? ABSENT_SAM : CONTROL_FASTQ);
		final WorkerState worker=ownedOutputWorker(snapshot);
		write(work.resolve("observed.tsv"), "exit\tgnuTimeout\touterTimeout\tprimaryStack\townedOutputWorker\townedOutputNonDaemon\toutputExists\toutputRegular\toutputSize\tinputUnchanged\n"+
			exit+"\t"+gnuTimeout+"\t"+outerTimeout+"\t"+(!primary.isEmpty())+"\t"+worker.present+"\t"+
			worker.nonDaemon+"\t"+outputExists+"\t"+outputRegular+"\t"+outputSize+"\t"+inputUnchanged+"\n");
		check(stdout.length()==0, name+": child unexpectedly wrote stdout.");
		check(inputUnchanged, name+": child must not modify its input bytes.");
		if(absentSam){
			check(noFastqData, name+": absent-base assertion must not write FASTQ data; exists="+outputExists+" regular="+outputRegular+" size="+outputSize);
			check(!hasSuccessSummary(stderr), name+": absent-base assertion printed a success summary.");
			check(primary.contains("java.lang.AssertionError") && primary.contains("jgi.AddAdapters.addAdapter"),
				name+": primary stack must show the AddAdapters.addAdapter assertion.");
			check(stderr.contains("java.lang.AssertionError") && stderr.contains("jgi.AddAdapters.addAdapter"),
				name+": child stderr must show the rethrown main-thread assertion.");
			check(!stderr.contains("SamStreamer worker failed") && !stderr.contains("KillSwitch"),
				name+": public trigger must not be hidden behind a worker KillSwitch path.");
			if(baseline){
				check(gnuTimeout && !outerTimeout, name+": baseline must leave an owned worker alive for GNU timeout without parent cleanup.");
				check(worker.present, name+": baseline snapshot must retain a live ReadStreamByteWriter.");
				check(worker.nonDaemon, name+": baseline retained ReadStreamByteWriter must be non-daemon.");
			}else{
				check(!gnuTimeout && !outerTimeout, name+": fixed child must terminate naturally.");
				check(exit==1, name+": fixed child must exit1 after rethrowing the primary assertion, not a signal exit: "+exit);
				check(!worker.present, name+": fixed snapshot must not retain a live owned output worker.");
			}
		}else{
			check(!gnuTimeout && !outerTimeout, name+": control must terminate naturally.");
			check(exit==0, name+": control wrong Java exit: "+exit);
			check(matches(output, CONTROL_OUTPUT), name+": control generated unexpected FASTQ.");
			check(!stderr.contains("AssertionError") && !stderr.contains("terminated in an error state"),
				name+": control must have clean diagnostics.");
		}
		cases++;
	}

	/**
	 * Child mode owns the AddAdapters driver on this JVM main thread, records the exact
	 * escaping failure and all still-live thread stacks after return from AddAdapters,
	 * then rethrows so the parent can distinguish natural termination from timeout.
	 */
	private static void child(final String kind, final String input, final String output, final String work) throws Throwable{
		require(("absent".equals(kind) || "control".equals(kind)) && input!=null && output!=null && work!=null,
			"Child requires kind=absent|control plus input/output/work.");
		final Path root=Paths.get(work).toAbsolutePath().normalize();
		final String[] addAdaptersArgs={
			"mode=write", "in="+input, "out="+output, "literal=C", "rate=1", "adderrors=f", "ow=t", "app=f",
			"t=1", "simd=f", "int=f", "qin=33", "qout=33", "arc=f", "right=t", "reads=-1"
		};
		try{
			AddAdapters.main(addAdaptersArgs);
			write(root.resolve("child-returned.txt"), "AddAdapters.main returned normally for "+kind+"\n");
		}catch(Throwable t){
			write(root.resolve("primary-stack.txt"), stack(t));
			write(root.resolve("thread-snapshot.tsv"), snapshot());
			throw t;
		}
	}

	/** Launches an owned child through GNU timeout; this parent never inspects unrelated JVMs. */
	private int run(final Path work, final String kind, final Path input, final Path output) throws Exception{
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/timeout", "--verbose",
			"--signal=TERM", "--kill-after=2s", "8s", "bash",
			tree.resolve("dev/addadaptersgenerationcleanuptest.sh").toString(),
			"mode=child", "kind="+kind, "input="+input, "output="+output, "work="+work));
		Files.write(work.resolve("command.argv"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(tree.toFile())
			.redirectOutput(work.resolve("stdout").toFile()).redirectError(work.resolve("stderr").toFile());
		final String inherited=builder.environment().get("JAVA_TOOL_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", (inherited==null ? "" : inherited+" ")+"-Dfile.encoding=UTF-8");
		final Process process=builder.start();
		Throwable primary=null;
		boolean exceeded=false;
		try{
			process.getOutputStream().close();
			exceeded=!process.waitFor(15, TimeUnit.SECONDS);
		}catch(Exception e){primary=e; throw e;}
		catch(Error e){primary=e; throw e;}
		finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){if(primary==null){throw e;} primary.addSuppressed(e);}
			catch(Error e){if(primary==null){throw e;} primary.addSuppressed(e);}
		}
		final int exit=process.exitValue();
		write(work.resolve("exit.txt"), exit+"\n");
		final String diagnostics=readIfExists(work.resolve("stderr"));
		if(exceeded){write(work.resolve("outer-timeout.txt"), "Parent reaped the owned child after the outer 15s fallback.\n");}
		if(exit==124 && diagnostics.contains("timeout: sending signal TERM")){
			write(work.resolve("gnu-timeout.txt"), "GNU timeout sent TERM after the expected 8s deadline and exited 124.\n");
		}
		children++;
		return exit;
	}

	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 5)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){throw new AssertionError("Owned AddAdapters cleanup child survived forced cleanup.");}
	}

	/** Interrupts cannot cancel cleanup; restores the flag after the owned process is reaped. */
	private static boolean waitOwnedChild(final Process process, final int seconds){
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(seconds);
		boolean interrupted=false;
		try{
			while(process.isAlive()){
				final long remaining=deadline-System.nanoTime();
				if(remaining<=0){return false;}
				try{return process.waitFor(remaining, TimeUnit.NANOSECONDS);}
				catch(InterruptedException e){interrupted=true;}
			}
			return true;
		}finally{if(interrupted){Thread.currentThread().interrupt();}}
	}

	private static String snapshot(){
		final List<Map.Entry<Thread, StackTraceElement[]>> entries=
			new ArrayList<Map.Entry<Thread, StackTraceElement[]>>(Thread.getAllStackTraces().entrySet());
		Collections.sort(entries, new Comparator<Map.Entry<Thread, StackTraceElement[]>>(){
			@Override
			public int compare(final Map.Entry<Thread, StackTraceElement[]> a, final Map.Entry<Thread, StackTraceElement[]> b){
				final int name=a.getKey().getName().compareTo(b.getKey().getName());
				return name==0 ? Long.compare(a.getKey().getId(), b.getKey().getId()) : name;
			}
		});
		final StringBuilder sb=new StringBuilder(4096);
		for(final Map.Entry<Thread, StackTraceElement[]> entry : entries){
			final Thread thread=entry.getKey();
			sb.append("thread\t").append(thread.getId()).append('\t').append(thread.getName()).append('\t')
				.append(thread.getClass().getName()).append('\t').append(thread.getState()).append('\t')
				.append(thread.isDaemon()).append('\n');
			for(final StackTraceElement frame : entry.getValue()){
				sb.append("frame\t").append(frame.getClassName()).append('\t').append(frame.getMethodName()).append('\t')
					.append(frame.getFileName()).append('\t').append(frame.getLineNumber()).append('\n');
			}
		}
		return sb.toString();
	}

	private static WorkerState ownedOutputWorker(final String snapshot){
		boolean present=false, nonDaemon=false;
		for(final String line : snapshot.split("\n")){
			if(!line.startsWith("thread\t")){continue;}
			final String[] split=line.split("\t", -1);
			if(split.length>=6 && split[3].equals("stream.ReadStreamByteWriter") && !split[4].equals("TERMINATED")){
				present=true;
				nonDaemon|=split[5].equals("false");
			}
		}
		return new WorkerState(present, nonDaemon);
	}

	private static String stack(final Throwable t){
		final StringWriter out=new StringWriter();
		t.printStackTrace(new PrintWriter(out));
		return out.toString();
	}

	private static boolean matches(final Path path, final String expected) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static boolean saveBytes(final Path work, final Path path, final String name, final String expected) throws Exception{
		if(!Files.isRegularFile(path)){write(work.resolve(name+".missing"), "Input is no longer a regular file.\n"); return false;}
		final byte[] actual=Files.readAllBytes(path);
		Files.write(work.resolve(name), actual, StandardOpenOption.CREATE_NEW);
		return Arrays.equals(actual, expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static long size(final Path path) throws Exception{
		return Files.isRegularFile(path) ? Files.size(path) : 0;
	}
	private static boolean hasSuccessSummary(final String text){
		return text.contains("Adapters Added:") || text.contains("Valid Output:");
	}
	private static String readIfExists(final Path path) throws Exception{
		return Files.isRegularFile(path) ? new String(Files.readAllBytes(path), StandardCharsets.US_ASCII) : "";
	}
	private static boolean safePath(final Path path){
		final String text=path.toString();
		return text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0 && text.indexOf('\r')<0;
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}
	}

	private static final String ABSENT_SAM="@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:1\nabsent0\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
	private static final String CONTROL_FASTQ="@present0\nG\n+\nI\n";
	private static final String CONTROL_OUTPUT="@0_1_0\nC\n+\nI\n";

	private static final class WorkerState{
		WorkerState(final boolean present_, final boolean nonDaemon_){present=present_; nonDaemon=nonDaemon_;}
		final boolean present, nonDaemon;
	}

	private final Path tree, root;
	private final boolean baseline;
	private int checks, failures, cases, children;
}
