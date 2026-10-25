package test.synth;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

/** Captures the existing-output ow=f startup failure without leaving a hung reader alive.
 * @author Jean
 */
public final class FuseByTaxaOutputPreflightTest{

	/** Runs the captured baseline or the corrected output-preflight contract. */
	public static void main(String[] args) throws Exception{
		check(args.length==2, "Expected <BBTools dir> and expect=baseline or expect=fixed");
		final boolean fixed=args[1].equals("expect=fixed");
		check(fixed || args[1].equals("expect=baseline"), "Unknown expectation: "+args[1]);
		final Path code=Paths.get(args[0]).toAbsolutePath();
		final Path root=Files.createTempDirectory("fusebytaxa-output-preflight-");
		final Result ea=run(code, root, true);
		if(fixed){
			check(!ea.timedOut, "-ea existing-output failure must return to its caller.");
			check(ea.exit!=0, "-ea ow=f existing output must fail before streaming: "+ea.log);
			check(ea.outputPreserved, "-ea ow=f must preserve the existing output.");
		}else{
			check(ea.timedOut, "Baseline -ea should strand the child after ByteStreamWriter asserts: "+ea.log);
			check(ea.logText.contains("AssertionError"), "Hung baseline must carry a visible assertion: "+ea.log);
			check(ea.logText.contains("Please add the flag ow"), "Assertion diagnostic must identify ow=f: "+ea.log);
			check(ea.outputPreserved, "-ea assertion should fire before opening and truncating output.");
		}

		final Result da=run(code, root, false);
		if(fixed){
			check(!da.timedOut, "-da existing-output failure must return to its caller.");
			check(da.exit!=0, "-da ow=f existing output must fail before streaming: "+da.log);
			check(da.outputPreserved, "-da ow=f must preserve the existing output.");
		}else{
			check(!da.timedOut, "Baseline -da should complete after the constructor assertion is disabled: "+da.log);
			check(da.exit==0, "Baseline -da currently overwrites the existing output: "+da.log);
			check(!da.outputPreserved, "Baseline -da must expose the missing runtime preflight.");
		}
		check(ea.inputPreserved && da.inputPreserved, "Refusing one output must never mutate the input FASTA.");
		System.out.println("FUSEBYTAXA_OUTPUT_PREFLIGHT_TEST_PASS expect="+(fixed ? "fixed" : "baseline"));
	}

	/** Launches only this test's child JVM through the wrapper and kills only that process on timeout. */
	private static Result run(final Path code, final Path root, final boolean assertions) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(assertions ? "ea" : "da"));
		final Path in=dir.resolve("input.fa");
		final Path out=dir.resolve("out.fa");
		final Path log=dir.resolve("child.log");
		final Path pgid=dir.resolve("child.pgid");
		final Path launcher=dir.resolve("launch.sh");
		final byte[] input=(">tid|7|a\nA\n>tid|7|b\nC\n>tid|7|c\nG\n"+
				">tid|7|d\nT\n>tid|7|e\nN\n").getBytes(StandardCharsets.US_ASCII);
		final byte[] sentinel="do not overwrite\n".getBytes(StandardCharsets.US_ASCII);
		Files.write(in, input);
		Files.write(out, sentinel);
		Files.write(launcher, ("#!/bin/bash\nset -u\nsetsid \"$@\" &\nchild=$!\n"+
				"printf '%s\\n' \"$child\" > \"$PGID_FILE\"\nwait \"$child\"\n").getBytes(StandardCharsets.US_ASCII));

		final List<String> command=new ArrayList<String>(Arrays.asList("bash",
				launcher.toString(), "bash", code.resolve("fusebytaxa.sh").toString(), "in="+in, "out="+out, "ow=f",
				"floor=0", "maxlen=1k", "pad=0", "t=4", "readbufferlen=1",
				"fastastreamer2=f", "usepigz=f", "useunpigz=f",
				"-Xmx256m", "-Xms256m", assertions ? "-ea" : "-da"));
		final ProcessBuilder pb=makeChild(command, log);
		pb.environment().put("PGID_FILE", pgid.toString());
		final Process child=pb.start();
		child.getOutputStream().close();
		final boolean finished=child.waitFor(5, TimeUnit.SECONDS);
		int exit=Integer.MIN_VALUE;
		if(finished){
			exit=child.exitValue();
		}else{
			terminateOwnedGroup(child, pgid, log);
		}
		final byte[] afterInput=Files.readAllBytes(in);
		final byte[] afterOutput=Files.readAllBytes(out);
		return new Result(!finished, exit, log, new String(Files.readAllBytes(log), StandardCharsets.UTF_8),
				Arrays.equals(input, afterInput), Arrays.equals(sentinel, afterOutput));
	}

	/** Kills only the process group started by this test case's setsid launcher. */
	private static void terminateOwnedGroup(final Process child, final Path pgid, final Path log) throws Exception{
		for(int i=0; i<20 && !Files.exists(pgid); i++){Thread.sleep(50);}
		check(Files.exists(pgid), "Watchdog could not find the owned child process group: "+log);
		final String group=new String(Files.readAllBytes(pgid), StandardCharsets.US_ASCII).trim();
		check(group.matches("[0-9]+"), "Invalid owned process group id "+group+" in "+pgid);
		new ProcessBuilder("kill", "-TERM", "--", "-"+group).start().waitFor();
		if(!waitForGroupExit(group)){
			new ProcessBuilder("kill", "-KILL", "--", "-"+group).start().waitFor();
			child.destroyForcibly();
			check(waitForGroupExit(group), "Owned FuseByTaxa group survived the watchdog: "+log);
		}
		child.waitFor(2, TimeUnit.SECONDS);
		check(!child.isAlive(), "Owned FuseByTaxa child did not die after the watchdog: "+log);
	}

	/** Waits briefly for the process group owned by this fixture to disappear. */
	private static boolean waitForGroupExit(final String group) throws Exception{
		for(int i=0; i<20; i++){
			if(!groupAlive(group)){return true;}
			Thread.sleep(100);
		}
		return !groupAlive(group);
	}

	/** Returns true if any process remains in the owned child process group. */
	private static boolean groupAlive(final String group) throws Exception{
		return new ProcessBuilder("kill", "-0", "--", "-"+group).start().waitFor()==0;
	}

	/** Test expectations stay fatal with assertions disabled in the parent JVM. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Pins wrapper subprocesses to this JVM's Java and a deterministic tiny runtime. */
	private static ProcessBuilder makeChild(final List<String> command, final Path log){
		final ProcessBuilder pb=new ProcessBuilder(command).redirectErrorStream(true)
				.redirectOutput(log.toFile());
		final String javaHome=System.getProperty("java.home");
		final String oldPath=pb.environment().get("PATH");
		pb.environment().put("JAVA_HOME", javaHome);
		pb.environment().put("PATH", Paths.get(javaHome, "bin")+File.pathSeparator+
				(oldPath==null ? "" : oldPath));
		return pb;
	}

	private static final class Result{
		Result(final boolean timedOut_, final int exit_, final Path log_, final String logText_,
				final boolean inputPreserved_, final boolean outputPreserved_){
			timedOut=timedOut_;
			exit=exit_;
			log=log_;
			logText=logText_;
			inputPreserved=inputPreserved_;
			outputPreserved=outputPreserved_;
		}
		final boolean timedOut, inputPreserved, outputPreserved;
		final int exit;
		final Path log;
		final String logText;
	}
}
