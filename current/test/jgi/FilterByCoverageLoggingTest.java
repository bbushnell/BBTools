package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import parse.Parse;

/** Captures FilterByCoverage missing-coverage logging failures, repair, and log collision coverage.
 * All cases invoke the real shell launcher in an owned child process and preserve its literal
 * command, exit code, stdout, stderr, log file and clean/dirty outputs for root review.
 * @author Jean
 */
public final class FilterByCoverageLoggingTest{

	/** Requires phase=baseline or fixed and a new outdir; does not alter user inputs. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, scope="logging";
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if("phase".equals(pair[0])){phase=pair[1];}
			else if("outdir".equals(pair[0])){directory=pair[1];}
			else if("scope".equals(pair[0])){scope=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		final boolean collision="collision".equals(scope);
		if(!collision && !"logging".equals(scope)){throw new IllegalArgumentException("Set scope=logging or scope=collision.");}
		check("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed.");
		final boolean fixed="fixed".equals(phase);
		check(directory!=null && directory.indexOf('=')<0,
				"Set a new outdir without equals signs; arguments preserve raw absolute paths.");
		final Path launcher=Paths.get("filterbycoverage.sh").toAbsolutePath().normalize();
		check(Files.isRegularFile(launcher), "Missing actual filter launcher: "+launcher);
		final Path root=Files.createDirectory(Paths.get(directory).toAbsolutePath());
		if(collision){
			logCollision(root, launcher, true, fixed);
			logCollision(root, launcher, false, fixed);
			if(fixed){
				fixedPreflightCollision(root, launcher, "log-clean-collision", true, false,
						"specified for multiple output streams");
				fixedPreflightCollision(root, launcher, "qfout-cov-collision", false, true,
						"specified for both input and output streams");
			}
			System.out.println("FILTER_BY_COVERAGE_LOG_COLLISION_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
		}else{
			missingCoverage(root, launcher, false, fixed);
			missingCoverage(root, launcher, true, fixed);
			coveredControl(root, launcher);
			System.out.println("FILTER_BY_COVERAGE_LOGGING_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
		}
	}

	/** Checks the baseline log failure and the repaired zero-valued row for missing after-coverage. */
	private static void missingCoverage(final Path root, final Path launcher,
			final boolean withCov0, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(withCov0 ? "missing-with-cov0" : "missing-no-cov0"));
		final Path input=work.resolve("missing.fa"), coverage=work.resolve("coverage.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa"), log=work.resolve("results.tsv");
		final String fasta=">absent\nACGTACGT\n";
		write(input, fasta);
		write(coverage, COVERAGE_HEADER+"\n"+coveredRow("other")+"\n");
		final Path coverage0=work.resolve("coverage0.tsv");
		final String expectedLog=LOG_HEADER+(withCov0 ? "\tavgFold0\treads0\tnormRatio" : "")+"\n"
				+(fixed ? "missing\tabsent\t1\t8\t0.00\t0\t0.00"
						+(withCov0 ? "\t0.00\t0\t0.00" : "")+"\n" : "");
		write(work.resolve("expected-log.tsv"), expectedLog);
		write(work.resolve("expected-clean.fa"), "");
		write(work.resolve("expected-dirty.fa"), fixed ? fasta : "");

		final ArrayList<String> arguments=baseArguments(input, coverage, clean, dirty, log);
		if(withCov0){
			write(coverage0, COVERAGE_HEADER+"\n"+coveredRow("absent")+"\n");
			arguments.add("cov0="+coverage0);
		}
		final int exit=run(launcher, arguments, work);

		check(exit==(fixed ? 0 : 1), "Missing coverage exit="+exit+"; inspect "+work);
		checkNoTimeout(work);
		check(Files.size(work.resolve("stdout"))==0, "Missing coverage wrote stdout: "+work);
		check("".equals(read(clean)), "Missing coverage clean output differs: "+work);
		check((fixed ? fasta : "").equals(read(dirty)), "Missing coverage dirty output differs: "+work);
		check(expectedLog.equals(read(log)), "Missing coverage log differs: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		if(fixed){
			check(!diagnostics.contains("IllegalFormatConversionException"),
					"Missing coverage retained the format exception: "+work);
			check(diagnostics.contains("\nTime:"), "Missing coverage lacked the final report: "+work);
			check(stat(diagnostics, "Reads In", 1) && stat(diagnostics, "Reads Out", 0)
					&& stat(diagnostics, "Reads Filtered", 1), "Missing coverage fixed counters differ: "+work);
		}else{
			check(diagnostics.contains("IllegalFormatConversionException"), "Missing coverage lacks format exception: "+work);
			check(diagnostics.contains("f != java.lang.Integer"), "Missing coverage lacks integer-to-float diagnostic: "+work);
			check(!diagnostics.contains("\nTime:"), "Missing coverage reached the final report: "+work);
		}
	}

	/** Covered input keeps the exact seven-column logging row and normal output accounting. */
	private static void coveredControl(final Path root, final Path launcher) throws Exception{
		final Path work=Files.createDirectory(root.resolve("covered-control"));
		final Path input=work.resolve("assembly.fa"), coverage=work.resolve("coverage.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa"), log=work.resolve("results.tsv");
		final String fasta=">covered\nACGTACGT\n";
		final String expectedLog=LOG_HEADER+"\nassembly\tcovered\t0\t8\t10.00\t5\t100.00\n";
		write(input, fasta);
		write(coverage, COVERAGE_HEADER+"\n"+coveredRow("covered")+"\n");
		write(work.resolve("expected-log.tsv"), expectedLog);
		write(work.resolve("expected-clean.fa"), fasta);
		write(work.resolve("expected-dirty.fa"), "");

		final int exit=run(launcher, baseArguments(input, coverage, clean, dirty, log), work);

		check(exit==0, "Covered logging control exit="+exit+"; inspect "+work);
		checkNoTimeout(work);
		check(Files.size(work.resolve("stdout"))==0, "Covered logging wrote stdout: "+work);
		check(fasta.equals(read(clean)), "Covered clean output differs: "+work);
		check("".equals(read(dirty)), "Covered dirty output differs: "+work);
		check(expectedLog.equals(read(log)), "Covered log differs: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(stat(diagnostics, "Reads In", 1), "Covered control Reads In differs: "+work);
		check(stat(diagnostics, "Reads Out", 1), "Covered control Reads Out differs: "+work);
		check(stat(diagnostics, "Reads Filtered", 0), "Covered control Reads Filtered differs: "+work);
		check(diagnostics.contains("\nTime:"), "Covered control lacks the final report: "+work);
	}

	/** Captures log=cov collisions and read/read cov0=cov aliasing. */
	private static void logCollision(final Path root, final Path launcher,
			final boolean collide, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(collide ? "log-cov-collision" : "log-distinct-control"));
		final Path input=work.resolve("assembly.fa"), coverage=work.resolve("coverage.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa");
		final Path log=collide ? coverage : work.resolve("results.tsv");
		final String fasta=">covered\nACGTACGT\n";
		final String originalCoverage=COVERAGE_HEADER+"\n"+coveredRow("covered")+"\n";
		final String expectedLog=LOG_HEADER+"\tavgFold0\treads0\tnormRatio\n"
				+"assembly\tcovered\t0\t8\t10.00\t5\t100.00\t10.00\t5\t1.00\n";
		final boolean rejected=collide && fixed;
		final String expectedCoverage=rejected || !collide ? originalCoverage : expectedLog;
		write(input, fasta);
		write(coverage, originalCoverage);
		write(work.resolve("original-coverage.tsv"), originalCoverage);
		write(work.resolve("expected-clean.fa"), rejected ? "" : fasta);
		write(work.resolve("expected-dirty.fa"), "");
		write(work.resolve("expected-coverage-after.tsv"), expectedCoverage);
		if(!rejected){write(work.resolve("expected-log.tsv"), expectedLog);}

		final ArrayList<String> arguments=baseArguments(input, coverage, clean, dirty, log, "t");
		arguments.add("cov0="+coverage);
		final int exit=run(launcher, arguments, work);

		check(exit==(rejected ? 1 : 0), "Log collision exit="+exit+"; inspect "+work);
		checkNoTimeout(work);
		check(Files.size(work.resolve("stdout"))==0, "Log collision wrote stdout: "+work);
		check(expectedCoverage.equals(read(coverage)), "Log collision coverage-after differs: "+work);
		if(rejected){
			check(!Files.exists(clean), "Rejected collision still wrote clean output: "+work);
			check(!Files.exists(dirty), "Rejected collision still wrote dirty output: "+work);
		}else{
			check(fasta.equals(read(clean)), "Log collision clean output differs: "+work);
			check("".equals(read(dirty)), "Log collision dirty output differs: "+work);
			check(expectedLog.equals(read(log)), "Log collision log differs: "+work);
		}
		final String diagnostics=read(work.resolve("stderr"));
		if(rejected){
			check(diagnostics.contains("specified for both input and output streams"),
					"Rejected collision lacks the clear path diagnostic: "+work);
			check(!diagnostics.contains("\nTime:"), "Rejected collision reached the final report: "+work);
		}else{
			check(stat(diagnostics, "Reads In", 1), "Log collision Reads In differs: "+work);
			check(stat(diagnostics, "Reads Out", 1), "Log collision Reads Out differs: "+work);
			check(stat(diagnostics, "Reads Filtered", 0), "Log collision Reads Filtered differs: "+work);
			check(diagnostics.contains("\nTime:"), "Log collision lacks the final report: "+work);
		}
	}

	/** Fixed-only cases for roles not covered by the captured baseline pair. */
	private static void fixedPreflightCollision(final Path root, final Path launcher, final String name,
			final boolean logClean, final boolean qfoutCoverage, final String diagnostic) throws Exception{
		final Path work=Files.createDirectory(root.resolve(name));
		final Path input=work.resolve("assembly.fa"), coverage=work.resolve("coverage.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa");
		final Path log=logClean ? clean : work.resolve("results.tsv");
		final String originalCoverage=COVERAGE_HEADER+"\n"+coveredRow("covered")+"\n";
		write(input, ">covered\nACGTACGT\n");
		write(coverage, originalCoverage);
		write(work.resolve("original-coverage.tsv"), originalCoverage);
		write(work.resolve("expected-coverage-after.tsv"), originalCoverage);
		final ArrayList<String> arguments=baseArguments(input, coverage, clean, dirty, log, "t");
		arguments.add("cov0="+coverage);
		if(qfoutCoverage){arguments.add("qfout="+coverage);}

		final int exit=run(launcher, arguments, work);

		check(exit==1, "Fixed preflight collision exit="+exit+"; inspect "+work);
		checkNoTimeout(work);
		check(Files.size(work.resolve("stdout"))==0, "Fixed preflight collision wrote stdout: "+work);
		check(originalCoverage.equals(read(coverage)), "Fixed preflight collision changed coverage: "+work);
		check(!Files.exists(clean), "Fixed preflight collision wrote clean output: "+work);
		check(!Files.exists(dirty), "Fixed preflight collision wrote dirty output: "+work);
		check(!Files.exists(log), "Fixed preflight collision wrote log output: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(diagnostics.contains(diagnostic), "Fixed preflight collision lacks the clear diagnostic: "+work);
		check(!diagnostics.contains("\nTime:"), "Fixed preflight collision reached the final report: "+work);
	}

	private static ArrayList<String> baseArguments(final Path input, final Path coverage, final Path clean,
			final Path dirty, final Path log){
		return baseArguments(input, coverage, clean, dirty, log, "f");
	}
	private static ArrayList<String> baseArguments(final Path input, final Path coverage, final Path clean,
			final Path dirty, final Path log, final String overwrite){
		return new ArrayList<String>(Arrays.asList("-Xmx128m", "-Xms32m", "in="+input, "cov="+coverage,
				"out="+clean, "outdirty="+dirty, "log="+log, "mincoverage=5", "minpercent=40",
				"minreads=0", "minlength=1", "trim=0", "fastawrap=1000", "t=1", "int=f", "ow="+overwrite,
				"simd=f"));
	}

	private static String coveredRow(final String id){
		return id+"\t10\t8\t0.5\t100\t8\t3\t2\t0.5\t10\t0.0\t0";
	}

	private static boolean stat(final String diagnostics, final String label, final long value){
		return Pattern.compile("(?m)^"+Pattern.quote(label)+":\\s+"+value+"(?:\\s|$)").matcher(diagnostics).find();
	}
	private static void checkNoTimeout(final Path work) throws Exception{
		check(!Files.exists(work.resolve("timeout.txt")), "Owned filter launcher timed out: "+work);
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}
	private static String read(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Supervises only this invocation's owned process group; saves the complete literal command and diagnostics. */
	private static int run(final Path launcher, final List<String> arguments, final Path work) throws Exception{
		final ArrayList<String> command=launcherCommand(launcher, arguments);
		Files.write(work.resolve("args.txt"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).directory(launcher.getParent().toFile())
				.redirectOutput(work.resolve("stdout").toFile()).redirectError(work.resolve("stderr").toFile()).start();
		Throwable primary=null;
		try{
			process.getOutputStream().close();
			if(!process.waitFor(60, TimeUnit.SECONDS)){
				Files.write(work.resolve("timeout.txt"), Collections.singletonList("Owned launcher exceeded60 seconds."),
						StandardCharsets.US_ASCII);
				throw new AssertionError("Finite FilterByCoverage logging baseline exceeded60 seconds: "+work);
			}
			final int exit=process.exitValue();
			write(work.resolve("exit.txt"), exit+"\n");
			return exit;
		}catch(Exception e){
			primary=e; throw e;
		}catch(Error e){
			primary=e; throw e;
		}finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
			catch(Error e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
		}
	}

	private static ArrayList<String> launcherCommand(final Path launcher, final List<String> arguments){
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/timeout", "--verbose",
				"--signal=TERM", "--kill-after=2s", "30s", "bash", launcher.toString()));
		command.addAll(arguments);
		return command;
	}

	/** Allows GNU timeout's two-second group KILL grace before forcing/reaping its exact owned supervisor. */
	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 5)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){
			throw new AssertionError("Owned FilterByCoverage launcher did not terminate after forced cleanup.");
		}
	}
	/** Keeps repeated interrupts from removing the cleanup deadline; restores interruption on return. */
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

	private static final String COVERAGE_HEADER="#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases"
			+"\tPlus_reads\tMinus_reads\tRead_GC\tMedian_fold\tStd_Dev\tUnder_min";
	private static final String LOG_HEADER="#assembly\tcontig\tcontam\tlength\tavgFold\treads\tpercentCovered";
	private static int checks;
}
