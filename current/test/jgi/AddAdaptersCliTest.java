package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import parse.Parse;

/**
 * Twelve fixed-only AddAdapters CLI cases with literal generation bytes and grading counts.
 * Generates private inputs, invokes the actual launcher in owned supervised children and
 * snapshots inputs only after child termination. Covers deliberately wrong trimmed lengths,
 * both header families and decoder failures without executing an incompatible baseline.
 * Percentage formatting, including existing zero-denominator values, is not redefined.
 * @author Sucrose
 */
public final class AddAdaptersCliTest{

	/** Requires a canonical {@code tree=...} and new {@code outdir=...}; accepts only {@code phase=fixed}. */
	public static void main(final String[] args) throws Exception{
		String sourceTree=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("tree")){sourceTree=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(!pair[0].equals("phase") || !pair[1].equals("fixed")){
				throw new IllegalArgumentException("Fixed-only CLI fixture: "+arg);
			}
		}
		require(sourceTree!=null && directory!=null, "Set tree and a new outdir.");
		final Path tree=Paths.get(sourceTree).toAbsolutePath().normalize();
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		require(safePath(tree) && safePath(root), "Paths must survive AddAdapters equals splitting and fit one evidence row.");
		require(Files.isRegularFile(tree.resolve("addadapters.sh")), "Missing actual AddAdapters launcher: "+tree);
		require(Files.isExecutable(Paths.get("/usr/bin/timeout")), "Owned child supervision requires GNU timeout.");
		Files.createDirectory(root);
		final AddAdaptersCliTest test=new AddAdaptersCliTest(tree, root);
		for(int i=0; i<CASES.length; i++){test.scenario(CASES[i], i==0);}
		test.check(test.cases==12 && test.children==12, "Frozen CLI matrix requires12 cases/12 terminal children: "+test.cases+"/"+test.children);
		if(test.failures>0){throw new AssertionError("AddAdapters CLI failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADD_ADAPTERS_CLI_PASS phase=fixed cases="+test.cases+" children="+test.children+" checks="+test.checks);
	}

	private AddAdaptersCliTest(final Path tree_, final Path root_){tree=tree_; root=root_;}

	/** All children process their complete tiny inputs; generation alone receives a data-output path. */
	private void scenario(final CliCase row, final boolean generation) throws Exception{
		final Path work=Files.createDirectory(root.resolve(row.name));
		final Path input1=work.resolve("input1.fq"), input2=work.resolve("input2.fq"), output=work.resolve("generated.fq");
		write(input1, row.first); write(work.resolve("input1-before.fq"), row.first);
		if(row.second!=null){write(input2, row.second); write(work.resolve("input2-before.fq"), row.second);}
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("-ea", "-Xmx256m", "-Xms32m", "in="+input1,
			"mode="+(generation ? "write" : "grade"), "reads=-1", "minlen=1", "int=f", "qin=33", "qout=33",
			"simd=f", "t=1", "ow=t", "app=f"));
		if(row.second!=null){args.add("in2="+input2);}
		if(row.custom!=null){args.add("parsecustom="+row.custom);}
		if(generation){
			args.addAll(Arrays.asList("out="+output, "literal=AC", "rate=1", "adderrors=f", "arc=f", "right=t"));
			write(work.resolve("expected-generation.fq"), GENERATED);
		}
		write(work.resolve("expected.tsv"), "totalReads\ttotalBases\tcorrectReads\tcorrectBases\trejectedHeader\n"+
			row.totalReads+"\t"+row.totalBases+"\t"+row.correctReads+"\t"+row.correctBases+"\t"+
			(row.badHeader==null ? "none" : row.badHeader)+"\n");
		final int exit=run(args, work);
		//run has reaped the exact owned supervisor before any input observation.
		final boolean same1=saveInput(work, input1, "input1-after.fq", row.first);
		final boolean same2=row.second==null || saveInput(work, input2, "input2-after.fq", row.second);
		final String diagnostics=read(work.resolve("stderr"));
		final long[] total=summary(diagnostics, "Total output:"), correct=summary(diagnostics, "Perfectly Correct (% of output):");
		write(work.resolve("observed.tsv"), "exit\ttotalReads\ttotalBases\tcorrectReads\tcorrectBases\tinput1Unchanged\tinput2Unchanged\n"+
			exit+"\t"+total[0]+"\t"+total[1]+"\t"+correct[0]+"\t"+correct[1]+"\t"+same1+"\t"+
			(row.second==null ? "NA" : Boolean.toString(same2))+"\n");
		if(Files.exists(work.resolve("timeout.txt"))){throw new AssertionError("Timeout/forced cleanup is never an accepted CLI rejection: "+work);}
		check(same1 && same2, row.name+": every original input byte must survive the child unchanged.");
		check(Files.size(work.resolve("stdout"))==0, row.name+": file/statistics mode unexpectedly wrote stdout.");
		if(row.badHeader!=null){
			check(exit==1, row.name+": decoder rejection must be a real Java exit1, observed="+exit);
			final Pattern primary=Pattern.compile("(?m)^(?:Exception in thread \\\"main\\\" )?java\\.lang\\.IllegalArgumentException: "+
				Pattern.quote("Cannot grade read header '"+row.badHeader+"':"));
			check(primary.matcher(diagnostics).find() && diagnostics.contains(row.reason),
				row.name+": primary decoder failure/header/reason missing; inspect stderr.");
			check(!hasGradingSummary(diagnostics), row.name+": malformed input printed a partial grading success summary.");
			check(!Files.exists(output), row.name+": grade failure created data output.");
		}else{
			check(exit==0, row.name+": expected successful fixed CLI exit, observed="+exit);
			check(!diagnostics.contains("Exception") && !diagnostics.contains("AssertionError"), row.name+": successful CLI contains an exception.");
			if(generation){
				check(matches(output, GENERATED), "Ordinary generation must emit exact0_1_0/A/I FASTQ bytes.");
				checkCount(diagnostics, "Adapters Added:", 1, 1, row.name);
				checkCount(diagnostics, "Valid Output:", 0, 0, row.name);
				check(!hasGradingSummary(diagnostics), "The generator must not claim a grading result.");
			}else{
				checkCount(diagnostics, "Total output:", row.totalReads, row.totalBases, row.name);
				checkCount(diagnostics, "Perfectly Correct (% of output):", row.correctReads, row.correctBases, row.name);
				check(!Files.exists(output), row.name+": grading unexpectedly created data output.");
			}
		}
		cases++;
	}

	/** Reads counts only; percentage fields deliberately retain the program's existing policy. */
	private static Pattern countPattern(final String label){
		return Pattern.compile("(?m)^"+Pattern.quote(label)+"[\\t ]+([0-9]+) reads[^\\r\\n]*?[\\t ]+([0-9]+) bases[^\\r\\n]*$");
	}
	private static long[] summary(final String text, final String label){
		final Matcher matcher=countPattern(label).matcher(text);
		return matcher.find() ? new long[]{Long.parseLong(matcher.group(1)), Long.parseLong(matcher.group(2))} : new long[]{-1, -1};
	}
	private void checkCount(final String text, final String label, final long reads, final long bases, final String name){
		final Matcher matcher=countPattern(label).matcher(text);
		final boolean found=matcher.find();
		check(found && Long.parseLong(matcher.group(1))==reads && Long.parseLong(matcher.group(2))==bases,
			name+": expected "+label+" "+reads+" reads/"+bases+" bases.");
		check(!found || !matcher.find(), name+": duplicate "+label+" summaries are not one completed result.");
	}
	private static boolean hasGradingSummary(final String text){
		return text.contains("Total output:") || text.contains("Perfectly Correct (% of output):") || text.contains("Incorrect (% of output):");
	}

	/** Uses only the supervisor handle returned by this launch; GNU timeout owns its isolated child group. */
	private int run(final List<String> args, final Path work) throws Exception{
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/timeout", "--verbose",
			"--signal=TERM", "--kill-after=2s", "30s", "bash", tree.resolve("addadapters.sh").toString()));
		command.addAll(args);
		Files.write(work.resolve("command.argv"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(tree.toFile())
			.redirectOutput(work.resolve("stdout").toFile()).redirectError(work.resolve("stderr").toFile());
		final String inherited=builder.environment().get("JAVA_TOOL_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", (inherited==null ? "" : inherited+" ")+"-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process process=builder.start();
		Throwable primary=null;
		boolean exceeded=false;
		try{
			process.getOutputStream().close();
			exceeded=!process.waitFor(45, TimeUnit.SECONDS);
		}catch(Exception e){primary=e; throw e;}
		catch(Error e){primary=e; throw e;}
		finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){if(primary==null){throw e;}else{primary.addSuppressed(e);}}
			catch(Error e){if(primary==null){throw e;}else{primary.addSuppressed(e);}}
		}
		final int exit=process.exitValue();
		write(work.resolve("exit.txt"), exit+"\n");
		if(exceeded || exit==124 || exit==137 || read(work.resolve("stderr")).contains("timeout: sending signal")){
			write(work.resolve("timeout.txt"), "Owned CLI child exceeded its deadline or was killed.\n");
		}
		children++;
		return exit;
	}
	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 5)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){throw new AssertionError("Owned AddAdapters CLI supervisor survived forced cleanup.");}
	}
	/** Interrupts cannot cancel the cleanup deadline; restores the flag after the owned process is reaped. */
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
	private static boolean saveInput(final Path work, final Path path, final String name, final String expected) throws Exception{
		if(!Files.isRegularFile(path)){write(work.resolve(name+".missing"), "Generated input is no longer a regular file.\n"); return false;}
		final byte[] actual=Files.readAllBytes(path);
		Files.write(work.resolve(name), actual, StandardOpenOption.CREATE_NEW);
		return Arrays.equals(actual, expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static boolean matches(final Path path, final String expected) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static String record(final String id, final String bases){
		final char[] quality=new char[bases.length()]; Arrays.fill(quality, 'I');
		return "@"+id+"\n"+bases+"\n+\n"+new String(quality)+"\n";
	}
	private static String read(final Path path) throws Exception{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static boolean safePath(final Path path){
		final String text=path.toString();
		return text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0 && text.indexOf('\r')<0;
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}}

	/** Each row owns literal inputs and an independent count oracle; null custom means the option is omitted. */
	private static final class CliCase{
		CliCase(final String name_, final String first_, final String second_, final String custom_,
			final long totalReads_, final long totalBases_, final long correctReads_, final long correctBases_){
			this(name_, first_, second_, custom_, totalReads_, totalBases_, correctReads_, correctBases_, null, null);
		}
		CliCase(final String name_, final String first_, final String badHeader_, final String reason_){
			this(name_, first_, null, null, -1, -1, -1, -1, badHeader_, reason_);
		}
		private CliCase(final String name_, final String first_, final String second_, final String custom_,
			final long totalReads_, final long totalBases_, final long correctReads_, final long correctBases_,
			final String badHeader_, final String reason_){
			name=name_; first=first_; second=second_; custom=custom_;
			totalReads=totalReads_; totalBases=totalBases_; correctReads=correctReads_; correctBases=correctBases_;
			badHeader=badHeader_; reason=reason_;
		}
		final String name, first, second, custom, badHeader, reason;
		final long totalReads, totalBases, correctReads, correctBases;
	}
	private static final String GENERATED="@0_1_0\nA\n+\nI\n";
	private static final String SYN="SYN_7_0_5_7_+_100_1_._ref&7_0_7_7_-_200_1_._ref";
	private static final String SYN_INPUT=record(SYN+" 1:", "GGGGGG")+record(SYN+" 2:", "GGGGGGG");
	private static final CliCase[] CASES={
		new CliCase("ordinary-generation", record("ordinary", "G"), null, null, -1, -1, -1, -1),
		new CliCase("legacy-mixed", record("7_6_2", "GG")+record("8_6_2", "GGG")+
			record("9_6_4", "GG")+record("10_6_6", "GGGGGG"), null, null, 4, 13, 2, 8),
		new CliCase("legacy-zero", record("7_6_0", "G"), null, "f", 1, 1, 1, 0),
		new CliCase("shared-suffix-pair", record("7_6_2_8_5 /1", "GG"), record("7_6_2_8_5 /2", "GGGGG"), null, 2, 7, 2, 7),
		new CliCase("shared-context-pair", record("7_6_2_8_5", "GG"), record("7_6_2_8_5", "GGGGG"), null, 2, 7, 2, 7),
		new CliCase("rename-pair", record("7_6_2 1:", "GG"), record("7_8_5 2:", "GGGGG"), null, 2, 7, 2, 7),
		new CliCase("orphan-second", record("7_6_2_8_5 /2", "GGGGG"), null, null, 1, 5, 1, 5),
		new CliCase("synthetic-default", SYN_INPUT, null, null, 2, 13, 2, 13),
		new CliCase("synthetic-false", SYN_INPUT, null, "f", 2, 13, 2, 13),
		new CliCase("synthetic-true", SYN_INPUT, null, "t", 2, 13, 2, 13),
		new CliCase("good-then-malformed", record("7_6_2", "GG")+record("8_6_x", "GG"), "8_6_x", "nondecimal character"),
		new CliCase("unpaired-ambiguous", record("7_6_2_8_5", "GG"), "7_6_2_8_5", "mate suffix or paired input")
	};
	private final Path tree, root;
	private int checks, failures, cases, children;
}
