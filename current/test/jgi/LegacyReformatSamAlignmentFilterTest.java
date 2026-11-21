package test.jgi;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

import parse.Parse;

/** Constructor-only legacy SAM breaklength rejection cases for original-alignment filters.
 * @author Jean
 */
public final class LegacyReformatSamAlignmentFilterTest{

	/** Runs fresh helper JVMs so legacy constructor statics cannot cross-contaminate cases. */
	public static void main(final String[] args) throws Exception{
		final Parsed parsed=parse(args);
		if(parsed.helper){
			construct(parsed.driver, parsed.driverArgs());
			return;
		}
		require("baseline".equals(parsed.phase) || "fixed".equals(parsed.phase),
			"Set phase=baseline|fixed.");
		require(parsed.directory!=null && !parsed.directory.isEmpty(), "Set a new outdir.");
		require(parsed.rootPath!=null && !parsed.rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/legacyreformatsamalignmentfiltertest.sh with assertions enabled.");
		final LegacyReformatSamAlignmentFilterTest test=new LegacyReformatSamAlignmentFilterTest(
			Paths.get(parsed.directory).toAbsolutePath().normalize(),
			Paths.get(parsed.rootPath).toAbsolutePath().normalize(), parsed.phase.equals("baseline"));
		for(final Case item : CASES){test.run(item);}
		test.check(test.cases==CASES.length, "Expected "+CASES.length+" cases, saw "+test.cases);
		test.check(test.children==CASES.length, "Expected "+CASES.length+" children, saw "+test.children);
		test.check(test.changed==7, "Expected 7 repaired rows, saw "+test.changed);
		test.check(test.controls==6, "Expected six unchanged controls, saw "+test.controls);
		test.check(test.sourceMatches==(test.baseline ? 13 : 6), "Unexpected source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==(test.baseline ? 6 : 13), "Unexpected corrected matches: "+test.correctedMatches);
		final String summary="LEGACY_REFORMAT_SAM_ALIGNMENT_FILTER_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+parsed.phase+" cases="+test.cases+" children="+test.children+
			" changed="+test.changed+" controls="+test.controls+
			" sourceMatches="+test.sourceMatches+" correctedMatches="+test.correctedMatches+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+test.out);}
		System.out.print(summary);
	}

	/** Binds the test root and initializes the evidence table. */
	private LegacyReformatSamAlignmentFilterTest(final Path out_, final Path repo_, final boolean baseline_)
			throws Exception{
		out=out_; repo=repo_; baseline=baseline_;
		launcher=repo.resolve("dev/legacyreformatsamalignmentfiltertest.sh");
		require(Files.isDirectory(repo.resolve("current")), "Missing current/ at "+repo);
		require(Files.isRegularFile(launcher), "Missing helper launcher at "+launcher);
		Files.createDirectory(out);
		Files.write(out.resolve("cases.tsv"), ("case\tdriver\tphase\texit\toutputExists"+
			"\tsourceMatches\tcorrectedMatches\tstderrMatches\tstdoutEmpty\tinputsUnchanged\tmatchPhase\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Launches one constructor-only helper JVM and checks exact fixed rejection or success. */
	private void run(final Case item) throws Exception{
		cases++;
		children++;
		if(item.source.changed(item.corrected)){changed++;}else{controls++;}
		final Path dir=out.resolve(item.name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.sam"), output=dir.resolve("output.sam");
		Files.write(input, SAM.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final byte[] original=Files.readAllBytes(input);
		check(!Files.exists(output), item.name+": output path existed before helper");
		final ArrayList<String> command=command(item, input, output);
		Files.write(dir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		Files.write(dir.resolve("expected.tsv"), expected(item).getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		final Path stdout=dir.resolve("stdout.txt"), stderr=dir.resolve("stderr.txt");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdout.toFile()).redirectError(stderr.toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		Files.write(dir.resolve("exit.txt"), (exit+"\n").getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		final byte[] stdoutBytes=Files.readAllBytes(stdout);
		final String stderrText=new String(Files.readAllBytes(stderr), StandardCharsets.UTF_8);
		final boolean outputExists=Files.exists(output);
		final boolean sourceMatch=item.source.matches(exit, outputExists, stderrText);
		final boolean correctedMatch=item.corrected.matches(exit, outputExists, stderrText);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean stderrMatches=(baseline ? item.source : item.corrected).stderrMatches(stderrText);
		final boolean stdoutEmpty=stdoutBytes.length==0;
		final boolean unchanged=Arrays.equals(original, Files.readAllBytes(input));
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;
		Files.write(out.resolve("cases.tsv"), (item.name+"\t"+item.driver.simpleName+"\t"+
			(baseline ? "baseline" : "fixed")+"\t"+exit+"\t"+outputExists+"\t"+sourceMatch+
			"\t"+correctedMatch+"\t"+stderrMatches+"\t"+stdoutEmpty+"\t"+unchanged+
			"\t"+phaseMatch+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
		check(phaseMatch, item.name+": wrong constructor result; stderr="+stderrText);
		check(stderrMatches, item.name+": stderr mismatch: "+stderrText);
		check(stdoutEmpty, item.name+": helper wrote to stdout");
		check(unchanged, item.name+": input bytes changed");
		check(!outputExists, item.name+": constructor created primary output");
		check(command.contains("breaklength="+(item.breakLength ? 3 : 0)) &&
			command.contains("int=f") && command.contains("t=1") &&
			command.contains("wt=1") && command.contains("simd=f") &&
			command.contains("qin=33") && command.contains("qout=33") &&
			!command.contains("sample=1") && !command.contains("samplebases=1"),
			item.name+": helper argv lost a safe constructor-only SAM-to-SAM pin");
	}

	/** Builds one helper JVM command. */
	private ArrayList<String> command(final Case item, final Path input, final Path output){
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("helper=t");
		command.add("driver="+item.driver.className);
		command.add("in="+input);
		command.add("out="+output);
		command.add("t=1");
		command.add("wt=1");
		command.add("int=f");
		command.add("ow=f");
		command.add("app=f");
		command.add("simd=f");
		command.add("minlength=1");
		command.add("breaklength="+(item.breakLength ? 3 : 0));
		command.add("qin=33");
		command.add("qout=33");
		for(final String flag : item.flags){command.add(flag);}
		return command;
	}

	/** Reflects the real public constructor and never calls driver main() or process(). */
	private static void construct(final String className, final String[] driverArgs) throws Exception{
		require(className!=null && className.startsWith("jgi.Reformat"), "Set driver=jgi.Reformat*.");
		final Class<?> type=Class.forName(className);
		final Constructor<?> constructor=type.getConstructor(String[].class);
		try{constructor.newInstance((Object)driverArgs);}
		catch(InvocationTargetException e){
			e.getCause().printStackTrace();
			System.exit(1);
		}
	}

	/** Parses parent and helper flags, leaving all driver flags in order. */
	private static Parsed parse(final String[] args){
		final Parsed parsed=new Parsed();
		final ArrayList<String> driverArgs=new ArrayList<String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){parsed.phase=pair[1];}
			else if(pair[0].equals("outdir")){parsed.directory=pair[1];}
			else if(pair[0].equals("bbtools")){parsed.rootPath=pair[1];}
			else if(pair[0].equals("helper")){parsed.helper=Parse.parseBoolean(pair[1]);}
			else if(pair[0].equals("driver")){parsed.driver=pair[1];}
			else{driverArgs.add(arg);}
		}
		parsed.driverArgs=driverArgs.toArray(new String[driverArgs.size()]);
		return parsed;
	}

	/** Expected constructor results for evidence transcripts. */
	private static String expected(final Case item){
		return "phase\tresult\nsource\t"+item.source+"\ncorrected\t"+item.corrected+"\n";
	}

	/** Shell-escaped argv transcript for the evidence directory. */
	private static String join(final ArrayList<String> command){
		final StringBuilder builder=new StringBuilder();
		for(final String arg : command){
			if(builder.length()>0){builder.append(' ');}
			builder.append('\'').append(arg.replace("'", "'\\''")).append('\'');
		}
		return builder.append('\n').toString();
	}

	/** Records one assertion while retaining later row evidence. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	/** Hard construction-time requirement. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static Case reject(final Driver driver, final String branch, final String flag){
		return new Case(driver.simpleName+"-"+branch, driver, true, new String[]{flag},
			Success.INSTANCE, new Reject(driver));
	}

	private static Case control(final Driver driver, final String name, final boolean breakLength,
			final String[] flags){
		return new Case(driver.simpleName+"-"+name, driver, breakLength, flags,
			Success.INSTANCE, Success.INSTANCE);
	}

	/** One helper-only constructor case. */
	private static final class Case{
		Case(final String name_, final Driver driver_, final boolean breakLength_,
				final String[] flags_, final Expect source_, final Expect corrected_){
			name=name_;
			driver=driver_;
			breakLength=breakLength_;
			flags=flags_;
			source=source_;
			corrected=corrected_;
		}
		final String name;
		final Driver driver;
		final boolean breakLength;
		final String[] flags;
		final Expect source;
		final Expect corrected;
	}

	/** Legacy driver identity and exact fixed-constructor stack frame. */
	private static final class Driver{
		Driver(final String simpleName_, final String className_){
			simpleName=simpleName_;
			className=className_;
		}
		final String simpleName;
		final String className;
	}

	/** Parsed flags for parent or helper mode. */
	private static final class Parsed{
		String[] driverArgs(){return driverArgs;}
		boolean helper;
		String phase;
		String directory;
		String rootPath;
		String driver;
		String[] driverArgs;
	}

	/** Expected constructor result for one phase. */
	private interface Expect{
		boolean matches(int exit, boolean outputExists, String stderr);
		boolean stderrMatches(String stderr);
		boolean changed(Expect other);
	}

	/** Baseline or fixed control constructor success. */
	private static final class Success implements Expect{
		public boolean matches(final int exit, final boolean outputExists, final String stderr){
			return exit==0 && !outputExists && stderrMatches(stderr);
		}
		public boolean stderrMatches(final String stderr){
			return !stderr.contains("Exception in thread") && !stderr.contains("Unknown parameter") &&
				!stderr.contains("java.lang.") && !stderr.contains("ERROR") && !stderr.contains("Error");
		}
		public boolean changed(final Expect other){return other!=this && !toString().equals(other.toString());}
		@Override
		public String toString(){return "constructor-success-no-primary-output";}
		static final Success INSTANCE=new Success();
	}

	/** Future constructor rejection before the output writer creates the SAM destination. */
	private static final class Reject implements Expect{
		Reject(final Driver driver_){driver=driver_;}
		public boolean matches(final int exit, final boolean outputExists, final String stderr){
			return exit==1 && !outputExists && stderrMatches(stderr);
		}
		public boolean stderrMatches(final String stderr){
			return stderr.contains("java.lang.IllegalArgumentException: "+DIAGNOSTIC+"\n") &&
				stderr.contains("\tat "+driver.className+".<init>("+driver.simpleName+".java:");
		}
		public boolean changed(final Expect other){return other!=this;}
		@Override
		public String toString(){return "fixed-preflight-rejection";}
		final Driver driver;
	}

	/*--------------------------------------------------------------*/

	private final Path out;
	private final Path repo;
	private final Path launcher;
	private final boolean baseline;
	private int cases;
	private int children;
	private int changed;
	private int controls;
	private int sourceMatches;
	private int correctedMatches;
	private int checks;
	private int failures;

	private static final Driver REFORMAT2=new Driver("Reformat2", "jgi.Reformat2");
	private static final Driver REFORMAT_READS=new Driver("ReformatReads", "jgi.ReformatReads");
	private static final Driver REFORMAT_READS_OLD=new Driver("ReformatReadsOld", "jgi.ReformatReadsOld");
	private static final String DIAGNOSTIC=
		"breaklength with SAM/BAM input is incompatible with alignment-dependent filters.";
	private static final String SAM="@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\n"+
		"mapped\t0\tchr1\t1\t60\t8=\t*\t0\t0\tACGTACGT\tIIIIIIII\n";
	private static final Case[] CASES={
		reject(REFORMAT2, "minmapq", "minmapq=60"),
		reject(REFORMAT_READS, "minmapq", "minmapq=60"),
		reject(REFORMAT_READS, "idfilter", "idfilter=0.90"),
		reject(REFORMAT_READS, "mappedonly", "mappedonly=t"),
		reject(REFORMAT_READS_OLD, "minmapq", "minmapq=60"),
		reject(REFORMAT_READS_OLD, "idfilter", "idfilter=0.90"),
		reject(REFORMAT_READS_OLD, "mappedonly", "mappedonly=t"),
		control(REFORMAT2, "no-break-filter", false, new String[]{"minmapq=60"}),
		control(REFORMAT2, "breaklength-no-filter", true, new String[0]),
		control(REFORMAT_READS, "no-break-filter", false, new String[]{"minmapq=60"}),
		control(REFORMAT_READS, "breaklength-no-filter", true, new String[0]),
		control(REFORMAT_READS_OLD, "no-break-filter", false, new String[]{"minmapq=60"}),
		control(REFORMAT_READS_OLD, "breaklength-no-filter", true, new String[0])
	};
}
