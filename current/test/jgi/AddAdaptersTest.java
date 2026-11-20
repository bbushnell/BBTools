package test.jgi;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import jgi.AddAdapters;
import parse.Parse;
import stream.Read;

/**
 * Finite end-selection and exact input/output-path collision checks for deprecated AddAdapters.
 * Baseline records the six known selector mismatches and observes the two unguarded input2
 * collisions without assuming which asynchronous read/write ordering wins. Fixed mode requires
 * the intended selectors and rejection of all five collision cases. Grading is never invoked.
 * Inputs are generated only beneath a newly created output directory and are never deleted.
 * @author Sucrose
 */
public final class AddAdaptersTest{

	/** Requires {@code phase=baseline|fixed}, {@code tree=...}, and a new {@code outdir=...}. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, sourceTree=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("tree")){sourceTree=pair[1];}
			else{throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && sourceTree!=null, "Set tree and a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Collision baselines must exercise AddAdapters assertions; use -ea.");
		final Path tree=Paths.get(sourceTree).toAbsolutePath().normalize();
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		require(safePath(tree) && safePath(root), "Fixture paths cannot contain equals signs, tabs or newlines.");
		require(Files.isRegularFile(tree.resolve("addadapters.sh")), "Missing deprecated launcher: "+tree);
		require(Files.isExecutable(Paths.get("/usr/bin/timeout")), "Finite child supervision requires GNU timeout.");
		Files.createDirectory(root);
		final AddAdaptersTest test=new AddAdaptersTest(tree, root, phase.equals("baseline"));
		test.helpers();
		test.cliEnds();
		for(int i=0; i<7; i++){test.alias(i);}
		test.check(test.helperCases==14 && test.cliCases==3 && test.aliasCases==7,
			"Frozen case counts must be14 helpers/3 CLI selectors/7 alias-controls: "+
			test.helperCases+"/"+test.cliCases+"/"+test.aliasCases);
		test.check(test.selectorMismatches==(test.baseline ? 6 : 0),
			"Only the six source-confirmed helper selectors differ from intended behavior before repair: "+test.selectorMismatches);
		if(test.failures>0){throw new AssertionError("AddAdapters fixture failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADD_ADAPTERS_PASS phase="+phase+" helperCases=14 cliCases=3 aliasCases=7"+
			" selectorMismatches="+test.selectorMismatches+" checks="+test.checks);
	}

	private AddAdaptersTest(final Path tree_, final Path root_, final boolean baseline_){
		tree=tree_; root=root_; baseline=baseline_;
	}

	/** Calls the real private helper at loc2; the six-base adapter prevents random filler on either side. */
	private void helpers() throws Exception{
		final Path work=Files.createDirectory(root.resolve("helpers"));
		final Path input=work.resolve("input.fq");
		write(input, ONE);
		write(work.resolve("input-before.fq"), ONE);
		write(work.resolve("cases.tsv"), "case\tintendedBases\tintendedId\tactualBases\tactualId\tqualitiesUnchanged\tnumericId\n");
		final Method method=AddAdapters.class.getDeclaredMethod("addAdapter", Read.class, int.class);
		method.setAccessible(true);
		for(final EndCase row : ENDS){
			final ArrayList<String> args=baseArguments(input, work.resolve("unused-"+row.name+".fq"), "ACGTAC", "1");
			if(row.before!=null){args.add(row.before);}
			args.add(row.option);
			final AddAdapters app=new AddAdapters(args.toArray(new String[args.size()]));
			final byte[] quality={30, 30, 30, 30, 30, 30};
			final Read read=new Read("GGGGGG".getBytes(StandardCharsets.US_ASCII), quality, "original", 7);
			try{method.invoke(app, read, 2);}
			catch(InvocationTargetException e){throw new AssertionError("Actual addAdapter helper failed: "+row.name, e.getCause());}
			final String bases=new String(read.bases, StandardCharsets.US_ASCII);
			final boolean qualitySame=Arrays.equals(read.quality, new byte[]{30, 30, 30, 30, 30, 30});
			final String intendedBases=row.intendedRight ? "GGACGT" : "TACGGG";
			final String intendedId=row.intendedRight ? "6_2" : "6_3";
			append(work.resolve("cases.tsv"), row.name+"\t"+intendedBases+"\t"+intendedId+"\t"+
				bases+"\t"+read.id+"\t"+qualitySame+"\t"+read.numericID+"\n");
			final boolean expectedRight=baseline ? row.baselineRight : row.intendedRight;
			check(bases.equals(expectedRight ? "GGACGT" : "TACGGG") && read.id.equals(expectedRight ? "6_2" : "6_3"),
				row.name+": fixed loc2 must produce the literal phase-specific sequence/remaining length; observed "+bases+"/"+read.id);
			check(qualitySame && read.numericID==7, row.name+": adapter selection must preserve qualities and numericID7.");
			check(!Files.exists(work.resolve("unused-"+row.name+".fq")), row.name+": constructing/calling a helper must not open output.");
			if(!bases.equals(intendedBases) || !read.id.equals(intendedId)){selectorMismatches++;}
			helperCases++;
		}
		check(saveInput(work, input, "input-after.fq", ONE), "Helper constructors modified their generated input.");
	}

	/** One base forces loc0: first adapter base A is right, last adapter base C is left. */
	private void cliEnds() throws Exception{
		final String[] names={"cli-right", "cli-5prime", "cli-end-5prime"};
		final String[] options={"right=t", "5prime=t", "end=5prime"};
		for(int i=0; i<names.length; i++){
			final Path work=Files.createDirectory(root.resolve(names[i]));
			final Path input=work.resolve("input.fq"), output=work.resolve("output.fq");
			write(input, ONE); write(work.resolve("input-before.fq"), ONE);
			final ArrayList<String> args=baseArguments(input, output, "AC", "1");
			args.add(options[i]);
			final String base=i==0 || (i==2 && baseline) ? "A" : "C";
			final String expected=record("0_1_0", base, "I");
			write(work.resolve("expected.fq"), expected);
			write(work.resolve("intended.fq"), record("0_1_0", i==0 ? "A" : "C", "I"));
			final int exit=run(args, work);
			saveInput(work, input, "input-after.fq", ONE);
			checkSuccess(exit, work);
			check(matches(output, expected), names[i]+": expected exact one-base FASTQ including generated ID/quality/newline.");
			check(matches(input, ONE), names[i]+": generated input changed.");
			cliCases++;
		}
	}

	/** Uses only fresh generated twins; unsafe baseline target cases record writer evidence, not a guessed corruption. */
	private void alias(final int index) throws Exception{
		assert(index>=0 && index<7) : "The fixed seven-case alias matrix determines all paths and oracles: "+index;
		final String[] names={"out1-in2-paired", "out1-in2-interleaved", "out1-in1", "out2-in1", "out2-in2", "distinct-paired", "distinct-unpaired"};
		final Path work=Files.createDirectory(root.resolve(names[index]));
		final Path input1=work.resolve("input1.fq"), input2=work.resolve("input2.fq");
		final Path fresh1=work.resolve("output1.fq"), fresh2=work.resolve("output2.fq");
		write(input1, PAIR1); write(work.resolve("input1-before.fq"), PAIR1);
		if(index!=6){write(input2, PAIR2); write(work.resolve("input2-before.fq"), PAIR2);}
		final Path output1=index<2 ? input2 : (index==2 ? input1 : fresh1);
		final Path output2=index==1 || index==6 ? null : (index==3 ? input1 : (index==4 ? input2 : fresh2));
		final ArrayList<String> args=baseArguments(input1, output1, "AC", "0");
		if(index!=6){args.add("in2="+input2);}
		if(output2!=null){args.add("out2="+output2);}
		final int exit=run(args, work);
		final boolean same1=saveInput(work, input1, "input1-after.fq", PAIR1);
		final boolean same2=index==6 || saveInput(work, input2, "input2-after.fq", PAIR2);
		final String diagnostics=read(work.resolve("stderr"));
		final boolean collision=collisionDiagnostic(diagnostics);
		write(work.resolve("observations.tsv"), "exit\tinput1Unchanged\tinput2Unchanged\tfreshOut1Exists\tfreshOut2Exists\tcollisionDiagnostic\n"+
			exit+"\t"+same1+"\t"+(index==6 ? "NA" : same2)+"\t"+Files.exists(fresh1)+"\t"+Files.exists(fresh2)+"\t"+collision+"\n");
		noTimeout(work);
		if(index<2 && baseline){
			check(exit==0 || exit==1, names[index]+": baseline writer observation requires a normal Java exit, not signal/supervisor failure: "+exit);
			check(!collision, names[index]+": baseline unexpectedly rejected the missing input2 alias; inspect source/runtime identity.");
			check(same1, names[index]+": the unrelated first input changed.");
			check(!same2 || Files.exists(fresh2), names[index]+": no actual writer evidence; an unrelated failure cannot reproduce the missing guard.");
		}else if(index<5){
			check(exit==1 && collision, names[index]+": require explicit Java collision rejection, observed exit="+exit+"; inspect stderr.");
			check(same1 && same2, names[index]+": rejection must preserve both original input byte arrays.");
			check(!Files.exists(fresh1) && !Files.exists(fresh2), names[index]+": rejection created a fresh output.");
			check(!diagnostics.contains("Adapters Added:"), names[index]+": rejected input/output collision reached the success summary.");
		}else{
			checkSuccess(exit, work);
			final String first=record(index==6 ? "0_2_2" : "0_2_2_2_2 /1", "GG", "II");
			write(work.resolve("expected1.fq"), first);
			check(matches(output1, first), names[index]+": first output must preserve GG/II with the literal generated name.");
			if(index==5){
				final String second=record("0_2_2_2_2 /2", "TT", "II");
				write(work.resolve("expected2.fq"), second);
				check(matches(output2, second), names[index]+": mate output must preserve TT/II and pair2 name.");
			}
			check(same1 && same2, names[index]+": successful distinct output modified an input.");
		}
		aliasCases++;
	}

	/** Explicitly disables synthetic-header parsing after AddAdapters enables it in its constructor. */
	private static ArrayList<String> baseArguments(final Path input, final Path output, final String literal, final String rate){
		return new ArrayList<String>(Arrays.asList("-Xmx256m", "-Xms32m", "in="+input, "out="+output,
			"mode=write", "literal="+literal, "rate="+rate, "adderrors=f", "arc=f", "reads=1", "int=f",
			"parsecustom=f", "qin=33", "qout=33", "simd=f", "t=1", "ow=t", "app=f"));
	}

	/** Records argv and owns the returned supervisor handle; no global process search or unrelated signal is used. */
	private int run(final List<String> args, final Path work) throws Exception{
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/timeout", "--verbose",
			"--signal=TERM", "--kill-after=2s", "30s", "bash", tree.resolve("addadapters.sh").toString()));
		command.addAll(args);
		Files.write(work.resolve("command.argv"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(tree.toFile())
			.redirectOutput(work.resolve("stdout").toFile()).redirectError(work.resolve("stderr").toFile());
		final String inherited=builder.environment().get("JAVA_TOOL_OPTIONS");
		final String options=(inherited==null ? "" : inherited+" ")+"-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8";
		builder.environment().put("JAVA_TOOL_OPTIONS", options);
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
			write(work.resolve("timeout.txt"), "Owned AddAdapters child exceeded its deadline or was killed.\n");
		}
		//Callers retain post-termination input snapshots before rejecting a timeout.
		return exit;
	}

	/** GNU timeout supervises only this launch's group; allow its TERM/KILL grace before reaping the exact handle. */
	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 5)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){throw new AssertionError("Owned AddAdapters supervisor survived forced cleanup.");}
	}

	/** Repeated interruption cannot remove the cleanup deadline; restore interruption only after cleanup. */
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

	private void checkSuccess(final int exit, final Path work) throws Exception{
		noTimeout(work);
		final String diagnostics=read(work.resolve("stderr"));
		check(exit==0, "Expected successful finite generation, exit="+exit+": "+work);
		check(Files.size(work.resolve("stdout"))==0, "File-output case unexpectedly wrote stdout: "+work);
		check(diagnostics.contains("Adapters Added:") && diagnostics.contains("Valid Output:") &&
			!diagnostics.contains("Exception") && !diagnostics.contains("AssertionError"),
			"Successful generation needs its real summary without an exception: "+work);
	}
	private static void noTimeout(final Path work){
		if(Files.exists(work.resolve("timeout.txt"))){
			throw new AssertionError("A timed-out or killed child is never an accepted baseline defect: "+work);
		}
	}

	/** Retains observed bytes before making any acceptance claim, including for intentionally unsafe baseline cases. */
	private static boolean saveInput(final Path work, final Path input, final String after, final String expected) throws Exception{
		if(!Files.isRegularFile(input)){
			write(work.resolve(after+".missing"), "Generated input is no longer a regular file.\n");
			return false;
		}
		final byte[] actual=Files.readAllBytes(input);
		Files.write(work.resolve(after), actual, StandardOpenOption.CREATE_NEW);
		return Arrays.equals(actual, expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static boolean collisionDiagnostic(final String text){
		return text.contains("Input file and output file have same name") ||
			text.contains("out1 and out2 have same name") || text.contains("specified for both input and output streams");
	}
	private static boolean matches(final Path path, final String expected) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII));
	}
	private static boolean safePath(final Path path){
		final String text=path.toString();
		return text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0 && text.indexOf('\r')<0;
	}
	private static String record(final String name, final String bases, final String qualities){return "@"+name+"\n"+bases+"\n+\n"+qualities+"\n";}
	private static String read(final Path path) throws Exception{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}
	}

	/** Literal baseline and intended selections, independent of the parser implementation under test. */
	private static final class EndCase{
		EndCase(final String name_, final String before_, final String option_, final boolean baselineRight_, final boolean intendedRight_){
			name=name_; before=before_; option=option_; baselineRight=baselineRight_; intendedRight=intendedRight_;
		}
		final String name, before, option;
		final boolean baselineRight, intendedRight;
	}
	private static final EndCase[] ENDS={
		new EndCase("right-true", null, "right=t", true, true),
		new EndCase("right-false", null, "right=f", false, false),
		new EndCase("5prime-true", null, "5prime=t", false, false),
		new EndCase("5prime-false", null, "5prime=f", true, true),
		new EndCase("end-3quote", "left=t", "end=3'", true, true),
		new EndCase("end-3prime", "left=t", "end=3prime", true, true),
		new EndCase("end-3-prime", "left=t", "end=3-prime", true, true),
		new EndCase("end-right", "left=t", "end=right", true, true),
		new EndCase("end-r", "left=t", "end=r", false, true),
		new EndCase("end-5quote", "right=t", "end=5'", true, false),
		new EndCase("end-5prime", "right=t", "end=5prime", true, false),
		new EndCase("end-5-prime", "right=t", "end=5-prime", true, false),
		new EndCase("end-left", "right=t", "end=left", true, false),
		new EndCase("end-l", "right=t", "end=l", true, false)
	};
	private static final String ONE="@one\nG\n+\nI\n";
	private static final String PAIR1="@pair/1\nGG\n+\nII\n", PAIR2="@pair/2\nTT\n+\nII\n";
	private final Path tree, root;
	private final boolean baseline;
	private int checks, failures, helperCases, cliCases, aliasCases, selectorMismatches;
}
