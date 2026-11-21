package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import parse.Parse;

/**
 * Finite BBMask empty-reference controls for repeat and entropy worker queues.
 * @author Jean
 */
public final class BBMaskEmptyReferenceTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask empty-reference option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmaskemptyreferencetest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final BBMaskEmptyReferenceTest test=new BBMaskEmptyReferenceTest(out, repo,
			phase.equals("baseline"));
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==6, "Expected six BBMask empty-reference cases: "+test.cases);
		test.check(test.children==6, "Expected six fresh bbmask.sh children: "+test.children);
		test.check(test.changed==3, "Expected three repaired rows: "+test.changed);
		test.check(test.controls==3, "Expected three unchanged controls: "+test.controls);
		final int expectedSourceMatches=(test.baseline ? 6 : 3);
		final int expectedCorrectedMatches=(test.baseline ? 3 : 6);
		test.check(test.sourceMatches==expectedSourceMatches,
			"Expected "+expectedSourceMatches+" source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==expectedCorrectedMatches,
			"Expected "+expectedCorrectedMatches+" corrected matches: "+test.correctedMatches);
		final String summary="BBMASK_EMPTY_REFERENCE_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" changed="+test.changed+" controls="+test.controls+
			" sourceMatches="+test.sourceMatches+
			" correctedMatches="+test.correctedMatches+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private BBMaskEmptyReferenceTest(final Path out_, final Path repo_, final boolean baseline_)
			throws IOException{
		out=out_; repo=repo_; baseline=baseline_;
		launcher=repo.resolve("bbmask.sh");
		require(Files.isRegularFile(launcher), "Missing bbmask.sh at "+launcher);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\tobservedCounter"+
			"\tobservedDenominator\toutputExists\toutputBytes\tsourceMatches"+
			"\tcorrectedMatches\tstderrMatches\tstdoutEmpty\tinputsUnchanged"+
			"\tmatchPhase\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Launches one bbmask.sh child and checks exact output or exact queue failure. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		if(row.source.changed(row.corrected)){changed++;}else{controls++;}
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final Path ref=caseDir.resolve("ref.fa");
		write(ref, row.inputBytes());
		final Path output=caseDir.resolve("out.fa");
		check(!Files.exists(output), row.name+": output path exists before child: "+output);
		final ArrayList<String> command=baseCommand(row, ref, output);
		Files.write(caseDir.resolve("expected.tsv"), expected(row).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		Files.write(caseDir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);

		final Path stdoutPath=caseDir.resolve("stdout.txt");
		final Path stderrPath=caseDir.resolve("stderr.txt");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdoutPath.toFile()).redirectError(stderrPath.toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		final byte[] stdout=Files.readAllBytes(stdoutPath);
		final byte[] stderr=Files.readAllBytes(stderrPath);
		final String errText=new String(stderr, StandardCharsets.UTF_8);
		Files.write(caseDir.resolve("exit.txt"), (Integer.toString(exit)+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);

		final boolean outputExists=Files.isRegularFile(output);
		final byte[] outputBytes=outputExists ? Files.readAllBytes(output) : null;
		final Counter counter=parseCounter(errText);
		final boolean sourceMatch=row.source.matches(exit, outputExists, outputBytes, counter, errText);
		final boolean correctedMatch=row.corrected.matches(exit, outputExists, outputBytes, counter, errText);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean stderrMatches=(baseline ? row.source : row.corrected).stderrMatches(errText);
		final boolean stdoutEmpty=stdout.length==0;
		final boolean unchanged=inputsUnchanged();
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;

		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t"+counter.numeratorText()+"\t"+counter.denominatorText()+
			"\t"+outputExists+"\t"+(outputBytes==null ? "" : outputBytes.length)+
			"\t"+sourceMatch+"\t"+correctedMatch+"\t"+stderrMatches+
			"\t"+stdoutEmpty+"\t"+unchanged+"\t"+phaseMatch+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);

		check(phaseMatch, row.name+": output/counter/exception mismatch: "+errText);
		check(stderrMatches, row.name+": stderr did not match oracle: "+errText);
		check(stdoutEmpty, row.name+": bbmask.sh wrote to stdout");
		check(unchanged, row.name+": a child modified literal input bytes");
		check(command.contains("-ea") && command.contains("t=1") &&
			command.contains("simd=f") && command.contains("split=f") &&
			command.contains("kr=1") && command.contains("ke=1") &&
			command.contains("window=4") && command.contains("minlen=4") &&
			command.contains("mincount=2") && command.contains("entropy=0.5") &&
			command.contains("ow=t") && command.contains("mr="+row.mr) &&
			command.contains("me="+row.me) && command.contains("entropymode="+row.entropyMode) &&
			command.contains("minscaf="+row.minscaf),
			row.name+": child argv lost a frozen CLI pin");
	}

	/** Launch arguments with all matrix parameters explicit. */
	private ArrayList<String> baseCommand(final Row row, final Path ref, final Path output){
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("t=1");
		command.add("simd=f");
		command.add("split=f");
		command.add("kr=1");
		command.add("ke=1");
		command.add("window=4");
		command.add("minlen=4");
		command.add("mincount=2");
		command.add("entropy=0.5");
		command.add("ow=t");
		command.add("mr="+row.mr);
		command.add("me="+row.me);
		command.add("entropymode="+row.entropyMode);
		command.add("minscaf="+row.minscaf);
		command.add("in="+ref);
		command.add("out="+output);
		return command;
	}

	/** Writes a literal input and records it for after-child preservation checks. */
	private void write(final Path path, final byte[] bytes) throws IOException{
		Files.write(path, bytes, StandardOpenOption.CREATE_NEW);
		original.put(path, bytes);
	}

	/** Expected source and corrected row for one case. */
	private static String expected(final Row row){
		return "phase\tresult\tcounter\tdenominator\tmethod\toutputBytes\n"+
			"source\t"+row.source+"\n"+
			"corrected\t"+row.corrected+"\n";
	}

	/** Extracts BBMask's final masked counter, if it reached the final report. */
	private static Counter parseCounter(final String stderr){
		final Matcher matcher=COUNTER.matcher(stderr);
		if(!matcher.find()){return Counter.absent();}
		return new Counter(Long.parseLong(matcher.group(1)), Long.parseLong(matcher.group(2)));
	}

	/** Requires every literal input path to retain its original bytes. */
	private boolean inputsUnchanged() throws IOException{
		for(final Path path : original.keySet()){
			if(!Arrays.equals(original.get(path), Files.readAllBytes(path))){return false;}
		}
		return true;
	}

	/** Records one assertion without hiding later evidence after a failure. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	/** Hard construction-time requirement. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Shell-escaped argv transcript for the evidence directory. */
	private static String join(final List<String> command){
		final StringBuilder builder=new StringBuilder();
		for(final String arg : command){
			if(builder.length()>0){builder.append(' ');}
			builder.append('\'').append(arg.replace("'", "'\\''")).append('\'');
		}
		return builder.append('\n').toString();
	}

	/** Exact expectation for one phase. */
	private static final class Expect{
		Expect(final long numerator_, final long denominator_, final byte[] outputBytes_){
			numerator=numerator_;
			denominator=denominator_;
			outputBytes=outputBytes_;
			failureMethod=null;
		}

		Expect(final String failureMethod_){
			numerator=-1;
			denominator=-1;
			outputBytes=null;
			failureMethod=failureMethod_;
		}

		boolean changed(final Expect other){
			if(failureMethod==null ? other.failureMethod!=null : !failureMethod.equals(other.failureMethod)){
				return true;
			}
			return numerator!=other.numerator || denominator!=other.denominator ||
				!Arrays.equals(outputBytes, other.outputBytes);
		}

		boolean matches(final int exit, final boolean outputExists, final byte[] observed,
				final Counter counter, final String stderr){
			if(failureMethod!=null){
				return exit==1 && !outputExists && counter.absent && stderrMatches(stderr);
			}
			return exit==0 && outputExists && Arrays.equals(outputBytes, observed) &&
				counter.numerator==numerator && counter.denominator==denominator &&
				stderrMatches(stderr);
		}

		boolean stderrMatches(final String stderr){
			if(failureMethod!=null){
				return stderr.contains("Exception in thread \"main\" java.lang.IllegalArgumentException") &&
					stderr.contains("java.util.concurrent.ArrayBlockingQueue.<init>") &&
					stderr.contains("jgi.BBMask."+failureMethod);
			}
			return !stderr.contains("Exception in thread") && !stderr.contains("java.lang.") &&
				!stderr.contains("ERROR") && !stderr.contains("Error");
		}

		@Override
		public String toString(){
			if(failureMethod!=null){return "failure\t\t\t"+failureMethod+"\t";}
			return "success\t"+numerator+"\t"+denominator+"\t\t"+outputBytes.length;
		}

		final long numerator;
		final long denominator;
		final byte[] outputBytes;
		final String failureMethod;
	}

	/** Final counter parsed from stderr. */
	private static final class Counter{
		Counter(final long numerator_, final long denominator_){
			numerator=numerator_;
			denominator=denominator_;
			absent=false;
		}

		static Counter absent(){return new Counter();}

		private Counter(){
			numerator=-1;
			denominator=-1;
			absent=true;
		}

		String numeratorText(){return absent ? "" : Long.toString(numerator);}
		String denominatorText(){return absent ? "" : Long.toString(denominator);}

		final long numerator;
		final long denominator;
		final boolean absent;
	}

	/** One direct child row. */
	private static final class Row{
		Row(final String name_, final byte[] input_, final String mr_, final String me_,
				final String entropyMode_, final int minscaf_, final Expect source_,
				final Expect corrected_){
			name=name_;
			input=input_;
			mr=mr_;
			me=me_;
			entropyMode=entropyMode_;
			minscaf=minscaf_;
			source=source_;
			corrected=corrected_;
		}

		byte[] inputBytes(){return input;}

		final String name;
		final byte[] input;
		final String mr;
		final String me;
		final String entropyMode;
		final int minscaf;
		final Expect source;
		final Expect corrected;
	}

	/*--------------------------------------------------------------*/

	/** Evidence output directory. */
	private final Path out;
	/** BBTools repository root containing bbmask.sh. */
	private final Path repo;
	/** Absolute path to bbmask.sh. */
	private final Path launcher;
	/** True when checking the uncorrected source phase. */
	private final boolean baseline;
	/** Original literal input bytes, keyed by path. */
	private final LinkedHashMap<Path, byte[]> original=new LinkedHashMap<Path, byte[]>();

	/** Cases run. */
	private int cases=0;
	/** Fresh bbmask.sh children run. */
	private int children=0;
	/** Rows expected to change after the repair. */
	private int changed=0;
	/** Rows expected to match in both phases. */
	private int controls=0;
	/** Rows matching the source oracle. */
	private int sourceMatches=0;
	/** Rows matching the corrected oracle. */
	private int correctedMatches=0;
	/** Checks evaluated. */
	private int checks=0;
	/** Checks failed. */
	private int failures=0;

	private static final byte[] EMPTY=new byte[0];
	private static final byte[] ACGT=">ref\nACGT\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] A8=">ref\nAAAAAAAA\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] N8=">ref\nNNNNNNNN\n".getBytes(StandardCharsets.US_ASCII);
	private static final Expect EMPTY_OK=new Expect(0, 0, EMPTY);
	private static final Expect N8_OK=new Expect(8, 8, N8);
	private static final Pattern COUNTER=Pattern.compile("Total Bases Masked:\\s*(\\d+)/(\\d+)");
	private static final Row[] ROWS={
		new Row("empty-repeat", EMPTY, "t", "f", "t", 1,
			new Expect("maskRepeats"), EMPTY_OK),
		new Row("empty-entropy", EMPTY, "f", "t", "t", 1,
			new Expect("maskLowEntropy"), EMPTY_OK),
		new Row("filtered-entropy", ACGT, "f", "t", "t", 5,
			new Expect("maskLowEntropy"), EMPTY_OK),
		new Row("empty-low-complexity-control", EMPTY, "f", "t", "f", 1,
			EMPTY_OK, EMPTY_OK),
		new Row("nonempty-repeat-control", A8, "t", "f", "t", 1,
			N8_OK, N8_OK),
		new Row("nonempty-entropy-control", A8, "f", "t", "t", 1,
			N8_OK, N8_OK)
	};
}
