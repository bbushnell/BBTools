package test.jgi;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import parse.Parse;

/**
 * Seven AssemblyStats3 child-process cases for Assembly FASTA blank-line
 * handling. Baseline rejects zero-length physical lines; fixed must skip them
 * without changing any nonempty sequence line, including a one-space line.
 * @author Jean
 */
public final class AssemblyBlankLinesTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown Assembly blank-line option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblyblanklinestest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final AssemblyBlankLinesTest test=new AssemblyBlankLinesTest(out, repo, phase.equals("baseline"));
		test.writeInputs();
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==7, "Expected seven Assembly blank-line cases: "+test.cases);
		test.check(test.children==7, "Expected one fresh child per case: "+test.children);
		test.check(test.childFailures==(test.baseline ? 5 : 0),
			"Unexpected child failure count: "+test.childFailures);
		test.check(test.successfulRows==(test.baseline ? 2 : 7),
			"Unexpected successful data rows: "+test.successfulRows);
		final String summary="ASSEMBLY_BLANK_LINES_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" childFailures="+test.childFailures+" successfulRows="+test.successfulRows+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private AssemblyBlankLinesTest(final Path out_, final Path repo_, final boolean baseline_) throws IOException{
		out=out_; repo=repo_; baseline=baseline_; inputs=out.resolve("inputs");
		launcher=repo.resolve("stats3.sh");
		require(Files.isRegularFile(launcher), "Missing stats3.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\texpectedFailure\toutput"+
			"\tstdoutBytes\tstderrBytes\tdataRows\toutputMatches\tstdoutEmpty"+
			"\tstderrMatches\tinputsUnchanged\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes every retained literal FASTA exactly once. */
	private void writeInputs() throws IOException{
		for(final Row row : ROWS){write(row.name+".fa", row.fasta);}
	}

	/** Writes one literal FASTA file and records its original bytes. */
	private void write(final String name, final String text) throws IOException{
		final Path path=inputs.resolve(name);
		final byte[] data=text.getBytes(StandardCharsets.US_ASCII);
		Files.write(path, data, StandardOpenOption.CREATE_NEW);
		original.put(path, data);
	}

	/** Launches stats3.sh once and checks the phase-specific blank-line oracle. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final Path input=inputs.resolve(row.name+".fa");
		final Path output=caseDir.resolve(row.name+".tsv");
		final boolean expectFailure=baseline && row.blankLine;
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("simd=f");
		command.add("t=1");
		command.add("in="+input);
		command.add("out="+output);
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
		Files.write(caseDir.resolve("exit.txt"), (Integer.toString(exit)+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		final boolean outputPresent=Files.isRegularFile(output);
		final byte[] outputBytes=outputPresent ? Files.readAllBytes(output) : new byte[0];
		final int dataRows=expectFailure ? 0 : 1;
		final String expected=expectFailure ? HEADER : expectedText(input, row);
		final boolean outputMatches=outputPresent && expected.equals(new String(outputBytes, StandardCharsets.UTF_8));
		final boolean stdoutEmpty=stdout.length==0;
		final boolean stderrMatches=errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
			(expectFailure ? errText.contains("java.lang.ArrayIndexOutOfBoundsException") &&
				errText.contains("jgi.Assembly.load") :
				!errText.contains("Exception") && !errText.contains("AssertionError") &&
				!errText.contains("Error reading assembly input"));
		final boolean unchanged=inputsUnchanged();
		if(exit!=0){childFailures++;}
		successfulRows+=dataRows;
		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t"+expectFailure+"\t"+output+"\t"+stdout.length+"\t"+
			stderr.length+"\t"+dataRows+"\t"+outputMatches+"\t"+stdoutEmpty+"\t"+
			stderrMatches+"\t"+unchanged+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);
		check(expectFailure ? exit!=0 : exit==0, row.name+": child exit was "+exit);
		check(outputPresent, row.name+": expected output missing: "+output);
		check(outputMatches, row.name+": output text mismatch at "+output);
		check(stdoutEmpty, row.name+": stats3 wrote primary table to stdout despite out=");
		check(stderrMatches, row.name+": stderr did not match the phase oracle: "+errText);
		check(command.contains("-ea") && command.contains("simd=f") && command.contains("t=1"),
			row.name+": child argv lost assertions, simd=f or t=1");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Builds the exact Stats3 table for one literal input. */
	private static String expectedText(final Path input, final Row row){
		return HEADER+input+"\t"+row.size+"\t"+row.contigs+"\t0\t"+row.maxContig+"\t0\t0\t0\t0\n";
	}

	/** Checks every retained input against the literal bytes written before children ran. */
	private boolean inputsUnchanged() throws IOException{
		for(final Map.Entry<Path, byte[]> entry : original.entrySet()){
			if(!Arrays.equals(entry.getValue(), Files.readAllBytes(entry.getKey()))){return false;}
		}
		return true;
	}

	/** Joins argv for stable evidence. */
	private static String join(final List<String> command){
		final StringBuilder sb=new StringBuilder();
		for(final String s : command){sb.append(s).append('\n');}
		return sb.toString();
	}

	/** Accumulates mismatches so every child leaves stdout, stderr and exit evidence. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}

	/** Throws on fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** One fresh Stats3 CLI case and its fixed output metrics. */
	private static final class Row{
		Row(final String name_, final String fasta_, final boolean blankLine_, final int size_,
				final int contigs_, final int maxContig_){
			name=name_; fasta=fasta_; blankLine=blankLine_;
			size=size_; contigs=contigs_; maxContig=maxContig_;
		}

		final String name, fasta;
		final boolean blankLine;
		final int size, contigs, maxContig;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new HashMap<Path, byte[]>();
	private int cases, children, childFailures, successfulRows, checks, failures;
	private final boolean baseline;

	private static final String HEADER=
		"fname\tsize\tcontigs\tgc\tmaxContig\t5kplus\t10kplus\t25kplus\t50kplus\n";

	private static final Row[] ROWS={
		new Row("no-blanks", ">a\nAA\n>b\nAAA\n", false, 5, 2, 3),
		new Row("leading-blank", "\n>a\nAA\n", true, 2, 1, 2),
		new Row("after-header", ">a\n\nAA\n", true, 2, 1, 2),
		new Row("inside-sequence", ">a\nAA\n\nAAA\n", true, 5, 1, 5),
		new Row("between-records", ">a\nAA\n\n>b\nAAA\n", true, 5, 2, 3),
		new Row("trailing-blank", ">a\nAA\n\n", true, 2, 1, 2),
		new Row("space-line", ">a\nAA\n \nA\n", false, 4, 1, 4)
	};
}
