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
 * Six AssemblyStats3 stats3.sh child-process cases for overwrite and append
 * policy forwarding. Baseline ignores explicit policy flags; fixed must still
 * create absent outputs with ow=f, reject an existing ow=f output before
 * changing it, and let app=t preserve the sentinel before a complete table.
 * @author Jean
 */
public final class AssemblyStats3OutputPolicyTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown Stats3 output-policy option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystats3outputpolicytest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final AssemblyStats3OutputPolicyTest test=new AssemblyStats3OutputPolicyTest(out, repo,
			phase.equals("baseline"));
		test.writeInput();
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==6, "Expected six Stats3 output-policy cases: "+test.cases);
		test.check(test.children==6, "Expected one fresh child per case: "+test.children);
		test.check(test.rejections==(test.baseline ? 0 : 1), "Unexpected rejection count: "+test.rejections);
		test.check(test.successfulRows==(test.baseline ? 6 : 5),
			"Unexpected successful data rows: "+test.successfulRows);
		final String summary="ASSEMBLY_STATS3_OUTPUT_POLICY_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" rejections="+test.rejections+" successfulRows="+test.successfulRows+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private AssemblyStats3OutputPolicyTest(final Path out_, final Path repo_,
			final boolean baseline_) throws IOException{
		out=out_; repo=repo_; baseline=baseline_; input=out.resolve("tiny.fa");
		launcher=repo.resolve("stats3.sh");
		require(Files.isRegularFile(launcher), "Missing stats3.sh at "+launcher);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\texpectedReject\texpectedText"+
			"\tstdoutBytes\tstderrBytes\toutputBytes\toutputMatches\tstdoutEmpty"+
			"\tstderrMatches\tinputsUnchanged\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes the retained four-A FASTA exactly once. */
	private void writeInput() throws IOException{
		final byte[] data=(">tiny\nAAAA\n").getBytes(StandardCharsets.US_ASCII);
		Files.write(input, data, StandardOpenOption.CREATE_NEW);
		original.put(input, data);
	}

	/** Launches stats3.sh once and checks the phase-specific output-policy oracle. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final Path output=caseDir.resolve("stats3.tsv");
		if(row.existing){Files.write(output, SENTINEL.getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);}

		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("simd=f");
		command.add("t=1");
		command.add("in="+input);
		command.add("out="+output);
		for(final String arg : row.extraArgs){command.add(arg);}
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

		final boolean reject=row.rejectsWhenFixed && !baseline;
		final String table=expectedTable(input);
		final String expected=row.appendWhenFixed && !baseline ? SENTINEL+table :
			(reject ? SENTINEL : table);
		final byte[] outputBytes=Files.isRegularFile(output) ? Files.readAllBytes(output) : new byte[0];
		final boolean outputMatches=expected.equals(new String(outputBytes, StandardCharsets.UTF_8));
		final boolean stdoutEmpty=stdout.length==0;
		final boolean stderrMatches=stderrMatches(errText, output, reject);
		final boolean unchanged=inputsUnchanged();
		if(reject){rejections++;}else{successfulRows++;}

		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t"+reject+"\t"+clean(expected)+"\t"+stdout.length+"\t"+stderr.length+
			"\t"+outputBytes.length+"\t"+outputMatches+"\t"+stdoutEmpty+"\t"+stderrMatches+
			"\t"+unchanged+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(reject ? exit!=0 : exit==0, row.name+": child exit was "+exit);
		check(Files.isRegularFile(output), row.name+": expected output file missing: "+output);
		check(outputMatches, row.name+": output text mismatch at "+output);
		check(stdoutEmpty, row.name+": stats3 wrote primary table to stdout despite out=");
		check(stderrMatches, row.name+": stderr did not match output policy oracle: "+errText);
		check(command.contains("-ea") && command.contains("simd=f") && command.contains("t=1"),
			row.name+": child argv lost assertions, simd=f or t=1");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Checks normal and ow=f rejection stderr without probing OS-specific failure text. */
	private static boolean stderrMatches(final String errText, final Path output, final boolean reject){
		final boolean common=errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
			!errText.contains("Unknown parameter") && !errText.contains("Error writing assembly statistics");
		if(reject){
			return common && errText.contains("java.lang.AssertionError") &&
				errText.contains(output.toString());
		}
		return common && !errText.contains("Exception") && !errText.contains("AssertionError");
	}

	/** Builds the exact Stats3 table for the retained all-A four-base input. */
	private static String expectedTable(final Path input){
		return "fname\tsize\tcontigs\tgc\tmaxContig\t5kplus\t10kplus\t25kplus\t50kplus\n"+
			input+"\t4\t1\t0\t4\t0\t0\t0\t0\n";
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

	private static String clean(final String value){
		return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
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

	/** One fresh stats3.sh CLI case and its fixed-phase output policy. */
	private static final class Row{
		Row(final String name_, final boolean existing_, final boolean appendWhenFixed_,
				final boolean rejectsWhenFixed_, final String... extraArgs_){
			name=name_;
			existing=existing_;
			appendWhenFixed=appendWhenFixed_;
			rejectsWhenFixed=rejectsWhenFixed_;
			extraArgs=extraArgs_;
		}
		final String name;
		final boolean existing, appendWhenFixed, rejectsWhenFixed;
		final String[] extraArgs;
	}

	private final Path out, repo, input, launcher;
	private final Map<Path, byte[]> original=new HashMap<Path, byte[]>();
	private int cases, children, rejections, successfulRows, checks, failures;
	private final boolean baseline;

	private static final String SENTINEL="SENTINEL\n";
	private static final Row[] ROWS={
		new Row("existing-default", true, false, false),
		new Row("existing-overwrite-true", true, false, false, "ow=t"),
		new Row("absent-overwrite-false", false, false, false, "ow=f"),
		new Row("existing-overwrite-false", true, false, true, "ow=f"),
		new Row("existing-append", true, true, false, "app=t"),
		new Row("existing-append-overwrite-false", true, true, false, "app=t", "ow=f")
	};
}
