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
 * Seven AssemblyStats3 stats3.sh child-process cases for assembly input
 * aliases. Baseline accepts input/input1/in1 through Parser but never copies
 * them into the local assembly list; fixed should make those aliases behave
 * like in= while preserving bare-file and existing-file-before-comma handling.
 * @author Jean
 */
public final class AssemblyStats3InputAliasTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown Stats3 input-alias option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystats3inputaliastest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final AssemblyStats3InputAliasTest test=new AssemblyStats3InputAliasTest(out, repo,
			phase.equals("baseline"));
		test.writeInputs();
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==7, "Expected seven Stats3 input-alias cases: "+test.cases);
		test.check(test.children==7, "Expected one fresh child per case: "+test.children);
		test.check(test.successfulRows==(test.baseline ? 5 : 9),
			"Unexpected successful data rows: "+test.successfulRows);
		final String summary="ASSEMBLY_STATS3_INPUT_ALIAS_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" successfulRows="+test.successfulRows+" checks="+test.checks+
			" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private AssemblyStats3InputAliasTest(final Path out_, final Path repo_,
			final boolean baseline_) throws IOException{
		out=out_; repo=repo_; baseline=baseline_; inputs=out.resolve("inputs");
		launcher=repo.resolve("stats3.sh");
		require(Files.isRegularFile(launcher), "Missing stats3.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\tstdoutBytes\tstderrBytes"+
			"\trows\tfilenames\toutputMatches\tstdoutEmpty\tstderrClean\tinputsUnchanged\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes retained single-contig FASTA files with literal, distinct lengths. */
	private void writeInputs() throws IOException{
		write("a.fa", 4);
		write("b.fa", 6);
		write("literal,comma.fa", 8);
	}

	/** Writes one single-contig all-A FASTA whose integral Stats3 GC prints as literal 0. */
	private void write(final String name, final int length) throws IOException{
		final char[] bases=new char[length];
		Arrays.fill(bases, 'A');
		final byte[] data=(">"+name+"\n"+new String(bases)+"\n").getBytes(StandardCharsets.US_ASCII);
		final Path path=inputs.resolve(name);
		Files.write(path, data, StandardOpenOption.CREATE_NEW);
		original.put(path, data);
		literalLengths.put(path.toString(), length);
	}

	/** Launches stats3.sh once and checks the phase-specific alias oracle. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final Path output=caseDir.resolve(row.name+".tsv");
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("simd=f");
		command.add("t=1");
		for(final String arg : row.args(inputs, output)){command.add(arg);}
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
		final List<String> expected=row.expectedInputs(inputs, baseline);
		final String expectedText=expectedText(expected);
		final boolean outputMatches=outputPresent && expectedText.equals(new String(outputBytes, StandardCharsets.UTF_8));
		final boolean stdoutEmpty=stdout.length==0;
		final boolean stderrClean=errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
			!errText.contains("Unknown parameter") && !errText.contains("Error reading assembly input") &&
			!errText.contains("Exception") && !errText.contains("AssertionError");
		final boolean unchanged=inputsUnchanged();
		successfulRows+=expected.size();
		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t"+stdout.length+"\t"+stderr.length+"\t"+expected.size()+
			"\t"+expected+"\t"+outputMatches+"\t"+stdoutEmpty+"\t"+stderrClean+
			"\t"+unchanged+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(exit==0, row.name+": child exit was "+exit);
		check(outputPresent, row.name+": expected output missing: "+output);
		check(outputMatches, row.name+": output text mismatch at "+output);
		check(stdoutEmpty, row.name+": stats3 wrote primary table to stdout despite out=");
		check(stderrClean, row.name+": stderr contains an unexpected diagnostic: "+errText);
		check(command.contains("-ea") && command.contains("simd=f") && command.contains("t=1"),
			row.name+": child argv lost assertions, simd=f or t=1");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Builds the exact Stats3 table for these all-A single-contig inputs. */
	private String expectedText(final List<String> names){
		final StringBuilder sb=new StringBuilder("fname\tsize\tcontigs\tgc\tmaxContig\t5kplus\t10kplus\t25kplus\t50kplus\n");
		for(final String name : names){
			final Integer length=literalLengths.get(name);
			require(length!=null, "No literal length registered for "+name);
			sb.append(name).append('\t').append(length.intValue()).append("\t1\t0\t")
				.append(length.intValue()).append("\t0\t0\t0\t0\n");
		}
		return sb.toString();
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

	/** One child-process case and its phase-dependent Stats3 alias oracle. */
	private static final class Row{
		Row(final String name_, final String[] args_, final String[] baseline_, final String[] fixed_){
			name=name_; args=args_; baseline=baseline_; fixed=fixed_;
		}
		String[] args(final Path inputs, final Path out){
			final String[] expanded=new String[args.length];
			for(int i=0; i<args.length; i++){expanded[i]=expand(args[i], inputs, out);}
			return expanded;
		}
		List<String> expectedInputs(final Path inputs, final boolean baselinePhase){
			final String[] source=baselinePhase ? baseline : fixed;
			final ArrayList<String> list=new ArrayList<String>();
			for(final String name : source){list.add(expand(name, inputs, null));}
			return list;
		}
		private static String expand(final String value, final Path inputs, final Path out){
			return value.replace("${A}", inputs.resolve("a.fa").toString())
				.replace("${B}", inputs.resolve("b.fa").toString())
				.replace("${COMMA}", inputs.resolve("literal,comma.fa").toString())
				.replace("${OUT}", out==null ? "" : out.toString());
		}
		final String name;
		final String[] args, baseline, fixed;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new HashMap<Path, byte[]>();
	private final Map<String, Integer> literalLengths=new HashMap<String, Integer>();
	private int cases, children, successfulRows, checks, failures;
	private final boolean baseline;

	private static final Row[] ROWS={
		new Row("in-a", new String[]{"in=${A}", "out=${OUT}"},
			new String[]{"${A}"}, new String[]{"${A}"}),
		new Row("bare-a", new String[]{"${A}", "out=${OUT}"},
			new String[]{"${A}"}, new String[]{"${A}"}),
		new Row("whole-comma", new String[]{"in=${COMMA}", "out=${OUT}"},
			new String[]{"${COMMA}"}, new String[]{"${COMMA}"}),
		new Row("comma-list", new String[]{"in=${A},${B}", "out=${OUT}"},
			new String[]{"${A}", "${B}"}, new String[]{"${A}", "${B}"}),
		new Row("input-a", new String[]{"input=${A}", "out=${OUT}"},
			new String[]{}, new String[]{"${A}"}),
		new Row("input1-b", new String[]{"input1=${B}", "out=${OUT}"},
			new String[]{}, new String[]{"${B}"}),
		new Row("in1-list", new String[]{"in1=${A},${B}", "out=${OUT}"},
			new String[]{}, new String[]{"${A}", "${B}"})
	};
}
