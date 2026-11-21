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
 * Eleven AssemblyStatsWrapper child-process cases for values containing the
 * first '=' in a path.  Each case enters through statswrapper.sh so the public
 * wrapper expansion and AssemblyStats2 preflight decide the outcome.  Baseline
 * must fail on four missing prefixes created by split("="); fixed must retain
 * complete in/ref values.  Existing comma-list expansion, whole existing paths,
 * and empty in=/ref= no-ops are retained controls.
 * @author Jean
 */
public final class AssemblyStatsWrapperFirstEqualsTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown Wrapper fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystatswrapperfirstequalstest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final AssemblyStatsWrapperFirstEqualsTest test=new AssemblyStatsWrapperFirstEqualsTest(out, repo,
				phase.equals("baseline"));
		test.writeInputs();
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==11, "Expected eleven Wrapper first-equals cases: "+test.cases);
		test.check(test.children==11, "Expected one fresh child per case: "+test.children);
		test.check(test.baselineFailures==(test.baseline ? 4 : 0),
			"Expected four baseline missing-prefix failures and zero fixed failures: "+test.baselineFailures);
		test.check(test.successfulRows==(test.baseline ? 8 : 13),
			"Unexpected successful data row count: "+test.successfulRows);
		test.check(test.checks==81, "Eleven cases with seven checks plus four summary checks: "+test.checks);
		final String summary="ASSEMBLY_WRAPPER_FIRST_EQUALS_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" baselineFailures="+test.baselineFailures+" successfulRows="+test.successfulRows+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private AssemblyStatsWrapperFirstEqualsTest(final Path out_, final Path repo_, final boolean baseline_) throws IOException{
		out=out_; repo=repo_; baseline=baseline_; inputs=out.resolve("inputs");
		launcher=repo.resolve("statswrapper.sh");
		require(Files.isRegularFile(launcher), "Missing statswrapper.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\texpectedExit\trows\texpectedRows"+
			"\tbaselineFailure\tfailureIndex\theaders\tfilenames\tlengths\tmissingPrefix\tstdoutBytes"+
			"\tstderrBytes\tinputsUnchanged\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes five retained FASTA files with literal, distinct scaffold lengths. */
	private void writeInputs() throws IOException{
		write("plain-alpha.fa", 4);
		write("plain-beta.fa", 6);
		write("eq=one.fa", 8);
		write("rhs=two.fa", 10);
		write("literal,comma.fa", 5);
	}

	/** Writes one single-contig FASTA whose scaffold and contig length equal length. */
	private void write(final String name, final int length) throws IOException{
		final char[] bases=new char[length];
		Arrays.fill(bases, 'A');
		final byte[] data=(">"+name+"\n"+new String(bases)+"\n").getBytes(StandardCharsets.US_ASCII);
		final Path path=inputs.resolve(name);
		Files.write(path, data, StandardOpenOption.CREATE_NEW);
		original.put(path, data);
		literalLengths.put(canonical(path), length);
	}

	/** Launches statswrapper.sh once and checks the phase-specific oracle. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("simd=f");
		command.add("t=1");
		for(final String arg : row.args(inputs)){command.add(arg);}
		Files.write(caseDir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
				StandardOpenOption.CREATE_NEW);
		final Path stdoutPath=caseDir.resolve("stdout.txt");
		final Path stderrPath=caseDir.resolve("stderr.txt");
		final Process process=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdoutPath.toFile()).redirectError(stderrPath.toFile()).start();
		final int exit=process.waitFor();
		final byte[] stdout=Files.readAllBytes(stdoutPath);
		final byte[] stderr=Files.readAllBytes(stderrPath);
		final String outText=new String(stdout, StandardCharsets.UTF_8);
		final String errText=new String(stderr, StandardCharsets.UTF_8);
		Files.write(caseDir.resolve("exit.txt"), (Integer.toString(exit)+"\n").getBytes(StandardCharsets.UTF_8),
				StandardOpenOption.CREATE_NEW);
		final boolean expectFailure=baseline && row.missingPrefix!=null;
		if(expectFailure){baselineFailures++;}
		final List<String> expected=row.expectedFilenames(inputs, expectFailure);
		final ParsedOutput parsed=parseOutput(outText);
		successfulRows+=parsed.filenames.size();
		final boolean unchanged=inputsUnchanged();
		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t"+(expectFailure ? 1 : 0)+"\t"+parsed.filenames.size()+"\t"+expected.size()+
			"\t"+expectFailure+"\t"+row.failureIndex+"\t"+parsed.headerCount+"\t"+parsed.filenames+
			"\t"+parsed.lengths+"\t"+row.expectedMissingPrefix(inputs)+"\t"+stdout.length+
			"\t"+stderr.length+"\t"+unchanged+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(exit==(expectFailure ? 1 : 0), row.name+": exit expected "+(expectFailure ? 1 : 0)+", saw "+exit);
		check(parsed.headerCount==(expected.isEmpty() ? 0 : 1), row.name+": header count expected "+
			(expected.isEmpty() ? 0 : 1)+", saw "+parsed.headerCount);
		check(parsed.filenames.equals(expected), row.name+": filenames expected "+expected+", saw "+parsed.filenames);
		check(lengthsMatch(parsed.filenames, parsed.lengths), row.name+": lengths mismatched filenames "+
			parsed.filenames+" / "+parsed.lengths);
		if(row.missingPrefix!=null){
			final String diagnostic="Input file does not appear to be valid: "+row.expectedMissingPrefix(inputs);
			check(!expectFailure || errText.contains(diagnostic), row.name+": missing exact diagnostic "+diagnostic);
		}else{
			check(!errText.contains("Input file does not appear to be valid"), row.name+": unexpected preflight failure");
		}
		check(!expected.isEmpty() || outText.isEmpty(), row.name+": empty-value control wrote stdout");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Parses the format=3 table written by statswrapper.sh. */
	private static ParsedOutput parseOutput(final String text){
		final ArrayList<String> filenames=new ArrayList<String>();
		final ArrayList<Integer> lengths=new ArrayList<Integer>();
		int headerCount=0, scafBp=-1, contigBp=-1, filename=-1;
		for(final String line : text.split("\\R")){
			if(line.isEmpty()){continue;}
			final String[] fields=line.split("\t", -1);
			if(fields[0].equals("n_scaffolds")){
				headerCount++;
				scafBp=index(fields, "scaf_bp");
				contigBp=index(fields, "contig_bp");
				filename=index(fields, "filename");
			}else{
				require(scafBp>=0 && contigBp>=0 && filename>=0, "Data row arrived before a complete header: "+line);
				final int scaffold=Integer.parseInt(fields[scafBp]);
				final int contig=Integer.parseInt(fields[contigBp]);
				require(scaffold==contig, "Literal inputs have no gaps: "+line);
				lengths.add(scaffold);
				filenames.add(fields[filename]);
			}
		}
		return new ParsedOutput(filenames, lengths, headerCount);
	}

	/** Finds one exact format=3 column. */
	private static int index(final String[] fields, final String name){
		for(int i=0; i<fields.length; i++){
			if(fields[i].equals(name)){return i;}
		}
		return -1;
	}

	/** Checks every retained input against the literal bytes written before children ran. */
	private boolean inputsUnchanged() throws IOException{
		for(final Map.Entry<Path, byte[]> entry : original.entrySet()){
			if(!Arrays.equals(entry.getValue(), Files.readAllBytes(entry.getKey()))){return false;}
		}
		return true;
	}

	/** Verifies each output row's literal scaf_bp against the row filename. */
	private boolean lengthsMatch(final List<String> filenames, final List<Integer> lengths){
		if(filenames.size()!=lengths.size()){return false;}
		for(int i=0; i<filenames.size(); i++){
			final Integer expected=literalLengths.get(filenames.get(i));
			if(expected==null || expected.intValue()!=lengths.get(i).intValue()){return false;}
		}
		return true;
	}

	/** Returns the same canonical spelling AssemblyStats2 appends as the filename column. */
	private static String canonical(final Path path) throws IOException{return path.toFile().getCanonicalPath();}

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

	/** One child-process case and its phase-dependent order oracle. */
	private static final class Row{
		Row(final String name_, final int failureIndex_, final String missingPrefix_,
				final String[] args_, final String... expected_){
			name=name_; failureIndex=failureIndex_; missingPrefix=missingPrefix_; args=args_; expected=expected_;
		}
		String[] args(final Path inputs){
			final String[] expanded=new String[args.length];
			for(int i=0; i<args.length; i++){expanded[i]=expand(args[i], inputs);}
			return expanded;
		}
		List<String> expectedFilenames(final Path inputs, final boolean expectFailure) throws IOException{
			final ArrayList<String> list=new ArrayList<String>();
			final int limit=expectFailure ? failureIndex : expected.length;
			for(int i=0; i<limit; i++){
				list.add(canonical(Paths.get(expand(expected[i], inputs))));
			}
			return list;
		}
		String expectedMissingPrefix(final Path inputs){
			return missingPrefix==null ? "" : expand(missingPrefix, inputs);
		}
		private static String expand(final String value, final Path inputs){
			return value.replace("${A}", inputs.resolve("plain-alpha.fa").toString())
				.replace("${B}", inputs.resolve("plain-beta.fa").toString())
				.replace("${E1}", inputs.resolve("eq=one.fa").toString())
				.replace("${E1P}", inputs.resolve("eq").toString())
				.replace("${E2}", inputs.resolve("rhs=two.fa").toString())
				.replace("${E2P}", inputs.resolve("rhs").toString())
				.replace("${C}", inputs.resolve("literal,comma.fa").toString());
		}
		final String name, missingPrefix;
		final int failureIndex;
		final String[] args, expected;
	}

	/** Parsed output columns under test. */
	private static final class ParsedOutput{
		ParsedOutput(final ArrayList<String> filenames_, final ArrayList<Integer> lengths_, final int headerCount_){
			filenames=filenames_; lengths=lengths_; headerCount=headerCount_;
		}
		final ArrayList<String> filenames;
		final ArrayList<Integer> lengths;
		final int headerCount;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new HashMap<Path, byte[]>();
	private final Map<String, Integer> literalLengths=new HashMap<String, Integer>();
	private int cases, children, baselineFailures, successfulRows, checks, failures;
	private final boolean baseline;

	private static final Row[] ROWS={
		new Row("in-single", -1, null, new String[]{"in=${A}"}, "${A}"),
		new Row("ref-single", -1, null, new String[]{"ref=${B}"}, "${B}"),
		new Row("in-list", -1, null, new String[]{"in=${A},${B}"}, "${A}", "${B}"),
		new Row("ref-list", -1, null, new String[]{"ref=${B},${A}"}, "${B}", "${A}"),
		new Row("in-whole-comma", -1, null, new String[]{"in=${C}"}, "${C}"),
		new Row("in-equals-single", 0, "${E1P}", new String[]{"in=${E1}"}, "${E1}"),
		new Row("ref-equals-single", 0, "${E2P}", new String[]{"ref=${E2}"}, "${E2}"),
		new Row("in-equals-first-list", 0, "${E1P}", new String[]{"in=${E1},${B}"}, "${E1}", "${B}"),
		new Row("in-equals-second-list", 1, "${E2P}", new String[]{"in=${A},${E2}"}, "${A}", "${E2}"),
		new Row("in-empty", -1, null, new String[]{"in="}),
		new Row("ref-empty", -1, null, new String[]{"ref="})
	};
}
