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
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import jgi.BBMask;
import parse.Parse;
import stream.SamLine;
import structures.CoverageArray2;
import structures.IntList;

/**
 * Finite BBMask SAM coordinate controls for padded direct masks and coverage
 * ranges with soft clips, hard clips and deletion gaps.
 * @author Jean
 */
public final class BBMaskCoordinatesTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask coordinates option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmaskcoordinatestest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final BBMaskCoordinatesTest test=new BBMaskCoordinatesTest(out, repo, phase.equals("baseline"));
		test.writeInputs();
		for(final Row row : ROWS){test.runCoverage(row);}
		for(final Row row : ROWS){
			if(row.directCli){test.runDirectCli(row);}
		}
		for(final CoverageCliRow row : COVERAGE_CLI_ROWS){test.runCoverageCli(row);}
		test.check(test.coverageRows==11, "Expected eleven public increment rows: "+test.coverageRows);
		test.check(test.directChildren==10, "Expected ten direct-mask children: "+test.directChildren);
		final int expectedCoverageChildren=(test.baseline ? 0 : 1);
		test.check(test.coverageChildren==expectedCoverageChildren,
			"Expected "+expectedCoverageChildren+" coverage-threshold children: "+test.coverageChildren);
		final int expectedSkipped=(test.baseline ? 1 : 0);
		test.check(test.skipped==expectedSkipped, "Unexpected skipped coverage children: "+test.skipped);
		final int expectedSourceMatches=(test.baseline ? 21 : 13);
		final int expectedCorrectedMatches=(test.baseline ? 13 : 22);
		test.check(test.sourceMatches==expectedSourceMatches,
			"Expected "+expectedSourceMatches+" source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==expectedCorrectedMatches,
			"Expected "+expectedCorrectedMatches+" corrected matches: "+test.correctedMatches);
		final String summary="BBMASK_COORDINATES_"+(test.failures==0 ? "PASS" : "FAIL")+
		" phase="+phase+" helperRows="+test.coverageRows+
			" directChildren="+test.directChildren+" coverageChildren="+
			test.coverageChildren+" skipped="+test.skipped+
			" sourceMatches="+test.sourceMatches+" correctedMatches="+
			test.correctedMatches+" checks="+test.checks+
			" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates evidence tables and binds the script under test. */
	private BBMaskCoordinatesTest(final Path out_, final Path repo_, final boolean baseline_)
			throws IOException{
		out=out_; repo=repo_; inputs=out.resolve("inputs"); baseline=baseline_;
		launcher=repo.resolve("bbmask.sh");
		require(Files.isRegularFile(launcher), "Missing bbmask.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("coverage.tsv"), ("case\tphase\tcigar\tpos"+
			"\treflen\tsamPad\tincludeDels\tsuppliedMatch\tasserted\tpairs"+
			"\tassertionClass\tassertionMessage\tassertionTop\tassertionStack"+
			"\tsourceAsserted\tcorrectedAsserted\tsourcePairs\tcorrectedPairs"+
			"\tcoverage\tsourceCoverage\tcorrectedCoverage\tsourceMatches"+
			"\tcorrectedMatches\tmatchPhase\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		Files.write(out.resolve("direct.tsv"), ("case\tphase\texit\texpectedExit"+
			"\tcounter\tsourceCounter\tcorrectedCounter\tobservedFasta"+
			"\tsourceFasta\tcorrectedFasta\tsourceMatches\tcorrectedMatches"+
			"\tstderrMatches\tstdoutEmpty\tinputsUnchanged\tmatchPhase\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		Files.write(out.resolve("coverage-cli.tsv"), ("case\tphase\tskipped"+
			"\texit\texpectedExit\tcounter\tcorrectedCounter\tobservedFasta"+
			"\tcorrectedFasta\tcorrectedMatches\tstderrMatches\tstdoutEmpty"+
			"\tinputsUnchanged\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	/** Writes all literal FASTA and SAM inputs before any child runs. */
	private void writeInputs() throws IOException{
		for(final Row row : ROWS){
			if(row.directCli){writeSamInput(row);}
		}
	}

	/** Writes one reference/SAM pair for a row. */
	private void writeSamInput(final Row row) throws IOException{
		final Path caseInputs=inputs.resolve(row.name);
		Files.createDirectory(caseInputs);
		final byte[] fastaBytes=fasta(row.reflen).getBytes(StandardCharsets.US_ASCII);
		final byte[] samBytes=sam(row).getBytes(StandardCharsets.US_ASCII);
		final Path fasta=caseInputs.resolve("ref.fa");
		final Path sam=caseInputs.resolve("read.sam");
		Files.write(fasta, fastaBytes, StandardOpenOption.CREATE_NEW);
		Files.write(sam, samBytes, StandardOpenOption.CREATE_NEW);
		original.put(fasta, fastaBytes);
		original.put(sam, samBytes);
	}

	/** Calls public BBMask.increment() and checks exact ranges and coverage values. */
	private void runCoverage(final Row row) throws IOException{
		coverageRows++;
		final CoverageArray2 ca=new CoverageArray2(0, row.reflen);
		final IntList ranges=new IntList(8);
		AssertionError caught=null;
		try{
			BBMask.increment(ca, samLine(row), row.match(), row.reflen, ranges,
				row.includeDels, row.samPad);
		}catch(final AssertionError ae){
			caught=ae;
		}
		final boolean asserted=(caught!=null);
		final int[] observedPairs=toArray(ranges);
		final String observedCoverage=coverage(ca, row.reflen);
		final boolean sourceMatch=(assertionMatches(caught, row.sourceAsserted) &&
			Arrays.equals(observedPairs, row.sourcePairs) &&
			observedCoverage.equals(row.sourceCoverage));
		final boolean correctedMatch=(assertionMatches(caught, row.correctedAsserted) &&
			Arrays.equals(observedPairs, row.correctedPairs) &&
			observedCoverage.equals(row.correctedCoverage));
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;
		Files.write(out.resolve("coverage.tsv"), (row.name+"\t"+
			(baseline ? "baseline" : "fixed")+"\t"+row.cigar+"\t"+row.pos+
			"\t"+row.reflen+"\t"+row.samPad+"\t"+row.includeDels+"\t"+
			(row.suppliedMatch==null ? "" : row.suppliedMatch)+"\t"+asserted+
			"\t"+pairs(observedPairs)+"\t"+
			(caught==null ? "" : caught.getClass().getName())+"\t"+
			clean(caught==null || caught.getMessage()==null ? "" : caught.getMessage())+"\t"+
			clean(topFrame(caught))+"\t"+clean(stackTrace(caught))+"\t"+
			row.sourceAsserted+"\t"+row.correctedAsserted+"\t"+
			pairs(row.sourcePairs)+"\t"+pairs(row.correctedPairs)+"\t"+
			observedCoverage+"\t"+row.sourceCoverage+"\t"+
			row.correctedCoverage+"\t"+sourceMatch+"\t"+correctedMatch+
			"\t"+phaseMatch+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(phaseMatch, row.name+": public increment pairs "+pairs(observedPairs)+
			" coverage "+observedCoverage+" asserted="+asserted);
		check((observedPairs.length&1)==0, row.name+": odd number of range endpoints");
	}

	/** Launches one direct-mask bbmask.sh child and checks its exact FASTA and counter. */
	private void runDirectCli(final Row row) throws Exception{
		directChildren++;
		final Path caseDir=out.resolve("direct-"+row.name);
		Files.createDirectory(caseDir);
		writeDirectExpected(caseDir, row);
		final Path output=caseDir.resolve("out.fa");
		final ArrayList<String> command=baseCommand(row, output);
		Files.write(caseDir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		final ChildResult result=runChild(command, caseDir, output);
		final boolean sourceMatch=row.sourceDirect.matches(result.records, result.counter);
		final boolean correctedMatch=row.correctedDirect.matches(result.records, result.counter);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean stderrMatches=stderrMatches(result.stderr);
		final boolean stdoutEmpty=result.stdout.length==0;
		final boolean unchanged=inputsUnchanged();
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;
		Files.write(out.resolve("direct.tsv"), (row.name+"\t"+
			(baseline ? "baseline" : "fixed")+"\t"+result.exit+"\t0\t"+
			result.counter+"\t"+row.sourceDirect.counter+"\t"+
			row.correctedDirect.counter+"\t"+clean(serialize(result.records))+"\t"+
			clean(row.sourceDirect.serialize())+"\t"+
			clean(row.correctedDirect.serialize())+"\t"+sourceMatch+"\t"+
			correctedMatch+"\t"+stderrMatches+"\t"+stdoutEmpty+"\t"+
			unchanged+"\t"+phaseMatch+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);
		check(result.exit==0, row.name+": direct child exit was "+result.exit);
		check(phaseMatch, row.name+": direct FASTA/counter mismatch");
		check(stderrMatches, row.name+": direct stderr did not match oracle: "+result.stderr);
		check(stdoutEmpty, row.name+": direct child wrote to stdout");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Launches fixed-only coverage CLI checks that assert on the baseline source. */
	private void runCoverageCli(final CoverageCliRow row) throws Exception{
		if(baseline){
			skipped++;
			final Path caseDir=out.resolve("coverage-"+row.row.name);
			Files.createDirectory(caseDir);
			writeCoverageExpected(caseDir, row);
			Files.write(out.resolve("coverage-cli.tsv"), (row.row.name+"\tbaseline\ttrue\t"+
				"\t0\t\t"+row.corrected.counter+"\t\t"+
				clean(row.corrected.serialize())+"\tfalse\tfalse\tfalse\ttrue\n")
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
			Files.write(caseDir.resolve("SKIPPED.txt"),
				("Baseline intentionally does not execute this asserting coverage path.\n")
				.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			return;
		}
		coverageChildren++;
		final Path caseDir=out.resolve("coverage-"+row.row.name);
		Files.createDirectory(caseDir);
		writeCoverageExpected(caseDir, row);
		final Path output=caseDir.resolve("out.fa");
		final ArrayList<String> command=baseCommand(row.row, output);
		command.add("mincov=1");
		Files.write(caseDir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		final ChildResult result=runChild(command, caseDir, output);
		final boolean correctedMatch=row.corrected.matches(result.records, result.counter);
		if(correctedMatch){correctedMatches++;}
		final boolean stderrMatches=stderrMatches(result.stderr);
		final boolean stdoutEmpty=result.stdout.length==0;
		final boolean unchanged=inputsUnchanged();
		Files.write(out.resolve("coverage-cli.tsv"), (row.row.name+"\tfixed\tfalse\t"+
			result.exit+"\t0\t"+result.counter+"\t"+row.corrected.counter+"\t"+
			clean(serialize(result.records))+"\t"+clean(row.corrected.serialize())+
			"\t"+correctedMatch+"\t"+stderrMatches+"\t"+stdoutEmpty+"\t"+
			unchanged+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(result.exit==0, row.row.name+": coverage child exit was "+result.exit);
		check(correctedMatch, row.row.name+": coverage FASTA/counter mismatch");
		check(stderrMatches, row.row.name+": coverage stderr did not match oracle: "+result.stderr);
		check(stdoutEmpty, row.row.name+": coverage child wrote to stdout");
		check(unchanged, row.row.name+": a child modified literal input bytes");
	}

	/** Returns the common bbmask.sh command for a row. */
	private ArrayList<String> baseCommand(final Row row, final Path output){
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("t=1");
		command.add("simd=f");
		command.add("mr=f");
		command.add("me=f");
		command.add("minscaf=1");
		command.add("ow=t");
		command.add("sampad="+row.samPad);
		command.add("in="+input(row.name, "ref.fa"));
		command.add("out="+output);
		command.add("sam="+input(row.name, "read.sam"));
		return command;
	}

	/** Runs a child process and parses its standard artifacts. */
	private ChildResult runChild(final ArrayList<String> command, final Path caseDir,
			final Path output) throws Exception{
		final Path stdoutPath=caseDir.resolve("stdout.txt");
		final Path stderrPath=caseDir.resolve("stderr.txt");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdoutPath.toFile()).redirectError(stderrPath.toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		Files.write(caseDir.resolve("exit.txt"), (Integer.toString(exit)+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		final byte[] stdout=Files.readAllBytes(stdoutPath);
		final String stderr=new String(Files.readAllBytes(stderrPath), StandardCharsets.UTF_8);
		final ArrayList<FastaRecord> records=Files.isRegularFile(output) ?
			parseFasta(output) : new ArrayList<FastaRecord>();
		return new ChildResult(exit, stdout, stderr, parseCounter(stderr), records);
	}

	/** Writes both the baseline source oracle and corrected oracle for one direct child. */
	private static void writeDirectExpected(final Path caseDir, final Row row) throws IOException{
		final String text="sourceCounter\tsourceFasta\tcorrectedCounter\tcorrectedFasta\n"+
			row.sourceDirect.counter+"\t"+clean(row.sourceDirect.serialize())+"\t"+
			row.correctedDirect.counter+"\t"+clean(row.correctedDirect.serialize())+"\n";
		Files.write(caseDir.resolve("expected.tsv"), text.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	/** Writes the exact fixed-only oracle for one coverage child. */
	private static void writeCoverageExpected(final Path caseDir, final CoverageCliRow row)
			throws IOException{
		final String text="runInBaseline\tcorrectedCounter\tcorrectedFasta\nfalse\t"+
			row.corrected.counter+"\t"+clean(row.corrected.serialize())+"\n";
		Files.write(caseDir.resolve("expected.tsv"), text.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	private Path input(final String rowName, final String fileName){
		return inputs.resolve(rowName).resolve(fileName);
	}

	/** Builds a minimal parsed record for public increment() calls. */
	private static SamLine samLine(final Row row){
		final SamLine sl=new SamLine();
		sl.qname=row.name;
		sl.flag=0;
		sl.pos=row.pos;
		sl.mapq=60;
		sl.cigar=row.cigar;
		sl.seq=bases(queryLength(row.cigar)).getBytes(StandardCharsets.US_ASCII);
		return sl;
	}

	/** Builds one single-read SAM file for bbmask.sh children. */
	private static String sam(final Row row){
		final String seq=bases(queryLength(row.cigar));
		final char[] quals=new char[seq.length()];
		Arrays.fill(quals, 'I');
		return "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:"+row.reflen+"\n"+
			row.name+"\t0\tref\t"+row.pos+"\t60\t"+row.cigar+"\t*\t0\t0\t"+
			seq+"\t"+new String(quals)+"\n";
	}

	/** Returns A bases of the requested length. */
	private static String fasta(final int len){return ">ref\n"+bases(len)+"\n";}

	/** Returns A bases of the requested length. */
	private static String bases(final int len){
		final char[] array=new char[len];
		Arrays.fill(array, 'A');
		return new String(array);
	}

	/** Counts query bases consumed by a CIGAR, excluding H and D/N. */
	private static int queryLength(final String cigar){
		int len=0, current=0;
		for(int i=0; i<cigar.length(); i++){
			final char c=cigar.charAt(i);
			if(c>='0' && c<='9'){current=(current*10)+(c-'0');}
			else{
				if(c=='M' || c=='=' || c=='X' || c=='I' || c=='S'){len+=current;}
				current=0;
			}
		}
		return len;
	}

	/** Parses the short FASTA output generated by these cases. */
	private static ArrayList<FastaRecord> parseFasta(final Path output) throws IOException{
		final ArrayList<FastaRecord> records=new ArrayList<FastaRecord>();
		String id=null;
		final StringBuilder bases=new StringBuilder();
		for(final String line : Files.readAllLines(output, StandardCharsets.US_ASCII)){
			if(line.startsWith(">")){
				if(id!=null){records.add(new FastaRecord(id, bases.toString()));}
				id=line.substring(1);
				bases.setLength(0);
			}else{bases.append(line.trim());}
		}
		if(id!=null){records.add(new FastaRecord(id, bases.toString()));}
		return records;
	}

	/** Extracts the end-of-run Total Bases Masked value from stderr. */
	private static long parseCounter(final String errText){
		final Matcher m=COUNTER_PATTERN.matcher(errText);
		return m.find() ? Long.parseLong(m.group(1)) : Long.MIN_VALUE;
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
		for(final String part : command){sb.append(part).append('\n');}
		return sb.toString();
	}

	/** Serializes FASTA records into the compact table form. */
	private static String serialize(final List<FastaRecord> records){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<records.size(); i++){
			if(i>0){sb.append(';');}
			final FastaRecord record=records.get(i);
			sb.append(record.id).append('=').append(record.bases);
		}
		return sb.toString();
	}

	/** Escapes strings for a one-line TSV artifact. */
	private static String clean(final String text){return text.replace('\t', ' ').replace('\n', ';');}

	private static boolean assertionMatches(final AssertionError caught,
			final boolean expected){
		if(caught==null){return !expected;}
		if(!expected || caught.getClass()!=AssertionError.class){return false;}
		final StackTraceElement[] stack=caught.getStackTrace();
		return stack.length>0 && "jgi.BBMask".equals(stack[0].getClassName()) &&
			"fillRanges".equals(stack[0].getMethodName());
	}

	private static String topFrame(final AssertionError caught){
		if(caught==null){return "";}
		final StackTraceElement[] stack=caught.getStackTrace();
		return stack.length<1 ? "" : frame(stack[0]);
	}

	private static String stackTrace(final AssertionError caught){
		if(caught==null){return "";}
		final StringBuilder sb=new StringBuilder();
		for(final StackTraceElement element : caught.getStackTrace()){
			if(sb.length()>0){sb.append(';');}
			sb.append(frame(element));
		}
		return sb.toString();
	}

	private static String frame(final StackTraceElement element){
		return element.getClassName()+"."+element.getMethodName()+"("+
			element.getFileName()+":"+element.getLineNumber()+")";
	}

	/** Records a failed invariant without aborting subsequent checks. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println(message);}
	}

	/** Fails fast for invalid harness setup. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Copies IntList contents to a compact array. */
	private static int[] toArray(final IntList list){
		final int[] array=new int[list.size];
		for(int i=0; i<array.length; i++){array[i]=list.get(i);}
		return array;
	}

	/** Serializes adjacent range endpoints as half-open pairs. */
	private static String pairs(final int[] array){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<array.length; i+=2){
			if(i>0){sb.append(',');}
			sb.append('[').append(array[i]).append(',').append(array[i+1]).append(')');
		}
		return sb.toString();
	}

	/** Serializes a coverage array as exact comma-delimited values. */
	private static String coverage(final CoverageArray2 ca, final int len){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<len; i++){
			if(i>0){sb.append(',');}
			sb.append(ca.get(i));
		}
		return sb.toString();
	}

	private static boolean stderrMatches(final String errText){
		return errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
			!errText.contains("Exception") && !errText.contains("AssertionError") &&
			!errText.contains("Unknown parameter") && !errText.contains("Error reading");
	}

	/** One output FASTA record. */
	private static final class FastaRecord{
		FastaRecord(final String id_, final String bases_){id=id_; bases=bases_;}
		@Override
		public boolean equals(final Object other){
			if(this==other){return true;}
			if(!(other instanceof FastaRecord)){return false;}
			final FastaRecord f=(FastaRecord)other;
			return id.equals(f.id) && bases.equals(f.bases);
		}
		@Override
		public int hashCode(){return id.hashCode()*31+bases.hashCode();}
		final String id, bases;
	}

	/** Exact FASTA/counter oracle for one child. */
	private static final class Expect{
		Expect(final long counter_, final FastaRecord record_){counter=counter_; record=record_;}
		boolean matches(final ArrayList<FastaRecord> records, final long observedCounter){
			return counter==observedCounter && records.size()==1 && record.equals(records.get(0));
		}
		String serialize(){return BBMaskCoordinatesTest.serialize(Arrays.asList(record));}
		final long counter;
		final FastaRecord record;
	}

	/** One child process result. */
	private static final class ChildResult{
		ChildResult(final int exit_, final byte[] stdout_, final String stderr_,
				final long counter_, final ArrayList<FastaRecord> records_){
			exit=exit_; stdout=stdout_; stderr=stderr_; counter=counter_; records=records_;
		}
		final int exit;
		final byte[] stdout;
		final String stderr;
		final long counter;
		final ArrayList<FastaRecord> records;
	}

	/** One finite coordinate row. */
	private static final class Row{
		Row(final String name_, final int pos_, final String cigar_, final int reflen_,
				final int samPad_, final boolean includeDels_, final String suppliedMatch_,
				final boolean directCli_,
				final Expect sourceDirect_, final Expect correctedDirect_,
				final boolean sourceAsserted_, final boolean correctedAsserted_,
				final int[] sourcePairs_, final int[] correctedPairs_,
				final String sourceCoverage_, final String correctedCoverage_){
			name=name_; pos=pos_; cigar=cigar_; reflen=reflen_; samPad=samPad_;
			includeDels=includeDels_; suppliedMatch=suppliedMatch_; directCli=directCli_;
			sourceDirect=sourceDirect_; correctedDirect=correctedDirect_;
			sourceAsserted=sourceAsserted_; correctedAsserted=correctedAsserted_;
			sourcePairs=sourcePairs_; correctedPairs=correctedPairs_;
			sourceCoverage=sourceCoverage_; correctedCoverage=correctedCoverage_;
		}
		byte[] match(){
			return suppliedMatch==null ? null : suppliedMatch.getBytes(StandardCharsets.US_ASCII);
		}
		final String name, cigar, suppliedMatch;
		final int pos, reflen, samPad;
		final boolean includeDels, sourceAsserted, correctedAsserted, directCli;
		final Expect sourceDirect, correctedDirect;
		final int[] sourcePairs, correctedPairs;
		final String sourceCoverage, correctedCoverage;
	}

	/** One fixed-only coverage child row. */
	private static final class CoverageCliRow{
		CoverageCliRow(final Row row_, final Expect corrected_){
			row=row_; corrected=corrected_;
		}
		final Row row;
		final Expect corrected;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new LinkedHashMap<Path, byte[]>();
	private int coverageRows, directChildren, coverageChildren, skipped;
	private int sourceMatches, correctedMatches, checks, failures;
	private final boolean baseline;

	private static final Pattern COUNTER_PATTERN=Pattern.compile("Total Bases Masked:\\s+(\\d+)/");
	private static final String ZERO12="0,0,0,0,0,0,0,0,0,0,0,0";
	private static final String ZERO20="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0";

	private static final Row ROW_RIGHT_SOFT_TRUNCATED=new Row(
		"right-soft-truncated", 9, "4M2S", 12, 0, true, null, true,
		new Expect(4, new FastaRecord("ref", "AAAAAAAANNNN")),
		new Expect(4, new FastaRecord("ref", "AAAAAAAANNNN")),
		true, false, new int[] {8, 12}, new int[] {8, 12},
		ZERO12, "0,0,0,0,0,0,0,0,1,1,1,1");

	private static final Row[] ROWS={
		new Row("pad0-middle", 10, "5M", 20, 0, true, null, true,
			new Expect(5, new FastaRecord("ref", "AAAAAAAAANNNNNAAAAAA")),
			new Expect(5, new FastaRecord("ref", "AAAAAAAAANNNNNAAAAAA")),
			false, false, new int[] {9, 14}, new int[] {9, 14},
			"0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0",
			"0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0"),
		new Row("pad2-middle", 10, "5M", 20, 2, true, null, true,
			new Expect(7, new FastaRecord("ref", "AAAAAAANNNNNNNAAAAAA")),
			new Expect(9, new FastaRecord("ref", "AAAAAAANNNNNNNNNAAAA")),
			false, false, new int[] {7, 16}, new int[] {7, 16},
			"0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0",
			"0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0"),
		new Row("left-boundary-pad", 2, "4M", 20, 2, true, null, true,
			new Expect(6, new FastaRecord("ref", "NNNNNNAAAAAAAAAAAAAA")),
			new Expect(7, new FastaRecord("ref", "NNNNNNNAAAAAAAAAAAAA")),
			false, false, new int[] {0, 7}, new int[] {0, 7},
			"1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0",
			"1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0"),
		new Row("left-soft-unpadded", 1, "2S4M", 20, 0, true, null, true,
			new Expect(6, new FastaRecord("ref", "NNNNNNAAAAAAAAAAAAAA")),
			new Expect(4, new FastaRecord("ref", "NNNNAAAAAAAAAAAAAAAA")),
			false, false, new int[] {2, 6}, new int[] {0, 4},
			"0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
			"1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"),
		new Row("left-soft-padded", 1, "2S4M", 20, 2, true, null, true,
			new Expect(8, new FastaRecord("ref", "NNNNNNNNAAAAAAAAAAAA")),
			new Expect(6, new FastaRecord("ref", "NNNNNNAAAAAAAAAAAAAA")),
			false, false, new int[] {0, 8}, new int[] {0, 6},
			"1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0",
			"1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0"),
		new Row("right-soft-inbounds", 5, "2M2S", 20, 0, true, null, true,
			new Expect(4, new FastaRecord("ref", "AAAANNNNAAAAAAAAAAAA")),
			new Expect(4, new FastaRecord("ref", "AAAANNNNAAAAAAAAAAAA")),
			false, false, new int[] {4, 6}, new int[] {4, 6},
			"0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
			"0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0"),
		ROW_RIGHT_SOFT_TRUNCATED,
		new Row("hard-clip-control", 10, "2H4M2H", 20, 0, true, null, true,
			new Expect(4, new FastaRecord("ref", "AAAAAAAAANNNNAAAAAAA")),
			new Expect(4, new FastaRecord("ref", "AAAAAAAAANNNNAAAAAAA")),
			false, false, new int[] {9, 13}, new int[] {9, 13},
			"0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0",
			"0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0"),
		new Row("excluded-deletion-pad", 10, "2M2D2M", 20, 1, false, null, true,
			new Expect(7, new FastaRecord("ref", "AAAAAAAANNNNNNNAAAAA")),
			new Expect(8, new FastaRecord("ref", "AAAAAAAANNNNNNNNAAAA")),
			false, false, new int[] {8, 12, 12, 16}, new int[] {8, 12, 12, 16},
			"0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0",
			"0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0"),
		new Row("supplied-excluded-deletion-pad", 10, "2M2D2M", 20, 1, false, "mmDDmm", false,
			new Expect(0, new FastaRecord("ref", "")),
			new Expect(0, new FastaRecord("ref", "")),
			false, false, new int[] {8, 12, 12, 16}, new int[] {8, 12, 12, 16},
			"0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0",
			"0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0"),
		new Row("right-boundary-pad", 18, "3M", 20, 2, true, null, true,
			new Expect(5, new FastaRecord("ref", "AAAAAAAAAAAAAAANNNNN")),
			new Expect(5, new FastaRecord("ref", "AAAAAAAAAAAAAAANNNNN")),
			false, false, new int[] {15, 20}, new int[] {15, 20},
			"0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1",
			"0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1")
	};

	private static final CoverageCliRow[] COVERAGE_CLI_ROWS={
		new CoverageCliRow(ROW_RIGHT_SOFT_TRUNCATED,
			new Expect(8, new FastaRecord("ref", "NNNNNNNNAAAA")))
	};
}
