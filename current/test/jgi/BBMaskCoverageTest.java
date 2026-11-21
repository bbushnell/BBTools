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

import parse.Parse;

/**
 * Fresh BBMask SAM coverage-threshold child-process cases for one-file controls
 * and fixed-only multi-file aggregate coverage checks.
 * @author Jean
 */
public final class BBMaskCoverageTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask coverage option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmaskcoveragetest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final BBMaskCoverageTest test=new BBMaskCoverageTest(out, repo, phase.equals("baseline"));
		test.writeInputs();
		for(final Row row : ROWS){test.run(row);}
		final int expectedCases=(test.baseline ? 3 : 7);
		test.check(test.cases==expectedCases, "Expected "+expectedCases+" BBMask coverage cases: "+test.cases);
		test.check(test.children==expectedCases, "Expected "+expectedCases+" fresh bbmask.sh children: "+test.children);
		final int expectedSkipped=(test.baseline ? 4 : 0);
		test.check(test.skipped==expectedSkipped, "Unexpected skipped cases: "+test.skipped);
		final String summary="BBMASK_COVERAGE_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" skipped="+test.skipped+" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private BBMaskCoverageTest(final Path out_, final Path repo_, final boolean baseline_)
			throws IOException{
		out=out_; repo=repo_; inputs=out.resolve("inputs"); baseline=baseline_;
		launcher=repo.resolve("bbmask.sh");
		require(Files.isRegularFile(launcher), "Missing bbmask.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\texpectedExit\tcounter"+
			"\texpectedCounter\toutputExists\tstdoutBytes\tstderrBytes\tobservedFasta"+
			"\toutputMatches\tcounterMatches\tstderrMatches\tstdoutEmpty\tinputsUnchanged\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes all literal FASTA and SAM inputs before any child runs. */
	private void writeInputs() throws IOException{
		for(final Row row : ROWS){
			final Path caseInputs=inputs.resolve(row.name);
			Files.createDirectory(caseInputs);
			final Path fasta=caseInputs.resolve("ref.fa");
			final byte[] fastaBytes=FASTA.getBytes(StandardCharsets.US_ASCII);
			Files.write(fasta, fastaBytes, StandardOpenOption.CREATE_NEW);
			original.put(fasta, fastaBytes);
			for(final SamFile sam : row.sams){
				final Path samPath=caseInputs.resolve(sam.name);
				final byte[] samBytes=sam.text().getBytes(StandardCharsets.US_ASCII);
				Files.write(samPath, samBytes, StandardOpenOption.CREATE_NEW);
				original.put(samPath, samBytes);
			}
		}
	}

	/** Launches one safe bbmask.sh child and checks its exact FASTA and counter. */
	private void run(final Row row) throws Exception{
		if(baseline && !row.baseline){
			skipped++;
			writeSkipped(row);
			return;
		}
		cases++;
		children++;
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		writeExpected(caseDir, row);
		final Path output=caseDir.resolve("out.fa");

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
		command.add("in="+input(row.name, "ref.fa"));
		command.add("out="+output);
		if(row.mincov>-1){command.add("mincov="+row.mincov);}
		if(row.maxcov>-1){command.add("maxcov="+row.maxcov);}
		command.add("sam="+samList(row));
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

		final boolean outputExists=Files.isRegularFile(output);
		final ArrayList<FastaRecord> records=outputExists ? parseFasta(output) : new ArrayList<FastaRecord>();
		final long observedCounter=parseCounter(errText);
		final boolean outputMatches=row.expected.recordsEqual(records);
		final boolean counterMatches=row.expected.counterMatches(observedCounter);
		final boolean stderrMatches=row.expected.stderrMatches(errText);
		final boolean stdoutEmpty=stdout.length==0;
		final boolean unchanged=inputsUnchanged();

		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t0\t"+(observedCounter==NO_COUNTER ? "" : Long.toString(observedCounter))+
			"\t"+row.expected.counter+"\t"+outputExists+"\t"+stdout.length+"\t"+
			stderr.length+"\t"+clean(serialize(records))+"\t"+outputMatches+"\t"+
			counterMatches+"\t"+stderrMatches+"\t"+stdoutEmpty+"\t"+unchanged+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);

		check(exit==0, row.name+": child exit was "+exit);
		check(outputExists, row.name+": no FASTA output was written");
		check(outputMatches, row.name+": output FASTA mismatch: "+serialize(records));
		check(counterMatches, row.name+": Total Bases Masked mismatch: "+observedCounter);
		check(stderrMatches, row.name+": stderr did not match oracle: "+errText);
		check(stdoutEmpty, row.name+": bbmask.sh wrote to stdout");
		check(command.contains("-ea") && command.contains("t=1") &&
			command.contains("simd=f") && command.contains("mr=f") &&
			command.contains("me=f") && command.contains("minscaf=1") &&
			!command.contains("split=t"), row.name+": child argv lost a frozen CLI pin");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Records that the unsafe source-derived affected path is deliberately not run. */
	private void writeSkipped(final Row row) throws IOException{
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		writeExpected(caseDir, row);
		Files.write(caseDir.resolve("SKIPPED.txt"),
			("Baseline intentionally does not execute this affected multi-SAM threshold path.\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes the exact fixed oracle for one child. */
	private static void writeExpected(final Path caseDir, final Row row) throws IOException{
		final String text="runInBaseline\tcounter\tfasta\tdescription\n"+
			row.baseline+"\t"+row.expected.counter+"\t"+
			clean(row.expected.serialize())+"\t"+row.description+"\n";
		Files.write(caseDir.resolve("expected.tsv"), text.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
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
			}else{
				bases.append(line.trim());
			}
		}
		if(id!=null){records.add(new FastaRecord(id, bases.toString()));}
		return records;
	}

	/** Extracts the end-of-run Total Bases Masked value from stderr. */
	private static long parseCounter(final String errText){
		final Matcher m=COUNTER_PATTERN.matcher(errText);
		return m.find() ? Long.parseLong(m.group(1)) : NO_COUNTER;
	}

	/** Checks every retained input against the literal bytes written before children ran. */
	private boolean inputsUnchanged() throws IOException{
		for(final Map.Entry<Path, byte[]> entry : original.entrySet()){
			if(!Arrays.equals(entry.getValue(), Files.readAllBytes(entry.getKey()))){return false;}
		}
		return true;
	}

	/** Returns a comma-delimited sam= value in the row's exact order. */
	private String samList(final Row row){
		final ArrayList<String> paths=new ArrayList<String>(row.sams.length);
		for(final SamFile sam : row.sams){paths.add(input(row.name, sam.name).toString());}
		return join(paths, ",");
	}

	private Path input(final String rowName, final String fileName){return inputs.resolve(rowName).resolve(fileName);}

	/** Joins argv for stable evidence. */
	private static String join(final List<String> command){return join(command, " ")+"\n";}

	/** Joins values with a delimiter. */
	private static String join(final List<String> values, final String delimiter){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<values.size(); i++){
			if(i>0){sb.append(delimiter);}
			sb.append(values.get(i));
		}
		return sb.toString();
	}

	/** Serializes FASTA records into the compact cases.tsv form. */
	private static String serialize(final List<FastaRecord> records){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<records.size(); i++){
			final FastaRecord record=records.get(i);
			if(i>0){sb.append(';');}
			sb.append(record.id).append('=').append(record.bases);
		}
		return sb.toString();
	}

	/** Escapes strings for a one-line TSV artifact. */
	private static String clean(final String text){
		return text.replace('\t', ' ').replace('\n', ';');
	}

	/** Records a failed invariant without aborting subsequent checks. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){
			failures++;
			System.err.println(message);
		}
	}

	/** Fails fast for invalid harness setup. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
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

	/** One literal SAM file. */
	private static final class SamFile{

		SamFile(final String name_, final SamRecord... records_){name=name_; records=records_;}

		String text(){
			final StringBuilder sb=new StringBuilder(SAM_HEADER);
			for(final SamRecord record : records){sb.append(record.text());}
			return sb.toString();
		}

		final String name;
		final SamRecord[] records;
	}

	/** One SAM alignment to the four-base ref. */
	private static final class SamRecord{

		SamRecord(final String read_, final int pos_, final String seq_){read=read_; pos=pos_; seq=seq_;}

		String text(){
			final char[] quals=new char[seq.length()];
			Arrays.fill(quals, 'I');
			return read+"\t0\tref\t"+pos+"\t60\t"+seq.length()+"M\t*\t0\t0\t"+seq+"\t"+
				new String(quals)+"\n";
		}

		final String read, seq;
		final int pos;
	}

	/** Exact oracle for one successful child. */
	private static final class Expect{

		Expect(final long counter_, final FastaRecord... records_){counter=counter_; records=records_;}

		boolean recordsEqual(final ArrayList<FastaRecord> observed){
			return Arrays.equals(records, observed.toArray(new FastaRecord[observed.size()]));
		}

		boolean counterMatches(final long observed){return counter==observed;}

		boolean stderrMatches(final String errText){
			return errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
				!errText.contains("Exception") && !errText.contains("AssertionError") &&
				!errText.contains("Unknown parameter") && !errText.contains("Error reading");
		}

		String serialize(){return BBMaskCoverageTest.serialize(Arrays.asList(records));}

		final long counter;
		final FastaRecord[] records;
	}

	/** One child-process case. */
	private static final class Row{

		Row(final String name_, final boolean baseline_, final int mincov_, final int maxcov_,
				final Expect expected_, final String description_, final SamFile... sams_){
			name=name_; baseline=baseline_; mincov=mincov_; maxcov=maxcov_;
			expected=expected_; description=description_; sams=sams_;
		}

		final String name, description;
		final boolean baseline;
		final int mincov, maxcov;
		final Expect expected;
		final SamFile[] sams;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new LinkedHashMap<Path, byte[]>();
	private int cases, children, skipped, checks, failures;
	private final boolean baseline;

	private static final long NO_COUNTER=Long.MIN_VALUE;
	private static final Pattern COUNTER_PATTERN=Pattern.compile("Total Bases Masked:\\s+(\\d+)/");

	private static final String FASTA=">ref\nACGT\n";
	private static final String SAM_HEADER="@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:4\n";
	private static final FastaRecord MINCOV_OUT=new FastaRecord("ref", "ANNT");
	private static final FastaRecord MAXCOV_OUT=new FastaRecord("ref", "NCGT");
	private static final FastaRecord DIRECT_OUT=new FastaRecord("ref", "NCGN");
	private static final Expect MINCOV_EXPECT=new Expect(2, MINCOV_OUT);
	private static final Expect MAXCOV_EXPECT=new Expect(1, MAXCOV_OUT);
	private static final Expect DIRECT_EXPECT=new Expect(2, DIRECT_OUT);

	private static final SamRecord POS1_A=new SamRecord("pos1", 1, "A");
	private static final SamRecord POS1_B=new SamRecord("pos1b", 1, "A");
	private static final SamRecord POS4_T=new SamRecord("pos4", 4, "T");

	private static final Row[] ROWS={
		new Row("single-mincov", true, 1, -1, MINCOV_EXPECT,
			"single combined SAM covers endpoints; mincov masks uncovered middle",
			new SamFile("combined.sam", POS1_A, POS4_T)),
		new Row("single-maxcov", true, -1, 1, MAXCOV_EXPECT,
			"single combined SAM covers position one twice; maxcov masks duplicate position",
			new SamFile("combined.sam", POS1_A, POS1_B)),
		new Row("direct-two-file", true, -1, -1, DIRECT_EXPECT,
			"two SAM files without thresholds direct-mask both covered endpoints",
			new SamFile("one.sam", POS1_A), new SamFile("two.sam", POS4_T)),
		new Row("mincov-forward", false, 1, -1, MINCOV_EXPECT,
			"fixed-only aggregate mincov with forward file order",
			new SamFile("one.sam", POS1_A), new SamFile("two.sam", POS4_T)),
		new Row("mincov-reverse", false, 1, -1, MINCOV_EXPECT,
			"fixed-only aggregate mincov with reverse file order",
			new SamFile("two.sam", POS4_T), new SamFile("one.sam", POS1_A)),
		new Row("maxcov-forward", false, -1, 1, MAXCOV_EXPECT,
			"fixed-only aggregate maxcov with duplicate position in forward file order",
			new SamFile("one.sam", POS1_A), new SamFile("two.sam", POS1_B)),
		new Row("maxcov-reverse", false, -1, 1, MAXCOV_EXPECT,
			"fixed-only aggregate maxcov with duplicate position in reverse file order",
			new SamFile("two.sam", POS1_B), new SamFile("one.sam", POS1_A)),
	};
}
