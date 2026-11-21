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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import parse.Parse;

/**
 * Ten fresh BBMask split-mode child-process cases for existing and generated
 * FASTA header collisions, split-piece numbering, and split counter parity.
 * @author Jean
 */
public final class BBMaskSplitTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask split option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmasksplittest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final BBMaskSplitTest test=new BBMaskSplitTest(out, repo, phase.equals("baseline"));
		test.writeInputs();
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==10, "Expected ten BBMask split cases: "+test.cases);
		test.check(test.children==10, "Expected one fresh bbmask.sh child per case: "+test.children);
		test.check(test.sourceMatches==(test.baseline ? 10 : 4),
			"Unexpected source-derived matches: "+test.sourceMatches);
		test.check(test.correctedMatches==(test.baseline ? 4 : 10),
			"Unexpected corrected-oracle matches: "+test.correctedMatches);
		final String summary="BBMASK_SPLIT_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" sourceMatches="+test.sourceMatches+" correctedMatches="+test.correctedMatches+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private BBMaskSplitTest(final Path out_, final Path repo_, final boolean baseline_) throws IOException{
		out=out_; repo=repo_; inputs=out.resolve("inputs"); baseline=baseline_;
		launcher=repo.resolve("bbmask.sh");
		require(Files.isRegularFile(launcher), "Missing bbmask.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\texpectedExit\tcounter"+
			"\texpectedCounter\toutputExists\tstdoutBytes\tstderrBytes\tobservedFasta"+
			"\tsourceMatches\tcorrectedMatches\toutputMatches\tcounterMatches"+
			"\tstderrMatches\tstdoutEmpty\tinputsUnchanged\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	/** Writes all literal FASTA and optional SAM inputs before any child runs. */
	private void writeInputs() throws IOException{
		for(final Row row : ROWS){
			final Path fasta=input(row.inputName);
			final byte[] fastaBytes=row.fastaText().getBytes(StandardCharsets.US_ASCII);
			Files.write(fasta, fastaBytes, StandardOpenOption.CREATE_NEW);
			original.put(fasta, fastaBytes);
			if(row.sam!=null){
				final Path sam=sam(row.name);
				final byte[] samBytes=row.sam.text().getBytes(StandardCharsets.US_ASCII);
				Files.write(sam, samBytes, StandardOpenOption.CREATE_NEW);
				original.put(sam, samBytes);
			}
		}
	}

	/** Launches one bbmask.sh child and checks baseline or corrected split behavior. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		final Expect source=row.source;
		final Expect corrected=row.corrected;
		final Expect expected=baseline ? source : corrected;
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final Path output=caseDir.resolve("out.fa");
		writeExpected(caseDir, source, corrected);

		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("t=1");
		command.add("simd=f");
		command.add("split=t");
		command.add("mr=f");
		command.add("me=f");
		command.add("minscaf=1");
		command.add("ow=t");
		command.add("in="+input(row.inputName));
		command.add("out="+output);
		if(row.lowercase){command.add("lc=t");}
		if(row.sam!=null){command.add("sam="+sam(row.name));}
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
		final boolean sourceMatch=source.matches(exit, outputExists, records, observedCounter, errText);
		final boolean correctedMatch=corrected.matches(exit, outputExists, records, observedCounter, errText);
		final boolean outputMatches=expected.recordsEqual(records);
		final boolean counterMatches=expected.counterMatches(observedCounter);
		final boolean stderrMatches=expected.stderrMatches(errText);
		final boolean outputFileMatches=expected.outputFileMatches(outputExists);
		final boolean stdoutEmpty=stdout.length==0;
		final boolean unchanged=inputsUnchanged();
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}

		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t"+(expected.exitZero ? 0 : 1)+"\t"+
			(observedCounter==NO_COUNTER ? "" : Long.toString(observedCounter))+"\t"+
			(expected.counter==NO_COUNTER ? "" : Long.toString(expected.counter))+"\t"+
			outputExists+"\t"+stdout.length+"\t"+stderr.length+"\t"+clean(serialize(records))+"\t"+
			sourceMatch+"\t"+correctedMatch+"\t"+outputMatches+"\t"+counterMatches+
			"\t"+stderrMatches+"\t"+stdoutEmpty+"\t"+unchanged+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);

		check((exit==0)==expected.exitZero, row.name+": child exit was "+exit);
		check(outputFileMatches, row.name+": output file existence did not match oracle");
		check(outputMatches, row.name+": output FASTA mismatch: "+serialize(records));
		check(counterMatches, row.name+": Total Bases Masked mismatch: "+observedCounter);
		check(stderrMatches, row.name+": stderr did not match source-derived oracle: "+errText);
		check(stdoutEmpty, row.name+": bbmask.sh wrote to stdout");
		check(command.contains("-ea") && command.contains("t=1") &&
			command.contains("simd=f") && command.contains("split=t") &&
			command.contains("mr=f") && command.contains("me=f") &&
			command.contains("minscaf=1"), row.name+": child argv lost a frozen CLI pin");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Writes the current-source and corrected oracles before launching each child. */
	private static void writeExpected(final Path caseDir, final Expect source,
			final Expect corrected) throws IOException{
		final String text="oracle\texit\tcounter\tfasta\tfailure\n"+
			"source\t"+(source.exitZero ? 0 : 1)+"\t"+
			(source.counter==NO_COUNTER ? "" : Long.toString(source.counter))+"\t"+
			clean(source.serialize())+"\t"+source.failureToken+"\n"+
			"corrected\t"+(corrected.exitZero ? 0 : 1)+"\t"+
			(corrected.counter==NO_COUNTER ? "" : Long.toString(corrected.counter))+"\t"+
			clean(corrected.serialize())+"\t"+corrected.failureToken+"\n";
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

	/** Joins argv for stable evidence. */
	private static String join(final List<String> command){
		final StringBuilder sb=new StringBuilder();
		for(final String s : command){sb.append(s).append('\n');}
		return sb.toString();
	}

	private static String serialize(final List<FastaRecord> records){
		final StringBuilder sb=new StringBuilder();
		for(final FastaRecord record : records){
			if(sb.length()>0){sb.append(';');}
			sb.append(record.id).append('=').append(record.bases);
		}
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

	private Path input(final String name){return inputs.resolve(name+".fa");}

	private Path sam(final String name){return inputs.resolve(name+".sam");}

	/** One input FASTA record. */
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

	/** Optional SAM alignment that masks a tail range. */
	private static final class SamCase{

		SamCase(final int pos_, final String cigar_, final String seq_){
			pos=pos_; cigar=cigar_; seq=seq_;
		}

		String text(){
			final char[] quals=new char[seq.length()];
			Arrays.fill(quals, 'I');
			return "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:4\n"+
				"r1\t0\tref\t"+pos+"\t60\t"+cigar+"\t*\t0\t0\t"+seq+"\t"+
				new String(quals)+"\n";
		}

		final int pos;
		final String cigar, seq;
	}

	/** Current-source or corrected oracle for one case. */
	private static final class Expect{

		Expect(final long counter_, final String failureToken_, final FastaRecord... records_){
			exitZero=(failureToken_==null);
			counter=counter_;
			failureToken=failureToken_==null ? "" : failureToken_;
			records=records_;
		}

		boolean matches(final int exit, final boolean outputExists, final ArrayList<FastaRecord> observed,
				final long observedCounter, final String errText){
			return (exit==0)==exitZero && outputFileMatches(outputExists) && recordsEqual(observed) &&
				counterMatches(observedCounter) && stderrMatches(errText);
		}

		boolean recordsEqual(final ArrayList<FastaRecord> observed){
			return Arrays.equals(records, observed.toArray(new FastaRecord[observed.size()]));
		}

		boolean counterMatches(final long observed){
			return counter==observed;
		}

		boolean outputFileMatches(final boolean exists){
			return exists==exitZero;
		}

		boolean stderrMatches(final String errText){
			final boolean common=errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
				!errText.contains("Unknown parameter") && !errText.contains("Error reading");
			if(!common){return false;}
			if(exitZero){return !errText.contains("Exception") && !errText.contains("AssertionError");}
			return errText.contains("AssertionError") && errText.contains(failureToken);
		}

		String serialize(){return BBMaskSplitTest.serialize(Arrays.asList(records));}

		final boolean exitZero;
		final long counter;
		final String failureToken;
		final FastaRecord[] records;
	}

	/** One bbmask.sh child-process case. */
	private static final class Row{

		Row(final String name_, final String inputName_, final boolean lowercase_,
				final SamCase sam_, final Expect source_, final Expect corrected_,
				final FastaRecord... inputRecords_){
			name=name_; inputName=inputName_; lowercase=lowercase_; sam=sam_;
			source=source_; corrected=corrected_; inputRecords=inputRecords_;
		}

		String fastaText(){
			final StringBuilder sb=new StringBuilder();
			for(final FastaRecord record : inputRecords){
				sb.append('>').append(record.id).append('\n').append(record.bases).append('\n');
			}
			return sb.toString();
		}

		final String name, inputName;
		final boolean lowercase;
		final SamCase sam;
		final Expect source, corrected;
		final FastaRecord[] inputRecords;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new HashMap<Path, byte[]>();
	private int cases, children, sourceMatches, correctedMatches, checks, failures;
	private final boolean baseline;

	private static final long NO_COUNTER=Long.MIN_VALUE;
	private static final Pattern COUNTER_PATTERN=Pattern.compile("Total Bases Masked:\\s+(\\d+)/");

	private static final FastaRecord S=new FastaRecord("s", "ACGT");
	private static final FastaRecord S_0_A=new FastaRecord("s_0", "A");
	private static final FastaRecord S_0_AC=new FastaRecord("s_0", "AC");
	private static final FastaRecord S_0_TT=new FastaRecord("s_0", "TT");
	private static final FastaRecord S_1_AC=new FastaRecord("s_1", "AC");
	private static final FastaRecord S_1_GT=new FastaRecord("s_1", "GT");
	private static final FastaRecord S_1_T=new FastaRecord("s_1", "T");
	private static final FastaRecord S_2_GT=new FastaRecord("s_2", "GT");
	private static final FastaRecord REF_0_A=new FastaRecord("ref_0", "A");
	private static final FastaRecord REF_0_AC=new FastaRecord("ref_0", "AC");

	private static final Expect DUPLICATE=new Expect(NO_COUNTER, "Duplicate id");
	private static final Expect COUNTER_2=new Expect(NO_COUNTER, "0, 2, 0, 2, 0");
	private static final Expect COUNTER_1=new Expect(NO_COUNTER, "0, 1, 0, 1, 0");
	private static final Expect EMPTY=new Expect(0, null);

	private static final Row[] ROWS={
		new Row("unsplit", "unsplit", false, null, new Expect(0, null, S),
			new Expect(0, null, S), S),
		new Row("leading", "leading", false, null, new Expect(0, null, S_0_AC),
			new Expect(0, null, S_0_AC), new FastaRecord("s", "NNAC")),
		new Row("trailing", "trailing", false, null, new Expect(0, null, S_0_AC),
			new Expect(0, null, S_0_AC), new FastaRecord("s", "ACNN")),
		new Row("internal", "internal", false, null, DUPLICATE,
			new Expect(0, null, S_0_AC, S_1_GT), new FastaRecord("s", "ACNNGT")),
		new Row("all-n", "all-n", false, null, EMPTY, EMPTY, new FastaRecord("s", "NNNN")),
		new Row("existing-before", "existing-before", false, null, DUPLICATE,
			new Expect(0, null, S_0_TT, S_1_AC, S_2_GT),
			new FastaRecord("s_0", "TT"), new FastaRecord("s", "ACNNGT")),
		new Row("existing-after", "existing-after", false, null,
			new Expect(0, null, S_0_TT), new Expect(0, null, S_1_AC, S_0_TT),
			new FastaRecord("s", "ACNN"), new FastaRecord("s_0", "TT")),
		new Row("sam-tail", "sam-tail", false, new SamCase(3, "2M", "GT"), COUNTER_2,
			new Expect(2, null, REF_0_AC), new FastaRecord("ref", "ACGT")),
		new Row("n-plus-sam-tail", "n-plus-sam-tail", false, new SamCase(4, "1M", "T"),
			COUNTER_1, new Expect(1, null, REF_0_A), new FastaRecord("ref", "ANNT")),
		new Row("lowercase-preexisting", "lowercase-preexisting", true, null, DUPLICATE,
			new Expect(0, null, S_0_A, S_1_T), new FastaRecord("s", "AcgtT")),
	};
}
