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
 * Finite BBMask constructor parsing controls for named filenames containing
 * additional '=' characters.
 * @author Jean
 */
public final class BBMaskFirstEqualsTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask first-equals option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmaskfirstequalstest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final BBMaskFirstEqualsTest test=new BBMaskFirstEqualsTest(out, repo,
			phase.equals("baseline"));
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==5, "Expected five BBMask first-equals cases: "+test.cases);
		test.check(test.children==5, "Expected five fresh bbmask.sh children: "+test.children);
		test.check(test.changed==3, "Expected three behavioral rows: "+test.changed);
		test.check(test.unchanged==2, "Expected two unchanged controls: "+test.unchanged);
		final int expectedSourceMatches=(test.baseline ? 5 : 2);
		final int expectedCorrectedMatches=(test.baseline ? 2 : 5);
		test.check(test.sourceMatches==expectedSourceMatches,
			"Expected "+expectedSourceMatches+" source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==expectedCorrectedMatches,
			"Expected "+expectedCorrectedMatches+" corrected matches: "+test.correctedMatches);
		final String summary="BBMASK_FIRST_EQUALS_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" changed="+test.changed+" unchanged="+test.unchanged+
			" sourceMatches="+test.sourceMatches+
			" correctedMatches="+test.correctedMatches+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates the evidence table and binds the script under test. */
	private BBMaskFirstEqualsTest(final Path out_, final Path repo_, final boolean baseline_)
			throws IOException{
		out=out_; repo=repo_; baseline=baseline_;
		launcher=repo.resolve("bbmask.sh");
		require(Files.isRegularFile(launcher), "Missing bbmask.sh at "+launcher);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\texpectedExit\tcounter"+
			"\tsourceCounter\tcorrectedCounter\tobservedDestination\tobservedFasta"+
			"\tsourceFasta\tcorrectedFasta\tintendedExists\ttruncatedExists"+
			"\tsourceIntendedExists\tsourceTruncatedExists\tcorrectedIntendedExists"+
			"\tcorrectedTruncatedExists\tsourceMatches\tcorrectedMatches"+
			"\tstderrMatches\tstdoutEmpty\tinputsUnchanged\tmatchPhase\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Launches one bbmask.sh child and checks its exact FASTA, counter and destination. */
	private void run(final Row row) throws Exception{
		cases++;
		children++;
		if(row.source.changed(row.corrected)){changed++;}else{unchanged++;}
		final Path caseDir=out.resolve(row.name);
		Files.createDirectory(caseDir);
		final CaseFiles files=writeInputs(caseDir, row);
		final ArrayList<String> command=baseCommand(row, files);
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

		final boolean intendedExists=Files.isRegularFile(files.intendedOut);
		final boolean truncatedExists=files.truncatedOut!=null && Files.isRegularFile(files.truncatedOut);
		final Path observedPath=intendedExists ? files.intendedOut :
			(truncatedExists ? files.truncatedOut : files.intendedOut);
		final ArrayList<FastaRecord> records=Files.isRegularFile(observedPath) ?
			parseFasta(observedPath) : new ArrayList<FastaRecord>();
		final long observedCounter=parseCounter(errText);
		final boolean sourceMatch=row.source.matches(records, observedCounter,
			intendedExists, truncatedExists);
		final boolean correctedMatch=row.corrected.matches(records, observedCounter,
			intendedExists, truncatedExists);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean stderrMatches=row.corrected.stderrMatches(errText);
		final boolean stdoutEmpty=stdout.length==0;
		final boolean unchanged=inputsUnchanged();
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;

		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+
			(baseline ? "baseline" : "fixed")+"\t"+exit+"\t0\t"+
			(observedCounter==NO_COUNTER ? "" : Long.toString(observedCounter))+"\t"+
			row.source.counter+"\t"+row.corrected.counter+"\t"+
			clean(out.relativize(observedPath).toString())+"\t"+clean(serialize(records))+"\t"+
			clean(row.source.serialize())+"\t"+clean(row.corrected.serialize())+"\t"+
			intendedExists+"\t"+truncatedExists+"\t"+row.source.intendedExists+
			"\t"+row.source.truncatedExists+"\t"+row.corrected.intendedExists+
			"\t"+row.corrected.truncatedExists+"\t"+sourceMatch+"\t"+
			correctedMatch+"\t"+stderrMatches+"\t"+stdoutEmpty+"\t"+
			unchanged+"\t"+phaseMatch+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);

		check(exit==0, row.name+": child exit was "+exit);
		check(phaseMatch, row.name+": FASTA/counter/destination mismatch");
		check(stderrMatches, row.name+": stderr did not match oracle: "+errText);
		check(stdoutEmpty, row.name+": bbmask.sh wrote to stdout");
		check(command.contains("-ea") && command.contains("t=1") &&
			command.contains("simd=f") && command.contains("mr=f") &&
			command.contains("me=f") && command.contains("minscaf=1") &&
			command.contains("sampad=0") && !command.contains("split=t"),
			row.name+": child argv lost a frozen CLI pin");
		check(unchanged, row.name+": a child modified literal input bytes");
	}

	/** Writes the literal references and SAMs needed for one child case. */
	private CaseFiles writeInputs(final Path caseDir, final Row row) throws IOException{
		Path intendedRef=caseDir.resolve("ref.fa");
		final Path intendedSam=caseDir.resolve("read.sam");
		final Path intendedOut, truncatedOut;
		if(row.kind==Kind.REFERENCE_EQUALS){
			final Path refDir=caseDir.resolve("refs.fa=dir");
			Files.createDirectory(refDir);
			intendedRef=refDir.resolve("ref.fa");
			write(caseDir.resolve("refs.fa"), ">ref\nTTTT\n");
			write(intendedRef, FASTA);
		}else{
			write(intendedRef, FASTA);
		}
		if(row.kind==Kind.SAM_EQUALS){
			final Path samDir=caseDir.resolve("reads.sam=dir");
			Files.createDirectory(samDir);
			final Path withEquals=samDir.resolve("read.sam");
			write(caseDir.resolve("reads.sam"), sam(POS4));
			write(withEquals, sam(POS1));
			return new CaseFiles(intendedRef, withEquals, caseDir.resolve("result.fa"), null);
		}
		if(row.kind==Kind.OUTPUT_EQUALS){
			final Path outDir=caseDir.resolve("out.fa=dir");
			Files.createDirectory(outDir);
			intendedOut=outDir.resolve("result.fa");
			truncatedOut=caseDir.resolve("out.fa");
		}else{
			intendedOut=caseDir.resolve("result.fa");
			truncatedOut=null;
		}
		return new CaseFiles(intendedRef, intendedSam, intendedOut, truncatedOut);
	}

	/** Launch arguments preserving each case's intentional '=' placement. */
	private ArrayList<String> baseCommand(final Row row, final CaseFiles files){
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		if(row.kind==Kind.POSITIONAL){command.add(files.ref.toString());}
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("t=1");
		command.add("simd=f");
		command.add("mr=f");
		command.add("me=f");
		command.add("minscaf=1");
		command.add("sampad=0");
		command.add("ow=t");
		if(row.kind!=Kind.POSITIONAL){command.add("in="+files.ref);}
		if(row.kind==Kind.SAM_EQUALS){command.add("sam="+files.sam);}
		command.add("out="+files.intendedOut);
		return command;
	}

	/** Writes a literal input and records it for after-child preservation checks. */
	private void write(final Path path, final String text) throws IOException{
		final byte[] bytes=text.getBytes(StandardCharsets.US_ASCII);
		Files.write(path, bytes, StandardOpenOption.CREATE_NEW);
		original.put(path, bytes);
	}

	/** Expected source and corrected row for one case. */
	private static String expected(final Row row){
		return "phase\tcounter\tfasta\tintendedExists\ttruncatedExists\n"+
			"source\t"+row.source.counter+"\t"+clean(row.source.serialize())+
			"\t"+row.source.intendedExists+"\t"+row.source.truncatedExists+"\n"+
			"corrected\t"+row.corrected.counter+"\t"+clean(row.corrected.serialize())+
			"\t"+row.corrected.intendedExists+"\t"+row.corrected.truncatedExists+"\n";
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

	/** Joins argv for stable evidence. */
	private static String join(final List<String> command){
		final StringBuilder sb=new StringBuilder();
		for(final String s : command){
			if(sb.length()>0){sb.append(' ');}
			sb.append(s);
		}
		return sb.append('\n').toString();
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

	/** Expected output for one source or corrected phase. */
	private static final class Expect{

		Expect(final long counter_, final boolean intendedExists_,
				final boolean truncatedExists_, final FastaRecord... records_){
			counter=counter_; intendedExists=intendedExists_;
			truncatedExists=truncatedExists_; records=records_;
		}

		boolean matches(final ArrayList<FastaRecord> observed, final long observedCounter,
				final boolean observedIntended, final boolean observedTruncated){
			return counter==observedCounter && intendedExists==observedIntended &&
				truncatedExists==observedTruncated &&
				Arrays.equals(records, observed.toArray(new FastaRecord[observed.size()]));
		}

		boolean changed(final Expect other){
			return counter!=other.counter || intendedExists!=other.intendedExists ||
				truncatedExists!=other.truncatedExists || !Arrays.equals(records, other.records);
		}

		boolean stderrMatches(final String errText){
			return errText.contains("Picked up JAVA_TOOL_OPTIONS") &&
				!errText.contains("Exception") && !errText.contains("AssertionError") &&
				!errText.contains("Unknown parameter") && !errText.contains("Error reading");
		}

		String serialize(){return BBMaskFirstEqualsTest.serialize(Arrays.asList(records));}

		final long counter;
		final boolean intendedExists, truncatedExists;
		final FastaRecord[] records;
	}

	/** Paths used by one child. */
	private static final class CaseFiles{

		CaseFiles(final Path ref_, final Path sam_, final Path intendedOut_,
				final Path truncatedOut_){
			ref=ref_; sam=sam_; intendedOut=intendedOut_; truncatedOut=truncatedOut_;
		}

		final Path ref, sam, intendedOut, truncatedOut;
	}

	/** One named or positional child-process case. */
	private static final class Row{

		Row(final String name_, final Kind kind_, final Expect source_,
				final Expect corrected_){
			name=name_; kind=kind_; source=source_; corrected=corrected_;
		}

		final String name;
		final Kind kind;
		final Expect source, corrected;
	}

	private enum Kind { REFERENCE_EQUALS, SAM_EQUALS, OUTPUT_EQUALS, NAMED, POSITIONAL }

	private final Path out, repo, launcher;
	private final Map<Path, byte[]> original=new LinkedHashMap<Path, byte[]>();
	private int cases, children, changed, unchanged, sourceMatches, correctedMatches, checks, failures;
	private final boolean baseline;

	private static final long NO_COUNTER=Long.MIN_VALUE;
	private static final Pattern COUNTER_PATTERN=Pattern.compile("Total Bases Masked:\\s+(\\d+)/");
	private static final String FASTA=">ref\nACGT\n";
	private static final String SAM_HEADER="@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:4\n";
	private static final FastaRecord ACGT=new FastaRecord("ref", "ACGT");
	private static final FastaRecord TTTT=new FastaRecord("ref", "TTTT");
	private static final FastaRecord ACGN=new FastaRecord("ref", "ACGN");
	private static final FastaRecord NCGT=new FastaRecord("ref", "NCGT");
	private static final Expect EXPECT_ACGT=new Expect(0, true, false, ACGT);
	private static final SamRecord POS1=new SamRecord("pos1", 1);
	private static final SamRecord POS4=new SamRecord("pos4", 4);

	private static final Row[] ROWS={
		new Row("reference-equals", Kind.REFERENCE_EQUALS,
			new Expect(0, true, false, TTTT), EXPECT_ACGT),
		new Row("sam-equals", Kind.SAM_EQUALS,
			new Expect(1, true, false, ACGN), new Expect(1, true, false, NCGT)),
		new Row("output-equals", Kind.OUTPUT_EQUALS,
			new Expect(0, false, true, ACGT), EXPECT_ACGT),
		new Row("named-control", Kind.NAMED, EXPECT_ACGT, EXPECT_ACGT),
		new Row("positional-control", Kind.POSITIONAL, EXPECT_ACGT, EXPECT_ACGT),
	};

	/** One SAM alignment to the four-base ref. */
	private static final class SamRecord{

		SamRecord(final String read_, final int pos_){read=read_; pos=pos_;}

		String text(){return read+"\t0\tref\t"+pos+"\t60\t1M\t*\t0\t0\tA\tI\n";}

		final String read;
		final int pos;
	}

	/** One literal SAM file. */
	private static String sam(final SamRecord record){
		return SAM_HEADER+record.text();
	}
}
