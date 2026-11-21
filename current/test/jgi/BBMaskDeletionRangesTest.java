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
import structures.IntList;

/**
 * Finite BBMask deletion-coverage range checks for raw half-open fillRanges()
 * pairs and valid 2M2D2M CLI threshold cases.
 * @author Jean
 */
public final class BBMaskDeletionRangesTest{

	/** Requires phase=baseline|fixed, a new outdir and the BBTools root. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown BBMask deletion-ranges option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/bbmaskdeletionrangestest.sh with assertions enabled.");
		final Path out=Paths.get(directory).toAbsolutePath().normalize();
		final Path repo=Paths.get(rootPath).toAbsolutePath().normalize();
		Files.createDirectory(out);
		final BBMaskDeletionRangesTest test=new BBMaskDeletionRangesTest(out, repo,
			phase.equals("baseline"));
		test.writeInputs();
		for(final DirectRow row : DIRECT_ROWS){test.runDirect(row);}
		for(final CliRow row : CLI_ROWS){test.runCli(row);}
		test.check(test.directCases==14, "Expected fourteen direct fillRanges cases: "+test.directCases);
		test.check(test.cliCases==4, "Expected four BBMask CLI cases: "+test.cliCases);
		test.check(test.children==4, "Expected four fresh bbmask.sh children: "+test.children);
		test.check(test.changedDirect==5, "Expected five changed direct rows: "+test.changedDirect);
		test.check(test.changedCli==2, "Expected two changed CLI rows: "+test.changedCli);
		test.check(test.unchanged==11, "Expected eleven unchanged controls: "+test.unchanged);
		final int expectedSourceMatches=(test.baseline ? 18 : 11);
		final int expectedCorrectedMatches=(test.baseline ? 11 : 18);
		test.check(test.sourceMatches==expectedSourceMatches,
			"Expected "+expectedSourceMatches+" source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==expectedCorrectedMatches,
			"Expected "+expectedCorrectedMatches+" corrected matches: "+test.correctedMatches);
		final String summary="BBMASK_DELETION_RANGES_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+(test.directCases+test.cliCases)+
			" direct="+test.directCases+" children="+test.children+
			" changedDirect="+test.changedDirect+" changedCli="+test.changedCli+
			" unchanged="+test.unchanged+" sourceMatches="+test.sourceMatches+
			" correctedMatches="+test.correctedMatches+" checks="+test.checks+
			" failures="+test.failures+"\n";
		Files.write(out.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+out);}
		System.out.print(summary);
	}

	/** Creates evidence tables and binds the script under test. */
	private BBMaskDeletionRangesTest(final Path out_, final Path repo_, final boolean baseline_)
			throws IOException{
		out=out_; repo=repo_; inputs=out.resolve("inputs"); baseline=baseline_;
		launcher=repo.resolve("bbmask.sh");
		require(Files.isRegularFile(launcher), "Missing bbmask.sh at "+launcher);
		Files.createDirectory(inputs);
		Files.write(out.resolve("direct.tsv"), ("case\tphase\tincludeDels\tmatch\tstart"+
			"\tstop\trawPairs\tsourceRawPairs\tcorrectedRawPairs\tcovered"+
			"\tsourceCovered\tcorrectedCovered\tsourceMatches\tcorrectedMatches\tmatchPhase\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		Files.write(out.resolve("cli.tsv"), ("case\tphase\tdelcov\texit\texpectedExit"+
			"\tcounter\tsourceCounter\tcorrectedCounter\tobservedFasta"+
			"\tsourceFasta\tcorrectedFasta\tsourceMatches\tcorrectedMatches"+
			"\tstderrMatches\tstdoutEmpty\tinputsUnchanged\tmatchPhase\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes all literal FASTA and SAM inputs before any child runs. */
	private void writeInputs() throws IOException{
		for(final CliRow row : CLI_ROWS){
			final Path caseInputs=inputs.resolve(row.name);
			Files.createDirectory(caseInputs);
			final Path fasta=caseInputs.resolve("ref.fa");
			final Path sam=caseInputs.resolve("read.sam");
			final byte[] fastaBytes=FASTA.getBytes(StandardCharsets.US_ASCII);
			final byte[] samBytes=SAM.getBytes(StandardCharsets.US_ASCII);
			Files.write(fasta, fastaBytes, StandardOpenOption.CREATE_NEW);
			Files.write(sam, samBytes, StandardOpenOption.CREATE_NEW);
			original.put(fasta, fastaBytes);
			original.put(sam, samBytes);
		}
	}

	/** Calls public BBMask.fillRanges() and checks raw ranges plus covered positions. */
	private void runDirect(final DirectRow row) throws IOException{
		directCases++;
		if(row.changed()){changedDirect++;}
		else{unchanged++;}
		final IntList ranges=new IntList(8);
		BBMask.fillRanges(row.match.getBytes(StandardCharsets.US_ASCII), row.start,
			row.stop, ranges, row.includeDels);
		final int[] observed=toArray(ranges);
		final String observedCovered=covered(observed);
		final boolean sourceMatch=Arrays.equals(observed, row.sourcePairs) &&
			observedCovered.equals(row.sourceCovered);
		final boolean correctedMatch=Arrays.equals(observed, row.correctedPairs) &&
			observedCovered.equals(row.correctedCovered);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;
		Files.write(out.resolve("direct.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+row.includeDels+"\t"+row.match+"\t"+row.start+"\t"+row.stop+
			"\t"+pairs(observed)+"\t"+pairs(row.sourcePairs)+"\t"+
			pairs(row.correctedPairs)+"\t"+observedCovered+"\t"+
			row.sourceCovered+"\t"+row.correctedCovered+"\t"+
			sourceMatch+"\t"+correctedMatch+"\t"+phaseMatch+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check((observed.length&1)==0, row.name+": odd number of raw range endpoints");
		check(phaseMatch, row.name+": raw pairs were "+pairs(observed)+
			", covered positions were "+observedCovered);
	}

	/** Launches one valid bbmask.sh child and checks the exact threshold output. */
	private void runCli(final CliRow row) throws Exception{
		cliCases++;
		children++;
		if(row.changed()){changedCli++;}
		else{unchanged++;}
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
		command.add("delcov="+(row.delcov ? "t" : "f"));
		if(row.mincov>-1){command.add("mincov="+row.mincov);}
		if(row.maxcov>-1){command.add("maxcov="+row.maxcov);}
		command.add("in="+input(row.name, "ref.fa"));
		command.add("out="+output);
		command.add("sam="+input(row.name, "read.sam"));
		Files.write(caseDir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		final Path stdoutPath=caseDir.resolve("stdout.txt");
		final Path stderrPath=caseDir.resolve("stderr.txt");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdoutPath.toFile()).redirectError(stderrPath.toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		Files.write(caseDir.resolve("exit.txt"), (Integer.toString(exit)+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		final byte[] stdout=Files.readAllBytes(stdoutPath);
		final String errText=new String(Files.readAllBytes(stderrPath), StandardCharsets.UTF_8);
		final ArrayList<FastaRecord> records=Files.isRegularFile(output) ?
			parseFasta(output) : new ArrayList<FastaRecord>();
		final long counter=parseCounter(errText);
		final boolean sourceMatch=row.sourceExpect.matches(records, counter);
		final boolean correctedMatch=row.correctedExpect.matches(records, counter);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean stderrMatches=stderrMatches(errText);
		final boolean stdoutEmpty=stdout.length==0;
		final boolean unchanged=inputsUnchanged();
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;
		Files.write(out.resolve("cli.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+row.delcov+"\t"+exit+"\t0\t"+counter+"\t"+
			row.sourceExpect.counter+"\t"+row.correctedExpect.counter+"\t"+
			clean(serialize(records))+"\t"+clean(row.sourceExpect.serialize())+"\t"+
			clean(row.correctedExpect.serialize())+"\t"+sourceMatch+"\t"+
			correctedMatch+"\t"+stderrMatches+"\t"+stdoutEmpty+"\t"+unchanged+
			"\t"+phaseMatch+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);
		check(exit==0, row.name+": child exit was "+exit);
		check(phaseMatch, row.name+": FASTA/counter mismatch: "+serialize(records)+" "+counter);
		check(stderrMatches, row.name+": stderr did not match oracle: "+errText);
		check(stdoutEmpty, row.name+": bbmask.sh wrote to stdout");
		check(unchanged, row.name+": a child modified literal input bytes");
		check(command.contains("delcov="+(row.delcov ? "t" : "f")) &&
			command.contains("-ea") && command.contains("t=1") &&
			command.contains("simd=f") && command.contains("minscaf=1") &&
			!command.contains("split=t"),
			row.name+": child argv lost a frozen CLI pin");
	}

	/** Writes both the baseline source oracle and corrected oracle for one child. */
	private static void writeExpected(final Path caseDir, final CliRow row) throws IOException{
		final String text="delcov\tsourceCounter\tsourceFasta\tcorrectedCounter\tcorrectedFasta\n"+
			row.delcov+"\t"+row.sourceExpect.counter+"\t"+
			clean(row.sourceExpect.serialize())+"\t"+row.correctedExpect.counter+
			"\t"+clean(row.correctedExpect.serialize())+"\n";
		Files.write(caseDir.resolve("expected.tsv"), text.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	private Path input(final String rowName, final String fileName){return inputs.resolve(rowName).resolve(fileName);}

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

	/** Records a failed invariant without aborting subsequent checks. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println(message);}
	}

	/** Fails fast for invalid harness setup. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static int[] toArray(final IntList list){
		final int[] array=new int[list.size];
		for(int i=0; i<array.length; i++){array[i]=list.get(i);}
		return array;
	}

	private static String pairs(final int[] array){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<array.length; i+=2){
			if(i>0){sb.append(',');}
			sb.append('[').append(array[i]).append(',').append(array[i+1]).append(')');
		}
		return sb.toString();
	}

	private static String covered(final int[] array){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<array.length; i+=2){
			for(int j=array[i]; j<array[i+1]; j++){
				if(sb.length()>0){sb.append(',');}
				sb.append(j);
			}
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

	/** Exact CLI oracle for one phase. */
	private static final class Expect{
		Expect(final long counter_, final FastaRecord record_){counter=counter_; record=record_;}
		boolean matches(final ArrayList<FastaRecord> records, final long observedCounter){
			return counter==observedCounter && records.size()==1 && record.equals(records.get(0));
		}
		String serialize(){return BBMaskDeletionRangesTest.serialize(Arrays.asList(record));}
		final long counter;
		final FastaRecord record;
	}

	/** One direct fillRanges case. */
	private static final class DirectRow{
		DirectRow(final String name_, final String match_, final int start_, final int stop_,
				final boolean includeDels_, final int[] sourcePairs_, final int[] correctedPairs_,
				final String sourceCovered_, final String correctedCovered_){
			name=name_; match=match_; start=start_; stop=stop_; includeDels=includeDels_;
			sourcePairs=sourcePairs_; correctedPairs=correctedPairs_;
			sourceCovered=sourceCovered_; correctedCovered=correctedCovered_;
		}
		boolean changed(){return !Arrays.equals(sourcePairs, correctedPairs) ||
			!sourceCovered.equals(correctedCovered);}
		final String name, match;
		final int start, stop;
		final boolean includeDels;
		final int[] sourcePairs, correctedPairs;
		final String sourceCovered, correctedCovered;
	}

	/** One valid 2M2D2M child-process case. */
	private static final class CliRow{
		CliRow(final String name_, final boolean delcov_, final int mincov_, final int maxcov_,
				final Expect sourceExpect_, final Expect correctedExpect_){
			name=name_; delcov=delcov_; mincov=mincov_; maxcov=maxcov_;
			sourceExpect=sourceExpect_; correctedExpect=correctedExpect_;
		}
		boolean changed(){return !sourceExpect.record.equals(correctedExpect.record) ||
			sourceExpect.counter!=correctedExpect.counter;}
		final String name;
		final boolean delcov;
		final int mincov, maxcov;
		final Expect sourceExpect, correctedExpect;
	}

	private final Path out, repo, inputs, launcher;
	private final Map<Path, byte[]> original=new LinkedHashMap<Path, byte[]>();
	private int directCases, cliCases, children, changedDirect, changedCli, unchanged;
	private int sourceMatches, correctedMatches, checks, failures;
	private final boolean baseline;

	private static final Pattern COUNTER_PATTERN=Pattern.compile("Total Bases Masked:\\s+(\\d+)/");
	private static final String FASTA=">ref\nACGTAC\n";
	private static final String SAM="@HD\tVN:1.6\tSO:unsorted\n"+
		"@SQ\tSN:ref\tLN:6\n"+
		"read\t0\tref\t1\t60\t2M2D2M\t*\t0\t0\tACAC\tIIII\n";
	private static final FastaRecord SOURCE_MIN=new FastaRecord("ref", "ACNNAC");
	private static final FastaRecord SOURCE_MAX=new FastaRecord("ref", "NNGTNN");
	private static final FastaRecord FIXED_MIN=new FastaRecord("ref", "ACGTAC");
	private static final FastaRecord FIXED_MAX=new FastaRecord("ref", "NNNNNN");

	private static final DirectRow[] DIRECT_ROWS={
		new DirectRow("internal-true", "mmDDmm", 0, 6, true,
			new int[] {0, 2, 4, 6}, new int[] {0, 2, 2, 4, 4, 6},
			"0,1,4,5", "0,1,2,3,4,5"),
		new DirectRow("internal-false", "mmDDmm", 0, 6, false,
			new int[] {0, 2, 4, 6}, new int[] {0, 2, 4, 6},
			"0,1,4,5", "0,1,4,5"),
		new DirectRow("terminal-true", "mmDD", 0, 4, true,
			new int[] {0, 2, 2, 4}, new int[] {0, 2, 2, 4},
			"0,1,2,3", "0,1,2,3"),
		new DirectRow("terminal-false", "mmDD", 0, 4, false,
			new int[] {0, 2}, new int[] {0, 2},
			"0,1", "0,1"),
		new DirectRow("nodeletion-true", "mmmm", 0, 4, true,
			new int[] {0, 4}, new int[] {0, 4},
			"0,1,2,3", "0,1,2,3"),
		new DirectRow("nodeletion-false", "mmmm", 0, 4, false,
			new int[] {0, 4}, new int[] {0, 4},
			"0,1,2,3", "0,1,2,3"),
		new DirectRow("leading-true", "DDmm", 0, 4, true,
			new int[] {2, 4}, new int[] {0, 2, 2, 4},
			"2,3", "0,1,2,3"),
		new DirectRow("leading-false", "DDmm", 0, 4, false,
			new int[] {2, 4}, new int[] {2, 4},
			"2,3", "2,3"),
		new DirectRow("clip-true", "CCDDmm", 0, 6, true,
			new int[] {0, 2, 4, 6}, new int[] {2, 4, 4, 6},
			"0,1,4,5", "2,3,4,5"),
		new DirectRow("clip-false", "CCDDmm", 0, 6, false,
			new int[] {4, 6}, new int[] {4, 6},
			"4,5", "4,5"),
		new DirectRow("insertion-true", "IIDDmm", 0, 4, true,
			new int[] {0, 0, 2, 4}, new int[] {0, 2, 2, 4},
			"2,3", "0,1,2,3"),
		new DirectRow("insertion-false", "IIDDmm", 0, 4, false,
			new int[] {2, 4}, new int[] {2, 4},
			"2,3", "2,3"),
		new DirectRow("offset-internal-true", "mmDDmm", 10, 16, true,
			new int[] {10, 12, 14, 16}, new int[] {10, 12, 12, 14, 14, 16},
			"10,11,14,15", "10,11,12,13,14,15"),
		new DirectRow("offset-internal-false", "mmDDmm", 10, 16, false,
			new int[] {10, 12, 14, 16}, new int[] {10, 12, 14, 16},
			"10,11,14,15", "10,11,14,15")
	};

	private static final CliRow[] CLI_ROWS={
		new CliRow("mincov-delcov-true", true, 1, -1,
			new Expect(2, SOURCE_MIN), new Expect(0, FIXED_MIN)),
		new CliRow("mincov-delcov-false", false, 1, -1,
			new Expect(2, SOURCE_MIN), new Expect(2, SOURCE_MIN)),
		new CliRow("maxcov-delcov-true", true, -1, 0,
			new Expect(4, SOURCE_MAX), new Expect(6, FIXED_MAX)),
		new CliRow("maxcov-delcov-false", false, -1, 0,
			new Expect(4, SOURCE_MAX), new Expect(4, SOURCE_MAX))
	};
}
