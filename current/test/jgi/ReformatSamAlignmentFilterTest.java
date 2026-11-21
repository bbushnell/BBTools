package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import parse.Parse;

/** Finite SAM breaklength rejection cases for original-alignment filters.
 * @author Jean
 */
public final class ReformatSamAlignmentFilterTest{

	/** Runs public ReformatStreamer children against source and corrected oracles. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, rootPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){rootPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown Reformat alignment-filter option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(rootPath!=null && !rootPath.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/reformatsamalignmentfiltertest.sh with assertions enabled.");
		final ReformatSamAlignmentFilterTest test=new ReformatSamAlignmentFilterTest(
			Paths.get(directory).toAbsolutePath().normalize(),
			Paths.get(rootPath).toAbsolutePath().normalize(), phase.equals("baseline"));
		for(final Case item : CASES){
			if(test.baseline && item.fixedOnly){continue;}
			test.run(item);
		}
		final int expectedCases=CASES.length-(test.baseline ? FIXED_ONLY_CASES : 0);
		test.check(test.cases==expectedCases, "Expected "+expectedCases+" cases, saw "+test.cases);
		test.check(test.children==expectedCases, "Expected "+expectedCases+" children, saw "+test.children);
		test.check(test.changed==(test.baseline ? 21 : 23), "Expected 21/23 repaired rows, saw "+test.changed);
		test.check(test.controls==2, "Expected two unchanged controls, saw "+test.controls);
		test.check(test.sourceMatches==(test.baseline ? 23 : 2), "Unexpected source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==(test.baseline ? 2 : 25), "Unexpected corrected matches: "+test.correctedMatches);
		final String summary="REFORMAT_SAM_ALIGNMENT_FILTER_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" changed="+test.changed+" controls="+test.controls+
			" sourceMatches="+test.sourceMatches+" correctedMatches="+test.correctedMatches+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+test.out);}
		System.out.print(summary);
	}

	/** Binds the launcher and initializes the evidence table. */
	private ReformatSamAlignmentFilterTest(final Path out_, final Path repo_, final boolean baseline_)
			throws Exception{
		out=out_; repo=repo_; baseline=baseline_;
		launcher=repo.resolve("reformat.sh");
		require(Files.isRegularFile(launcher), "Missing reformat.sh at "+launcher);
		Files.createDirectory(out);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\tdataRows\toutputExists"+
			"\tinputReads\tinputBases\toutputReads\toutputBases"+
			"\tsourceMatches\tcorrectedMatches\tstderrMatches\tstdoutEmpty\tinputsUnchanged\tmatchPhase\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Launches one fresh public reformat.sh child and checks exact rows or fixed preflight rejection. */
	private void run(final Case item) throws Exception{
		cases++;
		children++;
		if(item.source.changed(item.corrected)){changed++;}else{controls++;}
		final Path dir=out.resolve(item.name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.sam"), output=dir.resolve("output"+item.outputExt);
		Files.write(input, item.sam().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final byte[] original=Files.readAllBytes(input);
		check(!Files.exists(output), item.name+": output path existed before child");
		final ArrayList<String> command=command(item, input, output);
		Files.write(dir.resolve("argv.txt"), join(command).getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		Files.write(dir.resolve("expected.tsv"), expected(item).getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		final Path stdout=dir.resolve("stdout.txt"), stderr=dir.resolve("stderr.txt");
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(stdout.toFile()).redirectError(stderr.toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		Files.write(dir.resolve("exit.txt"), (exit+"\n").getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		final byte[] stdoutBytes=Files.readAllBytes(stdout);
		final String stderrText=new String(Files.readAllBytes(stderr), StandardCharsets.UTF_8);
		final ArrayList<String> rows=readDataRows(output);
		final boolean outputExists=Files.isRegularFile(output);
		final Counts counts=Counts.parse(stderrText);
		final boolean sourceMatch=item.source.matches(exit, outputExists, rows, stderrText, counts);
		final boolean correctedMatch=item.corrected.matches(exit, outputExists, rows, stderrText, counts);
		if(sourceMatch){sourceMatches++;}
		if(correctedMatch){correctedMatches++;}
		final boolean stderrMatches=(baseline ? item.source : item.corrected).stderrMatches(stderrText);
		final boolean stdoutEmpty=stdoutBytes.length==0;
		final boolean unchanged=Arrays.equals(original, Files.readAllBytes(input));
		final boolean phaseMatch=baseline ? sourceMatch : correctedMatch;
		Files.write(out.resolve("cases.tsv"), (item.name+"\t"+(baseline ? "baseline" : "fixed")+
			"\t"+exit+"\t"+rows.size()+"\t"+outputExists+"\t"+counts.tsv()+"\t"+sourceMatch+
			"\t"+correctedMatch+"\t"+stderrMatches+"\t"+stdoutEmpty+
			"\t"+unchanged+"\t"+phaseMatch+"\n").getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.APPEND);
		check(phaseMatch, item.name+": wrong rows/rejection; stderr="+stderrText);
		check(stderrMatches, item.name+": stderr mismatch: "+stderrText);
		check(stdoutEmpty, item.name+": Reformat wrote to stdout");
		check(unchanged, item.name+": input bytes changed");
		check(command.contains("breaklength="+(item.breakLength ? 3 : 0)) &&
			command.contains("int=f") && command.contains("t=1") && command.contains("wt=1") &&
			command.contains("tin=1") && command.contains("tout=1") && command.contains("simd=f") &&
			command.contains("qin=33") && command.contains("qout=33") &&
			!command.contains("sample=1") && !command.contains("samplebases=1"),
			item.name+": child argv lost a safe SAM breaklength pin");
	}

	/** Builds the exact public launcher command for one row. */
	private ArrayList<String> command(final Case item, final Path input, final Path output){
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher.toString());
		command.add("-Xmx128m");
		command.add("-ea");
		command.add("in="+input);
		command.add("out="+output);
		command.add("t=1");
		command.add("wt=1");
		command.add("tin=1");
		command.add("tout=1");
		command.add("int=f");
		command.add("ow=f");
		command.add("app=f");
		command.add("simd=f");
		command.add("minlength=1");
		command.add("breaklength="+(item.breakLength ? 3 : 0));
		command.add("qin=33");
		command.add("qout=33");
		for(final String flag : item.flags){command.add(flag);}
		return command;
	}

	/** Reads non-header SAM rows; absent fixed-rejection outputs read as zero rows. */
	private static ArrayList<String> readDataRows(final Path output) throws Exception{
		final ArrayList<String> rows=new ArrayList<String>();
		if(!Files.isRegularFile(output)){return rows;}
		for(final String line : Files.readAllLines(output, StandardCharsets.US_ASCII)){
			if(!line.isEmpty() && line.charAt(0)!='@'){rows.add(line);}
		}
		return rows;
	}

	/** Expected rows for evidence transcripts. */
	private static String expected(final Case item){
		return "phase\tresult\trows\nsource\t"+item.source+"\ncorrected\t"+item.corrected+"\n";
	}

	/** Shell-escaped argv transcript for the evidence directory. */
	private static String join(final ArrayList<String> command){
		final StringBuilder builder=new StringBuilder();
		for(final String arg : command){
			if(builder.length()>0){builder.append(' ');}
			builder.append('\'').append(arg.replace("'", "'\\''")).append('\'');
		}
		return builder.append('\n').toString();
	}

	/** Records one assertion while retaining later row evidence. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	/** Hard construction-time requirement. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static Case reject(final String name, final Sam sam, final String flag){
		return new Case(name, true, sam, new String[]{flag}, sam.rows(), Reject.INSTANCE);
	}

	private static Case fixedOnlyReject(final String name, final String outputExt){
		return new Case(name, true, MAPPED, new String[]{"minmapq=60"}, NotRun.INSTANCE,
			Reject.INSTANCE, outputExt, true);
	}

	private static Case control(final String name, final Sam sam, final String[] flags){
		return new Case(name, false, sam, flags, sam.rows(), sam.rows());
	}

	/** One exact SAM input row. */
	private static final class Sam{
		Sam(final String name_, final int flag_, final String rname_, final int pos_, final int mapq_,
				final String cigar_, final String rnext_, final int pnext_, final int tlen_,
				final String bases_, final String qual_){
			name=name_;
			row=name_+'\t'+flag_+'\t'+rname_+'\t'+pos_+'\t'+mapq_+'\t'+cigar_+
				'\t'+rnext_+'\t'+pnext_+'\t'+tlen_+'\t'+bases_+'\t'+qual_;
		}
		String input(){return "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\n"+row+'\n';}
		Rows rows(){return new Rows(row);}
		final String name;
		final String row;
	}

	/** One public launcher child row. */
	private static final class Case{
		Case(final String name_, final boolean breakLength_, final Sam sam_, final String[] flags_,
				final Expect source_, final Expect corrected_){
			this(name_, breakLength_, sam_, flags_, source_, corrected_, ".sam", false);
		}
		Case(final String name_, final boolean breakLength_, final Sam sam_, final String[] flags_,
				final Expect source_, final Expect corrected_, final String outputExt_,
				final boolean fixedOnly_){
			name=name_;
			breakLength=breakLength_;
			sam=sam_;
			flags=flags_;
			source=source_;
			corrected=corrected_;
			outputExt=outputExt_;
			fixedOnly=fixedOnly_;
		}
		String sam(){return sam.input();}
		final String name;
		final boolean breakLength;
		final Sam sam;
		final String[] flags;
		final Expect source;
		final Expect corrected;
		final String outputExt;
		final boolean fixedOnly;
	}

	/** Expected result for one phase. */
	private interface Expect{
		boolean matches(int exit, boolean outputExists, ArrayList<String> rows, String stderr,
			Counts counts);
		boolean stderrMatches(String stderr);
		boolean changed(Expect other);
	}

	/** Exact successful data rows. */
	private static final class Rows implements Expect{
		Rows(final String... rows_){this(1, 8, rows_);}
		Rows(final int inputReads_, final int inputBases_, final String... rows_){
			rows=rows_;
			int bases=0;
			for(final String row : rows){
				final String[] fields=row.split("\t");
				require(fields.length==11, "Bad SAM oracle row: "+row);
				bases+=fields[9].length();
			}
			counts=new Counts(inputReads_, inputBases_, rows.length, bases);
		}
		public boolean matches(final int exit, final boolean outputExists, final ArrayList<String> observed,
				final String stderr, final Counts observedCounts){
			return exit==0 && outputExists && observed.equals(Arrays.asList(rows)) &&
				counts.equals(observedCounts) && stderrMatches(stderr);
		}
		public boolean stderrMatches(final String stderr){
			return !stderr.contains("Exception in thread") && !stderr.contains("Unknown parameter") &&
				!stderr.contains("java.lang.") && !stderr.contains("ERROR") && !stderr.contains("Error");
		}
		public boolean changed(final Expect other){return other!=this && !toString().equals(other.toString());}
		@Override
		public String toString(){return "success\t"+Arrays.toString(rows);}
		final String[] rows;
		final Counts counts;
	}

	/** Future preflight rejection before the output writer creates the SAM destination. */
	private static final class Reject implements Expect{
		public boolean matches(final int exit, final boolean outputExists, final ArrayList<String> rows,
				final String stderr, final Counts counts){
			return exit==1 && !outputExists && rows.isEmpty() && counts.isAbsent() &&
				!Counts.hasRawCounter(stderr) && stderrMatches(stderr);
		}
		public boolean stderrMatches(final String stderr){
			return stderr.contains("Exception in thread \"main\" java.lang.IllegalArgumentException: "+
				"breaklength with SAM/BAM input is incompatible with alignment-dependent filters.\n") &&
				stderr.contains("\tat jgi.ReformatStreamer.<init>(ReformatStreamer.java:");
		}
		public boolean changed(final Expect other){return other!=this;}
		@Override
		public String toString(){return "fixed-preflight-rejection";}
		static final Reject INSTANCE=new Reject();
	}

	/** Source path deliberately skipped because it is unsafe to replay. */
	private static final class NotRun implements Expect{
		public boolean matches(final int exit, final boolean outputExists, final ArrayList<String> rows,
				final String stderr, final Counts counts){
			return false;
		}
		public boolean stderrMatches(final String stderr){return false;}
		public boolean changed(final Expect other){return other!=this;}
		@Override
		public String toString(){return "not run; unsafe baseline path";}
		static final NotRun INSTANCE=new NotRun();
	}

	/** Parsed ReformatStreamer completion counters from stderr. */
	private static final class Counts{
		Counts(final long inputReads_, final long inputBases_, final long outputReads_,
				final long outputBases_){
			inputReads=inputReads_;
			inputBases=inputBases_;
			outputReads=outputReads_;
			outputBases=outputBases_;
		}
		static Counts parse(final String stderr){
			final long[] input=parseLine(stderr, "Input");
			final long[] output=parseLine(stderr, "Output");
			return new Counts(input[0], input[1], output[0], output[1]);
		}
		static boolean hasRawCounter(final String stderr){
			return Pattern.compile("(?m)^(?:Input|Output):").matcher(stderr).find();
		}
		private static long[] parseLine(final String stderr, final String label){
			final Matcher matcher=Pattern.compile("(?m)^"+label+
				":\\s*([0-9][0-9,]*)\\s+reads?(?:\\s+\\([0-9]+\\.[0-9][0-9]%\\))?"+
				"\\s+([0-9][0-9,]*)\\s+bases?(?:\\s+\\([0-9]+\\.[0-9][0-9]%\\))?\\s*$")
				.matcher(stderr);
			if(!matcher.find()){return ABSENT;}
			final long[] counts=new long[] {Long.parseLong(matcher.group(1).replace(",", "")),
				Long.parseLong(matcher.group(2).replace(",", ""))};
			if(matcher.find()){return ABSENT;}
			return counts;
		}
		boolean isAbsent(){
			return inputReads<0 && inputBases<0 && outputReads<0 && outputBases<0;
		}
		String tsv(){
			return format(inputReads)+"\t"+format(inputBases)+"\t"+
				format(outputReads)+"\t"+format(outputBases);
		}
		private static String format(final long value){return value<0 ? "absent" : Long.toString(value);}
		@Override
		public boolean equals(final Object other){
			if(!(other instanceof Counts)){return false;}
			final Counts counts=(Counts)other;
			return inputReads==counts.inputReads && inputBases==counts.inputBases &&
				outputReads==counts.outputReads && outputBases==counts.outputBases;
		}
		@Override
		public int hashCode(){
			long hash=inputReads;
			hash=hash*31+inputBases;
			hash=hash*31+outputReads;
			hash=hash*31+outputBases;
			return (int)(hash^(hash>>>32));
		}
		private static final long[] ABSENT=new long[] {-1, -1};
		final long inputReads;
		final long inputBases;
		final long outputReads;
		final long outputBases;
	}

	/*--------------------------------------------------------------*/

	private final Path out;
	private final Path repo;
	private final Path launcher;
	private final boolean baseline;
	private int cases;
	private int children;
	private int changed;
	private int controls;
	private int sourceMatches;
	private int correctedMatches;
	private int checks;
	private int failures;

	private static final Sam MAPPED=new Sam("mapped", 0, "chr1", 1, 60, "8=", "*", 0, 0,
		"ACGTACGT", "IIIIIIII");
	private static final Sam REQUIRED=new Sam("required", 16, "chr1", 1, 60, "8=", "*", 0, 0,
		"ACGTACGT", "IIIIIIII");
	private static final Sam UNMAPPED=new Sam("unmapped", 4, "*", 0, 0, "*", "*", 0, 0,
		"ACGTACGT", "IIIIIIII");
	private static final Sam PAIRED=new Sam("paired", 67, "chr1", 1, 60, "8=", "=", 1, 8,
		"ACGTACGT", "IIIIIIII");
	private static final Sam MAXID=new Sam("maxid", 0, "chr1", 1, 60, "8X", "*", 0, 0,
		"ACGTACGT", "IIIIIIII");
	private static final Rows SPLIT=new Rows(3, 8,
		"mapped_1\t4\t*\t0\t0\t*\t*\t0\t0\tACG\tIII",
		"mapped_2\t4\t*\t0\t0\t*\t*\t0\t0\tTAC\tIII",
		"mapped_3\t4\t*\t0\t0\t*\t*\t0\t0\tGT\tII");
	private static final Case[] CASES={
		new Case("legal-split", true, MAPPED, new String[0], MAPPED.rows(), SPLIT),
		control("no-break-combined", PAIRED, new String[]{"minmapq=60", "maxmapq=60",
			"filterbits=512", "requiredbits=2", "mappedonly=t", "primaryonly=t", "pairedonly=t",
			"idfilter=1", "subfilter=0", "minsubs=0", "clipfilter=0", "delfilter=0",
			"insfilter=0", "indelfilter=0", "dellenfilter=9", "inslenfilter=9", "editfilter=0"}),
		control("no-break-unmapped", UNMAPPED, new String[]{"unmappedonly=t"}),
		reject("minmapq", MAPPED, "minmapq=60"),
		fixedOnlyReject("minmapq-fastq", ".fq"),
		fixedOnlyReject("minmapq-fasta", ".fa"),
		reject("maxmapq", MAPPED, "maxmapq=60"),
		reject("filterbits", MAPPED, "filterbits=4"),
		reject("requiredbits", REQUIRED, "requiredbits=16"),
		reject("mappedonly", MAPPED, "mappedonly=t"),
		reject("unmappedonly", UNMAPPED, "unmappedonly=t"),
		reject("primaryonly", MAPPED, "primaryonly=t"),
		reject("pairedonly", PAIRED, "pairedonly=t"),
		reject("unpairedonly", MAPPED, "unpairedonly=t"),
		reject("idfilter", MAPPED, "idfilter=1"),
		reject("maxid", MAXID, "maxid=0.5"),
		reject("subfilter", MAPPED, "subfilter=0"),
		reject("minsubs", MAPPED, "minsubs=0"),
		reject("clipfilter", MAPPED, "clipfilter=0"),
		reject("delfilter", MAPPED, "delfilter=0"),
		reject("insfilter", MAPPED, "insfilter=0"),
		reject("indelfilter", MAPPED, "indelfilter=0"),
		reject("dellenfilter", MAPPED, "dellenfilter=9"),
		reject("inslenfilter", MAPPED, "inslenfilter=9"),
		reject("editfilter", MAPPED, "editfilter=0")
	};
	private static final int FIXED_ONLY_CASES=2;
}
