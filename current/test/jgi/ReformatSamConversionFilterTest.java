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

/** Finite SAM-to-FASTQ identity/edit oracles without splitting or ambiguous M operations.
 * @author Sucrose
 */
public final class ReformatSamConversionFilterTest{

	/** Runs six fresh public-launcher children in the selected phase. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, repoPath=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){repoPath=pair[1];}
			else{throw new IllegalArgumentException("Unknown option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && repoPath!=null, "Set outdir and bbtools.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/reformatsamconversionfiltertest.sh with assertions enabled.");
		final ReformatSamConversionFilterTest test=new ReformatSamConversionFilterTest(
			Paths.get(directory).toAbsolutePath().normalize(), Paths.get(repoPath).toAbsolutePath().normalize(),
			phase.equals("baseline"));
		for(final Case row : ROWS){test.run(row);}
		test.check(test.cases==6, "Expected six cases: "+test.cases);
		test.check(test.children==6, "Expected six children: "+test.children);
		test.check(test.changed==2, "Expected two changed rows: "+test.changed);
		test.check(test.controls==4, "Expected four controls: "+test.controls);
		test.check(test.sourceMatches==(test.baseline ? 6 : 4), "Source matches: "+test.sourceMatches);
		test.check(test.correctedMatches==(test.baseline ? 4 : 6), "Corrected matches: "+test.correctedMatches);
		final String summary="REFORMAT_SAM_CONVERSION_FILTER_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" children="+test.children+
			" changed="+test.changed+" controls="+test.controls+" sourceMatches="+test.sourceMatches+
			" correctedMatches="+test.correctedMatches+" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+test.out);}
		System.out.print(summary);
	}

	/** Binds a new result directory and the real reformat.sh entry point. */
	private ReformatSamConversionFilterTest(final Path out_, final Path repo_, final boolean baseline_) throws Exception{
		out=out_; repo=repo_; baseline=baseline_;
		require(Files.isRegularFile(repo.resolve("reformat.sh")), "Missing reformat.sh");
		Files.createDirectory(out);
		Files.write(out.resolve("cases.tsv"), ("case\tphase\texit\toutputBytes\tsourceMatches\tcorrectedMatches"+
			"\tinputReads\tinputBases\toutputReads\toutputBases\tinputsUnchanged\n").getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
	}

	/** Saves and checks exact FASTQ output, diagnostics, counters and original SAM bytes. */
	private void run(final Case row) throws Exception{
		cases++; children++;
		if(row.sourceReads!=row.correctedReads){changed++;}else{controls++;}
		final Path work=out.resolve(row.name);
		Files.createDirectory(work);
		final Path input=work.resolve("input.sam"), output=work.resolve("output.fq");
		final byte[] original=("@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\n"+
			row.readName+"\t0\tchr1\t1\t60\t"+row.cigar+"\t*\t0\t0\tACGTACGT\tIIIIIIII\n")
			.getBytes(StandardCharsets.US_ASCII);
		Files.write(input, original, StandardOpenOption.CREATE_NEW);
		check(!Files.exists(output), "Output already exists: "+work);
		final byte[] source=expected(row, row.sourceReads), corrected=expected(row, row.correctedReads);
		Files.write(work.resolve("source.fq"), source, StandardOpenOption.CREATE_NEW);
		Files.write(work.resolve("corrected.fq"), corrected, StandardOpenOption.CREATE_NEW);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(repo.resolve("reformat.sh").toString());
		command.addAll(Arrays.asList("-Xmx128m", "-ea", "in="+input, "out="+output, "t=1", "wt=1",
			"tin=1", "tout=1", "int=f", "ow=f", "app=f", "simd=f", "minlength=1", "breaklength=0",
			"qin=33", "qout=33", row.filter));
		if(row.forceParse){command.add("forceparse=t");}
		final StringBuilder argv=new StringBuilder();
		for(final String arg : command){argv.append('\'').append(arg.replace("'", "'\\''")).append("' ");}
		Files.write(work.resolve("argv.txt"), argv.append('\n').toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(work.resolve("stdout.txt").toFile()).redirectError(work.resolve("stderr.txt").toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		Files.write(work.resolve("exit.txt"), (exit+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final boolean exists=Files.isRegularFile(output);
		final byte[] actual=exists ? Files.readAllBytes(output) : null;
		final String stderr=new String(Files.readAllBytes(work.resolve("stderr.txt")), StandardCharsets.UTF_8);
		final long[] counts=counts(stderr);
		final boolean sourceMatches=exists && Arrays.equals(actual, source) && correctCounts(counts, row.sourceReads);
		final boolean correctedMatches=exists && Arrays.equals(actual, corrected) && correctCounts(counts, row.correctedReads);
		if(sourceMatches){this.sourceMatches++;}
		if(correctedMatches){this.correctedMatches++;}
		final boolean unchanged=Arrays.equals(original, Files.readAllBytes(input));
		Files.write(out.resolve("cases.tsv"), (row.name+"\t"+(baseline ? "baseline" : "fixed")+"\t"+exit+
			"\t"+(actual==null ? "absent" : actual.length)+"\t"+sourceMatches+"\t"+correctedMatches+
			"\t"+counts[0]+"\t"+counts[1]+"\t"+counts[2]+"\t"+counts[3]+"\t"+unchanged+"\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
		check(exit==0, "Child failed: "+work);
		check(Files.size(work.resolve("stdout.txt"))==0, "Child wrote stdout: "+work);
		check(!stderr.contains("Exception") && !stderr.contains("java.lang.") && !stderr.contains("Unknown parameter") &&
			!stderr.contains("ERROR") && !stderr.contains("Error"), "Unexpected error diagnostic: "+work);
		check(baseline ? sourceMatches : correctedMatches, "Exact output or counters differ: "+work);
		check(unchanged, "Input changed: "+work);
		check(command.contains("breaklength=0") && command.contains("tin=1") && command.contains("tout=1") &&
			command.contains("wt=1") && command.contains("simd=f") && command.contains("qin=33") &&
			command.contains("qout=33") && command.contains("forceparse=t")==row.forceParse,
			"Lost safe conversion flags: "+work);
	}

	/** Both phases use literal four-line FASTQ or an existing empty output file. */
	private static byte[] expected(final Case row, final int reads){
		return (reads==0 ? "" : "@"+row.readName+"\nACGTACGT\n+\nIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	}

	/** Reads unique complete Input/Output lines, rejecting absent, duplicate or malformed statistics. */
	private static long[] counts(final String stderr){
		final long[] in=countLine(stderr, "Input"), out=countLine(stderr, "Output");
		return new long[]{in[0], in[1], out[0], out[1]};
	}
	private static long[] countLine(final String stderr, final String label){
		final Matcher matcher=Pattern.compile("(?m)^"+label+
			":\\s*([0-9]+)\\s+reads?(?:\\s+\\([0-9]+\\.[0-9]{2}%\\))?"+
			"\\s+([0-9]+)\\s+bases?(?:\\s+\\([0-9]+\\.[0-9]{2}%\\))?\\s*$").matcher(stderr);
		if(!matcher.find()){return new long[]{-1, -1};}
		final long a=Long.parseLong(matcher.group(1)), b=Long.parseLong(matcher.group(2));
		if(matcher.find()){return new long[]{-1, -1};}
		return new long[]{a, b};
	}
	private static boolean correctCounts(final long[] counts, final int outputReads){
		return counts[0]==1 && counts[1]==8 && counts[2]==outputReads && counts[3]==outputReads*8;
	}
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** One independent filter/parse/record combination, with literal expected retained-read counts. */
	private static final class Case{
		Case(final String name_, final String cigar_, final String filter_, final boolean forceParse_,
				final int sourceReads_, final int correctedReads_){
			name=name_; cigar=cigar_; filter=filter_; forceParse=forceParse_;
			sourceReads=sourceReads_; correctedReads=correctedReads_;
			readName=cigar.equals("8X") ? "bad" : "good";
		}
		final String name, cigar, filter, readName;
		final boolean forceParse;
		final int sourceReads, correctedReads;
	}
	private final Path out, repo;
	private final boolean baseline;
	private int cases, children, changed, controls, sourceMatches, correctedMatches, checks, failures;
	private static final Case[] ROWS={
		new Case("identity-bad-default", "8X", "idfilter=0.90", false, 1, 0),
		new Case("identity-bad-forced", "8X", "idfilter=0.90", true, 0, 0),
		new Case("identity-good-default", "8=", "idfilter=0.90", false, 1, 1),
		new Case("edit-bad-default", "8X", "editfilter=0", false, 1, 0),
		new Case("edit-bad-forced", "8X", "editfilter=0", true, 0, 0),
		new Case("edit-good-default", "8=", "editfilter=0", false, 1, 1)
	};
}
