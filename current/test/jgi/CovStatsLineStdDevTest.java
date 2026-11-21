package test.jgi;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import jgi.CovStatsLine;
import parse.Parse;

/** Finite literal oracles for unavailable merged coverage standard deviations.
 * Uses only the pre-repair public API so the same fixture compiles in both phases.
 * @author Sucrose
 */
public final class CovStatsLineStdDevTest{

	/** Requires phase=baseline|fixed, a new outdir and the actual BBTools root. */
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
		require(assertions, "Run with assertions through dev/covstatslinestddevtest.sh.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		require(root.toString().indexOf('=')<0, "Legacy merger requires an outdir without '='.");
		Files.createDirectory(root);
		final CovStatsLineStdDevTest test=new CovStatsLineStdDevTest(root,
			Paths.get(repoPath).toAbsolutePath().normalize(), phase.equals("fixed"));
		test.numeric();
		test.omittedNumeric();
		test.unavailable();
		test.add(false);
		test.add(true);
		test.rejectedAdd();
		for(final String token : new String[]{"NaN", "Infinity", "-Infinity"}){test.nonfinite(token);}
		test.merge("merge-positive", HEADER, FIRST+"\n"+SECOND+"\n", MERGED);
		test.merge("merge-singleton", HEADER, FIRST+"\n", SINGLETON);
		test.merge("merge-reordered", REORDERED_HEADER, REORDERED_FIRST+"\n"+REORDERED_SECOND+"\n",
			"otuA\tSTD\t10.0000\t2000\t0.4000\t90.0000\t1800\t14\t6\tNA\t1200\t0.2500");
		test.merge("merge-omitted", CORE, CORE_FIRST+"\n"+CORE_SECOND+"\n",
			"otuA\t10.0000\t2000\t0.4000\t90.0000\t1800\t14\t6");
		test.filter(true);
		test.filter(false);
		test.check(test.cases==15, "Expected fifteen cases: "+test.cases);
		test.check(test.direct==9, "Expected nine direct cases: "+test.direct);
		test.check(test.children==6, "Expected six real CLI children: "+test.children);
		test.check(test.inputsUnchanged(), "A literal input changed.");
		final String summary="COV_STATS_STD_DEV_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" direct="+test.direct+" children="+test.children+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+root);}
		System.out.print(summary);
	}

	/** Binds the output and real shell launchers; no production compilation occurs here. */
	private CovStatsLineStdDevTest(final Path root_, final Path repo_, final boolean fixed_){
		root=root_; repo=repo_; fixed=fixed_;
		require(Files.isRegularFile(repo.resolve("mergeOTUs.sh")), "Missing mergeOTUs.sh");
		require(Files.isRegularFile(repo.resolve("filterbycoverage.sh")), "Missing filterbycoverage.sh");
	}

	/** Known singleton values and the existing numeric formatting stay intact. */
	private void numeric() throws Exception{
		final Path work=begin("direct-numeric", true);
		input(work.resolve("input.tsv"), HEADER+"\n"+FIRST+"\n");
		CovStatsLine.initializeHeader(HEADER);
		final CovStatsLine row=new CovStatsLine(FIRST);
		save(work.resolve("observed.tsv"), row.toString()+"\n");
		save(work.resolve("expected.tsv"), DIRECT_SINGLETON+"\n");
		check(DIRECT_SINGLETON.equals(row.toString()), "Numeric singleton output: "+work);
		check(bits(row.stdDev)==bits(1.5), "Numeric singleton storage: "+work);
	}

	/** An omitted deviation column retains the known-zero singleton default in fixed-order output. */
	private void omittedNumeric() throws Exception{
		final Path work=begin("direct-omitted-numeric", true);
		input(work.resolve("input.tsv"), CORE+"\n"+CORE_FIRST+"\n");
		CovStatsLine.initializeHeader(CORE);
		final CovStatsLine row=new CovStatsLine(CORE_FIRST);
		final String expected="ctgA otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3\t0\t0\t0.0000\t0.0000";
		save(work.resolve("observed.tsv"), row.toString()+"\n");
		save(work.resolve("expected.tsv"), expected+"\n");
		check(expected.equals(row.toString()), "Omitted deviation lost known-zero serialization: "+work);
		check(bits(row.stdDev)==bits(0.0), "Omitted deviation lost positive-zero storage: "+work);
	}

	/** Literal NA currently fails synchronously; the repair must round-trip it without numeric NaN. */
	private void unavailable() throws Exception{
		final Path work=begin("direct-na", true);
		final String line=FIRST.substring(0, FIRST.lastIndexOf('\t')+1)+"NA";
		input(work.resolve("input.tsv"), HEADER+"\n"+line+"\n");
		CovStatsLine.initializeHeader(HEADER);
		CovStatsLine row=null;
		NumberFormatException failure=null;
		try{row=new CovStatsLine(line);}catch(NumberFormatException e){failure=e;}
		save(work.resolve("exception.txt"), stack(failure));
		if(row!=null){save(work.resolve("observed.tsv"), row.toString()+"\n");}
		final String expected=DIRECT_SINGLETON.substring(0, DIRECT_SINGLETON.lastIndexOf('\t')+1)+"NA";
		save(work.resolve("expected.tsv"), fixed ? expected+"\n" : "NumberFormatException: NA\n");
		if(fixed){
			check(failure==null, "Literal NA rejected: "+work);
			check(row!=null && expected.equals(row.toString()), "Literal NA output: "+work);
			check(row!=null && bits(row.stdDev)==bits(0.0), "Unavailable storage must be zero: "+work);
		}else{
			check(failure!=null && failure.getClass()==NumberFormatException.class, "Expected exact NA parse failure: "+work);
			check(failure!=null && stack(failure).contains("NA") && stack(failure).contains("jgi.CovStatsLine.parseFinite"),
				"NA failure must come from coverage scalar parsing: "+work);
			check(row==null, "NA baseline unexpectedly constructed a row: "+work);
		}
	}

	/** Successful positive and zero-length merges change only the intended availability contract. */
	private void add(final boolean zero) throws Exception{
		final Path work=begin(zero ? "direct-add-zero" : "direct-add-positive", true);
		final String first=zero ? ZERO_FIRST : FIRST, second=zero ? ZERO_SECOND : SECOND;
		input(work.resolve("input.tsv"), HEADER+"\n"+first+"\n"+second+"\n");
		CovStatsLine.initializeHeader(HEADER);
		final CovStatsLine row=new CovStatsLine(first), other=new CovStatsLine(second);
		final String otherBefore=state(other);
		RuntimeException failure=null;
		try{row.add(other);}catch(RuntimeException e){failure=e;}
		final String expected=(zero ? DIRECT_ZERO : DIRECT_MERGED).replace("STD", fixed ? "NA" : "1.5000");
		save(work.resolve("expected.tsv"), expected+"\n");
		save(work.resolve("observed.tsv"), row.toString()+"\n");
		save(work.resolve("exception.txt"), stack(failure));
		save(work.resolve("source-before.txt"), otherBefore+"\n");
		save(work.resolve("source-after.txt"), state(other)+"\n");
		check(failure==null, "Finite add failed: "+work);
		check(expected.equals(row.toString()), "Merged complete row differs: "+work);
		check(otherBefore.equals(state(other)), "Add mutated source row: "+work);
		check(bits(row.stdDev)==bits(fixed ? 0.0 : 1.5), "Merged numeric storage differs: "+work);
	}

	/** A checked length overflow must leave numeric singleton deviation and every field unchanged. */
	private void rejectedAdd() throws Exception{
		final Path work=begin("direct-rejected-add", true);
		final String first="limit otu\t1.25\t9223372036854775807\t0.4\t0\t1\t2\t3\t8\t0\t0.25\t1.5";
		final String second="extra otu\t2\t1\t0.6\t100\t1\t4\t5\t9\t1\t0.75\t2.5";
		final String expected="limit otu\t1.2500\t9223372036854775807\t0.4000\t0.0000\t1\t2\t3\t8\t0\t0.2500\t1.5000";
		input(work.resolve("input.tsv"), HEADER+"\n"+first+"\n"+second+"\n");
		CovStatsLine.initializeHeader(HEADER);
		final CovStatsLine row=new CovStatsLine(first), other=new CovStatsLine(second);
		final String before=state(row), otherBefore=state(other);
		ArithmeticException failure=null;
		try{row.add(other);}catch(ArithmeticException e){failure=e;}
		save(work.resolve("receiver-before.txt"), before+"\n");
		save(work.resolve("receiver-after.txt"), state(row)+"\n");
		save(work.resolve("source-before.txt"), otherBefore+"\n");
		save(work.resolve("source-after.txt"), state(other)+"\n");
		save(work.resolve("expected.tsv"), expected+"\n");
		save(work.resolve("observed.tsv"), row.toString()+"\n");
		save(work.resolve("exception.txt"), stack(failure));
		check(failure!=null && failure.getClass()==ArithmeticException.class, "Missing checked overflow: "+work);
		check(before.equals(state(row)), "Rejected add mutated receiver: "+work);
		check(expected.equals(row.toString()), "Rejected add changed singleton serialization: "+work);
		check(otherBefore.equals(state(other)), "Rejected add mutated source: "+work);
	}

	/** Accepting literal NA must not admit IEEE nonfinite scalar tokens. */
	private void nonfinite(final String token) throws Exception{
		final Path work=begin("direct-reject-"+token, true);
		final String line=FIRST.substring(0, FIRST.lastIndexOf('\t')+1)+token;
		input(work.resolve("input.tsv"), HEADER+"\n"+line+"\n");
		CovStatsLine.initializeHeader(HEADER);
		NumberFormatException failure=null;
		try{new CovStatsLine(line);}catch(NumberFormatException e){failure=e;}
		save(work.resolve("exception.txt"), stack(failure));
		check(failure!=null && failure.getClass()==NumberFormatException.class, "Nonfinite token accepted: "+work);
		check(failure!=null && failure.getMessage().contains("Nonfinite coverage value") &&
			failure.getMessage().contains("Std_Dev") && failure.getMessage().contains("ctgA otuA") &&
			failure.getMessage().contains(token), "Nonfinite diagnostic lost column/id/token: "+work);
	}

	/** Real merger retains complete headers, column placement, optional omissions and singleton values. */
	private void merge(final String name, final String header, final String rows, final String expectedRow) throws Exception{
		final Path work=begin(name, false);
		final Path in=work.resolve("input.tsv"), out=work.resolve("output.tsv");
		input(in, header+"\n"+rows);
		final String expected=header+"\n"+expectedRow.replace("STD", fixed ? "NA" : "1.5000")+"\n";
		save(work.resolve("expected.tsv"), expected);
		//The legacy main reads only its first two arguments; keep in/out before JVM flags.
		final int exit=run("mergeOTUs.sh", work, "in="+in, "out="+out, "-Xmx128m", "-Xms32m", "-ea", "simd=f");
		final String stderr=read(work.resolve("stderr.txt"));
		check(exit==0, "Merger failed: "+work);
		check(Files.size(work.resolve("stdout.txt"))==0, "Merger wrote stdout: "+work);
		check(clean(stderr), "Merger error diagnostic: "+work);
		check(Files.isRegularFile(out) && expected.equals(read(out)), "Merger complete output differs: "+work);
		check(inputsUnchanged(), "Merger modified input bytes: "+work);
	}

	/** Coverage NA is parsed before input/output pipelines; numeric and repaired NA retain the same FASTA. */
	private void filter(final boolean na) throws Exception{
		final Path work=begin(na ? "filter-na" : "filter-numeric", false);
		final Path in=work.resolve("input.fa"), cov=work.resolve("coverage.tsv");
		final Path cleanOut=work.resolve("clean.fa"), dirtyOut=work.resolve("dirty.fa");
		final StringBuilder bases=new StringBuilder(1000);
		for(int i=0; i<250; i++){bases.append("ACGT");}
		final String fasta=">ctgA\n"+bases+"\n";
		String row=FIRST.replace("ctgA otuA", "ctgA");
		if(na){row=row.substring(0, row.lastIndexOf('\t')+1)+"NA";}
		input(in, fasta);
		input(cov, HEADER+"\n"+row+"\n");
		final boolean success=!na || fixed;
		save(work.resolve("expected-clean.fa"), success ? fasta : "ABSENT\n");
		save(work.resolve("expected-dirty.fa"), success ? "" : "ABSENT\n");
		final int exit=run("filterbycoverage.sh", work, "-Xmx128m", "-Xms32m", "-ea", "in="+in, "cov="+cov,
			"out="+cleanOut, "outdirty="+dirtyOut, "basesundermin=-1", "mincoverage=5", "minpercent=40",
			"minreads=0", "minlength=1", "trim=0", "fastawrap=1000", "t=1", "int=f", "ow=f", "simd=f");
		final String stderr=read(work.resolve("stderr.txt"));
		check(exit==(success ? 0 : 1), "Filter exit differs: "+work);
		check(Files.size(work.resolve("stdout.txt"))==0, "Filter wrote stdout: "+work);
		check(success ? Files.isRegularFile(cleanOut) && fasta.equals(read(cleanOut)) : !Files.exists(cleanOut),
			"Filter clean destination differs: "+work);
		check(success ? Files.isRegularFile(dirtyOut) && Files.size(dirtyOut)==0 : !Files.exists(dirtyOut),
			"Filter dirty destination differs: "+work);
		check(success ? clean(stderr) && stat(stderr, "Reads In", 1) && stat(stderr, "Reads Out", 1) &&
			stat(stderr, "Reads Filtered", 0) : stderr.contains("Exception in thread \"main\" java.lang.NumberFormatException") &&
			stderr.contains("NA") && stderr.contains("jgi.CovStatsLine.parseFinite") &&
			stderr.contains("jgi.FilterByCoverage.process") && !stderr.contains("\nTime:"),
			"Filter diagnostics differ: "+work);
		check(inputsUnchanged(), "Filter modified literal inputs: "+work);
	}

	/** Creates one fresh case directory, keeping direct cases separate from child processes. */
	private Path begin(final String name, final boolean directCase) throws Exception{
		cases++;
		if(directCase){direct++;}else{children++;}
		return Files.createDirectory(root.resolve(name));
	}

	/** Runs a finite actual launcher, saving exact argv and separate diagnostics without a timeout oracle. */
	private int run(final String launcher, final Path work, final String... args) throws Exception{
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(repo.resolve(launcher).toString());
		command.addAll(Arrays.asList(args));
		final StringBuilder argv=new StringBuilder();
		for(final String arg : command){argv.append('\'').append(arg.replace("'", "'\\''")).append("' ");}
		save(work.resolve("argv.txt"), argv.append('\n').toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(work.resolve("stdout.txt").toFile()).redirectError(work.resolve("stderr.txt").toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		save(work.resolve("exit.txt"), exit+"\n");
		return exit;
	}

	/** Saves original input bytes for preservation checks across every later case. */
	private void input(final Path path, final String value) throws Exception{
		final byte[] bytes=value.getBytes(StandardCharsets.UTF_8);
		Files.write(path, bytes, StandardOpenOption.CREATE_NEW);
		original.put(path, bytes);
	}
	private boolean inputsUnchanged() throws Exception{
		for(final Map.Entry<Path, byte[]> entry : original.entrySet()){
			if(!Arrays.equals(entry.getValue(), Files.readAllBytes(entry.getKey()))){return false;}
		}
		return true;
	}
	private static void save(final Path path, final String value) throws Exception{
		Files.write(path, value.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}
	private static String read(final Path path) throws Exception{return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);}
	private static long bits(final double value){return Double.doubleToRawLongBits(value);}

	/** Captures every pre-repair public field plus serialization, including median availability and signed zeros. */
	private static String state(final CovStatsLine row){
		return row.id+"\t"+row.length+"\t"+row.coveredBases+"\t"+row.plusReads+"\t"+row.minusReads+
			"\t"+row.median+"\t"+row.medianAvailable+"\t"+row.underMin+"\t"+bits(row.avgFold)+
			"\t"+bits(row.refGC)+"\t"+bits(row.readGC)+"\t"+bits(row.stdDev)+"\n"+row.toString();
	}
	private static String stack(final Throwable failure){
		if(failure==null){return "none\n";}
		final StringWriter text=new StringWriter();
		failure.printStackTrace(new PrintWriter(text));
		return text.toString();
	}
	private static boolean clean(final String text){
		return !text.contains("Exception in thread") && !text.contains("java.lang.") &&
			!text.contains("ERROR") && !text.contains("Error");
	}
	private static boolean stat(final String text, final String label, final int value){
		return java.util.regex.Pattern.compile("(?m)^"+java.util.regex.Pattern.quote(label)+":\\s*"+value+"(?:\\s|$)")
			.matcher(text).find();
	}
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private final Path root, repo;
	private final boolean fixed;
	private final LinkedHashMap<Path, byte[]> original=new LinkedHashMap<Path, byte[]>();
	private int cases, direct, children, checks, failures;
	private static final String CORE="#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads";
	private static final String HEADER=CORE+"\tMedian_fold\tUnder_min\tRead_GC\tStd_Dev";
	private static final String FIRST="ctgA otuA\t10\t1000\t0.4\t90\t900\t7\t3\t8\t600\t0.25\t1.5";
	private static final String SECOND="ctgB otuA\t10\t1000\t0.4\t90\t900\t7\t3\t8\t600\t0.25\t2.5";
	private static final String ZERO_FIRST="zeroA otuA\t3\t0\t0.4\t0\t0\t0\t0\t8\t0\t0.25\t1.5";
	private static final String ZERO_SECOND="zeroB otuA\t9\t0\t0.6\t0\t0\t0\t0\t8\t0\t0.75\t2.5";
	private static final String DIRECT_SINGLETON="ctgA otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3\t8\t600\t0.2500\t1.5000";
	private static final String SINGLETON="otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3\t8\t600\t0.2500\t1.5000";
	private static final String MERGED="otuA\t10.0000\t2000\t0.4000\t90.0000\t1800\t14\t6\tNA\t1200\t0.2500\tSTD";
	private static final String DIRECT_MERGED="ctgA otuA\t10.0000\t2000\t0.4000\t90.0000\t1800\t14\t6\tNA\t1200\t0.2500\tSTD";
	private static final String DIRECT_ZERO="zeroA otuA\t0.0000\t0\t0.0000\t0.0000\t0\t0\t0\tNA\t0\t0.0000\tSTD";
	private static final String CORE_FIRST="ctgA otuA\t10\t1000\t0.4\t90\t900\t7\t3";
	private static final String CORE_SECOND="ctgB otuA\t10\t1000\t0.4\t90\t900\t7\t3";
	private static final String REORDERED_HEADER="#ID\tStd_Dev\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tMedian_fold\tUnder_min\tRead_GC";
	private static final String REORDERED_FIRST="ctgA otuA\t1.5\t10\t1000\t0.4\t90\t900\t7\t3\t8\t600\t0.25";
	private static final String REORDERED_SECOND="ctgB otuA\t2.5\t10\t1000\t0.4\t90\t900\t7\t3\t8\t600\t0.25";
}
