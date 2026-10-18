package test.ml;

import ml.EnrichFailures;
import ml.RankingEval;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Stream;

/**
 * Validates fixed dim27/dim29 ranking-vector users on supported and malformed
 * vector widths.
 * <p>
 * The accepted width-29 and width-33 cases check query-boundary detection,
 * enrichment duplication and full text reports. Width-28 cases run in child
 * JVMs because the underlying tools reject the vector layout through
 * {@code KillSwitch.assertDie}, preventing an ordinary {@code finally} block
 * from cleaning the failing case; the parent removes the child temp directories.
 * Run this fixture with assertions enabled and serialize invocations that
 * capture shared {@code System.err}.
 *
 * @author Sucrose
 */
public class RankingFixedIndexTest{

	/** Runs all fixture modes, or one named mode for child-JVM rejection checks. */
	public static void main(String[] args) throws Exception{
		if(args.length==1){
			run(args[0]);
			return;
		}
		assert(args.length==0) : "This fixture accepts only one optional case name";
		run("enrich29");
		run("enrich33");
		run("ranking29");
		run("ranking33");
		runRejectJVM("enrich28");
		runRejectJVM("ranking28");
		System.err.println("RankingFixedIndexTest PASS");
	}

	/** Runs one case so malformed-width crashes can be isolated in child JVMs. */
	private static void run(String name) throws Exception{
		final Path dir=Files.createTempDirectory("ranking-fixed-index-");
		try{
			if(name.equals("enrich28")){
				runEnrichRaw(dir, 28);
			}else if(name.equals("ranking28")){
				runRankingRaw(dir, 28);
			}else if(name.equals("enrich29")){
				checkEnrich(dir, 29);
			}else if(name.equals("enrich33")){
				checkEnrich(dir, 33);
			}else if(name.equals("ranking29")){
				checkRanking(dir, 29);
			}else if(name.equals("ranking33")){
				checkRanking(dir, 33);
			}else{
				throw new IllegalArgumentException("Unknown fixture mode: "+name);
			}
			System.err.println("RankingFixedIndexTest PASS "+name);
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path->{
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	/** The width-28 cases must reject before x[28] access and may terminate the VM. */
	private static void runRejectJVM(String name) throws Exception{
		final Path tmp=Files.createTempDirectory("ranking-fixed-index-jvm-");
		try{
			final ArrayList<String> command=new ArrayList<String>();
			command.add(javaBin());
			command.add("-ea");
			command.add("-Xmx256m");
			command.add("-Djava.io.tmpdir="+tmp);
			if(!System.getProperty("java.specification.version", "8").equals("1.8")){
				command.add("--add-modules");
				command.add("jdk.incubator.vector");
			}
			command.add("-cp");
			command.add(System.getProperty("java.class.path"));
			command.add(RankingFixedIndexTest.class.getName());
			command.add(name);

			final Process process=new ProcessBuilder(command).redirectErrorStream(true).start();
			final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
			final InputStream is=process.getInputStream();
			final byte[] temp=new byte[4096];
			for(int len; (len=is.read(temp))>=0; ){buffer.write(temp, 0, len);}
			final int code=process.waitFor();
			final String output=new String(buffer.toByteArray(), StandardCharsets.UTF_8);
			assert(code!=0) : name+" accepted #dims 28";
			assert(output.contains("#dims must declare at least 29 inputs")) :
				name+" rejected with the wrong diagnostic: "+output;
			System.err.println("RankingFixedIndexTest PASS "+name);
		}finally{
			try(Stream<Path> paths=Files.walk(tmp)){
				paths.sorted(Comparator.reverseOrder()).forEach(path->{
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	/** Checks that accepted widths preserve query boundaries and exact enrichment output. */
	private static void checkEnrich(Path dir, int width) throws Exception{
		final Path out=runEnrich(dir, width);
		final String actual=new String(Files.readAllBytes(out), StandardCharsets.US_ASCII);
		final String expected=vectorText(width, "#dims\t"+width+"\t1\n", 3);
		assert(actual.equals(expected)) : "EnrichFailures width "+width+" output drifted:\n"+
			"Expected:\n"+expected+"Actual:\n"+actual;
	}

	/** Checks exact ranking reports for accepted widths and both scoring modes. */
	private static void checkRanking(Path dir, int width) throws Exception{
		final String actual=runRanking(dir, width).replace(File.separatorChar, '/');
		final Path net=dir.resolve("nets").resolve("zero"+width+".bbnet");
		final String expected=rankingReport(width, net).replace(File.separatorChar, '/');
		assert(actual.equals(expected)) : "RankingEval width "+width+" report drifted:\n"+
			"Expected:\n"+expected+"Actual:\n"+actual;
	}

	/** Runs EnrichFailures and returns its output file after checking its exact report. */
	private static Path runEnrich(Path dir, int width) throws Exception{
		final Path out=dir.resolve("enriched"+width+".tsv");
		final Path vectors=writeVectors(dir, width);
		final Path net=writeNet(dir, width);
		final String report=withErr(new ThrowingRunnable(){
			@Override
			public void run() throws Exception{
				EnrichFailures.main(new String[]{"vectors="+vectors, "net="+net, "out="+out, "copies=2"});
			}
		});
		final String expected="Queries: 2  failures: 1 (50.0%)  threshold: 0.1\n"+
			"Output rows: 7 (3 normal + 4 enriched, 2x copies of failures)\n";
		assert(report.equals(expected)) : "EnrichFailures width "+width+" report drifted:\n"+
			"Expected:\n"+expected+"Actual:\n"+report;
		return out;
	}

	/** Runs EnrichFailures without redirecting diagnostics for child-JVM rejection cases. */
	private static void runEnrichRaw(Path dir, int width) throws Exception{
		final Path vectors=writeVectors(dir, width);
		final Path net=writeNet(dir, width);
		EnrichFailures.main(new String[]{"vectors="+vectors, "net="+net,
			"out="+dir.resolve("enriched"+width+".tsv"), "copies=2"});
	}

	/** Runs RankingEval and returns its stderr report. */
	private static String runRanking(Path dir, int width) throws Exception{
		final Path vectors=writeVectors(dir, width);
		final Path netDir=Files.createDirectories(dir.resolve("nets"));
		writeNet(netDir, width);
		return withErr(new ThrowingRunnable(){
			@Override
			public void run() throws Exception{
				RankingEval.main(new String[]{"vectors="+vectors, "nets="+netDir});
			}
		});
	}

	/** Runs RankingEval without redirecting diagnostics for child-JVM rejection cases. */
	private static void runRankingRaw(Path dir, int width) throws Exception{
		final Path vectors=writeVectors(dir, width);
		final Path netDir=Files.createDirectories(dir.resolve("nets"));
		writeNet(netDir, width);
		RankingEval.main(new String[]{"vectors="+vectors, "nets="+netDir});
	}

	/** Writes a valid ranking vector with two queries and the fixed fields populated. */
	private static Path writeVectors(Path dir, int width) throws Exception{
		final Path path=dir.resolve("vectors"+width+".tsv");
		Files.write(path, Arrays.asList(
			"#dims\t"+width+"\t1",
			row(width, 0.8f, 0.2f, 0.2f),
			row(width, 0.1f, 0.2f, 1.0f),
			row(width, 0.5f, 0.3f, 0.8f)
		), StandardCharsets.US_ASCII);
		return path;
	}

	/** Emits one deterministic zero-output linear net whose input width matches the vector file. */
	private static Path writeNet(Path dir, int width) throws Exception{
		final Path path=dir.resolve("zero"+width+".bbnet");
		final StringBuilder cell=new StringBuilder();
		cell.append('C').append(width+1).append(" LINEAR 0");
		for(int i=0; i<width; i++){cell.append(" 0");}
		Files.write(path, Arrays.asList(
			"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
			"#seed 1", "#layers 2", "#dims "+width+" 1", "#edges "+width,
			"#coding decimal", "##layer 1", cell.toString()
		), StandardCharsets.US_ASCII);
		return path;
	}

	/** Builds a row with dim27 and dim29 set when present and the label in column #dims. */
	private static String row(int width, float dim27, float dim29, float label){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<width; i++){
			if(i>0){sb.append('\t');}
			final float value=(i==26 ? dim27 : i==28 ? dim29 : 0);
			sb.append(Float.toString(value));
		}
		return sb.append('\t').append(Float.toString(label)).toString();
	}

	/** Returns expected EnrichFailures output for one duplicated failure query and one control query. */
	private static String vectorText(int width, String header, int failureRepetitions){
		final String row1=row(width, 0.8f, 0.2f, 0.2f);
		final String row2=row(width, 0.1f, 0.2f, 1.0f);
		final String row3=row(width, 0.5f, 0.3f, 0.8f);
		final StringBuilder sb=new StringBuilder(header);
		for(int i=0; i<failureRepetitions; i++){
			sb.append(row1).append('\n');
			sb.append(row2).append('\n');
		}
		sb.append(row3).append('\n');
		return sb.toString();
	}

	/** Expected RankingEval report. */
	private static String rankingReport(int width, Path net){
		return "Loaded 3 hits, 2 queries from "+net.getParent().getParent().resolve("vectors"+width+".tsv")+"\n"+
			"\n=== NET SIZES ===\n"+
			String.format("%-24s %8d bytes  (%,.1f KB)", "zero"+width, net.toFile().length(), net.toFile().length()/1024.0)+"\n"+
			"\n=== TOP-1 SELECTION (which hit does each ranker pick?) ===\n"+
			String.format("%-24s %10s %10s %12s", "ranker", "exact%", ">=genus%", "meanDepth")+"\n"+
			String.format("%-24s %10.2f %10.2f %12.4f", "composite(incumbent)", 0.0, 50.0, 0.5)+"\n"+
			String.format("%-24s %10.2f %10.2f %12.4f", "zero"+width, 0.0, 50.0, 0.5)+"\n"+
			String.format("%-24s %10.2f %10.2f %12.4f", "ORACLE", 50.0, 100.0, 0.9)+"\n"+
			"\n=== POSITION-WEIGHTED SCORE (10/i weights, normalized by oracle) ===\n"+
			String.format("%-24s %10s %12s %12s", "ranker", "nWDCG", "rawWDCG", "oracleWDCG")+"\n"+
			String.format("%-24s %10.4f %12.4f %12.4f", "composite(incumbent)", 15.0/19, 7.5, 9.5)+"\n"+
			String.format("%-24s %10.4f %12.4f %12.4f", "zero"+width, 15.0/19, 7.5, 9.5)+"\n"+
			"\n=== WITHIN-QUERY CONCORDANCE vs LCA depth (0.5=random) ===\n"+
			String.format("%-24s %10s %10s %8s %10s", "ranker", "withTies", "noTies", "tie%", "pairs")+"\n"+
			String.format("%-24s %10.4f %10.4f %7.1f%% %10d", "composite", 0.0, 0.0, 0.0, 1)+"\n"+
			String.format("%-24s %10.4f %10.4f %7.1f%% %10d", "zero"+width, 0.5, 0.0, 100.0, 1)+"\n"+
			"\n=== COMPOSITE GRADE (10*meanTop1Depth + concordance, higher=better) ===\n"+
			String.format("%-24s %8s %10s", "ranker", "KB", "grade")+"\n"+
			String.format("%-24s %8s %10.4f", "composite", "-", 5.0)+"\n"+
			String.format("%-24s %7.1f %10.4f", "zero"+width, net.toFile().length()/1024.0, 5.5)+"\n";
	}

	/** Captures stderr while a tool writes its report. */
	private static String withErr(ThrowingRunnable runnable) throws Exception{
		final PrintStream old=System.err;
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(buffer, true, "UTF-8");
		try{
			System.setErr(ps);
			runnable.run();
		}finally{
			System.setErr(old);
			ps.close();
		}
		return new String(buffer.toByteArray(), StandardCharsets.UTF_8);
	}

	/** Locates the current JVM executable for child malformed-vector cases. */
	private static String javaBin(){
		return System.getProperty("java.home")+File.separator+"bin"+File.separator+"java";
	}

	private static interface ThrowingRunnable{
		void run() throws Exception;
	}
}
