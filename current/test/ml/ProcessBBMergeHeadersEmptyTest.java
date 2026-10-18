package test.ml;

import ml.ProcessBBMergeHeaders;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

import shared.Timer;

/**
 * Preserves ProcessBBMergeHeaders' accepted-row and header-only normal EOF cases.
 * @author Sucrose
 */
public class ProcessBBMergeHeadersEmptyTest{

	/** Runs one ordinary normal-EOF case. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Usage: ProcessBBMergeHeadersEmptyTest <case> <artifact-dir>";
		final String name=args[0];
		final Path root=Paths.get(args[1]);
		assert(Files.isDirectory(root)) :
			"ProcessBBMergeHeadersEmptyTest runner must create the artifact root: "+root;
		final Path dir=root.resolve(name);
		if(!(name.equals("empty") || name.equals("all-filtered") || name.equals("accepted"))){
			throw new IllegalArgumentException("Unknown ProcessBBMergeHeaders empty case: "+name);
		}
		Files.createDirectory(dir);
		if(name.equals("empty")){
			empty(dir);
		}else if(name.equals("all-filtered")){
			allFiltered(dir);
		}else{
			accepted(dir);
		}
		System.err.println("ProcessBBMergeHeadersEmptyTest PASS "+name);
	}

	/** Empty input should still write the synthetic table header in .header output. */
	private static void empty(final Path dir) throws Exception{
		final Path in=dir.resolve("empty.fq");
		Files.write(in, Arrays.asList(), StandardCharsets.US_ASCII);
		runCase(dir, in, "empty.header", HEADER_ONLY, 0);
	}

	/** Nonempty input without accepted BBMerge rows already writes exactly the header. */
	private static void allFiltered(final Path dir) throws Exception{
		final Path in=dir.resolve("all-filtered.fq");
		Files.write(in, Arrays.asList(
			"@plain-read-1",
			"ACGT",
			"+",
			"IIII",
			"@plain-read-2",
			"TGCA",
			"+",
			"JJJJ"
		), StandardCharsets.US_ASCII);
		runCase(dir, in, "all-filtered.header", HEADER_ONLY, 2);
	}

	/** Accepted and rejected rows should keep one synthetic header and one metric row. */
	private static void accepted(final Path dir) throws Exception{
		final Path in=dir.resolve("accepted.fq");
		Files.write(in, Arrays.asList(
			"@"+ACCEPTED_BBMERGE_ID,
			"ACGT",
			"+",
			"IIII",
			"@plain-read",
			"TGCA",
			"+",
			"JJJJ"
		), StandardCharsets.US_ASCII);
		runCase(dir, in, "accepted.header", ACCEPTED_LINES, 2);
	}

	/** Runs the tool and checks the literal .header rows and summary count. */
	private static void runCase(final Path dir, final Path in, final String outName,
			final List<String> expectedLines, final long expectedReads) throws Exception{
		final Path out=dir.resolve(outName);
		final String report=withErr(new ThrowingRunnable(){
			@Override
			public void run() throws Exception{
				new ProcessBBMergeHeaders(new String[]{"in="+in, "out="+out}).process(new Timer());
			}
		});
		Files.write(dir.resolve("report.stderr"), report.getBytes(StandardCharsets.US_ASCII));
		Files.write(dir.resolve("report-normalized.stderr"),
			normalizeReport(report).getBytes(StandardCharsets.US_ASCII));

		final List<String> lines=Files.readAllLines(out, StandardCharsets.US_ASCII);
		assert(lines.equals(expectedLines)) : outName+" .header rows drifted:\n"+
			"Expected:\n"+expectedLines+"\nActual:\n"+lines;
		assert(report.contains("Reads Processed:    "+expectedReads+" \t")) :
			outName+" must count "+expectedReads+" normal EOF input reads exactly:\n"+report;
	}

	/** Removes only elapsed time and rate variability from the report. */
	private static String normalizeReport(final String report){
		final StringBuilder sb=new StringBuilder();
		for(String line : report.split("\n", -1)){
			if(line.startsWith("Time:")){
				line="Time:                         \t<elapsed>";
			}else if(line.startsWith("Reads Processed:")){
				final int tab=line.indexOf('\t');
				assert(tab>=0) :
					"ProcessBBMergeHeaders read summary keeps a tab before the rate: "+line;
				line=line.substring(0, tab+1)+"<rate> reads/sec";
			}
			if(sb.length()>0){sb.append('\n');}
			sb.append(line);
		}
		return sb.toString();
	}

	/** Captures stderr while ProcessBBMergeHeaders writes its report. */
	private static String withErr(final ThrowingRunnable runnable) throws Exception{
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

	/** Runnable that can throw checked exceptions. */
	private static interface ThrowingRunnable{
		void run() throws Exception;
	}

	private static final String BBMERGE_HEADER=
		"#Correct\tminOverlap\tbestOverlap\tbestBadInt\tsecondBestOverlap\tsecondBestBadInt\t"+
		"expectedErrors\tbestExpectedErrors\tbestRatio\tbestBad\tsecondBestRatio\tsecondBestBad\t"+
		"probability";

	private static final String ACCEPTED_BBMERGE_ID=
		"insert=202 synthetic mo=14_r1ee=5.2728_r2ee=3.4856"+
		"_bi=202_bo=98_bb=5.3063_br=0.0598_bbi=6_sbi=270_sbo=30"+
		"_sbb=12.4775_sbr=0.4343_sbbi=14_be=6.5990_pr=0.0007";

	private static final String ACCEPTED_ROW=
		"1\t14\t98\t6\t30\t14\t8.7584\t6.5990\t0.0598\t5.3063\t0.4343\t12.4775\t0.00070000";

	private static final List<String> HEADER_ONLY=Arrays.asList(BBMERGE_HEADER);

	private static final List<String> ACCEPTED_LINES=Arrays.asList(BBMERGE_HEADER, ACCEPTED_ROW);

}
