package test.ml;

import ml.ErrorTaper;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * Exercises ErrorTaper accepted-row edges independently of score ordering.
 *
 * @author Sucrose
 */
public class ErrorTaperAcceptedRowsTest{

	/** Runs one accepted-row case. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Usage: ErrorTaperAcceptedRowsTest <case> <artifact-dir>";
		final String name=args[0];
		final Path dir=Paths.get(args[1]);
		assert(Files.isDirectory(dir)) : "Artifact directory must already exist: "+dir;
		if(name.equals("empty-no-header")){
			emptyNoHeader(dir);
		}else if(name.equals("empty")){
			empty(dir);
		}else if(name.equals("malformed-only")){
			malformedOnly(dir);
		}else if(name.equals("singleton")){
			singleton(dir);
		}else if(name.equals("ordinary2+control")){
			ordinaryControl(dir);
		}else{
			throw new IllegalArgumentException("Unknown ErrorTaper accepted-row case: "+name);
		}
		System.err.println("ErrorTaperAcceptedRowsTest PASS "+name);
	}

	/** Zero accepted rows with no header emit nothing and retain the default width. */
	private static void emptyNoHeader(final Path dir) throws Exception{
		final Path vectors=dir.resolve("empty-no-header.tsv");
		final Path net=dir.resolve("zero.bbnet");
		final Path out=dir.resolve("empty-no-header.out.tsv");
		Files.write(vectors, new byte[0]);
		writeNet(net);

		final String report=taper(vectors, net, out, "0.10");
		final String actualOut=read(out);
		assert(actualOut.length()==0) : "Expected no output from no-header zero rows, found:\n"+
			actualOut;
		final String expectedReport="Loaded 0 rows, NI=33\n"+
			"Kept 0 of 0 (0.0%). floor=0.10. Error range: NA - NA (median NA)\n";
		assert(report.equals(expectedReport)) : "No-header empty report drifted:\n"+
			"Expected:\n"+expectedReport+"Actual:\n"+report;
	}

	/** Zero accepted rows preserve #dims when it is the retained header. */
	private static void empty(final Path dir) throws Exception{
		final Path vectors=dir.resolve("empty.tsv");
		final Path net=dir.resolve("zero.bbnet");
		final Path out=dir.resolve("empty.out.tsv");
		Files.write(vectors, Arrays.asList("#dims\t1\t1"), StandardCharsets.US_ASCII);
		writeNet(net);

		final String report=taper(vectors, net, out, "0.10");
		assertHeaderOnlyOutput(out, "#dims\t1\t1\n");
		final String expectedReport="Loaded 0 rows, NI=1\n"+
			"Kept 0 of 0 (0.0%). floor=0.10. Error range: NA - NA (median NA)\n";
		assert(report.equals(expectedReport)) : "Empty accepted-row report drifted:\n"+
			"Expected:\n"+expectedReport+"Actual:\n"+report;
	}

	/** Comment and short rows still count as zero accepted rows. */
	private static void malformedOnly(final Path dir) throws Exception{
		final Path vectors=dir.resolve("malformed.tsv");
		final Path net=dir.resolve("zero.bbnet");
		final Path out=dir.resolve("malformed.out.tsv");
		Files.write(vectors, Arrays.asList("#dims\t1\t1", "0.0", "#ignored", ""),
			StandardCharsets.US_ASCII);
		writeNet(net);

		final String report=taper(vectors, net, out, "0.10");
		assertHeaderOnlyOutput(out, "#ignored\n");
		final String expectedReport="Loaded 0 rows, NI=1\n"+
			"Kept 0 of 0 (0.0%). floor=0.10. Error range: NA - NA (median NA)\n";
		assert(report.equals(expectedReport)) : "Malformed-only report drifted:\n"+
			"Expected:\n"+expectedReport+"Actual:\n"+report;
	}

	/** One accepted row is kept instead of being rejected by a NaN probability. */
	private static void singleton(final Path dir) throws Exception{
		final Path vectors=dir.resolve("singleton.tsv");
		final Path net=dir.resolve("zero.bbnet");
		final Path out=dir.resolve("singleton.out.tsv");
		Files.write(vectors, Arrays.asList("#dims\t1\t1", "0.0\t0.5"),
			StandardCharsets.US_ASCII);
		writeNet(net);

		final String report=taper(vectors, net, out, "0.10");
		final String expectedOut="#dims\t1\t1\n0.0\t0.5\n";
		final String actualOut=read(out);
		assert(actualOut.equals(expectedOut)) : "Singleton output drifted:\n"+
			"Expected:\n"+expectedOut+"Actual:\n"+actualOut;
		final String expectedReport="Loaded 1 rows, NI=1\n"+
			"Kept 1 of 1 (100.0%). floor=0.10. Error range: 0.5000 - 0.5000 (median 0.5000)\n";
		assert(report.equals(expectedReport)) : "Singleton report drifted:\n"+
			"Expected:\n"+expectedReport+"Actual:\n"+report;
	}

	/** The ordinary two-or-more path still sorts retained rows by descending error. */
	private static void ordinaryControl(final Path dir) throws Exception{
		final Path vectors=dir.resolve("ordinary.tsv");
		final Path net=dir.resolve("zero.bbnet");
		final Path out=dir.resolve("ordinary.out.tsv");
		Files.write(vectors, Arrays.asList(
			"#dims\t1\t1",
			"0.0\t0.2",
			"0.0\t0.7",
			"0.0\t0.4"
		), StandardCharsets.US_ASCII);
		writeNet(net);

		final String report=taper(vectors, net, out, "1");
		final String expectedOut="#dims\t1\t1\n"+
			"0.0\t0.7\n"+
			"0.0\t0.4\n"+
			"0.0\t0.2\n";
		final String actualOut=read(out);
		assert(actualOut.equals(expectedOut)) : "Ordinary output drifted:\n"+
			"Expected:\n"+expectedOut+"Actual:\n"+actualOut;
		final String expectedReport=
			"Loaded 3 rows, NI=1\n"+
			"Kept 3 of 3 (100.0%). floor=1.00. Error range: 0.2000 - 0.7000 (median 0.4000)\n";
		assert(report.equals(expectedReport)) : "Ordinary report drifted:\n"+
			"Expected:\n"+expectedReport+"Actual:\n"+report;
	}

	/** Runs ErrorTaper and captures its stderr report. */
	private static String taper(final Path vectors, final Path net, final Path out,
		final String floor) throws Exception{
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream old=System.err;
		final PrintStream ps=new PrintStream(buffer, true, "UTF-8");
		try{
			System.setErr(ps);
			ErrorTaper.main(new String[]{
				"vectors="+vectors, "net="+net, "out="+out, "floor="+floor, "seed=42"});
		}finally{
			System.setErr(old);
			ps.close();
			Files.write(vectors.getParent().resolve("report.stderr"), buffer.toByteArray());
		}
		return new String(buffer.toByteArray(), StandardCharsets.UTF_8);
	}

	/** Writes a one-input zero-output linear net. */
	private static void writeNet(final Path path) throws Exception{
		Files.write(path, Arrays.asList(
			"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
			"#seed 1", "#layers 2", "#dims 1 1", "#edges 1", "#coding decimal",
			"##layer 1", "C2 LINEAR 0 0"
		), StandardCharsets.US_ASCII);
	}

	/** Verifies that zero accepted rows preserve only the retained header. */
	private static void assertHeaderOnlyOutput(final Path out, final String expectedOut)
		throws Exception{
		final String actualOut=read(out);
		assert(actualOut.equals(expectedOut)) :
			"Expected header-only output from zero accepted rows:\n"+
			"Expected:\n"+expectedOut+"Actual:\n"+
			actualOut;
	}

	/** Reads a US-ASCII text file. */
	private static String read(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
}
