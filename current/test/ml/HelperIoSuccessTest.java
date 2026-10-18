package test.ml;

import ml.ErrorTaper;
import ml.ReduceColumns;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * Preserves ordinary helper output while I/O completion handling is tightened.
 * Each named case runs the real helper main with small literal inputs and
 * checks exact retained output under the artifact directory supplied by the
 * caller.
 *
 * @author Sucrose
 */
public class HelperIoSuccessTest{

	/** Runs one success case. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Usage: HelperIoSuccessTest <case> <artifact-dir>";
		final String name=args[0];
		final Path dir=Paths.get(args[1]);
		assert(Files.isDirectory(dir)) : "Artifact directory must already exist: "+dir;
		if(name.equals("taper")){
			taper(dir);
		}else if(name.equals("reduce-explicit")){
			reduceExplicit(dir);
		}else if(name.equals("reduce-open")){
			reduceOpen(dir);
		}else{
			throw new IllegalArgumentException("Unknown helper I/O case: "+name);
		}
		System.err.println("HelperIoSuccessTest PASS "+name);
	}

	/** Runs ErrorTaper with three retained rows sorted by distinct absolute error. */
	private static void taper(final Path dir) throws Exception{
		final Path vectors=dir.resolve("vectors.tsv");
		final Path net=dir.resolve("zero.bbnet");
		final Path out=dir.resolve("tapered.tsv");
		writeVectors(vectors);
		writeNet(net);

		final String report=withErr(new ThrowingRunnable(){
			@Override
			public void run() throws Exception{
				ErrorTaper.main(new String[]{
					"vectors="+vectors, "net="+net, "out="+out, "floor=1", "seed=42"});
			}
		});
		final String expectedOut=
			"#dims\t33\t1\n"+
			vectorRow(0.75f)+"\n"+
			vectorRow(0.5f)+"\n"+
			vectorRow(0.1f)+"\n";
		final String actualOut=read(out);
		assert(actualOut.equals(expectedOut)) : "ErrorTaper output drifted:\n"+
			"Expected:\n"+expectedOut+"Actual:\n"+actualOut;
		final String expectedReport=
			"Loaded 3 rows, NI=33\n"+
			"Kept 3 of 3 (100.0%). floor=1.00. Error range: 0.1000 - 0.7500 (median 0.5000)\n";
		assert(report.equals(expectedReport)) : "ErrorTaper report drifted:\n"+
			"Expected:\n"+expectedReport+"Actual:\n"+report;
	}

	/** Runs ReduceColumns with an explicit range plus singleton. */
	private static void reduceExplicit(final Path dir) throws Exception{
		final Path in=writeReduceInput(dir);
		final Path out=dir.resolve("explicit.tsv");
		ReduceColumns.main(new String[]{in.toString(), out.toString(), "0-1", "3"});
		final String expected="#dims\t2\t1\n"+
			"a\tb\td\n"+
			"1\t2\t4\n";
		final String actual=read(out);
		assert(actual.equals(expected)) : "ReduceColumns explicit output drifted:\n"+
			"Expected:\n"+expected+"Actual:\n"+actual;
	}

	/** Runs ReduceColumns with the open-ended suffix range. */
	private static void reduceOpen(final Path dir) throws Exception{
		final Path in=writeReduceInput(dir);
		final Path out=dir.resolve("open.tsv");
		ReduceColumns.main(new String[]{in.toString(), out.toString(), "2+"});
		final String expected="#dims\t2\t1\n"+
			"c\td\te\n"+
			"3\t4\t5\n";
		final String actual=read(out);
		assert(actual.equals(expected)) : "ReduceColumns open output drifted:\n"+
			"Expected:\n"+expected+"Actual:\n"+actual;
	}

	/** Writes three vectors whose zero-net absolute errors are distinct. */
	private static void writeVectors(final Path path) throws Exception{
		Files.write(path, Arrays.asList(
			"#dims\t33\t1",
			vectorRow(0.1f),
			vectorRow(0.75f),
			vectorRow(0.5f)
		), StandardCharsets.US_ASCII);
	}

	/** Writes a two-data-row column matrix with a skipped comment line. */
	private static Path writeReduceInput(final Path dir) throws Exception{
		final Path path=dir.resolve("matrix.tsv");
		Files.write(path, Arrays.asList(
			"#comment\tignored",
			"a\tb\tc\td\te",
			"1\t2\t3\t4\t5"
		), StandardCharsets.US_ASCII);
		return path;
	}

	/** Emits one deterministic zero-output linear net with 33 inputs. */
	private static void writeNet(final Path path) throws Exception{
		final StringBuilder cell=new StringBuilder();
		cell.append("C34 LINEAR 0");
		for(int i=0; i<33; i++){cell.append(" 0");}
		Files.write(path, Arrays.asList(
			"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
			"#seed 1", "#layers 2", "#dims 33 1", "#edges 33", "#coding decimal",
			"##layer 1", cell.toString()
		), StandardCharsets.US_ASCII);
	}

	/** Builds one 33-input zero feature row with the supplied label. */
	private static String vectorRow(final float label){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<33; i++){
			if(i>0){sb.append('\t');}
			sb.append("0.0");
		}
		return sb.append('\t').append(Float.toString(label)).toString();
	}

	/** Reads a US-ASCII text file. */
	private static String read(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Captures stderr while a tool writes its report. */
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
}
