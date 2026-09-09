package ml;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Stream;

/**
 * Fixture for {@link BBNetCompare}: one PASS case (comment lines ignored, exact PASS line
 * incl. the worst cell) and every fatal condition, each asserted by its own message:
 * row count mismatch, row ID mismatch, width mismatch (either file), non-finite value
 * (either file), over tolerance, invalid tolerance, zero expected rows, no value columns,
 * wrong argument count.
 *
 * <p>Run with assertions on: {@code java -ea ml.BBNetCompareTest}</p>
 *
 * @author UMP45 (2026-09-08)
 */
public class BBNetCompareTest {

	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("bbnet-compare-test-");
		try{
			final Path e=w(dir, "e.tsv", "#row\tpred0\tpred1", "0\t0.10000000\t0.50000000", "1\t0.20000000\t0.60000000");
			final Path a=w(dir, "a.tsv", "# different comment", "0\t0.10000000\t0.50000000", "1\t0.20000000\t0.60040000");

			final String pass=capture(new String[]{e.toString(), a.toString(), "0.001"});
			// The worst cell is |0.60040000-0.60000000| computed in double exactly as the comparator does;
			// its decimal text is whatever Double.toString gives for that double, not a hand-typed literal.
			final double worst=Math.abs(Double.parseDouble("0.60040000")-Double.parseDouble("0.60000000"));
			final String expectedPass="PASS rows=2 outputs=2 max_abs_diff="+worst+" worst_row=1 worst_col=2 tol=0.001\n";
			check(pass.equals(expectedPass), "unexpected PASS line: "+pass+" (expected "+expectedPass+")");
			System.err.println("BBNetCompareTest (1) PASS: exact PASS line with worst cell");

			reject("over tolerance", e, a, "0.0001", "FAIL: max|diff|=");
			reject("row count mismatch", e, w(dir, "a1.tsv", "0\t0.1\t0.5"), "0.001", "row count mismatch");
			reject("row ID mismatch", e, w(dir, "a2.tsv", "0\t0.1\t0.5", "9\t0.2\t0.6"), "0.001", "row ID mismatch");
			reject("actual width mismatch", e, w(dir, "a3.tsv", "0\t0.1\t0.5", "1\t0.2"), "0.001", "width 2 != 3");
			reject("expected width mismatch", w(dir, "e1.tsv", "0\t0.1\t0.5", "1\t0.2"), a, "0.001", "width 2 != 3");
			reject("non-finite actual", e, w(dir, "a4.tsv", "0\t0.1\tNaN", "1\t0.2\t0.6"), "0.001", "non-finite value NaN");
			reject("non-finite expected", w(dir, "e2.tsv", "0\tInfinity\t0.5", "1\t0.2\t0.6"), a, "0.001", "non-finite value Infinity");
			reject("invalid tolerance 0", e, a, "0", "invalid tolerance");
			reject("invalid tolerance NaN", e, a, "NaN", "invalid tolerance");
			reject("zero expected rows", w(dir, "e3.tsv", "# only a comment"), a, "0.001", "zero data rows");
			reject("no value columns", w(dir, "e4.tsv", "0", "1"), w(dir, "a5.tsv", "0", "1"), "0.001", "no value columns");
			rejectArgs("wrong argument count", new String[]{e.toString(), a.toString()}, "Usage");
			System.err.println("BBNetCompareTest (2) PASS: 12 fatal conditions each rejected for the stated reason");

			System.err.println("BBNetCompareTest: PASS");
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}catch(Exception ex){throw new RuntimeException(ex);}
				});
			}
		}
	}

	private static Path w(Path dir, String name, String... lines) throws Exception{
		final Path p=dir.resolve(name);
		Files.write(p, Arrays.asList(lines), StandardCharsets.US_ASCII);
		return p;
	}

	private static String capture(String[] args) throws Exception{
		final PrintStream old=System.out;
		final ByteArrayOutputStream bo=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(bo, true, "UTF-8");
		System.setOut(ps);
		try{BBNetCompare.main(args);}
		finally{System.setOut(old); ps.flush();}
		return new String(bo.toByteArray(), StandardCharsets.UTF_8);
	}

	private static void reject(String what, Path e, Path a, String tol, String fragment) throws Exception{
		rejectArgs(what, new String[]{e.toString(), a.toString(), tol}, fragment);
	}

	private static void rejectArgs(String what, String[] args, String fragment) throws Exception{
		final PrintStream old=System.out;
		System.setOut(new PrintStream(new ByteArrayOutputStream(), true, "UTF-8"));
		try{
			BBNetCompare.main(args);
		}catch(RuntimeException ex){
			check(ex.getMessage()!=null && ex.getMessage().contains(fragment),
				what+": wrong rejection message: "+ex.getMessage()+" (expected '"+fragment+"')");
			return;
		}finally{
			System.setOut(old);
		}
		throw new AssertionError(what+": was accepted");
	}

	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
}
