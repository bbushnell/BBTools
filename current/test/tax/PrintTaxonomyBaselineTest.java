package test.tax;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;

import shared.Timer;
import tax.PrintTaxonomy;
import tax.TaxTree;
import tax.TaxTreeText;

/** Finite baselines for PrintTaxonomy null output and column newline behavior.
 * @author Sucrose */
public final class PrintTaxonomyBaselineTest{

	/** Runs the captured baseline or expected fixed expectations.
	 * @param args Optional {@code expect=baseline} or {@code expect=fixed}
	 * @throws Exception If setup, processing or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		checkNullOutput(baseline);
		checkHealthyNameOutput();
		checkHealthyTextColumn();
		checkReadColumnNewline(baseline);
		System.out.println("PRINT_TAXONOMY_BASELINE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private PrintTaxonomyBaselineTest(){}

	/** Verifies out=null reaches the baseline null TextStreamWriter or the fixed rejection.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If fixture construction or invocation fails */
	private static void checkNullOutput(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-null-");
		try{
			final Path tree=writeTree(root);
			final Throwable failure=constructAndProcess(new String[]{
					"tree="+tree,
					"out=null",
					"missing"});
			if(baseline){
				check(failure instanceof NullPointerException,
						"out=null baseline must dereference the missing primary TextStreamWriter");
				check(failure.getStackTrace().length>0 &&
						failure.getStackTrace()[0].getClassName().equals("tax.PrintTaxonomy") &&
						failure.getStackTrace()[0].getMethodName().equals("printTaxonomy"),
						"out=null baseline must fail at the PrintTaxonomy writer dereference: "+
						failure);
			}else{
				check(failure instanceof IllegalArgumentException,
						"fixed out=null must reject the missing primary TextStreamWriter clearly: "+failure);
				check(failure.getMessage().contains("out"),
						"fixed out=null rejection must guide the output argument: "+failure.getMessage());
			}
		}finally{deleteTempTree(root);}
	}

	/** Verifies command-line name output works with the real String[] constructor.
	 * @throws Exception If setup, output reading or cleanup fails */
	private static void checkHealthyNameOutput() throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-name-");
		try{
			final Path tree=writeTree(root);
			final Path out=root.resolve("out.txt");
			final Throwable failure=process(new PrintTaxonomy(new String[]{
					"tree="+tree,
					"out="+out,
					"missing"}));
			check(failure==null, "healthy name-output control must not fail: "+failure);
			check(text(out).equals("\nmissing:\nCould not find node.\n"),
					"healthy name-output control must use the command-line name path");
		}finally{deleteTempTree(root);}
	}

	/** Verifies tabular TextFile column translation prints exactly one newline per row.
	 * @throws Exception If setup, output reading or cleanup fails */
	private static void checkHealthyTextColumn() throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-text-");
		try{
			final Path tree=writeTree(root);
			final Path in=root.resolve("input.tsv");
			final Path out=root.resolve("out.txt");
			Files.write(in, ("missing\tleft\nabsent\tright\n").getBytes(StandardCharsets.UTF_8));
			final Throwable failure=process(new PrintTaxonomy(new String[]{
					"tree="+tree,
					"in="+in,
					"out="+out,
					"column=0"}));
			check(failure==null, "healthy TextFile column control must not fail: "+failure);
			check(text(out).equals("missing_***NOT_FOUND***\tleft\n"+
					"absent_***NOT_FOUND***\tright\n"),
					"healthy TextFile column control must retain one newline per row");
		}finally{deleteTempTree(root);}
	}

	/** Verifies read column translation writes one extra blank line until fixed.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If setup, output reading or cleanup fails */
	private static void checkReadColumnNewline(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-read-");
		try{
			final Path tree=writeTree(root);
			final Path in=root.resolve("input.fq");
			final Path out=root.resolve("out.txt");
			Files.write(in, ("@missing\tleft\nA\n+\n!\n"+
					"@absent\tright\nC\n+\n!\n").getBytes(StandardCharsets.UTF_8));
			final Throwable failure=process(new PrintTaxonomy(new String[]{
					"tree="+tree,
					"in="+in,
					"out="+out,
					"column=0"}));
			check(failure==null, "read column newline case must not fail: "+failure);
			final String fixed="missing_***NOT_FOUND***\tleft\n"+
					"absent_***NOT_FOUND***\tright\n";
			final String expected=baseline ? fixed.replace("\n", "\n\n") : fixed;
			check(text(out).equals(expected),
					"read column newline "+(baseline ? "baseline" : "fixed")+
					" output mismatch");
		}finally{deleteTempTree(root);}
	}

	/** Writes the finite taxonomy tree used by every case.
	 * @param root Temporary fixture root
	 * @return Serialized tree path
	 * @throws Exception If tree construction or writing fails */
	private static Path writeTree(final Path root) throws Exception{
		final Path tree=root.resolve("tree.tsv");
		TaxTreeText.write(TaxTreeContractTest.tree(TaxTree.GENUS_E), tree.toString(), true);
		return tree;
	}

	/** Invokes the package-private PrintTaxonomy process method.
	 * @param pt Fixture instance to run
	 * @return Thrown failure, or null on success
	 * @throws Exception If reflection itself fails */
	private static Throwable process(final PrintTaxonomy pt) throws Exception{
		final Method process=PrintTaxonomy.class.getDeclaredMethod("process", Timer.class);
		process.setAccessible(true);
		try{
			process.invoke(pt, new Timer());
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	/** Constructs and runs PrintTaxonomy, preserving constructor-time rejections.
	 * @param args Command-line arguments
	 * @return Thrown failure, or null on success
	 * @throws Exception If reflection itself fails */
	private static Throwable constructAndProcess(final String[] args) throws Exception{
		try{
			return process(new PrintTaxonomy(args));
		}catch(IllegalArgumentException e){
			return e;
		}
	}

	/** Reads an output file as UTF-8 text.
	 * @param path Output path
	 * @return File contents
	 * @throws Exception If reading fails */
	private static String text(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
