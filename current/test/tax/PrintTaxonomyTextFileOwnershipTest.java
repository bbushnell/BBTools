package test.tax;

import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import tax.PrintTaxonomy;
import tax.TaxTree;
import tax.TaxTreeText;

/** Finite baselines for PrintTaxonomy TextFile ownership.
 * @author Sucrose */
public final class PrintTaxonomyTextFileOwnershipTest{

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
		checkNormalCompletion(baseline);
		checkBodyFailureClosesInput(baseline);
		checkReportedCloseError(baseline);
		System.out.println("PRINT_TAXONOMY_TEXTFILE_OWNERSHIP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private PrintTaxonomyTextFileOwnershipTest(){}

	/** Verifies that processFile owns and closes its TextFile on a normal exit.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If setup, invocation or cleanup fails */
	private static void checkNormalCompletion(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-textfile-normal-");
		ScriptedTextFile tf=null;
		try{
			final Path tree=writeTree(root);
			final Path out=root.resolve("out.txt");
			final PrintTaxonomy pt=new PrintTaxonomy(new String[]{
					"tree="+tree,
					"out="+out,
					"column=0"});
			tf=new ScriptedTextFile(root.resolve("input.tsv"), "missing\tleft");
			final TextStreamWriter tsw=new TextStreamWriter(out.toString(), true, false, false);
			tsw.start();
			final boolean[] writerError=new boolean[1];
			final Throwable failure=processFileAndCloseWriter(pt, tf, tsw, writerError);
			check(failure==null, "normal TextFile processing must not fail: "+failure);
			check(!writerError[0], "normal TextFile output must close cleanly");
			check(text(out).equals("missing_***NOT_FOUND***\tleft\n"),
					"normal TextFile processing must preserve the translated row");
			check(tf.isOpen()==baseline,
					"TextFile normal completion open state must match expectation");
			check(!pt.errorState, "normal TextFile close must not set errorState");
		}finally{
			if(tf!=null && tf.isOpen()){tf.close();}
			deleteTempTree(root);
		}
	}

	/** Verifies that a body failure remains primary while processFile still closes TextFile.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If setup, invocation or cleanup fails */
	private static void checkBodyFailureClosesInput(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-textfile-failure-");
		ScriptedTextFile tf=null;
		try{
			final Path tree=writeTree(root);
			final Path out=root.resolve("out.txt");
			final PrintTaxonomy pt=new PrintTaxonomy(new String[]{
					"tree="+tree,
					"out="+out,
					"column=0"});
			tf=new ScriptedTextFile(root.resolve("input.tsv"), "missing\tleft");
			final Throwable failure=processFile(pt, tf, null);
			check(failure instanceof NullPointerException,
					"TextFile body failure must preserve the original writer NPE: "+failure);
			check(topFrame(failure, "tax.PrintTaxonomy", "processFile"),
					"TextFile body failure must originate in processFile: "+failure);
			check(tf.isOpen()==baseline,
					"TextFile exceptional completion open state must match expectation");
			check(!pt.errorState, "a normal close after body failure must not set errorState");
		}finally{
			if(tf!=null && tf.isOpen()){tf.close();}
			deleteTempTree(root);
		}
	}

	/** Verifies that a reported TextFile close error reaches PrintTaxonomy.errorState.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If setup, invocation or cleanup fails */
	private static void checkReportedCloseError(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-textfile-close-");
		ScriptedTextFile tf=null;
		InputStream original=null;
		try{
			final Path tree=writeTree(root);
			final Path out=root.resolve("out.txt");
			final PrintTaxonomy pt=new PrintTaxonomy(new String[]{
					"tree="+tree,
					"out="+out,
					"column=0"});
			tf=new ScriptedTextFile(root.resolve("input.tsv"), "missing\tleft");
			original=tf.is;
			final CloseErrorInputStream closeError=new CloseErrorInputStream();
			tf.is=closeError;
			final TextStreamWriter tsw=new TextStreamWriter(out.toString(), true, false, false);
			tsw.start();
			final boolean[] writerError=new boolean[1];
			final Throwable failure=processFileAndCloseWriter(pt, tf, tsw, writerError);
			check(failure==null, "reported-close TextFile processing body must not fail: "+failure);
			check(!writerError[0], "reported-close TextFile output must close cleanly");
			check(text(out).equals("missing_***NOT_FOUND***\tleft\n"),
					"reported-close TextFile processing must preserve the translated row");
			check(closeError.closed!=baseline,
					"artificial raw stream close state must match expectation");
			check(tf.isOpen()==baseline,
					"reported-close TextFile open state must match expectation");
			check(pt.errorState!=baseline,
					"reported close status must reach PrintTaxonomy.errorState only after the fix");
			if(baseline){
				tf.is=original;
				original=null;
			}
		}finally{
			if(tf!=null && tf.isOpen()){
				if(original!=null){tf.is=original;}
				tf.close();
			}
			deleteTempTree(root);
		}
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

	/** Invokes the package-private PrintTaxonomy.processFile method.
	 * @param pt Fixture instance to run
	 * @param tf Text input owned by processFile after the fix
	 * @param tsw Output writer, or null to force a body failure
	 * @return Thrown failure, or null on success
	 * @throws Exception If reflection itself fails */
	private static Throwable processFile(final PrintTaxonomy pt, final TextFile tf,
			final TextStreamWriter tsw) throws Exception{
		final Method processFile=PrintTaxonomy.class.getDeclaredMethod("processFile",
				TextFile.class, TextStreamWriter.class);
		processFile.setAccessible(true);
		try{
			processFile.invoke(pt, tf, tsw);
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	/** Invokes processFile and closes the writer while preserving reflection failures.
	 * @param pt Fixture instance to run
	 * @param tf Text input owned by processFile after the fix
	 * @param tsw Output writer to close after invocation
	 * @param writerError Single-value holder for the reported close status
	 * @return Thrown failure, or null on success
	 * @throws Exception If reflection itself fails */
	private static Throwable processFileAndCloseWriter(final PrintTaxonomy pt,
			final TextFile tf, final TextStreamWriter tsw, final boolean[] writerError)
			throws Exception{
		Throwable reflectionFailure=null;
		try{
			return processFile(pt, tf, tsw);
		}catch(Exception | Error e){
			reflectionFailure=e;
			throw e;
		}finally{
			try{
				writerError[0]=tsw.poisonAndWait();
			}catch(RuntimeException | Error e){
				if(reflectionFailure==null){throw e;}
				reflectionFailure.addSuppressed(e);
			}
		}
	}

	/** Checks a reflected failure's top frame.
	 * @param failure Throwable to inspect
	 * @param className Required class name
	 * @param methodName Required method name
	 * @return Whether the top frame matches */
	private static boolean topFrame(final Throwable failure, final String className,
			final String methodName){
		return failure.getStackTrace().length>0 &&
				failure.getStackTrace()[0].getClassName().equals(className) &&
				failure.getStackTrace()[0].getMethodName().equals(methodName);
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

	/** TextFile with finite scripted logical lines and a real final close path. */
	private static final class ScriptedTextFile extends TextFile{

		/** Opens a real backing file and installs a finite logical script.
		 * @param backing Real file path used for final TextFile.close()
		 * @param lines Logical lines to emit before EOF
		 * @throws IOException If the backing file cannot be created */
		ScriptedTextFile(final Path backing, final String... lines) throws IOException{
			super(createBacking(backing));
			this.lines=lines;
		}

		/** Returns the next scripted line or EOF without touching the real reader.
		 * @return Next scripted line, or null at EOF */
		@Override
		public String nextLine(){return next<lines.length ? lines[next++] : null;}

		/** Finite logical-line script. */
		private final String[] lines;
		/** Next script index. */
		private int next=0;
	}

	/** Creates a real empty file before the TextFile constructor opens it.
	 * @param backing Owned backing path
	 * @return String form of backing
	 * @throws IOException If the file cannot be created */
	private static String createBacking(final Path backing) throws IOException{
		Files.write(backing, new byte[0]);
		return backing.toString();
	}

	/** InputStream that reports a bounded close error after normal wrappers close. */
	private static final class CloseErrorInputStream extends InputStream{

		/** Constructs an unclosed artificial raw stream. */
		private CloseErrorInputStream(){}

		/** Reports EOF if accidentally read.
		 * @return Always -1 to report EOF */
		@Override
		public int read(){return -1;}

		/** Reports one artificial raw-stream close error. */
		@Override
		public void close() throws IOException{
			closed=true;
			throw new IOException("artificial TextFile raw close failure");
		}

		/** Whether the artificial raw stream was closed. */
		private boolean closed=false;
	}
}
