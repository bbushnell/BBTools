package test.tax;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Shared;
import tax.PrintTaxonomy;
import tax.TaxNode;
import tax.TaxTree;
import tax.TaxTreeText;

/** Finite baselines for PrintTaxonomy CRIS and outer-writer cleanup after processing failures.
 * @author Sucrose */
public final class PrintTaxonomyCrisOwnershipTest{

	/** Marker printed and written only after the intended processReads NPE is observed. */
	private static final String CHILD_PROCESS_READS_NPE="CHILD_PROCESS_READS_NPE";
	/** Marker printed and written only after the intended parse failure is observed. */
	private static final String CHILD_OUTER_WRITER_SENTINEL="CHILD_OUTER_WRITER_SENTINEL";
	/** Child mode for the exceptional CRIS cleanup baseline. */
	private static final String CHILD_CRIS_MODE="child-cris";
	/** Child mode for the reported CRIS close-error case. */
	private static final String CHILD_CRIS_CLOSE_MODE="child-cris-close";
	/** Child mode for the exceptional outer-writer cleanup baseline. */
	private static final String CHILD_WRITER_MODE="child-writer";
	/** Marker printed and written only after the close-error CRIS child passes. */
	private static final String CHILD_CRIS_CLOSE="CHILD_CRIS_CLOSE";
	/** Nonexistent path that ReadWrite maps to System.in for the close-error case. */
	private static final String STDIN_CRIS_CLOSE_NAME="stdin.printtaxonomy-close.fastq.gz";
	/** Exact sentinel thrown after the real outer TextStreamWriter starts. */
	private static final RuntimeException OUTER_WRITER_SENTINEL=
			new RuntimeException("artificial parseNodeFromHeader failure");
	/** Milliseconds to wait for child startup or fixed child completion on a loaded login node. */
	private static final long CHILD_WAIT_MILLIS=30000;
	/** Explicit CRIS output buffer count for the exceptional child case. */
	private static final int HANG_BUFFERS=2;
	/** Explicit read buffer length for the exceptional child case. */
	private static final int HANG_BUFFER_LEN=4;
	/** Number of reads large enough to exhaust the explicit CRIS and raw producer queues. */
	private static final int HANG_READS=2000;

	/** Runs the captured baseline or expected fixed expectations.
	 * @param args Optional {@code expect=baseline} or {@code expect=fixed}, or child-mode args
	 * @throws Exception If setup, processing or cleanup fails */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && args[0].equals(CHILD_CRIS_MODE)){
			childCrisMain(args);
			return;
		}
		if(args.length>0 && args[0].equals(CHILD_CRIS_CLOSE_MODE)){
			childCrisCloseMain(args);
			return;
		}
		if(args.length>0 && args[0].equals(CHILD_WRITER_MODE)){
			childWriterMain(args);
			return;
		}

		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		checkHealthyReadColumn();
		checkExceptionalReadCleanup(baseline);
		checkReportedCrisCloseError(baseline);
		checkOuterWriterCleanup(baseline);
		System.out.println("PRINT_TAXONOMY_CRIS_OWNERSHIP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private PrintTaxonomyCrisOwnershipTest(){}

	/** Verifies normal read-column translation remains finite and exact.
	 * @throws Exception If setup, processing, output reading or cleanup fails */
	private static void checkHealthyReadColumn() throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-cris-healthy-");
		try{
			final Path tree=writeTree(root);
			final Path in=root.resolve("healthy.fq");
			final Path out=root.resolve("out.txt");
			writeFastq(in, 2);
			final Throwable failure=process(new PrintTaxonomy(new String[]{
					"tree="+tree,
					"in="+in,
					"out="+out,
					"column=0"}));
			check(failure==null, "healthy read-column control must not fail: "+failure);
			check(text(out).equals("missing0_***NOT_FOUND***\tleft0\n"+
					"missing1_***NOT_FOUND***\tleft1\n"),
					"healthy read-column control must translate exact read IDs");
		}finally{deleteTempTree(root);}
	}

	/** Verifies that an exceptional read exit closes CRIS instead of leaking live threads.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If child setup, execution or cleanup fails */
	private static void checkExceptionalReadCleanup(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-cris-exception-");
		final Path log=root.resolve("child.log");
		Process child=null;
		try{
			final Path marker=root.resolve("npe.marker");
			child=javaProcess(assertionsEnabled(),
					PrintTaxonomyCrisOwnershipTest.class.getName(),
					CHILD_CRIS_MODE, root.toString(), marker.toString())
					.redirectErrorStream(true)
					.redirectOutput(log.toFile())
					.start();
			child.getOutputStream().close();

			check(waitForMarker(marker, CHILD_WAIT_MILLIS),
					"child must prove processReads hit the intended NPE before liveness checks: "+
					text(log));
			System.out.println(CHILD_PROCESS_READS_NPE+": "+text(marker));

			if(baseline){
				Thread.sleep(500);
				check(child.isAlive(),
						"baseline child should still be alive after leaking CRIS threads: "+
						text(log));
				System.out.println("PRINT_TAXONOMY_CRIS_BASELINE_CHILD_ALIVE");
			}else{
				check(child.waitFor(CHILD_WAIT_MILLIS, TimeUnit.MILLISECONDS),
						"fixed child should exit after processReads closes CRIS: "+text(log));
				check(child.exitValue()==0,
						"fixed child should preserve the processReads NPE and exit cleanly: "+
						text(log));
				System.out.println("PRINT_TAXONOMY_CRIS_FIXED_CHILD_EXIT code="+
						child.exitValue());
			}
		}finally{
			terminateChild(child, log);
			deleteTempTree(root);
		}
	}

	/** Verifies that a reported CRIS close error reaches PrintTaxonomy.errorState.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If child setup, execution or cleanup fails */
	private static void checkReportedCrisCloseError(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-cris-close-");
		final Path log=root.resolve("close.log");
		Process child=null;
		try{
			final Path marker=root.resolve("close.marker");
			child=javaProcess(assertionsEnabled(),
					PrintTaxonomyCrisOwnershipTest.class.getName(),
					CHILD_CRIS_CLOSE_MODE, root.toString(), marker.toString(),
					Boolean.toString(baseline))
					.redirectErrorStream(true)
					.redirectOutput(log.toFile())
					.start();
			child.getOutputStream().close();

			check(child.waitFor(CHILD_WAIT_MILLIS, TimeUnit.MILLISECONDS),
					"CRIS close-error child should complete promptly: "+text(log));
			check(child.exitValue()==0,
					"CRIS close-error child should preserve output and close status: "+
					text(log));
			check(Files.exists(marker),
					"CRIS close-error child should write its pass marker: "+text(log));
			System.out.println(CHILD_CRIS_CLOSE+": "+text(marker));
		}finally{
			terminateChild(child, log);
			deleteTempTree(root);
		}
	}

	/** Verifies that an exceptional text path closes the outer TextStreamWriter.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception If child setup, execution or cleanup fails */
	private static void checkOuterWriterCleanup(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("print-taxonomy-writer-exception-");
		final Path log=root.resolve("writer.log");
		Process child=null;
		try{
			final Path marker=root.resolve("writer.marker");
			child=javaProcess(assertionsEnabled(),
					PrintTaxonomyCrisOwnershipTest.class.getName(),
					CHILD_WRITER_MODE, root.toString(), marker.toString())
					.redirectErrorStream(true)
					.redirectOutput(log.toFile())
					.start();
			child.getOutputStream().close();

			check(waitForMarker(marker, CHILD_WAIT_MILLIS),
					"child must prove process caught the parse sentinel before liveness checks: "+
					text(log));
			System.out.println(CHILD_OUTER_WRITER_SENTINEL+": "+text(marker));

			if(baseline){
				Thread.sleep(500);
				check(child.isAlive(),
						"baseline child should still be alive after leaking the writer: "+
						text(log));
				System.out.println("PRINT_TAXONOMY_WRITER_BASELINE_CHILD_ALIVE");
			}else{
				check(child.waitFor(CHILD_WAIT_MILLIS, TimeUnit.MILLISECONDS),
						"fixed child should exit after closing the writer: "+text(log));
				check(child.exitValue()==0,
						"fixed child should preserve the parse sentinel and exit cleanly: "+
						text(log));
				System.out.println("PRINT_TAXONOMY_WRITER_FIXED_CHILD_EXIT code="+
						child.exitValue());
			}
		}finally{
			terminateChild(child, log);
			deleteTempTree(root);
		}
	}

	/** Runs the exceptional processReads child case.
	 * @param args Child-mode arguments: child mode, temp root, marker path
	 * @throws Exception If setup or reflection fails */
	private static void childCrisMain(final String[] args) throws Exception{
		check(args.length==3, "child mode requires root and marker arguments");
		final Path root=Paths.get(args[1]);
		final Path marker=Paths.get(args[2]);
		final Path tree=writeTree(root);
		final Path in=root.resolve("exception.fq");
		final Path out=root.resolve("out.txt");
		writeFastq(in, HANG_READS);
		Shared.setBuffers(HANG_BUFFERS);
		Shared.setBufferLen(HANG_BUFFER_LEN);

		final Throwable failure=processReads(new PrintTaxonomy(new String[]{
				"tree="+tree,
				"in="+in,
				"out="+out,
				"column=0"}), null);
		check(failure instanceof NullPointerException,
				"read body failure must preserve the original writer NPE: "+failure);
		check(topFrame(failure, "tax.PrintTaxonomy", "processReads"),
				"read body failure must originate in processReads: "+failure);
		System.err.println(CHILD_PROCESS_READS_NPE);
		Files.write(marker, CHILD_PROCESS_READS_NPE.getBytes(StandardCharsets.UTF_8));
	}

	/** Runs the reported CRIS close-error child case.
	 * @param args Child-mode arguments: child mode, temp root, marker path, baseline flag
	 * @throws Exception If setup or reflection fails */
	private static void childCrisCloseMain(final String[] args) throws Exception{
		check(args.length==4, "CRIS close child mode requires root, marker and baseline");
		final Path root=Paths.get(args[1]);
		final Path marker=Paths.get(args[2]);
		final boolean baseline=Boolean.parseBoolean(args[3]);
		final Path tree=writeTree(root);
		final Path out=root.resolve("close-out.txt");

		check(!Files.exists(Paths.get(STDIN_CRIS_CLOSE_NAME)),
				"CRIS close-error stdin.* shadow path should not exist: "+
				STDIN_CRIS_CLOSE_NAME);

		final InputStream oldIn=System.in;
		final boolean oldUseUnbgzip=ReadWrite.USE_UNBGZIP;
		final boolean oldUseUnpigz=ReadWrite.USE_UNPIGZ;
		final boolean oldUseGunzip=ReadWrite.USE_GUNZIP;
		final boolean oldAllowNativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF;
		final boolean oldForceBf1=ByteFile.FORCE_MODE_BF1;
		final boolean oldForceBf2=ByteFile.FORCE_MODE_BF2;
		final boolean oldForceBf3=ByteFile.FORCE_MODE_BF3;
		final boolean oldForceBf4=ByteFile.FORCE_MODE_BF4;
		final CloseErrorInputStream closeError=new CloseErrorInputStream(gzipFastq(
				"@missing0\tleft0\nA\n+\n!\n"));
		try{
			ReadWrite.USE_UNBGZIP=false;
			ReadWrite.USE_UNPIGZ=false;
			ReadWrite.USE_GUNZIP=false;
			ReadWrite.ALLOW_NATIVE_BGZF=false;
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF3=false;
			ByteFile.FORCE_MODE_BF4=false;
			System.setIn(closeError);

			final PrintTaxonomy pt=new PrintTaxonomy(new String[]{
					"tree="+tree,
					"in="+STDIN_CRIS_CLOSE_NAME,
					"out="+out,
					"column=0"});
			final TextStreamWriter tsw=new TextStreamWriter(out.toString(), true, false, false);
			tsw.start();
			final boolean[] writerError=new boolean[1];
			final Throwable failure=processReadsAndCloseWriter(pt, tsw, writerError);
			check(failure==null,
					"reported-close processReads body must not fail: "+failure);
			check(!writerError[0], "reported-close output must close cleanly");
			check(text(out).equals("missing0_***NOT_FOUND***\tleft0\n"),
					"reported-close processReads must preserve translated read output");
			check(closeError.closeCount()>0,
					"controlled CRIS input should have been closed; found "+
					closeError.closeCount());
			check(pt.errorState!=baseline,
					"reported CRIS close status must reach PrintTaxonomy.errorState "+
					"only after the fix");
			final String markerText=CHILD_CRIS_CLOSE+" close_count="+
					closeError.closeCount()+" errorState="+pt.errorState;
			System.err.println(markerText);
			Files.write(marker, markerText.getBytes(StandardCharsets.UTF_8));
		}finally{
			System.setIn(oldIn);
			ReadWrite.USE_UNBGZIP=oldUseUnbgzip;
			ReadWrite.USE_UNPIGZ=oldUseUnpigz;
			ReadWrite.USE_GUNZIP=oldUseGunzip;
			ReadWrite.ALLOW_NATIVE_BGZF=oldAllowNativeBgzf;
			ByteFile.FORCE_MODE_BF1=oldForceBf1;
			ByteFile.FORCE_MODE_BF2=oldForceBf2;
			ByteFile.FORCE_MODE_BF3=oldForceBf3;
			ByteFile.FORCE_MODE_BF4=oldForceBf4;
		}
	}

	/** Runs the exceptional outer-writer child case.
	 * @param args Child-mode arguments: child mode, temp root, marker path
	 * @throws Exception If setup or reflection fails */
	private static void childWriterMain(final String[] args) throws Exception{
		check(args.length==3, "writer child mode requires root and marker arguments");
		final Path root=Paths.get(args[1]);
		final Path marker=Paths.get(args[2]);
		final Path tree=writeTree(root);
		final Path in=root.resolve("writer.txt");
		final Path out=root.resolve("out.txt");
		Files.write(in, "missing\n".getBytes(StandardCharsets.UTF_8));

		final Throwable failure=process(new ParseFailingPrintTaxonomy(new String[]{
				"tree="+tree,
				"in="+in,
				"out="+out}));
		check(failure==OUTER_WRITER_SENTINEL,
				"outer writer child must preserve the exact parse sentinel: "+failure);
		System.err.println(CHILD_OUTER_WRITER_SENTINEL);
		Files.write(marker, CHILD_OUTER_WRITER_SENTINEL.getBytes(StandardCharsets.UTF_8));
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

	/** Writes a temporary FASTQ with tabular read identifiers.
	 * @param path FASTQ path
	 * @param reads Number of single-base reads to write
	 * @throws Exception If writing fails */
	private static void writeFastq(final Path path, final int reads) throws Exception{
		final StringBuilder sb=new StringBuilder(reads*32);
		for(int i=0; i<reads; i++){
			sb.append("@missing").append(i).append('\t').append("left").append(i).append('\n');
			sb.append("A\n+\n!\n");
		}
		Files.write(path, sb.toString().getBytes(StandardCharsets.UTF_8));
	}

	/** Gzips a UTF-8 FASTQ payload with the Java gzip implementation.
	 * @param text FASTQ text to compress
	 * @return gzip-compressed bytes
	 * @throws IOException If compression fails */
	private static byte[] gzipFastq(final String text) throws IOException{
		final ByteArrayOutputStream baos=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(baos)){
			gzip.write(text.getBytes(StandardCharsets.UTF_8));
		}
		return baos.toByteArray();
	}

	/** Invokes the package-private PrintTaxonomy.process method.
	 * @param pt Fixture instance to run
	 * @return Thrown failure, or null on success
	 * @throws Exception If reflection itself fails */
	private static Throwable process(final PrintTaxonomy pt) throws Exception{
		final Method process=PrintTaxonomy.class.getDeclaredMethod("process",
				shared.Timer.class);
		process.setAccessible(true);
		try{
			process.invoke(pt, new shared.Timer());
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	/** Invokes the package-private PrintTaxonomy.processReads method.
	 * @param pt Fixture instance to run
	 * @param tsw Output writer, or null to force a body failure
	 * @return Thrown failure, or null on success
	 * @throws Exception If reflection itself fails */
	private static Throwable processReads(final PrintTaxonomy pt, final TextStreamWriter tsw)
			throws Exception{
		final Method processReads=PrintTaxonomy.class.getDeclaredMethod("processReads",
				TextStreamWriter.class);
		processReads.setAccessible(true);
		try{
			processReads.invoke(pt, tsw);
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	/** Invokes processReads and closes the writer while preserving reflection failures.
	 * @param pt Fixture instance to run
	 * @param tsw Output writer to close after invocation
	 * @param writerError Single-value holder for the reported close status
	 * @return Thrown failure, or null on success
	 * @throws Exception If reflection itself fails */
	private static Throwable processReadsAndCloseWriter(final PrintTaxonomy pt,
			final TextStreamWriter tsw, final boolean[] writerError) throws Exception{
		Throwable reflectionFailure=null;
		try{
			return processReads(pt, tsw);
		}catch(Exception | Error e){
			reflectionFailure=e;
			throw e;
		}finally{
			try{
				writerError[0]=tsw.poisonAndWait();
			}catch(RuntimeException | Error e){
				if(reflectionFailure==null){throw e;}
				if(reflectionFailure!=e){reflectionFailure.addSuppressed(e);}
			}
		}
	}

	/** Creates a bounded child JVM with the requested assertion state.
	 * @param assertions true to enable child assertions, false to disable them
	 * @param className fully qualified class name to run
	 * @param args arguments passed after the class name
	 * @return configured child process builder */
	private static ProcessBuilder javaProcess(final boolean assertions, final String className,
			final String... args){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+"/bin/java");
		command.add("-Xmx128m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add(assertions ? "-ea" : "-da");
		if(!"1.8".equals(System.getProperty("java.specification.version"))){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(className);
		command.addAll(Arrays.asList(args));
		return new ProcessBuilder(command);
	}

	/** Polls until a child writes its intended-failure marker.
	 * @param marker Marker path written after the intended child NPE
	 * @param timeoutMillis Maximum wait in milliseconds
	 * @return true if the marker appeared before the timeout
	 * @throws Exception If interrupted while sleeping */
	private static boolean waitForMarker(final Path marker, final long timeoutMillis)
			throws Exception{
		final long deadline=System.nanoTime()+TimeUnit.MILLISECONDS.toNanos(timeoutMillis);
		while(System.nanoTime()<deadline){
			if(Files.exists(marker)){return true;}
			Thread.sleep(50);
		}
		return Files.exists(marker);
	}

	/** Terminates an intentionally leaked baseline child with a bounded forced wait.
	 * @param child Child process, or null if launch failed
	 * @param log Combined stdout/stderr log to report if termination fails
	 * @throws Exception If the child does not terminate after the forced destroy */
	private static void terminateChild(final Process child, final Path log) throws Exception{
		if(child==null || !child.isAlive()){return;}
		child.destroyForcibly();
		check(child.waitFor(CHILD_WAIT_MILLIS, TimeUnit.MILLISECONDS),
				"child should terminate after forced destroy: "+text(log));
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

	/** Returns whether the current JVM runs with assertions enabled.
	 * @return true if assertions are enabled */
	private static boolean assertionsEnabled(){
		boolean enabled=false;
		assert(enabled=true);
		return enabled;
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

	/** ByteArrayInputStream wrapper that reports a real IOException from close(). */
	private static final class CloseErrorInputStream extends InputStream{

		/** Creates a stream over one in-memory gzipped FASTQ.
		 * @param data bytes consumed by CRIS before close */
		CloseErrorInputStream(final byte[] data){
			delegate=new ByteArrayInputStream(data);
		}

		@Override
		public int read() throws IOException{
			return delegate.read();
		}

		@Override
		public int read(final byte[] b, final int off, final int len) throws IOException{
			return delegate.read(b, off, len);
		}

		@Override
		public void close() throws IOException{
			closeCount++;
			delegate.close();
			throw new IOException("artificial CRIS raw close failure");
		}

		/** Reports how often cleanup reached the controlled CRIS input.
		 * @return number of close attempts observed */
		int closeCount(){return closeCount;}

		/** Backing stream containing the controlled gzipped FASTQ. */
		private final ByteArrayInputStream delegate;
		/** Number of close calls from CRIS's reader chain. */
		private int closeCount;
	}

	/** PrintTaxonomy variant that fails after the outer writer is started and used. */
	private static final class ParseFailingPrintTaxonomy extends PrintTaxonomy{

		/** Creates a PrintTaxonomy instance with the real constructor.
		 * @param args Real PrintTaxonomy arguments */
		ParseFailingPrintTaxonomy(final String[] args){super(args);}

		/** Throws a sentinel when taxonomy resolution reaches header parsing.
		 * @param header Header routed through the real text-input process
		 * @return Never returns normally */
		@Override
		public TaxNode parseNodeFromHeader(final String header){
			throw OUTER_WRITER_SENTINEL;
		}
	}

}
