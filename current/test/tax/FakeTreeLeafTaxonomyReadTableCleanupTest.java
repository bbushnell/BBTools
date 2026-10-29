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
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPOutputStream;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import tax.FakeTreeLeafTaxonomy;

/** Baseline/fixed fixture for FakeTreeLeafTaxonomy.readTable exceptional input cleanup.
 * Baseline mode requires skipped close failures; fixed mode requires one close,
 * the same primary failures, and suppressed close-status failures.
 * @author Jean */
public final class FakeTreeLeafTaxonomyReadTableCleanupTest{

	/** Runs parent or child checks for readTable cleanup.
	 * @param args Optional expect=baseline|fixed or child-mode arguments
	 * @throws Exception if fixture setup, child execution or cleanup fails */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && CHILD_EMPTY.equals(args[0])){
			if(checkEmptyChild(Boolean.parseBoolean(args[1]), Paths.get(args[2]))){
				System.exit(0);
			}
			return;
		}
		if(args.length>0 && CHILD_UNCHECKED.equals(args[0])){
			checkUncheckedChild(Boolean.parseBoolean(args[1]), Paths.get(args[2]));
			return;
		}

		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		boolean success=false;
		final Path root=Files.createTempDirectory("faketree-readtable-cleanup-");
		try{
			checkHealthyTable(root);
			runParentCase(root, CHILD_EMPTY, baseline, false);
			runParentCase(root, CHILD_EMPTY, baseline, true);
			runParentCase(root, CHILD_UNCHECKED, baseline, false);
			runParentCase(root, CHILD_UNCHECKED, baseline, true);
			success=true;
		}finally{
			if(success){deleteTree(root);}
			else{System.err.println("Retained FakeTree readTable child logs in "+root);}
		}
		System.out.println("FAKETREE_READTABLE_CLEANUP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private FakeTreeLeafTaxonomyReadTableCleanupTest(){}

	/** Verifies normal table parsing still returns a row with the expected data.
	 * @param root Parent temporary directory for the healthy table
	 * @throws Exception if reflection or table parsing fails */
	private static void checkHealthyTable(final Path root) throws Exception{
		final Path table=root.resolve("healthy.tsv");
		Files.write(table, ("fake_taxid\tvalue\n17\talpha\n").getBytes(StandardCharsets.UTF_8));

		final List<?> rows=invokeReadTable(fileFormat(table.toString()));
		check(rows.size()==1, "Healthy table should produce one row, not "+rows.size());
		check("17".equals(rowValue(rows.get(0), "fake_taxid")),
				"Healthy row should retain fake_taxid=17");
		check("alpha".equals(rowValue(rows.get(0), "value")),
				"Healthy row should retain value=alpha");
	}

	/** Runs one bounded child JVM and verifies its marker.
	 * @param root Parent-owned temporary directory for child logs
	 * @param childMode Child case selector
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @param assertions Whether to enable child assertions
	 * @throws Exception if the child cannot be launched, read or verified */
	private static void runParentCase(final Path root, final String childMode,
			final boolean baseline, final boolean assertions) throws Exception{
		final String label=childMode+(assertions ? "-ea" : "-da");
		final Path log=root.resolve(label+".log");
		final Path marker=root.resolve(label+".marker");
		Process child=null;
		final long started=System.nanoTime();
		try{
			final ProcessBuilder pb=javaProcess(assertions, childMode,
					Boolean.toString(baseline), marker.toString());
			pb.redirectErrorStream(true);
			pb.redirectOutput(log.toFile());
			child=pb.start();
			child.getOutputStream().close();
			if(!child.waitFor(CHILD_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
				terminateChild(child, label, log);
				throw new AssertionError("Timed out waiting for "+label+" child: "+
						readAscii(log));
			}

			final String output=readAscii(log);
			final long elapsedMs=TimeUnit.NANOSECONDS.toMillis(System.nanoTime()-started);
			System.out.println(label+" child exit="+child.exitValue()+
					" elapsed_ms="+elapsedMs+" output='"+compact(output)+"'");
			check(child.exitValue()==0,
					"Child "+label+" failed with exit "+child.exitValue()+": "+output);
			check(Files.exists(marker), "Child "+label+" did not write "+marker+
					": "+output);
			final String receipt=readAscii(marker);
			check(output.contains(receipt),
					"Child "+label+" did not print its marker before exit: "+output);
		}finally{
			if(child!=null && child.isAlive()){terminateChild(child, label, log);}
		}
	}

	/** Exercises the empty-table pre-close escape with a BF4 Java-gzip stdin reader.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @param marker Path written after the oracle has passed
	 * @return True only for the baseline child that must exit explicitly after receipt
	 * @throws Exception if the oracle fails */
	private static boolean checkEmptyChild(final boolean baseline, final Path marker)
			throws Exception{
		check(!Files.exists(Paths.get(STDIN_EMPTY_NAME)),
				"Empty stdin.*.gz shadow path should not exist: "+STDIN_EMPTY_NAME);
		final CloseErrorInputStream controlledIn=new CloseErrorInputStream(gzip(""), false);
		final ReaderFlags flags=new ReaderFlags();
		try{
			forceJavaGzipReader(4, controlledIn);
			final RuntimeException failure=expectFailure(fileFormat(STDIN_EMPTY_NAME));
			check(("Empty table: "+STDIN_EMPTY_NAME).equals(failure.getMessage()),
					"Empty table should remain primary: "+failure);
			assertSuppressedClose(failure, baseline, STDIN_EMPTY_NAME);
			check(controlledIn.closeCount()==(baseline ? 0 : 1),
					"Unexpected empty-table close count: "+controlledIn.closeCount());

			final String receipt="case=empty close_count="+controlledIn.closeCount()+
					" explicit_exit="+baseline;
			writeReceipt(marker, receipt);
			return baseline;
		}finally{
			flags.restore();
		}
	}

	/** Exercises unchecked read failure after Java-gzip reader acquisition.
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @param marker Path written after the oracle has passed
	 * @throws Exception if the oracle fails */
	private static void checkUncheckedChild(final boolean baseline, final Path marker)
			throws Exception{
		check(!Files.exists(Paths.get(STDIN_UNCHECKED_NAME)),
				"Unchecked stdin.*.gz shadow path should not exist: "+STDIN_UNCHECKED_NAME);
		final CloseErrorInputStream controlledIn=new CloseErrorInputStream(
				gzip("fake_taxid\n17\n"), true);
		final ReaderFlags flags=new ReaderFlags();
		try{
			forceJavaGzipReader(1, controlledIn);
			final RuntimeException failure=expectFailure(fileFormat(STDIN_UNCHECKED_NAME));
			check(failure==READ_SENTINEL,
					"Unchecked read failure should preserve the exact sentinel: "+failure);
			check(controlledIn.sawByteFile1DataRead(),
					"Unchecked sentinel should come from ByteFile1 fillBuffer");
			assertSuppressedClose(failure, baseline, STDIN_UNCHECKED_NAME);
			check(controlledIn.closeCount()==(baseline ? 0 : 1),
					"Unexpected unchecked close count: "+controlledIn.closeCount());

			writeReceipt(marker, "case=unchecked close_count="+controlledIn.closeCount());
		}finally{
			flags.restore();
		}
	}

	/** Configures a child to use Java gzip and one concrete ByteFile implementation.
	 * @param bfMode ByteFile implementation to force, 1 or 4
	 * @param in Replacement standard input */
	private static void forceJavaGzipReader(final int bfMode, final InputStream in){
		System.setIn(in);
		ReadWrite.USE_UNBGZIP=false;
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_GUNZIP=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		ByteFile.FORCE_MODE_BF1=bfMode==1;
		ByteFile.FORCE_MODE_BF2=false;
		ByteFile.FORCE_MODE_BF3=false;
		ByteFile.FORCE_MODE_BF4=bfMode==4;
	}

	/** Requires one fixed-mode suppressed close failure and none in baseline mode.
	 * @param failure Primary RuntimeException from readTable
	 * @param baseline Whether baseline mode should suppress nothing
	 * @param name Input name expected in the close-status message */
	private static void assertSuppressedClose(final RuntimeException failure,
			final boolean baseline, final String name){
		final Throwable[] suppressed=failure.getSuppressed();
		if(baseline){
			check(suppressed.length==0,
					"Baseline should not suppress a skipped close failure");
		}else{
			check(suppressed.length==1,
					"Fixed mode should suppress one close failure");
			check(suppressed[0] instanceof RuntimeException,
					"Suppressed close failure should be RuntimeException: "+suppressed[0]);
			check(("Error reading "+name).equals(suppressed[0].getMessage()),
					"Suppressed close failure should preserve readTable wording: "+
					suppressed[0]);
		}
	}

	/** Invokes readTable and requires a RuntimeException.
	 * @param ff Input FileFormat
	 * @return RuntimeException thrown by readTable
	 * @throws Exception if reflection fails or readTable does not throw */
	private static RuntimeException expectFailure(final FileFormat ff) throws Exception{
		try{
			invokeReadTable(ff);
			throw new AssertionError("readTable should have thrown for "+ff.name());
		}catch(RuntimeException expected){
			return expected;
		}
	}

	/** Invokes the private readTable helper.
	 * @param ff Input FileFormat
	 * @return Parsed rows
	 * @throws Exception if reflection fails */
	private static List<?> invokeReadTable(final FileFormat ff) throws Exception{
		final Method method=FakeTreeLeafTaxonomy.class.getDeclaredMethod("readTable",
				FileFormat.class);
		method.setAccessible(true);
		try{
			return (List<?>) method.invoke(null, ff);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof RuntimeException){throw (RuntimeException) cause;}
			if(cause instanceof Error){throw (Error) cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Reads one value from a private FakeTreeLeafTaxonomy Row.
	 * @param row Row object returned by readTable
	 * @param key Column key
	 * @return Stored value for key
	 * @throws Exception if reflection fails */
	private static String rowValue(final Object row, final String key) throws Exception{
		final Method method=row.getClass().getDeclaredMethod("get", String.class);
		method.setAccessible(true);
		return (String) method.invoke(row, key);
	}

	/** Builds an input FileFormat matching FakeTreeLeafTaxonomy's constructor.
	 * @param name Input name
	 * @return Text FileFormat for name */
	private static FileFormat fileFormat(final String name){
		return FileFormat.testInput(name, FileFormat.TEXT, null, true, true);
	}

	/** Writes and prints one child receipt.
	 * @param marker Destination marker path
	 * @param receipt Exact receipt text
	 * @throws IOException if the marker cannot be written */
	private static void writeReceipt(final Path marker, final String receipt) throws IOException{
		System.out.println(receipt);
		Files.write(marker, receipt.getBytes(StandardCharsets.UTF_8));
	}

	/** Creates a java command for a child-mode invocation of this fixture.
	 * @param assertions Whether to enable child assertions
	 * @param args Child-mode arguments
	 * @return Configured ProcessBuilder */
	private static ProcessBuilder javaProcess(final boolean assertions, final String... args){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add("-Xmx256m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		command.add(assertions ? "-ea" : "-da");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(FakeTreeLeafTaxonomyReadTableCleanupTest.class.getName());
		command.addAll(Arrays.asList(args));
		return new ProcessBuilder(command);
	}

	/** Adds the Vector incubator module when this JDK exposes it.
	 * @param command Child JVM command to extend */
	private static void addVectorModule(final List<String> command){
		try{
			final Class<?> moduleLayer=Class.forName("java.lang.ModuleLayer");
			final Object boot=moduleLayer.getMethod("boot").invoke(null);
			final Object optional=moduleLayer.getMethod("findModule",
					String.class).invoke(boot, "jdk.incubator.vector");
			final boolean present=((Boolean)optional.getClass()
					.getMethod("isPresent").invoke(optional)).booleanValue();
			if(present){
				command.add("--add-modules");
				command.add("jdk.incubator.vector");
			}
		}catch(ClassNotFoundException e){
			return;
		}catch(Exception e){
			throw new RuntimeException("Could not probe jdk.incubator.vector", e);
		}
	}

	/** Terminates an overdue child and waits briefly for process cleanup.
	 * @param child Child process to terminate
	 * @param label Child label for diagnostics
	 * @param log Child output log
	 * @throws Exception if the timeout diagnostic cannot be read */
	private static void terminateChild(final Process child, final String label,
			final Path log) throws Exception{
		child.destroyForcibly();
		if(child.waitFor(5, TimeUnit.SECONDS)){return;}
		throw new AssertionError("Could not terminate "+label+" child; log="+readAscii(log));
	}

	/** Gzips one UTF-8 table payload.
	 * @param text Table text to compress
	 * @return gzip-compressed bytes
	 * @throws IOException if compression fails */
	private static byte[] gzip(final String text) throws IOException{
		final ByteArrayOutputStream baos=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(baos)){
			gzip.write(text.getBytes(StandardCharsets.UTF_8));
		}
		return baos.toByteArray();
	}

	/** Reads one log or marker as ASCII.
	 * @param path Text path
	 * @return Complete text, or empty for an absent file
	 * @throws IOException if the file cannot be read */
	private static String readAscii(final Path path) throws IOException{
		return Files.exists(path) ? new String(Files.readAllBytes(path),
				StandardCharsets.US_ASCII) : "";
	}

	/** Compacts child output to one line for parent diagnostics.
	 * @param text Raw child output
	 * @return One-line output */
	private static String compact(final String text){
		return text.replace('\n', '|').replace('\r', '~');
	}

	/** Deletes a temporary tree in reverse path order.
	 * @param root Tree root
	 * @throws IOException if a path cannot be deleted */
	private static void deleteTree(final Path root) throws IOException{
		try(Stream<Path> stream=Files.walk(root)){
			final List<Path> paths=stream.sorted(Comparator.reverseOrder())
					.collect(Collectors.<Path>toList());
			for(Path path : paths){Files.deleteIfExists(path);}
		}
	}

	/** Test expectations stay fatal with assertions disabled.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Restores process-wide input-reader flags after one child case. */
	private static final class ReaderFlags{
		/** Saves all globals this fixture mutates. */
		ReaderFlags(){
			in=System.in;
			useUnbgzip=ReadWrite.USE_UNBGZIP;
			useUnpigz=ReadWrite.USE_UNPIGZ;
			useGunzip=ReadWrite.USE_GUNZIP;
			allowNativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF;
			forceBf1=ByteFile.FORCE_MODE_BF1;
			forceBf2=ByteFile.FORCE_MODE_BF2;
			forceBf3=ByteFile.FORCE_MODE_BF3;
			forceBf4=ByteFile.FORCE_MODE_BF4;
		}

		/** Restores saved globals. */
		void restore(){
			System.setIn(in);
			ReadWrite.USE_UNBGZIP=useUnbgzip;
			ReadWrite.USE_UNPIGZ=useUnpigz;
			ReadWrite.USE_GUNZIP=useGunzip;
			ReadWrite.ALLOW_NATIVE_BGZF=allowNativeBgzf;
			ByteFile.FORCE_MODE_BF1=forceBf1;
			ByteFile.FORCE_MODE_BF2=forceBf2;
			ByteFile.FORCE_MODE_BF3=forceBf3;
			ByteFile.FORCE_MODE_BF4=forceBf4;
		}

		/** Original standard input. */
		private final InputStream in;
		/** Original native BGZF Java fallback gate. */
		private final boolean useUnbgzip;
		/** Original pigz input gate. */
		private final boolean useUnpigz;
		/** Original external gunzip input gate. */
		private final boolean useGunzip;
		/** Original native BGZF input permission. */
		private final boolean allowNativeBgzf;
		/** Original BF1 force flag. */
		private final boolean forceBf1;
		/** Original BF2 force flag. */
		private final boolean forceBf2;
		/** Original BF3 force flag. */
		private final boolean forceBf3;
		/** Original BF4 force flag. */
		private final boolean forceBf4;
	}

	/** ByteArrayInputStream wrapper with optional post-acquisition read failure. */
	private static final class CloseErrorInputStream extends InputStream{

		/** Wraps gzip-compressed bytes with a controllable read and close failure.
		 * @param data gzip-compressed input
		 * @param throwAfterGzip Whether array reads after GZIP construction should fail */
		CloseErrorInputStream(final byte[] data, final boolean throwAfterGzip){
			delegate=new ByteArrayInputStream(data);
			throwAfterGzip_=throwAfterGzip;
		}

		@Override
		public int read() throws IOException{
			noteGzipAcquisition();
			return delegate.read();
		}

		@Override
		public int read(final byte[] b, final int off, final int len) throws IOException{
			noteGzipAcquisition();
			if(throwAfterGzip_ && byteFile1DataRead()){
				check(sawByteFile1GzipAcquire,
						"Unchecked ByteFile1 data read fired before ByteFile1 GZIP acquisition");
				sawByteFile1DataRead=true;
				throw READ_SENTINEL;
			}
			return delegate.read(b, off, len);
		}

		@Override
		public void close() throws IOException{
			closeCount++;
			delegate.close();
			throw new IOException("TAX-034 controlled close failure");
		}

		/** Reports close attempts from the Java GZIP wrapper.
		 * @return number of close calls */
		int closeCount(){return closeCount;}

		/** Reports whether the sentinel was thrown from a ByteFile1 data read.
		 * @return true if ByteFile1.fillBuffer attempted a post-acquisition read */
		boolean sawByteFile1DataRead(){return sawByteFile1DataRead;}

		/** Records gzip construction by ByteFile1.open after FileFormat probing.
		 * FileFormat's gzip probes must not satisfy the acquisition proof. */
		private void noteGzipAcquisition(){
			if(acquiringGzip() && byteFile1FactoryAcquire()){
				sawByteFile1GzipAcquire=true;
			}
		}

		/** Detects reads performed during GZIPInputStream header acquisition.
		 * @return True if the current stack is constructing a Java GZIPInputStream */
		private static boolean acquiringGzip(){
			for(StackTraceElement element : Thread.currentThread().getStackTrace()){
				if("java.util.zip.GZIPInputStream".equals(element.getClassName())){
					final String method=element.getMethodName();
					if("<init>".equals(method) || "readHeader".equals(method)){
						return true;
					}
				}
			}
			return false;
		}

		/** Detects reader acquisition inside ByteFile1's constructor path.
		 * @return True if the current stack is inside ByteFile1.open */
		private static boolean byteFile1FactoryAcquire(){
			for(StackTraceElement element : Thread.currentThread().getStackTrace()){
				if("fileIO.ByteFile1".equals(element.getClassName()) &&
						"open".equals(element.getMethodName())){
					return true;
				}
			}
			return false;
		}

		/** Detects ByteFile1 data reads after the factory has acquired the reader.
		 * @return True if the current stack is reading data through ByteFile1 */
		private static boolean byteFile1DataRead(){
			for(StackTraceElement element : Thread.currentThread().getStackTrace()){
				if("fileIO.ByteFile1".equals(element.getClassName()) &&
						"fillBuffer".equals(element.getMethodName())){
					return true;
				}
			}
			return false;
		}

		/** Gzipped data source. */
		private final ByteArrayInputStream delegate;
		/** Whether the first post-acquisition bulk read should throw READ_SENTINEL. */
		private final boolean throwAfterGzip_;
		/** Whether ByteFile1.open constructed Java GZIP before a sentinel was thrown. */
		private boolean sawByteFile1GzipAcquire;
		/** Whether READ_SENTINEL was thrown from ByteFile1.fillBuffer. */
		private boolean sawByteFile1DataRead;
		/** Number of close calls from the Java GZIP wrapper. */
		private int closeCount;
	}

	/** Child mode for the empty-table BF4 case. */
	private static final String CHILD_EMPTY="child-empty";
	/** Child mode for the unchecked-read BF1 case. */
	private static final String CHILD_UNCHECKED="child-unchecked";
	/** Nonexistent Java-gzip stdin path for the empty-table case. */
	private static final String STDIN_EMPTY_NAME="stdin.faketree-empty.tsv.gz";
	/** Nonexistent Java-gzip stdin path for the unchecked-read case. */
	private static final String STDIN_UNCHECKED_NAME="stdin.faketree-unchecked.tsv.gz";
	/** Milliseconds to wait for child completion before forced termination. */
	private static final long CHILD_TIMEOUT_MS=30000;
	/** Exact unchecked sentinel thrown from ByteFile1 after GZIPInputStream construction. */
	private static final RuntimeException READ_SENTINEL=
			new RuntimeException("TAX-034 controlled read failure");
}
