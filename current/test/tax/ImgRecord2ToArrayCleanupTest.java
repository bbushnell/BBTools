package test.tax;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;

import tax.ImgRecord2;

/** Finite child-JVM descriptor cleanup check for ImgRecord2.toArray.
 * @author Sucrose */
public final class ImgRecord2ToArrayCleanupTest{

	/** Runs parent or child checks for ImgRecord2.toArray cleanup.
	 * @param args Optional child-mode or baseline/fixed expectation arguments
	 * @throws Exception if fixture setup, child execution, or cleanup fails */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			runChildCase(args[1], Paths.get(args[2]), Boolean.parseBoolean(args[3]));
			return;
		}

		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("imgrecord2-toarray-cleanup-");
		try{
			final Path valid=root.resolve("valid.tsv");
			final Path skipped=root.resolve("skipped.tsv");
			final Path malformed=root.resolve("malformed.tsv");
			Files.write(valid, VALID_LINES, StandardCharsets.UTF_8);
			Files.write(skipped, SKIPPED_LINES, StandardCharsets.UTF_8);
			Files.write(malformed, MALFORMED_LINES, StandardCharsets.UTF_8);

			runParentCase(root, "valid", valid, baseline, false);
			runParentCase(root, "valid", valid, baseline, true);
			runParentCase(root, "skipped", skipped, baseline, false);
			runParentCase(root, "skipped", skipped, baseline, true);
			runParentCase(root, "malformed", malformed, baseline, false);
			runParentCase(root, "malformed", malformed, baseline, true);
			runParentCase(root, "close-error", Paths.get(STDIN_CLOSE_ERROR_NAME),
					baseline, false);
			runParentCase(root, "close-error", Paths.get(STDIN_CLOSE_ERROR_NAME),
					baseline, true);
		}finally{
			deleteTree(root);
		}
		System.out.println("IMGRECORD2_TOARRAY_CLEANUP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private ImgRecord2ToArrayCleanupTest(){}

	/** Runs one bounded child JVM and validates its per-case oracle.
	 * @param root Parent-owned temporary directory for the child log
	 * @param caseName Child case selector
	 * @param input Child input file path
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @param assertions Whether to enable child assertions
	 * @throws Exception if the child cannot be launched, read, or verified */
	private static void runParentCase(final Path root, final String caseName,
			final Path input, final boolean baseline, final boolean assertions) throws Exception{
		final String label=caseName+(assertions ? "-ea" : "-da");
		final Path log=root.resolve(label+".log");
		Process child=null;
		final long started=System.nanoTime();
		try{
			final ProcessBuilder pb=javaProcess(ImgRecord2ToArrayCleanupTest.class.getName(),
					assertions, "child", caseName, input.toString(), Boolean.toString(baseline));
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
			check(output.contains("case="+caseName),
					"Child "+label+" did not report its case marker: "+output);
		}finally{
			if(child!=null && child.isAlive()){terminateChild(child, label, log);}
		}
	}

	/** Runs one descriptor-counting toArray case inside a fresh child JVM.
	 * @param caseName Case selector
	 * @param input Input file path to parse
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception if parsing or descriptor inspection diverges from the oracle */
	private static void runChildCase(final String caseName, final Path input,
			final boolean baseline) throws Exception{
		final boolean oldStoreName=ImgRecord2.storeName;
		try{
			ImgRecord2.storeName=true;
			if("close-error".equals(caseName)){
				checkCloseErrorStatus(input, baseline);
				return;
			}

			final int before=countOpenDescriptors(input);
			if("valid".equals(caseName)){
				final ImgRecord2[] records=ImgRecord2.toArray(input.toString());
				check(records.length==2, "Valid input should emit two records");
				checkRecord(records[0], 11, 101, "Genome Alpha");
				checkRecord(records[1], 22, 202, "Genome Beta");
				checkNoLeak(caseName, input, before);
			}else if("skipped".equals(caseName)){
				final ImgRecord2[] records=ImgRecord2.toArray(input.toString());
				check(records.length==1, "Non-digit malformed row should be skipped");
				checkRecord(records[0], 33, 303, "Genome Gamma");
				checkNoLeak(caseName, input, before);
			}else if("malformed".equals(caseName)){
				checkMalformedFailure(input, before, baseline);
			}else{
				throw new IllegalArgumentException("Unknown child case: "+caseName);
			}
		}finally{
			ImgRecord2.storeName=oldStoreName;
		}
	}

	/** Requires a malformed digit-starting row to keep the primary NumberFormatException.
	 * @param input Input file with a bad taxID in the second record
	 * @param before Target-file descriptor count before parsing
	 * @param baseline Whether the parser is expected to leak the input descriptor
	 * @throws Exception if descriptor inspection fails */
	private static void checkMalformedFailure(final Path input, final int before,
			final boolean baseline) throws Exception{
		try{
			ImgRecord2.toArray(input.toString());
			throw new AssertionError("Malformed taxID should throw NumberFormatException");
		}catch(NumberFormatException expected){
			check("For input string: \"bad\"".equals(expected.getMessage()),
					"Original NumberFormatException should remain primary: "+expected);
		}
		final int after=countOpenDescriptors(input);
		if(baseline){
			check(after>before,
					"Baseline should leak the malformed input; before="+before+
					" after="+after);
		}else{
			check(after==before,
					"Fixed parser should close malformed input; before="+before+
					" after="+after);
		}
		System.out.println("case=malformed descriptor_delta="+(after-before));
	}

	/** Requires the TextFile true close-error status to become a rejection.
	 * @param input Nonexistent stdin shadow path selected by ReadWrite's stdin handling
	 * @param baseline Whether the parser is expected to ignore the close-error status
	 * @throws Exception if the controlled stdin path is shadowed or not fully closed */
	private static void checkCloseErrorStatus(final Path input, final boolean baseline)
			throws Exception{
		check(STDIN_CLOSE_ERROR_NAME.equals(input.toString()),
				"Close-error case must use the stdin.* shadow path: "+input);
		check(!Files.exists(input),
				"Close-error stdin.* shadow path should not exist: "+input);

		final InputStream oldIn=System.in;
		final CloseErrorInputStream controlledIn=new CloseErrorInputStream(
				CLOSE_ERROR_LINE.getBytes(StandardCharsets.UTF_8));
		System.setIn(controlledIn);
		try{
			if(baseline){
				final ImgRecord2[] records=ImgRecord2.toArray(STDIN_CLOSE_ERROR_NAME);
				check(records.length==1,
						"Baseline close-error input should emit one record");
				checkRecord(records[0], 66, 606, "Genome Zeta");
			}else{
				try{
					ImgRecord2.toArray(STDIN_CLOSE_ERROR_NAME);
					throw new AssertionError(
							"Fixed close-error input should throw RuntimeException");
				}catch(RuntimeException expected){
					final String message="Error closing "+STDIN_CLOSE_ERROR_NAME;
					check(message.equals(expected.getMessage()),
							"Close-error message should name the input; expected "
							+message+" but found "+expected);
				}
			}
		}finally{
			System.setIn(oldIn);
		}

		check(controlledIn.eofReached(),
				"Controlled stdin should have been consumed to EOF");
		check(controlledIn.closeCount()>0,
				"Controlled stdin should have been closed; found "
				+controlledIn.closeCount());
		System.out.println("case=close-error eof=true close_count="
				+controlledIn.closeCount());
	}

	/** ByteArrayInputStream wrapper that reports a real IOException from close(). */
	private static final class CloseErrorInputStream extends InputStream{

		/** Creates a stream over one in-memory fixture payload.
		 * @param data bytes consumed by TextFile before close */
		CloseErrorInputStream(final byte[] data){
			delegate=new ByteArrayInputStream(data);
		}

		@Override
		public int read() throws IOException{
			final int value=delegate.read();
			if(value<0){eof=true;}
			return value;
		}

		@Override
		public int read(final byte[] b, final int off, final int len) throws IOException{
			final int count=delegate.read(b, off, len);
			if(count<0){eof=true;}
			return count;
		}

		@Override
		public void close() throws IOException{
			closeCount++;
			delegate.close();
			throw new IOException(CLOSE_ERROR_MESSAGE);
		}

		/** Reports whether the reader consumed the full input.
		 * @return whether TextFile read through EOF before closing */
		boolean eofReached(){return eof;}

		/** Reports how often cleanup reached the controlled input.
		 * @return number of close attempts observed */
		int closeCount(){return closeCount;}

		/** Backing stream containing the controlled IMG row. */
		private final ByteArrayInputStream delegate;
		/** Whether a read observed EOF. */
		private boolean eof;
		/** Number of close calls from TextFile's reader chain. */
		private int closeCount;
	}

	/** Verifies one parsed record's exact fields.
	 * @param record Parsed ImgRecord2 instance
	 * @param imgID Expected IMG identifier
	 * @param taxID Expected NCBI taxonomy identifier
	 * @param name Expected retained name */
	private static void checkRecord(final ImgRecord2 record, final long imgID,
			final int taxID, final String name){
		check(record.imgID==imgID,
				"Expected imgID "+imgID+" but found "+record.imgID);
		check(record.taxID==taxID,
				"Expected taxID "+taxID+" but found "+record.taxID);
		check(name.equals(record.name),
				"Expected name "+name+" but found "+record.name);
	}

	/** Requires the parsed input to have no retained target-file descriptors.
	 * @param caseName Case name for diagnostics
	 * @param input Input file to inspect
	 * @param before Target-file descriptor count before parsing
	 * @throws IOException if descriptor inspection fails */
	private static void checkNoLeak(final String caseName, final Path input, final int before)
			throws IOException{
		final int after=countOpenDescriptors(input);
		check(after==before,
				caseName+" case should not retain input descriptors; before="+before+
				" after="+after);
		System.out.println("case="+caseName+" descriptor_delta="+(after-before));
	}

	/** Counts this JVM's open descriptors for exactly one fixture path.
	 * @param target Fixture file to count
	 * @return Number of /proc/self/fd entries pointing to target
	 * @throws IOException if descriptor inspection is unavailable */
	private static int countOpenDescriptors(final Path target) throws IOException{
		check(Files.isDirectory(FD_ROOT),
				"ImgRecord2ToArrayCleanupTest requires /proc/self/fd");
		final String canonical=target.toRealPath().toString();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : fds){
				try{
					if(Files.readSymbolicLink(fd).toString().equals(canonical)){count++;}
				}catch(IOException ignored){
					// Descriptor disappeared between directory listing and readlink.
				}
			}
		}
		return count;
	}

	/** Creates a bounded child JVM command.
	 * @param className Fully qualified class name to run
	 * @param assertions Whether to enable child assertions
	 * @param args Arguments passed after the class name
	 * @return configured process builder */
	private static ProcessBuilder javaProcess(final String className, final boolean assertions,
			final String... args){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+File.separator+"bin"+
				File.separator+"java");
		command.add(assertions ? "-ea" : "-da");
		command.add("-Xmx128m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(className);
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

	/** Terminates a child JVM with bounded graceful and forced waits.
	 * @param child Child process to terminate
	 * @param label Child case label for diagnostics
	 * @param log Child log retained for diagnostics
	 * @throws Exception if the child cannot be terminated or its log cannot be read */
	private static void terminateChild(final Process child, final String label,
			final Path log) throws Exception{
		child.destroy();
		if(child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){return;}

		child.destroyForcibly();
		if(child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){return;}

		throw new AssertionError("Could not terminate "+label+" child: "
				+readAscii(log));
	}

	/** Reads one small ASCII child log.
	 * @param path Child log path
	 * @return Full child log
	 * @throws IOException if the log cannot be read */
	private static String readAscii(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Compacts child output into one retained parent line.
	 * @param text Full child output
	 * @return one-line output for parent diagnostics */
	private static String compact(final String text){
		return text.replace('\n', '|').replace('\r', '|');
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Owned temporary directory root
	 * @throws Exception if the tree cannot be walked or deleted */
	private static void deleteTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps checks active with assertions disabled.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Time to wait for each bounded child JVM. */
	private static final long CHILD_TIMEOUT_MS=30000;

	/** Time to wait for graceful child termination after a timeout or parent failure. */
	private static final long DESTROY_TIMEOUT_MS=5000;

	/** Linux descriptor directory for this JVM only. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");

	/** Fully valid two-record IMG table with exact expected field values. */
	private static final java.util.List<String> VALID_LINES=Arrays.asList(
			"11\t101\tGenome Alpha",
			"22\t202\tGenome Beta");

	/** Non-digit malformed row followed by a valid row. */
	private static final java.util.List<String> SKIPPED_LINES=Arrays.asList(
			"bad\ttax\tname",
			"33\t303\tGenome Gamma");

	/** Valid first row followed by a digit-starting row with a malformed taxID. */
	private static final java.util.List<String> MALFORMED_LINES=Arrays.asList(
			"44\t404\tGenome Delta",
			"55\tbad\tGenome Epsilon");

	/** Nonexistent path that ReadWrite maps to System.in only when unshadowed. */
	private static final String STDIN_CLOSE_ERROR_NAME="stdin.imgrecord2-close-error";

	/** Valid IMG payload supplied through controlled stdin. */
	private static final String CLOSE_ERROR_LINE="66\t606\tGenome Zeta\n";

	/** Sentinel text emitted by the controlled stdin close failure. */
	private static final String CLOSE_ERROR_MESSAGE="TAX-049 controlled close failure";
}
