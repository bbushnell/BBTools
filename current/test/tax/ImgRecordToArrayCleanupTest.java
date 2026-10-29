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

import tax.ImgRecord;

/** Finite child-JVM cleanup check for ImgRecord.toArray baseline and fixed cases.
 * @author Sucrose */
public final class ImgRecordToArrayCleanupTest{

	/** Runs parent or child checks for ImgRecord.toArray cleanup.
	 * @param args Optional expect=baseline|fixed or child-mode arguments
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

		boolean success=false;
		final Path root=Files.createTempDirectory("imgrecord-toarray-cleanup-");
		try{
			final Path normal=root.resolve("normal.tsv");
			final Path malformed=root.resolve("malformed.tsv");
			Files.write(normal, NORMAL_LINES, StandardCharsets.UTF_8);
			Files.write(malformed, MALFORMED_LINES, StandardCharsets.UTF_8);

			runParentCase(root, "normal", normal, baseline, false);
			runParentCase(root, "normal", normal, baseline, true);
			runParentCase(root, "hq", normal, baseline, false);
			runParentCase(root, "hq", normal, baseline, true);
			runParentCase(root, "malformed", malformed, baseline, false);
			runParentCase(root, "malformed", malformed, baseline, true);
			runParentCase(root, "close", Paths.get(STDIN_CLOSE_ERROR_NAME),
					baseline, false);
			runParentCase(root, "close", Paths.get(STDIN_CLOSE_ERROR_NAME),
					baseline, true);
			if(!baseline){
				runParentCase(root, "combined", Paths.get(STDIN_COMBINED_ERROR_NAME),
						baseline, false);
				runParentCase(root, "combined", Paths.get(STDIN_COMBINED_ERROR_NAME),
						baseline, true);
			}
			success=true;
		}finally{
			if(success){deleteTree(root);}
			else{System.err.println("Retained ImgRecord child logs in "+root);}
		}
		System.out.println("IMGRECORD_TOARRAY_CLEANUP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private ImgRecordToArrayCleanupTest(){}

	/** Runs one bounded child JVM and validates its per-mode oracle.
	 * @param root Parent-owned temporary directory for child logs
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
			final ProcessBuilder pb=javaProcess(ImgRecordToArrayCleanupTest.class.getName(),
					assertions, "child", caseName, input.toString(),
					Boolean.toString(baseline));
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
			check(output.contains("case="),
					"Child "+label+" did not report an oracle marker: "+output);
		}finally{
			if(child!=null && child.isAlive()){terminateChild(child, label, log);}
		}
	}

	/** Runs one toArray case inside a fresh child JVM.
	 * @param caseName Case selector
	 * @param input Input file path to parse
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception if parsing or descriptor inspection diverges from the oracle */
	private static void runChildCase(final String caseName, final Path input,
			final boolean baseline) throws Exception{
		final boolean oldStoreName=ImgRecord.storeName;
		try{
			ImgRecord.storeName=true;
			if("normal".equals(caseName)){
				checkNormalRecords(input);
			}else if("hq".equals(caseName)){
				checkHighQualityRecords(input);
			}else if("malformed".equals(caseName)){
				checkMalformedFailure(input, baseline);
			}else if("close".equals(caseName)){
				checkCloseErrorStatus(input, baseline);
			}else if("combined".equals(caseName)){
				checkCombinedMalformedClose(input, baseline);
			}else{
				throw new IllegalArgumentException("Unknown child case: "+caseName);
			}
		}finally{
			ImgRecord.storeName=oldStoreName;
		}
	}

	/** Requires all three valid digit-starting IMG rows and no retained descriptor.
	 * @param input Input table with mixed HQ values
	 * @throws IOException if descriptor inspection fails */
	private static void checkNormalRecords(final Path input) throws IOException{
		final int before=countOpenDescriptors(input);
		final ImgRecord[] records=ImgRecord.toArray(input.toString(), false);
		check(records.length==3, "Normal input should emit three records");
		checkRecord(records[0], 11, "Genome Alpha", 101, true, false, 0, true);
		checkRecord(records[1], 22, "Genome Beta", 202, true, false, 2, false);
		checkRecord(records[2], 33, "Genome Gamma", -1, false, true, 1, true);
		checkNoLeak("normal", input, before);
	}

	/** Requires the high-quality filter to retain only records whose column 7 is yes.
	 * @param input Input table with mixed HQ values
	 * @throws IOException if descriptor inspection fails */
	private static void checkHighQualityRecords(final Path input) throws IOException{
		final int before=countOpenDescriptors(input);
		final ImgRecord[] records=ImgRecord.toArray(input.toString(), true);
		check(records.length==2, "HQ input should emit two records");
		checkRecord(records[0], 11, "Genome Alpha", 101, true, false, 0, true);
		checkRecord(records[1], 33, "Genome Gamma", -1, false, true, 1, true);
		checkNoLeak("hq", input, before);
	}

	/** Requires malformed nonempty taxID parsing to retain its exception semantics.
	 * @param input Input file with a bad taxID in the second record
	 * @param baseline Whether the parser is expected to leak the input descriptor
	 * @throws Exception if descriptor inspection fails */
	private static void checkMalformedFailure(final Path input, final boolean baseline)
			throws Exception{
		final int before=countOpenDescriptors(input);
		try{
			ImgRecord.toArray(input.toString(), false);
			throw new AssertionError("Malformed taxID should throw RuntimeException");
		}catch(RuntimeException expected){
			check(RuntimeException.class.equals(expected.getClass()),
					"Malformed taxID should throw a bare RuntimeException: "+expected);
			check(expected.getMessage()==null,
					"Malformed taxID RuntimeException should have null message: "+expected);
			check(expected.getCause()==null,
					"Malformed taxID RuntimeException should have null cause: "+expected);
			checkStackContains(expected, "tax.ImgRecord", "<init>");
			checkStackContains(expected, "tax.ImgRecord", "toArray");
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

	/** Requires baseline and fixed handling of a reported TextFile close failure.
	 * @param input Nonexistent stdin shadow path selected by ReadWrite's stdin handling
	 * @param baseline Whether the parser is expected to ignore the close-error status
	 * @throws Exception if the controlled stdin path is shadowed or not fully closed */
	private static void checkCloseErrorStatus(final Path input, final boolean baseline)
			throws Exception{
		check(STDIN_CLOSE_ERROR_NAME.equals(input.toString()),
				"Close-error case must use the stdin.* shadow path: "+input);
		check(!Files.exists(input), "Close-error stdin.* shadow path should not exist: "+input);

		final InputStream oldIn=System.in;
		final CloseErrorInputStream controlledIn=new CloseErrorInputStream(
				CLOSE_ERROR_LINE.getBytes(StandardCharsets.UTF_8));
		System.setIn(controlledIn);
		try{
			if(baseline){
				final ImgRecord[] records=ImgRecord.toArray(STDIN_CLOSE_ERROR_NAME, false);
				check(records.length==1, "Baseline close-error input should emit one record");
				checkRecord(records[0], 66, "Genome Zeta", 606, true, false, 0, true);
			}else{
				try{
					ImgRecord.toArray(STDIN_CLOSE_ERROR_NAME, false);
					throw new AssertionError("Fixed close-error input should throw RuntimeException");
				}catch(RuntimeException expected){
					final String message="Error closing "+STDIN_CLOSE_ERROR_NAME;
					check(message.equals(expected.getMessage()),
							"Close-error message should name the input; expected "+
							message+" but found "+expected);
				}
			}
		}finally{
			System.setIn(oldIn);
		}

		check(controlledIn.eofReached(), "Controlled stdin should have been consumed to EOF");
		check(controlledIn.closeCount()>0,
				"Controlled stdin should have been closed; found "+controlledIn.closeCount());
		System.out.println("case=close eof=true close_count="+controlledIn.closeCount());
	}

	/** Requires a fixed parser to suppress a close failure behind a parse failure.
	 * @param input Nonexistent stdin shadow path selected by ReadWrite's stdin handling
	 * @param baseline Whether the parser is expected to ignore the close-error status
	 * @throws Exception if the primary or suppressed failure diverges from the oracle */
	private static void checkCombinedMalformedClose(final Path input, final boolean baseline)
			throws Exception{
		check(!baseline, "Combined malformed-close child should be fixed-only");
		check(STDIN_COMBINED_ERROR_NAME.equals(input.toString()),
				"Combined case must use the stdin.* shadow path: "+input);
		check(!Files.exists(input), "Combined stdin.* shadow path should not exist: "+input);

		final InputStream oldIn=System.in;
		final CloseErrorInputStream controlledIn=new CloseErrorInputStream(
				COMBINED_ERROR_LINE.getBytes(StandardCharsets.UTF_8));
		System.setIn(controlledIn);
		try{
			try{
				ImgRecord.toArray(STDIN_COMBINED_ERROR_NAME, false);
				throw new AssertionError(
						"Combined malformed-close input should throw RuntimeException");
			}catch(RuntimeException expected){
				check(RuntimeException.class.equals(expected.getClass()),
						"Malformed taxID should remain primary: "+expected);
				check(expected.getMessage()==null,
						"Primary RuntimeException should have null message: "+expected);
				check(expected.getCause()==null,
						"Primary RuntimeException should have null cause: "+expected);
				checkStackContains(expected, "tax.ImgRecord", "<init>");
				checkStackContains(expected, "tax.ImgRecord", "toArray");

				final Throwable[] suppressed=expected.getSuppressed();
				check(suppressed.length==1,
						"Combined case should have one suppressed close failure: "+
						Arrays.toString(suppressed));
				final String message="Error closing "+STDIN_COMBINED_ERROR_NAME;
				check(RuntimeException.class.equals(suppressed[0].getClass()),
						"Suppressed close failure should be RuntimeException: "+suppressed[0]);
				check(message.equals(suppressed[0].getMessage()),
						"Suppressed close failure should name the input; expected "+
						message+" but found "+suppressed[0]);
			}
		}finally{
			System.setIn(oldIn);
		}

		check(controlledIn.closeCount()==1,
				"Combined stdin should be closed once; found "+controlledIn.closeCount());
		System.out.println("case=combined suppressed=1 close_count="+
				controlledIn.closeCount());
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
	 * @param record Parsed ImgRecord instance
	 * @param imgID Expected IMG identifier
	 * @param name Expected retained name
	 * @param taxID Expected NCBI taxonomy identifier
	 * @param isPublic Expected public flag
	 * @param obsolete Expected obsolete flag
	 * @param genomeType Expected genome-type index
	 * @param highQuality Expected high-quality flag */
	private static void checkRecord(final ImgRecord record, final long imgID,
			final String name, final int taxID, final boolean isPublic,
			final boolean obsolete, final int genomeType, final boolean highQuality){
		check(record.imgID==imgID,
				"Expected imgID "+imgID+" but found "+record.imgID);
		check(name.equals(record.name),
				"Expected name "+name+" but found "+record.name);
		check(record.taxID==taxID,
				"Expected taxID "+taxID+" but found "+record.taxID);
		check(record.isPublic==isPublic,
				"Expected public flag "+isPublic+" but found "+record.isPublic);
		check(record.obsolete==obsolete,
				"Expected obsolete flag "+obsolete+" but found "+record.obsolete);
		check(record.genomeType==genomeType,
				"Expected genomeType "+genomeType+" but found "+record.genomeType);
		check(record.highQuality==highQuality,
				"Expected highQuality "+highQuality+" but found "+record.highQuality);
	}

	/** Requires one frame to be present in the thrown exception's stack.
	 * @param thrown Exception to inspect
	 * @param className Expected stack-frame class
	 * @param methodName Expected stack-frame method */
	private static void checkStackContains(final Throwable thrown, final String className,
			final String methodName){
		for(StackTraceElement ste : thrown.getStackTrace()){
			if(className.equals(ste.getClassName()) && methodName.equals(ste.getMethodName())){
				return;
			}
		}
		throw new AssertionError("Expected stack to contain "+className+'.'+methodName+
				" but found "+Arrays.toString(thrown.getStackTrace()));
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
		check(Files.isDirectory(FD_ROOT), "ImgRecordToArrayCleanupTest requires /proc/self/fd");
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

		throw new AssertionError("Could not terminate "+label+" child: "+readAscii(log));
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
	 * @return One-line output for parent diagnostics */
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

	/** IMG-format table with a skipped header, mixed HQ values, and exact metadata. */
	private static final List<String> NORMAL_LINES=Arrays.asList(
			"img_id\tname\ttax_id\tpublic\tobsolete\tgenome_type\tunused\thq",
			"11\tGenome Alpha\t101\tyes\tno\tisolate\tunused\tyes",
			"22\tGenome Beta\t202\tyes\tno\tmetagenome\tunused\tno",
			"33\tGenome Gamma\t\tno\tyes\tsingle_cell\tunused\tyes");

	/** Valid first row followed by a digit-starting row with a malformed taxID. */
	private static final List<String> MALFORMED_LINES=Arrays.asList(
			"44\tGenome Delta\t404\tyes\tno\tisolate\tunused\tyes",
			"55\tGenome Epsilon\tbad\tyes\tno\tisolate\tunused\tyes");

	/** Nonexistent path that ReadWrite maps to System.in only when unshadowed. */
	private static final String STDIN_CLOSE_ERROR_NAME="stdin.imgrecord-close-error";

	/** Valid IMG payload supplied through controlled stdin. */
	private static final String CLOSE_ERROR_LINE=
			"66\tGenome Zeta\t606\tyes\tno\tisolate\tunused\tyes\n";

	/** Nonexistent path for combined malformed-record and close-failure checks. */
	private static final String STDIN_COMBINED_ERROR_NAME="stdin.imgrecord-combined-error";

	/** Malformed IMG payload supplied through controlled stdin. */
	private static final String COMBINED_ERROR_LINE=
			"77\tGenome Eta\tbad\tyes\tno\tisolate\tunused\tyes\n";

	/** Sentinel text emitted by the controlled stdin close failure. */
	private static final String CLOSE_ERROR_MESSAGE="TAX-070 controlled close failure";
}
