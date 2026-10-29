package test.tax;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Timer;
import tax.ShrinkAccession;

/** Bounded actual-entry cleanup baseline for ShrinkAccession.process failures. */
public final class ShrinkAccessionCleanupTest{

	/** Runs the parent harness or one leak-counting child.
	 * @param args Optional expect=baseline|fixed; child mode is internal
	 * @throws Exception If fixture construction, child launch or cleanup fails */
	public static void main(final String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			try{
				runChild(parseCase(args), parseExpectation(args, "expect=", "fixed", true),
						parseExpectation(args, "startup=", null, true), parseRoot(args));
			}catch(Throwable failure){
				failure.printStackTrace();
				System.exit(1);
			}
			return;
		}
		runParent(parseExpectation(args, "expect=", "fixed", false),
				parseExpectation(args, "startup=", null, false));
	}

	/** Prevents construction of this standalone fixture. */
	private ShrinkAccessionCleanupTest(){}

	/** Runs one explicit ByteFile1 child in the parent JVM's assertion mode.
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @param startupExpectation Expected writer-start cleanup behavior, or null to skip
	 * @throws Exception If a child fails, hangs or cannot be cleaned up */
	private static void runParent(final String expectation,
			final String startupExpectation) throws Exception{
		boolean success=false;
		final Path root=Files.createTempDirectory("shrink-cleanup-");
		try{
			checkChild(root, "cleanup", expectation, null);
			if(startupExpectation!=null){
				checkChild(root, "writer-before", expectation, startupExpectation);
				checkChild(root, "writer-after", expectation, startupExpectation);
			}
			success=true;
		}finally{
			if(success){deleteTree(root);}
			else{System.err.println("Retained ShrinkAccession child logs in "+root);}
		}
		System.out.println("SHRINK_ACCESSION_CLEANUP_TEST_PASS expect="+expectation+
				" startup="+startupExpectation);
	}

	/** Launches one bounded ShrinkAccession child.
	 * @param root Parent-owned temporary tree
	 * @param caseName Child case selector
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @param startupExpectation Expected writer-start cleanup behavior
	 * @throws Exception If the child fails or times out */
	private static void checkChild(final Path root, final String caseName,
			final String expectation, final String startupExpectation) throws Exception{
		final String label=caseName+"-bf1-"+(assertionsEnabled() ? "ea" : "da");
		final Path childRoot=root.resolve(label);
		final Path stdout=root.resolve(label+".stdout");
		final Path stderr=root.resolve(label+".stderr");
		Files.createDirectories(childRoot);

		Process child=null;
		Throwable primary=null;
		final long started=System.nanoTime();
		try{
			child=new ProcessBuilder(childCommand(caseName, expectation,
					startupExpectation, childRoot))
					.redirectOutput(stdout.toFile())
					.redirectError(stderr.toFile())
					.start();
			child.getOutputStream().close();
			if(!waitForExit(child, CHILD_TIMEOUT_MS)){
				throw new AssertionError(label+" child timed out\nstdout\n"+
						readAscii(stdout)+"\nstderr\n"+readAscii(stderr));
			}
			final String childOut=readAscii(stdout);
			final String childErr=readAscii(stderr);
			System.out.println(label+" child exit="+child.exitValue()+
					" elapsed_ms="+
					TimeUnit.NANOSECONDS.toMillis(System.nanoTime()-started)+
					" stdout='"+compact(childOut)+"' stderr='"+compact(childErr)+"'");
			check(child.exitValue()==0, label+" child exit "+child.exitValue()+
					"\nstdout\n"+childOut+"\nstderr\n"+childErr);
			check(childOut.contains("case="+caseName),
					label+" child did not report retained counts: "+childOut);
		}catch(final Throwable failure){
			primary=failure;
			throw failure;
		}finally{
			stopChild(child, primary);
		}
	}

	/** Builds a child JVM command with the same assertion mode as the parent.
	 * @param caseName Child case selector
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @param startupExpectation Expected writer-start cleanup behavior, or null to skip
	 * @param root Child artifact directory
	 * @return Child command */
	private static List<String> childCommand(final String caseName,
			final String expectation, final String startupExpectation, final Path root){
		final List<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+File.separator+"bin"+
				File.separator+"java");
		command.add(assertionsEnabled() ? "-ea" : "-da");
		command.add("-Xmx512m");
		command.add("-Xms512m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(ShrinkAccessionCleanupTest.class.getName());
		command.add("child");
		command.add("case="+caseName);
		command.add("expect="+expectation);
		if(startupExpectation!=null){command.add("startup="+startupExpectation);}
		command.add("root="+root.toString());
		return command;
	}

	/** Runs one selected child case.
	 * @param caseName Child case selector
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @param startupExpectation Expected writer-start cleanup behavior
	 * @param root Child artifact directory
	 * @throws Throwable If the selected case diverges */
	private static void runChild(final String caseName, final String expectation,
			final String startupExpectation, final Path root) throws Throwable{
		if("cleanup".equals(caseName)){
			checkCleanupFailure(expectation, root);
		}else if("writer-before".equals(caseName)){
			checkWriterStartFailure(startupExpectation, StartMode.BEFORE_SUPER, root);
		}else if("writer-after".equals(caseName)){
			checkWriterStartFailure(startupExpectation, StartMode.AFTER_SUPER, root);
		}else{
			throw new IllegalArgumentException("Unknown child case: "+caseName);
		}
	}

	/** Runs the real failing ShrinkAccession invocation and inspects cleanup state.
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @param root Child artifact directory
	 * @throws Exception If the parse failure or cleanup state diverges */
	private static void checkCleanupFailure(final String expectation, final Path root) throws Exception{
		Files.createDirectories(root);
		final Path in=root.resolve("assembly.tsv");
		final Path out=root.resolve("assembly.out");
		Files.write(in, MALFORMED_ASSEMBLY.getBytes(StandardCharsets.US_ASCII));

		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=false;
		ByteFile.FORCE_MODE_BF3=false;
		ByteFile.FORCE_MODE_BF4=false;

		try{
			ShrinkAccession.main(new String[] {"in="+in.toString(),
					"out="+out.toString(), "asm=t", "ow=t"});
			throw new AssertionError("Truncated assembly row should fail");
		}catch(Throwable expected){
			checkExpectedFailure(expected);
			final int readers=countFds(in, O_RDONLY);
			final int writers=countFds(out, O_WRONLY);
			final int liveWriters="baseline".equals(expectation) ?
					awaitLiveByteStreamWritersInProcessJobs() : countLiveByteStreamWriters();
			if("baseline".equals(expectation)){
				check(readers==1, "Baseline should leak one reader FD, found "+readers);
				check(writers==1, "Baseline should leak one writer FD, found "+writers);
				check(liveWriters==1,
						"Baseline should retain one live ByteStreamWriter, found "+liveWriters);
			}else{
				check(readers==0, "Fixed path should close all reader FDs, found "+readers);
				check(writers==0, "Fixed path should close all writer FDs, found "+writers);
				check(liveWriters==0,
						"Fixed path should stop ByteStreamWriter, found "+liveWriters);
			}
			System.out.println("case=cleanup bf1=true expect="+expectation+
					" readers="+readers+" writers="+writers+
					" live_bsw="+liveWriters+" failure="+
					expected.getClass().getSimpleName()+" assertions="+assertionsEnabled());
			System.exit(0);
		}
	}

	/** Injects one writer.start failure and inspects directly opened resources.
	 * @param expectation Expected baseline or fixed writer-start behavior
	 * @param mode Whether to fail before or after {@code super.start()}
	 * @param root Child artifact directory
	 * @throws Throwable If startup or cleanup state diverges */
	private static void checkWriterStartFailure(final String expectation,
			final StartMode mode, final Path root) throws Throwable{
		check(expectation!=null, "writer startup cases require startup=baseline|fixed");
		Files.createDirectories(root);
		final Path in=root.resolve(mode.caseName+".tsv");
		final Path out=root.resolve(mode.caseName+".out");
		Files.write(in, MALFORMED_ASSEMBLY.getBytes(StandardCharsets.US_ASCII));

		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=false;
		ByteFile.FORCE_MODE_BF3=false;
		ByteFile.FORCE_MODE_BF4=false;

		lastFailingWriter=null;
		try{
			runInjectedShrink(in, out, mode);
			throw new AssertionError("Injected writer.start should fail");
		}catch(Throwable expected){
			check(expected==mode.sentinel,
					"ShrinkAccession must preserve writer.start sentinel identity");
			final int readers=countFds(in, O_RDONLY);
			final int writers=countFds(out, O_WRONLY);
			check(lastFailingWriter!=null, "Injected writer was not retained");
			final int liveWriters=mode==StartMode.AFTER_SUPER &&
					"baseline".equals(expectation) ?
					awaitLiveByteStreamWritersInProcessJobs() :
					countLiveByteStreamWriters();
			check(readers==0, "Writer-start failure should close the reader, found "+readers);
			if("baseline".equals(expectation)){
				check(writers==1, "Baseline should leak one writer FD, found "+writers);
				check(liveWriters==(mode==StartMode.AFTER_SUPER ? 1 : 0),
						"Baseline live writer count was "+liveWriters+" for "+mode);
			}else{
				check(writers==0, "Fixed path should close all writer FDs, found "+writers);
				check(liveWriters==0,
						"Fixed path should stop ByteStreamWriter, found "+liveWriters);
			}
			System.out.println("case="+mode.caseName+" expect="+expectation+
					" readers="+readers+" writers="+writers+
					" live_bsw="+liveWriters+" primary=exact");
			System.exit(0);
		}
	}

	/** Runs ShrinkAccession with an injected writer-start failure.
	 * @param in Input assembly-summary path
	 * @param out Output table path
	 * @param mode Writer start failure mode
	 * @throws Throwable If the invocation fails */
	private static void runInjectedShrink(final Path in, final Path out,
			final StartMode mode) throws Throwable{
		final InjectedShrinkAccession shrink=new InjectedShrinkAccession(new String[] {
				"in="+in.toString(), "out="+out.toString(), "asm=t", "ow=t"}, mode);
		final Method process=ShrinkAccession.class.getDeclaredMethod("process", Timer.class);
		process.setAccessible(true);
		try{
			process.invoke(shrink, new Timer());
		}catch(InvocationTargetException e){
			throw e.getCause();
		}
	}

	/** Requires the current assertion-mode-specific LineParser1 failure.
	 * @param thrown Throwable from the malformed assembly row */
	private static void checkExpectedFailure(final Throwable thrown){
		if(assertionsEnabled()){
			check(thrown instanceof AssertionError,
					"-ea should throw AssertionError, not "+thrown);
		}else{
			check(thrown instanceof ArrayIndexOutOfBoundsException,
					"-da should throw ArrayIndexOutOfBoundsException, not "+thrown);
		}
		final StackTraceElement[] stack=thrown.getStackTrace();
		check(stack.length>0, "Expected a nonempty failure stack");
		check("parse.LineParser1".equals(stack[0].getClassName()) &&
				"parseByte".equals(stack[0].getMethodName()),
				"Expected parse.LineParser1.parseByte origin, not "+
				Arrays.toString(stack));
	}

	/** Counts process descriptors pointing at one file with a requested access mode.
	 * @param file File path to match
	 * @param mode O_ACCMODE bits to match
	 * @return Matching descriptor count
	 * @throws Exception If fdinfo cannot be read */
	private static int countFds(final Path file, final int mode) throws Exception{
		final Path real=file.toRealPath();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : fds){
				final Path target;
				try{
					target=Files.readSymbolicLink(fd);
				}catch(Exception e){
					continue;
				}
				if(target.equals(real) && (fdFlags(fd)&O_ACCMODE)==mode){
					count++;
				}
			}
		}
		return count;
	}

	/** Reads one /proc descriptor flags field.
	 * @param fd Descriptor path under /proc/self/fd
	 * @return Parsed octal flags
	 * @throws Exception If fdinfo cannot be read */
	private static int fdFlags(final Path fd) throws Exception{
		final Path info=Paths.get("/proc/self/fdinfo",
				fd.getFileName().toString());
		for(String line : Files.readAllLines(info, StandardCharsets.US_ASCII)){
			if(line.startsWith("flags:")){
				return Integer.parseInt(line.substring(6).trim(), 8);
			}
		}
		throw new IllegalStateException("No flags line in "+info);
	}

	/** Counts all live ByteStreamWriter thread instances.
	 * @return Number of matching live workers */
	private static int countLiveByteStreamWriters(){
		int count=0;
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter && thread.isAlive()){
				count++;
			}
		}
		return count;
	}

	/** Counts live ByteStreamWriter workers blocked or active in processJobs.
	 * @return Number of matching live workers */
	private static int countLiveByteStreamWritersInProcessJobs(){
		int count=0;
		for(Map.Entry<Thread, StackTraceElement[]> entry :
				Thread.getAllStackTraces().entrySet()){
			if(entry.getKey()==Thread.currentThread() || !entry.getKey().isAlive()){
				continue;
			}
			for(StackTraceElement frame : entry.getValue()){
				if("fileIO.ByteStreamWriter".equals(frame.getClassName()) &&
						"processJobs".equals(frame.getMethodName())){
					count++;
					break;
				}
			}
		}
		return count;
	}

	/** Waits briefly for a just-started baseline writer to block in processJobs.
	 * @return Number of observed live ByteStreamWriter workers
	 * @throws InterruptedException If interrupted while waiting */
	private static int awaitLiveByteStreamWritersInProcessJobs() throws InterruptedException{
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(2);
		int count=countLiveByteStreamWritersInProcessJobs();
		while(count<1 && System.nanoTime()<deadline){
			Thread.sleep(5);
			count=countLiveByteStreamWritersInProcessJobs();
		}
		return count;
	}

	/** Parses one optional baseline-or-fixed flag.
	 * @param args Command-line arguments
	 * @param prefix Flag prefix to find
	 * @param defaultValue Value returned when no flag is present
	 * @param childMode Whether child-only flags are accepted
	 * @return baseline, fixed, or {@code defaultValue} */
	private static String parseExpectation(final String[] args, final String prefix,
			final String defaultValue, final boolean childMode){
		String expectation=defaultValue;
		boolean seen=false;
		for(String arg : args){
			if(startsWithName(arg, prefix)){
				check(!seen, "duplicate argument "+prefix);
				seen=true;
				expectation=arg.substring(prefix.length());
			}else{
				checkAllowedArgument(arg, childMode);
			}
		}
		if(expectation!=null){
			check("baseline".equals(expectation) || "fixed".equals(expectation),
					prefix+" must be baseline or fixed, not "+expectation);
		}
		return expectation;
	}

	/** Parses the child case selector.
	 * @param args Child command-line arguments
	 * @return Selected child case */
	private static String parseCase(final String[] args){
		final String value=parseValue(args, "case=");
		check(value!=null, "child requires case=...");
		return value;
	}

	/** Parses the child artifact directory.
	 * @param args Child command-line arguments
	 * @return Child artifact directory */
	private static Path parseRoot(final String[] args){
		final String value=parseValue(args, "root=");
		check(value!=null, "child requires root=...");
		return Paths.get(value);
	}

	/** Returns an optional child argument value by prefix.
	 * @param args Command-line arguments
	 * @param prefix Argument prefix
	 * @return Argument value, or null */
	private static String parseValue(final String[] args, final String prefix){
		String value=null;
		boolean seen=false;
		for(String arg : args){
			if(startsWithName(arg, prefix)){
				check(!seen, "duplicate argument "+prefix);
				seen=true;
				value=arg.substring(prefix.length());
			}else{
				checkAllowedArgument(arg, true);
			}
		}
		return value;
	}

	/** Requires one argument to be allowed in this JVM mode.
	 * @param arg Argument to check
	 * @param childMode Whether child-only flags are accepted */
	private static void checkAllowedArgument(final String arg, final boolean childMode){
		if("child".equalsIgnoreCase(arg)){check(childMode, "child flag is child-only");}
		else if(startsWithName(arg, "case=") || startsWithName(arg, "root=")){
			check(childMode, "Child-only argument in parent mode: "+arg);
		}else if(!startsWithName(arg, "expect=") &&
				!startsWithName(arg, "startup=")){
			throw new IllegalArgumentException("Unknown argument "+arg);
		}
	}

	/** Tests an argument name case-insensitively.
	 * @param arg Full argument
	 * @param prefix Expected name plus equals sign
	 * @return true if the argument starts with the expected name */
	private static boolean startsWithName(final String arg, final String prefix){
		return arg.regionMatches(true, 0, prefix, 0, prefix.length());
	}

	/** Builds a literal assembly_summary header line.
	 * @return Header line */
	private static String headerLine(){
		final String[] fields=filledLine();
		fields[0]="#assembly_accession";
		fields[5]="taxid";
		fields[23]="assembly_type";
		fields[25]="genome_size";
		fields[26]="genome_size_ungapped";
		fields[28]="replicon_count";
		fields[29]="scaffold_count";
		fields[30]="contig_count";
		fields[34]="total_gene_count";
		fields[35]="protein_coding_gene_count";
		fields[36]="non_coding_gene_count";
		return join(fields);
	}

	/** Builds the ignored fields in a literal 37-column assembly_summary row.
	 * @return Filled field array */
	private static String[] filledLine(){
		final String[] fields=new String[37];
		Arrays.fill(fields, "na");
		return fields;
	}

	/** Joins one assembly_summary row.
	 * @param fields Fields to join
	 * @return Tab-delimited row with trailing newline */
	private static String join(final String[] fields){
		final StringBuilder sb=new StringBuilder(256);
		for(int i=0; i<fields.length; i++){
			if(i>0){sb.append('\t');}
			sb.append(fields[i]);
		}
		return sb.append('\n').toString();
	}

	/** Checks whether Java assertions are enabled.
	 * @return True if assertions are enabled */
	private static boolean assertionsEnabled(){
		boolean enabled=false;
		assert enabled=true;
		return enabled;
	}

	/** Adds the incubator vector module on JDKs where it is present.
	 * @param command Command under construction */
	private static void addVectorModule(final List<String> command){
		try{
			Class.forName("jdk.incubator.vector.ByteVector");
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}catch(final ClassNotFoundException e){
			// Java 8: no module flag needed or accepted.
		}
	}

	/** Stops a child without hiding an earlier test failure.
	 * @param child Child process, possibly null
	 * @param primary Earlier failure to retain and annotate, if any */
	private static void stopChild(final Process child, final Throwable primary){
		if(child==null || !child.isAlive()){return;}
		boolean interrupted=false;
		Throwable failure=null;
		child.destroy();
		if(!waitForExit(child, DESTROY_TIMEOUT_MS)){
			child.destroyForcibly();
			if(!waitForExit(child, FORCE_DESTROY_TIMEOUT_MS)){
				failure=new AssertionError(
						"child stayed alive after destroyForcibly");
			}
		}
		if(Thread.interrupted()){
			interrupted=true;
			failure=addOrWrap(failure,
					new InterruptedException("interrupted while stopping child"));
		}
		if(interrupted){Thread.currentThread().interrupt();}
		if(failure!=null){
			if(primary!=null){primary.addSuppressed(failure);}
			else if(failure instanceof Error){throw (Error)failure;}
			else if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			else{throw new AssertionError(failure);}
		}
	}

	/** Waits for child exit while retaining any interrupt request.
	 * @param child Child process
	 * @param timeoutMs Wall-time bound in milliseconds
	 * @return True if the child exits before the bound */
	private static boolean waitForExit(final Process child, final long timeoutMs){
		final long deadline=System.nanoTime()+
				TimeUnit.MILLISECONDS.toNanos(timeoutMs);
		boolean interrupted=false;
		while(child.isAlive()){
			final long remaining=deadline-System.nanoTime();
			if(remaining<=0){break;}
			try{
				if(child.waitFor(remaining, TimeUnit.NANOSECONDS)){break;}
			}catch(final InterruptedException e){
				interrupted=true;
			}
		}
		if(interrupted){Thread.currentThread().interrupt();}
		return !child.isAlive();
	}

	/** Appends a suppressed failure to an existing failure.
	 * @param primary Existing failure, or null
	 * @param suppressed Failure to suppress
	 * @return The retained failure */
	private static Throwable addOrWrap(final Throwable primary,
			final Throwable suppressed){
		if(primary!=null){primary.addSuppressed(suppressed); return primary;}
		return suppressed;
	}

	/** Reads one US-ASCII file, returning empty text when absent.
	 * @param path File to read
	 * @return Decoded file contents
	 * @throws Exception If the file exists but cannot be read */
	private static String readAscii(final Path path) throws Exception{
		return Files.exists(path) ?
				new String(Files.readAllBytes(path), StandardCharsets.US_ASCII) : "";
	}

	/** Collapses process output for a parent-visible one-line child report.
	 * @param text Text to compact
	 * @return Single-line compact text */
	private static String compact(final String text){
		return text.trim().replaceAll("\\s+", " ");
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTree(final Path root) throws Exception{
		if(!Files.exists(root)){return;}
		try(final Stream<Path> paths=Files.walk(root)){
			for(final Path path : (Iterable<Path>)paths.sorted(
					Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under both assertion modes.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** ShrinkAccession subclass with one injectable writer factory. */
	private static final class InjectedShrinkAccession extends ShrinkAccession{
		/** Creates an injected ShrinkAccession.
		 * @param args BBTools arguments
		 * @param mode_ Writer start failure mode */
		InjectedShrinkAccession(final String[] args, final StartMode mode_){
			super(args);
			mode=mode_;
		}
		@Override
		protected ByteStreamWriter makeWriter(final FileFormat ff){
			lastFailingWriter=new StartFailingByteStreamWriter(ff, mode);
			return lastFailingWriter;
		}
		/** Writer start failure mode. */
		final StartMode mode;
	}

	/** Writer that fails before or after starting its worker. */
	private static final class StartFailingByteStreamWriter extends ByteStreamWriter{
		/** Creates an injected writer.
		 * @param ff Output descriptor
		 * @param mode_ Failure mode */
		StartFailingByteStreamWriter(final FileFormat ff, final StartMode mode_){
			super(ff);
			mode=mode_;
		}
		@Override
		public synchronized void start(){
			if(mode==StartMode.BEFORE_SUPER){throw mode.sentinel;}
			super.start();
			throw mode.sentinel;
		}
		/** Failure mode. */
		final StartMode mode;
	}

	/** Writer start failure sites. */
	private enum StartMode{
		/** Fail before starting the writer worker. */
		BEFORE_SUPER("writer-before", new RuntimeException("SHRINK_START_BEFORE")),
		/** Fail after starting the writer worker. */
		AFTER_SUPER("writer-after", new RuntimeException("SHRINK_START_AFTER"));
		/** Creates a failure mode.
		 * @param caseName_ Child case name
		 * @param sentinel_ Exact primary expected from writer.start() */
		StartMode(final String caseName_, final RuntimeException sentinel_){
			caseName=caseName_;
			sentinel=sentinel_;
		}
		/** Child case name. */
		final String caseName;
		/** Exact primary expected from writer.start(). */
		final RuntimeException sentinel;
	}

	/** Keeps an injected writer reachable until descriptor-count oracles finish. */
	private static StartFailingByteStreamWriter lastFailingWriter;

	/** Parent child-process wall-time bound. */
	private static final long CHILD_TIMEOUT_MS=15000;
	/** Time to wait for graceful child destruction before kill. */
	private static final long DESTROY_TIMEOUT_MS=1000;
	/** Time to wait for forcible child termination before reporting failure. */
	private static final long FORCE_DESTROY_TIMEOUT_MS=1000;

	/** /proc descriptor directory for this JVM only. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Access-mode mask in Linux fdinfo flags. */
	private static final int O_ACCMODE=3;
	/** Linux read-only descriptor access mode. */
	private static final int O_RDONLY=0;
	/** Linux write-only descriptor access mode. */
	private static final int O_WRONLY=1;

	/** Valid header followed by one non-comment row missing all metric fields. */
	private static final String MALFORMED_ASSEMBLY=headerLine()+"x\n";
}
