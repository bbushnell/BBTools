package test.fileIO;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;

/** Linux plain-backend regression fixture for ByteStreamWriter pre-start close.
 * @author Brian Bushnell, Sucrose */
public final class ByteStreamWriterPrestartCloseTest{

	/** Runs healthy and pre-start close children behind a finite timeout.
	 * @param args Optional expect=baseline or expect=fixed; child mode is internal
	 * @throws Exception If fixture setup or a child process fails */
	public static void main(final String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			check(args.length==4, "child requires child, case=..., expect=..., root=...");
			runChild(parseExpectation(args[2]), parseCase(args[1]), parseRoot(args[3]));
			return;
		}

		check(args.length<=1, "parent requires zero args or expect=...");
		final String expectation=(args.length==0 ? "fixed" : parseExpectation(args[0]));
		final Path root=Files.createTempDirectory("bsw-prestart-close-");
		try{
			runBounded("healthy", expectation, root);
			runBounded("prestart", expectation, root);
			if("fixed".equals(expectation)){
				runBounded("ordered", expectation, root);
				runBounded("forceprint", expectation, root);
			}
			deleteTree(root);
			System.out.println("BSW_PRESTART_CLOSE_TEST_PASS expect="+expectation);
		}catch(final Throwable t){
			System.err.println("BSW_PRESTART_CLOSE_TEST_ROOT "+root);
			throw t;
		}
	}

	/** Prevents construction of this standalone fixture. */
	private ByteStreamWriterPrestartCloseTest(){}

	/** Runs one child JVM with a bounded wall time.
	 * @param label Child case label
	 * @param expectation baseline or fixed
	 * @param root Parent artifact directory
	 * @throws Exception If the child hangs, fails or omits its receipt */
	private static void runBounded(final String label, final String expectation,
			final Path root) throws Exception{
		final Path childRoot=root.resolve(label);
		final Path stdout=root.resolve(label+".stdout");
		final Path stderr=root.resolve(label+".stderr");
		Files.createDirectories(childRoot);

		Process child=null;
		Throwable primary=null;
		final long started=System.nanoTime();
		try{
			child=new ProcessBuilder(childCommand(label, expectation, childRoot))
					.redirectOutput(stdout.toFile())
					.redirectError(stderr.toFile())
					.start();
			child.getOutputStream().close();
			if(!child.waitFor(CHILD_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
				final AssertionError failure=new AssertionError(label+
						" child timed out\nstdout\n"+
						readAscii(stdout)+"\nstderr\n"+readAscii(stderr));
				stopChild(child, failure);
				throw failure;
			}

			final String childOut=readAscii(stdout);
			final String childErr=readAscii(stderr);
			System.out.println(label+" child exit="+child.exitValue()+
					" elapsed_ms="+
					TimeUnit.NANOSECONDS.toMillis(System.nanoTime()-started)+
					" stdout='"+compact(childOut)+"' stderr='"+compact(childErr)+"'");
			check(child.exitValue()==0, label+" child exit "+child.exitValue()+
					"\nstdout\n"+childOut+"\nstderr\n"+childErr);
			check(childOut.contains("BSW_"+label.toUpperCase()+"_RECEIPT"),
					label+" child did not report receipt: "+childOut);
		}catch(final Throwable t){
			primary=t;
			throw t;
		}finally{
			stopChild(child, primary);
		}
	}

	/** Builds a child JVM command with the same assertion mode as the parent.
	 * @param label Child case label
	 * @param expectation baseline or fixed
	 * @param root Child artifact directory
	 * @return Child command */
	private static List<String> childCommand(final String label,
			final String expectation, final Path root){
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
		command.add(ByteStreamWriterPrestartCloseTest.class.getName());
		command.add("child");
		command.add("case="+label);
		command.add("expect="+expectation);
		command.add("root="+root.toString());
		return command;
	}

	/** Dispatches a child case.
	 * @param expectation baseline or fixed
	 * @param label Child case label
	 * @param root Child artifact directory
	 * @throws Exception If the case fails */
	private static void runChild(final String expectation, final String label,
			final Path root) throws Exception{
		check(Files.isDirectory(FD_ROOT),
				"ByteStreamWriterPrestartCloseTest requires /proc/self/fd");
		Files.createDirectories(root);
		if("healthy".equals(label)){runHealthy(root);}
		else if("prestart".equals(label)){runPrestart(expectation, root);}
		else if("ordered".equals(label)){runOrdered(root);}
		else if("forceprint".equals(label)){runForcePrint(root);}
		else{throw new IllegalArgumentException("Unknown child case: "+label);}
	}

	/** Runs a healthy start/write/close control over a plain file.
	 * @param root Child artifact directory
	 * @throws Exception If writing or descriptor closure diverges */
	private static void runHealthy(final Path root) throws Exception{
		final Path out=root.resolve("healthy.txt");
		createClosedFile(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.toString(), true,
				false, false);
		check(countWriterFds(out)==1,
				"Constructor should open one plain writer FD before start");
		writer.start();
		writer.println("healthy");
		final boolean error=writer.poisonAndWait();
		check(!error, "Healthy writer should not report errorState");
		check(countWriterFds(out)==0,
				"Healthy poisonAndWait should close the plain writer FD");
		check("healthy\n".equals(readAscii(out)),
				"Healthy writer emitted unexpected output: "+readAscii(out));
		System.out.println("BSW_HEALTHY_RECEIPT writer_fds="+
				countWriterFds(out)+" error="+error);
	}

	/** Runs the pre-start close API check over a constructor-opened plain file.
	 * @param expectation baseline or fixed
	 * @param root Child artifact directory
	 * @throws Exception If the baseline or fixed contract diverges */
	private static void runPrestart(final String expectation, final Path root)
			throws Exception{
		final Path out=root.resolve("prestart.txt");
		createClosedFile(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.toString(), true,
				false, false);
		check(countWriterFds(out)==1,
				"Constructor should open one plain writer FD before closeBeforeStart");
		check(writer.getState()==Thread.State.NEW,
				"New writer should be NEW before closeBeforeStart, not "+
				writer.getState());
		final Method closeBeforeStart=findCloseBeforeStart();

		if("baseline".equals(expectation)){
			check(closeBeforeStart==null,
					"Baseline should not expose closeBeforeStart yet");
			writer.start();
			final boolean error=writer.poisonAndWait();
			check(!error, "Baseline cleanup writer should close cleanly");
			check(countWriterFds(out)==0,
					"Baseline cleanup should close the plain writer FD after start");
			System.out.println("BSW_PRESTART_RECEIPT expect=baseline api=absent "+
					"constructor_fds=1 writer_fds="+countWriterFds(out));
			return;
		}

		check(closeBeforeStart!=null, "Fixed writer must expose closeBeforeStart");
		invokeCloseBeforeStart(closeBeforeStart, writer);
		check(writer.getState()==Thread.State.NEW,
				"closeBeforeStart should leave writer NEW, not "+
				writer.getState());
		check(countWriterFds(out)==0,
				"closeBeforeStart should close the constructor-opened writer FD");
		invokeCloseBeforeStart(closeBeforeStart, writer);
		check(writer.getState()==Thread.State.NEW,
				"Second closeBeforeStart should leave writer NEW, not "+
				writer.getState());
		check(countWriterFds(out)==0,
				"Second closeBeforeStart should be idempotent");

		final String startMessage=expectPreclosedStartRejection(writer);

		final boolean poisonError=writer.poisonAndWait();
		writer.waitForFinish();
		check(!poisonError, "Preclosed poisonAndWait should not report errorState");
		check(countWriterFds(out)==0,
				"Preclosed poisonAndWait/waitForFinish should not reopen output");
		System.out.println("BSW_PRESTART_RECEIPT expect=fixed api=present "+
				"constructor_fds=1 writer_fds="+countWriterFds(out)+
				" idempotent=true start_rejected=true"+
				" start_message='"+compact(startMessage)+"'"+
				" poison_error="+poisonError+
				" poison_wait_returned=true wait_returned=true");
	}

	/** Runs pre-start close over the plain ordered-output path.
	 * @param root Child artifact directory
	 * @throws Exception If the ordered pre-start contract diverges */
	private static void runOrdered(final Path root) throws Exception{
		final Path out=root.resolve("ordered.txt");
		createClosedFile(out);
		final ByteStreamWriter writer=newPlainWriter(out, true);
		check(countWriterFds(out)==1,
				"Ordered constructor should open one plain writer FD");
		check(writer.getState()==Thread.State.NEW,
				"New ordered writer should be NEW before closeBeforeStart, not "+
				writer.getState());
		final Method closeBeforeStart=findCloseBeforeStart();

		check(closeBeforeStart!=null,
				"Fixed ordered writer must expose closeBeforeStart");
		invokeCloseBeforeStart(closeBeforeStart, writer);
		check(writer.getState()==Thread.State.NEW,
				"Ordered closeBeforeStart should leave writer NEW, not "+
				writer.getState());
		check(countWriterFds(out)==0,
				"Ordered closeBeforeStart should close the plain writer FD");
		invokeCloseBeforeStart(closeBeforeStart, writer);
		check(countWriterFds(out)==0,
				"Second ordered closeBeforeStart should be idempotent");
		check(writer.getState()==Thread.State.NEW,
				"Second ordered closeBeforeStart should leave writer NEW, not "+
				writer.getState());

		final String startMessage=expectPreclosedStartRejection(writer);
		final boolean poisonError=writer.poisonAndWait();
		writer.waitForFinish();
		check(!poisonError,
				"Preclosed ordered poisonAndWait should not report errorState");
		check(countWriterFds(out)==0,
				"Preclosed ordered shutdown should not reopen output");
		check(Files.size(out)==0,
				"Preclosed ordered writer should not emit bytes: "+
				Files.size(out));
		System.out.println("BSW_ORDERED_RECEIPT writer_fds="+
				countWriterFds(out)+" empty_bytes="+Files.size(out)+
				" idempotent=true start_message='"+compact(startMessage)+"'"+
				" poison_error="+poisonError);
	}

	/** Runs a pre-start direct-byte write followed by closeBeforeStart.
	 * @param root Child artifact directory
	 * @throws Exception If direct bytes are not preserved exactly */
	private static void runForcePrint(final Path root) throws Exception{
		final Path out=root.resolve("forceprint.txt");
		final byte[] expected="prestart-header\n".getBytes(
				StandardCharsets.US_ASCII);
		createClosedFile(out);
		final ByteStreamWriter writer=newPlainWriter(out, false);
		writer.forcePrint(expected);
		check(countWriterFds(out)==1,
				"forcePrint before closeBeforeStart should keep one writer FD");
		check(writer.getState()==Thread.State.NEW,
				"New forcePrint writer should be NEW before closeBeforeStart, not "+
				writer.getState());
		final Method closeBeforeStart=findCloseBeforeStart();

		check(closeBeforeStart!=null,
				"Fixed forcePrint writer must expose closeBeforeStart");
		invokeCloseBeforeStart(closeBeforeStart, writer);
		check(countWriterFds(out)==0,
				"forcePrint closeBeforeStart should close the plain writer FD");
		check(Arrays.equals(expected, Files.readAllBytes(out)),
				"closeBeforeStart should flush exact pre-start forcePrint bytes");
		invokeCloseBeforeStart(closeBeforeStart, writer);
		check(countWriterFds(out)==0,
				"Second forcePrint closeBeforeStart should be idempotent");
		check(writer.getState()==Thread.State.NEW,
				"Second forcePrint closeBeforeStart should leave writer NEW, not "+
				writer.getState());
		check(Arrays.equals(expected, Files.readAllBytes(out)),
				"Second closeBeforeStart should preserve pre-start bytes");

		final String startMessage=expectPreclosedStartRejection(writer);
		final boolean poisonError=writer.poisonAndWait();
		writer.waitForFinish();
		check(!poisonError,
				"Preclosed forcePrint poisonAndWait should not report errorState");
		check(countWriterFds(out)==0,
				"Preclosed forcePrint shutdown should not reopen output");
		check(Arrays.equals(expected, Files.readAllBytes(out)),
				"Preclosed forcePrint shutdown should not truncate bytes");
		System.out.println("BSW_FORCEPRINT_RECEIPT writer_fds="+
				countWriterFds(out)+" preserved_bytes="+expected.length+
				" idempotent=true start_message='"+compact(startMessage)+"'"+
				" poison_error="+poisonError);
	}

	/** Creates and closes the plain output file before writer construction.
	 * @param out Output file path
	 * @throws Exception If the file cannot be created or descriptor state diverges */
	private static void createClosedFile(final Path out) throws Exception{
		Files.write(out, new byte[0]);
		check(Files.isRegularFile(out), "Output should be a regular file: "+out);
		check(countWriterFds(out)==0,
				"Closed plain file should have zero writer FDs before construction");
	}

	/** Opens a normal plain ByteStreamWriter without starting it.
	 * @param out Output file path
	 * @param ordered Whether to activate ordered job submission
	 * @return Newly constructed writer */
	private static ByteStreamWriter newPlainWriter(final Path out,
			final boolean ordered){
		final FileFormat ff=FileFormat.testOutput(out.toString(),
				FileFormat.TEXT, 0, 0, false, true, false, ordered);
		return new ByteStreamWriter(ff);
	}

	/** Looks up the fixed pre-start close API without binding the baseline.
	 * @return closeBeforeStart method, or null on the old API */
	private static Method findCloseBeforeStart(){
		try{return ByteStreamWriter.class.getMethod("closeBeforeStart");}
		catch(final NoSuchMethodException e){return null;}
	}

	/** Invokes closeBeforeStart and unwraps its primary runtime failure.
	 * @param method closeBeforeStart reflection handle
	 * @param writer writer to close
	 * @throws Exception If invocation fails */
	private static void invokeCloseBeforeStart(final Method method,
			final ByteStreamWriter writer) throws Exception{
		try{
			method.invoke(writer);
		}catch(final InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			throw e;
		}
	}

	/** Verifies that start() rejects a writer closed before startup.
	 * @param writer Preclosed writer
	 * @return The nonempty rejection message */
	private static String expectPreclosedStartRejection(
			final ByteStreamWriter writer){
		boolean startRejected=false;
		String startMessage=null;
		try{
			writer.start();
		}catch(final IllegalStateException expected){
			startRejected=true;
			startMessage=expected.getMessage();
		}
		check(startRejected, "start() should reject a preclosed writer");
		check(startMessage!=null && startMessage.length()>0,
				"start() rejection should explain why the writer is preclosed");
		return startMessage;
	}

	/** Counts write-only process descriptors pointing at one file.
	 * @param file File path to match
	 * @return Matching writer descriptor count
	 * @throws Exception If fdinfo cannot be read */
	private static int countWriterFds(final Path file) throws Exception{
		final Path real=file.toRealPath();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : fds){
				final Path target;
				try{
					target=Files.readSymbolicLink(fd);
				}catch(final Exception e){
					continue;
				}
				if(target.equals(real) && (fdFlags(fd)&O_ACCMODE)==O_WRONLY){
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

	/** Parses one BBTools-style expectation flag.
	 * @param arg Command-line argument
	 * @return baseline or fixed */
	private static String parseExpectation(final String arg){
		final String prefix="expect=";
		check(arg.startsWith(prefix), "unknown expectation argument "+arg);
		final String expectation=arg.substring(prefix.length());
		check("baseline".equals(expectation) || "fixed".equals(expectation),
				"expectation must be baseline or fixed, not "+expectation);
		return expectation;
	}

	/** Parses the child case flag.
	 * @param arg Child case argument
	 * @return Child case label */
	private static String parseCase(final String arg){
		final String prefix="case=";
		check(arg.startsWith(prefix), "unknown child argument "+arg);
		return arg.substring(prefix.length());
	}

	/** Parses the child artifact directory.
	 * @param arg Child root argument
	 * @return Child artifact directory */
	private static Path parseRoot(final String arg){
		final String prefix="root=";
		check(arg.startsWith(prefix), "unknown child argument "+arg);
		return Paths.get(arg.substring(prefix.length()));
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

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root owned temporary directory root
	 * @throws IOException if the tree cannot be walked or deleted */
	private static void deleteTree(final Path root) throws IOException{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(
					java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Reads an ASCII file.
	 * @param path File to read
	 * @return File contents
	 * @throws IOException If the file cannot be read */
	private static String readAscii(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Makes child output compact enough for a parent receipt.
	 * @param text Text to compact
	 * @return One-line compacted text */
	private static String compact(final String text){
		return text.trim().replace('\n', '|').replace('\r', '|');
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

	/** Test checks stay active under -da.
	 * @param condition condition required for test success
	 * @param message diagnostic message for failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Reports whether Java assertions are active in this JVM.
	 * @return True if Java assertions are active */
	private static boolean assertionsEnabled(){
		boolean enabled=false;
		assert enabled=true;
		return enabled;
	}

	/** Child process wall-time bound. */
	private static final long CHILD_TIMEOUT_MS=15000;
	/** Time to wait for graceful child destruction before kill. */
	private static final long DESTROY_TIMEOUT_MS=1000;
	/** Time to wait for forcible child destruction before failing. */
	private static final long FORCE_DESTROY_TIMEOUT_MS=1000;

	/** /proc descriptor directory for this JVM only. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Access-mode mask in Linux fdinfo flags. */
	private static final int O_ACCMODE=3;
	/** Linux write-only descriptor access mode. */
	private static final int O_WRONLY=1;
}
