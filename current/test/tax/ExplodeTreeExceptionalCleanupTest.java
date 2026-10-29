package test.tax;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import structures.ByteBuilder;
import tax.ExplodeTree;
import tax.TaxTree;
import tax.TaxTreeText;

/** Baseline/fixed checks for ExplodeTree exceptional cleanup paths.
 * Baseline mode records skipped cleanup; fixed mode requires directory failures,
 * active writers and input readers to be closed or reported.
 * @author Sucrose */
public final class ExplodeTreeExceptionalCleanupTest{

	/** Runs parent or child checks for directory and processInner cleanup failures.
	 * @param args Optional {@code expect=baseline|fixed} or child-mode arguments
	 * @throws Exception if setup, child execution, or cleanup fails */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			final Boolean startupBaseline=parseStartupBaseline(args);
			runChildCase(args[1], Boolean.parseBoolean(args[2]), Paths.get(args[3]),
					startupBaseline);
			return;
		}

		boolean baseline=false;
		Boolean startupBaseline=null;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.equals("startup=baseline")){startupBaseline=Boolean.TRUE;}
			else if(arg.equals("startup=fixed")){startupBaseline=Boolean.FALSE;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		boolean success=false;
		final Path root=Files.createTempDirectory("explodetree-exceptional-parent-");
		try{
			runParentCase(root, "healthy", baseline, false, null);
			runParentCase(root, "healthy", baseline, true, null);
			runParentCase(root, "mkdirs", baseline, false, null);
			runParentCase(root, "mkdirs", baseline, true, null);
			runParentCase(root, "name", baseline, false, null);
			runParentCase(root, "name", baseline, true, null);
			runParentCase(root, "exceptional", baseline, false, null);
			runParentCase(root, "exceptional", baseline, true, null);
			if(startupBaseline!=null){
				runParentCase(root, "writer-before", baseline, false, startupBaseline);
				runParentCase(root, "writer-before", baseline, true, startupBaseline);
				runParentCase(root, "writer-after", baseline, false, startupBaseline);
				runParentCase(root, "writer-after", baseline, true, startupBaseline);
			}
			success=true;
		}finally{
			if(success){deleteTree(root);}
			else{System.err.println("Retained ExplodeTree child logs in "+root);}
		}
		System.out.println("EXPLODETREE_EXCEPTIONAL_CLEANUP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private ExplodeTreeExceptionalCleanupTest(){}

	/** Parses the optional child-side writer startup expectation.
	 * @param args Child-mode command-line arguments
	 * @return Expected startup baseline state, or null when absent
	 * @throws IllegalArgumentException if child-mode arguments are malformed */
	private static Boolean parseStartupBaseline(final String[] args){
		if(args.length==4){return null;}
		if(args.length!=5){
			throw new IllegalArgumentException(
					"Child mode requires 4 or 5 arguments, not "+args.length);
		}
		if("-".equals(args[4])){return null;}
		if("true".equals(args[4])){return Boolean.TRUE;}
		if("false".equals(args[4])){return Boolean.FALSE;}
		throw new IllegalArgumentException(
				"Unknown startup expectation: "+args[4]);
	}

	/** Runs one bounded child JVM and validates its receipt.
	 * @param root Parent-owned temporary directory for child logs
	 * @param caseName Child case selector
	 * @param baseline Whether the child should assert baseline or fixed behavior
	 * @param assertions Whether to enable assertions in the child JVM
	 * @param startupBaseline Expected writer-start behavior, or null to skip
	 * @throws Exception if the child cannot be launched, read, or verified */
	private static void runParentCase(final Path root, final String caseName,
			final boolean baseline, final boolean assertions,
			final Boolean startupBaseline) throws Exception{
		final String label=caseName+(assertions ? "-ea" : "-da");
		final Path log=root.resolve(label+".log");
		final Path childRoot=root.resolve(label+"-root");
		Process child=null;
		final long started=System.nanoTime();
		try{
			final ProcessBuilder pb=javaProcess(assertions,
					ExplodeTreeExceptionalCleanupTest.class.getName(), "child", caseName,
					Boolean.toString(baseline), childRoot.toString(),
					startupBaseline==null ? "-" : startupBaseline.toString());
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
			deleteTree(childRoot);
		}finally{
			if(child!=null && child.isAlive()){terminateChild(child, label, log);}
		}
	}

	/** Dispatches one child-side fixture case.
	 * @param caseName Case selector
	 * @param baseline Whether to assert baseline or fixed behavior
	 * @param root Parent-owned temporary directory assigned to this child
	 * @param startupBaseline Expected writer-start behavior, or null for non-startup cases
	 * @throws Exception if setup, tested code or cleanup fails */
	private static void runChildCase(final String caseName, final boolean baseline,
			final Path root, final Boolean startupBaseline) throws Exception{
		if("healthy".equals(caseName)){
			checkHealthyProcess(root);
		}else if("mkdirs".equals(caseName)){
			checkMkdirsFailure(root, baseline);
		}else if("name".equals(caseName)){
			checkNameWriteFailure(root, baseline);
		}else if("exceptional".equals(caseName)){
			checkExceptionalProcess(root, baseline);
		}else if("writer-before".equals(caseName)){
			check(startupBaseline!=null, "writer-before requires startup expectation");
			forceJvmGzip();
			checkWriterStartFailure(root, startupBaseline.booleanValue(),
					StartMode.BEFORE_SUPER);
		}else if("writer-after".equals(caseName)){
			check(startupBaseline!=null, "writer-after requires startup expectation");
			forceJvmGzip();
			checkWriterStartFailure(root, startupBaseline.booleanValue(),
					StartMode.AFTER_SUPER);
		}else{
			throw new IllegalArgumentException("Unknown child case: "+caseName);
		}
	}

	/** Requires the unchanged one-record success path to still close and write output.
	 * @param root Parent-owned child directory for this run
	 * @throws Exception if setup, output reading, or cleanup fails */
	private static void checkHealthyProcess(final Path root) throws Exception{
		Files.createDirectories(root);
		final Path input=Files.createTempFile(root, "input", ".fa");
		Files.write(input, new byte[0]);
		lastFailingWriter=null;
		final ThrowingByteFile reader=new ThrowingByteFile(input, false);
		final InjectedExplodeTree tree=new InjectedExplodeTree(input, root,
				TaxTreeContractTest.tree(TaxTree.GENUS_E), reader);
		Throwable failure=null;
		try{tree.process(new Timer());}
		catch(Throwable t){failure=t;}

		check(failure==null, "Healthy ExplodeTree process must not fail: "+failure);
		check(reader.closes==1, "Healthy ExplodeTree must close the reader once");
		check(EXPECTED_FASTA.equals(readOnlyOutput(root)),
				"Healthy ExplodeTree must preserve the one-record output");
		System.out.println("case=healthy close_count="+reader.closes);
	}

	/** Requires directory-creation failure to be visible only after the repair.
	 * @param root Parent-owned child directory for this run
	 * @param baseline Whether to assert the captured baseline
	 * @throws Exception if setup or cleanup fails */
	private static void checkMkdirsFailure(final Path root, final boolean baseline)
			throws Exception{
		Files.createDirectories(root);
		final Path rootFile=Files.createTempFile(root, "regular-root", ".tmp");
		final DirectoryExplodeTree tree=new DirectoryExplodeTree(singleNodeTree("Life"));
		tree.makeDirectoryTree(rootFile.toString(), false);
		check(!Files.exists(rootFile.resolve("1")),
				"Regular-file root must prevent creation of the TaxID 1 directory");
		if(baseline){
			check(!tree.errorState, "Baseline mkdirs failure must leave errorState false");
		}else{
			check(tree.errorState, "Fixed mkdirs failure must set errorState");
		}
		System.out.println("case=mkdirs errorState="+tree.errorState);
	}

	/** Requires .name write failure to be distinct from directory creation.
	 * @param root Parent-owned child directory for this run
	 * @param baseline Whether to assert the captured baseline
	 * @throws Exception if setup or cleanup fails */
	private static void checkNameWriteFailure(final Path root, final boolean baseline)
			throws Exception{
		Files.createDirectories(root);
		final PrintStream originalErr=System.err;
		final ByteArrayOutputStream printed=new ByteArrayOutputStream();
		try{
			System.setErr(new PrintStream(printed, true, "UTF-8"));
			final DirectoryExplodeTree tree=new DirectoryExplodeTree(singleNodeTree(longName()));
			tree.makeDirectoryTree(root.toString(), true);
			final Path nodeDir=root.resolve("1");
			check(Files.isDirectory(nodeDir),
					"Name-write case must create the TaxID 1 directory before failing");
			check(!Files.exists(nodeDir.resolve(longName()+".name")),
					"Overlong .name target must not be created");
			if(baseline){
				check(!tree.errorState, "Baseline .name failure must leave errorState false");
				check(printed.toString("UTF-8").contains("FileNotFoundException"),
						"Baseline .name failure must print and swallow the write failure");
			}else{
				check(tree.errorState, "Fixed .name failure must set errorState");
			}
			System.out.println("case=name errorState="+tree.errorState+
					" stderr_bytes="+printed.size());
		}finally{
			System.setErr(originalErr);
		}
	}

	/** Requires exceptional input to close the reader and stop the active writer when fixed.
	 * @param root Parent-owned child directory for this run
	 * @param baseline Whether to assert the captured baseline
	 * @throws Exception if setup or cleanup fails */
	private static void checkExceptionalProcess(final Path root, final boolean baseline)
			throws Exception{
		Files.createDirectories(root);
		final Path input=Files.createTempFile(root, "input", ".fa");
		Files.write(input, new byte[0]);
		final HashSet<Long> before=liveByteStreamWriters();
		final ThrowingByteFile reader=new ThrowingByteFile(input, true);
		final InjectedExplodeTree tree=new InjectedExplodeTree(input, root,
				TaxTreeContractTest.tree(TaxTree.GENUS_E), reader);

		Throwable failure=null;
		try{tree.process(new Timer());}
		catch(Throwable t){failure=t;}

		check(failure==ThrowingByteFile.SENTINEL,
				"ExplodeTree must preserve the primary nextLine sentinel identity: "+
				failure);
		if(baseline){
			final String thread=onlyNewLiveByteStreamWriter(before);
			check(reader.closes==0,
					"Baseline exceptional process must skip ByteFile.close");
			System.out.println("case=exceptional primary=sentinel close_count="+
					reader.closes+" before_writers="+before.size()+
					" leaked_writer="+thread);
			System.out.flush();
			System.exit(0);
		}else{
			check(reader.closes==1,
					"Fixed exceptional process must close the reader exactly once");
			check(noNewLiveByteStreamWriters(before),
					"Fixed exceptional process must poison and join the active writer");
			System.out.println("case=exceptional primary=sentinel close_count="+
					reader.closes+" before_writers="+before.size());
		}
	}

	/** Forces gzip output to use the in-JVM GZIPOutputStream branch. */
	private static void forceJvmGzip(){
		ReadWrite.USE_BGZF=false;
		ReadWrite.USE_BGZIP=false;
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_GZIP=false;
		ReadWrite.USE_UNPIGZ=false;
	}

	/** Requires writer-start failure cleanup to match the requested baseline.
	 * @param root Parent-owned child directory for this run
	 * @param baseline Whether to assert the captured baseline
	 * @param mode Writer start failure mode
	 * @throws Exception if setup or cleanup inspection fails */
	private static void checkWriterStartFailure(final Path root, final boolean baseline,
			final StartMode mode) throws Exception{
		Files.createDirectories(root);
		final Path input=Files.createTempFile(root, "input", ".fa");
		Files.write(input, new byte[0]);
		final HashSet<Long> before=liveByteStreamWriters();
		final ThrowingByteFile reader=new ThrowingByteFile(input, false);
		final InjectedExplodeTree tree=new InjectedExplodeTree(input, root,
				TaxTreeContractTest.tree(TaxTree.GENUS_E), reader, mode);

		Throwable failure=null;
		try{tree.process(new Timer());}
		catch(Throwable t){failure=t;}

		check(failure==mode.sentinel,
				"ExplodeTree must preserve writer.start sentinel identity: "+failure);
		final int writerFds=countWriteFdsUnder(root);
		check(lastFailingWriter!=null, "Injected writer was not retained");
		if(baseline){
			check(reader.closes==1,
					"Baseline writer-start failure must close the reader once");
			check(writerFds==1, "Baseline should leak one writer FD, found "+writerFds);
			if(mode==StartMode.AFTER_SUPER){
				final String thread=onlyNewLiveByteStreamWriter(before);
				System.out.println("case="+mode.caseName+" primary=exact close_count="+
						reader.closes+" writer_fds="+writerFds+
						" leaked_writer="+thread);
				System.out.flush();
				System.exit(0);
			}
			check(noNewLiveByteStreamWriters(before),
					"Before-super baseline must not start a writer thread");
		}else{
			check(reader.closes==1,
					"Fixed writer-start failure must close the reader once");
			check(writerFds==0, "Fixed path should close all writer FDs, found "+writerFds);
			check(noNewLiveByteStreamWriters(before),
					"Fixed path should leave no new live ByteStreamWriter");
		}
		System.out.println("case="+mode.caseName+" primary=exact close_count="+
				reader.closes+" writer_fds="+writerFds);
	}

	/** Creates a bounded child JVM command.
	 * @param assertions Whether to enable child assertions
	 * @param className Fully qualified class name to run
	 * @param args Arguments passed after the class name
	 * @return configured process builder */
	private static ProcessBuilder javaProcess(final boolean assertions,
			final String className, final String... args){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+"/bin/java");
		command.add(assertions ? "-ea" : "-da");
		command.add("-Xmx128m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(className);
		for(String arg : args){command.add(arg);}
		return new ProcessBuilder(command);
	}

	/** Adds the Vector incubator module when this JDK exposes it.
	 * @param command Child JVM command to extend */
	private static void addVectorModule(final ArrayList<String> command){
		try{
			final Class<?> moduleLayer=Class.forName("java.lang.ModuleLayer");
			final Object boot=moduleLayer.getMethod("boot").invoke(null);
			final Object optional=moduleLayer.getMethod("findModule", String.class)
					.invoke(boot, "jdk.incubator.vector");
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

	/** Reads the single emitted FASTA gzip payload from the directory tree.
	 * @param root Temporary output root
	 * @return Decompressed FASTA text emitted by ExplodeTree
	 * @throws Exception if the output cannot be found or read */
	private static String readOnlyOutput(final Path root) throws Exception{
		Path out=null;
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths::iterator){
				if(path.getFileName().toString().equals("4.fa.gz")){
					check(out==null, "Expected exactly one output file named 4.fa.gz");
					out=path;
				}
			}
		}
		check(out!=null, "ExplodeTree must write taxid 4 output under the output tree");
		try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(out));
				ByteArrayOutputStream buffer=new ByteArrayOutputStream()){
			byte[] temp=new byte[256];
			for(int len=in.read(temp); len>=0; len=in.read(temp)){
				buffer.write(temp, 0, len);
			}
			return new String(buffer.toByteArray(), StandardCharsets.UTF_8);
		}
	}

	/** Builds a one-node taxonomy for directory setup checks.
	 * @param name Literal root node name
	 * @return Loaded tree containing only TaxID 1
	 * @throws Exception if the tree cannot be serialized, loaded, or cleaned up */
	private static TaxTree singleNodeTree(final String name) throws Exception{
		final ByteBuilder bb=new ByteBuilder();
		bb.append("##TaxTree2\t1\n##Encoding\tA48\n##Tree\t2\t1\n");
		bb.append("##Config\t0\t0\t0\t0\t0\n##MergedMap\t0\t0\n");
		bb.append("##LevelExtended\tTID\tPID\tNumChildren\t");
		bb.append("MinParentLevelExtended\tMaxChildLevelExtended\tFlag\tName\n");
		bb.appendA48(TaxTree.LIFE_E).tab().appendA48(1).tab().appendA48(1)
				.tab().appendA48(0).tab().appendA48(TaxTree.LIFE_E).tab()
				.appendA48(0).tab().appendA48(TaxTree.LIFE_E).tab().append(name).nl();
		final Path file=Files.createTempFile("explodetree-one-node-", ".tsv");
		try{
			Files.write(file, bb.toBytes());
			return TaxTreeText.load(file.toString(), 1);
		}finally{
			Files.delete(file);
		}
	}

	/** Finds the live ByteStreamWriter thread ids visible to this JVM.
	 * @return Set of live ByteStreamWriter ids */
	private static HashSet<Long> liveByteStreamWriters(){
		final HashSet<Long> ids=new HashSet<Long>();
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter && thread.isAlive()){
				ids.add(Long.valueOf(thread.getId()));
			}
		}
		return ids;
	}

	/** Requires exactly one new ByteStreamWriter thread.
	 * @param before Live ByteStreamWriter ids before ExplodeTree started
	 * @return Diagnostic identity of the newly live writer */
	private static String onlyNewLiveByteStreamWriter(final HashSet<Long> before){
		Thread only=null;
		int count=0;
		for(Map.Entry<Thread, StackTraceElement[]> entry : Thread.getAllStackTraces().entrySet()){
			final Thread thread=entry.getKey();
			if(thread instanceof ByteStreamWriter && thread.isAlive() &&
					!before.contains(Long.valueOf(thread.getId()))){
				only=thread;
				count++;
			}
		}
		check(count==1, "Baseline exceptional process should leak one new ByteStreamWriter, not "+count);
		return only.getName()+"#"+only.getId()+":"+only.getState();
	}

	/** Checks that no extra ByteStreamWriter thread remains live.
	 * @param before Live ByteStreamWriter ids before ExplodeTree started
	 * @return true when no new ByteStreamWriter remains */
	private static boolean noNewLiveByteStreamWriters(final HashSet<Long> before){
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter && thread.isAlive() &&
					!before.contains(Long.valueOf(thread.getId()))){return false;}
		}
		return true;
	}

	/** Deletes an owned temporary tree after the tested writers have closed.
	 * @param root Temporary directory whose descendants are owned by this fixture
	 * @throws Exception if directory traversal or deletion fails */
	private static void deleteTree(final Path root) throws Exception{
		if(!Files.exists(root)){return;}
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Builds a filename that exceeds the usual filesystem component limit.
	 * @return Overlong alphanumeric filename stem */
	private static String longName(){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<300; i++){sb.append('A');}
		return sb.toString();
	}

	/** Reads an ASCII text file.
	 * @param path Path to read
	 * @return File text decoded as ASCII
	 * @throws Exception if reading fails */
	private static String readAscii(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Compacts child output into one parent log line.
	 * @param text Text to compact
	 * @return Single-line representation of the text */
	private static String compact(final String text){
		return text.replace('\n', '|').replace('\r', '|');
	}

	/** Keeps checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Counts write-only descriptors to files beneath one root.
	 * @param root Root directory
	 * @return Matching descriptor count
	 * @throws Exception if fdinfo cannot be read */
	private static int countWriteFdsUnder(final Path root) throws Exception{
		final Path real=root.toRealPath();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : fds){
				final Path target;
				try{
					target=Files.readSymbolicLink(fd);
				}catch(Exception e){
					continue;
				}
				if(target.startsWith(real) && (fdFlags(fd)&O_ACCMODE)==O_WRONLY){
					count++;
				}
			}
		}
		return count;
	}

	/** Reads one /proc descriptor flags field.
	 * @param fd Descriptor path under /proc/self/fd
	 * @return Parsed octal flags
	 * @throws Exception if fdinfo cannot be read */
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

	/** ExplodeTree variant that returns the injected finite reader instead of opening input. */
	private static final class InjectedExplodeTree extends ExplodeTree{
		/** Creates an ExplodeTree wired to a finite injected reader.
		 * @param input Placeholder input file path for FileFormat validation
		 * @param out Temporary output root
		 * @param tree_ Injected taxonomy tree
		 * @param reader_ Injected input reader */
		InjectedExplodeTree(final Path input, final Path out, final TaxTree tree_,
				final ByteFile reader_){
			this(input, out, tree_, reader_, null);
		}
		/** Creates an ExplodeTree wired to injected input and output.
		 * @param input Placeholder input file path for FileFormat validation
		 * @param out Temporary output root
		 * @param tree_ Injected taxonomy tree
		 * @param reader_ Injected input reader
		 * @param mode_ Optional writer start failure mode */
		InjectedExplodeTree(final Path input, final Path out, final TaxTree tree_,
				final ByteFile reader_, final StartMode mode_){
			super(FileFormat.testInput(input.toString(), FileFormat.FASTA, null, true, true),
					tree_, out.toString()+File.separator, "", true);
			reader=reader_;
			mode=mode_;
		}
		@Override
		protected ByteFile makeByteFile(final FileFormat ff){return reader;}
		@Override
		protected ByteStreamWriter makeWriter(final FileFormat ff){
			if(mode==null){return super.makeWriter(ff);}
			lastFailingWriter=new StartFailingByteStreamWriter(ff, mode);
			return lastFailingWriter;
		}
		/** Reader returned by {@link #makeByteFile(FileFormat)}. */
		final ByteFile reader;
		/** Optional writer start failure mode. */
		final StartMode mode;
	}

	/** ExplodeTree variant exposing directory construction with an injected tree. */
	private static final class DirectoryExplodeTree extends ExplodeTree{
		/** Creates a directory-only ExplodeTree with no input.
		 * @param tree_ Injected taxonomy tree */
		DirectoryExplodeTree(final TaxTree tree_){super(null, tree_, null, "", true);}
	}

	/** Minimal ByteFile that can throw after starting one output writer. */
	private static final class ThrowingByteFile extends ByteFile{
		/** Creates a finite or throwing reader over {@link ExplodeTreeExceptionalCleanupTest#LINES}.
		 * @param input Placeholder file path for the ByteFile base class
		 * @param throwAfterRecord_ Whether to throw after the first FASTA record */
		ThrowingByteFile(final Path input, final boolean throwAfterRecord_){
			super(FileFormat.testInput(input.toString(), FileFormat.FASTA, null, true, true));
			throwAfterRecord=throwAfterRecord_;
		}
		@Override
		public void reset(){next=0; open=true;}
		@Override
		public InputStream is(){return null;}
		@Override
		public long lineNum(){return next;}
		@Override
		public boolean close(){closes++; open=false; return false;}
		@Override
		public void pushBack(final byte[] line){throw new AssertionError("ExplodeTree does not use pushBack");}
		@Override
		public boolean isOpen(){return open;}
		@Override
		public byte[] nextLine(){
			if(!open){return null;}
			if(next<LINES.length){return LINES[next++].getBytes(StandardCharsets.UTF_8);}
			if(throwAfterRecord){throw SENTINEL;}
			return null;
		}
		/** Marker exception thrown after starting one writer. */
		static final RuntimeException SENTINEL=new RuntimeException("EXPLODETREE_SENTINEL");
		/** Whether to throw after the first FASTA record. */
		final boolean throwAfterRecord;
		/** Next input line index and number of close attempts. */
		int next=0, closes=0;
		/** Whether the reader can still return lines. */
		boolean open=true;
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
		BEFORE_SUPER("writer-before", new RuntimeException("EXPLODE_START_BEFORE")),
		/** Fail after starting the writer worker. */
		AFTER_SUPER("writer-after", new RuntimeException("EXPLODE_START_AFTER"));
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

	/** Injected FASTA lines read by ThrowingByteFile. */
	private static final String[] LINES={">tid|4 read0", "ACGT"};
	/** Expected decompressed output for the injected FASTA. */
	private static final String EXPECTED_FASTA=">tid|4 read0\nACGT\n";
	/** /proc descriptor directory for this JVM only. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Access-mode mask in Linux fdinfo flags. */
	private static final int O_ACCMODE=3;
	/** Linux write-only descriptor access mode. */
	private static final int O_WRONLY=1;
	/** Maximum time to wait for one child JVM. */
	private static final long CHILD_TIMEOUT_MS=30000;
	/** Time to wait for graceful child termination after a timeout or parent failure. */
	private static final long DESTROY_TIMEOUT_MS=5000;
}
