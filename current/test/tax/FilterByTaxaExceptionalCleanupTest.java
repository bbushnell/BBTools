package test.tax;

import java.io.File;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

import shared.Timer;
import tax.FilterByTaxa;
import tax.TaxFilter;

/** Bounded exceptional-cleanup baseline for FilterByTaxa after input and output startup.
 * @author Sucrose */
public final class FilterByTaxaExceptionalCleanupTest{

	/** Runs the parent harness or one leak-counting child.
	 * @param args Optional expect=baseline|fixed; child mode is internal
	 * @throws Exception If fixture construction, child launch or cleanup fails */
	public static void main(final String[] args) throws Exception{
		if(args.length>0 && "child".equals(args[0])){
			try{
				runChild(args);
				System.exit(0);
			}catch(Throwable failure){
				failure.printStackTrace();
				System.exit(1);
			}
			return;
		}
		runParent(parseExpectation(args));
	}

	/** Prevents construction of this standalone fixture. */
	private FilterByTaxaExceptionalCleanupTest(){}

	/** Runs cleanup and healthy-output children for t=1/t=4 and assertions on/off.
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @throws Exception If any child fails or times out */
	private static void runParent(final String expectation) throws Exception{
		final Path root=Files.createTempDirectory("filter-by-taxa-cleanup-");
		try{
			for(final int threads : THREAD_CASES){
				for(final boolean assertions : ASSERTION_CASES){
					checkChild(root, expectation, threads, assertions, "cleanup");
					checkChild(root, expectation, threads, assertions, "pass");
				}
			}
		}finally{
			deleteTree(root);
		}
		System.out.println("FILTER_BY_TAXA_EXCEPTIONAL_CLEANUP_TEST_PASS expect="+expectation);
	}

	/** Launches one bounded child JVM.
	 * @param root Parent-owned temporary tree
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @param threads Shared.threads value requested through t=
	 * @param assertions Whether to enable JVM assertions in the child
	 * @param mode cleanup for the throwing path, pass for the healthy control
	 * @throws Exception If the child fails or times out */
	private static void checkChild(final Path root, final String expectation, final int threads,
			final boolean assertions, final String mode) throws Exception{
		final String label=mode+"-t"+threads+"-"+(assertions ? "ea" : "da");
		final Path childRoot=root.resolve(label);
		final Path stdout=root.resolve(label+".stdout");
		final Path stderr=root.resolve(label+".stderr");
		Files.createDirectories(childRoot);

		Process child=null;
		try{
			child=new ProcessBuilder(childCommand(expectation, childRoot, threads, assertions, mode))
					.redirectOutput(stdout.toFile())
					.redirectError(stderr.toFile())
					.start();
			child.getOutputStream().close();
			if(!child.waitFor(CHILD_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
				terminateChild(child, label);
				throw new AssertionError(label+" timed out after "+CHILD_TIMEOUT_MS+" ms"+
						"\nstdout:\n"+textIfExists(stdout)+
						"\nstderr:\n"+textIfExists(stderr));
			}
			final String childOut=text(stdout), childErr=text(stderr);
			if(child.exitValue()!=0){
				throw new AssertionError(label+" exited "+child.exitValue()+
						"\nstdout:\n"+childOut+"\nstderr:\n"+childErr);
			}
			check(childOut.contains("case="+mode), label+" child did not report case: "+childOut);
			check(childOut.contains("threads="+threads), label+" child did not report threads: "+childOut);
			check(childOut.contains("assertions="+assertions),
					label+" child did not report assertion mode: "+childOut);
			System.out.print(childOut);
		}finally{
			if(child!=null && child.isAlive()){
				terminateChild(child, label);
			}
		}
	}

	/** Stops a child process without an unbounded forced wait.
	 * @param child Child process to stop
	 * @param label Child case label for diagnostics
	 * @throws Exception If the child cannot be stopped within the bounded waits */
	private static void terminateChild(final Process child, final String label) throws Exception{
		child.destroy();
		if(!child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
			child.destroyForcibly();
			if(!child.waitFor(DESTROY_TIMEOUT_MS, TimeUnit.MILLISECONDS)){
				throw new AssertionError(label+" ignored destroyForcibly for "+
						DESTROY_TIMEOUT_MS+" ms");
			}
		}
	}

	/** Builds a child JVM command.
	 * @param expectation Expected baseline or fixed cleanup behavior
	 * @param root Child artifact directory
	 * @param threads Shared.threads value to request through t=
	 * @param assertions Whether to enable JVM assertions
	 * @param mode cleanup or pass
	 * @return Child command */
	private static List<String> childCommand(final String expectation, final Path root,
			final int threads, final boolean assertions, final String mode){
		final List<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+File.separator+"bin"+File.separator+"java");
		command.add(assertions ? "-ea" : "-da");
		command.add("-Xmx512m");
		command.add("-Xms128m");
		command.add("-XX:ActiveProcessorCount="+threads);
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(FilterByTaxaExceptionalCleanupTest.class.getName());
		command.add("child");
		command.add("mode="+mode);
		command.add("expect="+expectation);
		command.add("root="+root);
		command.add("threads="+threads);
		return command;
	}

	/** Runs one child mode and prints its receipt.
	 * @param args Child command-line arguments
	 * @throws Exception If setup, processing or cleanup-state checks fail */
	private static void runChild(final String[] args) throws Exception{
		final ChildArgs parsed=parseChildArgs(args);
		Files.createDirectories(parsed.root);
		final Path in=parsed.root.resolve("tiny.fa");
		final Path out=parsed.root.resolve("owned.fa");
		writeTinyFasta(in);
		if("cleanup".equals(parsed.mode)){
			checkExceptionalCleanup(in, out, parsed);
		}else if("pass".equals(parsed.mode)){
			checkHealthyOutput(in, out, parsed);
		}else{
			throw new IllegalArgumentException("Unknown child mode: "+parsed.mode);
		}
	}

	/** Installs a throwing filter and verifies resources after the sentinel escapes.
	 * @param in Owned FASTA input
	 * @param out Owned FASTA output
	 * @param parsed Child settings
	 * @throws Exception If the sentinel or cleanup state diverges */
	private static void checkExceptionalCleanup(final Path in, final Path out,
			final ChildArgs parsed) throws Exception{
		final FilterByTaxa filter=newFilter(in, out, parsed.threads);
		check(reflectedFilter(filter).tree()==null, "No-tree constructor loaded a TaxTree");
		setFilter(filter, new ThrowingTaxFilter());

		try{
			filter.process(new Timer());
			throw new AssertionError("ThrowingTaxFilter should have thrown after startup");
		}catch(SentinelException expected){
			final int inputFds=countFds(in);
			final int outputFds=countFds(out);
			final int liveWriters=("baseline".equals(parsed.expectation) && parsed.threads==4) ?
					awaitLiveFastqWriterST2Workers() : countLiveFastqWriterST2Workers();
			if("baseline".equals(parsed.expectation)){
				check(outputFds==1, "Baseline should retain one owned output FD, found "+
						outputFds+" inputFds="+inputFds);
				check(liveWriters==(parsed.threads==4 ? 1 : 0),
						"Baseline t="+parsed.threads+" should have "+
						(parsed.threads==4 ? "one live" : "no")+
						" FastqWriterST2 worker, found "+liveWriters);
			}else{
				check(inputFds==0, "Fixed path should close owned input FDs, found "+inputFds);
				check(outputFds==0, "Fixed path should close owned output FDs, found "+outputFds);
				check(liveWriters==0,
						"Fixed path should stop FastqWriterST2 workers, found "+liveWriters);
			}
			System.out.println("case=cleanup expect="+parsed.expectation+
					" threads="+parsed.threads+" assertions="+assertionsEnabled()+
					" input_fds="+inputFds+" output_fds="+outputFds+
					" fastqst2_workers="+liveWriters+
					" failure="+expected.getClass().getSimpleName());
		}
	}

	/** Verifies that the same tree-free constructor can write a tiny normal FASTA.
	 * @param in Owned FASTA input
	 * @param out Owned FASTA output
	 * @param parsed Child settings
	 * @throws Exception If pass-through output or cleanup diverges */
	private static void checkHealthyOutput(final Path in, final Path out, final ChildArgs parsed)
			throws Exception{
		final FilterByTaxa filter=newFilter(in, out, parsed.threads);
		check(reflectedFilter(filter).tree()==null, "No-tree constructor loaded a TaxTree");
		filter.process(new Timer());
		final String output=text(out);
		check(TINY_FASTA.equals(output), "Healthy output mismatch: "+output);
		check(countFds(in)==0, "Healthy control should close owned input");
		check(countFds(out)==0, "Healthy control should close owned output");
		check(countLiveFastqWriterST2Workers()==0,
				"Healthy control should stop FastqWriterST2 workers");
		System.out.println("case=pass expect="+parsed.expectation+
				" threads="+parsed.threads+" assertions="+assertionsEnabled());
	}

	/** Builds a FilterByTaxa instance without real taxonomy resources.
	 * @param in Owned FASTA input
	 * @param out Owned FASTA output
	 * @param threads Shared thread count requested through Parser
	 * @return Configured fixture */
	private static FilterByTaxa newFilter(final Path in, final Path out, final int threads){
		return new FilterByTaxa(new String[] {
				"in="+in,
				"out="+out,
				"ow=t",
				"t="+threads});
	}

	/** Replaces the private final filter after construction but before process startup.
	 * @param target FilterByTaxa instance
	 * @param replacement Replacement filter
	 * @throws Exception If reflection fails */
	private static void setFilter(final FilterByTaxa target, final TaxFilter replacement)
			throws Exception{
		final Field field=FilterByTaxa.class.getDeclaredField("filter");
		field.setAccessible(true);
		field.set(target, replacement);
		check(reflectedFilter(target)==replacement, "Reflected filter replacement failed");
	}

	/** Reads the private filter field.
	 * @param target FilterByTaxa instance
	 * @return Current TaxFilter
	 * @throws Exception If reflection fails */
	private static TaxFilter reflectedFilter(final FilterByTaxa target) throws Exception{
		final Field field=FilterByTaxa.class.getDeclaredField("filter");
		field.setAccessible(true);
		return (TaxFilter)field.get(target);
	}

	/** Writes the finite two-record FASTA fixture.
	 * @param path Destination path
	 * @throws Exception If writing fails */
	private static void writeTinyFasta(final Path path) throws Exception{
		Files.write(path, TINY_FASTA.getBytes(StandardCharsets.US_ASCII));
	}

	/** Counts process descriptors pointing at one owned file.
	 * @param file File path to match
	 * @return Matching descriptor count
	 * @throws Exception If fdinfo cannot be read */
	private static int countFds(final Path file) throws Exception{
		if(!Files.exists(file)){return 0;}
		final Path real=file.toRealPath();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(final Path fd : fds){
				final Path target;
				try{
					target=Files.readSymbolicLink(fd);
				}catch(Exception e){
					continue;
				}
				if(target.equals(real)){count++;}
			}
		}
		return count;
	}

	/** Counts live worker threads in FastqWriterST2's private writer loop.
	 * @return Number of matching worker stacks */
	private static int countLiveFastqWriterST2Workers(){
		int count=0;
		for(final Map.Entry<Thread, StackTraceElement[]> entry :
				Thread.getAllStackTraces().entrySet()){
			if(entry.getKey()==Thread.currentThread() || !entry.getKey().isAlive()){
				continue;
			}
			for(final StackTraceElement frame : entry.getValue()){
				if("stream.FastqWriterST2$WriterRunnable".equals(frame.getClassName())){
					count++;
					break;
				}
			}
		}
		return count;
	}

	/** Waits briefly for a just-started baseline writer to reach WriterRunnable.
	 * @return Number of observed live workers
	 * @throws InterruptedException If interrupted while waiting */
	private static int awaitLiveFastqWriterST2Workers() throws InterruptedException{
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(2);
		int count=countLiveFastqWriterST2Workers();
		while(count<1 && System.nanoTime()<deadline){
			Thread.sleep(5);
			count=countLiveFastqWriterST2Workers();
		}
		return count;
	}

	/** Parses the BBTools-style expectation flag.
	 * @param args Command-line arguments
	 * @return baseline or fixed */
	private static String parseExpectation(final String[] args){
		if(args.length<1){return "fixed";}
		check(args.length==1, "unknown extra argument "+(args.length>1 ? args[1] : ""));
		final String arg=args[0];
		final String prefix="expect=";
		check(arg.regionMatches(true, 0, prefix, 0, prefix.length()),
				"unknown argument "+arg);
		final String expectation=arg.substring(prefix.length());
		check("baseline".equals(expectation) || "fixed".equals(expectation),
				"expectation must be baseline or fixed, not "+expectation);
		return expectation;
	}

	/** Parses child mode, expectation, root and thread-count settings.
	 * @param args Child arguments
	 * @return Parsed settings */
	private static ChildArgs parseChildArgs(final String[] args){
		check(args.length==5, "child requires mode, expect, root and threads");
		return new ChildArgs(parsePrefixed(args[1], "mode="),
				parsePrefixed(args[2], "expect="),
				Paths.get(parsePrefixed(args[3], "root=")),
				Integer.parseInt(parsePrefixed(args[4], "threads=")));
	}

	/** Parses one required prefixed argument.
	 * @param arg Complete argument
	 * @param prefix Required prefix
	 * @return Suffix */
	private static String parsePrefixed(final String arg, final String prefix){
		check(arg.startsWith(prefix), "Expected "+prefix+" argument, not "+arg);
		return arg.substring(prefix.length());
	}

	/** Adds the Vector API module on JDKs new enough to provide it.
	 * @param command Mutable command */
	private static void addVectorModule(final List<String> command){
		if(javaFeatureVersion()>=16){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
	}

	/** Returns the Java feature version.
	 * @return Feature version number */
	private static int javaFeatureVersion(){
		final String spec=System.getProperty("java.specification.version");
		return spec.startsWith("1.") ? Integer.parseInt(spec.substring(2)) :
			Integer.parseInt(spec);
	}

	/** Reads a small text file as US-ASCII.
	 * @param path File to read
	 * @return File contents
	 * @throws Exception If reading fails */
	private static String text(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Reads a small diagnostic file if it already exists.
	 * @param path File to read
	 * @return File contents, or a marker if the child did not create the file
	 * @throws Exception If an existing file cannot be read */
	private static String textIfExists(final Path path) throws Exception{
		return Files.exists(path) ? text(path) : "<missing>";
	}

	/** Recursively deletes a temporary tree.
	 * @param root Root directory
	 * @throws Exception If deletion fails */
	private static void deleteTree(final Path root) throws Exception{
		if(!Files.exists(root)){return;}
		final ArrayList<Path> paths=new ArrayList<Path>();
		try(final Stream<Path> stream=Files.walk(root)){
			stream.forEach(paths::add);
		}
		paths.sort((a, b)->b.compareTo(a));
		for(final Path path : paths){Files.deleteIfExists(path);}
	}

	/** Returns whether child JVM assertions are enabled.
	 * @return True if assertions are enabled */
	private static boolean assertionsEnabled(){
		boolean assertions=false;
		assert assertions=true;
		return assertions;
	}

	/** Throws an AssertionError independent of JVM assertion status.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** TaxFilter test seam that throws only after FilterByTaxa has started I/O. */
	private static final class ThrowingTaxFilter extends TaxFilter{
		/** Builds a no-tree filter used only as a post-construction test seam. */
		ThrowingTaxFilter(){super(null, false);}

		/** Always throws the TAX059 sentinel instead of applying production filter logic.
		 * @param name Read name that reached the post-startup filter point
		 * @return Never returns normally */
		@Override
		public boolean passesFilter(final String name){
			throw new SentinelException(name);
		}
	}

	/** Distinct process-body failure injected by ThrowingTaxFilter. */
	private static final class SentinelException extends RuntimeException{
		/** Builds a sentinel carrying the read name that triggered the failure.
		 * @param name Read name passed to TaxFilter */
		SentinelException(final String name){super("TAX059 sentinel for "+name);}
	}

	/** Child settings. */
	private static final class ChildArgs{
		/** Captures one child process mode.
		 * @param mode_ cleanup or pass
		 * @param expectation_ baseline or fixed
		 * @param root_ Child temporary directory
		 * @param threads_ Shared thread count requested through t= */
		ChildArgs(final String mode_, final String expectation_, final Path root_,
				final int threads_){
			mode=mode_;
			expectation=expectation_;
			root=root_;
			threads=threads_;
			check("baseline".equals(expectation) || "fixed".equals(expectation),
					"expectation must be baseline or fixed, not "+expectation);
			check(threads==1 || threads==4, "threads must be 1 or 4, not "+threads);
		}

		/** Child case to run: cleanup or pass. */
		final String mode;
		/** Expected FilterByTaxa cleanup behavior: baseline or fixed. */
		final String expectation;
		/** Directory for this child's owned inputs and outputs. */
		final Path root;
		/** Shared thread count to request through t=. */
		final int threads;
	}

	/** Shared.threads test values: one caller-driven writer and one worker-backed writer. */
	private static final int[] THREAD_CASES=new int[] {1, 4};
	/** JVM assertion modes that should not affect cleanup behavior. */
	private static final boolean[] ASSERTION_CASES=new boolean[] {true, false};
	/** Per-child startup and run timeout, in milliseconds. */
	private static final long CHILD_TIMEOUT_MS=15000;
	/** Bounded graceful or forced child-destruction wait, in milliseconds. */
	private static final long DESTROY_TIMEOUT_MS=500;
	/** Linux descriptor directory used to count owned open files. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Tiny two-record FASTA that the no-taxonomy FilterByTaxa can pass through. */
	private static final String TINY_FASTA=">taxa_keep_1\nACGT\n>taxa_keep_2\nTTAA\n";
}
