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
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.TimeUnit;

import tax.IDTree;

/** Bounded checks for IDTree input ownership and reported close failures.
 * @author Shinobu */
public final class IDTreeInputCleanupTest{

	/** Runs isolated cases or one child case.
	 * @param args Optional expect=baseline|fixed; child, case and root are internal
	 * @throws Exception If fixture setup, launch or validation fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false, child=false;
		String caseName=null;
		Path childRoot=null;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.equals("child=t")){child=true;}
			else if(arg.startsWith("case=")){caseName=arg.substring(5);}
			else if(arg.startsWith("root=")){childRoot=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(child){
			check(caseName!=null && childRoot!=null, "Child needs case and root");
			runChild(caseName, childRoot, baseline);
			return;
		}
		final Path root=Files.createTempDirectory("idtree-input-cleanup-");
		try{
			for(String name : new String[]{"normal", "limited", "malformed", "close-error"}){
				runParentCase(root, name, baseline, true);
				runParentCase(root, name, baseline, false);
			}
		}finally{
			deleteTree(root);
		}
		System.out.println("IDTREE_INPUT_CLEANUP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private IDTreeInputCleanupTest(){}

	/** Runs one bounded child and retains its diagnostic receipt in parent output.
	 * @param root Parent-owned temporary directory
	 * @param name Case name
	 * @param baseline Whether to require the original failure behavior
	 * @param assertions Whether the child uses assertions
	 * @throws Exception If the child cannot start, finish or pass */
	private static void runParentCase(final Path root, final String name,
			final boolean baseline, final boolean assertions) throws Exception{
		final String label=name+(assertions ? "-ea" : "-da");
		final Path childRoot=Files.createDirectory(root.resolve(label));
		final Path log=root.resolve(label+".log");
		final List<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+File.separator+"bin"+File.separator+"java");
		command.add(assertions ? "-ea" : "-da");
		command.add("-Xmx512m");
		command.add("-Xms512m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(IDTreeInputCleanupTest.class.getName());
		command.add("child=t");
		command.add("case="+name);
		command.add("root="+childRoot);
		command.add(baseline ? "expect=baseline" : "expect=fixed");
		Process process=null;
		final long started=System.nanoTime();
		try{
			process=new ProcessBuilder(command).redirectErrorStream(true)
					.redirectOutput(log.toFile()).start();
			process.getOutputStream().close();
			if(!process.waitFor(15, TimeUnit.SECONDS)){
				stop(process);
				throw new AssertionError("Timed out: "+label+" "+read(log));
			}
			final String output=read(log);
			final long millis=TimeUnit.NANOSECONDS.toMillis(System.nanoTime()-started);
			System.out.println(label+" exit="+process.exitValue()+" elapsed_ms="+millis+
					" output='"+output.replace('\n', ' ').replace('\r', ' ')+"'");
			check(process.exitValue()==0, label+" failed: "+output);
			check(output.contains("IDTREE_CASE_PASS case="+name), "Missing receipt: "+output);
		}finally{
			if(process!=null && process.isAlive()){stop(process);}
		}
	}

	/** Exercises the real public entry point in a fresh child JVM.
	 * @param name Case name
	 * @param root Owned child directory
	 * @param baseline Whether to require the original failure behavior
	 * @throws Exception If input preparation or a contract check fails */
	private static void runChild(final String name, final Path root,
			final boolean baseline) throws Exception{
		if(name.equals("close-error")){
			checkCloseError(root, baseline);
			return;
		}
		check(name.equals("normal") || name.equals("limited") || name.equals("malformed"),
				"Unknown case: "+name);
		final Path input=root.resolve("matrix.tsv");
		final Path output=root.resolve("tree.nwk");
		Files.write(input, (name.equals("normal") ? "A\n" : "A\nB\tbad\n")
				.getBytes(StandardCharsets.US_ASCII));
		check(countDescriptors(input)==0, "Fixture input should not already be open");
		final String[] command=arguments(input.toString(), output, name.equals("limited"));
		if(name.equals("malformed")){
			try{
				IDTree.main(command);
				throw new AssertionError("Malformed identity should fail");
			}catch(NumberFormatException expected){
				check("For input string: \"bad\"".equals(expected.getMessage()),
						"Unexpected primary parse failure: "+expected);
			}
			check(!Files.exists(output), "Failed parsing must not write a tree");
		}else{
			IDTree.main(command);
			check("A;".equals(read(output)), "Single retained leaf must be A;");
		}
		final int descriptors=countDescriptors(input);
		final int expected=name.equals("malformed") && baseline ? 1 : 0;
		check(descriptors==expected, name+" reader descriptors "+descriptors+" != "+expected);
		System.out.println("IDTREE_CASE_PASS case="+name+" read_fds="+descriptors);
	}

	/** Forces a real TextFile.close boolean error after valid stdin data.
	 * The Reader wrappers close this substituted stdin; ReadWrite reports the
	 * injected IOException through its boolean result rather than throwing it.
	 * @param root Owned child directory
	 * @param baseline Whether IDTree should currently ignore the close error
	 * @throws Exception If the close path or failure outcome is wrong */
	private static void checkCloseError(final Path root, final boolean baseline)
			throws Exception{
		final Path output=root.resolve("tree.nwk");
		final String inputName="stdin.idtree-close-fixture";
		check(!Files.exists(Paths.get(inputName)), "A physical file must not shadow the stdin fixture");
		final InputStream original=System.in;
		final CloseFailureInput input=new CloseFailureInput();
		RuntimeException failure=null;
		try{
			System.setIn(input);
			try{IDTree.main(arguments(inputName, output, false));}
			catch(RuntimeException e){failure=e;}
		}finally{
			System.setIn(original);
		}
		check(input.closeCalls>0, "Reader wrapper must reach the injected close failure");
		if(baseline){
			check(failure==null, "Baseline should ignore reported close failure: "+failure);
			check("A;".equals(read(output)), "Baseline should still write its single leaf");
		}else{
			check(failure!=null && !(failure instanceof NumberFormatException),
					"Close failure must reject otherwise valid input");
			check(failure.getMessage()!=null && failure.getMessage().contains("stdin"),
					"Close diagnostic must identify input: "+failure);
			check(!Files.exists(output), "Reported input close error must precede tree output");
		}
		System.out.println("IDTREE_CASE_PASS case=close-error close_calls="+input.closeCalls+
				" rejected="+(failure!=null));
	}

	/** Builds real IDTree arguments with quiet diagnostics and deterministic input.
	 * @param input Input path or stdin
	 * @param output Output Newick path
	 * @param limited Whether to retain only the first row
	 * @return Command-line arguments */
	private static String[] arguments(final String input, final Path output, final boolean limited){
		return new String[]{"in="+input, "out="+output, "outstream=null", "ow=t", "t=1",
				"unpigz=f", "gunzip=f", "unbgzip=f", "reads="+(limited ? 1 : -1)};
	}

	/** Counts descriptors for exactly the owned regular input file.
	 * @param input Input path
	 * @return Open descriptor count
	 * @throws IOException If procfs cannot be inspected */
	private static int countDescriptors(final Path input) throws IOException{
		final Path real=input.toRealPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){
			for(Path descriptor : descriptors){
				try{if(Files.readSymbolicLink(descriptor).equals(real)){count++;}}
				catch(IOException ignored){/* A descriptor may disappear during the scan. */}
			}
		}
		return count;
	}

	/** Enables the optional vector module only on runtimes that provide it.
	 * @param command Child command to extend
	 * @throws Exception If a present module system cannot be queried */
	private static void addVectorModule(final List<String> command) throws Exception{
		final Class<?> finderClass;
		try{finderClass=Class.forName("java.lang.module.ModuleFinder");}
		catch(ClassNotFoundException java8){return;}
		final Object finder=finderClass.getMethod("ofSystem").invoke(null);
		final Optional<?> module=(Optional<?>)finderClass.getMethod("find", String.class)
				.invoke(finder, "jdk.incubator.vector");
		if(module.isPresent()){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
	}

	/** Stops only the child owned by this fixture, with bounded waits.
	 * @param process Owned child
	 * @throws Exception If interruption or failed termination prevents cleanup */
	private static void stop(final Process process) throws Exception{
		if(!process.isAlive()){return;}
		process.destroy();
		if(!process.waitFor(1, TimeUnit.SECONDS)){
			process.destroyForcibly();
			check(process.waitFor(3, TimeUnit.SECONDS), "Owned child did not terminate after forced stop");
		}
	}

	/** Reads a small owned fixture or log.
	 * @param path File to read
	 * @return ASCII contents
	 * @throws IOException If reading fails */
	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Deletes only the temporary fixture tree after its children have stopped.
	 * @param root Owned temporary root
	 * @throws IOException If listing or deletion fails */
	private static void deleteTree(final Path root) throws IOException{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active with assertions disabled.
	 * @param condition Required condition
	 * @param message Diagnostic on failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Finite valid input whose close reports an IOException. */
	private static final class CloseFailureInput extends ByteArrayInputStream{
		/** Supplies one valid leaf without any physical stdin dependency. */
		CloseFailureInput(){super("A\n".getBytes(StandardCharsets.US_ASCII));}
		/** Reports the injected backend error on every actual close attempt.
		 * @throws IOException Always, to exercise TextFile's reported close status */
		@Override
		public void close() throws IOException{
			closeCalls++;
			throw new IOException("Injected IDTree input close failure");
		}
		/** Number of actual backend close attempts. */
		int closeCalls=0;
	}
}
