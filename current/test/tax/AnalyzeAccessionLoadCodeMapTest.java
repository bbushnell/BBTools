package test.tax;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.TimeUnit;

import tax.AnalyzeAccession;
import tax.AnalyzeAccession_ST;

/** Finite loadCodeMap ownership oracles for AnalyzeAccession variants.
 * @author Jean */
public final class AnalyzeAccessionLoadCodeMapTest{

	/** Runs regular-file and reported-close-error loader checks.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction, child launch or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		boolean child=false;
		String childMode=null;
		boolean st=false;
		Path input=null;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.startsWith("child=")){
				child=true;
				childMode=arg.substring(6);
			}else if(arg.equals("variant=mt")){st=false;}
			else if(arg.equals("variant=st")){st=true;}
			else if(arg.startsWith("input=")){input=Paths.get(arg.substring(6));}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		if(child){
			if(childMode.equals("regular")){checkRegularChild(st, input, baseline);}
			else if(childMode.equals("read-error")){checkReadErrorChild(st, baseline);}
			else if(childMode.equals("close-error")){checkCloseErrorChild(st, baseline);}
			else if(childMode.equals("read-close-error")){checkReadCloseErrorChild(st);}
			else{throw new IllegalArgumentException("Unknown child mode: "+childMode);}
			return;
		}

		final Path root=Files.createTempDirectory("analyze-accession-codemap-");
		try{
			for(boolean assertions : new boolean[] {true, false}){
				for(boolean variantSt : new boolean[] {false, true}){
					checkRegular(root, assertions, variantSt, baseline);
					checkReadError(root, assertions, variantSt, baseline);
					checkCloseError(root, assertions, variantSt, baseline);
					if(!baseline){checkReadCloseError(root, assertions, variantSt);}
				}
			}
		}finally{
			deleteTempTree(root);
		}
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private AnalyzeAccessionLoadCodeMapTest(){}

	/** Runs a regular-file ownership child.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param baseline True to expect the current open reader
	 * @throws Exception If the child launch or file IO fails */
	private static void checkRegular(Path root, boolean assertions, boolean st,
			boolean baseline) throws Exception{
		final String label=variant(st)+"-regular-"+(assertions ? "ea" : "da");
		final Path input=root.resolve(label+".txt");
		Files.write(input, MAP_TEXT.getBytes(StandardCharsets.US_ASCII));

		final RunResult result=runSelf(root, label, assertions, st, input, "regular",
				baseline);
		check(result.exit==0, label+" child exited "+result.exit+
				" with stderr\n"+result.stderr);
		check(result.stdout.contains("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS"),
				label+" child did not report pass\n"+result.stdout+
				"\nstderr\n"+result.stderr);
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_PARENT_CHILD_PASS case="+
				label+" "+result.stdout.trim());
	}

	/** Runs a read-error child.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param baseline True to expect the current skipped close
	 * @throws Exception If the child launch or file IO fails */
	private static void checkReadError(Path root, boolean assertions, boolean st,
			boolean baseline) throws Exception{
		final String label=variant(st)+"-read-error-"+(assertions ? "ea" : "da");
		final RunResult result=runSelf(root, label, assertions, st, null, "read-error",
				baseline);
		check(result.exit==0, label+" child exited "+result.exit+
				" with stderr\n"+result.stderr);
		check(result.stdout.contains("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS"),
				label+" child did not report pass\n"+result.stdout+
				"\nstderr\n"+result.stderr);
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_PARENT_CHILD_PASS case="+
				label+" "+result.stdout.trim());
	}

	/** Runs a reported-close-error child.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param baseline True to expect the current ignored close error
	 * @throws Exception If the child launch or file IO fails */
	private static void checkCloseError(Path root, boolean assertions, boolean st,
			boolean baseline) throws Exception{
		final String label=variant(st)+"-close-error-"+(assertions ? "ea" : "da");
		final RunResult result=runSelf(root, label, assertions, st, null, "close-error",
				baseline);
		check(result.exit==0, label+" child exited "+result.exit+
				" with stderr\n"+result.stderr);
		check(result.stdout.contains("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS"),
				label+" child did not report pass\n"+result.stdout+
				"\nstderr\n"+result.stderr);
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_PARENT_CHILD_PASS case="+
				label+" "+result.stdout.trim());
	}

	/** Loads a tiny regular pattern file in a child JVM and counts live FDs.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param input Pattern file path
	 * @param baseline True to expect the current open reader
	 * @throws Exception If loading, digitization or procfs inspection fails */
	private static void checkRegularChild(boolean st, Path input, boolean baseline)
			throws Exception{
		check(input!=null, "Missing input path");
		check(countDescriptors(input)==0, "Input descriptor open before loadCodeMap");
		checkMap(load(st, input.toString()), st);
		checkDigitize(st);

		final int descriptors=countDescriptors(input);
		final int expected=baseline ? 1 : 0;
		check(descriptors==expected,
				"Expected "+expected+" live pattern descriptor(s); found "+descriptors);
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS variant="+
				variant(st)+" case=regular read_fds="+descriptors+
				" expect="+(baseline ? "baseline" : "fixed"));
	}

	/** Forces a TextFile read error after finite valid stdin data.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param baseline True to expect the current skipped close
	 * @throws Exception If loader behavior differs from the selected expectation */
	private static void checkReadErrorChild(boolean st, boolean baseline)
			throws Exception{
		final String inputName="stdin.analyze-accession-codemap-read-"+variant(st);
		check(!Files.exists(Paths.get(inputName)),
				"A physical file must not shadow the stdin fixture");
		final InputStream original=System.in;
		final ReadFailureInput input=new ReadFailureInput();
		RuntimeException failure=null;
		try{
			System.setIn(input);
			try{load(st, inputName);}
			catch(RuntimeException e){failure=e;}
		}finally{
			System.setIn(original);
		}

		check(System.in==original, "System.in was not restored after read failure");
		check(failure!=null, "Injected read failure did not reject "+variant(st));
		check(rootCause(failure) instanceof IOException,
				"Original read error was not preserved: "+failure);
		check(rootCause(failure).getMessage()!=null &&
				rootCause(failure).getMessage().contains("read failure"),
				"Unexpected read-error cause: "+rootCause(failure));
		check(input.closeCalls==(baseline ? 0 : 1),
				"Unexpected read-error close count "+input.closeCalls);
		checkNoCodeMap(st);
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS variant="+
				variant(st)+" case=read-error close_calls="+input.closeCalls+
				" expect="+(baseline ? "baseline" : "fixed"));
	}

	/** Runs a fixed read-plus-close-error child.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @throws Exception If the child launch or file IO fails */
	private static void checkReadCloseError(Path root, boolean assertions, boolean st)
			throws Exception{
		final String label=variant(st)+"-read-close-error-"+
				(assertions ? "ea" : "da");
		final RunResult result=runSelf(root, label, assertions, st, null,
				"read-close-error", false);
		check(result.exit==0, label+" child exited "+result.exit+
				" with stderr\n"+result.stderr);
		check(result.stdout.contains("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS"),
				label+" child did not report pass\n"+result.stdout+
				"\nstderr\n"+result.stderr);
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_PARENT_CHILD_PASS case="+
				label+" "+result.stdout.trim());
	}

	/** Forces a TextFile close error after otherwise valid stdin data.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param baseline True to expect the current ignored close error
	 * @throws Exception If loader behavior differs from the selected expectation */
	private static void checkCloseErrorChild(boolean st, boolean baseline)
			throws Exception{
		final String inputName="stdin.analyze-accession-codemap-"+variant(st);
		check(!Files.exists(Paths.get(inputName)),
				"A physical file must not shadow the stdin fixture");
		final InputStream original=System.in;
		final CloseFailureInput input=new CloseFailureInput();
		RuntimeException failure=null;
		try{
			System.setIn(input);
			try{
				checkMap(load(st, inputName), st);
				checkDigitize(st);
			}catch(RuntimeException e){
				failure=e;
			}
		}finally{
			System.setIn(original);
		}

		check(System.in==original, "System.in was not restored after close failure");
		if(baseline){
			check(input.closeCalls==0,
					"Baseline should not reach close; calls="+input.closeCalls);
			check(failure==null, "Baseline should ignore close failure: "+failure);
		}else{
			check(input.closeCalls>0, "Fixed loader did not close injected input");
			check(failure!=null,
					"Fixed loader should reject reported close failure");
			check(failure.getMessage()!=null && failure.getMessage().contains("stdin"),
					"Close diagnostic must identify input: "+failure);
			checkNoCodeMap(st);
		}
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS variant="+
				variant(st)+" case=close-error close_calls="+input.closeCalls+
				" rejected="+(failure!=null)+" expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Forces read and close failures in the same fixed loader call.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @throws Exception If loader behavior differs from the fixed expectation */
	private static void checkReadCloseErrorChild(boolean st) throws Exception{
		final String inputName="stdin.analyze-accession-codemap-read-close-"+variant(st);
		check(!Files.exists(Paths.get(inputName)),
				"A physical file must not shadow the stdin fixture");
		final InputStream original=System.in;
		final ReadFailureInput input=new ReadFailureInput(true);
		RuntimeException failure=null;
		try{
			System.setIn(input);
			try{load(st, inputName);}
			catch(RuntimeException e){failure=e;}
		}finally{
			System.setIn(original);
		}

		check(System.in==original, "System.in was not restored after read failure");
		check(failure!=null, "Combined read and close failure did not reject "+
				variant(st));
		check(rootCause(failure) instanceof IOException,
				"Original read error was not preserved: "+failure);
		check(rootCause(failure).getMessage()!=null &&
				rootCause(failure).getMessage().contains("read failure"),
				"Unexpected read-error cause: "+rootCause(failure));
		check(input.closeCalls==1,
				"Fixed loader did not close injected input once: "+input.closeCalls);
		check(failure.getSuppressed().length==1,
				"Expected one suppressed close diagnostic; found "+
				failure.getSuppressed().length);
		check(failure.getSuppressed()[0].getMessage()!=null &&
				failure.getSuppressed()[0].getMessage().contains(inputName),
				"Suppressed close diagnostic did not identify input: "+
				failure.getSuppressed()[0]);
		checkNoCodeMap(st);
		System.out.println("ANALYZE_ACCESSION_LOAD_CODE_MAP_CASE_PASS variant="+
				variant(st)+" case=read-close-error close_calls="+
				input.closeCalls+" suppressed="+failure.getSuppressed().length+
				" expect=fixed");
	}

	/** Runs this fixture in one isolated child mode.
	 * @param root Temporary artifact directory
	 * @param label Case label for captured stdout/stderr
	 * @param assertions True for -ea, false for -da
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param input Regular input path, or null for stdin injection
	 * @param childMode Child mode to run
	 * @param baseline True to expect current baseline behavior
	 * @return Child exit status and captured output
	 * @throws Exception If process launch, wait or file IO fails */
	private static RunResult runSelf(Path root, String label, boolean assertions,
			boolean st, Path input, String childMode, boolean baseline)
			throws Exception{
		final ArrayList<String> command=baseCommand(assertions);
		command.add(AnalyzeAccessionLoadCodeMapTest.class.getName());
		command.add(baseline ? "expect=baseline" : "expect=fixed");
		command.add("child="+childMode);
		command.add(st ? "variant=st" : "variant=mt");
		if(input!=null){command.add("input="+input.toString());}
		return run(root, label, command);
	}

	/** Loads one variant's pattern table.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @param input Pattern input path
	 * @return Pattern-to-index map */
	private static HashMap<String, Integer> load(boolean st, String input){
		return st ? AnalyzeAccession_ST.loadCodeMap(input) :
				AnalyzeAccession.loadCodeMap(input);
	}

	/** Checks the tiny literal pattern map.
	 * @param map Loaded pattern map
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession */
	private static void checkMap(HashMap<String, Integer> map, boolean st){
		check(map.size()==1, variant(st)+" map size="+map.size());
		check(map.containsKey(PATTERN), variant(st)+" missing "+PATTERN);
		check(map.get(PATTERN).intValue()==0,
				variant(st)+" unexpected pattern index "+map.get(PATTERN));
	}

	/** Checks String and byte-array digitization against the literal map.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession */
	private static void checkDigitize(boolean st){
		final byte[] bytes=ACCESSION.getBytes(StandardCharsets.US_ASCII);
		final byte[] missing=MISSING.getBytes(StandardCharsets.US_ASCII);
		if(st){
			check(AnalyzeAccession_ST.digitize(ACCESSION)==EXPECTED_KEY,
					"ST String digitize failed");
			check(AnalyzeAccession_ST.digitize(bytes)==EXPECTED_KEY,
					"ST byte digitize failed");
			check(AnalyzeAccession_ST.digitize(MISSING)==-2,
					"ST missing String pattern should return -2");
			check(AnalyzeAccession_ST.digitize(missing)==-2,
					"ST missing byte pattern should return -2");
		}else{
			check(AnalyzeAccession.digitize(ACCESSION)==EXPECTED_KEY,
					"MT String digitize failed");
			check(AnalyzeAccession.digitize(bytes)==EXPECTED_KEY,
					"MT byte digitize failed");
			check(AnalyzeAccession.digitize(MISSING)==-2,
					"MT missing String pattern should return -2");
			check(AnalyzeAccession.digitize(missing)==-2,
					"MT missing byte pattern should return -2");
		}
	}

	/** Checks that no static codeMap has been published.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession */
	private static void checkNoCodeMap(boolean st){
		if(st){
			check(AnalyzeAccession_ST.codeMap==null, "ST codeMap was published");
		}else{
			check(AnalyzeAccession.codeMap==null, "MT codeMap was published");
		}
	}

	/** Returns the final cause for a wrapped failure.
	 * @param failure Throwable to inspect
	 * @return Deepest cause, or failure itself when it has no cause */
	private static Throwable rootCause(Throwable failure){
		Throwable root=failure;
		while(root.getCause()!=null){root=root.getCause();}
		return root;
	}

	/** Builds a Java command with the active classpath.
	 * @param assertions True for -ea, false for -da
	 * @return Mutable command list
	 * @throws Exception If the vector module query fails */
	private static ArrayList<String> baseCommand(boolean assertions)
			throws Exception{
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+"/bin/java");
		command.add("-Xmx128m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		command.add("-XX:ActiveProcessorCount=1");
		command.add(assertions ? "-ea" : "-da");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		return command;
	}

	/** Runs one child process with bounded cleanup.
	 * @param root Temporary artifact directory
	 * @param label Case label for captured stdout/stderr
	 * @param command Complete child command
	 * @return Child exit status and captured output
	 * @throws Exception If process launch, wait or file IO fails */
	private static RunResult run(Path root, String label, List<String> command)
			throws Exception{
		final Path stdout=root.resolve(label+".stdout");
		final Path stderr=root.resolve(label+".stderr");
		final Process child=new ProcessBuilder(command)
				.redirectOutput(stdout.toFile())
				.redirectError(stderr.toFile())
				.start();
		try{
			child.getOutputStream().close();
			check(child.waitFor(15, TimeUnit.SECONDS),
					label+" child did not finish within its bounded wait");
			return new RunResult(child.exitValue(), readAscii(stdout), readAscii(stderr));
		}finally{
			stop(child);
		}
	}

	/** Stops only the child owned by this fixture, with bounded waits.
	 * @param process Owned child
	 * @throws Exception If interruption or failed termination prevents cleanup */
	private static void stop(Process process) throws Exception{
		if(!process.isAlive()){return;}
		process.destroy();
		if(!process.waitFor(1, TimeUnit.SECONDS)){
			process.destroyForcibly();
			check(process.waitFor(1, TimeUnit.SECONDS),
					"Child process ignored forcible termination");
		}
	}

	/** Counts descriptors for exactly the owned regular pattern file.
	 * @param input Input path
	 * @return Open descriptor count
	 * @throws IOException If procfs cannot be inspected */
	private static int countDescriptors(Path input) throws IOException{
		final Path real=input.toRealPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(FD_ROOT)){
			for(Path descriptor : descriptors){
				try{
					if(Files.readSymbolicLink(descriptor).equals(real)){count++;}
				}catch(IOException ignored){
					/* A descriptor may disappear during the scan. */
				}
			}
		}
		return count;
	}

	/** Enables the optional vector module only on runtimes that provide it.
	 * @param command Child command to extend
	 * @throws Exception If a present module system cannot be queried */
	private static void addVectorModule(List<String> command) throws Exception{
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

	/** Returns the short variant name for messages.
	 * @param st True for AnalyzeAccession_ST, false for AnalyzeAccession
	 * @return st or mt */
	private static String variant(boolean st){return st ? "st" : "mt";}

	/** Reads one US-ASCII file, returning empty text when absent.
	 * @param path File to read
	 * @return Decoded file contents
	 * @throws IOException If the file exists but cannot be read */
	private static String readAscii(Path path) throws IOException{
		return Files.exists(path) ?
				new String(Files.readAllBytes(path), StandardCharsets.US_ASCII) : "";
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws IOException If walking or deletion fails */
	private static void deleteTempTree(Path root) throws IOException{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Tiny pattern table with one comment and one encodable pattern row. */
	private static final String MAP_TEXT="#comment\nLLLDDD\tignored\n";
	/** Encodable pattern loaded from {@link #MAP_TEXT}. */
	private static final String PATTERN="LLLDDD";
	/** Accession matching {@link #PATTERN}. */
	private static final String ACCESSION="ABC123";
	/** Accession whose pattern is absent from {@link #MAP_TEXT}. */
	private static final String MISSING="ABC";
	/** Expected compact key for {@link #ACCESSION}. */
	private static final long EXPECTED_KEY=112496;
	/** Linux procfs descriptor root. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");

	/** Finite valid pattern input whose close reports an IOException. */
	private static final class CloseFailureInput extends ByteArrayInputStream{

		/** Supplies the literal pattern table without a physical stdin file. */
		CloseFailureInput(){super(MAP_TEXT.getBytes(StandardCharsets.US_ASCII));}

		/** Reports the injected backend error on every actual close attempt.
		 * @throws IOException Always, to exercise TextFile's reported close status */
		@Override
		public void close() throws IOException{
			closeCalls++;
			throw new IOException("Injected AnalyzeAccession loadCodeMap close failure");
		}

		/** Number of actual backend close attempts. */
		int closeCalls=0;
	}

	/** Finite valid pattern input that fails on the first read past its prefix. */
	private static final class ReadFailureInput extends InputStream{

		/** Supplies a finite prefix before the injected read error. */
		private ReadFailureInput(){this(false);}

		/** Supplies a finite prefix and optional injected close failure.
		 * @param failOnClose_ True to report a close failure */
		private ReadFailureInput(boolean failOnClose_){
			failOnClose=failOnClose_;
		}

		/** Reads as many prefix bytes as possible, then reports the injected failure.
		 * @param b Destination buffer
		 * @param off Destination offset
		 * @param len Maximum bytes to read
		 * @return Prefix bytes copied
		 * @throws IOException After the finite prefix has been consumed */
		@Override
		public int read(byte[] b, int off, int len) throws IOException{
			if(len==0){return 0;}
			if(pos>=prefix.length){throw new IOException("Injected read failure");}
			final int copy=Math.min(len, prefix.length-pos);
			System.arraycopy(prefix, pos, b, off, copy);
			pos+=copy;
			return copy;
		}

		/** Reads one prefix byte, then reports the injected failure.
		 * @return Next unsigned prefix byte
		 * @throws IOException After the finite prefix has been consumed */
		@Override
		public int read() throws IOException{
			if(pos>=prefix.length){throw new IOException("Injected read failure");}
			return prefix[pos++]&0xFF;
		}

		/** Records and optionally fails TextFile's injected source close.
		 * @throws IOException When this instance should report a close error */
		@Override
		public void close() throws IOException{
			closeCalls++;
			if(failOnClose){
				throw new IOException("Injected AnalyzeAccession read-close failure");
			}
		}

		/** Valid prefix before the injected read failure. */
		final byte[] prefix=MAP_TEXT.getBytes(StandardCharsets.US_ASCII);
		/** Whether close should report a backend failure. */
		final boolean failOnClose;
		/** Next prefix byte to serve. */
		int pos=0;
		/** Number of actual backend close attempts. */
		int closeCalls=0;
	}

	/** Child process result. */
	private static final class RunResult{

		/** Captures process status and output.
		 * @param exit_ Child process exit code
		 * @param stdout_ Captured standard output
		 * @param stderr_ Captured standard error */
		RunResult(int exit_, String stdout_, String stderr_){
			exit=exit_;
			stdout=stdout_;
			stderr=stderr_;
		}

		/** Child process exit code. */
		final int exit;
		/** Captured standard output. */
		final String stdout;
		/** Captured standard error. */
		final String stderr;
	}
}
