package test.tax;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.TimeUnit;

import tax.AnalyzeAccession_ST;

/** Finite header and exceptional-reader oracles for AnalyzeAccession_ST.
 * @author Jean */
public final class AnalyzeAccessionSTInputTest{

	/** Runs ST header, MT control, and ST exceptional-reader expectations.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction, child launch or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		boolean child=false;
		Path input=null, output=null;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.equals("child=reader-leak")){child=true;}
			else if(arg.startsWith("input=")){input=Paths.get(arg.substring(6));}
			else if(arg.startsWith("output=")){output=Paths.get(arg.substring(7));}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		if(child){
			checkReaderLeakChild(input, output, baseline);
			return;
		}

		final Path root=Files.createTempDirectory("analyze-accession-st-");
		try{
			for(boolean assertions : new boolean[] {true, false}){
				checkHeader(root, assertions, "st", "tax.AnalyzeAccession_ST",
						baseline ? EXPECTED_BASELINE : EXPECTED_FIXED);
				checkHeader(root, assertions, "mt", "tax.AnalyzeAccession",
						EXPECTED_FIXED);
				checkReaderLeak(root, assertions, baseline);
			}
		}finally{
			deleteTempTree(root);
		}
		System.out.println("ANALYZE_ACCESSION_ST_INPUT_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private AnalyzeAccessionSTInputTest(){}

	/** Runs one header-counting child and checks exact pattern output.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param prefix Case prefix
	 * @param mainClass Child main class
	 * @param expected Exact expected output table
	 * @throws Exception If the child launch or file IO fails */
	private static void checkHeader(Path root, boolean assertions, String prefix,
			String mainClass, String expected) throws Exception{
		final String label=prefix+"-header-"+(assertions ? "ea" : "da");
		final Path input=root.resolve(label+".tsv");
		final Path output=root.resolve(label+".out");
		Files.write(input, HEALTHY.getBytes(StandardCharsets.US_ASCII));

		RunResult result=runAnalyze(root, label, mainClass, assertions, input, output);
		check(result.exit==0, label+" child exited "+result.exit+
				" with stderr\n"+result.stderr);
		final String actual=readAscii(output);
		check(expected.equals(actual),
				label+" produced unexpected output\n"+actual);
		System.out.println("ANALYZE_ACCESSION_ST_HEADER_CASE_PASS case="+label+
				" bytes="+actual.length());
	}

	/** Runs an ST child that catches a non-ASCII row failure before process exit.
	 * @param root Temporary artifact directory
	 * @param assertions True for -ea, false for -da
	 * @param baseline True to expect the current open reader
	 * @throws Exception If the child launch or file IO fails */
	private static void checkReaderLeak(Path root, boolean assertions,
			boolean baseline) throws Exception{
		final String label="reader-leak-"+(assertions ? "ea" : "da");
		final Path input=root.resolve(label+".tsv");
		final Path output=root.resolve(label+".out");
		Files.write(input, poisonInput());

		RunResult result=runSelf(root, label, assertions, input, output, baseline);
		check(result.exit==0, label+" child exited "+result.exit+
				" with stderr\n"+result.stderr);
		check(result.stdout.contains("ANALYZE_ACCESSION_ST_READER_LEAK_CASE_PASS"),
				label+" child did not report its reader-leak pass\n"+result.stdout+
				"\nstderr\n"+result.stderr);
		System.out.println("ANALYZE_ACCESSION_ST_PARENT_CHILD_PASS case="+label+
				" "+result.stdout.trim());
	}

	/** Runs the exceptional-reader case inside a child JVM to inspect live FDs.
	 * @param input Input file that will throw from AnalyzeAccession_ST.increment
	 * @param output Output file requested from AnalyzeAccession_ST
	 * @param baseline True to expect the current open reader
	 * @throws Exception If the oracle fails or procfs cannot be read */
	private static void checkReaderLeakChild(Path input, Path output,
			boolean baseline) throws Exception{
		check(input!=null, "Missing input path");
		check(output!=null, "Missing output path");
		check(countDescriptors(input)==0, "Input descriptor open before ST child call");

		Throwable failure=null;
		try{
			AnalyzeAccession_ST.main(analyzeArgs(input, output));
		}catch(Throwable t){
			failure=t;
		}
		check(failure instanceof ArrayIndexOutOfBoundsException,
				"Original increment failure was not preserved: "+failure);
		check(hasFrame(failure, "tax.AnalyzeAccession_ST", "increment"),
				"Failure did not originate in AnalyzeAccession_ST.increment: "+
				stackSummary(failure));

		final int descriptors=countDescriptors(input);
		final int expected=baseline ? 1 : 0;
		check(descriptors==expected,
				"Expected "+expected+" live input descriptor(s); found "+descriptors);
		check(!Files.exists(output),
				"Non-ASCII input should fail before pattern output is written");

		System.out.println("ANALYZE_ACCESSION_ST_READER_LEAK_CASE_PASS read_fds="+
				descriptors+" expect="+(baseline ? "baseline" : "fixed"));
	}

	/** Runs AnalyzeAccession or AnalyzeAccession_ST in a child JVM.
	 * @param root Temporary artifact directory
	 * @param label Case label for captured stdout/stderr
	 * @param mainClass Child class
	 * @param assertions True for -ea, false for -da
	 * @param input Input TSV file
	 * @param output Requested pattern summary path
	 * @return Child exit status and captured output
	 * @throws Exception If process launch, wait or file IO fails */
	private static RunResult runAnalyze(Path root, String label, String mainClass,
			boolean assertions, Path input, Path output) throws Exception{
		final ArrayList<String> command=baseCommand(assertions);
		command.add(mainClass);
		command.addAll(Arrays.asList(analyzeArgs(input, output)));
		return run(root, label, command);
	}

	/** Runs this fixture in its descriptor-inspecting child mode.
	 * @param root Temporary artifact directory
	 * @param label Case label for captured stdout/stderr
	 * @param assertions True for -ea, false for -da
	 * @param input Input TSV file
	 * @param output Requested pattern summary path
	 * @param baseline True to expect the current open reader
	 * @return Child exit status and captured output
	 * @throws Exception If process launch, wait or file IO fails */
	private static RunResult runSelf(Path root, String label, boolean assertions,
			Path input, Path output, boolean baseline) throws Exception{
		final ArrayList<String> command=baseCommand(assertions);
		command.add(AnalyzeAccessionSTInputTest.class.getName());
		command.add(baseline ? "expect=baseline" : "expect=fixed");
		command.add("child=reader-leak");
		command.add("input="+input.toString());
		command.add("output="+output.toString());
		return run(root, label, command);
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
		child.getOutputStream().close();
		try{
			check(child.waitFor(15, TimeUnit.SECONDS),
					label+" child did not finish within its bounded wait");
			return new RunResult(child.exitValue(), readAscii(stdout), readAscii(stderr));
		}finally{
			stop(child);
		}
	}

	/** Builds shared AnalyzeAccession command-line arguments.
	 * @param input Input file
	 * @param output Requested pattern output
	 * @return Arguments forcing tiny regular-file BF1 execution */
	private static String[] analyzeArgs(Path input, Path output){
		return new String[] {"in="+input.toString(), "out="+output.toString(),
				"outstream=null", "ow=t", "t=1", "bf1=t", "unpigz=f",
				"gunzip=f", "unbgzip=f", "simd=f"};
	}

	/** Builds an accession2taxid input with a non-ASCII row after the header.
	 * @return Input bytes whose poison row fails AnalyzeAccession_ST.increment */
	private static byte[] poisonInput(){
		final byte[] prefix=HEALTHY.getBytes(StandardCharsets.US_ASCII);
		final byte[] suffix="\tBAD.1\t8\t0\n".getBytes(StandardCharsets.US_ASCII);
		final byte[] result=new byte[prefix.length+1+suffix.length];
		System.arraycopy(prefix, 0, result, 0, prefix.length);
		result[prefix.length]=(byte)0x80;
		System.arraycopy(suffix, 0, result, prefix.length+1, suffix.length);
		return result;
	}

	/** Counts descriptors for exactly the owned regular input file.
	 * @param input Input path
	 * @return Open descriptor count
	 * @throws IOException If procfs cannot be inspected */
	private static int countDescriptors(final Path input) throws IOException{
		final Path real=input.toRealPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(
				Paths.get("/proc/self/fd"))){
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
			check(process.waitFor(3, TimeUnit.SECONDS),
					"Owned child did not terminate after forced stop");
		}
	}

	/** Tests whether a throwable stack includes an exact source frame.
	 * @param failure Throwable to inspect
	 * @param className Expected class name
	 * @param methodName Expected method name
	 * @return True if the stack includes the expected frame */
	private static boolean hasFrame(Throwable failure, String className,
			String methodName){
		for(StackTraceElement element : failure.getStackTrace()){
			if(className.equals(element.getClassName())
					&& methodName.equals(element.getMethodName())){
				return true;
			}
		}
		return false;
	}

	/** Builds a compact stack summary for failed source-frame checks.
	 * @param failure Throwable to summarize
	 * @return First few stack frames, or null for a missing throwable */
	private static String stackSummary(Throwable failure){
		if(failure==null){return "null";}
		StringBuilder sb=new StringBuilder(failure.toString());
		StackTraceElement[] stack=failure.getStackTrace();
		for(int i=0; i<stack.length && i<4; i++){
			sb.append(" <- ").append(stack[i].getClassName()).append('.')
					.append(stack[i].getMethodName());
		}
		return sb.toString();
	}

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

	/** Accession2taxid input with one header row and one valid accession row. */
	private static final String HEALTHY=
			"accession\taccession.version\ttaxid\tgi\n"+
			"ABC123\tABC123.1\t2\t3\n";
	/** Fixed pattern table produced by {@link #HEALTHY}. */
	private static final String EXPECTED_FIXED=
			"#Pattern\tCount\tCombos\tBits\n"+
			"LLLDDD\t1\t17576000\t24.07\n";
	/** Current AnalyzeAccession_ST table when it counts the header token. */
	private static final String EXPECTED_BASELINE=
			"#Pattern\tCount\tCombos\tBits\n"+
			"LLLLLLLLL\t1\t5429503678976\t42.30\n"+
			"LLLDDD\t1\t17576000\t24.07\n";

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
