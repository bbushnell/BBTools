package test.sort;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

/** Checks MergeSorted paired output descriptors and required input validation.
 * @author Jean
 */
public final class MergeSortedPairTest{

	/** Takes baseline/fixed, a prepared runtime directory, and a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		assert(args.length==3) : "Expected phase, runtime, and fresh fixture directory";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown phase: "+args[0];
		final Path runtime=Paths.get(args[1]).toAbsolutePath(), root=Paths.get(args[2]).toAbsolutePath();
		Files.createDirectory(root);
		System.out.println("case\texit\texact\tclear_error");
		runPairedMerge(runtime, root, fixed);
		runMissingInput(runtime, root, fixed);
	}

	/** Merges one sorted interleaved input into two ordered mate files. */
	private static void runPairedMerge(Path runtime, Path root, boolean fixed) throws Exception{
		final String name="paired-merge";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("sorted.fq"), one=dir.resolve("one.fq"), two=dir.resolve("two.fq");
		write(input, A1+A2+B1+B2);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "t=1", "qin=33", "qout=33",
				"overwrite=f", "int=t", "in="+input, "out="+one, "out2="+two));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessResult result=run(command, dir);
		final boolean exact=equalsFile(one, A1+B1) && equalsFile(two, A2+B2);
		final boolean descriptorAssertion=result.stderr.contains("WriterFactory") && result.stderr.contains("AssertionError");
		if(fixed){assert(result.status==0 && exact) : "Wrong paired MergeSorted output";}
		else{assert(result.status!=0 && descriptorAssertion) : "Expected paired Writer descriptor rejection";}
		System.out.println(name+'\t'+result.status+'\t'+exact+'\t'+descriptorAssertion);
	}

	/** Ensures a run with no input fails with the intended input-specific message. */
	private static void runMissingInput(Path runtime, Path root, boolean fixed) throws Exception{
		final String name="missing-input";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "t=1", "overwrite=f", "out="+dir.resolve("out.fq")));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessResult result=run(command, dir);
		final boolean clear=result.stderr.contains("at least one input file");
		assert(result.status!=0) : "Missing input unexpectedly succeeded";
		assert(clear==fixed) : "Unexpected missing-input failure mode: "+result.stderr;
		System.out.println(name+'\t'+result.status+"\tfalse\t"+clear);
	}

	/** Runs a BBTools shell command in a fixture-specific directory. */
	private static ProcessResult run(ArrayList<String> command, Path dir) throws Exception{
		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.directory(dir.toFile());
		builder.environment().put("PATH", System.getProperty("java.home")+"/bin:"+builder.environment().get("PATH"));
		builder.redirectOutput(dir.resolve("stdout.txt").toFile());
		builder.redirectError(dir.resolve("stderr.txt").toFile());
		final int status=builder.start().waitFor();
		final String stderr=new String(Files.readAllBytes(dir.resolve("stderr.txt")), StandardCharsets.UTF_8);
		return new ProcessResult(status, stderr);
	}

	/** Compares exact record order, identifiers, bases and qualities. */
	private static boolean equalsFile(Path file, String expected) throws Exception{
		return Files.exists(file) && Arrays.equals(Files.readAllBytes(file), expected.getBytes(StandardCharsets.US_ASCII));
	}

	private static void write(Path file, String text) throws Exception{
		Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	private static final String A1="@alpha/1\nACGT\n+\nIIII\n";
	private static final String A2="@alpha/2\nTGCA\n+\nHHHH\n";
	private static final String B1="@beta/1\nGATC\n+\nGGGG\n";
	private static final String B2="@beta/2\nCTAG\n+\nFFFF\n";

	/** Completed subprocess status and captured stderr. */
	private static final class ProcessResult{
		ProcessResult(int status_, String stderr_){
			status=status_;
			stderr=stderr_;
		}
		final int status;
		final String stderr;
	}
}
