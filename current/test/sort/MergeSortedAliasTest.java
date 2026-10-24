package test.sort;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

/** Checks MergeSorted input/output alias and duplicate-input guards.
 * @author Jean
 */
public final class MergeSortedAliasTest{

	/** Takes baseline/followup-baseline/fixed, a prepared runtime directory, and a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		assert(args.length==3) : "Expected phase, runtime, and fresh fixture directory";
		final boolean fixed=args[0].equals("fixed");
		final boolean followupBaseline=args[0].equals("followup-baseline");
		assert(fixed || followupBaseline || args[0].equals("baseline")) : "Unknown phase: "+args[0];
		final Path runtime=Paths.get(args[1]).toAbsolutePath(), root=Paths.get(args[2]).toAbsolutePath();
		Files.createDirectory(root);
		System.out.println("case\texit\texact\tclear_error");
		runInputOutputAlias(runtime, root, fixed || followupBaseline);
		runDuplicateInput(runtime, root);
		runMaxFilesOne(runtime, root, fixed);
		runQualityOutput(runtime, root, fixed);
		runSeparateMateInput(runtime, root, fixed);
	}

	/** Rejects input/output aliasing before the writer can truncate the original input bytes. */
	private static void runInputOutputAlias(Path runtime, Path root, boolean fixed) throws Exception{
		final String name="input-output-alias";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("original.fq");
		write(input, ORIGINAL);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "t=1", "qin=33", "qout=33",
				"delete=f", "in="+input, "out="+input));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessResult result=run(command, dir);
		final boolean retained=equalsFile(input, ORIGINAL);
		final boolean clear=result.stderr.contains("specified multiple times");
		if(fixed){assert(result.status!=0 && retained && clear) : "Input/output alias was not safely rejected";}
		else{assert(!retained) : "Expected current MergeSorted to clobber aliased input/output";}
		System.out.println(name+'\t'+result.status+'\t'+retained+'\t'+clear);
	}

	/** Keeps existing duplicate-input rejection clear and output-free. */
	private static void runDuplicateInput(Path runtime, Path root) throws Exception{
		final String name="duplicate-input";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		write(input, ORIGINAL);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "t=1", "qin=33", "qout=33",
				"delete=f", "in="+input+","+input, "out="+output));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessResult result=run(command, dir);
		final boolean duplicated=equalsFile(output, DUPLICATED);
		final boolean clear=duplicateFileDiagnostic(result.stderr);
		assert(result.status!=0 && !Files.exists(output) && clear) :
			"Duplicate input was not rejected without writing output";
		System.out.println(name+'\t'+result.status+'\t'+duplicated+'\t'+clear);
	}

	/** Rejects maxfiles=1 before recursive merging can stall under disabled assertions. */
	private static void runMaxFilesOne(Path runtime, Path root, boolean fixed) throws Exception{
		final String name="maxfiles-one";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		write(input, ORIGINAL);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "-da", "t=1", "qin=33", "qout=33",
				"delete=f", "maxfiles=1", "in="+input, "out="+output));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessResult result=run(command, dir);
		final boolean exact=equalsFile(output, ORIGINAL);
		final boolean clear=result.stderr.contains("maxfiles must be at least 2");
		if(fixed){assert(result.status!=0 && !Files.exists(output) && clear) : "maxfiles=1 was not rejected";}
		else{assert(result.status==0 && exact) : "Expected current MergeSorted to accept maxfiles=1 for one input";}
		System.out.println(name+'\t'+result.status+'\t'+exact+'\t'+clear);
	}

	/** Rejects separate quality sidecar outputs because MergeSorted only writes FASTQ output streams. */
	private static void runQualityOutput(Path runtime, Path root, boolean fixed) throws Exception{
		final String name="quality-output";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		final Path qout=dir.resolve("out.qual");
		write(input, ORIGINAL);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "t=1", "qin=33", "qout=33",
				"delete=f", "in="+input, "out="+output, "qfout="+qout));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessResult result=run(command, dir);
		final boolean exact=equalsFile(output, ORIGINAL) && !Files.exists(qout);
		final boolean clear=result.stderr.contains("Separate quality output");
		if(fixed){assert(result.status!=0 && !Files.exists(output) && clear) : "qfout was not rejected";}
		else{assert(result.status==0 && exact) : "Expected current MergeSorted to ignore qfout";}
		System.out.println(name+'\t'+result.status+'\t'+exact+'\t'+clear);
	}

	/** Rejects in2 because MergeSorted consumes only the comma-delimited in list. */
	private static void runSeparateMateInput(Path runtime, Path root, boolean fixed) throws Exception{
		final String name="separate-mate-input";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input1=dir.resolve("input1.fq"), input2=dir.resolve("input2.fq");
		final Path output=dir.resolve("out.fq");
		write(input1, A);
		write(input2, B);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "t=1", "qin=33", "qout=33",
				"delete=f", "in="+input1, "in2="+input2, "out="+output));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessResult result=run(command, dir);
		final boolean exact=equalsFile(output, A);
		final boolean clear=result.stderr.contains("Separate in2 input is not supported");
		if(fixed){assert(result.status!=0 && !Files.exists(output) && clear) : "in2 was not rejected";}
		else{assert(result.status==0 && exact) : "Expected current MergeSorted to ignore in2";}
		System.out.println(name+'\t'+result.status+'\t'+exact+'\t'+clear);
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

	/** Accepts both current input-side and repaired general duplicate-file diagnostics. */
	private static boolean duplicateFileDiagnostic(String stderr){
		return stderr.contains("specified for multiple input streams") ||
			stderr.contains("specified multiple times");
	}

	private static void write(Path file, String text) throws Exception{
		Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	private static final String A="@alpha\nACGT\n+\nIIII\n";
	private static final String B="@beta\nTGCA\n+\nHHHH\n";
	private static final String ORIGINAL=A+B;
	private static final String DUPLICATED=A+A+B+B;

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
