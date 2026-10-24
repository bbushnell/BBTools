package test.sort;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import sort.ReadComparatorRandom;
import sort.SortByName;
import stream.FASTQ;

/** Focused regression coverage for random sort modes whose keys are not serialized to spills.
 * @author Jean
 */
public final class SortRandomModeGuardTest{

	/** Runs baseline|fixed against a prepared runtime directory and a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		if(args.length>0 && args[0].equals("static-merge")){
			runStaticMergeChild(args);
			return;
		}
		assert(args.length==3) : "Expected phase, runtime, and fresh fixture directory";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown phase: "+args[0];
		final Path runtime=Paths.get(args[1]).toAbsolutePath(), root=Paths.get(args[2]).toAbsolutePath();
		Files.createDirectory(root);
		System.out.println("case\texit\tcomplete\tclear_error\tinput_exact");
		runSortByNameRandomSpill(runtime, root, fixed, true,
				"sortbyname-random-temp-ea", "random=t");
		runSortByNameRandomSpill(runtime, root, fixed, false,
				"sortbyname-shuffle-temp-da", "shuffle=t");
		runSortByNameRandomSpill(runtime, root, fixed, true,
				"sortbyname-random-then-name", "random=t", "name=t");
		runSortByNameRandomSpill(runtime, root, fixed, false,
				"sortbyname-name-then-shuffle", "name=t", "shuffle=t");
		runMergeSortedRandom(runtime, root, fixed, true, "mergesorted-random-ea",
				"random=t");
		runMergeSortedRandom(runtime, root, fixed, false, "mergesorted-shuffle-da",
				"shuffle=t");
		runMergeSortedRandom(runtime, root, fixed, true, "mergesorted-random-then-name",
				"random=t", "name=t");
		runMergeSortedRandom(runtime, root, fixed, false, "mergesorted-name-then-shuffle",
				"name=t", "shuffle=t");
		runStaticRandomMerge(root, fixed, true);
		runStaticRandomMerge(root, fixed, false);
		runSortByNameInMemoryRandom(runtime, root, true);
		runSortByNameInMemoryRandom(runtime, root, false);
	}

	/** Rejects random SortByName only when temp spills could reread unserialized random keys. */
	private static void runSortByNameRandomSpill(Path runtime, Path root, boolean fixed,
			boolean assertions, String name, String... sortFlags) throws Exception{
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		write(input, INPUT);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("sortbyname.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", assertions ? "-ea" : "-da",
				"t=1", "qin=33", "qout=33", "delete=f", "allowtemp=t",
				"in="+input, "out="+output));
		command.addAll(Arrays.asList(sortFlags));
		final ProcessResult result=run(command, dir);
		final boolean complete=sameRecords(output, INPUT);
		final boolean clear=randomModeDiagnostic(result.stderr);
		final boolean inputExact=equalsFile(input, INPUT);
		final boolean shouldReject=isRandomLast(sortFlags);
		if(fixed && shouldReject){assert(result.status!=0 && !Files.exists(output) &&
				inputExact && clear) : name+" was not rejected before output";}
		else{assert(result.status==0 && complete) : "Expected current random SortByName to run";}
		System.out.println(name+'\t'+result.status+'\t'+complete+'\t'+clear+'\t'+inputExact);
	}

	/** Rejects MergeSorted random mode because every merge reads materialized records. */
	private static void runMergeSortedRandom(Path runtime, Path root, boolean fixed,
			boolean assertions, String name, String... sortFlags) throws Exception{
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		write(input, INPUT);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("mergesorted.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", assertions ? "-ea" : "-da",
				"t=1", "qin=33", "qout=33", "in="+input, "out="+output));
		command.addAll(Arrays.asList(sortFlags));
		final ProcessResult result=run(command, dir);
		final boolean complete=sameRecords(output, INPUT);
		final boolean clear=randomModeDiagnostic(result.stderr);
		final boolean inputExact=equalsFile(input, INPUT);
		final boolean shouldReject=isRandomLast(sortFlags);
		if(fixed && shouldReject){assert(result.status!=0 && !Files.exists(output) &&
				inputExact && clear) : name+" was not rejected before output";}
		else{assert(result.status==0 && complete) : "Expected current random MergeSorted to run";}
		System.out.println(name+'\t'+result.status+'\t'+complete+'\t'+clear+'\t'+inputExact);
	}

	/** Rejects direct public static mergeAndDump access after setComparator(Random). */
	private static void runStaticRandomMerge(Path root, boolean fixed, boolean assertions)
			throws Exception{
		final String name="sortbyname-static-random-merge-"+(assertions ? "ea" : "da");
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		write(input, INPUT);
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.addAll(Arrays.asList(assertions ? "-ea" : "-da", "-Xmx128m", "-Xms128m",
				"-XX:+ExitOnOutOfMemoryError"));
		if(Shared.javaVersion>=16){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
		command.addAll(Arrays.asList("-cp", System.getProperty("java.class.path"),
				SortRandomModeGuardTest.class.getName(), "static-merge", input.toString(),
				output.toString()));
		final ProcessResult result=run(command, dir);
		final boolean complete=sameRecords(output, INPUT);
		final boolean clear=randomModeDiagnostic(result.stderr);
		final boolean inputExact=equalsFile(input, INPUT);
		if(fixed){assert(result.status!=0 && !Files.exists(output) && inputExact && clear) :
			"Static random merge was not rejected before output";}
		else{assert(result.status==0 && complete && !Files.exists(input)) :
			"Expected current static random merge to write and delete input";}
		System.out.println(name+'\t'+result.status+'\t'+complete+'\t'+clear+'\t'+inputExact);
	}

	/** Child entry point that directly bypasses the command constructors. */
	private static void runStaticMergeChild(String[] args){
		assert(args.length==3) : "Expected static-merge, input, output";
		Shared.setThreads(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT=33;
		final ArrayList<String> names=new ArrayList<String>();
		names.add(args[1]);
		SortByName.setComparator(ReadComparatorRandom.comparator);
		final FileFormat ffout=FileFormat.testOutput(args[2], FileFormat.FASTQ, null, true,
				true, false, false);
		final boolean error=SortByName.mergeAndDump(names, ffout, null, true, false, false,
				System.err, 1000);
		if(error){throw new RuntimeException("Static random merge returned error state");}
	}

	/** Keeps bounded in-memory random sorting available. */
	private static void runSortByNameInMemoryRandom(Path runtime, Path root, boolean assertions)
			throws Exception{
		final String name="sortbyname-random-memory-"+(assertions ? "ea" : "da");
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		write(input, INPUT);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("sortbyname.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", assertions ? "-ea" : "-da",
				"t=1", "qin=33", "qout=33", "delete=f", "random=t", "allowtemp=f",
					"in="+input, "out="+output));
		final ProcessResult result=run(command, dir);
		final boolean complete=sameRecords(output, INPUT);
		assert(result.status==0 && complete) : "In-memory random SortByName was not retained";
		System.out.println(name+'\t'+result.status+'\t'+complete+"\tfalse\t"+equalsFile(input, INPUT));
	}

	/** Runs a BBTools shell command in a fixture-specific directory. */
	private static ProcessResult run(ArrayList<String> command, Path dir) throws Exception{
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8,
				StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.directory(dir.toFile());
		builder.environment().put("JAVA_HOME", System.getProperty("java.home"));
		builder.environment().put("PATH", System.getProperty("java.home")+"/bin:"+
				builder.environment().get("PATH"));
		builder.redirectOutput(dir.resolve("stdout.txt").toFile());
		builder.redirectError(dir.resolve("stderr.txt").toFile());
		final int status=builder.start().waitFor();
		final String stderr=new String(Files.readAllBytes(dir.resolve("stderr.txt")),
				StandardCharsets.UTF_8);
		return new ProcessResult(status, stderr);
	}

	/** Checks exact FASTQ record identity independent of output order. */
	private static boolean sameRecords(Path file, String expected) throws Exception{
		if(!Files.exists(file)){return false;}
		final ArrayList<String> a=records(new String(Files.readAllBytes(file),
				StandardCharsets.US_ASCII));
		final ArrayList<String> b=records(expected);
		Collections.sort(a);
		Collections.sort(b);
		return a.equals(b);
	}

	/** Checks exact file bytes. */
	private static boolean equalsFile(Path file, String expected) throws Exception{
		return Files.exists(file) && Arrays.equals(Files.readAllBytes(file),
				expected.getBytes(StandardCharsets.US_ASCII));
	}

	/** Splits a tiny FASTQ fixture into full record strings. */
	private static ArrayList<String> records(String text){
		final String[] lines=text.split("\n", -1);
		assert(lines.length%4==1 && lines[lines.length-1].isEmpty()) :
			"FASTQ fixture should have complete 4-line records";
		final ArrayList<String> list=new ArrayList<String>();
		for(int i=0; i<lines.length-1; i+=4){
			list.add(lines[i]+'\n'+lines[i+1]+'\n'+lines[i+2]+'\n'+lines[i+3]+'\n');
		}
		return list;
	}

	/** Accepts the fixed message only if it directs users to the supported shuffler. */
	private static boolean randomModeDiagnostic(String stderr){
		final String lower=stderr.toLowerCase();
		return lower.contains("random") && lower.contains("shuffle2.sh");
	}

	/** Returns whether the final comparator flag in this case selects random ordering. */
	private static boolean isRandomLast(String... sortFlags){
		for(int i=sortFlags.length-1; i>=0; i--){
			final String flag=sortFlags[i].toLowerCase();
			if(flag.startsWith("random=") || flag.startsWith("shuffle=")){return true;}
			if(flag.startsWith("name=")){return false;}
		}
		throw new AssertionError("Case lacks a relevant comparator flag");
	}

	private static void write(Path file, String text) throws Exception{
		Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	private static final String A="@alpha\nACGT\n+\nIIII\n";
	private static final String B="@beta\nTGCA\n+\nHHHH\n";
	private static final String C="@gamma\nGATC\n+\nJJJJ\n";
	private static final String INPUT=A+B+C;

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
