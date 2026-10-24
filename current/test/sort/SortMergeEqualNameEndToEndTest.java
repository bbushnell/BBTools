package test.sort;

import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import fileIO.FileFormat;
import shared.Shared;
import sort.ReadComparatorName;
import sort.SortByName;
import stream.FASTQ;

/** End-to-end equal-name merge boundary that leaves an equal queued pivot.
 * @author Jean
 */
public final class SortMergeEqualNameEndToEndTest{

	/** Runs expect=baseline|fixed or child in a fresh test JVM. */
	public static void main(String[] args) throws Exception{
		if(args.length==2 && args[0].equals("child")){
			runChild(Paths.get(args[1]).toAbsolutePath());
			return;
		}
		if(args.length!=2 || (!args[0].equals("expect=baseline") && !args[0].equals("expect=fixed"))){
			throw new IllegalArgumentException("Expected expect=baseline|expect=fixed and a fresh fixture directory");
		}
		final boolean fixed=args[0].equals("expect=fixed");
		final Path root=Paths.get(args[1]).toAbsolutePath();
		Files.createDirectory(root);
		final Path childDir=root.resolve("child");

		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add("-ea");
		command.add("-Xmx128m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		if(Shared.javaVersion>=16){
			command.add("--add-modules");
			command.add("jdk.incubator.vector");
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(SortMergeEqualNameEndToEndTest.class.getName());
		command.add("child");
		command.add(childDir.toString());

		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.redirectOutput(root.resolve("child.stdout").toFile());
		builder.redirectError(root.resolve("child.stderr").toFile());
		final Process process=builder.start();
		final boolean finished=process.waitFor(10, TimeUnit.SECONDS);
		if(!finished){
			process.destroyForcibly();
			process.waitFor();
			throw new AssertionError("Child JVM timed out.");
		}
		final int exit=process.exitValue();
		final String stderr=new String(Files.readAllBytes(root.resolve("child.stderr")),
			StandardCharsets.UTF_8);
		final boolean exact=exit==0 && exactFastq(childDir.resolve("out.fq"), 100002);
		if(fixed){
			assert(exact) : "Expected fixed child to finish with exact FASTQ output; exit="+exit+
				" stderr="+stderr;
		}else{
			assert(exit==BASELINE_ASSERTION_EXIT && stderr.contains("AssertionError") &&
					stderr.contains("SortByName.indexOfLowestAbovePivot") &&
					stderr.contains("SortByName.mergeAndDump")) :
				"Expected child to exit from the known equal-pivot assertion; exit="+exit+
				" stderr="+stderr;
		}
		System.out.println("equal-name-end-to-end	"+exit+'	'+exact);
	}

	/** Runs the merge in a disposable child JVM. */
	private static void runChild(Path root) throws Exception{
		Files.createDirectory(root);
		Shared.setThreads(1);
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.FORCE_INTERLEAVED=false;
		final Path left=root.resolve("left.fq"), right=root.resolve("right.fq");
		final Path out=root.resolve("out.fq");
		write(left, 100001);
		write(right, 1);
		final ArrayList<String> names=new ArrayList<String>();
		names.add(left.toString());
		names.add(right.toString());

		SortByName.setComparator(ReadComparatorName.ascending);
		final FileFormat ffout=FileFormat.testOutput(out.toString(), FileFormat.FASTQ,
			null, true, false, false, false);
		try{
			final boolean errorState=SortByName.mergeAndDump(names, ffout, null, false,
				false, false, System.err, 4);
			if(errorState){throw new AssertionError("Merge reported an error state.");}
		}catch(AssertionError e){
			e.printStackTrace(System.err);
			System.exit(BASELINE_ASSERTION_EXIT);
		}
	}

	/** Requires count exact one-base records with the shared equal name. */
	private static boolean exactFastq(Path file, int count) throws Exception{
		if(!Files.exists(file)){return false;}
		final List<String> lines=Files.readAllLines(file, StandardCharsets.US_ASCII);
		if(lines.size()!=count*4){return false;}
		for(int i=0; i<count; i++){
			final int offset=i*4;
			if(!lines.get(offset).equals("@same") ||
					!lines.get(offset+1).equals("A") ||
					!lines.get(offset+2).equals("+") ||
					!lines.get(offset+3).equals("I")){
				return false;
			}
		}
		return true;
	}

	/** Writes count identical-name one-base FASTQ records. */
	private static void write(Path file, int count) throws Exception{
		try(BufferedWriter writer=Files.newBufferedWriter(file, StandardCharsets.US_ASCII)){
			for(int i=0; i<count; i++){
				writer.write("@same");
				writer.newLine();
				writer.write("A");
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write("I");
				writer.newLine();
			}
		}
	}

	private static final int BASELINE_ASSERTION_EXIT=7;
}
