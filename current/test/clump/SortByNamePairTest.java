package test.clump;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

/** Checks the paired writer contract of clump's shared name-sort dependency.
 * Exercises direct writes and a forced small spill, with single-output controls.
 * @author Shinobu
 */
public final class SortByNamePairTest{

	/** Takes baseline/fixed, a prepared runtime directory, and a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		assert(args.length==3) : "Expected phase, runtime, and fresh fixture directory";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown phase: "+args[0];
		final Path runtime=Paths.get(args[1]).toAbsolutePath(), root=Paths.get(args[2]).toAbsolutePath();
		Files.createDirectory(root);
		System.out.println("case\texit\texact\tspill_observed\tinput_unchanged");
		for(boolean spill : new boolean[] {false, true}){
			for(boolean paired : new boolean[] {false, true}){
				run(runtime, root, paired, spill, fixed);
			}
		}
	}

	/** Runs a normal CLI invocation in its own directory and checks literal FASTQ output. */
	private static void run(Path runtime, Path root, boolean paired, boolean spill,
			boolean fixed) throws Exception{
		final String name=(paired ? "paired" : "single")+(spill ? "-spill" : "-direct");
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path left=dir.resolve("left.fq"), right=dir.resolve("right.fq");
		final Path one=dir.resolve("one.fq"), two=dir.resolve("two.fq");
		write(left, B1+A1); write(right, B2+A2);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("sortbyname.sh").toString());
		command.addAll(Arrays.asList("-Xmx256m", "-Xms32m", "t=1", "qin=33", "qout=33",
				"overwrite=f", "in="+left, "out="+one));
		if(paired){command.add("in2="+right); command.add("out2="+two);}
		else{command.add("interleaved=f");}
		if(spill){command.add("memmult=0.000001");}
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.directory(dir.toFile());
		builder.environment().put("PATH", System.getProperty("java.home")+"/bin:"+builder.environment().get("PATH"));
		builder.redirectOutput(dir.resolve("stdout.txt").toFile());
		builder.redirectError(dir.resolve("stderr.txt").toFile());
		System.err.println("Running "+name+"; logs in "+dir);
		final int status=builder.start().waitFor();
		final String stderr=new String(Files.readAllBytes(dir.resolve("stderr.txt")), StandardCharsets.UTF_8);
		final boolean spilled=stderr.contains("Merging 2 files.");
		final boolean exact=equalsFile(one, A1+B1) && (!paired || equalsFile(two, A2+B2));
		final boolean unchanged=equalsFile(left, B1+A1) && equalsFile(right, B2+A2);
		assert(spilled==spill) : "Unexpected sort path: "+name;
		assert(unchanged) : "Input changed: "+name;
		if(fixed || !paired){assert(status==0 && exact) : "Wrong sorted output: "+name;}
		else{
			assert(status!=0 && stderr.contains("WriterFactory") && stderr.contains("AssertionError"))
					: "Expected the old paired Writer descriptor rejection: "+name;
		}
		System.out.println(name+'\t'+status+'\t'+exact+'\t'+spilled+'\t'+unchanged);
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
}
