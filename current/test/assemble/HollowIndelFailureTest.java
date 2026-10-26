package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Exercises writer shutdown through real launcher processes, not in-process exceptions.
 * GNU timeout kills the complete child process group if a regression leaves a writer alive.
 * @author Brian Bushnell, Fischl */
public final class HollowIndelFailureTest {

	/** Accepts a BBTools installation and optional old mode to reproduce the original hang. */
	public static void main(final String[] args) throws Exception{
		check(args.length==1 || args.length==2, "Specify BBTools path and optionally old.");
		final Path code=Paths.get(args[0]).toAbsolutePath();
		final Path root=Files.createTempDirectory("hollow-indel-failure-");
		final boolean old=args.length==2 && args[1].equals("old");
		final String bases=repeat('A', 250);
		final Path shortSeed=write(root, "short.fa", ">ref_0-249\n"+bases.substring(1)+"\n");
		run(code, root, "short", old, false, "Delete mode expects 250bp FASTA shreds.",
				"mode=delete", "in="+shortSeed);
		if(old){System.out.println("HOLLOW_INDEL_OLD_HANG_REPRODUCED dir="+root); return;}
		final Path wrongQ=write(root, "wrongq.fq", "@read\n"+bases+"\n+\n"+repeat('4', 250)+"\n");
		run(code, root, "quality", false, false, "All original qualities must be Q20.",
				"mode=readindels", "in="+wrongQ, "qin=33");
		final Path ref=write(root, "ref.fa", ">ref\n"+repeat('A', 500)+"\n");
		final Path seed=write(root, "seed.fa", ">ref_0-249\n"+bases+"\n");
		final Path changed=write(root, "changed.fa", ">ref_0-249\nC"+bases.substring(1)+"\n");
		run(code, root, "prefix", false, false, "Original seed was changed.",
				"mode=grade", "ref="+ref, "shreds="+seed, "in="+changed);
		final Path unknown=write(root, "unknown.fa", ">unknown\n"+bases+"\n");
		run(code, root, "id", false, false, "Unknown or duplicate output ID:",
				"mode=grade", "ref="+ref, "shreds="+seed, "in="+unknown);
		run(code, root, "good", false, true, "HOLLOW_INDEL_GRADE_PASS",
				"mode=grade", "ref="+ref, "shreds="+seed, "in="+seed);
		System.out.println("HOLLOW_INDEL_FAILURE_TEST_PASS cases=5 dir="+root);
	}

	/** Checks both timely nonzero exit and the original validation diagnostic. */
	private static void run(final Path code, final Path root, final String name,
			final boolean hangExpected, final boolean success, final String diagnostic,
			final String... args) throws Exception{
		assert(!hangExpected || !success) : "Only the frozen buggy build is expected to time out.";
		final Path log=root.resolve(name+".log");
		final List<String> command=new ArrayList<String>(Arrays.asList("timeout", "--kill-after=5s", "20s",
				"bash", code.resolve("dev/hollowindeltruth.sh").toString()));
		command.addAll(Arrays.asList(args));
		command.addAll(Arrays.asList("out="+root.resolve(name+".out"), "details="+root.resolve(name+".tsv"),
				"hist="+root.resolve(name+".hist"), "t=1", "interleaved=f", "-Xmx256m", "-Xms32m"));
		final Process child=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile()).start();
		child.getOutputStream().close();
		final int status=child.waitFor();
		final String text=new String(Files.readAllBytes(log), StandardCharsets.UTF_8);
		check(text.contains(diagnostic), "Missing primary diagnostic; inspect "+log);
		if(hangExpected){check(status==124, "Old build did not reproduce the writer hang: "+status);}
		else{
			check(status!=124 && status!=137, "Child timed out instead of closing writers: "+log);
			check((status==0)==success, "Unexpected exit status "+status+": "+log);
			if(!success){check(!text.contains("_PASS"), "Invalid input reported success: "+log);}
		}
	}

	/** Creates deterministic tiny ASCII fixtures in a unique temporary directory. */
	private static Path write(final Path root, final String name, final String text) throws Exception{
		assert(root!=null && name!=null) : "Fixture paths must belong to the test's temporary root.";
		return Files.write(root.resolve(name), text.getBytes(StandardCharsets.US_ASCII));
	}

	/** Builds small repetitive fixtures without relying on Java11 String.repeat. */
	private static String repeat(final char c, final int length){
		final char[] array=new char[length];
		Arrays.fill(array, c);
		return new String(array);
	}

	/** Retains failure diagnostics even if JVM assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
