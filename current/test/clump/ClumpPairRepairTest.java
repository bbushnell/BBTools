package test.clump;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** End-to-end grouped mate repair checks through clumpify.sh.
 * Baseline mode records existing behavior; fixed mode requires literal paired outputs
 * or a clear unmatched-read error with the original input retained.
 * @author Shinobu
 */
public final class ClumpPairRepairTest{

	/** Takes baseline/fixed, runtime directory, and a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		assert(args.length==3) : "Expected phase, runtime, and fresh fixture directory";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown phase: "+args[0];
		final Path runtime=Paths.get(args[1]), root=Paths.get(args[2]);
		Files.createDirectory(root);
		System.out.println("case\texit\texact_pairs\tclear_unmatched_error\tinput_retained");
		for(String version : new String[] {"v1", "v2", "v3"}){
			run(runtime, root, version, "complete", COMPLETE, fixed, false);
			run(runtime, root, version, "leading", ORPHAN_A+COMPLETE, fixed, true);
			run(runtime, root, version, "trailing", COMPLETE+ORPHAN_Z, fixed, true);
		}
		if(fixed){
			run(runtime, root, "v3", "boundary", COMPLETE, true, false);
		}
	}

	/** Saves the command and raw process outputs; each case has an isolated input and JVM. */
	private static void run(Path runtime, Path root, String version, String kind,
			String input, boolean fixed, boolean orphan) throws Exception{
		final String name=version+"-"+kind;
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path in=dir.resolve("input.fq"), one=dir.resolve("out1.fq"), two=dir.resolve("out2.fq");
		final byte[] original=input.getBytes(StandardCharsets.US_ASCII);
		Files.write(in, original, StandardOpenOption.CREATE_NEW);
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("clumpify.sh").toString());
		command.addAll(Arrays.asList("-Xmx512m", "-Xms64m", "t=1", "k=7", "seed=1", "hashes=4",
				"border=0", "qin=33", "qout=33", "rcomp=f", "ecco=f", "correct=f", "dedupe=f",
				"condense=f", "quantize=f", "reorder=f", "overwrite=f", "groups=3", "unpair=t",
				"repair=t", "namesort=t", "interleaved=f", "compresstemp=f", version+"=t",
				"deleteinput="+orphan, "in="+in, "out="+one, "out2="+two));
		if(kind.equals("boundary")){command.add("readbufferlength=3");}
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.environment().put("PATH", System.getProperty("java.home")+"/bin:"+builder.environment().get("PATH"));
		builder.redirectOutput(dir.resolve("stdout.txt").toFile());
		builder.redirectError(dir.resolve("stderr.txt").toFile());
		System.err.println("Running "+name+"; logs in "+dir);
		final int status=builder.start().waitFor();
		final boolean exact=equalsFile(one, LEFT) && equalsFile(two, RIGHT);
		final String error=new String(Files.readAllBytes(dir.resolve("stderr.txt")), StandardCharsets.UTF_8);
		final boolean clear=error.contains("Unmatched read") && error.contains("two mate-output files");
		final boolean retained=Files.exists(in) && Arrays.equals(original, Files.readAllBytes(in));
		System.out.println(name+'\t'+status+'\t'+exact+'\t'+clear+'\t'+retained);
		if(fixed){
			if(orphan){
				assert(status!=0 && clear) : "Missing clear unmatched-read rejection: "+name;
			}else{
				assert(status==0 && exact) : "Complete pairs were not restored exactly: "+name;
			}
			assert(retained) : "Original input changed or was deleted: "+name;
		}
	}

	/** Compares literal FASTQ records including order, names, bases and qualities. */
	private static boolean equalsFile(Path file, String expected) throws Exception{
		if(!Files.exists(file)){return false;}
		final List<String> actual=Files.readAllLines(file, StandardCharsets.US_ASCII);
		return actual.equals(Arrays.asList(expected.split("\n")));
	}

	private static final String A1="@alpha/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n";
	private static final String A2="@alpha/2\nTTTTCCCCAAAA\n+\nGGGGGGGGGGGG\n";
	private static final String B1="@beta/1\nTTGCAACCGGTT\n+\nHHHHHHHHHHHH\n";
	private static final String B2="@beta/2\nAAAACCCCTTTT\n+\nFFFFFFFFFFFF\n";
	private static final String COMPLETE=B2+A1+B1+A2;
	private static final String LEFT=A1+B1, RIGHT=A2+B2;
	private static final String ORPHAN_A="@aardvark/1\nACGTCGTA\n+\nIIIIIIII\n";
	private static final String ORPHAN_Z="@zulu/2\nGTCACAGT\n+\nHHHHHHHH\n";
}
