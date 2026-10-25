package test.synth;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

/** Captures RandomReadsMG variance as whole-contig depth variation in both wave modes.
 * @author Jean
 */
public final class RandomReadsMGVarianceTest{

	private static final int CONTIGS=4;
	private static final int CONTIG_LENGTH=5200;
	private static final int READ_LENGTH=100;
	private static final int BASELINE_READS=52;
	private static final int EXPLICIT_READS=7;

	/** Runs the captured baseline or the corrected whole-contig variance contract. */
	public static void main(String[] args) throws Exception{
		check(args.length==3 || args.length==4,
				"Expected <BBTools dir> expect=baseline|expect=fixed <results dir> [baseline dir]");
		final boolean fixed=args[1].equals("expect=fixed");
		check(fixed || args[1].equals("expect=baseline"), "Unknown expectation: "+args[1]);
		final Path code=Paths.get(args[0]).toAbsolutePath();
		final Path root=Paths.get(args[2]).toAbsolutePath();
		final Path baselineRoot=(args.length==4 ? Paths.get(args[3]).toAbsolutePath() : null);
		Files.createDirectories(root);
		final Path ref=root.resolve("ref.fa");
		writeReference(ref);

		final Result noWave0=run(code, root, ref, "nowave-0", false, "0", "-1");
		final Result noWave0Repeat=run(code, root, ref, "nowave-0-repeat", false, "0", "-1");
		final Result noWaveV=run(code, root, ref, "nowave-v", false, "0.5", "-1");
		final Result noWaveExplicit=run(code, root, ref, "nowave-rpc7-v", false, "0.5", "7");
		final Result sine0=run(code, root, ref, "sine-0", true, "0", "-1");
		final Result sine0Repeat=run(code, root, ref, "sine-0-repeat", true, "0", "-1");
		final Result sineV=run(code, root, ref, "sine-v", true, "0.5", "-1");
		final Result sineExplicit=run(code, root, ref, "sine-rpc7-v", true, "0.5", "7");

		check(allEqual(noWave0.counts, BASELINE_READS),
				"variance=0 without sinewave must generate equal contig targets: "
						+Arrays.toString(noWave0.counts));
		check(allEqual(sine0.counts, BASELINE_READS),
				"variance=0 with sinewave must generate equal contig targets: "
						+Arrays.toString(sine0.counts));
		check(Arrays.equals(noWave0.bytes, noWave0Repeat.bytes),
				"variance=0 without sinewave must be deterministic for the pinned seed.");
		check(Arrays.equals(sine0.bytes, sine0Repeat.bytes),
				"variance=0 with sinewave must be deterministic for the pinned seed.");
		check(allEqual(noWaveExplicit.counts, EXPLICIT_READS),
				"Explicit readspercontig must ignore nonzero variance without sinewave: "
						+Arrays.toString(noWaveExplicit.counts));
		check(allEqual(sineExplicit.counts, EXPLICIT_READS),
				"Explicit readspercontig must ignore nonzero variance with sinewave: "
						+Arrays.toString(sineExplicit.counts));

		if(fixed){
			check(allWithin(noWaveV.counts, 26, 78),
					"Fixed non-sinewave counts escaped the variance bounds: "
							+Arrays.toString(noWaveV.counts));
			check(allWithin(sineV.counts, 26, 78),
					"Fixed sinewave counts escaped the variance bounds: "+Arrays.toString(sineV.counts));
			check(!allEqual(noWaveV.counts, BASELINE_READS),
					"Fixed non-sinewave variance must alter whole-contig targets.");
			check(!allEqual(sineV.counts, BASELINE_READS),
					"Fixed sinewave variance must alter whole-contig targets.");
			if(baselineRoot!=null){
				check(Arrays.equals(noWave0.bytes, readBaseline(baselineRoot, "nowave-0")),
						"Fixed non-sinewave variance=0 output changed from the Dori baseline.");
				check(Arrays.equals(sine0.bytes, readBaseline(baselineRoot, "sine-0")),
						"Fixed sinewave variance=0 output changed from the Dori baseline.");
			}
		}else{
			check(allEqual(noWaveV.counts, BASELINE_READS),
					"Baseline non-sinewave variance should not alter contig targets: "
							+Arrays.toString(noWaveV.counts));
			check(!Arrays.equals(noWave0.bytes, noWaveV.bytes),
					"Baseline non-sinewave variance should still perturb per-read positions.");
			check(Arrays.equals(sine0.bytes, sineV.bytes),
					"Baseline sinewave variance should currently be ignored by CoverageModel.");
		}
		System.out.println("RANDOMREADSMG_VARIANCE_TEST_PASS expect="+(fixed ? "fixed" : "baseline"));
	}

	/** Generates a deterministic equal-length four-contig FASTA. */
	private static void writeReference(final Path ref) throws Exception{
		final StringBuilder sb=new StringBuilder(CONTIGS*(CONTIG_LENGTH+20));
		for(int i=0; i<CONTIGS; i++){
			sb.append(">contig").append(i).append('\n');
			for(int j=0; j<CONTIG_LENGTH; j++){sb.append("ACGT".charAt(j&3));}
			sb.append('\n');
		}
		Files.write(ref, sb.toString().getBytes(StandardCharsets.US_ASCII));
	}

	/** Runs randomreadsmg.sh once with deterministic IO and returns exact output facts. */
	private static Result run(final Path code, final Path root, final Path ref,
			final String name, final boolean sinewave, final String variance,
			final String readsPerContig) throws Exception{
		final Path dir=root.resolve(name);
		Files.createDirectories(dir);
		final Path out=dir.resolve("out.fq");
		final Path log=dir.resolve("child.log");
		final Path pgid=dir.resolve("child.pgid");
		final Path launcher=dir.resolve("launch.sh");
		Files.deleteIfExists(pgid);
		Files.write(launcher, ("#!/bin/bash\nset -u\nsetsid \"$@\" &\nchild=$!\n"+
				"printf '%s\\n' \"$child\" > \"$PGID_FILE\"\nwait \"$child\"\n")
				.getBytes(StandardCharsets.US_ASCII));

		final List<String> command=new ArrayList<String>(Arrays.asList("bash",
				launcher.toString(), "bash", code.resolve("randomreadsmg.sh").toString(),
				"in="+ref, "out="+out, "depth=1", "variance="+variance, "seed=1", "tree=f",
				"paired=f", "length="+READ_LENGTH, "adderrors=f", "qrange=0",
				"readspercontig="+readsPerContig, "t=1", "sinewave="+(sinewave ? "t" : "f"),
				"ow=t", "usepigz=f", "useunpigz=f", "eoom",
				"-Xmx256m", "-Xms128m", "-ea"));
		if(sinewave){command.addAll(Arrays.asList("waves=1", "waveamp=0.7"));}
		final ProcessBuilder pb=makeChild(command, log);
		pb.environment().put("PGID_FILE", pgid.toString());
		final Process child=pb.start();
		child.getOutputStream().close();
		final boolean finished=child.waitFor(10, TimeUnit.SECONDS);
		if(!finished){
			terminateOwnedGroup(child, pgid, log);
			throw new AssertionError("randomreadsmg.sh timed out; log="+log);
		}
		check(child.exitValue()==0, "randomreadsmg.sh failed with exit "+child.exitValue()+"; log="+log);
		final byte[] bytes=Files.readAllBytes(out);
		final int[] counts=countContigs(bytes);
		System.out.println(name+"\t"+Arrays.toString(counts));
		return new Result(bytes, counts);
	}

	/** Reads the Dori baseline output for an exact cross-revision variance=0 check. */
	private static byte[] readBaseline(final Path baselineRoot, final String name) throws Exception{
		return Files.readAllBytes(baselineRoot.resolve(name).resolve("out.fq"));
	}

	/** Kills only the process group started by this test case's setsid launcher. */
	private static void terminateOwnedGroup(final Process child, final Path pgid, final Path log) throws Exception{
		for(int i=0; i<20 && !Files.exists(pgid); i++){Thread.sleep(50);}
		check(Files.exists(pgid), "Watchdog could not find the owned child process group: "+log);
		final String group=new String(Files.readAllBytes(pgid), StandardCharsets.US_ASCII).trim();
		check(group.matches("[0-9]+"), "Invalid owned process group id "+group+" in "+pgid);
		new ProcessBuilder("kill", "-TERM", "--", "-"+group).start().waitFor();
		if(!waitForGroupExit(group)){
			new ProcessBuilder("kill", "-KILL", "--", "-"+group).start().waitFor();
			child.destroyForcibly();
			check(waitForGroupExit(group), "Owned RandomReadsMG group survived the watchdog: "+log);
		}
		child.waitFor(2, TimeUnit.SECONDS);
		check(!child.isAlive(), "Owned RandomReadsMG child did not die after the watchdog: "+log);
	}

	/** Waits briefly for the process group owned by this fixture to disappear. */
	private static boolean waitForGroupExit(final String group) throws Exception{
		for(int i=0; i<20; i++){
			if(!groupAlive(group)){return true;}
			Thread.sleep(100);
		}
		return !groupAlive(group);
	}

	/** Returns true if any process remains in the owned child process group. */
	private static boolean groupAlive(final String group) throws Exception{
		return new ProcessBuilder("kill", "-0", "--", "-"+group).start().waitFor()==0;
	}

	/** Counts FASTQ headers by RandomReadsMG's c_N source-contig tag. */
	private static int[] countContigs(final byte[] fastq){
		final int[] counts=new int[CONTIGS];
		final String[] lines=new String(fastq, StandardCharsets.US_ASCII).split("\n");
		check(lines.length%4==0, "FASTQ line count must be divisible by four: "+lines.length);
		for(int i=0; i<lines.length; i+=4){
			final String header=lines[i];
			final int start=header.indexOf("_c_")+3;
			check(start>=3, "Missing c_ contig tag in header: "+header);
			final int end=header.indexOf("_s_", start);
			check(end>start, "Missing s_ strand tag after c_ in header: "+header);
			final int contig=Integer.parseInt(header.substring(start, end));
			check(contig>=0 && contig<CONTIGS, "Unexpected contig "+contig+" in header: "+header);
			counts[contig]++;
		}
		return counts;
	}

	private static boolean allEqual(final int[] counts, final int expected){
		for(final int count : counts){if(count!=expected){return false;}}
		return true;
	}

	private static boolean allWithin(final int[] counts, final int min, final int max){
		for(final int count : counts){if(count<min || count>max){return false;}}
		return true;
	}

	/** Test expectations stay fatal with assertions disabled in the parent JVM. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Pins wrapper subprocesses to this JVM's Java and a deterministic tiny runtime. */
	private static ProcessBuilder makeChild(final List<String> command, final Path log){
		final ProcessBuilder pb=new ProcessBuilder(command).redirectErrorStream(true)
				.redirectOutput(log.toFile());
		final String javaHome=System.getProperty("java.home");
		final String oldPath=pb.environment().get("PATH");
		pb.environment().put("JAVA_HOME", javaHome);
		pb.environment().put("PATH", Paths.get(javaHome, "bin")+File.pathSeparator+
				(oldPath==null ? "" : oldPath));
		return pb;
	}

	private static final class Result{
		Result(final byte[] bytes_, final int[] counts_){
			bytes=bytes_;
			counts=counts_;
		}
		final byte[] bytes;
		final int[] counts;
	}
}
