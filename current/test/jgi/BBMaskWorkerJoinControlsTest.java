package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import parse.Parse;

/** Ordinary two-worker CLI controls for BBMask's repeat, entropy and SAM joins.
 * @author Sucrose
 */
public final class BBMaskWorkerJoinControlsTest{

	/** Runs the same three literal controls before and after an interruption repair. */
	public static void main(final String[] args) throws Exception{
		String directory=null, root=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("bbtools")){root=pair[1];}
			else{throw new IllegalArgumentException("Unknown control option: "+arg);}
		}
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		require(root!=null && !root.isEmpty(), "Set bbtools=<repo root>.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmaskworkerjoincontrolstest.sh with assertions enabled.");
		final BBMaskWorkerJoinControlsTest test=new BBMaskWorkerJoinControlsTest(
			Paths.get(directory).toAbsolutePath().normalize(), Paths.get(root).toAbsolutePath().normalize());
		for(final String family : FAMILIES){test.run(family);}
		test.check(test.cases==3, "Expected three caller families: "+test.cases);
		test.check(test.children==3, "Expected three fresh bbmask.sh children: "+test.children);
		final String summary="BBMASK_WORKER_JOIN_CONTROLS_"+(test.failures==0 ? "PASS" : "FAIL")+
			" cases="+test.cases+" children="+test.children+" checks="+test.checks+
			" failures="+test.failures+"\n";
		Files.write(test.out.resolve("summary.txt"), summary.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+test.out);}
		System.out.print(summary);
	}

	/** Binds a new evidence directory and the real public launcher. */
	private BBMaskWorkerJoinControlsTest(final Path out_, final Path repo_) throws Exception{
		out=out_; repo=repo_;
		require(Files.isRegularFile(repo.resolve("bbmask.sh")), "Missing bbmask.sh");
		Files.createDirectory(out);
		Files.write(out.resolve("cases.tsv"), ("family\texit\toutputBytes\tmasked\tbases"+
			"\tstdoutEmpty\tstderrClean\tinputsUnchanged\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Checks exact two-record output and counters from one uninterrupted public child. */
	private void run(final String family) throws Exception{
		cases++; children++;
		final Path work=out.resolve(family);
		Files.createDirectory(work);
		final Path reference=work.resolve("ref.fa"), output=work.resolve("out.fa");
		final Path sam=family.equals("sam") ? work.resolve("reads.sam") : null;
		Files.write(reference, REFERENCE, StandardOpenOption.CREATE_NEW);
		if(sam!=null){Files.write(sam, SAM, StandardOpenOption.CREATE_NEW);}
		Files.write(work.resolve("expected.fa"), EXPECTED, StandardOpenOption.CREATE_NEW);
		check(!Files.exists(output), family+": output exists before child");
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(repo.resolve("bbmask.sh").toString());
		command.addAll(Arrays.asList("-Xmx128m", "-ea", "t=2", "simd=f", "split=f", "minscaf=1",
			"ow=t", "app=f", "kr=1", "ke=2", "window=4", "entropy=0.5", "mincount=2", "minlen=4",
			"entropymode=t", "mr="+(family.equals("repeat") ? "t" : "f"),
			"me="+(family.equals("entropy") ? "t" : "f"), "in="+reference, "out="+output));
		if(sam!=null){
			command.addAll(Arrays.asList("sam="+sam, "sampad=0", "mincov=-1", "maxcov=-1"));
		}
		final StringBuilder argv=new StringBuilder();
		for(final String arg : command){argv.append('\'').append(arg.replace("'", "'\\''")).append("' ");}
		Files.write(work.resolve("argv.txt"), argv.append('\n').toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command).directory(repo.toFile())
			.redirectOutput(work.resolve("stdout.txt").toFile()).redirectError(work.resolve("stderr.txt").toFile());
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=2 -Dfile.encoding=UTF-8");
		final int exit=builder.start().waitFor();
		Files.write(work.resolve("exit.txt"), (exit+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final byte[] actual=Files.isRegularFile(output) ? Files.readAllBytes(output) : null;
		final String stderr=new String(Files.readAllBytes(work.resolve("stderr.txt")), StandardCharsets.UTF_8);
		final long[] counts=counts(stderr);
		final boolean stdoutEmpty=Files.size(work.resolve("stdout.txt"))==0;
		final boolean clean=!stderr.contains("Exception") && !stderr.contains("java.lang.") &&
			!stderr.contains("ERROR") && !stderr.contains("Error") && !stderr.contains("Unknown parameter");
		final boolean unchanged=Files.isRegularFile(reference) && Arrays.equals(REFERENCE, Files.readAllBytes(reference)) &&
			(sam==null || (Files.isRegularFile(sam) && Arrays.equals(SAM, Files.readAllBytes(sam))));
		Files.write(out.resolve("cases.tsv"), (family+"\t"+exit+"\t"+(actual==null ? "absent" : actual.length)+
			"\t"+counts[0]+"\t"+counts[1]+"\t"+stdoutEmpty+"\t"+clean+"\t"+unchanged+"\n")
			.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
		check(exit==0, family+": child exit="+exit);
		check(stdoutEmpty, family+": child wrote stdout");
		check(clean, family+": unexpected diagnostic: "+stderr);
		check(actual!=null && Arrays.equals(EXPECTED, actual), family+": exact two-record FASTA differs");
		check(counts[0]==16 && counts[1]==16, family+": expected unique complete masked16/16 report: "+stderr);
		check(unchanged, family+": input bytes changed");
		check(command.contains("-ea") && command.contains("t=2") && command.contains("simd=f") &&
			command.contains("split=f") && command.contains("minscaf=1") && command.contains("kr=1") &&
			command.contains("ke=2") && command.contains("window=4") && command.contains("mincount=2") &&
			command.contains("minlen=4") && command.contains("entropy=0.5") && command.contains("entropymode=t") &&
			command.contains("mr="+(family.equals("repeat") ? "t" : "f")) &&
			command.contains("me="+(family.equals("entropy") ? "t" : "f")) &&
			(sam==null || (command.contains("sampad=0") && command.contains("mincov=-1") && command.contains("maxcov=-1"))),
			family+": lost frozen caller settings");
	}

	/** Returns counts only for one complete 100-percent line and one raw report marker. */
	private static long[] counts(final String stderr){
		final Matcher raw=RAW_COUNTER.matcher(stderr);
		int markers=0;
		while(raw.find()){markers++;}
		if(markers!=1){return new long[]{-1, -1};}
		final Matcher matcher=COUNTER.matcher(stderr);
		if(!matcher.find()){return new long[]{-1, -1};}
		final long masked=Long.parseLong(matcher.group(1)), bases=Long.parseLong(matcher.group(2));
		if(matcher.find()){return new long[]{-1, -1};}
		return new long[]{masked, bases};
	}

	/** Accumulates evidence failures while allowing the remaining controls to finish. */
	private void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}
	}

	/** Rejects missing setup before any child starts. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private final Path out, repo;
	private int cases, children, checks, failures;
	private static final String[] FAMILIES={"repeat", "entropy", "sam"};
	private static final byte[] REFERENCE=(">ref1\nAAAAAAAA\n>ref2\nAAAAAAAA\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] EXPECTED=(">ref1\nNNNNNNNN\n>ref2\nNNNNNNNN\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] SAM=("@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref1\tLN:8\n@SQ\tSN:ref2\tLN:8\n"+
		"one\t0\tref1\t1\t60\t8=\t*\t0\t0\tAAAAAAAA\tIIIIIIII\n"+
		"two\t0\tref2\t1\t60\t8=\t*\t0\t0\tAAAAAAAA\tIIIIIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final Pattern COUNTER=Pattern.compile("(?m)^Total Bases Masked:[ \\t]+([0-9]+)/([0-9]+)[ \\t]+100\\.000%$");
	private static final Pattern RAW_COUNTER=Pattern.compile("(?m)^Total Bases Masked:");
}
