package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import dna.AminoAcid;

/** Runs isolated Tadpole extension processes with exactly known reference continuations.
 * @author Brian Bushnell, Fischl */
public final class HollowExtensionTest {

	/** Accepts the installation path; keeps fixtures and child logs for diagnosis. */
	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify the BBTools installation.");
		final HollowExtensionTest test=new HollowExtensionTest(Paths.get(args[0]).toAbsolutePath());
		test.run();
		System.out.println("HOLLOW_EXTENSION_TEST_PASS checks="+checks+" dir="+test.root);
	}

	/** Isolated paths prevent one case from consuming another's outputs. */
	private HollowExtensionTest(final Path code_) throws Exception{
		code=code_;
		root=Files.createTempDirectory("hollow-extension-");
	}

	/** Tests middle-error rescue, default-off parity, true branches and invalid modes. */
	private void run() throws Exception{
		for(int span : new int[]{96, 124, 128}){
			final String reference=randomBases(span+40, 42000+span);
			final String original=reference.substring(0, span);
			final int middle=span/2;
			final char alternate=original.charAt(middle)=='A' ? 'C' : 'A';
			final String query=original.substring(0, middle)+alternate+original.substring(middle+1);
			final Path reads=write("reads"+span+".fa", ">f\n"+query+"\n>r\n"+reverse(query)+"\n");
			final Path evidence=write("extra"+span+".fa", ">a\n"+reference+"\n>b\n"+reference+"\n>c\n"+reference+"\n");
			final List<String> off=runCase("off"+span, reads, evidence, span, true, "hollow=f");
			check(off.equals(runCase("default"+span, reads, evidence, span, true)), "Default and explicit solid modes differ.");
			check(sequences(off).equals(Arrays.asList(query, reverse(query))), "Solid mode must not invent support across the error.");
			final List<String> on=runCase("on"+span, reads, evidence, span, true, "hollow=t");
			final String extended=query+reference.substring(span, span+16);
			check(sequences(on).equals(Arrays.asList(extended, reverse(extended))), "Hollow mode must add only the known continuation in both orientations.");
			check(on.equals(runCase("threads"+span, reads, evidence, span, true, "hollow=t", "t=4")), "Thread count changed extension output.");
			final Path clean=write("clean"+span+".fa", ">f\n"+original+"\n>r\n"+reverse(original)+"\n");
			check(runCase("cleanSolid"+span, clean, evidence, span, true).equals(
					runCase("cleanHollow"+span, clean, evidence, span, true, "hollow=t")),
					"Unambiguous clean extension should match solid mode.");
		}
		final int span=128;
		final String ref=randomBases(span+40, 42128), seed=ref.substring(0, span);
		final String alternate=seed+(ref.charAt(span)=='A' ? 'C' : 'A')+ref.substring(span+1);
		final Path reads=write("branchReads.fa", ">f\n"+seed+"\n");
		final Path extra=write("branchExtra.fa", ">a\n"+ref+"\n>b\n"+ref+"\n>c\n"+ref+
				"\n>d\n"+alternate+"\n>e\n"+alternate+"\n>f\n"+alternate+"\n");
		check(sequences(runCase("branch", reads, extra, span, true, "hollow=t")).equals(Arrays.asList(seed)),
				"Balanced supported branches must still stop.");
		for(String invalid : new String[]{"mode=correct", "k=64", "packed=f",
				"hashkmers=explicit", "hashkmers=fixed", "prefilter=t", "wash=t", "ecc=t", "ecco=t", "merge=t"}){
			runCase("invalid"+(++invalidCases), reads, extra, span, false, "hollow=t", invalid);
		}
	}

	/** Executes only normal launchers, with small resources and deterministic output order. */
	private List<String> runCase(final String name, final Path reads, final Path extra,
			final int span, final boolean success, final String... options) throws Exception{
		assert(reads!=null && extra!=null) : "Every comparison uses the same named evidence files.";
		final Path output=root.resolve(name+".fa"), log=root.resolve(name+".log");
		final List<String> command=new ArrayList<String>(Arrays.asList("bash", code.resolve("tadpole.sh").toString(),
				"in="+reads, "extra="+extra, "out="+output, "k="+span, "mode=extend",
				"el=16", "er=16", "mcs=2", "mce=2", "minprob=0", "ordered=t",
				"prealloc=f", "prefilter=f", "initialsize=1000", "t=1", "-Xmx512m", "-Xms32m", "-eoom"));
		command.addAll(Arrays.asList(options));
		final Process child=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile()).start();
		child.getOutputStream().close();
		final int status=child.waitFor();
		check((status==0)==success, "Unexpected status "+status+"; inspect "+log);
		if(!success){return null;}
		return Files.readAllLines(output, StandardCharsets.US_ASCII);
	}

	/** Reads only small fixture FASTAs; preserves sequence order and ignores headers. */
	private static List<String> sequences(final List<String> lines){
		final List<String> result=new ArrayList<String>();
		StringBuilder sequence=null;
		for(String line : lines){
			if(line.startsWith(">")){
				if(sequence!=null){result.add(sequence.toString());}
				sequence=new StringBuilder();
			}else{sequence.append(line);}
		}
		if(sequence!=null){result.add(sequence.toString());}
		return result;
	}

	/** Creates a fresh fixture; no production reads are modified. */
	private Path write(final String name, final String content) throws Exception{
		return Files.write(root.resolve(name), content.getBytes(StandardCharsets.US_ASCII));
	}

	/** Uses independently generated DNA, not assembler output, as continuation truth. */
	private static String randomBases(final int length, final long seed){
		final Random random=new Random(seed);
		final StringBuilder bases=new StringBuilder(length);
		for(int i=0; i<length; i++){bases.append("ACGT".charAt(random.nextInt(4)));}
		return bases.toString();
	}

	/** Converts whole fixtures to the opposite biological orientation. */
	private static String reverse(final String sequence){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		AminoAcid.reverseComplementBasesInPlace(bases);
		return new String(bases, StandardCharsets.US_ASCII);
	}

	/** Fails loudly on any scientific or execution expectation mismatch. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private final Path code, root;
	private int invalidCases;
	private static int checks;
}
