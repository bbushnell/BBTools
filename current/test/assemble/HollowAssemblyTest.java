package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import assemble.Contig;
import assemble.Tadpole2;
import ukmer.Kmer;

/** Tests Hollow auxiliary assembly routing and full-sequence endpoint identity.
 * Uses existing read fixtures, never generates a new read corpus.
 * @author Fischl */
public final class HollowAssemblyTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Takes the installation path and retains small child logs for failed cases. */
	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Specify the BBTools installation.");
		final HollowAssemblyTest test=new HollowAssemblyTest(Paths.get(args[0]).toAbsolutePath());
		test.endpoints();
		test.routing();
		System.out.println("HOLLOW_ASSEMBLY_TEST_PASS checks="+checks+" dir="+test.root);
	}

	/** All outputs are isolated from the existing read fixture and the installation. */
	private HollowAssemblyTest(final Path code_) throws Exception{
		code=code_;
		root=Files.createTempDirectory("hollow-assembly-");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Tests             ----------------*/
	/*--------------------------------------------------------------*/

	/** Arm matches, including palindromic arms, never substitute for full endpoint sequence. */
	private void endpoints(){
		Kmer.PACKED=true;
		for(int k : new int[]{65, 96, 124, 128}){
			final byte[] bases=SEQUENCE.getBytes(StandardCharsets.US_ASCII);
			final Contig contig=new Contig(bases, 0);
			final Kmer incoming=new Kmer(k, true), scratch=new Kmer(k, true);
			load(incoming, bases, 0);
			check(Tadpole2.incomingHollowEnd(contig, incoming, scratch)==0, "Exact left endpoint rejected at K"+k);
			final long hash=incoming.xor(), hash2=incoming.xor2();
			final byte[] changed=Arrays.copyOf(bases, k);
			changed[32]=(byte)(changed[32]=='A' ? 'C' : 'A');
			incoming.clear();
			for(byte b : changed){incoming.addRight(b);}
			check(hash==incoming.xor() && hash2==incoming.xor2(), "Fixture did not preserve both Hollow hashes.");
			check(Tadpole2.incomingHollowEnd(contig, incoming, scratch)<0, "Masked interior alias accepted at K"+k);
			load(incoming, bases, bases.length-k);
			incoming.rcomp();
			check(Tadpole2.incomingHollowEnd(contig, incoming, scratch)==2, "Reverse incoming right endpoint rejected.");
			incoming.rcomp();
			check(Tadpole2.incomingHollowEnd(contig, incoming, scratch)<0, "Outward endpoint accepted as an incoming join.");
			final byte[] tied=Arrays.copyOf(bases, k);
			Arrays.fill(tied, 0, 32, (byte)'A');
			Arrays.fill(tied, k-32, k, (byte)'T');
			final Contig shortContig=new Contig(tied, 1);
			load(incoming, tied, 0);
			check(incoming.corePalindrome(), "Fixture should have reverse-complement-identical arms.");
			check(Tadpole2.incomingHollowEnd(shortContig, incoming, scratch)==0, "Tied arms lost forward interior.");
			incoming.rcomp();
			check(Tadpole2.incomingHollowEnd(shortContig, incoming, scratch)==2, "Tied arms lost reverse interior.");
		}
	}

	/** Loads fixture sequence directly; the test does not use Contig's package-private loaders. */
	private static void load(final Kmer key, final byte[] bases, final int start){
		assert(start>=0 && start+key.kbig<=bases.length) : "Fixture endpoint must contain a full word.";
		key.clear();
		for(int i=start; i<start+key.kbig; i++){key.addRight(bases[i]);}
	}

	/** Exercises single-K, ordered, legacy, fusion evidence, washing, and final graph routing. */
	private void routing() throws Exception{
		final Path reads=code.resolve("testdata/crossk_long_bridge_x.fa");
		check(Files.isRegularFile(reads), "Missing existing read fixture.");
		run("single", reads, true, "k=96", "hollow=t");
		run("ordered", reads, true, "k=64,96,32", "hollow=t", "wash=t", "t=4");
		run("fusion", reads, true, "k=124,96,32", "hollow=t", "fusepath=t", "wash=t");
		run("legacy", reads, true, "k=124,96,32", "korder=legacy", "hollow=t", "wash=t");
		run("neural", reads, true, "k=124,96,32", "assemblydefaults=t", "hollow=t", "graphcover=t");
		run("short", reads, true, "k=64,32", "hollow=t");
		run("selected", reads, true, "k=64,124,96,32", "hollow=t", "hollowk=96");
		final String selected=new String(Files.readAllBytes(root.resolve("selected.log")), StandardCharsets.UTF_8);
		check(selected.contains("at k=96 (Hollow") && !selected.contains("at k=124 (Hollow"),
				"Only the selected K96 bridge may load Hollow counts.");
		run("selectedinitial", reads, true, "k=124,96,32", "hollow=t", "hollowk=96", "fusepath=t");
		final String initial=new String(Files.readAllBytes(root.resolve("selectedinitial.log")), StandardCharsets.UTF_8);
		check(!initial.contains("Loading Hollow initial graph"), "Unselected initial graph must remain solid.");
		run("deferredgraph", reads, true, "k=124,96,32", "hollow=t", "hollowk=124",
				"assemblydefaults=t", "pop=t");
		final String deferred=new String(Files.readAllBytes(root.resolve("deferredgraph.log")), StandardCharsets.UTF_8);
		final int handoff=deferred.indexOf("Loading Hollow initial graph");
		check(handoff>=0 && deferred.indexOf("Making contig graph")>handoff,
				"Constructor pop defaults must not run a solid graph before the Hollow handoff.");
		run("none", reads, true, "k=64,96,32", "hollow=t", "hollowk=none", "graphcover=t");
		run("matched", reads, true, "k=64,96,32", "hollow=f", "graphcover=t");
		check(Arrays.equals(Files.readAllBytes(root.resolve("none.fa")),
				Files.readAllBytes(root.resolve("matched.fa"))), "Empty Hollow selection changed matched solid output.");
		run("absent", reads, false, "k=64,96,32", "hollow=t", "hollowk=128");
		run("disabled", reads, false, "k=96", "hollow=f", "hollowk=96");
		run("disabledall", reads, false, "k=96", "hollow=f", "hollowk=all");
		run("implicitall", reads, false, "k=96", "hollowk=all");
		run("missing", reads, false, "k=96", "hollow=t", "hollowk");
		run("solid", reads, true, "k=64,96,32", "hollow=f");
		run("default", reads, true, "k=64,96,32");
		check(Arrays.equals(Files.readAllBytes(root.resolve("solid.fa")),
				Files.readAllBytes(root.resolve("default.fa"))), "Default solid assembly changed.");
		run("fixed", reads, false, "k=64,96,32", "hollow=t", "hashkmers=fixed");
		run("rescue", reads, false, "k=64,96,32", "hollow=t", "hollowrescue=t");
		run("graph", reads, true, "mode=graph", "contigs="+root.resolve("single.fa"), "k=96",
				"hollow=t", "out=null", "gfa="+root.resolve("graph.gfa"));
		check(Files.size(root.resolve("graph.gfa"))>0, "Graph mode did not emit its requested GFA.");
	}

	/** Runs public launchers only; process exit and output/log assertions are separate checks. */
	private void run(final String name, final Path reads, final boolean success,
			final String... options) throws Exception{
		assert(name!=null && reads!=null) : "Each test must have an isolated output and an existing input.";
		final Path output=root.resolve(name+".fa"), log=root.resolve(name+".log");
		final List<String> command=new ArrayList<String>(Arrays.asList("bash", code.resolve("tadpole.sh").toString(),
				"in="+reads, "out="+output, "assemblydefaults=f", "mode=contig", "mcs=1", "mce=1",
				"mincontig=1", "sweeplen=0", "prefilter=f", "prealloc=f", "initialsize=1000",
				"pop=f", "showstats=f", "minprob=0", "t=1", "-Xmx512m", "-Xms32m", "-eoom"));
		command.addAll(Arrays.asList(options));
		final Process child=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile()).start();
		child.getOutputStream().close();
		final int status=child.waitFor();
		check((status==0)==success, "Unexpected status "+status+"; inspect "+log);
		if(success && !name.equals("graph")){
			check(Files.size(output)>0, "Assembly was empty: "+name);
			if(!name.equals("solid") && !name.equals("default") && !name.equals("short") && !name.equals("matched")){
				final String messages=new String(Files.readAllBytes(log), StandardCharsets.UTF_8);
				check(messages.contains("Hollow"), "No Hollow phase reported: "+name);
			}
		}
	}

	/** Fails even if someone accidentally launches the test without assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Path code, root;
	private static int checks;
	private static final String SEQUENCE="CATACTCGCCTCGGCATATGTAGCCCGGAATCATCACACTCCTCAGGGAGCTCCTACTAGCTTGGAGCTCATCGCCGTTCACTCGGACTAACTCTGCAGCAAATGACGCCCGACCACGGTCCCAACTGACATCAGCCGCTTTGATGCGAAGCCAATTTCATTCCTCCTCTATCCTAAGCATCGATCTCTCCAACATCTCGACTCTTGGAG";
}
