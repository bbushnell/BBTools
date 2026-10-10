package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import structures.ByteBuilder;

/** Boundary, topology and native serialization checks for the D242 experiment.
 * @author Nilou
 */
public final class HbmRareStateProbeTest {

	/** Runs tiny native fixtures; no biological corpus or accepted model is modified. */
	public static void main(final String[] args) throws Exception{
		require(!HbmRareStateProbe.rare(1, 200), "exact0.5% must survive");
		require(HbmRareStateProbe.rare(1, 201), "strictly below0.5% must be removed");
		require(!HbmRareStateProbe.rare(1, 100, 100), "exact1% must survive");
		require(HbmRareStateProbe.rare(1, 101, 100), "strictly below1% must be removed");
		require(!HbmRareStateProbe.rare(0, 1000), "zero is not a removed state");
		require(!HbmRareStateProbe.rare(Integer.MAX_VALUE, 2L*Integer.MAX_VALUE), "wide cutoff arithmetic");
		final AAGraph graph=fixture();
		final HbmBundleBuilder.FamilyInput family=new HbmBundleBuilder.FamilyInput("fam0", graph.pivot, graph);
		HbmBundleBuilder.validateGraph(family);
		final byte[] original=HbmBundleBuilder.serializeBlock(family);
		final HbmRareStateProbe.Stats stats=HbmRareStateProbe.prune(graph);
		require(graph.ref[0].count[0]==1 && graph.ref[0].count[2]==0, "rare pivot protected, rare alternative removed");
		require(graph.ref[0].count[21]==3, "X must remain outside the real-residue pruning denominator");
		require(graph.ref[0].countSum==1002 && graph.ref[0].weightSum==1002, "all22 count sums recomputed");
		require(graph.del[0].countSum==0 && graph.del[1].countSum==0, "rare deletion occupancy removed");
		require(graph.ref[0].insEdge.countSum==10 && graph.ref[0].insEdge.insEdge==null, "rare insertion suffix removed");
		require(graph.del[0].insEdge!=null && graph.del[0].insEdge.countSum==10, "frequent insertion survives rare DEL parent");
		require(graph.ref[1].insEdge!=null && graph.ref[1].insEdge.countSum==5, "exact0.5% insertion survives");
		require(stats.insNodes==2 && stats.insMass==6, "five plus one insertion observations removed");
		require(stats.protectedPivots==1 && stats.delStates==2, "protected and removed events counted separately");
		HbmBundleBuilder.validateGraph(family);
		final byte[] candidate=HbmBundleBuilder.serializeBlock(family);
		require(candidate.length==original.length-2*HbmBundleFormat.INS_RECORD_BYTES, "only removed insertion records shrink raw size");
		final AAGraph restored=HbmRareStateProbe.reconstruct(candidate, graph.pivot);
		final ArrayList<HbmBundleBuilder.FamilyInput> families=new ArrayList<HbmBundleBuilder.FamilyInput>();
		families.add(new HbmBundleBuilder.FamilyInput("fam0", graph.pivot, restored));
		final byte[][] provenance=new byte[HbmBundleFormat.PROVENANCE_COUNT][];
		for(int i=0; i<provenance.length; i++){provenance[i]=HbmMemberIndexFormat.sha256(new byte[]{(byte)(i+1)}, 0, 1);}
		final Path root=Files.createTempDirectory("hbm-rare-state-test-");
		final Path bundle=root.resolve("candidate.mqhb");
		HbmBundleBuilder.build(bundle, families, provenance);
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(bundle, Arrays.asList("fam0"), name->graph.pivot, provenance);
		loaded.assertStructuralMatch(0, graph);
		final Path second=root.resolve("identity.mqhb");
		HbmBundleBuilder.build(second, families, provenance);
		require(Arrays.equals(HbmRareStateProbe.digest(bundle), HbmRareStateProbe.digest(second)), "native identity serialization must be exact");
		final ByteBuilder fasta=new ByteBuilder(); fasta.append(">fam0\n");
		for(byte residue:graph.pivot){fasta.append(AminoAcid.numberToAcid[residue]);} fasta.nl();
		final Path reference=root.resolve("consensus.faa"), manifest=root.resolve("provenance.tsv");
		Files.write(reference, fasta.toBytes());
		final String[] keys={"format_spec", "loader_source", "loader_class", "aagraph", "aagraphnode", "aagraphscorer",
			"glocalaminolinear", "blosum62", "builder_source", "builder_class", "roster", "consensus_ref", "source_corpus",
			"cluster_membership", "member_policy", "member_manifest"};
		final ByteBuilder provenanceText=new ByteBuilder();
		for(int i=0; i<keys.length; i++){
			final String hex=HbmMemberIndexFormat.toHexLower(provenance[i]);
			provenanceText.append(keys[i]).append("_sha80\t").append(hex.substring(hex.length()-20)).nl();
		}
		Files.write(manifest, provenanceText.toBytes());
		final Path originalFile=root.resolve("original.mqhb");
		families.set(0, new HbmBundleBuilder.FamilyInput("fam0", graph.pivot, HbmRareStateProbe.reconstruct(original, graph.pivot)));
		HbmBundleBuilder.build(originalFile, families, provenance);
		HbmRareStateProbe.run(originalFile, reference.toString(), manifest.toString(), root.resolve("identity-run"), false);
		HbmRareStateProbe.run(originalFile, reference.toString(), manifest.toString(), root.resolve("prune-run"), true);
		require(Arrays.equals(HbmRareStateProbe.digest(bundle), HbmRareStateProbe.digest(root.resolve("prune-run/candidate.mqhb"))),
			"full native file transform must match independently inspected graph");
		final AAGraph exact05=boundaryFixture(200), below05=boundaryFixture(201);
		HbmRareStateProbe.prune(exact05);
		HbmRareStateProbe.prune(below05);
		require(exact05.ref[0].count[0]==1, "default0.5% exact boundary parity");
		require(below05.ref[0].count[0]==0, "default0.5% below-boundary removal");
		final AAGraph exact01=boundaryFixture(100), below01=boundaryFixture(101);
		HbmRareStateProbe.prune(exact01, 100);
		HbmRareStateProbe.prune(below01, 100);
		require(exact01.ref[0].count[0]==1, "1% exact boundary survives");
		require(below01.ref[0].count[0]==0, "1% below-boundary removal");
		final AAGraph bad=fixture(); bad.ref[0].insEdge.insEdge.countSum=11;
		try{
			HbmRareStateProbe.prune(bad); throw new AssertionError("Nonmonotonic insertion counts were accepted");
		}catch(IllegalArgumentException expected){require(expected.getMessage().contains("prefix"), "diagnostic prefix-count failure");}
		final AAGraph changedDenominator=fixture();
		// The scorer uses position1 (1001 observations), not anchor0 (1004).
		set(changedDenominator.ref[0], 0, 2000);
		changedDenominator.ref[0].insEdge=ins(1, 10, 0);
		HbmRareStateProbe.prune(changedDenominator);
		require(changedDenominator.ref[0].insEdge!=null, "insertion denominator is current rpos, not anchor depth");
		System.out.println("HbmRareStateProbeTest PASS cutoff,pivot,X,DEL,INS suffix,position,monotonicity,native roundtrip");
		System.out.println("fixture_directory\t"+root);
	}

	private HbmRareStateProbeTest(){}

	/** Synthetic original counts make every removal decision independently inspectable. */
	private static AAGraph fixture(){
		final byte[] pivot={0, 1, 2};// Explicit native residue indices match the hand-written histogram fixture.
		final AAGraph graph=new AAGraph(pivot, 0);
		set(graph.ref[0], 0, 1); add(graph.ref[0], 1, 998); add(graph.ref[0], 2, 1); add(graph.ref[0], 21, 3);
		set(graph.ref[1], 1, 996); add(graph.ref[1], 2, 4);
		set(graph.ref[2], 2, 1000);
		graph.del[0].countSum=graph.del[0].weightSum=1;
		graph.del[1].countSum=graph.del[1].weightSum=1;
		graph.ref[0].insEdge=ins(1, 10, 0);
		graph.ref[0].insEdge.insEdge=ins(1, 5, 0);
		graph.ref[0].insEdge.insEdge.insEdge=ins(1, 1, 0);
		graph.del[0].insEdge=ins(1, 10, 0);
		graph.ref[1].insEdge=ins(2, 5, 0);
		return graph;
	}

	/** One rare alternative plus a protected pivot gives exact cutoff fixtures. */
	private static AAGraph boundaryFixture(final int denominator){
		final AAGraph graph=new AAGraph(new byte[]{1}, 0);
		set(graph.ref[0], 1, denominator-1); add(graph.ref[0], 0, 1);
		return graph;
	}

	/** Creates one insertion node with canonical integer histogram/weight equality. */
	private static AAGraphNode ins(final int position, final int count, final int residue){
		final AAGraphNode node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, position);
		set(node, residue, count); return node;
	}

	/** Replaces fixture histogram contents without calling the production observation loop. */
	private static void set(final AAGraphNode node, final int residue, final int count){
		Arrays.fill(node.count, 0); Arrays.fill(node.weight, 0); node.countSum=node.weightSum=0; add(node, residue, count);
	}

	/** Adds deterministic fixture counts and matching aggregate weights. */
	private static void add(final AAGraphNode node, final int residue, final int count){
		node.count[residue]+=count; node.weight[residue]+=count; node.countSum+=count; node.weightSum+=count;
	}

	/** Checks run even if a caller accidentally disables Java assertions. */
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
