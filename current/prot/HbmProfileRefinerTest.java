package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

/** Independent column-mapping, trace, conservation, and serialized round-trip fixtures. */
public final class HbmProfileRefinerTest {
	public static void main(String[] args) throws Exception{
		final Path out=Paths.get(HmmComparisonData.required(HmmComparisonData.options(args),"out"));
		Files.createDirectory(out);
		final double[] bg=new double[20];Arrays.fill(bg,0.05);
		final Random random=new Random(20261008);
		for(int trial=0; trial<400; trial++){
			final byte[] pivot=sequence(random,1+random.nextInt(12)),query=sequence(random,1+random.nextInt(15));
			final AAGraph old=new AAGraph(pivot,2),newer=new AAGraph(pivot,2);
			old.weightByIdentity=newer.weightByIdentity=(trial%2==0);
			final AAAlignment a=GlocalAminoLinear.align(query,old.pivot,true);
			old.add(query,a);newer.addTrace(query,a.tStart,a.match);
			new HbmBundleLoader.Loaded(new String[]{"f"},new AAGraph[]{newer}).assertStructuralMatch(0,old);
			check(Arrays.equals(old.traverse(),newer.traverseProfile().consensus()),"Old and snapshot traversal differ");
		}
		final AAGraph insertion=graph("ACD","AWCD","mImm",4);
		final AAGraph.Traversal snap=insertion.traverseProfile();
		check(text(snap.consensus()).equals("AWCD"),"Insertion was not emitted at the expected column");
		check(snap.sourceType(1)==AAGraphNode.INS && snap.sourcePosition(1)==1,"Inserted-column provenance differs");
		check(snap.count(1,enc("W")[0])==4,"Carried insertion histogram gained a scaffold observation");
		final HbmPositionModel profile=HbmPositionModel.deriveTraversal(snap,bg);
		check(profile.scoreAt(enc("W")[0],1)==Math.round(128*Math.log(19.81)/Math.log(2)),"Rebased insertion score differs from direct count oracle");
		insertion.ref[0].insEdge.count[enc("W")[0]]=99;
		final byte[] altered=snap.consensus();altered[0]=enc("W")[0];
		check(snap.count(1,enc("W")[0])==4 && text(snap.consensus()).equals("AWCD"),"Traversal leaked graph or consensus aliases");
		check(text(graph("AWCD","ACD","mDmm",4).traverseProfile().consensus()).equals("ACD"),"Deleted column persisted in emitted frame");
		check(text(graph("ACD","AED","mmm",4).traverseProfile().consensus()).equals("AED"),"Substitution fixture failed");
		final AAGraph trimmed=new AAGraph(enc("AC"),2);
		for(int i=0; i<4; i++){trimmed.addTrace(enc("AC"),2,ops("mm"));}
		final AAGraph.Traversal trim=trimmed.traverseProfile();
		check(text(trim.consensus()).equals("AC") && trim.sourcePosition(0)==2 && trim.sourcePosition(1)==3,"Trim detached counts from consensus columns");
		final AAGraph edge=new AAGraph(enc("ACD"),0);
		check(edge.addTrace(enc("WACDW"),0,ops("ImmmI"))==2,"Terminal exclusion accounting differs from pad0 policy");
		boolean rejected=false;
		try{profile.requireConsensus(enc("ACDW"));}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected,"A same-length wrong consensus must fail binding");
		rejected=false;
		try{edge.addTrace(enc("ACD"),0,ops("mmZ"));}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && edge.ref[0].countSum==2,"Invalid trace mutated the graph before failure");
		final byte[] seed=enc("ACD");
		final AAGraph initial=AAGraphScorer.buildModel(seed,new byte[][]{seed,seed,seed});
		final byte[] member=enc("WACDW");final byte[][] members={member,member,member,member};
		rejected=false;
		try{HbmProfileRefiner.refine(initial,members,bg,0);}catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("padding");}
		check(rejected,"Insufficient growth padding must fail explicitly");
		final HbmProfileRefiner.Result refined=HbmProfileRefiner.refine(initial,members,bg,2);
		check(text(refined.graph.pivot).equals("WACDW"),"Shared terminal extensions failed to enter the refined consensus");
		check(refined.members==4 && refined.residues==20 && refined.excludedTerminalResidues==0,"Refinement conservation counters differ");
		for(AAGraphNode node : refined.graph.ref){check(node.countSum==5,"Final graph must contain four members plus one new scaffold");}
		final HbmProfileRefiner.Result repeat=HbmProfileRefiner.refine(initial,members,bg,2);
		new HbmBundleLoader.Loaded(new String[]{"f"},new AAGraph[]{repeat.graph}).assertStructuralMatch(0,refined.graph);
		final HbmProfileRefiner.Result unchanged=HbmProfileRefiner.refine(initial,new byte[][]{seed,seed,seed},bg,2);
		check(Arrays.equals(unchanged.graph.pivot,seed),"Identical-member control changed consensus");
		check(initial.ref[0].countSum==4,"Refinement modified its initial HBM");
		final HbmBundleLoader.Loaded initialLoaded=new HbmBundleLoader.Loaded(new String[]{"f"},new AAGraph[]{initial});
		final double[] detached=initialLoaded.positionBackground();
		final double saved=detached[0];detached[0]=0;
		check(initialLoaded.positionBackground()[0]==saved,"Background getter exposed mutable loaded state");
		final HbmProfileRefiner.Result forwarded=initialLoaded.refineFamily(0,members,bg,2);
		new HbmBundleLoader.Loaded(new String[]{"f"},new AAGraph[]{forwarded.graph}).assertStructuralMatch(0,refined.graph);
		initialLoaded.assertStructuralMatch(0,AAGraphScorer.buildModel(seed,new byte[][]{seed,seed,seed}));
		final byte[][] provenance=new byte[HbmBundleFormat.PROVENANCE_COUNT][32];
		for(byte[] row : provenance){Arrays.fill(row,(byte)1);}
		final Path bundle=out.resolve("refined.mqhb");
		HbmBundleBuilder.build(bundle,Collections.singletonList(new HbmBundleBuilder.FamilyInput("f",refined.graph.pivot,refined.graph)),provenance);
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(bundle,Collections.singletonList("f"),id->refined.graph.pivot,provenance);
		loaded.assertStructuralMatch(0,refined.graph);
		System.err.println("PROFILE_REFINER_TEST_PASS trace_parity=400 columns=PASS edges=PASS conservation=PASS roundtrip=PASS");
	}
	private static AAGraph graph(String pivot,String member,String path,int n){
		final AAGraph graph=new AAGraph(enc(pivot),0);
		for(int i=0; i<n; i++){graph.addTrace(enc(member),0,ops(path));}
		return graph;
	}
	private static byte[] sequence(Random random,int n){final byte[] a=new byte[n];for(int i=0; i<n; i++){a[i]=(byte)random.nextInt(20);}return a;}
	private static byte[] enc(String s){return Blosum62.encode(s.getBytes(StandardCharsets.US_ASCII),s);}
	private static byte[] ops(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static String text(byte[] a){return AAGraph.decode(a);}
	private static void check(boolean value,String reason){if(!value){throw new AssertionError(reason);}}
}
