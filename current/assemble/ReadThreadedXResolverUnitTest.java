package assemble;

import java.util.ArrayList;
import java.util.HashMap;

import structures.IntList;
import structures.LongList;

/** Assertion-based regressions for evidence-gated X-repeat detection. */
public class ReadThreadedXResolverUnitTest {

	public static void main(String[] args){
		int failures=0;
		failures+=run("closedXIsCandidate", ReadThreadedXResolverUnitTest::closedXIsCandidate);
		failures+=run("supportedPairingResolves", ReadThreadedXResolverUnitTest::supportedPairingResolves);
		failures+=run("reverseReadsCount", ReadThreadedXResolverUnitTest::reverseReadsCount);
		failures+=run("oneReadCountsOnce", ReadThreadedXResolverUnitTest::oneReadCountsOnce);
		failures+=run("matePairCountsOnce", ReadThreadedXResolverUnitTest::matePairCountsOnce);
		failures+=run("oneCrossingReadDeclines", ReadThreadedXResolverUnitTest::oneCrossingReadDeclines);
		failures+=run("noiseCannotOutvotePath", ReadThreadedXResolverUnitTest::noiseCannotOutvotePath);
		failures+=run("openNeighborDeclines", ReadThreadedXResolverUnitTest::openNeighborDeclines);
		failures+=run("supportedXResolves", ReadThreadedXResolverUnitTest::supportedXResolves);
		failures+=run("reverseOrientedXResolves", ReadThreadedXResolverUnitTest::reverseOrientedXResolves);
		failures+=run("invalidSpliceRestoresOrientation", ReadThreadedXResolverUnitTest::invalidSpliceRestoresOrientation);
		failures+=run("exteriorEdgesArePreserved", ReadThreadedXResolverUnitTest::exteriorEdgesArePreserved);
		failures+=run("overlappingXDeclinesAfterNeighborMerge", ReadThreadedXResolverUnitTest::overlappingXDeclinesAfterNeighborMerge);
		failures+=run("misindexedContigsCrashLoudly", ReadThreadedXResolverUnitTest::misindexedContigsCrashLoudly);
		System.out.println(failures==0 ? "ALL TESTS PASSED" : failures+" TEST(S) FAILED");
		if(failures>0){System.exit(1);}
	}

	private static void closedXIsCandidate(){
		Fixture f=fixture();
		final String before=signature(f);
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		check(r.findCandidates()==1, "Expected exactly one closed X candidate");
		check(r.candidate(2)!=null, "Middle contig was not indexed");
		check(r.pairing(2)<0, "Topology alone authorized a pairing");
		check(r.resolveSupported()==0 && signature(f).equals(before),
				"An unsupported X changed the graph");
	}

	private static void supportedPairingResolves(){
		Fixture f=fixture();
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		IntList path=new IntList();
		path.add(0); path.add(2); path.add(3);
		r.observePath(path); r.observePath(path);
		path.clear(); path.add(1); path.add(2); path.add(4);
		r.observePath(path); r.observePath(path);
		check(r.pairing(2)==0, "Independent diagonal paths did not authorize pairing");
		check(r.support(2, 0, 0)==2 && r.support(2, 1, 1)==2, "Support matrix is wrong");
	}

	private static void reverseReadsCount(){
		Fixture f=fixture();
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		r.observe(3, 2, 0); r.observe(3, 2, 0);
		r.observe(4, 2, 1); r.observe(4, 2, 1);
		check(r.pairing(2)==0, "Reverse-oriented reads were not normalized");
	}

	private static void oneReadCountsOnce(){
		Fixture f=fixture();
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		IntList path=new IntList();
		path.add(0); path.add(2); path.add(3); path.add(0); path.add(2); path.add(3);
		r.observePath(path, new LongList(2));
		check(r.support(2, 0, 0)==1, "One read contributed duplicate support");
	}

	private static void matePairCountsOnce(){
		Fixture f=fixture();
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		IntList path=new IntList(); path.add(0); path.add(2); path.add(3);
		LongList seen=new LongList(2);
		r.observePath(path, seen);
		r.observePath(path, seen);
		check(r.support(2, 0, 0)==1, "Two mates from one fragment contributed two votes");
	}

	private static void oneCrossingReadDeclines(){
		Fixture f=fixture();
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		for(int i=0; i<4; i++){r.observe(0, 2, 3); r.observe(1, 2, 4);}
		r.observe(0, 2, 4);
		check(r.pairing(2)<0, "Conflicting traversal evidence was ignored");
	}

	private static void noiseCannotOutvotePath(){
		Fixture f=fixture();
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 5);
		r.findCandidates();
		for(int i=0; i<2; i++){r.observe(0, 2, 3); r.observe(1, 2, 4);}
		for(int i=0; i<3; i++){r.observe(0, 2, 4);}
		check(r.pairing(2)<0, "Configured noise outvoted a proposed path");
	}

	private static void openNeighborDeclines(){
		Fixture f=fixture();
		Contig external=contig(5, "TTTTTAAAAA");
		f.contigs.add(external);
		connect(f, f.contigs.get(0), true, external, false, 7);
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		check(r.findCandidates()==0, "X with a second joined-side edge was accepted");
	}

	private static void supportedXResolves(){
		Fixture f=fixture();
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		for(int i=0; i<3; i++){r.observe(0, 2, 3); r.observe(1, 2, 4);}
		check(r.resolveSupported()==1, "Supported X was not resolved");
		check(sequence(f.contigs.get(0)).equals("AAAAACCCCCCGGGGGGTTTTT"),
				"First resolved path is wrong: "+sequence(f.contigs.get(0)));
		check(sequence(f.contigs.get(1)).equals("TTTTTCCCCCCGGGGGGAAAAA"),
				"Second resolved path is wrong: "+sequence(f.contigs.get(1)));
		check(f.contigs.get(2).used() && f.contigs.get(3).used() && f.contigs.get(4).used(),
				"Absorbed X nodes were not retired");
		check(!f.destMap.containsKey(2) && !f.destMap.containsKey(3) && !f.destMap.containsKey(4),
				"Resolved X retained destination-map entries");
	}

	private static void reverseOrientedXResolves(){
		Fixture f=fixture();
		f.contigs.get(0).flip(f.destMap.get(0));
		f.contigs.get(4).flip(f.destMap.get(4));
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		for(int i=0; i<3; i++){r.observe(0, 2, 3); r.observe(1, 2, 4);}
		check(r.resolveSupported()==1, "Reverse-oriented supported X was not resolved");
		check(sequence(f.contigs.get(0)).equals("AAAAACCCCCCGGGGGGTTTTT"),
				"Reverse-oriented first path is wrong: "+sequence(f.contigs.get(0)));
		check(sequence(f.contigs.get(1)).equals("TTTTTCCCCCCGGGGGGAAAAA"),
				"Reverse-oriented second path is wrong: "+sequence(f.contigs.get(1)));
	}

	private static void invalidSpliceRestoresOrientation(){
		Fixture f=fixture();
		f.contigs.get(0).flip(f.destMap.get(0));
		f.contigs.get(0).leftEdges.get(0).bases[0]='A';
		final String before=signature(f);
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		for(int i=0; i<3; i++){r.observe(0, 2, 3); r.observe(1, 2, 4);}
		check(r.resolveSupported()==0, "Sequence-inconsistent X was resolved");
		check(signature(f).equals(before), "Declined X did not restore orientation exactly");
	}

	private static void exteriorEdgesArePreserved(){
		Fixture f=fixture();
		Contig beforeLeft=contig(5, "ACGTAATTTT");
		Contig afterRight=contig(6, "TTTTTACGTA");
		f.contigs.add(beforeLeft); f.contigs.add(afterRight);
		connect(f, f.contigs.get(0), false, beforeLeft, true, 7);
		connect(f, f.contigs.get(3), true, afterRight, false, 9);
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		r.findCandidates();
		for(int i=0; i<3; i++){r.observe(0, 2, 3); r.observe(1, 2, 4);}
		check(r.resolveSupported()==1, "Nonterminal supported X was not resolved");
		Contig product=f.contigs.get(0);
		check(product.getLeftEdge(5, -1)!=null && product.getRightEdge(6, -1)!=null,
				"Resolved product lost an exterior edge");
		check(afterRight.getLeftEdge(0, -1)!=null,
				"Right exterior reciprocal was not redirected to the product");
		BubblePopper checker=new BubblePopper(f.contigs, f.destMap, K);
		for(Contig c : f.contigs){check(checker.validate(c), "Exterior graph validation failed");}
	}

	private static void overlappingXDeclinesAfterNeighborMerge(){
		Fixture f=fixture();
		Contig mid2=contig(5, "GGGGGCCCCC");
		Contig left2=contig(6, "AAAAAGGGGG");
		Contig right2a=contig(7, "CCCCCAAAAA");
		Contig right2b=contig(8, "CCCCCTTTTT");
		f.contigs.add(mid2); f.contigs.add(left2); f.contigs.add(right2a); f.contigs.add(right2b);
		connect(f, mid2, false, f.contigs.get(3), true, 7);
		connect(f, mid2, false, left2, true, 9);
		connect(f, mid2, true, right2a, false, 11);
		connect(f, mid2, true, right2b, false, 13);
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		check(r.findCandidates()==2, "Overlapping-X fixture did not expose two candidates");
		for(int i=0; i<3; i++){
			r.observe(0, 2, 3); r.observe(1, 2, 4);
			r.observe(3, 5, 7); r.observe(6, 5, 8);
		}
		check(r.resolveSupported()==1, "Overlapping X was resolved after its neighbor changed");
		BubblePopper checker=new BubblePopper(f.contigs, f.destMap, K);
		for(Contig c : f.contigs){check(checker.validate(c), "Overlapping-X graph validation failed");}
	}

	private static void misindexedContigsCrashLoudly(){
		Fixture f=fixture();
		f.contigs.get(0).id=99;
		ReadThreadedXResolver r=new ReadThreadedXResolver(f.contigs, f.destMap, K, 2, 0);
		boolean threw=false;
		try{r.findCandidates();}catch(IllegalStateException e){threw=true;}
		check(threw, "Misindexed graph did not crash before candidate lookup");
	}

	private static final int K=5;

	private static Fixture fixture(){
		Fixture f=new Fixture();
		f.contigs.add(contig(0, "AAAAACCCCC"));
		f.contigs.add(contig(1, "TTTTTCCCCC"));
		f.contigs.add(contig(2, "CCCCCGGGGG"));
		f.contigs.add(contig(3, "GGGGGTTTTT"));
		f.contigs.add(contig(4, "GGGGGAAAAA"));
		Contig mid=f.contigs.get(2);
		connect(f, mid, false, f.contigs.get(0), true, 11);
		connect(f, mid, false, f.contigs.get(1), true, 13);
		connect(f, mid, true, f.contigs.get(3), false, 17);
		connect(f, mid, true, f.contigs.get(4), false, 19);
		return f;
	}

	private static void connect(Fixture f, Contig a, boolean aRight, Contig b, boolean bRight, int depth){
		final int orientation=(aRight ? 1 : 0)|(bRight ? 2 : 0);
		final byte abBase=(aRight && !bRight ? b.bases[K-1] : (byte)'A');
		final Edge ab=new Edge(a.id, b.id, 1, orientation, depth, new byte[]{abBase});
		final int reverseOrientation=(bRight ? 1 : 0)|(aRight ? 2 : 0);
		final byte baBase=(bRight && !aRight ? a.bases[K-1] : (byte)'T');
		final Edge ba=new Edge(b.id, a.id, 1, reverseOrientation, depth+1, new byte[]{baBase});
		if(aRight){a.rightEdges=append(a.rightEdges, ab);}else{a.leftEdges=append(a.leftEdges, ab);}
		if(bRight){b.rightEdges=append(b.rightEdges, ba);}else{b.leftEdges=append(b.leftEdges, ba);}
		addDest(f.destMap, ab); addDest(f.destMap, ba);
	}

	private static ArrayList<Edge> append(ArrayList<Edge> list, Edge e){
		if(list==null){list=new ArrayList<Edge>();}
		list.add(e);
		return list;
	}

	private static void addDest(HashMap<Integer, ArrayList<Edge>> map, Edge e){
		ArrayList<Edge> list=map.get(e.destination);
		if(list==null){list=new ArrayList<Edge>(); map.put(e.destination, list);}
		list.add(e);
	}

	private static Contig contig(int id, String bases){
		Contig c=new Contig(bases.getBytes(), id);
		c.coverage=20;
		c.leftCode=c.rightCode=Tadpole.F_BRANCH;
		return c;
	}

	private static String sequence(Contig c){return new String(c.bases);}

	private static String signature(Fixture f){
		StringBuilder sb=new StringBuilder();
		for(Contig c : f.contigs){
			sb.append(c.id).append(':').append(sequence(c)).append(':').append(c.flipped()).append(':')
					.append(c.used()).append('|');
			appendEdges(sb, c.leftEdges); sb.append('|'); appendEdges(sb, c.rightEdges); sb.append('\n');
		}
		return sb.toString();
	}

	private static void appendEdges(StringBuilder sb, ArrayList<Edge> edges){
		if(edges==null){sb.append("null"); return;}
		for(Edge e : edges){
			sb.append(e.origin).append('>').append(e.destination).append(':').append(e.orientation)
					.append(':').append(e.depth).append(':').append(new String(e.bases)).append(';');
		}
	}

	private static int run(String name, Test test){
		try{test.run(); System.out.println("PASS: "+name); return 0;}
		catch(Throwable t){System.err.println("FAIL: "+name); t.printStackTrace(); return 1;}
	}

	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}
	private interface Test {void run();}
	private static final class Fixture {
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		final HashMap<Integer, ArrayList<Edge>> destMap=new HashMap<Integer, ArrayList<Edge>>();
	}
}
