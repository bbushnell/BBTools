package assemble;

import java.util.ArrayList;
import java.util.HashSet;

import dna.AminoAcid;

/** Exact graph regressions for topology-only simple-omnitig extraction. */
public class SimpleOmnitigExtractorUnitTest {

	public static void main(String[] args){
		int failures=0;
		failures+=run("isolatedContigUnchanged", SimpleOmnitigExtractorUnitTest::isolatedContigUnchanged);
		failures+=run("linearChainCompacts", SimpleOmnitigExtractorUnitTest::linearChainCompacts);
		failures+=run("yPrefixIsReplicated", SimpleOmnitigExtractorUnitTest::yPrefixIsReplicated);
		failures+=run("biologicalBubbleKeepsBothFullPaths", SimpleOmnitigExtractorUnitTest::biologicalBubbleKeepsBothFullPaths);
		failures+=run("xEmitsContextsWithoutPairing", SimpleOmnitigExtractorUnitTest::xEmitsContextsWithoutPairing);
		failures+=run("reverseOrientedJoin", SimpleOmnitigExtractorUnitTest::reverseOrientedJoin);
		failures+=run("longLeftEdgeConvention", SimpleOmnitigExtractorUnitTest::longLeftEdgeConvention);
		failures+=run("inexactHairpinIsSkipped", SimpleOmnitigExtractorUnitTest::inexactHairpinIsSkipped);
		failures+=run("nonreciprocalArcIsSkipped", SimpleOmnitigExtractorUnitTest::nonreciprocalArcIsSkipped);
		failures+=run("branchCycleStopsBeforeRepeat", SimpleOmnitigExtractorUnitTest::branchCycleStopsBeforeRepeat);
		failures+=run("pureCycleIsPreserved", SimpleOmnitigExtractorUnitTest::pureCycleIsPreserved);
		failures+=run("nonredundantLinearChain", SimpleOmnitigExtractorUnitTest::nonredundantLinearChain);
		failures+=run("nonredundantYChoosesOneArm", SimpleOmnitigExtractorUnitTest::nonredundantYChoosesOneArm);
		failures+=run("nonredundantXUsesTwoCopies", SimpleOmnitigExtractorUnitTest::nonredundantXUsesTwoCopies);
		failures+=run("nonredundantDoesNotPromoteShortDebris", SimpleOmnitigExtractorUnitTest::nonredundantDoesNotPromoteShortDebris);
		failures+=run("nonredundantKeepsAnchoredShortBridge", SimpleOmnitigExtractorUnitTest::nonredundantKeepsAnchoredShortBridge);
		System.out.println(failures==0 ? "ALL TESTS PASSED" : failures+" TEST(S) FAILED");
		if(failures>0){System.exit(1);}
	}

	private static int run(String name, Test test){
		try{
			test.run();
			System.out.println("PASS: "+name);
			return 0;
		}catch(Throwable t){
			System.err.println("FAIL: "+name);
			t.printStackTrace();
			return 1;
		}
	}

	private static void isolatedContigUnchanged(){
		final ArrayList<Contig> graph=list(contig(0, "AAAACCC"));
		final String before=signature(graph);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(graph, K);
		final ArrayList<Contig> out=extractor.extract();
		check(out.size()==1, "Expected one isolated output, got "+out.size());
		check(sequenceSet(out).contains(canonical("AAAACCC")), "Isolated sequence changed");
		check(signature(graph).equals(before), "Extractor changed the input graph");
	}

	private static void linearChainCompacts(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCATGG"), c=contig(2, "GGTTTT");
		connect(a, true, b, false, 11);
		connect(b, true, c, false, 13);
		final ArrayList<Contig> out=new SimpleOmnitigExtractor(list(a, b, c), K).extract();
		check(out.size()==1, "Linear chain produced "+out.size()+" outputs");
		check(sequenceSet(out).contains(canonical("AAAACCATGGTTTT")), "Wrong linear-chain sequence: "+sequenceSet(out));
	}

	private static void yPrefixIsReplicated(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCATGG"), c=contig(2, "CCTTTG");
		connect(a, true, b, false, 11);
		connect(a, true, c, false, 13);
		final ArrayList<Contig> graph=list(a, b, c);
		final String before=signature(graph);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(graph, K);
		final ArrayList<Contig> out=extractor.extract();
		check(out.size()==2, "Y produced "+out.size()+" outputs");
		final HashSet<String> sequences=sequenceSet(out);
		check(sequences.contains(canonical("AAAACCATGG")), "Missing first Y path: "+sequences);
		check(sequences.contains(canonical("AAAACCTTTG")), "Missing second Y path: "+sequences);
		check(extractor.extendedWalks()==2, "Wrong extended-walk count");
		check(signature(graph).equals(before), "Y extraction changed the input graph");
	}

	private static void biologicalBubbleKeepsBothFullPaths(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCATGG"), d=contig(2, "CCTCGG"), c=contig(3, "GGTTTT");
		connect(a, true, b, false, 11);
		connect(a, true, d, false, 13);
		connect(b, true, c, false, 17);
		connect(d, true, c, false, 19);
		final ArrayList<Contig> out=new SimpleOmnitigExtractor(list(a, b, d, c), K).extract();
		check(out.size()==2, "Bubble produced "+out.size()+" outputs");
		final HashSet<String> sequences=sequenceSet(out);
		check(sequences.contains(canonical("AAAACCATGGTTTT")), "Missing first bubble path: "+sequences);
		check(sequences.contains(canonical("AAAACCTCGGTTTT")), "Missing second bubble path: "+sequences);
	}

	private static void xEmitsContextsWithoutPairing(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "TTTTCC"), r=contig(2, "CCACGG");
		final Contig c=contig(3, "GGTTTT"), d=contig(4, "GGAAAA");
		connect(a, true, r, false, 11);
		connect(b, true, r, false, 13);
		connect(r, true, c, false, 17);
		connect(r, true, d, false, 19);
		final ArrayList<Contig> out=new SimpleOmnitigExtractor(list(a, b, r, c, d), K).extract();
		check(out.size()==4, "X produced "+out.size()+" outputs instead of four contexts");
		final HashSet<String> sequences=sequenceSet(out);
		check(sequences.contains(canonical("AAAACCACGG")), "Missing A-R context: "+sequences);
		check(sequences.contains(canonical("TTTTCCACGG")), "Missing B-R context: "+sequences);
		check(sequences.contains(canonical("CCACGGTTTT")), "Missing R-C context: "+sequences);
		check(sequences.contains(canonical("CCACGGAAAA")), "Missing R-D context: "+sequences);
		for(Contig x : out){check(x.length()<14, "X extractor invented a three-node pairing: "+new String(x.bases));}
	}

	private static void reverseOrientedJoin(){
		final Contig a=contig(0, "AAAACC");
		final Contig b=contig(1, "AATTGG");
		connect(a, true, b, true, 11);
		final ArrayList<Contig> out=new SimpleOmnitigExtractor(list(a, b), K).extract();
		check(out.size()==1, "Reverse join produced "+out.size()+" outputs");
		check(sequenceSet(out).contains(canonical("AAAACCAATT")), "Wrong reverse-oriented sequence: "+sequenceSet(out));
	}

	private static void longLeftEdgeConvention(){
		final Contig a=contig(0, "AAAAGCT"), b=contig(1, "CTATTT");
		final Edge ab=new Edge(a.id, b.id, 4, 1, 11, "GCTA".getBytes());
		/* Reverse traversal spells CAGC: complement the novel G, then orient the
		 * destination's stored right kmer GCT as AGC. */
		final Edge ba=new Edge(b.id, a.id, 4, 2, 11, "GGCT".getBytes());
		a.addRightEdge(ab);
		b.addLeftEdge(ba);
		final ArrayList<Contig> out=new SimpleOmnitigExtractor(list(a, b), K).extract();
		check(out.size()==1, "Long left-edge graph produced "+out.size()+" outputs");
		check(sequenceSet(out).contains(canonical("AAAAGCTGCTATTT")),
				"Wrong long left-edge sequence: "+sequenceSet(out));
	}

	private static void pureCycleIsPreserved(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCAAAA");
		connect(a, true, b, false, 11);
		connect(b, true, a, false, 13);
		final ArrayList<Contig> graph=list(a, b);
		final String before=signature(graph);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(graph, K);
		final ArrayList<Contig> out=extractor.extract();
		check(out.size()==2, "Pure cycle was arbitrarily linearized: "+out.size());
		check(extractor.preservedCycleContigs()==2, "Cycle-preservation count is wrong");
		check(signature(graph).equals(before), "Cycle extraction changed the input graph");
	}

	private static void inexactHairpinIsSkipped(){
		final Contig a=contig(0, "AAAACCC");
		a.addRightEdge(new Edge(a.id, a.id, 2, 3, 11, "CG".getBytes()));
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a), K);
		final ArrayList<Contig> out=extractor.extract();
		check(out.size()==1, "Inexact hairpin changed the output count");
		check(sequenceSet(out).contains(canonical("AAAACCC")), "Inexact hairpin changed the contig");
		check(extractor.inexactEdges()==1, "Inexact hairpin was not reported");
	}

	private static void nonreciprocalArcIsSkipped(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCATGG");
		a.addRightEdge(new Edge(a.id, b.id, 1, 1, 11, new byte[]{'A'}));
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a, b), K);
		final ArrayList<Contig> out=extractor.extract();
		check(out.size()==2, "Nonreciprocal arc joined two contigs");
		check(extractor.nonreciprocalEdges()==1, "Nonreciprocal arc was not reported");
	}

	private static void branchCycleStopsBeforeRepeat(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCAAAA"), c=contig(2, "CCGGGG");
		connect(a, true, b, false, 11);
		connect(b, true, a, false, 13);
		connect(a, true, c, false, 17);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a, b, c), K);
		final ArrayList<Contig> out=extractor.extract();
		check(extractor.truncatedCycleWalks()>0, "Branch-attached cycle was not detected");
		for(Contig x : out){check(x.length()<=10, "Cycle closure repeated a contig: "+new String(x.bases));}
	}

	private static void nonredundantLinearChain(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCATGG"), c=contig(2, "GGTTTT");
		connect(a, true, b, false, 11);
		connect(b, true, c, false, 13);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a, b, c), K);
		final ArrayList<Contig> out=extractor.extractNonredundant();
		check(out.size()==1, "Nonredundant linear chain produced "+out.size()+" outputs");
		check(sequenceSet(out).contains(canonical("AAAACCATGGTTTT")), "Wrong nonredundant chain");
		check(extractor.selectedJoins()==2, "Linear path cover selected the wrong number of joins");
	}

	private static void nonredundantYChoosesOneArm(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCATGG"), c=contig(2, "CCTTTG");
		connect(a, true, b, false, 11);
		connect(a, true, c, false, 13);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a, b, c), K);
		final ArrayList<Contig> out=extractor.extractNonredundant();
		final HashSet<String> sequences=sequenceSet(out);
		check(out.size()==2, "Nonredundant Y produced "+out.size()+" outputs");
		check(sequences.contains(canonical("AAAACCTTTG")), "Y did not choose the better-supported arm: "+sequences);
		check(sequences.contains(canonical("CCATGG")), "Y did not retain the unselected arm: "+sequences);
	}

	private static void nonredundantXUsesTwoCopies(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "TTTTCC"), r=contig(2, "CCACGG");
		final Contig c=contig(3, "GGTTTT"), d=contig(4, "GGAAAA");
		connect(a, true, r, false, 11);
		connect(b, true, r, false, 13);
		connect(r, true, c, false, 17);
		connect(r, true, d, false, 19);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a, b, r, c, d), K);
		final ArrayList<Contig> out=extractor.extractNonredundant();
		final HashSet<String> sequences=sequenceSet(out);
		check(out.size()==4, "Copy-aware X produced "+out.size()+" outputs instead of four");
		check(sequences.contains(canonical("CCACGGTTTT")), "Missing first selected X boundary: "+sequences);
		check(sequences.contains(canonical("CCACGGAAAA")), "Missing second selected X boundary: "+sequences);
		check(sequences.contains(canonical("AAAACC")) && sequences.contains(canonical("TTTTCC")),
				"Unselected X boundary was not retained: "+sequences);
		for(Contig x : out){check(x.length()<14, "Copy-aware X invented a full pairing: "+new String(x.bases));}
		check(extractor.xCentersDuplicated()==1 && extractor.xCopiesEmitted()==2,
				"X center was not represented exactly twice");
		check(extractor.xCenterBasesDuplicated()==6 && extractor.xNetBases()==2,
				"X copy accounting is wrong");
	}

	private static void nonredundantDoesNotPromoteShortDebris(){
		final Contig a=contig(0, "AAAACC"), b=contig(1, "CCATGG"), c=contig(2, "GGTTTT");
		connect(a, true, b, false, 11);
		connect(b, true, c, false, 13);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a, b, c), K, 7);
		final ArrayList<Contig> out=extractor.extractNonredundant();
		check(out.size()==3, "Sub-threshold chain was promoted into an output product");
		check(extractor.selectedJoins()==0 && extractor.smallOnlyJoinsDiscarded()==2,
				"Sub-threshold joins were not discarded exactly");
	}

	private static void nonredundantKeepsAnchoredShortBridge(){
		final Contig a=contig(0, "AAAACCCC"), b=contig(1, "CCATGG"), c=contig(2, "GGTTTTAA");
		connect(a, true, b, false, 11);
		connect(b, true, c, false, 13);
		final SimpleOmnitigExtractor extractor=new SimpleOmnitigExtractor(list(a, b, c), K, 8);
		final ArrayList<Contig> out=extractor.extractNonredundant();
		check(out.size()==1, "Anchored path did not retain its short internal bridge");
		check(extractor.selectedJoins()==2 && extractor.smallOnlyJoinsDiscarded()==0,
				"Anchored path joins were discarded");
	}

	private static void connect(final Contig a, final boolean aRight,
			final Contig b, final boolean bRight, final int depth){
		final int orientation=(aRight ? 1 : 0)|(bRight ? 2 : 0);
		final boolean bForward=!bRight;
		final byte abTraversalBase=orientedBase(b, bForward, K-1);
		final byte abBase=aRight ? abTraversalBase : AminoAcid.baseToComplementExtended[abTraversalBase];
		final Edge ab=new Edge(a.id, b.id, 1, orientation, depth, new byte[]{abBase});

		final int reverseOrientation=(bRight ? 1 : 0)|(aRight ? 2 : 0);
		final boolean aReverseForward=!aRight;
		final byte baTraversalBase=orientedBase(a, aReverseForward, K-1);
		final byte baBase=bRight ? baTraversalBase : AminoAcid.baseToComplementExtended[baTraversalBase];
		final Edge ba=new Edge(b.id, a.id, 1, reverseOrientation, depth, new byte[]{baBase});
		if(aRight){a.addRightEdge(ab);}else{a.addLeftEdge(ab);}
		if(bRight){b.addRightEdge(ba);}else{b.addLeftEdge(ba);}
	}

	private static byte orientedBase(final Contig c, final boolean forward, final int pos){
		return forward ? c.bases[pos] : AminoAcid.baseToComplementExtended[c.bases[c.length()-1-pos]];
	}

	private static Contig contig(final int id, final String bases){
		final Contig c=new Contig(bases.getBytes(), id);
		c.coverage=20+id;
		c.minCov=10+id;
		c.maxCov=30+id;
		c.leftCode=c.rightCode=Tadpole.F_BRANCH;
		return c;
	}

	private static ArrayList<Contig> list(Contig... contigs){
		final ArrayList<Contig> list=new ArrayList<Contig>();
		for(Contig c : contigs){list.add(c);}
		return list;
	}

	private static HashSet<String> sequenceSet(final ArrayList<Contig> contigs){
		final HashSet<String> set=new HashSet<String>();
		for(Contig c : contigs){set.add(new String(c.bases));}
		return set;
	}

	private static String canonical(final String sequence){
		final byte[] bases=sequence.getBytes();
		final byte[] reverse=new byte[bases.length];
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			reverse[i]=AminoAcid.baseToComplementExtended[bases[j]];
		}
		final String rc=new String(reverse);
		return sequence.compareTo(rc)<=0 ? sequence : rc;
	}

	private static String signature(final ArrayList<Contig> graph){
		final StringBuilder sb=new StringBuilder();
		for(Contig c : graph){
			sb.append(c.id).append(':').append(new String(c.bases)).append('|');
			appendEdges(sb, c.leftEdges); sb.append('|'); appendEdges(sb, c.rightEdges); sb.append('\n');
		}
		return sb.toString();
	}

	private static void appendEdges(final StringBuilder sb, final ArrayList<Edge> edges){
		if(edges==null){sb.append("null"); return;}
		for(Edge e : edges){
			sb.append(e.origin).append('>').append(e.destination).append(':').append(e.orientation)
					.append(':').append(e.depth).append(':').append(new String(e.bases)).append(';');
		}
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final int K=3;
	private interface Test {void run();}
}
