package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

/** Tests rejection-only use of short endpoints after ordinary fusion selection. @author Fischl */
public final class ShortFusionEvidenceTest {

	/** Keeps the established overlap tests, then exercises the optional evidence veto. */
	public static void main(final String[] args){
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Tests require -ea.");}
		CrossKTipOverlapperUnitTest.main(new String[0]);
		for(boolean graphK : new boolean[]{false, true}){
			for(int mask=0; mask<8; mask++){
				testOverlap(mask, graphK, "CCCGGGGAA", false, 0);
				testOverlap(mask, graphK, "CCCCGGGGAA", false, 0);
				testOverlap(mask, graphK, "CCGGGGAA", false, 1);
				testOverlap(mask, graphK, "AACCCGGGG", true, 0);
				testOverlap(mask, graphK, "CCCGGGG", false, 1);
			}
			testEligibility(graphK);
		}
		testRepeatedAnchors();
		testTrimRules();
		testNoNewPairs();
		testCycle();
		testDriverParsing();
		System.out.println("SHORT_FUSION_EVIDENCE_PASS checks="+checks);
	}

	/** Equal/longer alternatives veto; shorter or fully consumed evidence does not. */
	private static void testOverlap(final int mask, final boolean graphK, final String shortBases,
			final boolean shortRight, final int expected){
		assert(mask>=0 && mask<8) : "Orientation fixtures enumerate all three input inversions.";
		final ArrayList<Contig> baseline=fixture(mask, shortBases, shortRight);
		check(overlapper(baseline, false, graphK, -1, false).addEdges()==1, "Off mode must retain the original pair.");
		final ArrayList<Contig> guarded=fixture(mask, shortBases, shortRight);
		final byte[][] before=new byte[guarded.size()][];
		for(int i=0; i<before.length; i++){before[i]=guarded.get(i).bases.clone();}
		check(overlapper(guarded, true, graphK, -1, false).addEdges()==expected,
				"Unexpected evidence decision: "+mask+", "+graphK+", "+shortBases+", "+shortRight);
		check(guarded.get(2).leftEdgeCount()+guarded.get(2).rightEdgeCount()==0,
				"An evidence-only contig must never receive a join edge.");
		for(int i=0; i<before.length; i++){
			check(Arrays.equals(before[i], guarded.get(i).bases), "Evidence discovery mutated sequence bytes.");
		}
		if(expected==1){checkOnlyOriginalPair(guarded);}
	}

	/** Current topology and the optional original-dead-end requirement both apply. */
	private static void testEligibility(final boolean graphK){
		final ArrayList<Contig> ineligible=fixture(0, "CCCGGGGAA", false);
		ineligible.get(2).leftBridgeEndpoint=false;
		ineligible.get(2).leftCode=Tadpole.F_BRANCH;
		check(overlapper(ineligible, true, graphK, -1, false).addEdges()==1,
				"A topology-ineligible short end cannot become evidence under this policy.");
		final ArrayList<Contig> historical=fixture(0, "CCCGGGGAA", false);
		historical.get(2).leftFusionEndpoint=false;
		check(overlapper(historical, true, graphK, -1, true).addEdges()==1,
				"Historical dead-end restriction must apply to evidence as well as ordinary candidates.");
	}

	/** A third identical short anchor must remain visible, not erase the blocker bucket. */
	private static void testRepeatedAnchors(){
		final ArrayList<Contig> list=fixture(0, "CCCGGGGAA", false);
		list.add(contig(3, "CCCGGGGAA", true, false));
		list.add(contig(4, "CCCGGGGAA", true, false));
		check(overlapper(list, true, false, -1, false).addEdges()==0,
				"Repeated short anchors must not hide evidence behind an ambiguous-map sentinel.");
	}

	/** Evidence must respect ordinary trim admission and discarded-flank validation. */
	private static void testTrimRules(){
		for(boolean trim : new boolean[]{false, true}){
			final CrossKTipOverlapper overlap=overlapper(fixture(0, "TTCCCGGGGAA", false), true, false, -1, false);
			overlap.allowTrim=trim;
			check(overlap.addEdges()==(trim ? 0 : 1), "Short evidence bypassed the trim switch.");
		}
		check(overlapper(fixture(0, "TTCCCGGGGAA", false), true, false, 0, false).addEdges()==1,
				"Contradicted trimmed context must not veto an otherwise unique pair.");
		check(overlapper(fixture(0, "TTCCCGGGGAA", false), true, false, 2, false).addEdges()==0,
				"Allowed mismatches must preserve an eligible short alternative.");
	}

	/** Evidence cannot create a short-long, short-short, or replacement long-long pair. */
	private static void testNoNewPairs(){
		final ArrayList<Contig> list=fixture(0, "CCCGGGGAA", false);
		list.get(1).leftBridgeEndpoint=false;
		check(overlapper(list, true, false, -1, false).addEdges()==0,
				"Short evidence created a pair when ordinary reciprocal selection had none.");
		final ArrayList<Contig> shortOnly=new ArrayList<Contig>();
		shortOnly.add(contig(0, "AACCCGGGG", false, true));
		shortOnly.add(contig(1, "CCCGGGGAA", true, false));
		check(overlapper(shortOnly, true, false, -1, false).addEdges()==0,
				"Two short evidence contigs must not join each other.");
	}

	/** Removing evidence-conflicted pairs cannot open a previously rejected cycle. */
	private static void testCycle(){
		final ArrayList<Contig> list=new ArrayList<Contig>();
		list.add(contig(0, "GGATTAAC", true, true));
		list.add(contig(1, "AACGGCCT", true, true));
		list.add(contig(2, "CCTCCGGA", true, true));
		list.add(contig(3, "AACT", true, false));
		final CrossKTipOverlapper overlap=new CrossKTipOverlapper(list, 3, 3, false, 8);
		overlap.shortEvidence=true;
		check(overlap.addEdges()==0, "The frozen cycle guard must run before the evidence veto.");
		for(Contig c : list){check(c.leftEdgeCount()+c.rightEdgeCount()==0, "Cycle rejection left an edge.");}
	}

	/** The experimental spelling is consumed by the driver, not forwarded to Tadpole. */
	private static void testDriverParsing(){
		final String[] base={"in=unused.fa", "out=unused-out.fa", "k=96,124,64,32", "assemblek=96"};
		for(String flag : new String[]{"fuseshortevidence", "fuseshortevidence=t", "fuseshortevidence=f"}){
			final String[] args=Arrays.copyOf(base, base.length+1);
			args[base.length]=flag;
			new MultiKBranchRescue(args);
		}
		check(!new TadpoleMulti.Config(base).fuseShortEvidence, "Public/default configuration must remain off.");
	}

	/** Makes the same long-long pair plus one deliberately length-excluded alternative. */
	private static ArrayList<Contig> fixture(final int mask, final String shortBases, final boolean shortRight){
		assert(shortBases.length()<12) : "Evidence fixtures must stay below their explicit emission threshold.";
		final ArrayList<Contig> list=new ArrayList<Contig>();
		list.add(contig(0, "ATATATACCCCGGGG", false, true));
		list.add(contig(1, "CCCGGGGTTATATAT", true, false));
		list.add(contig(2, shortBases, !shortRight, shortRight));
		for(int i=0; i<list.size(); i++){if((mask&(1<<i))!=0){list.get(i).rcomp();}}
		return list;
	}

	/** Sets matching current/historical endpoint metadata on a nonempty DNA contig. */
	private static Contig contig(final int id, final String bases, final boolean left, final boolean right){
		assert(id>=0 && !bases.isEmpty()) : "Fixture contigs need stable indices and explicit sequence.";
		final Contig c=new Contig(bases.getBytes(StandardCharsets.US_ASCII), id);
		c.leftBridgeEndpoint=c.leftFusionEndpoint=left;
		c.rightBridgeEndpoint=c.rightFusionEndpoint=right;
		c.leftCode=left ? Tadpole.KEEP_GOING : Tadpole.F_BRANCH;
		c.rightCode=right ? Tadpole.KEEP_GOING : Tadpole.F_BRANCH;
		c.coverage=10;
		return c;
	}

	/** Only the evidence toggle differs between matched ordinary/experimental instances. */
	private static CrossKTipOverlapper overlapper(final ArrayList<Contig> list, final boolean enabled,
			final boolean graphK, final int mismatches, final boolean deadEnds){
		assert(list!=null) : "Each policy fixture supplies its own unmodified contig set.";
		final CrossKTipOverlapper result=new CrossKTipOverlapper(list, 5, 9, graphK, 12, mismatches, deadEnds);
		result.shortEvidence=enabled;
		return result;
	}

	/** A surviving result must still connect exactly the two originally retained contigs. */
	private static void checkOnlyOriginalPair(final ArrayList<Contig> list){
		assert(list.size()>=2) : "The retained pair occupies the first two fixture indices.";
		for(int i=0; i<list.size(); i++){
			final Contig c=list.get(i);
			for(ArrayList<Edge> edges : Arrays.asList(c.leftEdges, c.rightEdges)){
				if(edges!=null){for(Edge edge : edges){check(i<2 && edge.destination==1-i, "Evidence admitted a new pair.");}}
			}
		}
	}

	/** Counts a fixture requirement and fails loudly with its diagnostic reason. */
	private static void check(final boolean condition, final String message){
		checks++;
		assert(condition) : message;
	}

	private static int checks;
}
