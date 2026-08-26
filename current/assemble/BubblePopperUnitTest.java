package assemble;

import java.util.ArrayList;
import java.util.HashMap;

/** Assertion-based regression tests for direct merging and indirect error-bubble cleanup. */
public class BubblePopperUnitTest {

	private static final int K=5;
	/** Mirrors Tadpole's errorPath=1 defaults; update if those defaults are retuned. */
	private static final BubblePopper.ErrorClassifier DEFAULT_ERROR_CLASSIFIER=
			(high, low) -> low*16f<high || (low<=4 && high>=Math.max(3, low*2.6f));

	public static void main(String[] args){
		check(!BubblePopper.popIndirect, "Indirect graph cleanup must default off");
		check(!BubblePopper.unzipBubbles, "True-bubble unzipping must default off");
		BubblePopper.verbose=false;
		BubblePopper.popDirect=true;
		BubblePopper.popIndirect=true;
		BubblePopper.unzipBubbles=false;
		BubblePopper.crossKMerge=false;
		BubblePopper.debranch=false;
		BubblePopper.validateGraph=true;

		int failures=0;
		failures+=run("directMergeKmerOverlap", BubblePopperUnitTest::directMergeKmerOverlap);
		failures+=run("soapBubbleCollapses", BubblePopperUnitTest::soapBubbleCollapses);
		failures+=run("trueBubbleSurvives", BubblePopperUnitTest::trueBubbleSurvives);
		failures+=run("isolatedTrueBubbleUnzips", BubblePopperUnitTest::isolatedTrueBubbleUnzips);
		failures+=run("reverseOrientedTrueBubbleUnzips", BubblePopperUnitTest::reverseOrientedTrueBubbleUnzips);
		failures+=run("nonterminalTrueBubbleDeclinesUnzip", BubblePopperUnitTest::nonterminalTrueBubbleDeclinesUnzip);
		failures+=run("declinedOrientationIsRestored", BubblePopperUnitTest::declinedOrientationIsRestored);
		failures+=run("errorBubbleDeclinesUnzip", BubblePopperUnitTest::errorBubbleDeclinesUnzip);
		failures+=run("higherShortArmDeclinesUnzip", BubblePopperUnitTest::higherShortArmDeclinesUnzip);
		failures+=run("multiArmPrunesOnlyError", BubblePopperUnitTest::multiArmPrunesOnlyError);
		failures+=run("bothBoundariesMustAgree", BubblePopperUnitTest::bothBoundariesMustAgree);
		failures+=run("reverseOrientedSoapBubbleCollapses", BubblePopperUnitTest::reverseOrientedSoapBubbleCollapses);
		failures+=run("reverseOrientedTrueBubbleSurvives", BubblePopperUnitTest::reverseOrientedTrueBubbleSurvives);
		failures+=run("multibaseEdgesSpliceExactly", BubblePopperUnitTest::multibaseEdgesSpliceExactly);
		failures+=run("fourArmsPruneOnlyErrorsAndRemainStable", BubblePopperUnitTest::fourArmsPruneOnlyErrorsAndRemainStable);
		failures+=run("duplicateEntryEdgeDeclines", BubblePopperUnitTest::duplicateEntryEdgeDeclines);

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

	private static void directMergeKmerOverlap(){
		Contig left=contig(0, "AAAAACCCCC", 10, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig right=contig(1, "CCCCGTTTTT", 20, Tadpole.DEAD_END, Tadpole.DEAD_END);
		Edge forward=edge(0, 1, 1, 15, "G");
		Edge reverse=edge(1, 0, 2, 15, null);
		left.rightEdges=list(forward);
		right.leftEdges=list(reverse);
		Fixture f=fixture(left, right, forward, reverse);

		int expansions=new BubblePopper(f.contigs, f.destMap, K).expand(left);
		check(expansions==1, "Expected one direct merge, got "+expansions);
		check(sequence(left).equals("AAAAACCCCCGTTTTT"), "Incorrect direct merge: "+sequence(left));
		check(approx(left.coverage, 15), "Incorrect direct coverage: "+left.coverage);
		check(right.used(), "Destination was not retired");
	}

	private static void soapBubbleCollapses(){
		BubbleFixture f=twoArmBubble(75, 3, 75, 3);
		int expansions=popper(f).expand(f.left);
		check(expansions==1, "Expected one full bubble pop, got "+expansions);
		check(sequence(f.left).equals("GGGGGAAAAATTCCCCCGGGGG"), "Incorrect bubble sequence: "+sequence(f.left));
		double expected=(30d*6+75d*6+50d*6)/18d;
		check(approx(f.left.coverage, expected), "Middle coverage was not conserved: "+f.left.coverage+" vs "+expected);
		check(f.mids[0].used(), "Representative mid was not absorbed");
		check(f.mids[1].associate(), "Error alternate was not detached");
		check(f.right.used(), "Right flank was not absorbed");
	}

	private static void trueBubbleSurvives(){
		BubbleFixture f=twoArmBubble(48, 45, 48, 45);
		int expansions=popper(f).expand(f.left);
		check(expansions==0, "Comparable-depth bubble was modified");
		check(!f.left.used() && !f.right.used(), "Bubble flanks were retired");
		check(!f.mids[0].used() && !f.mids[0].associate(), "Representative arm changed");
		check(!f.mids[1].used() && !f.mids[1].associate(), "Alternate arm changed");
		check(f.left.rightEdgeCount()==2 && f.right.leftEdgeCount()==2, "Bubble edges changed");
	}

	private static void isolatedTrueBubbleUnzips(){
		BubbleFixture f=twoArmBubble(48, 45, 48, 45);
		final float highCoverage=f.mids[0].coverage, lowCoverage=f.mids[1].coverage;
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;
		int expansions=popper(f).expand(f.left);
		BubblePopper.unzipBubbles=false;
		BubblePopper.popIndirect=true;

		check(expansions==1, "Expected one isolated true-bubble unzip, got "+expansions);
		check(sequence(f.left).equals("GGGGGAAAAATTCCCCCGGGGG"),
				"Incorrect representative allele product: "+sequence(f.left));
		check(sequence(f.right).equals("GGGGGAAAAAGGCCCCCGGGGG"),
				"Incorrect alternate allele product: "+sequence(f.right));
		check(approx(f.left.coverage, (30d*6*48/93+48d*6+50d*6*48/93)/18d),
				"Incorrect representative product coverage: "+f.left.coverage);
		check(approx(f.right.coverage, (30d*6*45/93+45d*6+50d*6*45/93)/18d),
				"Incorrect alternate product coverage: "+f.right.coverage);
		check(approx(f.left.coverage*18+f.right.coverage*18, 30d*6+48d*6+45d*6+50d*6),
				"Unzipping did not conserve span-weighted coverage mass");
		check(f.mids[0].used() && f.mids[1].used(), "Original allele arms were not retired");
		check(f.mids[0].coverage==highCoverage && f.mids[1].coverage==lowCoverage,
				"Original allele-arm coverage was mutated");
		check(!f.left.used() && !f.right.used(), "A terminal allele product was retired");
		check(f.left.leftEdgeCount()==0 && f.left.rightEdgeCount()==0
				&& f.right.leftEdgeCount()==0 && f.right.rightEdgeCount()==0,
				"An isolated allele product retained graph edges");
		check(f.left.leftCode==Tadpole.DEAD_END && f.left.rightCode==Tadpole.DEAD_END
				&& f.right.leftCode==Tadpole.DEAD_END && f.right.rightCode==Tadpole.DEAD_END,
				"An unzipped allele product is not terminal");
	}

	private static void reverseOrientedTrueBubbleUnzips(){
		BubbleFixture f=twoArmBubble(48, 45, 48, 45);
		flipMidsAndRight(f);
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;
		int expansions=popper(f).expand(f.left);
		BubblePopper.unzipBubbles=false;
		BubblePopper.popIndirect=true;
		check(expansions==1, "Reverse-oriented true bubble was not unzipped");
		check(sequence(f.left).equals("GGGGGAAAAATTCCCCCGGGGG"),
				"Reverse-oriented representative product is wrong: "+sequence(f.left));
		check(sequence(f.right).equals("GGGGGAAAAAGGCCCCCGGGGG"),
				"Reverse-oriented alternate product is wrong: "+sequence(f.right));
	}

	private static void nonterminalTrueBubbleDeclinesUnzip(){
		BubbleFixture f=twoArmBubble(48, 45, 48, 45);
		addExternalLeft(f);
		final String leftBefore=sequence(f.left), rightBefore=sequence(f.right);
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;
		int expansions=popper(f).expand(f.left);
		BubblePopper.unzipBubbles=false;
		BubblePopper.popIndirect=true;
		check(expansions==0, "Nonterminal true bubble was unzipped");
		check(sequence(f.left).equals(leftBefore) && sequence(f.right).equals(rightBefore),
				"Declined nonterminal bubble sequence changed");
		check(!f.mids[0].used() && !f.mids[1].used(), "Declined nonterminal arms were retired");
		check(f.left.leftEdgeCount()==1 && f.left.rightEdgeCount()==2,
				"Declined nonterminal topology changed");
	}

	private static void declinedOrientationIsRestored(){
		BubblePopper.popDirect=false;
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;

		BubbleFixture innerFlipped=twoArmBubble(48, 45, 48, 45);
		addExternalLeft(innerFlipped);
		flipMidsAndRight(innerFlipped);
		String innerBefore=graphSignature(innerFlipped);
		int innerExpansions=popper(innerFlipped).expand(innerFlipped.left);
		check(innerExpansions==0, "Reverse-oriented nonterminal bubble was modified");
		check(graphSignature(innerFlipped).equals(innerBefore),
				"Declined candidate did not restore mid/destination orientation exactly");

		BubbleFixture centerFlipped=twoArmBubble(48, 45, 48, 45);
		addExternalLeft(centerFlipped);
		centerFlipped.left.flip(centerFlipped.destMap.get(centerFlipped.left.id));
		String centerBefore=graphSignature(centerFlipped);
		int centerExpansions=popper(centerFlipped).expand(centerFlipped.left);
		check(centerExpansions==0, "Left-facing nonterminal bubble was modified");
		check(graphSignature(centerFlipped).equals(centerBefore),
				"Declined left-side attempt did not restore center orientation exactly");

		BubblePopper.unzipBubbles=false;
		BubblePopper.popIndirect=true;
		BubblePopper.popDirect=true;
	}

	private static void errorBubbleDeclinesUnzip(){
		BubbleFixture f=twoArmBubble(75, 3, 75, 3);
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;
		int expansions=popper(f).expand(f.left);
		BubblePopper.unzipBubbles=false;
		BubblePopper.popIndirect=true;
		check(expansions==0, "Error bubble was incorrectly unzipped into two products");
		check(!f.mids[0].used() && !f.mids[1].used() && !f.mids[1].associate(),
				"Declined error bubble changed ownership");
		check(f.left.rightEdgeCount()==2 && f.right.leftEdgeCount()==2,
				"Declined error bubble topology changed");
	}

	/** The representative selector favors a sufficiently long arm; classification must still
	 * reject unzipping when the other, shorter arm is the high-depth biological path. */
	private static void higherShortArmDeclinesUnzip(){
		Contig left=contig(0, "GGGGGAAAAA", 50, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig shortHigh=contig(1, "AAAATTTT", 75, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig longLow=contig(2, "AAAAGGTTTT", 3, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(3, "TTTTTCCCCC", 50, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connect(left, shortHigh, right, 75, 75, "T", "T", edges);
		connect(left, longLow, right, 3, 3, "G", "T", edges);
		BubbleFixture f=bubbleFixture(left, new Contig[]{shortHigh, longLow}, right, edges);
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;
		int expansions=popper(f).expand(left);
		BubblePopper.unzipBubbles=false;
		BubblePopper.popIndirect=true;
		check(expansions==0, "High-depth short arm was misclassified as a true alternate");
		check(!shortHigh.used() && !longLow.used(), "Declined asymmetric arms were retired");
		check(left.rightEdgeCount()==2 && right.leftEdgeCount()==2,
				"Declined asymmetric bubble topology changed");
	}

	private static void multiArmPrunesOnlyError(){
		Contig left=contig(0, "GGGGGAAAAA", 80, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high=contig(1, "AAAATTCCCC", 75, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig error=contig(2, "AAAAGGCCCC", 3, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig real=contig(3, "AAAACACCCC", 60, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(4, "CCCCCGGGGG", 80, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connect(left, high, right, 75, 75, "T", "C", edges);
		connect(left, error, right, 3, 3, "G", "C", edges);
		connect(left, real, right, 60, 60, "C", "C", edges);
		BubbleFixture f=bubbleFixture(left, new Contig[]{high, error, real}, right, edges);

		int expansions=popper(f).expand(left);
		check(expansions==1, "Expected one selective cleanup, got "+expansions);
		check(error.associate() && error.leftEdges==null && error.rightEdges==null, "Error arm survived");
		check(!high.used() && !high.associate() && !real.used() && !real.associate(), "Real arm was retired");
		check(left.rightEdgeCount()==2 && right.leftEdgeCount()==2, "Remaining biological bubble was damaged");
		check(left.getRightEdge(error.id, -1)==null && right.getLeftEdge(error.id, -1)==null,
				"Incident error edges survived");
	}

	private static void bothBoundariesMustAgree(){
		BubbleFixture f=twoArmBubble(75, 3, 75, 60);
		int expansions=popper(f).expand(f.left);
		check(expansions==0, "One-sided error evidence incorrectly pruned an arm");
		check(!f.mids[1].used() && !f.mids[1].associate(), "Alternate arm changed");
	}

	private static void reverseOrientedSoapBubbleCollapses(){
		BubbleFixture f=twoArmBubble(75, 3, 75, 3);
		flipMidsAndRight(f);
		int expansions=popper(f).expand(f.left);
		check(expansions==1, "Reverse-oriented error bubble was not collapsed");
		check(sequence(f.left).equals("GGGGGAAAAATTCCCCCGGGGG"),
				"Reverse-oriented bubble produced the wrong sequence: "+sequence(f.left));
		check(f.mids[0].used() && f.mids[1].associate() && f.right.used(),
				"Reverse-oriented bubble retired the wrong nodes");
	}

	private static void reverseOrientedTrueBubbleSurvives(){
		BubbleFixture f=twoArmBubble(48, 45, 48, 45);
		flipMidsAndRight(f);
		int expansions=popper(f).expand(f.left);
		check(expansions==0, "Reverse-oriented biological bubble was modified");
		check(!f.left.used() && !f.right.used(), "Reverse-oriented bubble flanks were retired");
		for(Contig mid : f.mids){
			check(!mid.used() && !mid.associate(), "Reverse-oriented biological arm changed");
		}
		check(f.left.rightEdgeCount()==2, "Reverse-oriented biological bubble lost an entry edge");
	}

	private static void flipMidsAndRight(BubbleFixture fixture){
		for(Contig mid : fixture.mids){mid.flip(fixture.destMap.get(mid.id));}
		fixture.right.flip(fixture.destMap.get(fixture.right.id));
	}

	private static void multibaseEdgesSpliceExactly(){
		Contig left=contig(0, "AAAAACCCCC", 10, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high=contig(1, "CCGTATTTTT", 20, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig error=contig(2, "CCACAGTTTTT", 2, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(3, "TTTCGAAAAA", 30, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connect(left, high, right, 75, 75, "GTA", "CG", edges);
		connect(left, error, right, 2, 2, "ACA", "CG", edges);
		BubbleFixture f=bubbleFixture(left, new Contig[]{high, error}, right, edges);

		int expansions=popper(f).expand(left);
		check(expansions==1, "Expected multibase bubble pop, got "+expansions);
		check(sequence(left).equals("AAAAACCCCCGTATTTTTCGAAAAA"), "Incorrect multibase splice: "+sequence(left));
		double expected=(10d*6+20d*6+30d*6+75d*2+75d)/21d;
		check(approx(left.coverage, expected), "Incorrect multibase coverage: "+left.coverage+" vs "+expected);
	}

	/** Four arms (2 comparable survivors + 2 error arms): only the two errors detach, and a
	 * fresh popper on the reduced bubble finds it already stable (Noelle, 2026-08-26). */
	private static void fourArmsPruneOnlyErrorsAndRemainStable(){
		Contig left=contig(0, "GGGGGAAAAA", 80, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high1=contig(1, "AAAATTCCCC", 75, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig high2=contig(2, "AAAAGGCCCC", 60, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig error1=contig(3, "AAAACACCCC", 3, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig error2=contig(4, "AAAATACCCC", 2, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(5, "CCCCCGGGGG", 80, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connect(left, high1, right, 75, 75, "T", "C", edges);
		connect(left, high2, right, 60, 60, "G", "C", edges);
		connect(left, error1, right, 3, 3, "C", "C", edges);
		connect(left, error2, right, 2, 2, "A", "C", edges);
		BubbleFixture f=bubbleFixture(left, new Contig[]{high1, high2, error1, error2}, right, edges);

		int expansions=popper(f).expand(left);
		check(expansions==1, "Expected one selective cleanup pass, got "+expansions);
		check(error1.associate() && error1.leftEdges==null && error1.rightEdges==null, "Error arm 1 survived");
		check(error2.associate() && error2.leftEdges==null && error2.rightEdges==null, "Error arm 2 survived");
		check(!high1.used() && !high1.associate() && !high2.used() && !high2.associate(),
				"A comparable-depth survivor was retired");
		check(left.rightEdgeCount()==2 && right.leftEdgeCount()==2, "Surviving biological bubble was damaged");
		check(left.getRightEdge(error1.id, -1)==null && left.getRightEdge(error2.id, -1)==null
				&& right.getLeftEdge(error1.id, -1)==null && right.getLeftEdge(error2.id, -1)==null,
				"Incident error edges survived");

		//A fresh popper on the now-reduced two-arm bubble must find nothing left to prune.
		int stableExpansions=popper(f).expand(left);
		check(stableExpansions==0, "Surviving biological bubble was not stable on a second expand, got "+stableExpansions);
		check(!high1.used() && !high1.associate() && !high2.used() && !high2.associate(),
				"A second expand pass mutated the stable survivors");
		check(left.rightEdgeCount()==2 && right.leftEdgeCount()==2, "A second expand pass changed the stable bubble's edges");
	}

	/** Two edges from left to the same mid (ambiguous entry): expand must decline entirely --
	 * fetchMidNodes' duplicate-mid check aborts before any arm is classified or mutated. Built by
	 * direct list manipulation since Contig.addRightEdge would auto-merge a real duplicate at
	 * insertion time; this exercises the defensive path for a malformed graph. */
	private static void duplicateEntryEdgeDeclines(){
		Contig left=contig(0, "GGGGGAAAAA", 80, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high=contig(1, "AAAATTCCCC", 75, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig low=contig(2, "AAAAGGCCCC", 3, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(3, "CCCCCGGGGG", 80, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connect(left, high, right, 75, 75, "T", "C", edges);
		connect(left, low, right, 3, 3, "G", "C", edges);
		Edge duplicate=edge(left.id, high.id, 1, 75, "T");
		left.rightEdges=append(left.rightEdges, duplicate);
		edges.add(duplicate);
		BubbleFixture f=bubbleFixture(left, new Contig[]{high, low}, right, edges);

		int leftEdgesBefore=left.rightEdgeCount();
		int rightEdgesBefore=right.leftEdgeCount();
		int expansions=popper(f).expand(left);
		check(expansions==0, "Ambiguous duplicate entry edge was not declined, got "+expansions);
		check(!high.used() && !high.associate() && !low.used() && !low.associate(),
				"Ambiguous bubble mutated arm ownership");
		check(left.rightEdgeCount()==leftEdgesBefore && right.leftEdgeCount()==rightEdgesBefore,
				"Ambiguous bubble mutated edge counts");
	}

	private static BubbleFixture twoArmBubble(int entryHigh, int entryLow, int exitHigh, int exitLow){
		Contig left=contig(0, "GGGGGAAAAA", 30, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high=contig(1, "AAAATTCCCC", entryHigh, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig low=contig(2, "AAAAGGCCCC", entryLow, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(3, "CCCCCGGGGG", 50, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connect(left, high, right, entryHigh, exitHigh, "T", "C", edges);
		connect(left, low, right, entryLow, exitLow, "G", "C", edges);
		return bubbleFixture(left, new Contig[]{high, low}, right, edges);
	}

	private static void addExternalLeft(BubbleFixture f){
		Contig external=contig(f.contigs.size(), "TTTTTGGGG", 20, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Edge externalToLeft=edge(external.id, f.left.id, 1, 20, "G");
		Edge leftToExternal=edge(f.left.id, external.id, 2, 20, null);
		external.rightEdges=list(externalToLeft);
		f.left.leftEdges=append(f.left.leftEdges, leftToExternal);
		f.contigs.add(external);
		addDest(f.destMap, externalToLeft);
		addDest(f.destMap, leftToExternal);
	}

	/** Adds both directed representations of a physical entry and exit connection. */
	private static void connect(Contig left, Contig mid, Contig right, int entryDepth, int exitDepth,
			String entryBases, String exitBases, ArrayList<Edge> allEdges){
		Edge lm=edge(left.id, mid.id, 1, entryDepth, entryBases);
		Edge ml=edge(mid.id, left.id, 2, entryDepth, null);
		Edge mr=edge(mid.id, right.id, 1, exitDepth, exitBases);
		Edge rm=edge(right.id, mid.id, 2, exitDepth, null);
		left.rightEdges=append(left.rightEdges, lm);
		mid.leftEdges=append(mid.leftEdges, ml);
		mid.rightEdges=append(mid.rightEdges, mr);
		right.leftEdges=append(right.leftEdges, rm);
		allEdges.add(lm); allEdges.add(ml); allEdges.add(mr); allEdges.add(rm);
	}

	private static BubbleFixture bubbleFixture(Contig left, Contig[] mids, Contig right, ArrayList<Edge> edges){
		ArrayList<Contig> contigs=new ArrayList<Contig>();
		contigs.add(left);
		for(Contig c : mids){contigs.add(c);}
		contigs.add(right);
		HashMap<Integer, ArrayList<Edge>> destMap=new HashMap<Integer, ArrayList<Edge>>();
		for(Edge e : edges){addDest(destMap, e);}
		return new BubbleFixture(left, mids, right, contigs, destMap);
	}

	private static BubblePopper popper(BubbleFixture fixture){
		return new BubblePopper(fixture.contigs, fixture.destMap, K, DEFAULT_ERROR_CLASSIFIER);
	}

	private static Fixture fixture(Contig a, Contig b, Edge... edges){
		ArrayList<Contig> contigs=new ArrayList<Contig>(); contigs.add(a); contigs.add(b);
		HashMap<Integer, ArrayList<Edge>> destMap=new HashMap<Integer, ArrayList<Edge>>();
		for(Edge e : edges){addDest(destMap, e);}
		return new Fixture(contigs, destMap);
	}

	private static Contig contig(int id, String bases, float coverage, int leftCode, int rightCode){
		Contig c=new Contig(bases.getBytes(), id);
		c.coverage=coverage; c.leftCode=leftCode; c.rightCode=rightCode;
		return c;
	}

	private static Edge edge(int origin, int destination, int orientation, int depth, String bases){
		return new Edge(origin, destination, bases==null ? 1 : bases.length(), orientation, depth,
				bases==null ? null : bases.getBytes());
	}

	private static ArrayList<Edge> append(ArrayList<Edge> edges, Edge edge){
		if(edges==null){edges=new ArrayList<Edge>();}
		edges.add(edge);
		return edges;
	}

	private static ArrayList<Edge> list(Edge... edges){
		ArrayList<Edge> list=new ArrayList<Edge>();
		for(Edge e : edges){list.add(e);}
		return list;
	}

	private static void addDest(HashMap<Integer, ArrayList<Edge>> map, Edge edge){
		ArrayList<Edge> list=map.get(edge.destination);
		if(list==null){list=new ArrayList<Edge>(); map.put(edge.destination, list);}
		list.add(edge);
	}

	private static String sequence(Contig c){return new String(c.bases);}
	private static String graphSignature(Fixture f){
		StringBuilder sb=new StringBuilder();
		for(Contig c : f.contigs){
			sb.append(c.id).append(':').append(sequence(c)).append(':').append(c.leftCode).append(':')
					.append(c.rightCode).append(':').append(c.coverage).append(':').append(c.used()).append(':')
					.append(c.associate()).append('|');
			appendEdgeSignature(sb, c.leftEdges);
			sb.append('|');
			appendEdgeSignature(sb, c.rightEdges);
			sb.append('\n');
		}
		return sb.toString();
	}
	private static void appendEdgeSignature(StringBuilder sb, ArrayList<Edge> edges){
		if(edges==null){sb.append("null"); return;}
		for(Edge e : edges){
			sb.append(e.origin).append('>').append(e.destination).append(':').append(e.orientation).append(':')
					.append(e.length).append(':').append(e.depth).append(':')
					.append(e.bases==null ? "null" : new String(e.bases)).append(';');
		}
	}
	private static boolean approx(double a, double b){return Math.abs(a-b)<0.001;}
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}

	private interface Test {void run();}

	private static class Fixture {
		Fixture(ArrayList<Contig> contigs_, HashMap<Integer, ArrayList<Edge>> destMap_){
			contigs=contigs_; destMap=destMap_;
		}
		final ArrayList<Contig> contigs;
		final HashMap<Integer, ArrayList<Edge>> destMap;
	}

	private static final class BubbleFixture extends Fixture {
		BubbleFixture(Contig left_, Contig[] mids_, Contig right_, ArrayList<Contig> contigs_,
				HashMap<Integer, ArrayList<Edge>> destMap_){
			super(contigs_, destMap_); left=left_; mids=mids_; right=right_;
		}
		final Contig left;
		final Contig[] mids;
		final Contig right;
	}
}
