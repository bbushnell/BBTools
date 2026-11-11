package test.assemble;

import assemble.BubblePopper;
import assemble.Contig;
import assemble.Edge;
import assemble.Tadpole;

import java.util.ArrayList;
import java.util.HashMap;

/** Assertion-based regression tests for direct merging and indirect error-bubble cleanup. */
public class BubblePopperUnitTest{

	private static final int K=5;
	/** Mirrors Tadpole's errorPath=1 defaults; update if those defaults are retuned. */
	private static final BubblePopper.ErrorClassifier DEFAULT_ERROR_CLASSIFIER=
			(high, low) -> low*16f<high || (low<=4 && high>=Math.max(3, low*2.6f));

	/** Runs the serial graph-mutation fixture suite.
	 * Arguments are ignored. The method rewrites BubblePopper's process-global
	 * switches and leaves several changed; launch this harness in a fresh JVM.
	 * Each case reports PASS/FAIL, failures print their thrown diagnostics, and
	 * the process exits with status 1 if any case fails. */
	public static void main(String[] args){
		System.err.println("BUBBLEPOPPER_ASSERTIONS enabled="+BubblePopperUnitTest.class.desiredAssertionStatus());
		check(!BubblePopper.popIndirect, "Indirect graph cleanup must default off");
		check(BubblePopper.unzipBubbles, "True-bubble unzipping must default on");
		BubblePopper.verbose=false;
		BubblePopper.popDirect=true;
		BubblePopper.popIndirect=true;
		BubblePopper.unzipBubbles=false;
		BubblePopper.crossKMerge=false;
		BubblePopper.debranch=false;
		BubblePopper.validateGraph=true;

		int failures=0;
		failures+=run("directMergeKmerOverlap", BubblePopperUnitTest::directMergeKmerOverlap);
		failures+=run("directMergeRefreshesCachedScalars",
				BubblePopperUnitTest::directMergeRefreshesCachedScalars);
		failures+=run("directMergeCoverageGate", BubblePopperUnitTest::directMergeCoverageGate);
		failures+=run("directMergeLongEdgeCoverageDoesNotOverflow",
				BubblePopperUnitTest::directMergeLongEdgeCoverageDoesNotOverflow);
		failures+=run("directMergeWithoutRetainedInbound", BubblePopperUnitTest::directMergeWithoutRetainedInbound);
		failures+=run("loopMergeRemovesInboundMapEntries", BubblePopperUnitTest::loopMergeRemovesInboundMapEntries);
		failures+=run("debranchKeepsUnselectedSharedInbound", BubblePopperUnitTest::debranchKeepsUnselectedSharedInbound);
		failures+=run("soapBubbleCollapses", BubblePopperUnitTest::soapBubbleCollapses);
		failures+=run("popRefreshesCachedScalars", BubblePopperUnitTest::popRefreshesCachedScalars);
		failures+=run("trueBubbleSurvives", BubblePopperUnitTest::trueBubbleSurvives);
		failures+=run("isolatedTrueBubbleUnzips", BubblePopperUnitTest::isolatedTrueBubbleUnzips);
		failures+=run("reverseOrientedTrueBubbleUnzips", BubblePopperUnitTest::reverseOrientedTrueBubbleUnzips);
		failures+=run("trueBubbleLongEdgeCoverageDoesNotOverflow",
				BubblePopperUnitTest::trueBubbleLongEdgeCoverageDoesNotOverflow);
		failures+=run("nonterminalTrueBubbleDeclinesUnzip", BubblePopperUnitTest::nonterminalTrueBubbleDeclinesUnzip);
		failures+=run("declinedOrientationIsRestored", BubblePopperUnitTest::declinedOrientationIsRestored);
		failures+=run("errorBubbleDeclinesUnzip", BubblePopperUnitTest::errorBubbleDeclinesUnzip);
		failures+=run("higherShortArmDeclinesUnzip", BubblePopperUnitTest::higherShortArmDeclinesUnzip);
		failures+=run("multiArmPrunesOnlyError", BubblePopperUnitTest::multiArmPrunesOnlyError);
		failures+=run("bothBoundariesMustAgree", BubblePopperUnitTest::bothBoundariesMustAgree);
		failures+=run("reverseOrientedSoapBubbleCollapses", BubblePopperUnitTest::reverseOrientedSoapBubbleCollapses);
		failures+=run("reverseOrientedTrueBubbleSurvives", BubblePopperUnitTest::reverseOrientedTrueBubbleSurvives);
		failures+=run("multibaseEdgesSpliceExactly", BubblePopperUnitTest::multibaseEdgesSpliceExactly);
		failures+=run("soapBubbleLongEdgeCoverageDoesNotOverflow",
				BubblePopperUnitTest::soapBubbleLongEdgeCoverageDoesNotOverflow);
		failures+=run("fourArmsPruneOnlyErrorsAndRemainStable", BubblePopperUnitTest::fourArmsPruneOnlyErrorsAndRemainStable);
		failures+=run("duplicateEntryEdgeDeclines", BubblePopperUnitTest::duplicateEntryEdgeDeclines);

		System.out.println(failures==0 ? "ALL TESTS PASSED" : failures+" TEST(S) FAILED");
		if(failures>0){System.exit(1);}
	}

	/** Runs one case and converts an assertion or unexpected throw into a failure count. */
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

	/** A direct merge must discard scalar caches filled before {@code left.bases} is rebuilt. */
	private static void directMergeRefreshesCachedScalars(){
		Contig left=contig(0, "AAAAACCCCC", 10, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig right=contig(1, "CCCCGATTTT", 20, Tadpole.DEAD_END, Tadpole.DEAD_END);
		Edge forward=edge(0, 1, 1, 15, "G");
		Edge reverse=edge(1, 0, 2, 15, null);
		left.rightEdges=list(forward);
		right.leftEdges=list(reverse);
		Fixture f=fixture(left, right, forward, reverse);
		checkScalars(left, 1/2f, 1f, 1/2f, ",gc=0.500", ",hh=1", ",caga=0.500");

		int expansions=new BubblePopper(f.contigs, f.destMap, K).expand(left);

		check(expansions==1, "Expected one direct scalar-cache merge, got "+expansions);
		check(sequence(left).equals("AAAAACCCCCGATTTT"),
				"Incorrect scalar-cache direct merge: "+sequence(left));
		// AAAAA CCCCC GATTTT: GC=6/16; HH=(AA+CC+TT)/(AA+TT+AT+CC+CG)=11/13.
		// CAGA=(1+(CA+TG-GA-TC)/(AC+AG+CA+GA+TC+TG+CT+GT))/2=(1+(0+0-1-0)/2)/2=1/4.
		checkScalars(left, 6/16f, 11/13f, 1/4f, ",gc=0.375", ",hh=0.846", ",caga=0.250");
	}

	/** Covers inclusive thresholds, unavailable depths, both endpoint orientations and cross-K isolation. */
	private static void directMergeCoverageGate(){
		check(BubblePopper.directCoverageGate==null, "Direct coverage filtering must default off");
		for(int orientation=0; orientation<4; orientation++){
			checkDirectCoverage(20, 35, 1.75f, true, orientation, false);
			checkDirectCoverage(20, Math.nextUp(35f), 1.75f, false, orientation, false);
			checkDirectCoverage(20, 80, 0, true, orientation, false);
			checkDirectCoverage(80, 80, 1.75f, true, orientation, false);
			checkDirectCoverage(0, 20, 1.75f, false, orientation, false);
			checkDirectCoverage(20, 0, 1.75f, false, orientation, false);
			checkDirectCoverage(20, 80, 1.75f, true, orientation, true);
		}
		check(BubblePopper.directCoverageGate==null, "Test leaked a direct coverage policy");
	}

	/** Executes one actual graph merge; rejected candidates must retain topology and sequence. */
	private static void checkDirectCoverage(final float leftDepth, final float rightDepth,
			final float ratio, final boolean accepted, final int orientation, final boolean crossK){
		assert(orientation>=0 && orientation<4) : "The fixture encodes independent source/destination flips in two bits.";
		final Contig left=contig(0, "AAAAACCCCC", leftDepth, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		final Contig right=contig(1, "CCCCGTTTTT", rightDepth, Tadpole.DEAD_END, Tadpole.DEAD_END);
		final Edge forward=edge(0, 1, 1, 15, "G"), reverse=edge(1, 0, 2, 15, null);
		left.rightEdges=list(forward);
		right.leftEdges=list(reverse);
		final Fixture f=fixture(left, right, forward, reverse);
		if((orientation&1)!=0){left.flip(f.destMap.get(left.id));}
		if((orientation&2)!=0){right.flip(f.destMap.get(right.id));}
		final String before=graphSignature(f);
		final BubblePopper.CoverageGate previous=BubblePopper.directCoverageGate;
		final boolean oldCrossK=BubblePopper.crossKMerge, oldUnzip=BubblePopper.unzipBubbles;
		final BubblePopper.CoverageGate gate=ratio==0 ? null : new BubblePopper.CoverageGate(ratio);
		BubblePopper.directCoverageGate=gate;
		BubblePopper.crossKMerge=crossK;
		BubblePopper.unzipBubbles=true;
		try{
			final int expanded=new BubblePopper(f.contigs, f.destMap, K).expand(left);
			check(expanded==(accepted ? 1 : 0), "Unexpected direct coverage decision at orientation "+orientation);
			check(right.used()==accepted, "Coverage decision did not preserve destination ownership");
			if(accepted){check(sequence(left).equals("AAAAACCCCCGTTTTT"), "Coverage gate changed accepted sequence");}
			else{check(graphSignature(f).equals(before), "Rejected coverage candidate mutated the graph");}
			if(gate!=null){
				check(gate.evaluations==(crossK ? 0 : 1), "Final-graph coverage policy leaked into cross-K processing");
				check(gate.rejected==(accepted ? 0 : 1), "Coverage rejection count differs from the merge decision");
			}
		}finally{
			BubblePopper.directCoverageGate=previous;
			BubblePopper.crossKMerge=oldCrossK;
			BubblePopper.unzipBubbles=oldUnzip;
		}
	}

	/** Abstract long-edge component: the edge is finite, sequence-consistent, and overflows
	 * only the same-K span-weighting product, not the merged contig length. */
	private static void directMergeLongEdgeCoverageDoesNotOverflow(){
		final int length=65537, depth=32768;
		final int expectedLength=65552, expectedDenominator=65548;
		check((long)depth*(length-1)>Integer.MAX_VALUE, "Fixture must overflow an int product");
		Contig left=contig(0, "AAAAACCCCC", 10, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig right=contig(1, "CCCCGTTTTT", 20, Tadpole.DEAD_END, Tadpole.DEAD_END);
		Edge forward=edge(0, 1, 1, depth, basesEndingWith(length, "CCCCG"));
		left.rightEdges=list(forward);
		Fixture f=fixture(left, right, forward);

		int expansions=new BubblePopper(f.contigs, f.destMap, K).expand(left);
		check(expansions==1, "Expected one long-edge direct merge, got "+expansions);
		check(left.length()==expectedLength, "Unexpected long-edge direct length: "+left.length());
		final double expected=(10d*6+20d*6+(long)depth*(length-1))/expectedDenominator;
		check(approxFloat(left.coverage, expected),
				"Long-edge direct coverage overflowed: "+left.coverage+" vs "+expected);
	}

	/** A terminal direct merge consumes the only inbound edge and must not leave a null-valued
	 * destination-map entry for the surviving contig. */
	private static void directMergeWithoutRetainedInbound(){
		Contig left=contig(0, "AAAAACCCCC", 10, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig right=contig(1, "CCCCGTTTTT", 20, Tadpole.DEAD_END, Tadpole.DEAD_END);
		Edge forward=edge(0, 1, 1, 15, "G");
		left.rightEdges=list(forward);
		Fixture f=fixture(left, right, forward);

		int expansions=new BubblePopper(f.contigs, f.destMap, K).expand(left);
		check(expansions==1, "Expected one terminal direct merge, got "+expansions);
		check(!f.destMap.containsKey(left.id), "Terminal merge left a null-valued destination-map entry");
		check(right.used(), "Terminal destination was not retired");
	}

	/** A merge that closes a loop removes external inbound edges from their sources; those
	 * consumed edges must also disappear from the destination map. */
	private static void loopMergeRemovesInboundMapEntries(){
		Contig left=contig(0, "AAAAACCCCC", 10, Tadpole.LOOP, Tadpole.F_BRANCH);
		Contig right=contig(1, "CCCCGTTTTT", 20, Tadpole.LOOP, Tadpole.DEAD_END);
		Contig external=contig(2, "TTTTTAAAAA", 8, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Edge forward=edge(0, 1, 3, 15, "G");
		Edge externalInbound=edge(2, 1, 1, 4, "G");
		left.rightEdges=list(forward);
		external.rightEdges=list(externalInbound);
		ArrayList<Contig> contigs=new ArrayList<Contig>();
		contigs.add(left); contigs.add(right); contigs.add(external);
		HashMap<Integer, ArrayList<Edge>> destMap=new HashMap<Integer, ArrayList<Edge>>();
		addDest(destMap, forward); addDest(destMap, externalInbound);

		int expansions=new BubblePopper(contigs, destMap, K).expand(left);
		check(expansions==1, "Expected one loop-closing merge, got "+expansions);
		check(left.leftCode==Tadpole.LOOP && left.rightCode==Tadpole.LOOP,
				"Merged contig was not marked as a loop");
		check(left.leftEdgeCount()==0 && left.rightEdgeCount()==0,
				"Closed loop retained outbound edges");
		check(external.rightEdgeCount()==0, "Closed loop retained an external inbound edge");
		check(!destMap.containsKey(left.id), "Closed loop retained consumed destination-map entries");
		check(right.used(), "Loop destination was not retired");
	}

	/** Debranch truncation must remove only the selected inbound edge when another source
	 * enters the same terminal side of the truncated contig. */
	private static void debranchKeepsUnselectedSharedInbound(){
		Contig source=contig(0, "AAAAACCCCC", 30, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig dead=contig(1, "CCCCCGGGGG", 1, Tadpole.B_BRANCH, Tadpole.DEAD_END);
		Contig live=contig(2, "GGGGGTTTTT", 20, Tadpole.B_BRANCH, Tadpole.F_BRANCH);
		Contig tail=contig(3, "TTTTTAAAAA", 20, Tadpole.B_BRANCH, Tadpole.DEAD_END);
		Contig external=contig(4, "AAAACCCCCA", 15, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Edge sourceToDead=edge(0, 1, 1, 1, "G");
		Edge deadToSource=edge(1, 0, 2, 2, null);
		Edge sourceToLive=edge(0, 2, 1, 20, "G");
		Edge liveToSource=edge(2, 0, 2, 21, null);
		Edge liveToTail=edge(2, 3, 1, 20, "G");
		Edge tailToLive=edge(3, 2, 2, 21, null);
		Edge externalToDead=edge(4, 1, 1, 3, "G");
		Edge deadToExternal=edge(1, 4, 2, 4, null);
		source.rightEdges=list(sourceToDead, sourceToLive);
		dead.leftEdges=list(deadToSource, deadToExternal);
		live.leftEdges=list(liveToSource);
		live.rightEdges=list(liveToTail);
		tail.leftEdges=list(tailToLive);
		external.rightEdges=list(externalToDead);
		ArrayList<Contig> contigs=new ArrayList<Contig>();
		contigs.add(source); contigs.add(dead); contigs.add(live); contigs.add(tail); contigs.add(external);
		HashMap<Integer, ArrayList<Edge>> destMap=new HashMap<Integer, ArrayList<Edge>>();
		addDest(destMap, sourceToDead); addDest(destMap, deadToSource);
		addDest(destMap, sourceToLive); addDest(destMap, liveToSource);
		addDest(destMap, liveToTail); addDest(destMap, tailToLive);
		addDest(destMap, externalToDead); addDest(destMap, deadToExternal);

		BubblePopper.debranch=true;
		try{
			new BubblePopper(contigs, destMap, K).debranch(source);
		}finally{
			BubblePopper.debranch=false;
		}

		check(source.rightEdgeCount()==1 && source.rightEdges.get(0)==sourceToLive,
				"Debranch did not truncate only the dead source arm");
		check(external.rightEdgeCount()==1 && external.rightEdges.get(0)==externalToDead,
				"Debranch removed an unselected inbound source edge");
		ArrayList<Edge> inbound=destMap.get(dead.id);
		check(inbound!=null && inbound.size()==1 && inbound.get(0)==externalToDead,
				"Debranch removed an unselected destination-map inbound edge");
		check(dead.leftEdgeCount()==1 && dead.leftEdges.get(0)==deadToExternal,
				"Debranch removed the unselected reciprocal dead-arm edge");
		ArrayList<Edge> sourceInbound=destMap.get(source.id);
		check(sourceInbound!=null && sourceInbound.size()==1 && sourceInbound.get(0)==liveToSource,
				"Debranch removed the reciprocal live-arm destination-map edge");
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

	/** A destructive pop must discard scalar caches filled before {@code left.bases} absorbs the path. */
	private static void popRefreshesCachedScalars(){
		Contig left=contig(0, "GGGGGAAAAA", 30, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high=contig(1, "AAAATGTGCC", 75, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig low=contig(2, "AAAAGGTGCC", 3, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(3, "TGCCACATGG", 50, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connect(left, high, right, 75, 75, "T", "A", edges);
		connect(left, low, right, 3, 3, "G", "A", edges);
		BubbleFixture f=bubbleFixture(left, new Contig[]{high, low}, right, edges);
		checkScalars(left, 1/2f, 1f, 0f, ",gc=0.500", ",hh=1", ",caga=0");

		int expansions=popper(f).expand(f.left);

		check(expansions==1, "Expected one scalar-cache bubble pop, got "+expansions);
		check(sequence(f.left).equals("GGGGGAAAAATGTGCCACATGG"),
				"Incorrect scalar-cache bubble sequence: "+sequence(f.left));
		// GGGGG AAAAA TGTGCC A CATGG: GC=12/22; HH=(AA+CC+GG)/(AA+AT+CC+GC+GG)=10/13.
		// CAGA=(1+(CA+TG-GA-TC)/(AC+AG+CA+GA+TC+TG+CT+GT))/2=(1+(2+3-1-0)/8)/2=3/4.
		checkScalars(f.left, 12/22f, 10/13f, 3/4f, ",gc=0.545", ",hh=0.769", ",caga=0.750");
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
		f.left.leftFusionEndpoint=true;
		f.right.rightFusionEndpoint=false;
		checkReciprocalDepthsAsymmetric(f);
		final float highCoverage=f.mids[0].coverage, lowCoverage=f.mids[1].coverage;
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;
		int expansions=popper(f).expand(f.left);
		BubblePopper.unzipBubbles=false;
		BubblePopper.popIndirect=true;

		check(expansions==1, "Expected one isolated true-bubble unzip, got "+expansions);
		check(f.left.leftFusionEndpoint && f.right.leftFusionEndpoint
				&& !f.left.rightFusionEndpoint && !f.right.rightFusionEndpoint,
				"Unzipping manufactured fusion eligibility instead of inheriting original outer ends");
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

	/** Abstract isolated bubble whose long exact arms cover pathCoverage() without relying on
	 * a cross-K read path or an error-classified destructive pop. */
	private static void trueBubbleLongEdgeCoverageDoesNotOverflow(){
		final int length=65537, highDepth=32768, lowDepth=32767;
		final int expectedLength=131094, expectedDenominator=131090;
		check((long)highDepth*(length-1)>Integer.MAX_VALUE, "Fixture must overflow an int product");
		Contig left=contig(0, "GGGGGAAAAA", 80, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high=contig(1, "AAAATTCCCC", 75, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig low=contig(2, "AAAAGGCCCC", 72, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(3, "CCCCCGGGGG", 80, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connectExact(left, high, right, highDepth, highDepth,
				basesEndingWith(length, "AAAAT"), basesEndingWith(length, "CCCCC"), edges);
		connectExact(left, low, right, lowDepth, lowDepth,
				basesEndingWith(length, "AAAAG"), basesEndingWith(length, "CCCCC"), edges);
		BubbleFixture f=bubbleFixture(left, new Contig[]{high, low}, right, edges);

		int expansions;
		final boolean oldPopIndirect=BubblePopper.popIndirect, oldUnzip=BubblePopper.unzipBubbles;
		try{
			BubblePopper.popIndirect=false;
			BubblePopper.unzipBubbles=true;
			expansions=popper(f).expand(left);
		}finally{
			BubblePopper.unzipBubbles=oldUnzip;
			BubblePopper.popIndirect=oldPopIndirect;
		}

		check(expansions==1, "Expected one long-edge true-bubble unzip, got "+expansions);
		check(f.left.length()==expectedLength, "Unexpected long-edge true-bubble length: "+f.left.length());
		final double fraction=highDepth/(double)(highDepth+lowDepth);
		final double expected=(80d*6*fraction+75d*6+80d*6*fraction
				+2L*highDepth*(length-1))/expectedDenominator;
		check(approxFloat(f.left.coverage, expected),
				"Long-edge true-bubble coverage overflowed: "+f.left.coverage+" vs "+expected);
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

	/** Abstract error bubble with a long representative arm: it reaches pop() and keeps
	 * both edge products above the int limit while all retained sequence is finite. */
	private static void soapBubbleLongEdgeCoverageDoesNotOverflow(){
		final int length=65537, depth=32768;
		final int expectedLength=131094, expectedDenominator=131090;
		check(2L*depth*(length-1)>Integer.MAX_VALUE, "Fixture must overflow int products");
		Contig left=contig(0, "GGGGGAAAAA", 80, Tadpole.DEAD_END, Tadpole.F_BRANCH);
		Contig high=contig(1, "AAAATTCCCC", 75, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig error=contig(2, "AAAAGGCCCC", 1, Tadpole.B_BRANCH, Tadpole.B_BRANCH);
		Contig right=contig(3, "CCCCCGGGGG", 80, Tadpole.F_BRANCH, Tadpole.DEAD_END);
		ArrayList<Edge> edges=new ArrayList<Edge>();
		connectExact(left, high, right, depth, depth,
				basesEndingWith(length, "AAAAT"), basesEndingWith(length, "CCCCC"), edges);
		connect(left, error, right, 1, 1, "G", "C", edges);
		BubbleFixture f=bubbleFixture(left, new Contig[]{high, error}, right, edges);

		int expansions=popper(f).expand(left);
		check(expansions==1, "Expected one long-edge bubble pop, got "+expansions);
		check(left.length()==expectedLength, "Unexpected long-edge bubble length: "+left.length());
		final double expected=(80d*6+75d*6+80d*6+2L*depth*(length-1))
				/expectedDenominator;
		check(approxFloat(left.coverage, expected),
				"Long-edge bubble coverage overflowed: "+left.coverage+" vs "+expected);
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

	/** Adds both independently measured directed representations of each physical connection. */
	private static void connect(Contig left, Contig mid, Contig right, int entryDepth, int exitDepth,
			String entryBases, String exitBases, ArrayList<Edge> allEdges){
		Edge lm=edge(left.id, mid.id, 1, entryDepth, entryBases);
		Edge ml=edge(mid.id, left.id, 2, entryDepth+1, null);
		Edge mr=edge(mid.id, right.id, 1, exitDepth, exitBases);
		Edge rm=edge(right.id, mid.id, 2, exitDepth+1, null);
		left.rightEdges=append(left.rightEdges, lm);
		mid.leftEdges=append(mid.leftEdges, ml);
		mid.rightEdges=append(mid.rightEdges, mr);
		right.leftEdges=append(right.leftEdges, rm);
		allEdges.add(lm); allEdges.add(ml); allEdges.add(mr); allEdges.add(rm);
	}

	/** As above, but keeps reciprocal measurements at the same span while omitting
	 * unused reverse inserted bases for long abstract graph components. */
	private static void connectExact(Contig left, Contig mid, Contig right, int entryDepth,
			int exitDepth, String entryBases, String exitBases, ArrayList<Edge> allEdges){
		Edge lm=edge(left.id, mid.id, 1, entryDepth, entryBases);
		Edge ml=edgeLengthOnly(mid.id, left.id, 2, entryBases.length(), entryDepth+1);
		Edge mr=edge(mid.id, right.id, 1, exitDepth, exitBases);
		Edge rm=edgeLengthOnly(right.id, mid.id, 2, exitBases.length(), exitDepth+1);
		left.rightEdges=append(left.rightEdges, lm);
		mid.leftEdges=append(mid.leftEdges, ml);
		mid.rightEdges=append(mid.rightEdges, mr);
		right.leftEdges=append(right.leftEdges, rm);
		allEdges.add(lm); allEdges.add(ml); allEdges.add(mr); allEdges.add(rm);
	}

	/** Real graph edge depths are directional measurements and need not match reciprocals. */
	private static void checkReciprocalDepthsAsymmetric(BubbleFixture f){
		for(Contig mid : f.mids){
			Edge entry=f.left.getRightEdge(mid.id, -1);
			Edge back=mid.getLeftEdge(f.left.id, -1);
			Edge exit=mid.getRightEdge(f.right.id, -1);
			Edge reverseExit=f.right.getLeftEdge(mid.id, -1);
			check(entry!=null && back!=null && exit!=null && reverseExit!=null,
					"Fixture is missing a reciprocal bubble edge");
			check(entry.depth!=back.depth && exit.depth!=reverseExit.depth,
					"Fixture did not exercise asymmetric reciprocal depths");
		}
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

	/** Creates a reverse edge with a measured span but no stored inserted sequence. */
	private static Edge edgeLengthOnly(int origin, int destination, int orientation, int length, int depth){
		return new Edge(origin, destination, length, orientation, depth, null);
	}

	/** Builds a finite artificial edge that satisfies exact-k validation at the target. */
	private static String basesEndingWith(int length, String suffix){
		check(length>=suffix.length(), "Suffix exceeds requested edge length");
		StringBuilder sb=new StringBuilder(length);
		for(int i=suffix.length(); i<length; i++){sb.append('G');}
		sb.append(suffix);
		return sb.toString();
	}

	/** Allows final float storage and float-vs-double accumulation order, not int-scale wrap. */
	private static boolean approxFloat(float actual, double expected){
		return Math.abs(actual-expected)<=4d*Math.ulp((float)expected);
	}

	/** Checks stale caches through direct accessors and the header builder TadpoleMergeTrace logs. */
	private static void checkScalars(Contig c, float gc, float hh, float caga,
			String headerGc, String headerHh, String headerCaga){
		check(approxFloat(c.gc(), gc), "GC cache was stale: "+c.name());
		check(approxFloat(c.hh(), hh), "HH cache was stale: "+c.name());
		check(approxFloat(c.caga(), caga), "CAGA cache was stale: "+c.name());
		String header=c.name();
		check(headerHasField(header, headerGc) && headerHasField(header, headerHh)
				&& headerHasField(header, headerCaga),
				"Header used stale scalar caches: "+header);
	}

	/** Returns true only when a comma-prefixed expected token is a complete header field. */
	private static boolean headerHasField(final String header, final String token){
		final int start=header.indexOf(token);
		return start>=0 && (start+token.length()==header.length()
				|| header.charAt(start+token.length())==',');
	}

	/** Appends one edge to an existing or previously absent owned-edge list. */
	private static ArrayList<Edge> append(ArrayList<Edge> edges, Edge edge){
		if(edges==null){edges=new ArrayList<Edge>();}
		edges.add(edge);
		return edges;
	}

	/** Builds an owned-edge list in the supplied fixture order. */
	private static ArrayList<Edge> list(Edge... edges){
		ArrayList<Edge> list=new ArrayList<Edge>();
		for(Edge e : edges){list.add(e);}
		return list;
	}

	/** Indexes one inbound edge under its destination contig id. */
	private static void addDest(HashMap<Integer, ArrayList<Edge>> map, Edge edge){
		ArrayList<Edge> list=map.get(edge.destination);
		if(list==null){list=new ArrayList<Edge>(); map.put(edge.destination, list);}
		list.add(edge);
	}

	/** Converts fixture sequence bytes using the JVM default charset. */
	private static String sequence(Contig c){return new String(c.bases);}

	/** Serializes contig state and owned edges after destructive pop attempts. */
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
	/** Appends one nullable edge list in a deterministic topology signature. */
	private static void appendEdgeSignature(StringBuilder sb, ArrayList<Edge> edges){
		if(edges==null){sb.append("null"); return;}
		for(Edge e : edges){
			sb.append(e.origin).append('>').append(e.destination).append(':').append(e.orientation).append(':')
					.append(e.length).append(':').append(e.depth).append(':')
					.append(e.bases==null ? "null" : new String(e.bases)).append(';');
		}
	}
	/** Compares expected coverages after production float accumulation and storage. */
	private static boolean approx(double a, double b){return Math.abs(a-b)<0.001;}

	/** Throws a named assertion failure for command-line execution. */
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}

	/** Zero-argument test body. */
	private interface Test{void run();}

	/** Mutable graph fixture plus destination-map index consumed by BubblePopper. */
	private static class Fixture{
		Fixture(ArrayList<Contig> contigs_, HashMap<Integer, ArrayList<Edge>> destMap_){
			contigs=contigs_; destMap=destMap_;
		}
		final ArrayList<Contig> contigs;
		final HashMap<Integer, ArrayList<Edge>> destMap;
	}

	/** Bubble fixture that retains the left flank, middle arms and right flank by role. */
	private static final class BubbleFixture extends Fixture{
		BubbleFixture(Contig left_, Contig[] mids_, Contig right_, ArrayList<Contig> contigs_,
				HashMap<Integer, ArrayList<Edge>> destMap_){
			super(contigs_, destMap_); left=left_; mids=mids_; right=right_;
		}
		final Contig left;
		final Contig[] mids;
		final Contig right;
	}
}
