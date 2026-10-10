package assemble;

import java.util.ArrayList;
import java.util.HashMap;

import dna.AminoAcid;

/** Assertion-based tests for exact low-depth tip-overlap joining. */
public class CrossKTipOverlapperUnitTest {

	public static void main(String[] args){
		BubblePopper.verbose=false;
		BubblePopper.popDirect=true;
		BubblePopper.popIndirect=false;
		BubblePopper.crossKMerge=true;
		BubblePopper.crossKMaxDepthRatio=3;
		BubblePopper.validateGraph=true;

		int failures=0;
		failures+=run("uniqueReciprocalOverlap", CrossKTipOverlapperUnitTest::uniqueReciprocalOverlap);
		failures+=run("inwardTipOverlap", CrossKTipOverlapperUnitTest::inwardTipOverlap);
		failures+=run("branchTipIgnored", CrossKTipOverlapperUnitTest::branchTipIgnored);
		failures+=run("equalBestIsAmbiguous", CrossKTipOverlapperUnitTest::equalBestIsAmbiguous);
		failures+=run("discardedTipDoesNotBlockJoin", CrossKTipOverlapperUnitTest::discardedTipDoesNotBlockJoin);
		failures+=run("repeatedAnchorOnOneTipIsAmbiguous", CrossKTipOverlapperUnitTest::repeatedAnchorOnOneTipIsAmbiguous);
		failures+=run("reverseOrientedOverlap", CrossKTipOverlapperUnitTest::reverseOrientedOverlap);
		failures+=run("cyclicComponentRejected", CrossKTipOverlapperUnitTest::cyclicComponentRejected);
		failures+=run("graphKUnbranchedOverlap", CrossKTipOverlapperUnitTest::graphKUnbranchedOverlap);
		failures+=run("graphKBranchIgnored", CrossKTipOverlapperUnitTest::graphKBranchIgnored);
		failures+=run("selfOverlapIsAmbiguous", CrossKTipOverlapperUnitTest::selfOverlapIsAmbiguous);
		failures+=run("trimmedFlankOrientations", CrossKTipOverlapperUnitTest::trimmedFlankOrientations);
		failures+=run("trimmedFlankBudget", CrossKTipOverlapperUnitTest::trimmedFlankBudget);
		failures+=run("trimmedFlankCoverage", CrossKTipOverlapperUnitTest::trimmedFlankCoverage);
		failures+=run("guardedInwardOverlap", CrossKTipOverlapperUnitTest::guardedInwardOverlap);
		failures+=run("realRepeatTrimRejected", CrossKTipOverlapperUnitTest::realRepeatTrimRejected);
		failures+=run("changedContextDeclinesMerge", CrossKTipOverlapperUnitTest::changedContextDeclinesMerge);
		failures+=run("historicalDeadEnds", CrossKTipOverlapperUnitTest::historicalDeadEnds);
		failures+=run("fusionEndpointOrientations", CrossKTipOverlapperUnitTest::fusionEndpointOrientations);
		failures+=run("conflictingRepeatPlacement", CrossKTipOverlapperUnitTest::conflictingRepeatPlacement);
		failures+=run("unconflictedFusionPreserved", CrossKTipOverlapperUnitTest::unconflictedFusionPreserved);
		failures+=run("coverageRatioVeto", CrossKTipOverlapperUnitTest::coverageRatioVeto);
		BubblePopper.crossKMerge=false;
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

	private static void uniqueReciprocalOverlap(){
		Contig a=contig(0, "AAAACCCCGGGG", false, true);
		Contig b=contig(1, "CCCGGGGTTTT", true, false);
		ArrayList<Contig> contigs=list(a, b);
		int pairs=new CrossKTipOverlapper(contigs, 5, 9).addEdges();
		check(pairs==1, "Expected one overlap pair, got "+pairs);
		check(a.rightEdgeCount()==1 && b.leftEdgeCount()==1, "Missing reciprocal overlap edges");
		check(a.rightEdges.get(0).overlap==7, "Wrong overlap length: "+a.rightEdges.get(0).overlap);
		int merged=popper(contigs).expand(a);
		check(merged==1, "Expected one merge, got "+merged);
		check(new String(a.bases).equals("AAAACCCCGGGGTTTT"), "Wrong merged sequence: "+new String(a.bases));
	}

	private static void inwardTipOverlap(){
		Contig a=contig(0, "AAAACCCCGGGTT", false, true);
		Contig b=contig(1, "TTCCCCGGGAAAA", true, false);
		a.rightCode=Tadpole.F_BRANCH;
		b.leftCode=Tadpole.B_BRANCH;
		ArrayList<Contig> contigs=list(a, b);
		int pairs=new CrossKTipOverlapper(contigs, 5, 9).addEdges();
		check(pairs==1, "Expected one inward overlap pair, got "+pairs);
		Edge edge=a.rightEdges.get(0);
		check(edge.overlap==7 && edge.sourceTrim==2 && edge.destTrim==2,
				"Wrong inward overlap geometry: "+edge);
		int merged=popper(contigs).expand(a);
		check(merged==1, "Expected one inward merge, got "+merged);
		check(new String(a.bases).equals("AAAACCCCGGGAAAA"), "Wrong inward merged sequence: "+new String(a.bases));
	}

	private static void branchTipIgnored(){
		Contig a=contig(0, "AAAACCCCGGGG", false, false);
		Contig b=contig(1, "CCCGGGGTTTT", true, false);
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 5, 9).addEdges()==0, "Ineligible branch tip was joined");
	}

	private static void equalBestIsAmbiguous(){
		Contig a=contig(0, "AAAACCCCGGGG", false, true);
		Contig b=contig(1, "CCCGGGGTTTT", true, false);
		Contig c=contig(2, "CCCGGGGAAAA", true, false);
		ArrayList<Contig> contigs=list(a, b, c);
		check(new CrossKTipOverlapper(contigs, 5, 9).addEdges()==0, "Equal best overlaps were not rejected");
	}

	private static void discardedTipDoesNotBlockJoin(){
		Contig a=contig(0, "AAAACCCCGGGG", false, true);
		Contig b=contig(1, "CCCGGGGTTTT", true, false);
		Contig debris=contig(2, "CCCGGGGA", true, false);
		ArrayList<Contig> contigs=list(a, b, debris);
		check(new CrossKTipOverlapper(contigs, 5, 9, false, 10).addEdges()==1,
				"Discarded short tip blocked a retained-contig join");
	}

	private static void repeatedAnchorOnOneTipIsAmbiguous(){
		Contig a=contig(0, "CCCCCCCCAAAA", false, true);
		Contig b=contig(1, "AAAGGG", true, false);
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 3, 4).addEdges()==0,
				"A repeated anchor on one tip was treated as unique");
	}

	private static void reverseOrientedOverlap(){
		Contig a=contig(0, "AAAACCCCGGGG", false, true);
		Contig b=contig(1, "AAAACCCCGGG", false, true);
		ArrayList<Contig> contigs=list(a, b);
		int pairs=new CrossKTipOverlapper(contigs, 5, 9).addEdges();
		check(pairs==1, "Expected one reverse-oriented pair, got "+pairs);
		check(a.rightEdges.get(0).destRight(), "Destination orientation was not right-facing");
		int merged=popper(contigs).expand(a);
		check(merged==1, "Expected one reverse-oriented merge, got "+merged);
		check(new String(a.bases).equals("AAAACCCCGGGGTTTT"), "Wrong reverse-oriented product: "+new String(a.bases));
	}

	private static void cyclicComponentRejected(){
		Contig a=contig(0, "GGATTAAC", true, true);
		Contig b=contig(1, "AACGGCCT", true, true);
		Contig c=contig(2, "CCTCCGGA", true, true);
		ArrayList<Contig> contigs=list(a, b, c);
		check(new CrossKTipOverlapper(contigs, 3, 3).addEdges()==0, "Cyclic overlap component was not rejected");
		for(Contig x : contigs){check(x.leftEdgeCount()==0 && x.rightEdgeCount()==0, "Cycle left graph edges");}
	}

	private static void graphKUnbranchedOverlap(){
		Contig a=contig(0, "AAAACCCCGGGG", false, false);
		Contig b=contig(1, "CCCGGGGTTTT", false, false);
		a.rightCode=Tadpole.KEEP_GOING;
		b.leftCode=Tadpole.KEEP_GOING;
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 5, 9, true).addEdges()==1,
				"Unbranched graph-k overlap was not selected");
	}

	private static void graphKBranchIgnored(){
		Contig a=contig(0, "AAAACCCCGGGG", false, false);
		Contig b=contig(1, "CCCGGGGTTTT", false, false);
		a.rightCode=Tadpole.F_BRANCH;
		b.leftCode=Tadpole.KEEP_GOING;
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 5, 9, true).addEdges()==0,
				"Branched graph-k overlap was selected");
	}

	private static void selfOverlapIsAmbiguous(){
		Contig a=contig(0, "AAACCCAAA", false, false);
		Contig b=contig(1, "AAAGGG", false, false);
		a.leftCode=a.rightCode=Tadpole.KEEP_GOING;
		b.leftCode=Tadpole.KEEP_GOING;
		ArrayList<Contig> contigs=list(a, b);
		check(new CrossKTipOverlapper(contigs, 3, 3, true).addEdges()==0,
				"A terminal self-overlap was ignored when selecting an external join");
	}

	/** Covers all stored orientations and the reciprocal direction of the same splice. */
	private static void trimmedFlankOrientations(){
		for(int orientation=0; orientation<4; orientation++){
			final boolean ar=(orientation&1)!=0, br=(orientation&2)!=0;
			final Contig a=oriented("GGACGTCAGTA", ar);
			final Contig b=oriented("ACGTCAGTACC", br);
			check(CrossKTipOverlapper.compatibleTrimmedFlanks(a, ar, 2, b, br, 2, 5, 0),
					"Compatible flanks rejected in orientation "+orientation);
			final Contig changed=oriented("GGACGTCAGTG", ar);
			for(int allowance=0; allowance<=1; allowance++){
				final boolean forward=CrossKTipOverlapper.compatibleTrimmedFlanks(changed, ar, 2, b, br, 2, 5, allowance);
				final boolean reverse=CrossKTipOverlapper.compatibleTrimmedFlanks(b, !br, 2, changed, !ar, 2, 5, allowance);
				check(forward==(allowance==1) && reverse==forward,
						"Mismatch budget or reciprocal geometry differs in orientation "+orientation);
			}
		}
	}

	/** The allowance is shared across both flanks, and N/N consumes a mismatch. */
	private static void trimmedFlankBudget(){
		final Contig a=oriented("GGACGTCAGTG", false);
		final Contig b=oriented("ATGTCAGTACC", false);
		check(!CrossKTipOverlapper.compatibleTrimmedFlanks(a, false, 2, b, false, 2, 5, 1),
				"One mismatch on each side incorrectly counted as one total");
		check(CrossKTipOverlapper.compatibleTrimmedFlanks(a, false, 2, b, false, 2, 5, 2),
				"Two allowed flank mismatches rejected");
		final Contig unknownA=oriented("GGACGTCAGNA", false);
		final Contig unknownB=oriented("ACGTCAGNACC", false);
		check(!CrossKTipOverlapper.compatibleTrimmedFlanks(unknownA, false, 2, unknownB, false, 2, 5, 0),
				"Matching unknown bases were treated as positive evidence");
		check(CrossKTipOverlapper.compatibleTrimmedFlanks(unknownA, false, 2, unknownB, false, 2, 5, 1),
				"One unknown comparison must consume exactly one allowed disagreement");
	}

	/** Every discarded base must have a corresponding base in the other contig. */
	private static void trimmedFlankCoverage(){
		final Contig a=oriented("GGACGTCAGTAC", false);
		final Contig b=oriented("ACGTCAGTA", false);
		check(!CrossKTipOverlapper.compatibleTrimmedFlanks(a, false, 3, b, false, 2, 5, 10),
				"Uncovered source tail was silently ignored");
		final Contig shortA=oriented("GTCAGTA", false);
		check(!CrossKTipOverlapper.compatibleTrimmedFlanks(shortA, false, 2, b, false, 2, 5, 10),
				"Uncovered destination prefix was silently ignored");
	}

	/** The historical inward-tip fixture discards four conflicting bases. */
	private static void guardedInwardOverlap(){
		for(int allowance=0; allowance<=1; allowance++){
			final Contig a=contig(0, "AAAACCCCGGGTT", false, true);
			final Contig b=contig(1, "TTCCCCGGGAAAA", true, false);
			check(new CrossKTipOverlapper(list(a, b), 5, 9, false, 0, allowance).addEdges()==0,
					"Conflicting inward-tip fixture survived allowance "+allowance);
		}
		final Contig a=contig(0, "AAAACCCCGGGG", false, true);
		final Contig b=contig(1, "CCCGGGGTTTT", true, false);
		check(new CrossKTipOverlapper(list(a, b), 5, 9, false, 0, 0).addEdges()==1,
				"Strict validation changed an ordinary untrimmed exact fusion");
	}

	/** Uses the observed M.ruber 55bp repeat and both contradictory continuations. */
	private static void realRepeatTrimRejected(){
		final String anchor="ATTCAAGCCGACCGAAGGGAGTAGAAAAGCCTTTCGGTAGTATCGTTTAGGCTTG";
		final String source="GGTGAAGG"+anchor+"CCACAGTGAACGATACTACCGAAATGCGTATGAGAGACT";
		final String dest=anchor+"TCAAAGTGAACAATACTACCGAAATGCGTATGAAGTCCG"+"GCTGC";
		for(int allowance=-1; allowance<=1; allowance++){
			final Contig a=contig(0, source, false, true), b=contig(1, dest, true, false);
			final int pairs=new CrossKTipOverlapper(list(a, b), 32, 94, false, 0, allowance).addEdges();
			check(pairs==(allowance<0 ? 1 : 0), "Real repeat fixture has wrong fusion count: "+pairs);
		}
	}

	/** Stale anchors and newly contradictory tails must decline before mutating either contig. */
	private static void changedContextDeclinesMerge(){
		final int oldAllowance=BubblePopper.crossKMaxMismatches;
		BubblePopper.crossKMaxMismatches=1;
		try{
			for(int changed=0; changed<3; changed++){
				final Contig a=contig(0, "GGACGTCAGTG", false, true);
				final Contig b=contig(1, "ACGTCAGTACC", true, false);
				a.addRightEdge(new Edge(0, 1, 0, 1, 20, null, 5, 2, 2));
				b.addLeftEdge(new Edge(1, 0, 0, 2, 20, null, 5, 2, 2));
				if(changed==1){b.bases[1]='T';}//A second disagreement in the opposite flank.
				if(changed==2){b.bases[4]='A';}//The exact anchor itself is no longer valid.
				final BubblePopper popper=popper(list(a, b));
				final int merged=popper.expand(a);
				check(merged==(changed==0 ? 1 : 0), "Stale fusion was not rechecked: "+changed);
				if(changed>0){
					check(!b.used() && a.length()==11 && popper.crossKRejectedFlanks>0,
							"Declined stale edge changed sequence ownership");
				}
			}
		}finally{
			BubblePopper.crossKMaxMismatches=oldAllowance;
		}
	}

	/** Clearing/reclassifying the graph must not manufacture original dead-end evidence. */
	private static void historicalDeadEnds(){
		for(boolean graphK : new boolean[]{false, true}){
			for(int condition=0; condition<3; condition++){
				final Contig a=contig(0, "AAAACCCCGGGG", false, true);
				final Contig b=contig(1, "CCCGGGGTTTT", true, false);
				if(condition==1){a.rightCode=Tadpole.F_BRANCH;}
				if(condition==2){a.addRightEdge(new Edge(0, 1, 1, 1, 20, new byte[]{'A'}));}
				a.markFusionEndpoints();
				b.markFusionEndpoints();
				a.rightEdges=null;
				if(graphK){
					a.leftCode=b.rightCode=Tadpole.F_BRANCH;
					a.rightCode=b.leftCode=Tadpole.KEEP_GOING;
				}
				final int pairs=new CrossKTipOverlapper(list(a, b), 5, 9, graphK, 0, 0, true).addEdges();
				check(pairs==(condition==0 ? 1 : 0),
						"Historical dead-end gate failed: condition="+condition+", graphK="+graphK);
				check(a.rightBridgeEndpoint && b.leftBridgeEndpoint,
						"Fusion restriction altered independent bridge eligibility");
			}
		}
	}

	/** Both orientation APIs and actual reverse merges preserve only surviving outer flags. */
	private static void fusionEndpointOrientations(){
		final Contig c=contig(0, "ACGTTGCA", false, true);
		c.leftFusionEndpoint=true;
		c.rcomp();
		check(!c.leftFusionEndpoint && c.rightFusionEndpoint, "rcomp did not swap fusion endpoints");
		c.rcomp();
		c.flip(null);
		check(!c.leftFusionEndpoint && c.rightFusionEndpoint, "flip did not swap fusion endpoints");
		c.flip(null);
		check(c.leftFusionEndpoint && !c.rightFusionEndpoint, "Double flip lost endpoint state");
		for(boolean reverse : new boolean[]{false, true}){
			final Contig a=contig(0, "AAAACCCCGGGG", false, true);
			final Contig b=contig(1, "CCCGGGGTTTT", true, false);
			a.leftFusionEndpoint=false;
			a.rightFusionEndpoint=b.leftFusionEndpoint=b.rightFusionEndpoint=true;
			if(reverse){b.rcomp();}
			final ArrayList<Contig> contigs=list(a, b);
			check(new CrossKTipOverlapper(contigs, 5, 9, false, 0, 0, true).addEdges()==1,
					"Eligible oriented ends failed to select a fusion");
			check(popper(contigs).expand(a)==1, "Eligible oriented fusion failed to merge");
			check(!a.leftFusionEndpoint && a.rightFusionEndpoint,
					"Merged eligibility did not follow the surviving outer ends");
		}
	}

	/** Real Spirulina event28: a shorter exact anchor has a larger contradictory placement span. */
	private static void conflictingRepeatPlacement(){
		final String source="TTCCCCTTTTTAAGGGGGGAGCCGCTCAAAGTCCCCCTTAAAAATAGGGGAGCCGCTCAAAGTCCCCCTTTTTAAGGGGGGAGCCGCTCAAAGTCCCCCTTTTTAAGGGGGGAGCCGCTCAAAGTCCCCCTTTTTAAGGG";
		final String dest="CCGCTCAAAGTCCCCCTTTTTAAGGGGGGAGCCGCTCAAAGTCCCCCTTTTTAAGGGGGGAGCCGCTCAAAGTCCCCCTTTTTAAGGGGGGAGCTGCTCAAAGTCCCCCTTTTTAAGGGGGATTTAGGGGGATCGATCGC";
		for(int orientation=0; orientation<4; orientation++){
			for(boolean swap : new boolean[]{false, true}){
				for(boolean enabled : new boolean[]{false, true}){
					final Contig a=contig(0, source, false, true), b=contig(1, dest, true, false);
					if((orientation&1)!=0){a.rcomp();}
					if((orientation&2)!=0){b.rcomp();}
					if(swap){a.id=1; b.id=0;}
					final ArrayList<Contig> contigs=swap ? list(b, a) : list(a, b);
					final int pairs=new CrossKTipOverlapper(contigs, 64, 94, false, 0, 1, false, enabled).addEdges();
					check(pairs==(enabled ? 0 : 1),
							"Conflict veto changed with orientation/order: "+orientation+", "+swap+", "+enabled);
					if(!enabled){
						final Edge edge=(orientation&1)==0 ? a.rightEdges.get(0) : a.leftEdges.get(0);
						check(edge.overlap==88 && edge.sourceTrim==0 && edge.destTrim==0,
								"Real repeat fixture no longer reproduces the observed 88bp fusion");
					}
				}
			}
		}
	}

	/** Ordinary exact joins remain available in initial and final graph-K overlap discovery. */
	private static void unconflictedFusionPreserved(){
		for(boolean graphK : new boolean[]{false, true}){
			final Contig a=contig(0, "AAAACCCCGGGG", false, true);
			final Contig b=contig(1, "CCCGGGGTTTT", true, false);
			if(graphK){a.rightCode=b.leftCode=Tadpole.KEEP_GOING;}
			check(new CrossKTipOverlapper(list(a, b), 5, 9, graphK, 0, 1, false, true).addEdges()==1,
					"Conflict policy removed an unconflicted exact join");
		}
	}

	/** Checks off/default behavior, inclusive boundaries, missing depth, and both graph routes/strands. */
	private static void coverageRatioVeto(){
		check(CrossKTipOverlapper.compatibleCoverage(20, 80, 0), "Disabled veto changed legacy behavior");
		check(CrossKTipOverlapper.compatibleCoverage(300, 300, 1), "High absolute depth was treated as a mismatch");
		check(!CrossKTipOverlapper.compatibleCoverage(0, 0, 2), "Missing depths supplied false copy-number evidence");
		for(boolean graphK : new boolean[]{false, true}){
			for(int orientation=0; orientation<4; orientation++){
				for(float depth : new float[]{20, 35, Math.nextUp(35f), 40, 0}){
					final Contig a=contig(0, "AAAACCCCGGGG", false, true);
					final Contig b=contig(1, "CCCGGGGTTTT", true, false);
					b.coverage=depth;
					if(graphK){a.rightCode=b.leftCode=Tadpole.KEEP_GOING;}
					if((orientation&1)!=0){a.rcomp();}
					if((orientation&2)!=0){b.rcomp();}
					final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(list(a, b), 5, 9, graphK);
					overlapper.maxCoverageRatio=1.75f;
					final int expected=(depth>0 && depth<=35 ? 1 : 0);
					check(overlapper.addEdges()==expected, "Coverage veto changed with depth/strand/graph route: "+
							depth+", "+orientation+", "+graphK);
					if(expected==0){
						check(a.leftEdgeCount()+a.rightEdgeCount()+b.leftEdgeCount()+b.rightEdgeCount()==0,
								"A rejected depth-discontinuous pair left live graph edges");
					}
				}
			}
		}
	}

	/** Stores a requested oriented sequence forward or reverse-complemented. */
	private static Contig oriented(final String sequence, final boolean reverse){
		final byte[] bases=sequence.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		if(reverse){AminoAcid.reverseComplementBasesInPlace(bases);}
		return new Contig(bases);
	}
	private static BubblePopper popper(ArrayList<Contig> contigs){
		HashMap<Integer, ArrayList<Edge>> map=new HashMap<Integer, ArrayList<Edge>>();
		for(Contig c : contigs){
			if(c.leftEdges!=null){for(Edge e : c.leftEdges){add(map, e);}}
			if(c.rightEdges!=null){for(Edge e : c.rightEdges){add(map, e);}}
		}
		return new BubblePopper(contigs, map, 5);
	}

	private static void add(HashMap<Integer, ArrayList<Edge>> map, Edge e){
		ArrayList<Edge> list=map.get(e.destination);
		if(list==null){list=new ArrayList<Edge>(); map.put(e.destination, list);}
		list.add(e);
	}

	private static Contig contig(int id, String bases, boolean left, boolean right){
		Contig c=new Contig(bases.getBytes(), id);
		c.coverage=20;
		c.leftCode=c.rightCode=Tadpole.DEAD_END;
		c.leftBridgeEndpoint=left;
		c.rightBridgeEndpoint=right;
		return c;
	}

	private static ArrayList<Contig> list(Contig... contigs){
		ArrayList<Contig> list=new ArrayList<Contig>();
		for(Contig c : contigs){list.add(c);}
		return list;
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private interface Test {void run();}
}
