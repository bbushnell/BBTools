package assemble;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Executable spec for Tadpole's path-preserving DBG bubble unzip (Noelle, lead;
 * production: {@link BubblePopper#unzipTrueBubble}, gated by BubblePopper.unzipBubbles,
 * default off).
 *
 * Adapted 2026-08-26 (Noelle's phase-2 bounded ask) to exercise the REAL default-off
 * production path directly: every motif builds a genuine Contig/Edge graph, sets the
 * real BubblePopper static flags, and calls the real BubblePopper.expand() -- no more
 * "ADAPTER HOOK" placeholders. Expected sequences/coverage are still computed by an
 * INDEPENDENT oracle (spliceHaplotype/pathCoverageOracle) and compared against what
 * production actually produced, rather than trusting production's own output.
 *
 * Per Noelle's file-ownership grant this session, edits are scoped to this file only --
 * no production class or BubblePopperUnitTest.java is touched.
 *
 * @author Nepgear
 */
public class PathPreservingBubbleSimplifierSpec {

	/** Mirrors Tadpole's errorPath=1 defaults -- the same mirror BubblePopperUnitTest.java
	 * uses. Keep in sync if isError's tuning is ever retuned (flagged by Amber, 2026-08-26). */
	private static final BubblePopper.ErrorClassifier DEFAULT_ERROR_CLASSIFIER=
			(high, low) -> low*16f<high || (low<=4 && high>=Math.max(3, low*2.6f));

	public static void main(String[] args){
		int checks=0, failures=0;

		checks++; failures+=preserveTwoAlleles() ? 0 : 1;
		checks++; failures+=rejectExternalBranch() ? 0 : 1;
		checks++; failures+=rejectBackwardBranch() ? 0 : 1;
		checks++; failures+=rejectCycle() ? 0 : 1;
		checks++; failures+=rejectMalformedEdge() ? 0 : 1;
		checks++; failures+=preserveAlleleCoverage() ? 0 : 1;

		System.err.println("Motifs checked: "+checks);
		System.err.println("Failures: "+failures);
		if(failures==0){System.err.println("ALL PASS");}
		if(failures>0){throw new RuntimeException("PathPreservingBubbleSimplifierSpec: "+failures+" failure(s)");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Motif 1            ----------------*/
	/*--------------------------------------------------------------*/

	/** Clean, fully isolated bubble: comparable-depth arms must unzip into two terminal products. */
	static boolean preserveTwoAlleles(){
		final int k=5;
		// Comparable arm depths (25 vs 22) -- a TRUE bubble under Brian's clarified charter:
		// existing BubblePopper.pop() collapsing a LOW-depth alternate is wanted soap-bubble
		// cleanup, not a misfeature. Only comparable-depth alternates unzip.
		IsolatedBubble fx=buildIsolatedBubble(k, 20, 25, 22, 20);
		String hap1=spliceHaplotype(fx.L, fx.M1, fx.R, k);
		String hap2=spliceHaplotype(fx.L, fx.M2, fx.R, k);

		configureUnzip();
		int expansions=new BubblePopper(fx.allContigs, fx.destMap, k, DEFAULT_ERROR_CLASSIFIER).expand(fx.L);

		boolean pass=true;
		pass&=check("preserveTwoAlleles: exactly one unzip", 1, expansions);
		pass&=check("preserveTwoAlleles: representative arm retired", true, fx.M1.used());
		pass&=check("preserveTwoAlleles: alternate arm retired", true, fx.M2.used());
		// M1 (depth 25) is the higher-depth representative, so L holds its spliced product.
		pass&=check("preserveTwoAlleles: left product is the representative's haplotype", hap1, new String(fx.L.bases));
		pass&=check("preserveTwoAlleles: right product is the alternate's haplotype", hap2, new String(fx.R.bases));
		pass&=check("preserveTwoAlleles: both products present and distinct", true, !new String(fx.L.bases).equals(new String(fx.R.bases)));
		pass&=check("preserveTwoAlleles: left product allele base", 'C', (char)fx.L.bases[10]);
		pass&=check("preserveTwoAlleles: right product allele base", 'T', (char)fx.R.bases[10]);
		pass&=check("preserveTwoAlleles: both products are edge-free terminal contigs",
				true, fx.L.rightEdges==null && fx.R.leftEdges==null);
		return pass;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Motif 2            ----------------*/
	/*--------------------------------------------------------------*/

	/** Repeat: a mid node reached by a second, unrelated source is ambiguous -- must be declined. */
	static boolean rejectExternalBranch(){
		final int k=5;
		IsolatedBubble fx=buildIsolatedBubble(k, 20, 25, 22, 20);
		byte[] originalLBases=fx.L.bases.clone(), originalRBases=fx.R.bases.clone();

		// L2 independently reaches M1 too (same trailing (k-1)-mer "AAAA" as L).
		Contig L2=new Contig("CCCCCAAAAA".getBytes(), fx.nextId());
		Edge l2ToM1=new Edge(L2.id, fx.M1.id, 1, orient(true,false), 1, new byte[]{'C'});
		Edge m1ToL2=new Edge(fx.M1.id, L2.id, 1, orient(false,true), 1, null);
		L2.addRightEdge(l2ToM1);
		fx.M1.addLeftEdge(m1ToL2);
		fx.allContigs.add(L2);
		addInbound(fx.destMap, l2ToM1);
		addInbound(fx.destMap, m1ToL2);

		int m1LeftInboundSources=countDistinctInboundSources(fx.destMap, fx.M1.id, false);

		configureUnzip();
		int expansions=new BubblePopper(fx.allContigs, fx.destMap, k, DEFAULT_ERROR_CLASSIFIER).expand(fx.L);

		boolean pass=true;
		pass&=check("rejectExternalBranch: M1 has 2 independent left-inbound sources (ambiguous)", 2, m1LeftInboundSources);
		pass&=check("rejectExternalBranch: real expand() declined (0 expansions)", 0, expansions);
		pass&=check("rejectExternalBranch: L bases untouched", true, java.util.Arrays.equals(originalLBases, fx.L.bases));
		pass&=check("rejectExternalBranch: R bases untouched", true, java.util.Arrays.equals(originalRBases, fx.R.bases));
		pass&=check("rejectExternalBranch: neither arm retired", true, !fx.M1.used() && !fx.M1.associate() && !fx.M2.used() && !fx.M2.associate());
		return pass;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Motif 3            ----------------*/
	/*--------------------------------------------------------------*/

	/** Backward branch: B_BRANCH must never be treated as a poppable forward-bubble side. */
	static boolean rejectBackwardBranch(){
		final int k=5;
		IsolatedBubble fx=buildIsolatedBubble(k, 20, 25, 22, 20);
		fx.L.rightCode=Tadpole.B_BRANCH;
		byte[] originalLBases=fx.L.bases.clone();

		boolean pass=true;
		pass&=check("rejectBackwardBranch: rightForwardBranch() is false for B_BRANCH", false, fx.L.rightForwardBranch());
		pass&=check("rejectBackwardBranch: rightBackwardBranch() is true", true, fx.L.rightBackwardBranch());
		pass&=check("rejectBackwardBranch: isBranchCode() still true (shares BRANCH_BIT)", true, fx.L.rightBranch());

		configureUnzip();
		int expansions=new BubblePopper(fx.allContigs, fx.destMap, k, DEFAULT_ERROR_CLASSIFIER).expand(fx.L);
		pass&=check("rejectBackwardBranch: real expand() attempts zero expansions", 0, expansions);
		pass&=check("rejectBackwardBranch: L bases untouched", true, java.util.Arrays.equals(originalLBases, fx.L.bases));
		pass&=check("rejectBackwardBranch: neither arm retired", true, !fx.M1.used() && !fx.M2.used());
		return pass;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Motif 4            ----------------*/
	/*--------------------------------------------------------------*/

	/** Cycle/self-loop: must be detected and left structurally alone, never linearized through. */
	static boolean rejectCycle(){
		Contig c=new Contig("AAAAACCCCC".getBytes(), 0);
		byte[] originalBases=c.bases.clone();
		c.leftEdges=new ArrayList<Edge>();
		c.leftEdges.add(new Edge(c.id, c.id, 1, orient(false,true), 5, null));
		c.rightEdges=new ArrayList<Edge>();
		c.rightEdges.add(new Edge(c.id, c.id, 1, orient(true,false), 5, null));

		HashMap<Integer, ArrayList<Edge>> destMap=new HashMap<Integer, ArrayList<Edge>>();
		ArrayList<Edge> inbound=new ArrayList<Edge>();
		inbound.add(c.leftEdges.get(0));
		inbound.add(c.rightEdges.get(0));
		destMap.put(c.id, inbound);

		ArrayList<Contig> allContigs=new ArrayList<Contig>();
		allContigs.add(c);
		configureUnzip();
		BubblePopper bp=new BubblePopper(allContigs, destMap, 5, DEFAULT_ERROR_CLASSIFIER);

		boolean pass=true;
		pass&=check("rejectCycle: isLoop() detects the self-loop", true, bp.isLoop(c));

		int expansions=bp.expand(c);
		// With configureUnzip() (popDirect=false), and this contig never marked as a branch
		// (leftCode/rightCode default to 0, not F_BRANCH), neither the unzip nor the legacy
		// popDirect flip path is even attempted -- a genuinely untouched decline.
		pass&=check("rejectCycle: real expand() makes zero structural changes (no merge/pop)", 0, expansions);
		pass&=check("rejectCycle: bases byte-identical to the original (opt-in unzip path never flips here)",
				true, java.util.Arrays.equals(originalBases, c.bases));
		pass&=check("rejectCycle: edges structurally untouched (still self-loop)", true, bp.isLoop(c));
		return pass;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Motif 5            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Asymmetric/malformed edge: an arm's forward (entry) edge and its reciprocal (back) edge
	 * disagree in depth -- a real one-sided-connection malformation, not just a missing edge
	 * (that simpler case is already caught upstream by midNodesConcur's leftEdges-null check).
	 * Exercises unzipTrueBubble's own validBubbleArm reciprocity check
	 * (entry.depth!=back.depth -> decline) directly.
	 */
	static boolean rejectMalformedEdge(){
		final int k=5;
		IsolatedBubble fx=buildIsolatedBubble(k, 20, 25, 22, 20);
		byte[] originalLBases=fx.L.bases.clone(), originalRBases=fx.R.bases.clone();

		if(fx.M1.leftEdgeCount()!=1){throw new RuntimeException("Fixture sanity failed: M1 must have exactly one back edge before corruption");}
		// Corrupt M1's back edge (M1->L) depth so it no longer matches the entry edge (L->M1) --
		// same claimed connection, inconsistent reciprocal record.
		fx.M1.leftEdges.get(0).depth=24;

		configureUnzip();
		int expansions=new BubblePopper(fx.allContigs, fx.destMap, k, DEFAULT_ERROR_CLASSIFIER).expand(fx.L);

		boolean pass=true;
		pass&=check("rejectMalformedEdge: real expand() declined (0 expansions)", 0, expansions);
		pass&=check("rejectMalformedEdge: L bases untouched", true, java.util.Arrays.equals(originalLBases, fx.L.bases));
		pass&=check("rejectMalformedEdge: R bases untouched", true, java.util.Arrays.equals(originalRBases, fx.R.bases));
		pass&=check("rejectMalformedEdge: neither arm retired", true, !fx.M1.used() && !fx.M1.associate() && !fx.M2.used() && !fx.M2.associate());
		return pass;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Motif 6            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * True bubble, near-equal arms (25 vs 22): each output's coverage derives from its own arm's
	 * full coverage PLUS the shared L/R flanks apportioned by entry/exit depth fraction (matching
	 * BubblePopper.pathCoverage's actual formula in this tree exactly: leftFraction/rightFraction
	 * = this arm's entry/exit depth over the SUM of both arms' entry/exit depth -- since L and R
	 * are shared between both output products, their coverage must be split, not double-counted
	 * on both haplotypes). Computed here by an independent oracle and compared against what
	 * production actually wrote -- confirmed against production's real numbers, not assumed.
	 */
	static boolean preserveAlleleCoverage(){
		final int k=5;
		final float lCov=20, repCov=25, altCov=22, rCov=20;
		IsolatedBubble fx=buildIsolatedBubble(k, lCov, repCov, altCov, rCov);
		String hap1=spliceHaplotype(fx.L, fx.M1, fx.R, k);
		String hap2=spliceHaplotype(fx.L, fx.M2, fx.R, k);

		// Entry/exit depths here are (int)repCov=25 and (int)altCov=22 on both boundaries (see
		// buildIsolatedBubble), so entryDepthSum==exitDepthSum==47.
		final double entryDepthSum=repCov+altCov, exitDepthSum=repCov+altCov;
		final double repLeftFraction=repCov/entryDepthSum, repRightFraction=repCov/exitDepthSum;
		final double altLeftFraction=altCov/entryDepthSum, altRightFraction=altCov/exitDepthSum;
		// Edge terms are entry.depth*(entry.length-1)+exit.depth*(exit.length-1); both edges here
		// are single-base (length=1), so those terms vanish and only the flank+mid terms remain.
		final double repExpected=(lCov*(fx.L.length()-k+1)*repLeftFraction + repCov*(fx.M1.length()-k+1)
				+ rCov*(fx.R.length()-k+1)*repRightFraction)/(double)(hap1.length()-k+1);
		final double altExpected=(lCov*(fx.L.length()-k+1)*altLeftFraction + altCov*(fx.M2.length()-k+1)
				+ rCov*(fx.R.length()-k+1)*altRightFraction)/(double)(hap2.length()-k+1);

		configureUnzip();
		int expansions=new BubblePopper(fx.allContigs, fx.destMap, k, DEFAULT_ERROR_CLASSIFIER).expand(fx.L);

		boolean pass=true;
		pass&=check("preserveAlleleCoverage: exactly one unzip", 1, expansions);
		pass&=check("preserveAlleleCoverage: representative product coverage matches independent oracle",
				repExpected, fx.L.coverage, 0.01);
		pass&=check("preserveAlleleCoverage: alternate product coverage matches independent oracle",
				altExpected, fx.R.coverage, 0.01);
		pass&=check("preserveAlleleCoverage: the two products' coverage remain distinct, neither blended",
				true, fx.L.coverage!=fx.R.coverage);
		return pass;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Fixture construction     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Fully isolated 2-arm bubble satisfying BubblePopper.unzipTrueBubble's complete contract:
	 * L/R have no outer edges (DEAD_END on the outside, no third-party inbound), exactly 2 arms
	 * on the inside, each arm single-edge each side, all four junction edges reciprocal (matching
	 * length/depth on both directions) and carrying the correct overlap byte(s) so
	 * validForwardEdge's sequence check passes. Same L/M1/M2/R base sequences and orientation
	 * convention as the prior buildBubble() (verified this session against production's
	 * validForwardEdge: shared=kbig-edge.length overlap bytes, then edge.bases must equal the
	 * next edge.length bytes of the target).
	 */
	static IsolatedBubble buildIsolatedBubble(int k, float lCov, float repCov, float altCov, float rCov){
		IsolatedBubble fx=new IsolatedBubble();
		fx.L=new Contig("TTTTTAAAAA".getBytes(), fx.nextId());
		fx.M1=new Contig("AAAACGGGG".getBytes(), fx.nextId());
		fx.M2=new Contig("AAAATGGGG".getBytes(), fx.nextId());
		fx.R=new Contig("GGGGGCCCCC".getBytes(), fx.nextId());

		fx.L.coverage=lCov; fx.M1.coverage=repCov; fx.M2.coverage=altCov; fx.R.coverage=rCov;
		fx.L.leftCode=Tadpole.DEAD_END; fx.L.rightCode=Tadpole.F_BRANCH;
		fx.M1.leftCode=Tadpole.F_BRANCH; fx.M1.rightCode=Tadpole.F_BRANCH;
		fx.M2.leftCode=Tadpole.F_BRANCH; fx.M2.rightCode=Tadpole.F_BRANCH;
		fx.R.leftCode=Tadpole.F_BRANCH; fx.R.rightCode=Tadpole.DEAD_END;

		// Entry edges L->M1/L->M2: single divergent overlap byte, taken from each mid's own sequence.
		Edge lToM1=new Edge(fx.L.id, fx.M1.id, 1, orient(true,false), (int)repCov, new byte[]{'C'});
		Edge lToM2=new Edge(fx.L.id, fx.M2.id, 1, orient(true,false), (int)altCov, new byte[]{'T'});
		fx.L.addRightEdge(lToM1);
		fx.L.addRightEdge(lToM2);
		Edge m1Back=new Edge(fx.M1.id, fx.L.id, 1, orient(false,true), (int)repCov, null);
		Edge m2Back=new Edge(fx.M2.id, fx.L.id, 1, orient(false,true), (int)altCov, null);
		fx.M1.addLeftEdge(m1Back);
		fx.M2.addLeftEdge(m2Back);

		// Exit edges M1->R/M2->R: no divergence at this boundary in these sequences (both overlaps are 'G').
		Edge m1ToR=new Edge(fx.M1.id, fx.R.id, 1, orient(true,false), (int)repCov, new byte[]{'G'});
		Edge m2ToR=new Edge(fx.M2.id, fx.R.id, 1, orient(true,false), (int)altCov, new byte[]{'G'});
		fx.M1.addRightEdge(m1ToR);
		fx.M2.addRightEdge(m2ToR);
		Edge rToM1=new Edge(fx.R.id, fx.M1.id, 1, orient(false,true), (int)repCov, null);
		Edge rToM2=new Edge(fx.R.id, fx.M2.id, 1, orient(false,true), (int)altCov, null);
		fx.R.addLeftEdge(rToM1);
		fx.R.addLeftEdge(rToM2);

		fx.allContigs.add(fx.L); fx.allContigs.add(fx.M1); fx.allContigs.add(fx.M2); fx.allContigs.add(fx.R);
		addInbound(fx.destMap, lToM1); addInbound(fx.destMap, lToM2);
		addInbound(fx.destMap, m1Back); addInbound(fx.destMap, m2Back);
		addInbound(fx.destMap, m1ToR); addInbound(fx.destMap, m2ToR);
		addInbound(fx.destMap, rToM1); addInbound(fx.destMap, rToM2);
		return fx;
	}

	/** Sets the real BubblePopper static flags to isolate exactly the unzip-true-bubble path:
	 * no direct merging, no indirect error pruning, unzip on, graph validation on. */
	static void configureUnzip(){
		BubblePopper.verbose=false;
		BubblePopper.popDirect=false;
		BubblePopper.popIndirect=false;
		BubblePopper.unzipBubbles=true;
		BubblePopper.validateGraph=true;
	}

	static class IsolatedBubble {
		Contig L, M1, M2, R;
		final ArrayList<Contig> allContigs=new ArrayList<Contig>();
		final HashMap<Integer, ArrayList<Edge>> destMap=new HashMap<Integer, ArrayList<Edge>>();
		private int idCounter=0;
		int nextId(){return idCounter++;}
	}

	static void addInbound(HashMap<Integer, ArrayList<Edge>> destMap, Edge e){
		ArrayList<Edge> list=destMap.get(e.destination);
		if(list==null){list=new ArrayList<Edge>(); destMap.put(e.destination, list);}
		list.add(e);
	}

	static int countDistinctInboundSources(HashMap<Integer, ArrayList<Edge>> destMap, int id, boolean destRight){
		ArrayList<Edge> inbound=destMap.get(id);
		if(inbound==null){return 0;}
		java.util.HashSet<Integer> sources=new java.util.HashSet<Integer>();
		for(Edge e : inbound){if(e.destRight()==destRight){sources.add(e.origin);}}
		return sources.size();
	}

	static int orient(boolean sourceRight, boolean destRight){
		return (sourceRight?1:0)|(destRight?2:0);
	}

	static byte[] reverseComplement(byte[] bases){
		byte[] out=new byte[bases.length];
		for(int i=0; i<bases.length; i++){
			byte b=bases[bases.length-1-i];
			out[i]=(b=='A' ? (byte)'T' : b=='T' ? (byte)'A' : b=='C' ? (byte)'G' : b=='G' ? (byte)'C' : b);
		}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Oracle arithmetic        ----------------*/
	/*--------------------------------------------------------------*/

	/** Splices left+arm+right into the single haplotype-linear product, trimming k-1 off both
	 * ends of arm -- the same construction BubblePopper.makePath performs via
	 * left+entry.bases+mid[kbig:]+exit.bases+right[kbig:] (equivalent here since both edges are
	 * single-base and equal the overlap byte already implied by the k-1 trim). Independent of
	 * production: used as the expected value, not derived from it. */
	static String spliceHaplotype(Contig left, Contig arm, Contig right, int k){
		String armStr=new String(arm.bases);
		String armMiddle=armStr.substring(k-1, armStr.length()-(k-1));
		return new String(left.bases)+armMiddle+new String(right.bases);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Check helpers       ----------------*/
	/*--------------------------------------------------------------*/

	static boolean check(String label, Object expected, Object actual){
		boolean ok=expected==null ? actual==null : expected.equals(actual);
		if(!ok){System.err.println("FAIL: "+label+" -- expected="+expected+", actual="+actual);}
		return ok;
	}

	static boolean check(String label, double expected, double actual, double tolerance){
		boolean ok=Math.abs(expected-actual)<=tolerance;
		if(!ok){System.err.println("FAIL: "+label+" -- expected="+expected+", actual="+actual);}
		return ok;
	}
}
