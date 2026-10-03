package prot;

/**
 * Conditional HBM-style scorer/builder for an {@link AAGraph} protein consensus model, ported
 * directly from {@code consensus.BaseGraph}/{@code BaseNode}'s ACTIVE scoring semantics (a
 * 20-residue generalization of the DNA 4-base case; no new formula invented). Composition, not
 * modification: reads {@link AAGraph}/{@link AAGraphNode}'s existing public fields; neither class
 * is touched by this file.
 *
 * <p>Intended as a FINAL CONDITIONAL scorer, per Brian's framing: applied only to a candidate that
 * has already passed every other filter (score/identity/overlap/length), never scored
 * unconditionally. {@link #scoreIfSurvivor} enforces this at the API boundary with an explicit
 * boolean gate; this class does not itself decide what "every other filter" means.</p>
 *
 * <p><b>X_CODE / undefined-residue behavior</b> (the amino-acid analog of DNA's non-ACGT case):
 * {@link BaseNode#baseProb} in the nucleotide original returns the uniform prior 0.25 for any base
 * that is not one of the 4 real symbols, and {@code calcProbs}'s denominator sums ONLY the 4 real
 * bases -- N observations increment {@code countSum} but have no slot in {@code acgtCount}, so they
 * are invisible to the frequency computation entirely. {@link Blosum62#X_CODE} is verified from
 * source (not inferred) to be <b>21</b>, not 20 -- {@code assert(X_CODE==21)} in
 * {@code Blosum62.buildMatrix}; an earlier draft of this file said 20 in its javadoc (never in a
 * bounds check that mattered, since no valid encoded value is ever exactly 20, but wrong as
 * documentation and worth being exact about). {@link #baseProb(AAGraphNode,byte)} returns the
 * uniform prior {@code 1/REAL_RESIDUES} (0.05) for X_CODE(21) or any other out-of-range byte, with
 * the denominator summing ONLY {@code count[0..19]} -- explicitly EXCLUDING {@code count[X_CODE]}
 * (AAGraphNode, unlike BaseNode, DOES give X its own histogram slot) -- same fallback/exclusion
 * semantics as the nucleotide original, generalized to 20 real symbols instead of 4.</p>
 *
 * <p>Anchoring contract (Yoimiya's pin, 2026-09-13): the {@link AAGraph} passed to the scorer
 * must be built with {@code pad=0} directly on the family's frozen, unpadded consensus sequence
 * (the same {@code target.enc} {@code GlocalAminoLinear} already aligns against in
 * {@code FamilyCalibrationAssay}), never re-anchored onto {@link AAGraph#traverse()}'s own output.
 * {@link #buildModel} enforces this by construction.</p>
 *
 * <p><b>Terminal insertions with pad=0 (Yoimiya's review, 2026-09-13): characterized, not fixed.</b>
 * This is inherited {@code BaseGraph}/{@code AAGraph} behavior, not a new bug, but the "the whole
 * query is placed" framing above is narrower than it sounds: {@code AAGraph.add}'s own loop
 * (unchanged, not touched by this file) drops a LEADING insertion when {@code prevNode==null} --
 * with {@code pad>0} this is a rare overhang-past-the-padding edge case, but with {@code pad=0}
 * (this class's anchoring contract) a leading insertion is the NORMAL case whenever a training
 * member overhangs the consensus's N-terminus at all, so it is silently dropped from training
 * every time it occurs, not just rarely. The scorer below has the same shape: its loop guard
 * {@code rpos<graph.ref.length} means a TRAILING insertion occurring exactly at
 * {@code rpos==graph.ref.length} (past the last real column) is silently dropped from the score.
 * Neither is remediated here -- see {@code AAGraphScorerTest#testTrailingInsertionDroppedAtPad0}
 * for a concrete trace of the scoring-side drop, added for this review rather than a new terminal
 * penalty invented to paper over it.</p>
 *
 * @author Eru, 2026-09-13 (protein-HBM prototype, Brian's proposal via Yoimiya, D104-adjacent;
 *         isolated -- not wired into any production path)
 */
public final class AAGraphScorer {

	private AAGraphScorer(){}

	/** Number of REAL amino-acid symbols (encoded 0-19). Deliberately a different name from
	 * {@link AAGraphNode#NAA} (=22, {@code Blosum62.X_CODE+1}, the node's full histogram-array
	 * length including X's own slot) -- reusing that name for this different quantity is exactly
	 * the kind of confusion Yoimiya's review caught; keep them visibly distinct. */
	private static final int REAL_RESIDUES=20;

	/*--------------------------------------------------------------*/
	/*----------------           Building            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds one family's HBM: an {@link AAGraph} anchored with {@code pad=0} directly on the
	 * supplied frozen consensus, with each training member folded in via the SAME canonical
	 * detailed aligner used for final-gate scoring ({@link GlocalAminoLinear#align(byte[],byte[],boolean)},
	 * {@code recordPath=true}) -- never {@link AAAligner#alignGlocal}, so training and scoring
	 * never disagree on gap model.
	 *
	 * @param consensusEnc The family's frozen, unpadded consensus (encoded); becomes {@code
	 *        graph.pivot} verbatim since pad=0.
	 * @param trainingMembersEnc Encoded training members for this family (caller decides which
	 *        rows belong here -- see the open scientific question on leakage in the feasibility
	 *        report; this method takes no position on it).
	 * @return The trained {@link AAGraph}.
	 */
	public static AAGraph buildModel(final byte[] consensusEnc, final byte[][] trainingMembersEnc){
		final AAGraph graph=new AAGraph(consensusEnc, 0);
		for(final byte[] memberEnc : trainingMembersEnc){
			final AAAlignment aln=GlocalAminoLinear.align(memberEnc, consensusEnc, true);
			graph.add(memberEnc, aln);
		}
		return graph;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Ported node ops       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Column empirical frequency of residue {@code b}, from COUNT (not weight) -- direct port of
	 * {@code BaseNode.calcProbs}'s ACTIVE formula: {@code mult=1f/max(1,sum); acgtCount[i]*mult}
	 * (multiply by the precomputed reciprocal, not a direct divide -- preserving the source's exact
	 * float operation order/form, not just its mathematical value). BaseNode's denominator sums
	 * ONLY the 4 real-base slots (there is no N slot in {@code acgtCount} to accidentally include);
	 * the faithful 20-symbol port sums ONLY {@code count[0..19]} -- explicitly EXCLUDING
	 * {@code count[X_CODE]}, since (unlike BaseNode) {@link AAGraphNode}'s array DOES give X_CODE
	 * its own slot. Corrected from an earlier draft that divided by {@code node.countSum} (which
	 * includes X observations) -- Yoimiya's repro: pivot=A, one X-only training observation added
	 * at that column -> correct P(A)=1 (X is invisible to the real-residue frequency, exactly like
	 * an N in the DNA original), the countSum-based version gave 0.5. Returns the uniform prior
	 * {@code 1/REAL_RESIDUES} for X_CODE(21) or any other out-of-range byte.
	 */
	public static float baseProb(final AAGraphNode node, final byte b){
		if(b<0 || b>=REAL_RESIDUES){return 1f/REAL_RESIDUES;}
		final float mult=1f/Math.max(1, sumRealCounts(node));
		return node.count[b]*mult;
	}

	/** Sum of {@code count[0..19]} only -- excludes {@code count[X_CODE]} (index 21), the
	 * denominator {@link #baseProb} needs. Not cached (BaseNode caches via a separate
	 * {@code calcProbs()} pass; this recomputes per call -- a performance difference, not a
	 * semantic one, and not something this review asked to change). */
	private static int sumRealCounts(final AAGraphNode node){
		int sum=0;
		for(int i=0; i<REAL_RESIDUES; i++){sum+=node.count[i];}
		return sum;
	}

	/** {@code BaseNode.baseScore}'s ACTIVE code is a plain passthrough of {@code baseProb}; ported identically
	 * (the commented-out alternatives in {@code BaseNode.java} are abandoned experiments, not the real behavior). */
	public static float baseScore(final AAGraphNode node, final byte b){
		return baseProb(node, b);
	}

	/** Least-observed real residue at this column, by count -- direct 20-way port of {@code BaseNode.minBase}. */
	public static byte minBase(final AAGraphNode node){
		int idx=0, count=node.count[0];
		for(int i=1; i<REAL_RESIDUES; i++){if(node.count[i]<count){idx=i; count=node.count[i];}}
		return (byte)idx;
	}

	/** Most-observed real residue at this column, by count -- direct 20-way port of {@code BaseNode.maxBase}. */
	public static byte maxBase(final AAGraphNode node){
		int idx=0, count=node.count[0];
		for(int i=1; i<REAL_RESIDUES; i++){if(node.count[i]>count){idx=i; count=node.count[i];}}
		return (byte)idx;
	}

	/** Total REF+DEL observations at padded-pivot position {@code rpos} -- direct port of {@code BaseGraph.countSum}. */
	public static int countSum(final AAGraph graph, final int rpos){
		return graph.ref[rpos].countSum+graph.del[rpos].countSum;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Scoring            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Public entry point enforcing the "final conditional scorer" contract: computes a score ONLY
	 * when {@code passedAllOtherFilters} is true, returning null otherwise. This method does not
	 * itself define what "every other filter" means -- the caller owns that decision and must not
	 * call this for the bulk of candidates, only the rare survivors (Brian's framing).
	 *
	 * @param graph The family's HBM, built by {@link #buildModel} (pad=0 on the frozen consensus).
	 * @param queryEnc The query's encoded residues (same array {@code aln} was computed from).
	 * @param aln The query's alignment against {@code graph.pivot} (e.g.
	 *        {@code GlocalAminoLinear.align(queryEnc, consensusEnc, true)} -- must have {@code
	 *        recordPath=true}, i.e. a non-null {@code match}).
	 * @param passedAllOtherFilters Whether this candidate already survived every other filter.
	 * @return {@code {rawScore, relativeScore}}, or null if {@code passedAllOtherFilters} is false.
	 *         {@code relativeScore} may be NaN or infinite in a degenerate case -- see
	 *         {@link #score} and {@code calibrationreview.HbmBoundaryReview}'s all-X case.
	 *         Consumers requiring a finite score must check {@code Float.isFinite}.
	 */
	public static float[] scoreIfSurvivor(final AAGraph graph, final byte[] queryEnc,
			final AAAlignment aln, final boolean passedAllOtherFilters){
		if(!passedAllOtherFilters){return null;}
		return score(graph, queryEnc, aln);
	}

	/**
	 * Direct port of {@code BaseGraph.score(Read,boolean,boolean)}'s per-position walk to 20
	 * amino acids -- same recurrence, same {@code 'm'/'D'/'I'} branch formulas, same insertion
	 * chain-walk (a running {@code prevNode}, not a fresh lookup per column, so consecutive
	 * insertion columns correctly walk deeper into the trained {@code insEdge} chain and fall back
	 * to the last known node when the query's insertion run exceeds what training saw -- exactly
	 * {@code BaseGraph.score}'s real 'I' branch, not a simplified reconstruction). No local mode:
	 * the current q-global intent ({@link GlocalAminoLinear}'s geometry) means the whole query is
	 * placed except for the terminal-insertion dropping documented in the class javadoc, so this
	 * always computes the GLOBAL sum, matching {@code BaseGraph.score}'s {@code local=false} path
	 * exactly. Returns BOTH the raw sum and the relative (0-1, worst-to-best-possible-alignment-
	 * scaled) score together so a normalization choice does not block use of either -- exactly
	 * {@code BaseGraph.score}'s {@code relative=true} formula, computed alongside the raw value.
	 *
	 * <p>Degenerate case (Yoimiya's review): {@code BaseGraph.score} has no guard for
	 * {@code maxSum==minSum} either, so {@code (sum-minSum)/(maxSum-minSum)} is left to produce
	 * whatever IEEE 754 gives (NaN for 0/0) rather than this file inventing a threshold-policy
	 * value (an earlier draft silently returned 0f here, which was a real, undocumented deviation
	 * from a byte-for-byte port). A nonzero numerator over zero also produces infinity;
	 * consumers requiring a finite score must check {@code Float.isFinite(relative)} explicitly.
	 * PRIVATE, not public (Yoimiya's review): {@link #scoreIfSurvivor} is the only public entry
	 * point, so the "only rare survivors" gate cannot be bypassed by calling this directly (an
	 * earlier draft left this public, defeating its own claimed contract).
	 */
	private static float[] score(final AAGraph graph, final byte[] queryEnc, final AAAlignment aln){
		final byte[] match=aln.match;
		assert(match!=null) : "score() needs a path-recorded alignment (align(...,true)).";
		float sum=0, minSum=0, maxSum=0;
		int qpos=aln.qStart, rpos=aln.tStart;
		AAGraphNode prevNode=(rpos<=0 ? null : graph.ref[rpos-1]);
		for(int mpos=0; mpos<match.length && rpos<graph.ref.length; mpos++){
			final byte m=match[mpos];
			final AAGraphNode next;
			if(m=='m'){
				next=graph.ref[rpos];
				final byte b=queryEnc[qpos];
				final float baseScore=baseScore(next, b);
				final float minBaseScore=baseScore(next, minBase(next));
				final float maxBaseScore=baseScore(next, maxBase(next));
				final float stateProb=next.countSum/(float)Math.max(1, countSum(graph, rpos));
				sum+=baseScore*stateProb+0.1f*(stateProb-1f);
				minSum+=minBaseScore*stateProb+0.1f*(stateProb-1f);
				maxSum+=maxBaseScore*stateProb+0.1f*(stateProb-1f);
				qpos++; rpos++;
			}else if(m=='D'){
				next=graph.del[rpos];
				final float stateProb=next.countSum/(float)Math.max(1, countSum(graph, rpos));
				final float dScore=2*(stateProb-0.5f);
				sum+=dScore; minSum+=dScore; maxSum+=Math.max(0, dScore);
				rpos++;
			}else{//'I'
				if(prevNode==null){next=null;}
				else if(prevNode.insEdge!=null){next=prevNode.insEdge;}
				else if(prevNode.type==AAGraphNode.INS){next=prevNode;}
				else{next=null;}
				final byte b=queryEnc[qpos];
				final float baseScore, minBaseScore, maxBaseScore, stateProb;
				if(next==null){
					baseScore=-1; minBaseScore=-1; maxBaseScore=-1; stateProb=0;
				}else{
					baseScore=baseScore(next, b);
					minBaseScore=baseScore(next, minBase(next));
					maxBaseScore=baseScore(next, maxBase(next));
					stateProb=next.countSum/(float)Math.max(1, countSum(graph, rpos));
				}
				final float iScore=(stateProb-0.5f)+(0.25f*(baseScore-1));
				float minScore=(stateProb-0.5f)+(0.25f*(minBaseScore-1));
				float maxScore=(stateProb-0.5f)+(0.25f*(maxBaseScore-1));
				minScore=Math.min(iScore, minScore);
				maxScore=Math.max(0, Math.max(iScore, maxScore));
				sum+=iScore; minSum+=minScore; maxSum+=maxScore;
				qpos++;
			}
			prevNode=next;
		}
		//No degenerate-case guard, deliberately: BaseGraph.score has none either, so this produces
		//whatever IEEE 754 gives (NaN for the 0/0 case) rather than inventing a threshold-policy
		//value. See the class/method javadoc.
		final float relative=(sum-minSum)/(maxSum-minSum);
		return new float[]{sum, relative};
	}
}
