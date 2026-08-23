package prok;

import consensus.BaseGraph;
import idaligner.AlignmentStats;
import ml.CellNet;
import ml.CellNetParser;
import stream.Read;

/**
 * Inference wiring stub for the tRNA boundary-precision NN (Brian's idea #3,
 * trna repo plans/boundary_precision_ideas.md). NOT wired into TrnaCaller's live
 * calling path -- that file is under a separate hold (#1a). This class shows how
 * a trained net WOULD be queried once one exists: load it, build the same v3
 * feature vector (stem, ani, enrichmentProfile x3, isStop, tipFuzziness x3,
 * lengthRatio, contigGC -- 11 values, or 14 with the cross-boundary enrichment
 * x3 appended when TrnaBoundaryVectorGen.INCLUDE_CROSS_ENRICHMENT is on, queued
 * directive #2) TrnaBoundaryVectorGen emits per candidate, feed it forward.
 * Uses enrichmentProfile's INFERENCE overload (no leave-one-out) -- a real
 * query candidate has no ground truth and is generally not even in the
 * training corpus.
 *
 * <p>LEAKAGE FIX (Brian, Aug 21 2026): the v1 version of this class took ani as a
 * per-locus-constant PARAMETER and reused it across all 5 offset queries in
 * bestOffset -- the EXACT SAME leak Brian caught in TrnaBoundaryVectorGen (both
 * files inherited the same Aug 19 design assumption).
 *
 * <p>ONE-ALIGNMENT-PER-LOCUS REDESIGN (Brian, Aug 21 2026, same day): the first
 * leakage fix corrected ani/fuzziness to be recomputed per candidate, but did so
 * by REALIGNING for every one of the 10 candidates per locus (5 start offsets +
 * 5 stop offsets) -- correct, but 10 ScrabbleAligner calls per locus is not
 * acceptable in CallGenes' hot path (~0.05s/genome budget). The shipping design
 * instead aligns the current best-guess candidate ONCE (refineBoundaries below),
 * then approximates ani via O(1) tip-adjustment (TrnaBoundaryFeatures.tipAdjustAni)
 * -- a handful of byte comparisons plus one division per shift, not a realignment.
 * The ~1-2% ani approximation error is a deliberate, accepted train/inference gap
 * (TrnaBoundaryVectorGen's training vectors still brute-force realign for ani,
 * correct, already generated) -- see tipAdjustAni's javadoc for the full rationale.
 *
 * <p>FUZZINESS IS HYBRID, NOT O(1) (Noire, Aug 22 2026): the pure O(1) tip-adjusted
 * fuzziness (TrnaBoundaryFeatures.fuzzinessFromAlignPos) was validated BAD
 * (TipAdjustCheck.java: 37% mean error, up to 100%) -- its ungapped-shift
 * assumption breaks whenever the true realignment absorbs a shift as an indel,
 * common right at a boundary. bestOffsetTipAdjusted below instead does a real
 * cheap realignment per candidate for fuzziness specifically (tipFuzzinessFeature)
 * while keeping ani tip-adjusted -- "start with correct," Noire's call, with an
 * explicit fallback if Part 3's grading shows this too slow for the genome budget.
 *
 * <p>SEARCH RANGE IS ASYMMETRIC PER BOUNDARY (Brian, Aug 22 2026, from the
 * boundary-offset histogram on the 203-genome bench AND the full 2.01M-record
 * corpus): START is TrnaBoundaryVectorGen.START_OFFSETS (-3..+2), STOP is
 * STOP_OFFSETS (-1..+4) -- not a uniform +-2. Shared fields with
 * TrnaBoundaryVectorGen so inference can never query a candidate set training
 * didn't score.
 *
 * <p>DESIGN NOTE -- ani now measures identity-to-PIVOT, not identity-to-cleaned-
 * consensus: the base alignment (for both ani and fuzziness) is against
 * model.original (the BaseGraph's raw pivot sequence), not library[bestModel]
 * (the majority-vote-cleaned consensus TrnaBoundaryVectorGen's ani uses). This
 * is necessary for ONE alignment to serve both features (fuzziness's countSum
 * data lives only on the BaseGraph, keyed to model.original's coordinates) --
 * flagged to Noire/Brian as a deliberate choice, not a rediscovery of the
 * original leak: pivot and cleaned-consensus sequences for the same model are
 * typically near-identical (the pivot IS a real cluster member; the consensus is
 * a majority vote across the same cluster), so the semantic drift should be
 * small, but it is a genuine definitional difference from the training vectors'
 * ani and worth Brian's explicit sign-off before this ships.
 *
 * Intended usage pattern (once Dori has trained a real net and Noire/Brian decide
 * to integrate it), mirroring Orf.calcOrfScoreHybrid's existing net-query style:
 *   CellNet startNet = TrnaBoundaryScorer.load(path5);
 *   CellNet stopNet = TrnaBoundaryScorer.load(path3);  // same object as startNet if one shared net
 *   int[] offsets = TrnaBoundaryScorer.refineBoundaries(startNet, stopNet, window, s, e, model,
 *       startTable, stopTable, startInside, startOutside, stopInside, stopOutside, contigGC);
 *   // offsets == null: no confident model alignment, don't refine.
 *   // offsets[0]: best start offset (one of START_OFFSETS, -3..+2) by net confidence.
 *   // offsets[1]: best stop offset (one of STOP_OFFSETS, -1..+4) by net confidence.
 *
 * <p>TWO-NET DISPATCH (Noire's spec, Aug 22 2026, queued directive #2/#3 infrastructure):
 * refineBoundaries takes SEPARATE nets for the start and stop sweeps -- startNet scores every
 * START_OFFSETS candidate, stopNet scores every STOP_OFFSETS candidate, dispatched by isStop
 * exactly where the sequential refinement already knows which boundary it's evaluating. A
 * caller with one shared net (the current shipped configuration) passes the same CellNet object
 * for both parameters -- behaviorally identical to the old single-net signature. Built to run
 * Noire's zero-new-training hybrid experiment (ml.Trainer's net for 5' moves, RT-sigmoid's net
 * for 3' moves, both EXISTING 11-dim v3 nets, no retraining) as this plumbing's own verification,
 * before committing to training dedicated per-boundary (10-dim, isStop-dropped) nets.
 *
 * @author Neptune
 */
public class TrnaBoundaryScorer {

	public static CellNet load(String path){
		return CellNetParser.load(path);
	}

	/** Confidence-margin cutoff (Brian, Aug 22 2026, queued directive #4): only move a
	 * boundary if the best candidate beats staying at the current (offset=0) position by
	 * more than this margin -- "always pick highest" may move a boundary on a razor-thin,
	 * noise-level margin in palindrome-confusion cases where moving is actively harmful.
	 * DEFAULT 0f preserves the original "always move if anything beats current" behavior
	 * exactly (a strictly-positive comparison already requires beating, not tying,
	 * current) -- fully backward compatible until explicitly swept. Toggle via
	 * CallGenes' trnaboundarymargin= flag (not yet added -- add when running the sweep;
	 * this field exists now so the sweep is a flag change, not a code change, per Brian's
	 * "make it toggleable" instruction while the STOP_OFFSETS-fix retrain is running). */
	/** Split into independent START/STOP thresholds (Brian, Aug 22 2026, queued directive #1),
	 * swept per-boundary on the rt_a5_d0.95 champion (14-dim, 95.08% both-exact at 0/0):
	 * START has a real, measured peak at margin=0.10 (95.15% both-exact, a broad plateau
	 * 0.05-0.20 all within 0.01pp, declining beyond that) -- some fraction of START moves are
	 * confident-but-wrong, and a small margin filters them out. STOP is the OPPOSITE: strictly,
	 * monotonically declining from the very first nonzero value tested (94.52% at margin=0.05,
	 * falling to 90.65% at 0.50) -- STOP's moves are reliably correct, so any threshold only
	 * blocks legitimate improving moves. Defaults set to the measured optimum for each boundary
	 * independently (both values were measured with the OTHER boundary's margin held at 0, so
	 * this exact combination -- 0.10/0 -- was directly tested, not composed from an
	 * independence assumption). Full sweep: plans/margin_sweep_perboundary_20260822.md. */
	static float MARGIN_THRESHOLD_START=0.10f;
	static float MARGIN_THRESHOLD_STOP=0f;

	/**
	 * Scores one candidate boundary directly (brute-force, no tip-adjustment) --
	 * kept for direct single-candidate queries and for cross-checking
	 * refineBoundaries' tip-adjusted approximation against a real realignment.
	 * isStop selects which end of [s,e] is being evaluated (matches
	 * TrnaBoundaryVectorGen's varyStart convention: when isStop==false, s is the
	 * candidate and e is held fixed; when isStop==true, e is the candidate and s
	 * is fixed).
	 * @param modelConsensus The chosen model's consensus sequence (library entry)
	 *   for ani -- ani is recomputed from THIS candidate's own [s,e] span against
	 *   it, never reused from another span.
	 * @param model The chosen model's HBM BaseGraph (same model as modelConsensus,
	 *   parallel-indexed), for the fuzziness feature. Null falls back to
	 *   PENDING_DORI (whole capability off), matching TrnaBoundaryFeatures.
	 * @param table Real enrichment table for THIS boundary type (null falls back to
	 *   PENDING_DORI). Caller passes the start-table when isStop==false, the
	 *   stop-table when isStop==true -- these differ (see BoundaryType).
	 * @param insideCount @param outsideCount The window composition THIS table was
	 *   built with (must match training exactly).
	 * @param contigGC Full-source-contig GC fraction -- a genuine per-locus
	 *   constant, unlike ani/fuzziness.
	 */
	public static float score(CellNet net, byte[] window, int s, int e, boolean isStop,
			byte[] modelConsensus, BaseGraph model, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, float contigGC,
			TrnaBoundaryFeatures.NinemerTable otherTable, int otherInside, int otherOutside){
		final boolean useStart=!isStop;
		float stem=TrnaBoundaryFeatures.stemFeature(window, s, e);
		int boundaryPos=(isStop ? e : s);
		final TrnaBoundaryFeatures.BoundaryType type=(isStop
			? TrnaBoundaryFeatures.BoundaryType.STOP : TrnaBoundaryFeatures.BoundaryType.START);
		float[] prof=TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos, type, insideCount, outsideCount, table);
		final byte[] candSeq=java.util.Arrays.copyOfRange(window, s, e+1);
		float ani=TrnaBoundaryFeatures.aniFeature(candSeq, modelConsensus);
		float[] fuzz=TrnaBoundaryFeatures.tipFuzzinessFeature(candSeq, model, useStart);
		float lengthRatio=TrnaBoundaryFeatures.lengthRatioFeature(s, e);
		final float[] in=buildInput(stem, ani, prof, isStop, fuzz, lengthRatio, contigGC,
			window, (isStop ? s : e), otherTable, otherInside, otherOutside);
		net.applyInput(in);
		return net.feedForward();
	}

	/** Builds the input vector, appending the cross-boundary enrichment profile (queued
	 * directive #2) when TrnaBoundaryVectorGen.INCLUDE_CROSS_ENRICHMENT is on -- 11 or 14
	 * values, matching whatever dims the loaded net was trained with exactly. otherPos is
	 * the OPPOSITE boundary's CURRENT position (its own s or e, whichever this call didn't
	 * vary) -- computed fresh by every caller from the CURRENT s0/e0, so if the opposite
	 * boundary already moved (sequential refinement, queued directive #3), this reads its
	 * NEW position, not a stale one. Inference overload of enrichmentProfile (no ground
	 * truth, no leave-one-out) -- matches how THIS boundary's own profile is scored too. */
	private static float[] buildInput(float stem, float ani, float[] prof, boolean isStop, float[] fuzz,
			float lengthRatio, float contigGC, byte[] window, int otherPos,
			TrnaBoundaryFeatures.NinemerTable otherTable, int otherInside, int otherOutside){
		if(!TrnaBoundaryVectorGen.INCLUDE_CROSS_ENRICHMENT){
			return new float[]{stem, ani, prof[0], prof[1], prof[2], (isStop ? 1f : 0f), fuzz[0], fuzz[1], fuzz[2], lengthRatio, contigGC};
		}
		final TrnaBoundaryFeatures.BoundaryType otherType=(isStop
			? TrnaBoundaryFeatures.BoundaryType.START : TrnaBoundaryFeatures.BoundaryType.STOP);
		final float[] otherProf=TrnaBoundaryFeatures.enrichmentProfile(window, otherPos, otherType, otherInside, otherOutside, otherTable);
		return new float[]{stem, ani, prof[0], prof[1], prof[2], (isStop ? 1f : 0f), fuzz[0], fuzz[1], fuzz[2], lengthRatio, contigGC,
			otherProf[0], otherProf[1], otherProf[2]};
	}

	/**
	 * Shipping inference entry point (Brian, Aug 21 2026): ONE real alignment of
	 * the given best-guess candidate [s,e] against model, then refines both the
	 * start and stop boundaries.
	 *
	 * <p>SEQUENTIAL, NOT PARALLEL (Brian, Aug 22 2026): the two boundaries are not
	 * independent -- moving the stop changes the candidate's length, its acceptor-
	 * stem partner positions, and its alignment frame, all of which feed the
	 * start's own score. Refining both boundaries against the ORIGINAL unshifted
	 * opposite position (the naive parallel approach) can move the "wrong"
	 * boundary because its features were read at a now-stale opposite position.
	 * Fix: score both boundaries at their CURRENT (offset=0) position first; the
	 * one with LOWER net confidence there is refined FIRST; the SECOND boundary's
	 * sweep then uses the FIRST boundary's NEW position, so its stem/lengthRatio/
	 * fuzziness features (all recomputed fresh per candidate regardless) reflect
	 * reality, not the stale original. This is a 1D refiner applied twice with
	 * state updated between calls -- architecturally correct, not a hack; the net
	 * itself is unchanged, only the order and the feature recompute between calls.
	 * (Ani's O(1) tip-adjustment still anchors to the ONE original base alignment
	 * regardless of which boundary moved first -- an accepted secondary
	 * approximation on top of the primary one documented above, per Noire.)
	 * @param startNet Net queried for every START_OFFSETS candidate. Pass the same object as
	 *   stopNet for a single shared net (behaviorally identical to the pre-two-net signature).
	 * @param stopNet Net queried for every STOP_OFFSETS candidate.
	 * @return {bestStartOffset, bestStopOffset} (bestStartOffset one of
	 *   TrnaBoundaryVectorGen.START_OFFSETS, bestStopOffset one of STOP_OFFSETS),
	 *   or null if no confident model alignment exists for the base candidate
	 *   (mirrors TrnaBoundaryFeatures.alignToModelFrame's 0.2-identity floor).
	 */
	public static int[] refineBoundaries(CellNet startNet, CellNet stopNet, byte[] window, int s, int e, BaseGraph model,
			TrnaBoundaryFeatures.NinemerTable startTable, TrnaBoundaryFeatures.NinemerTable stopTable,
			int startInside, int startOutside, int stopInside, int stopOutside, float contigGC){
		if(model==null){return null;}
		final byte[] candidate=java.util.Arrays.copyOfRange(window, s, e+1);
		final AlignmentStats stats=new AlignmentStats(true);
		final Read r=TrnaBoundaryFeatures.alignToModelFrame(candidate, model, stats);
		if(r==null){return null;}
		final int baseMatches=stats.matches;
		final int baseAlnLen=stats.matches+stats.subs+stats.ins+stats.dels;
		final byte[] modelBases=model.original;

		final float startConf=scoreOffset(startNet, window, s, s, e, false, model, modelBases,
			startTable, startInside, startOutside, contigGC, r.start, baseMatches, baseAlnLen, 0,
			stopTable, stopInside, stopOutside);
		final float stopConf=scoreOffset(stopNet, window, e, s, e, true, model, modelBases,
			stopTable, stopInside, stopOutside, contigGC, r.stop, baseMatches, baseAlnLen, 0,
			startTable, startInside, startOutside);

		final int bestStartOffset, bestStopOffset;
		if(startConf<=stopConf){//start is the worse (or tied) boundary -- refine it first
			bestStartOffset=applyMargin(bestOffsetTipAdjusted(startNet, window, s, e, false, model, modelBases,
				startTable, startInside, startOutside, contigGC, r.start, baseMatches, baseAlnLen,
				stopTable, stopInside, stopOutside), MARGIN_THRESHOLD_START);
			bestStopOffset=applyMargin(bestOffsetTipAdjusted(stopNet, window, s+bestStartOffset, e, true, model, modelBases,
				stopTable, stopInside, stopOutside, contigGC, r.stop, baseMatches, baseAlnLen,
				startTable, startInside, startOutside), MARGIN_THRESHOLD_STOP);
		}else{
			bestStopOffset=applyMargin(bestOffsetTipAdjusted(stopNet, window, s, e, true, model, modelBases,
				stopTable, stopInside, stopOutside, contigGC, r.stop, baseMatches, baseAlnLen,
				startTable, startInside, startOutside), MARGIN_THRESHOLD_STOP);
			bestStartOffset=applyMargin(bestOffsetTipAdjusted(startNet, window, s, e+bestStopOffset, false, model, modelBases,
				startTable, startInside, startOutside, contigGC, r.start, baseMatches, baseAlnLen,
				stopTable, stopInside, stopOutside), MARGIN_THRESHOLD_START);
		}
		return new int[]{bestStartOffset, bestStopOffset};
	}

	/** Applies a MARGIN_THRESHOLD gate (queued directive #4, split per-boundary per directive
	 * #1) to one boundary's sweep result {bestOffset, bestScore, zeroScore}: only actually
	 * move the boundary if the best candidate beat staying at the current (offset=0) position
	 * by more than the given margin. zeroScore is read from the SAME sweep, evaluated against
	 * whatever the OTHER boundary's position was at sweep time (already-refined if this is the
	 * second boundary) -- not a stale pre-refinement value, so the comparison is always fair to
	 * the boundary's actual current reality. margin=0f (default for both MARGIN_THRESHOLD_START
	 * and MARGIN_THRESHOLD_STOP) reproduces the original "always move if anything beats
	 * current" behavior exactly.
	 * @param margin The caller passes MARGIN_THRESHOLD_START or MARGIN_THRESHOLD_STOP,
	 *   matching which boundary this sweep result is for -- never the other one. */
	private static int applyMargin(float[] sweepResult, float margin){
		final int bestOffset=(int)sweepResult[0];
		final float bestScore=sweepResult[1], zeroScore=sweepResult[2];
		return (bestScore-zeroScore>margin) ? bestOffset : 0;
	}

	/** One boundary's candidate sweep, given the ONE base alignment's stats and the
	 * CURRENT (possibly already-refined) position of the OTHER boundary (s0 for
	 * START sweeps, e0 for STOP sweeps -- whichever is held fixed here). See
	 * refineBoundaries' javadoc for the sequential-refinement design and the
	 * hybrid-fuzziness rationale.
	 * @return {bestOffset, bestScore, zeroScore} -- zeroScore is this same sweep's own
	 *   score AT offset=0 (0 is always one of the swept candidates), the correct
	 *   "don't move" baseline for applyMargin's comparison. */
	private static float[] bestOffsetTipAdjusted(CellNet net, byte[] window, int s0, int e0, boolean isStop,
			BaseGraph model, byte[] modelBases, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, float contigGC,
			int baseAlignPos, int baseMatches, int baseAlnLen,
			TrnaBoundaryFeatures.NinemerTable otherTable, int otherInside, int otherOutside){
		final int boundaryPos=(isStop ? e0 : s0);
		final int[] offsets=(isStop ? TrnaBoundaryVectorGen.STOP_OFFSETS : TrnaBoundaryVectorGen.START_OFFSETS);
		int bestOffset=0; float bestScore=-1, zeroScore=Float.NaN;
		for(int offset : offsets){
			final float sc=scoreOffset(net, window, boundaryPos, s0, e0, isStop, model, modelBases,
				table, insideCount, outsideCount, contigGC, baseAlignPos, baseMatches, baseAlnLen, offset,
				otherTable, otherInside, otherOutside);
			if(Float.isNaN(sc)){continue;}
			if(offset==0){zeroScore=sc;}
			if(sc>bestScore){bestScore=sc; bestOffset=offset;}
		}
		assert(!Float.isNaN(zeroScore)) : "offset=0 is always in-bounds for a valid base candidate -- "
			+"scoreOffset should never return NaN for the unshifted span (window.length="+window.length
			+", s0="+s0+", e0="+e0+", isStop="+isStop+")";
		return new float[]{bestOffset, bestScore, zeroScore};
	}

	/** Scores ONE candidate offset for one boundary: builds the full 11-feature vector
	 * (ani tip-adjusted from the one base alignment; stem/enrichment/fuzziness/
	 * lengthRatio all computed fresh from the ACTUAL current [s,e] -- which is why
	 * passing an already-refined opposite-boundary s0/e0 in correctly updates them)
	 * and feeds it forward. Returns Float.NaN if this offset is out of bounds or the
	 * tip-adjusted ani is degenerate -- callers must skip, not score, a NaN.
	 * @param boundaryPos The boundary's OWN current position (s0 for START, e0 for
	 *   STOP) -- the tip-adjustment's "before any shift" reference point.
	 * @param s0 @param e0 The candidate's CURRENT full span (one of these is
	 *   boundaryPos; the other is the opposite boundary's current, possibly
	 *   already-refined, position). */
	private static float scoreOffset(CellNet net, byte[] window, int boundaryPos, int s0, int e0, boolean isStop,
			BaseGraph model, byte[] modelBases, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, float contigGC,
			int baseAlignPos, int baseMatches, int baseAlnLen, int offset,
			TrnaBoundaryFeatures.NinemerTable otherTable, int otherInside, int otherOutside){
		final int growDirection=(isStop ? 1 : -1);
		final int s=s0+(isStop ? 0 : offset), e=e0+(isStop ? offset : 0);
		if(s<0 || e>=window.length || e-s<15){return Float.NaN;}
		final float[] aniAndPos=TrnaBoundaryFeatures.tipAdjustAni(window, boundaryPos, baseAlignPos,
			growDirection, offset, modelBases, baseMatches, baseAlnLen);
		final float ani=aniAndPos[0];
		if(Float.isNaN(ani)){return Float.NaN;}//degenerate alnLen at this offset -- skip
		final float stem=TrnaBoundaryFeatures.stemFeature(window, s, e);
		final int boundaryPosC=(isStop ? e : s);
		final TrnaBoundaryFeatures.BoundaryType type=(isStop
			? TrnaBoundaryFeatures.BoundaryType.STOP : TrnaBoundaryFeatures.BoundaryType.START);
		final float[] prof=TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPosC, type, insideCount, outsideCount, table);
		final byte[] candSeq=java.util.Arrays.copyOfRange(window, s, e+1);
		final float[] fuzz=TrnaBoundaryFeatures.tipFuzzinessFeature(candSeq, model, !isStop);
		final float lengthRatio=TrnaBoundaryFeatures.lengthRatioFeature(s, e);
		//Cross-boundary enrichment (queued directive #2/#3): otherPos is whichever of s/e is
		//NOT this boundary -- the OTHER boundary's CURRENT position, read fresh here every
		//call, so a sequential-refinement second pass automatically sees the first boundary's
		//NEW position (Brian's explicit ask: full recompute, not just stem, when the opposite
		//end has moved).
		final int otherPos=(isStop ? s : e);
		final float[] in=buildInput(stem, ani, prof, isStop, fuzz, lengthRatio, contigGC,
			window, otherPos, otherTable, otherInside, otherOutside);
		net.applyInput(in);
		return net.feedForward();
	}
}
