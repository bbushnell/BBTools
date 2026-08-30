package prok;

import consensus.BaseGraph;
import ml.CellNet;
import ml.CellNetParser;
import shared.KillSwitch;

/**
 * Inference wiring for the ncRNA boundary-precision NN (C3, Noire's spec
 * plans/c3_ncrnaboundaryscorer_spec.md) -- family-agnostic analog of
 * TrnaBoundaryScorer, wired into NcrnaScavenger.trimToAlignmentExtent instead of
 * TrnaCaller's tRNA-specific path.
 *
 * <p>SIMPLER than TrnaBoundaryScorer in two structural ways, both because ncRNA
 * boundary refinement runs only on already-accepted, rare ncRNA calls (not
 * tRNA's per-genome hot loop):
 * <ul>
 * <li>No stem feature (10 dims, not 11) -- the acceptor-stem palindrome is
 *     tRNA-specific; see NcrnaBoundaryVectorGen's javadoc for the full 10-dim
 *     feature list this class must reproduce bit-for-bit.</li>
 * <li>No tip-adjustment approximation -- ani and fuzziness are BOTH brute-force
 *     recomputed per candidate (TrnaBoundaryFeatures.aniFeature /
 *     tipFuzzinessFeature directly), matching NcrnaBoundaryVectorGen's training
	 *     vectors exactly (zero train/inference gap). Up to 12 realignments per
 *     locus (6 START + 6 STOP candidates) is affordable at this call
 *     frequency; TrnaBoundaryScorer's one-alignment-per-locus tip-adjustment
 *     machinery exists specifically to avoid that cost in tRNA's hot path and
 *     has no equivalent need here.</li>
 * </ul>
 *
 * <p>Also unlike TrnaBoundaryScorer: no cross-boundary-enrichment option (the
 * ncRNA feature vector is always exactly 10 dims, never 14) and no dedicated
 * per-boundary START/STOP dispatch net-cloning story here -- callers own net
 * lifecycle (this class only loads and scores).
 *
	 * <p>Candidate offsets are supplied by the owning NcrnaFamily. The same arrays
	 * are consumed by NcrnaBoundaryVectorGen, so each family retains exact
	 * train/inference parity while allowing different endpoint-error geometry.
 *
 * @author G11
 */
public class NcrnaBoundaryScorer {

	/** The ncRNA boundary feature vector is always exactly this many dims (ani,
	 * prof0-2, isStop, fuzz0-2, lengthRatio, contigGC) -- no stem, no optional
	 * cross-boundary-enrichment variant, unlike TrnaBoundaryScorer's 11-or-14. */
	public static final int NUM_FEATURES=10;

	/** Loads a boundary net and fails loud if its declared input dimension isn't
	 * exactly NUM_FEATURES (Citan's explicit requirement, fail-fast on dimension
	 * mismatch): a net trained at a different dim was not trained by
	 * NcrnaBoundaryVectorGen's current 10-dim format, or is corrupt, and must
	 * halt the whole process rather than silently score garbage. */
	public static CellNet load(String path){
		final CellNet net=CellNetParser.load(path);
		final int dims=net.numInputs();
		assert(dims==NUM_FEATURES) : KillSwitch.assertDie("Ncrna boundary net '"+path+"' has "
			+dims+" inputs; the ncRNA boundary feature vector is always "+NUM_FEATURES
			+" dims (ani, prof0-2, isStop, fuzz0-2, lengthRatio, contigGC -- no stem, no "
			+"cross-boundary enrichment) -- this net was not trained by NcrnaBoundaryVectorGen's "
			+"current format, or is corrupt.");
		return net;
	}

	/** Confidence-margin cutoff, mirrors TrnaBoundaryScorer's MARGIN_THRESHOLD_START/STOP
	 * (see there for the full rationale). Both default 0f -- unlike tRNA's swept optimum
	 * (START=0.10/STOP=0), no margin sweep has run for any ncRNA family yet, so "always
	 * move if anything beats current" (margin=0f, fully backward-compatible with "no
	 * margin gate at all") is the correct un-swept default. */
	static float MARGIN_THRESHOLD_START=0f;
	static float MARGIN_THRESHOLD_STOP=0f;

	/**
	 * Scores one candidate boundary directly. isStop selects which end of [s,e] is being
	 * evaluated (matches NcrnaBoundaryVectorGen's varyStart convention: when
	 * isStop==false, s is the candidate and e is held fixed; when isStop==true, e is the
	 * candidate and s is fixed).
	 * @param modelConsensus The chosen model's library consensus sequence -- ani is
	 *   recomputed from THIS candidate's own [s,e] span against it on every call
	 *   (brute-force TrnaBoundaryFeatures.aniFeature, not tipAdjustAni), matching
	 *   NcrnaBoundaryVectorGen's training-vector generation exactly.
	 * @param model The chosen model's HBM BaseGraph, for the fuzziness feature (same
	 *   model as modelConsensus, parallel-indexed). Null falls back to PENDING_DORI,
	 *   matching TrnaBoundaryFeatures.
	 * @param table Real enrichment table for THIS boundary type (start-table when
	 *   isStop==false, stop-table when isStop==true -- these differ, see BoundaryType).
	 * @param insideCount @param outsideCount The window composition THIS table was built
	 *   with (must match training exactly).
	 * @param contigGC Full-source-contig GC fraction -- a genuine per-locus constant.
	 * @param meanLen The family's mean/median model length (an NcrnaFamily-level
	 *   constant, NOT TrnaBoundaryFeatures.lengthRatioFeature's hardcoded tRNA
	 *   MEAN_TRNA_LEN=76f) -- families differ hugely in length (rnasep ~380bp vs
	 *   srp_small ~95bp), so this must be threaded in per family, never shared.
	 */
	public static float score(CellNet net, byte[] window, int s, int e, boolean isStop,
			byte[] modelConsensus, BaseGraph model, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, float contigGC, float meanLen){
		final boolean useStart=!isStop;
		final int boundaryPos=(isStop ? e : s);
		final TrnaBoundaryFeatures.BoundaryType type=(isStop
			? TrnaBoundaryFeatures.BoundaryType.STOP : TrnaBoundaryFeatures.BoundaryType.START);
		final float[] prof=TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos, type, insideCount, outsideCount, table);
		final byte[] candSeq=java.util.Arrays.copyOfRange(window, s, e+1);
		final float ani=TrnaBoundaryFeatures.aniFeature(candSeq, modelConsensus);
		final float[] fuzz=TrnaBoundaryFeatures.tipFuzzinessFeature(candSeq, model, useStart);
		final float lengthRatio=(e-s+1)/meanLen;
		final float[] in=buildInput(ani, prof, isStop, fuzz, lengthRatio, contigGC);
		net.applyInput(in);
		return net.feedForward();
	}

	/** Builds the fixed 10-dim input vector. Exact order per NcrnaBoundaryVectorGen:143-144
	 * (the format that trained the nets -- must match bit-for-bit): ani, prof0-2, isStop,
	 * fuzz0-2, lengthRatio, contigGC. No stem, no cross-boundary-enrichment slot. */
	private static float[] buildInput(float ani, float[] prof, boolean isStop, float[] fuzz,
			float lengthRatio, float contigGC){
		return new float[]{ani, prof[0], prof[1], prof[2], (isStop ? 1f : 0f), fuzz[0], fuzz[1], fuzz[2], lengthRatio, contigGC};
	}

	/**
	 * Shipping inference entry point. Mirrors TrnaBoundaryScorer.refineBoundaries'
	 * sequential-refinement design exactly (see there for the full rationale: score both
	 * boundaries at their current offset=0 position first; the one with LOWER net
	 * confidence there is refined FIRST; the SECOND boundary's sweep then uses the
	 * FIRST boundary's NEW position, so its features -- all recomputed fresh per
	 * candidate regardless here -- reflect reality, not a stale original).
	 *
	 * <p>SIMPLER than the tRNA version: no AlignmentStats/alignToModelFrame pre-pass and
	 * no tip-adjustment bookkeeping threaded through -- every candidate's score() call
	 * independently realigns for ani and fuzziness (see score's javadoc for why this is
	 * affordable here).
	 * @return {bestStartOffset, bestStopOffset}, each drawn from the corresponding
	 *   family-configured candidate array. Unlike
	 *   TrnaBoundaryScorer, never returns null -- there is no separate "confident base
	 *   alignment" pre-check here (no base alignment is computed at all); a null/
	 *   no-confident-model situation is the caller's responsibility to gate before
	 *   calling, matching NcrnaScavenger's existing model-selection guards elsewhere.
	 */
	public static int[] refineBoundaries(CellNet startNet, CellNet stopNet, byte[] window, int s, int e,
			byte[] modelConsensus, BaseGraph model,
			TrnaBoundaryFeatures.NinemerTable startTable, TrnaBoundaryFeatures.NinemerTable stopTable,
			int startInside, int startOutside, int stopInside, int stopOutside, float contigGC, float meanLen){
		return refineBoundaries(startNet, stopNet, window, s, e, modelConsensus, model, startTable, stopTable,
			startInside, startOutside, stopInside, stopOutside, contigGC, meanLen,
			NcrnaFamily.LEGACY_START_OFFSETS, NcrnaFamily.LEGACY_STOP_OFFSETS);
	}

	public static int[] refineBoundaries(CellNet startNet, CellNet stopNet, byte[] window, int s, int e,
			byte[] modelConsensus, BaseGraph model,
			TrnaBoundaryFeatures.NinemerTable startTable, TrnaBoundaryFeatures.NinemerTable stopTable,
			int startInside, int startOutside, int stopInside, int stopOutside, float contigGC, float meanLen,
			int[] startOffsets, int[] stopOffsets){
		final float startConf=score(startNet, window, s, e, false, modelConsensus, model,
			startTable, startInside, startOutside, contigGC, meanLen);
		final float stopConf=score(stopNet, window, s, e, true, modelConsensus, model,
			stopTable, stopInside, stopOutside, contigGC, meanLen);

		final int bestStartOffset, bestStopOffset;
		if(startConf<=stopConf){//start is the worse (or tied) boundary -- refine it first
			bestStartOffset=applyMargin(bestOffset(startNet, window, s, e, false, modelConsensus, model,
				startTable, startInside, startOutside, contigGC, meanLen, startOffsets), MARGIN_THRESHOLD_START);
			bestStopOffset=applyMargin(bestOffset(stopNet, window, s+bestStartOffset, e, true, modelConsensus, model,
				stopTable, stopInside, stopOutside, contigGC, meanLen, stopOffsets), MARGIN_THRESHOLD_STOP);
		}else{
			bestStopOffset=applyMargin(bestOffset(stopNet, window, s, e, true, modelConsensus, model,
				stopTable, stopInside, stopOutside, contigGC, meanLen, stopOffsets), MARGIN_THRESHOLD_STOP);
			bestStartOffset=applyMargin(bestOffset(startNet, window, s, e+bestStopOffset, false, modelConsensus, model,
				startTable, startInside, startOutside, contigGC, meanLen, startOffsets), MARGIN_THRESHOLD_START);
		}
		return new int[]{bestStartOffset, bestStopOffset};
	}

	/** Applies a MARGIN_THRESHOLD gate to one boundary's sweep result {bestOffset,
	 * bestScore, zeroScore} -- see TrnaBoundaryScorer.applyMargin's javadoc for the full
	 * rationale (identical logic here). margin=0f (both defaults) reproduces "always move
	 * if anything beats current" exactly. */
	private static int applyMargin(float[] sweepResult, float margin){
		final int bestOffset=(int)sweepResult[0];
		final float bestScore=sweepResult[1], zeroScore=sweepResult[2];
		return (bestScore-zeroScore>margin) ? bestOffset : 0;
	}

	/**
	 * One boundary's candidate sweep, given the CURRENT (possibly already-refined)
	 * position of the OTHER boundary (s0 for START sweeps, e0 for STOP sweeps --
	 * whichever is held fixed here). See refineBoundaries' javadoc for the sequential-
	 * refinement design.
	 * @return {bestOffset, bestScore, zeroScore} -- zeroScore is this same sweep's own
	 *   score AT offset=0 (0 is always one of the swept candidates for a valid base
	 *   candidate), the correct "don't move" baseline for applyMargin's comparison.
	 */
	private static float[] bestOffset(CellNet net, byte[] window, int s0, int e0, boolean isStop,
			byte[] modelConsensus, BaseGraph model, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, float contigGC, float meanLen, int[] offsets){
		int bestOff=0; float bestScore=-Float.MAX_VALUE, zeroScore=Float.NaN;
		for(int offset : offsets){
			final int s=s0+(isStop ? 0 : offset), e=e0+(isStop ? offset : 0);
			if(s<0 || e>=window.length || e-s<15){continue;}
			final float sc=score(net, window, s, e, isStop, modelConsensus, model, table, insideCount, outsideCount, contigGC, meanLen);
			if(offset==0){zeroScore=sc;}
			if(sc>bestScore){bestScore=sc; bestOff=offset;}
		}
		//Kept identical to TrnaBoundaryScorer.bestOffsetTipAdjusted's own assert (plain, not
		//assertDie -- this runs deep inside NcrnaScavenger's normal per-contig call path, same
		//risk profile as the tRNA equivalent, not a producer/consumer pipeline where a silent
		//worker-thread death would hang a separate consumer).
		assert(!Float.isNaN(zeroScore)) : "offset=0 is always in-bounds for a valid base candidate -- "
			+"scoreOffset should never return NaN for the unshifted span (window.length="+window.length
			+", s0="+s0+", e0="+e0+", isStop="+isStop+")";
		return new float[]{bestOff, bestScore, zeroScore};
	}
}
