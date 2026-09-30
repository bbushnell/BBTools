package prok;

import consensus.BaseGraph;
import ml.CellNet;
import ml.CellNetParser;

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
	 * <li>V1 has no tip-adjustment approximation -- ani and fuzziness are BOTH brute-force
 *     recomputed per candidate (TrnaBoundaryFeatures.aniFeature /
 *     tipFuzzinessFeature directly), matching NcrnaBoundaryVectorGen's training
	 *     vectors exactly (zero train/inference gap). Up to 12 realignments per
 *     locus (6 START + 6 STOP candidates) is affordable at this call
 *     frequency; TrnaBoundaryScorer's one-alignment-per-locus tip-adjustment
 *     machinery exists specifically to avoid that cost in tRNA's hot path and
 *     has no equivalent need here.</li>
	 * </ul>
	 *
	 * <p>The opt-in V2 feature layout removes every per-site alignment: ANI is the
	 * accepted locus's existing QuantumAligner identity and fuzziness is a pair of
	 * three-value constants loaded for the selected model's 5' and 3' ends. The
	 * per-site enrichment profile, endpoint flag, length ratio, and contig GC retain
	 * their V1 definitions. V1 remains the default until V2 is trained and selected.
	 *
	 * <p>V3 keeps only measured per-site/per-locus inputs: enrichment profile (3),
	 * endpoint flag, length ratio, contig GC, and the accepted locus's existing
	 * QuantumAligner identity, followed by three explicit zero spare slots. It has
	 * no HBM-constant dependency. QuantumAligner consumes the full consensus query,
	 * so a proposed consensus-overlap fraction would always be 1 and is deliberately
	 * omitted rather than emitted as a fake feature.
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
	/** Shared 28-input rRNA feature contract; no network is loaded or evaluated.
	 * Candidate identity must be measured on the same inclusive raw span. */
	public static boolean rrnaEndpointFeatures(RrnaPositionalKmerTable.Table table,byte[] bases,
			int start,int stop,int consensusLength,float contigGC,float candidateIdentity,float[] scratch,float[] output){
		return RrnaEndpointVector.fill(table,bases,start,stop,consensusLength,contigGC,candidateIdentity,scratch,output);
	}

	/** The ncRNA boundary feature vector is always exactly this many dims (ani,
	 * prof0-2, isStop, fuzz0-2, lengthRatio, contigGC) -- no stem, no optional
	 * cross-boundary-enrichment variant, unlike TrnaBoundaryScorer's 11-or-14. */
	public static final int NUM_FEATURES=10;
	static final int FEATURES_V1=1,FEATURES_V2=2,FEATURES_V3=3;

	static int parseFeatureVersion(String value){
		if(value.equalsIgnoreCase("v1") || value.equals("1") || value.equalsIgnoreCase("current")){return FEATURES_V1;}
		if(value.equalsIgnoreCase("v2") || value.equals("2")){return FEATURES_V2;}
		if(value.equalsIgnoreCase("v3") || value.equals("3")){return FEATURES_V3;}
		throw new IllegalArgumentException("boundaryfeatures must be v1, v2, or v3: "+value);
	}

	/** V1 alone performs candidate-site Scrabble alignment. V2/V3 are structurally
	 * alignment-free at endpoint sites and reuse the accepted whole-locus identity. */
	static boolean usesPerSiteScrabble(int featureVersion){return featureVersion==FEATURES_V1;}

	/** Loads a boundary net and fails loud if its declared input dimension isn't
	 * exactly NUM_FEATURES (Citan's explicit requirement, fail-fast on dimension
	 * mismatch): a net trained at a different dim was not trained by
	 * NcrnaBoundaryVectorGen's current 10-dim format, or is corrupt, and must
	 * halt the whole process rather than silently score garbage. */
	public static CellNet load(String path){
		final CellNet net=CellNetParser.load(path);
		final int dims=net.numInputs();
		if(dims!=NUM_FEATURES){throw new IllegalArgumentException("Ncrna boundary net '"+path+"' has "
			+dims+" inputs; expected "+NUM_FEATURES+" for the ncRNA boundary feature vector.");}
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
		final float modelPlateau=(model==null ? 0 : TrnaBoundaryFeatures.modelMaxCoverage(model));
		return score(net,window,s,e,isStop,modelConsensus,model,table,insideCount,outsideCount,contigGC,meanLen,modelPlateau);
	}

	private static float score(CellNet net, byte[] window, int s, int e, boolean isStop,
			byte[] modelConsensus, BaseGraph model, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, float contigGC, float meanLen, float modelPlateau){
		return score(net,window,s,e,isStop,modelConsensus,model,table,insideCount,outsideCount,
			contigGC,meanLen,modelPlateau,FEATURES_V1,Float.NaN,null,null,null);
	}

	private static float score(CellNet net, byte[] window, int s, int e, boolean isStop,
			byte[] modelConsensus, BaseGraph model, TrnaBoundaryFeatures.NinemerTable table,
			int insideCount, int outsideCount, float contigGC, float meanLen, float modelPlateau,
			int featureVersion, float locusAni, float[] startFuzz, float[] stopFuzz,
			TrnaBoundaryFeatures.NinemerTable familyTable){
		final boolean useStart=!isStop;
		final int boundaryPos=(isStop ? e : s);
		final TrnaBoundaryFeatures.BoundaryType type=(isStop
			? TrnaBoundaryFeatures.BoundaryType.STOP : TrnaBoundaryFeatures.BoundaryType.START);
		final float[] prof=TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos, type, insideCount, outsideCount, table);
		if(familyTable!=null){
			if(featureVersion!=FEATURES_V3){throw new IllegalArgumentException("Family/consensus table blending requires boundaryfeatures=v3");}
			final float[] familyProf=TrnaBoundaryFeatures.enrichmentProfile(window, boundaryPos, type, insideCount, outsideCount, familyTable);
			NcrnaBoundaryVectorGen.blendProfiles(prof,familyProf);
		}
		final float ani;
		final float[] fuzz;
		if(featureVersion==FEATURES_V2){
			if(!Float.isFinite(locusAni)){throw new IllegalArgumentException("boundaryfeatures=v2 requires finite locus ANI from QuantumAligner");}
			fuzz=(useStart ? startFuzz : stopFuzz);
			if(fuzz==null || fuzz.length!=3){throw new IllegalArgumentException("boundaryfeatures=v2 requires three HBM constants per model end");}
			ani=locusAni;
		}else if(featureVersion==FEATURES_V3){
			if(!Float.isFinite(locusAni)){throw new IllegalArgumentException("boundaryfeatures=v3 requires finite locus identity from QuantumAligner");}
			ani=locusAni;fuzz=null;
		}else if(usesPerSiteScrabble(featureVersion)){
			final byte[] candSeq=java.util.Arrays.copyOfRange(window, s, e+1);
			final float[] aniFuzz=TrnaBoundaryFeatures.aniAndFuzzinessFeature(candSeq, modelConsensus, model, useStart, modelPlateau);
			ani=aniFuzz[0]; fuzz=new float[]{aniFuzz[1],aniFuzz[2],aniFuzz[3]};
		}else{
			throw new IllegalArgumentException("Unknown ncRNA boundary feature version: "+featureVersion);
		}
		final float lengthRatio=(e-s+1)/meanLen;
		final float[] in=(featureVersion==FEATURES_V3
			? buildV3Input(prof,isStop,lengthRatio,contigGC,locusAni)
			: buildInput(ani,prof,isStop,fuzz,lengthRatio,contigGC));
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

	/** V3 order: profile[0..2], isStop, lengthRatio, contigGC, Quantum identity,
	 * then three reserved zero slots. Package-visible for exact-layout fixtures. */
	static float[] buildV3Input(float[] prof, boolean isStop, float lengthRatio, float contigGC, float locusIdentity){
		if(!Float.isFinite(locusIdentity)){throw new IllegalArgumentException("boundaryfeatures=v3 requires finite locus identity from QuantumAligner");}
		return new float[]{prof[0],prof[1],prof[2],(isStop ? 1f : 0f),lengthRatio,contigGC,locusIdentity,0f,0f,0f};
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
		return refineBoundaries(startNet, stopNet, window, s, e, modelConsensus, model,
			startTable, stopTable, startInside, startOutside, stopInside, stopOutside,
			contigGC, meanLen, startOffsets, stopOffsets,
			MARGIN_THRESHOLD_START, MARGIN_THRESHOLD_STOP);
	}

	/** Family-scoped refinement. Margins travel with the family so an R58
	 * experiment cannot alter another ncRNA family's boundary behavior. */
	public static int[] refineBoundaries(CellNet startNet, CellNet stopNet, byte[] window, int s, int e,
			byte[] modelConsensus, BaseGraph model,
			TrnaBoundaryFeatures.NinemerTable startTable, TrnaBoundaryFeatures.NinemerTable stopTable,
			int startInside, int startOutside, int stopInside, int stopOutside, float contigGC, float meanLen,
			int[] startOffsets, int[] stopOffsets, float marginStart, float marginStop){
		return refineBoundariesWithCounts(startNet, stopNet, window, s, e, modelConsensus, model,
			startTable, stopTable, startInside, startOutside, stopInside, stopOutside,
			contigGC, meanLen, startOffsets, stopOffsets, marginStart, marginStop, false, false).offsets();
	}

	/** Refinement result plus the measured number of sweep sites evaluated for each endpoint. */
	static final class RefinementResult{
		RefinementResult(int startOffset_, int stopOffset_, int startSitesScored_, int stopSitesScored_,
				float startChosenScore_, float stopChosenScore_){
			startOffset=startOffset_; stopOffset=stopOffset_;
			startSitesScored=startSitesScored_; stopSitesScored=stopSitesScored_;
			startChosenScore=startChosenScore_; stopChosenScore=stopChosenScore_;
		}
		int[] offsets(){return new int[]{startOffset, stopOffset};}
		final int startOffset, stopOffset, startSitesScored, stopSitesScored;
		final float startChosenScore, stopChosenScore;
	}

	/** Same exhaustive refinement with explicit contig-edge state for the G16 radius proof. */
	static RefinementResult refineBoundariesWithCounts(CellNet startNet, CellNet stopNet, byte[] window, int s, int e,
			byte[] modelConsensus, BaseGraph model,
			TrnaBoundaryFeatures.NinemerTable startTable, TrnaBoundaryFeatures.NinemerTable stopTable,
			int startInside, int startOutside, int stopInside, int stopOutside, float contigGC, float meanLen,
			int[] startOffsets, int[] stopOffsets, float marginStart, float marginStop,
			boolean touchesContigStart, boolean touchesContigEnd){
		return refineBoundariesWithCounts(startNet,stopNet,window,s,e,modelConsensus,model,startTable,stopTable,
			startInside,startOutside,stopInside,stopOutside,contigGC,meanLen,startOffsets,stopOffsets,
			marginStart,marginStop,touchesContigStart,touchesContigEnd,FEATURES_V1,Float.NaN,null,null);
	}

	/** V2 shares one Quantum identity across the locus and reads fixed HBM endpoint
	 * constants. V3 shares the identity but replaces those constants with zero spares. */
	static RefinementResult refineBoundariesWithCounts(CellNet startNet, CellNet stopNet, byte[] window, int s, int e,
			byte[] modelConsensus, BaseGraph model,
			TrnaBoundaryFeatures.NinemerTable startTable, TrnaBoundaryFeatures.NinemerTable stopTable,
			int startInside, int startOutside, int stopInside, int stopOutside, float contigGC, float meanLen,
			int[] startOffsets, int[] stopOffsets, float marginStart, float marginStop,
			boolean touchesContigStart, boolean touchesContigEnd, int featureVersion, float locusAni,
			float[] startFuzz, float[] stopFuzz){
		return refineBoundariesWithCounts(startNet,stopNet,window,s,e,modelConsensus,model,startTable,stopTable,
			startInside,startOutside,stopInside,stopOutside,contigGC,meanLen,startOffsets,stopOffsets,
			marginStart,marginStop,touchesContigStart,touchesContigEnd,featureVersion,locusAni,startFuzz,stopFuzz,null,null);
	}

	/** V3 inference with the same independently normalized 50/50 family/consensus
	 * profile blend used by NcrnaBoundaryVectorGen. Null family tables preserve the
	 * historical single-table path byte-for-byte. */
	static RefinementResult refineBoundariesWithCounts(CellNet startNet, CellNet stopNet, byte[] window, int s, int e,
			byte[] modelConsensus, BaseGraph model,
			TrnaBoundaryFeatures.NinemerTable startTable, TrnaBoundaryFeatures.NinemerTable stopTable,
			int startInside, int startOutside, int stopInside, int stopOutside, float contigGC, float meanLen,
			int[] startOffsets, int[] stopOffsets, float marginStart, float marginStop,
			boolean touchesContigStart, boolean touchesContigEnd, int featureVersion, float locusAni,
			float[] startFuzz, float[] stopFuzz,
			TrnaBoundaryFeatures.NinemerTable familyStartTable,
			TrnaBoundaryFeatures.NinemerTable familyStopTable){
		if(featureVersion!=FEATURES_V1 && featureVersion!=FEATURES_V2 && featureVersion!=FEATURES_V3){throw new IllegalArgumentException(
			"Unknown ncRNA boundary feature version: "+featureVersion);}
		if((familyStartTable==null)!=(familyStopTable==null)){throw new IllegalArgumentException("Family start/stop tables must be supplied together");}
		if(familyStartTable!=null && featureVersion!=FEATURES_V3){throw new IllegalArgumentException("Family/consensus table blending requires boundaryfeatures=v3");}
		final float modelPlateau=(featureVersion==FEATURES_V1 && model!=null ? TrnaBoundaryFeatures.modelMaxCoverage(model) : 0);
		final float startConf=score(startNet, window, s, e, false, modelConsensus, model,
			startTable, startInside, startOutside, contigGC, meanLen, modelPlateau,featureVersion,locusAni,startFuzz,stopFuzz,familyStartTable);
		final float stopConf=score(stopNet, window, s, e, true, modelConsensus, model,
			stopTable, stopInside, stopOutside, contigGC, meanLen, modelPlateau,featureVersion,locusAni,startFuzz,stopFuzz,familyStopTable);

		final float[] startSweep, stopSweep;
		if(startConf<=stopConf){//start is the worse (or tied) boundary -- refine it first
			startSweep=bestOffset(startNet, window, s, e, false, modelConsensus, model,
				startTable, startInside, startOutside, contigGC, meanLen, startOffsets, touchesContigStart, modelPlateau,
				featureVersion,locusAni,startFuzz,stopFuzz,familyStartTable);
			final int bestStartOffset=applyMargin(startSweep, marginStart);
			stopSweep=bestOffset(stopNet, window, s+bestStartOffset, e, true, modelConsensus, model,
				stopTable, stopInside, stopOutside, contigGC, meanLen, stopOffsets, touchesContigEnd, modelPlateau,
				featureVersion,locusAni,startFuzz,stopFuzz,familyStopTable);
		}else{
			stopSweep=bestOffset(stopNet, window, s, e, true, modelConsensus, model,
				stopTable, stopInside, stopOutside, contigGC, meanLen, stopOffsets, touchesContigEnd, modelPlateau,
				featureVersion,locusAni,startFuzz,stopFuzz,familyStopTable);
			final int bestStopOffset=applyMargin(stopSweep, marginStop);
			startSweep=bestOffset(startNet, window, s, e+bestStopOffset, false, modelConsensus, model,
				startTable, startInside, startOutside, contigGC, meanLen, startOffsets, touchesContigStart, modelPlateau,
				featureVersion,locusAni,startFuzz,stopFuzz,familyStartTable);
		}
		return new RefinementResult(applyMargin(startSweep, marginStart), applyMargin(stopSweep, marginStop),
			(int)startSweep[3], (int)stopSweep[3], selectedScore(startSweep, marginStart), selectedScore(stopSweep, marginStop));
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

	/** Score at the endpoint actually selected after the margin gate. */
	static float selectedScore(float[] sweepResult, float margin){
		return applyMargin(sweepResult, margin)==0 ? sweepResult[2] : sweepResult[1];
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
			int insideCount, int outsideCount, float contigGC, float meanLen, int[] offsets,
			boolean touchesRelevantContigEnd, float modelPlateau, int featureVersion, float locusAni,
			float[] startFuzz, float[] stopFuzz, TrnaBoundaryFeatures.NinemerTable familyTable){
		int bestOff=0, sitesScored=0; float bestScore=-Float.MAX_VALUE, zeroScore=Float.NaN;
		for(int offset : offsets){
			final int s=s0+(isStop ? 0 : offset), e=e0+(isStop ? offset : 0);
			if(s<0 || e>=window.length || e-s<15){continue;}
			final float sc=score(net, window, s, e, isStop, modelConsensus, model, table, insideCount, outsideCount,
				contigGC, meanLen, modelPlateau,featureVersion,locusAni,startFuzz,stopFuzz,familyTable);
			sitesScored++;
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
		assertCompleteCenteredSweep(offsets, sitesScored, touchesRelevantContigEnd,
			(isStop ? e0 : s0), isStop, window.length);
		return new float[]{bestOff, bestScore, zeroScore, sitesScored};
	}

	/** Brian via G11, 2026-09-22: every site in -r..+r must be scored.  This
	 * plain assert follows the assertions skill because the sweep is synchronous
	 * on the calling thread; failure cannot strand a producer/consumer pipeline.
	 * It catches the fixed-PAD=10 regression that silently collapsed radii
	 * 12/16/20/25. */
	static void assertCompleteCenteredSweep(int[] offsets, int sitesScored,
			boolean touchesRelevantContigEnd, int position, boolean isStop, int windowLength){
		final int radius=centeredRadius(offsets);
		assert(radius<0 || touchesRelevantContigEnd || sitesScored==2*radius+1) :
			"Incomplete boundary sweep: r="+radius+" sitesScored="+sitesScored+" position="+position
			+" endpoint="+(isStop ? "stop" : "start")+" windowLength="+windowLength;
	}

	/** Returns r only for a complete ordered {-r,...,0,...,+r} array; otherwise -1. */
	static int centeredRadius(int[] offsets){
		if(offsets==null || offsets.length<1 || (offsets.length&1)==0){return -1;}
		final int radius=offsets.length/2;
		for(int i=0; i<offsets.length; i++){if(offsets[i]!=i-radius){return -1;}}
		return radius;
	}
}
