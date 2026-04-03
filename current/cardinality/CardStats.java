package cardinality;

import shared.Tools;

/**
 * Concrete cardinality statistics holder.
 * <p>
 * Single constructor computes all estimates from an NLZ histogram (counts[])
 * and optionally a per-bucket packed-value array (buckets[]).
 * All fields are final and set in dependency order.
 * No toArray() — callers use named getters or buildLegacyArray() for migration.
 * <p>
 * Constructor phases:
 * <ol>
 *   <li>Phase 1 (counts-only): LC, LCMin, DLC, HLL, Mean, GMean, HMean, CF-corrected values</li>
 *   <li>Phase 2a (history): sbs, HC, LDLC — when histBits &gt; 0</li>
 *   <li>Phase 2b (mantissa): hmeanMantissa, hybridDDL — when mantissaBits &gt; 0</li>
 *   <li>Phase 3 (null buckets): remaining fields cloned from phase-1 equivalents</li>
 *   <li>Final: hybridDLL, hybridDDL — depend on lcForHybrid which may be set by phase 2a</li>
 * </ol>
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class CardStats extends AbstractCardStats {

	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a CardStats from NLZ histogram and optional per-bucket data.
	 *
	 * @param buckets_       char[numBuckets] packed values [nlz|hist|luck|mantissa], or null
	 * @param counts_        int[66] NLZ histogram: [0]=empties, [k+1]=buckets at absNlz==k. MANDATORY.
	 * @param nlzBits_       bits for NLZ in each bucket value (e.g. 6 for PLL16c)
	 * @param histBits_      history bit width (0 when absent)
	 * @param luckBits_      luck bit width (0 when absent)
	 * @param mantissaBits_  mantissa bit width (0 when absent)
	 * @param numBuckets_    total bucket count
	 * @param microIndex_    64-bit micro-index bitmap
	 * @param added_         total elements added (for clamping)
	 * @param cfMatrix_      per-class CF matrix from resource file
	 * @param cfBuckets_     bucket count used when building cfMatrix
	 * @param mantissaOff_   fractional NLZ offset for mantissa classes (0.5 for DDL2/DDL8, 0.0 for DDL)
	 */
	CardStats(final char[] buckets_, final int[] counts_,
			final int nlzBits_, final int histBits_, final int luckBits_, final int mantissaBits_,
			final int numBuckets_, final long microIndex_, final long added_,
			final float[][] cfMatrix_, final int cfBuckets_,
			final double mantissaOff_){

		assert(counts_!=null && counts_.length==66);
		assert((histBits_|luckBits_|mantissaBits_)==0 || buckets_!=null) :
			"Per-bucket array required when sub-NLZ bits are present";

		/*--- Store configuration ---*/
		numBuckets=numBuckets_;
		final long effectiveMicro=USE_MICRO_INDEX ? microIndex_ : 0;
		microBits=(int)Long.bitCount(effectiveMicro);
		added=added_;
		nlzBits=nlzBits_;
		histBits=histBits_;
		luckBits=luckBits_;
		mantissaBits=mantissaBits_;
		counts=counts_;
		cfMatrix=cfMatrix_;
		cfBuckets=cfBuckets_;

		/*--- Phase 1: counts-only estimates (all classes) ---*/

		// Derive fundamental sums from NLZ histogram (counts[0]=empties, counts[k+1]=absNlz k)
		final double[] sums=sumsFromCounts(counts);
		difSum=sums[0];
		hllSumFilled=sums[1];
		gSum=sums[2];
		filled=(int)sums[3];
		assert filled+counts[0]==numBuckets :
			"counts sum ("+filled+"+"+counts[0]+") != numBuckets ("+numBuckets+")";

		// Empty buckets: counts[0], adjusted by microIndex for improved LC at low card
		final int microFilled=USE_MICRO_FOR_LC ? microBits : 0;
		V=numBuckets-Math.max(filled, microFilled);

		// Bias correction constants
		alpha_m=0.7213/(1.0+1.079/numBuckets);
		// NOTE: float cast matches CardinalityStats line 144 exactly — do not change to double
		correction=(filled+numBuckets)/(float)(numBuckets+numBuckets);

		// MicroIndex estimate: LC over 64 virtual buckets
		microEst=microEstimate(microBits);

		// LC floor from history bits and micro count
		final int lcFloor=(CardinalityTracker.USE_HISTORY_FOR_LC ? microFilled : 0);

		// 1. LC raw — basic linear counting
		lcRawF=Math.max(lcRaw(V, numBuckets), lcFloor);

		// 2. LCMin — tier-compensated LC (walks counts to find first tier with empties)
		lcMinF=lcMin(counts, V, numBuckets, lcFloor);

		// 3. DLC — primary CF-free estimator; info-power weighted blend (mode 2 = 0-param theory)
		//    dlcPure is the raw tier blend; dlcRaw blends it with lcMin at low occupancy.
		//    Must be set before any CF lookups since dlcRaw is used as the cardinality seed.
		dlcPureF=dlcPure(counts, V, numBuckets, DLC_INFO_POWER);
		dlcRawF=dlcBlendWithLcMin(dlcPureF, lcMinF, V, numBuckets, DLC_BLEND_LO, DLC_BLEND_HI);

		// 4. DLC variants
		dlcBestF=dlcBest(counts, V, numBuckets, lcMinF);
		dlc3bF=dlcBlend3(counts, V, numBuckets, lcMinF);
		dlcLowestF=lcMinF; // tier-compensated LC at actual floor; no separate computation needed

		// 5. HLL — all-buckets harmonic mean with LC fallback
		hllRawF=hllEstimate(hllSumFilled, V, numBuckets, alpha_m);

		// 6. HLL with history-corrected constant
		hllHistoryF=hllHistory(hllSumFilled, V, numBuckets, alpha_m);

		// 7. Mean, GMean, HMean — raw (pre-CF) estimates from counts-derived sums
		meanRawF=meanEstimate(difSum, filled, numBuckets);
		gmeanRawF=gmeanEstimate(gSum, filled, numBuckets);
		hmeanRawF=hmeanEstimate(hllSumFilled, filled, alpha_m);

		// 8. CF-corrected estimates — use DLC raw as cardinality seed for table lookup
		//    keyScale adjusts for table built at different bucket count than current
		final double keyScale=(CorrectionFactor.v1Buckets>0 ?
				(double)CorrectionFactor.v1Buckets/numBuckets : 1.0);
		meanCF=meanRawF*cardCF(meanRawF, CorrectionFactor.MEAN, keyScale);
		hmeanCF=hmeanRawF*cardCF(hmeanRawF, CorrectionFactor.HMEAN, keyScale);
		// NOTE: hmeanMCF uses the HMEANM CF row (index 3), not HMEAN (index 2),
		// because the mantissa-corrected harmonic mean has different bias characteristics.
		gmeanCF=gmeanRawF*cardCF(gmeanRawF, CorrectionFactor.GMEAN, keyScale);
		hllCF=hllRawF*cardCF(hllRawF, CorrectionFactor.HLL, keyScale);

		// DLC CF-corrected (DLC is already CF-free but gets table correction for residual bias)
		dlcCF=dlcRawF*cardCF(dlcRawF, CorrectionFactor.DLC, keyScale);
		dlc3bCF=dlc3bF*cardCF(dlc3bF, CorrectionFactor.DLC3B, keyScale);
		dlcBestCF=dlcBestF*cardCF(dlcBestF, CorrectionFactor.DLCBEST, keyScale);

		// HybridDLC and HybDLC50 are computed after phase 2b
		// because mantissa classes need mantissa-corrected mean for blending

		/*--- Phase 2a: history-corrected estimates (when histBits > 0) ---*/

		// Default lcForHybrid to lcMin; may be overridden by sbs below
		double lcForHybridTmp=lcMinF;
		double corrDifSum=difSum, corrGSum=gSum;
		boolean hasHistCorr=false;

		byte[] sbsStatesArr=null;
		if(histBits_>0 && buckets_!=null){
			hasHistCorr=true;
			final int histMask=(1<<histBits_)-1;
			final double termCF=StateTable.terminalCF(histBits_, 0);
			double cDif=0, cGSum=0;
			sbsStatesArr=new byte[numBuckets];
			for(int i=0; i<numBuckets; i++){
				final int val=buckets_[i];
				if(val==0){
					sbsStatesArr[i]=-1; // empty bucket
					continue;
				}
				final int absNlz=(val>>>histBits_)-1;
				final int histPattern=val&histMask;
				final int nlzBin=Math.min(absNlz, histBits_+1);

				// Per-bucket dif: must match old UDLL6 formula exactly
				final long dif=(absNlz==0 ? Long.MAX_VALUE : (absNlz<64 ? 1L<<(63-absNlz) : 1L));

				// Per-state correction for Mean/GMean only
				final double cf=StateTable.historyOffset(nlzBin, histBits_, histPattern);
				final double tierMult=Math.pow(2.0, -(cf+StateTable.CF_OFFSET))/termCF;

				cDif+=dif*tierMult;
				cGSum+=Math.log(Math.max(1, dif*tierMult));

				// LC history state index for sbs computation
				sbsStatesArr[i]=(byte)CorrectionFactor.sbsStateIndex(nlzBin, histPattern, histBits_);
			}
			corrDifSum=cDif;
			corrGSum=cGSum;
		}
		hasHistoryCorrection=hasHistCorr;

		// History-corrected Mean, GMean
		if(hasHistCorr){
			final int divH=Math.max(filled, 1);
			final double meanH=corrDifSum/divH;
			meanCorrRawF=2*(Long.MAX_VALUE/Tools.max(1.0, meanH))*divH*correction;
			meanCorrCF=meanCorrRawF*cardCF(meanCorrRawF, CorrectionFactor.MEAN, keyScale);
			final double gmeanH=Math.exp(corrGSum/divH);
			gmeanCorrRawF=2*(Long.MAX_VALUE/gmeanH)*divH*correction;
			gmeanCorrCF=gmeanCorrRawF*cardCF(gmeanCorrRawF, CorrectionFactor.GMEAN, keyScale);
		}else{
			meanCorrRawF=meanRawF;
			meanCorrCF=meanCF;
			gmeanCorrRawF=gmeanRawF;
			gmeanCorrCF=gmeanCF;
		}

		// SBS: history-aware LC from per-bucket state indices
		sbsF=computeSbs(sbsStatesArr, filled, numBuckets, lcRawF);
		sbsMultF=computeSbsMult(sbsStatesArr, filled, numBuckets, lcRawF);

		// If sbs table is available and enabled, use sbs for hybrid blending
		if(CardinalityTracker.USE_SBS_IN_HYBRID
				&& CorrectionFactor.SBS_CF_TABLE!=null
				&& histBits_>0 && sbsStatesArr!=null){
			lcForHybridTmp=sbsF;
		}
		lcForHybridF=lcForHybridTmp;

		// DSBS: DLC with sbs as fallback, extended blend zone
		dlcSbsF=dlcInfoPow(counts, V, numBuckets, sbsF, DLC_INFO_POWER, DLCSBS_BLEND_LO, DLCSBS_BLEND_HI);

		/*--- Phase 2b: mantissa estimates (when mantissaBits > 0) ---*/

		if(mantissaBits_>0 && buckets_!=null){
			// Phase 2b: scan buckets for mantissa-corrected sums.
			// The mantissa provides fractional NLZ precision beyond integer tier resolution.
			// All DDL classes use inverted mantissa: higher stored = smaller hash within tier.
			final int mantissaMask=(1<<mantissaBits_)-1;
			final double mantissaScale=(double)(1<<mantissaBits_);
			double hllSumM=0, difSumM=0, gSumM=0;
			int filledM=0;
			for(int i=0; i<numBuckets; i++){
				final int val=buckets_[i];
				if(val==0){continue;} // empty or phantom bucket
				filledM++;
				// Packed format: ((absNlz+1) << mantissaBits) | invMant
				// The +1 ensures val is never 0 for filled buckets
				final int absNlz=(val>>>mantissaBits_)-1;
				final int invMant=val&mantissaMask;
				// Fractional NLZ: exponent = -absNlz + offset - invMant/mantissaScale
				// offset is class-dependent: 0.5 for DDL2/DDL8, 0.0 for DDL
				hllSumM+=Math.pow(2.0, -absNlz+mantissaOff_-invMant/mantissaScale);
				// Mantissa-precise dif: reconstruct with uninverted mantissa bits
				final int lowbits=(~val)&mantissaMask;
				final long mantissa=(1L<<mantissaBits_)|lowbits;
				final int shift=WORDLEN-absNlz-mantissaBits_-1;
				final long dif;
				if(absNlz==0){dif=Long.MAX_VALUE;}
				else if(shift>=0){dif=mantissa<<shift;}
				else{dif=1L;}
				difSumM+=dif;
				gSumM+=Math.log(Math.max(1, dif));
			}
			hllSumFilledM=hllSumM;
			difSumM_=difSumM;
			gSumM_=gSumM;
			// Mantissa-corrected harmonic mean
			hmeanMRawF=(filledM==0 ? 0 : 2*alpha_m*(double)filledM*(double)filledM/hllSumM);
			hmeanMCF=hmeanMRawF*cardCF(hmeanMRawF, CorrectionFactor.HMEANM, keyScale);
			// Mantissa-corrected Mean and GMean (more precise than counts-only)
			final int divM=Math.max(filledM, 1);
			final double meanM=difSumM/divM;
			final float correctionM=(filledM+numBuckets)/(float)(numBuckets+numBuckets);
			meanMRawF=2*(Long.MAX_VALUE/Tools.max(1.0, meanM))*divM*correctionM;
			meanMCF_=meanMRawF*cardCF(meanMRawF, CorrectionFactor.MEAN, keyScale);
			final double gmeanM=Math.exp(gSumM/divM);
			gmeanMRawF=2*(Long.MAX_VALUE/gmeanM)*divM*correctionM;
			gmeanMCF_=gmeanMRawF*cardCF(gmeanMRawF, CorrectionFactor.GMEAN, keyScale);
		}else{
			// No mantissa: HMeanM = HMean, no CF correction applied.
			// The old code intentionally skips CF when hasMantissa==false because
			// the HMEANM CF was trained on mantissa-corrected data and would be wrong here.
			hllSumFilledM=hllSumFilled;
			difSumM_=difSum;
			gSumM_=gSum;
			hmeanMRawF=hmeanRawF;
			hmeanMCF=hmeanRawF; // raw, no CF
			meanMRawF=meanRawF;
			meanMCF_=meanCF;
			gmeanMRawF=gmeanRawF;
			gmeanMCF_=gmeanCF;
		}

		// hasMantissa: true when mantissa-corrected sums differ from plain sums
		// Currently false until phase 2b is implemented with actual mantissa correction
		hasMantissa=(mantissaBits_>0 && buckets_!=null);

		/*--- Phase 2c: luck estimates (when luckBits > 0) ---*/

		// TODO: luck implementation

		// Standard hybrids use mantissa-corrected or plain Mean (no history).
		// History-corrected hybrids are computed separately for +H columns.
		final double meanForHybrid=(hasMantissa ? meanMCF_ : meanCF);
		final double meanRawForHybrid=(hasMantissa ? meanMRawF : meanRawF);

		// HybridDLC: DLC → Mean blend (uses raw meanEst, CF applied to whole result)
		final double hybDLCraw=hybridDLC(dlcRawF, meanRawForHybrid, numBuckets);
		hybDLCcf=hybDLCraw*cardCF(hybDLCraw, CorrectionFactor.HYBDLC, keyScale);

		// HybDLC50: 50/50 DLC+Mean at high cardinality
		hybDLC50CF=hybDLC50Internal(meanRawForHybrid, keyScale);

		/*--- DlcThreshHybrid (uses DLC for zone detection, LC+Mean for values) ---*/
		dThHybF=dlcThreshHybridInternal(lcMinF, dlcRawF, meanForHybrid, hmeanMCF, numBuckets, hasMantissa);

		/*--- Final phase: standard hybrids (plain Mean, no history correction) ---*/
		hybridDLLF=hybridDLL(lcMinF, lcMinF, meanForHybrid, numBuckets);
		hybridDDLF=hybridDDL(lcMinF, lcMinF, meanForHybrid, hmeanMCF, numBuckets);

		/*--- History-corrected hybrids (use sbs and history-corrected Mean) ---*/
		final double meanForHybridH=(hasMantissa ? meanMCF_ : meanCorrCF);
		hybridDLLHistF=hybridDLL(lcForHybridF, lcMinF, meanForHybridH, numBuckets);
		hybridDDLHistF=hybridDDL(lcForHybridF, lcMinF, meanForHybridH, hmeanMCF, numBuckets);

		/*--- Hybrid+2: SBS → Mean+H blend, DLC as zone detector ---*/
		hybridPlus2F=hybridDLL(lcForHybridF, dlcRawF, meanForHybridH, numBuckets);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Private CF Helpers         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Cardinality-indexed CF lookup using DLC raw as seed.
	 * Handles v1, v3, and legacy table versions.
	 * Matches the logic of old CardinalityStats.cf() exactly.
	 */
	private double cardCF(final double rawEst, final int type, final double keyScale){
		// Formula mode for Mean CF: bypasses table entirely
		if(CorrectionFactor.USE_MEAN_CF_FORMULA && type==CorrectionFactor.MEAN){
			return CorrectionFactor.meanCfFormula(dlcRawF);
		}
		// v3+ tables: use dlcRawF as seed (the primary DLC estimator)
		if(CorrectionFactor.tableVersion>=3){
			if(!CorrectionFactor.USE_CORRECTION || CorrectionFactor.v1Matrix==null
					|| type==CorrectionFactor.LINEAR || type>=CorrectionFactor.v1Matrix.length){return 1;}
			final int iters=(dlcRawF*keyScale>MIN_SEED_CF_MULT*CorrectionFactor.v1Buckets ? DEFAULT_CF_ITERS : 1);
			return CorrectionFactor.getCF(dlcRawF, rawEst,
					CorrectionFactor.v1Matrix[type], CorrectionFactor.v1Keys,
					iters, DEFAULT_CF_DIF, keyScale);
		}
		// v1+ tables: use dlc3bF as seed (legacy DLC3B-indexed)
		if(CorrectionFactor.tableVersion>=1){
			if(!CorrectionFactor.USE_CORRECTION || CorrectionFactor.v1Matrix==null
					|| type==CorrectionFactor.LINEAR || type>=CorrectionFactor.v1Matrix.length){return 1;}
			final int iters=(dlc3bF*keyScale>MIN_SEED_CF_MULT*CorrectionFactor.v1Buckets ? DEFAULT_CF_ITERS : 1);
			return CorrectionFactor.getCF(dlc3bF, rawEst,
					CorrectionFactor.v1Matrix[type], CorrectionFactor.v1Keys,
					iters, DEFAULT_CF_DIF, keyScale);
		}
		// Legacy: occupancy-indexed + optional cardinality table
		if(CorrectionFactor.lastCardMatrix!=null && CardinalityTracker.USE_CARD_CF
				&& CorrectionFactor.USE_CORRECTION){
			return CorrectionFactor.getCF(cfMatrix, cfBuckets,
					CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys,
					filled, numBuckets, meanRawF, type);
		}
		return CorrectionFactor.getCF(cfMatrix, cfBuckets, filled, numBuckets, type);
	}

	/**
	 * HybDLC50: like hybridDLC but blends to 50/50 DLC+Mean at high cardinality.
	 * Matches old CardinalityStats.hybDLC50() exactly.
	 */
	private double hybDLC50Internal(final double meanRaw, final double keyScale){
		final double high=0.5*dlcRawF+0.5*meanRaw;
		final double hb0=0.20*numBuckets, hb1=5.0*numBuckets;
		final double raw;
		if(dlcRawF<=hb0){raw=dlcRawF;}
		else if(dlcRawF<hb1){
			final double t=Math.log(dlcRawF/hb0)/Math.log(hb1/hb0);
			raw=(1-t)*dlcRawF+t*high;
		}else{raw=high;}
		return raw*cardCF(raw, CorrectionFactor.HYBDLC50, keyScale);
	}

	/**
	 * DlcThreshHybrid: uses DLC for zone detection, LC+Mean for blended values.
	 * Matches old CardinalityStats.dlcThreshHybrid() exactly.
	 */
	private static double dlcThreshHybridInternal(final double lcForHybrid, final double dlc,
			final double meanCF, final double hmeanMCF, final int B, final boolean hasMant){
		final double hb0=0.20*B, hb1=5.0*B;
		if(dlc<=hb0){return lcForHybrid;}
		if(!hasMant){
			// DLL types: simple LC→Mean(CF) blend
			if(dlc<hb1){
				final double t=Math.log(dlc/hb0)/Math.log(hb1/hb0);
				return (1-t)*lcForHybrid+t*meanCF;
			}
			return meanCF;
		}
		// DDL types: LC→Mean(CF)→HMeanM(CF) blend
		final double hbMid1=1.0*B, hbMid2=2.5*B;
		final double t=Math.log(dlc/hb0)/Math.log(hb1/hb0);
		if(dlc<=hbMid1){
			return (1-t)*lcForHybrid+t*meanCF;
		}else if(dlc<=hbMid2){
			final double mix=(hbMid2-dlc)/(hbMid2-hbMid1);
			final double blended=meanCF*mix+hmeanMCF*(1-mix);
			return (1-t)*lcForHybrid+t*blended;
		}else if(dlc<hb1){
			return (1-t)*lcForHybrid+t*hmeanMCF;
		}
		return hmeanMCF;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Package-Private CF       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * CF lookup for callers that need to apply CF externally (e.g. DLL4m.cardinality()).
	 * Package-private; only for migration. Callers should eventually use pre-corrected getters.
	 */
	double cf(final double rawEst, final int type){
		final double ks=(CorrectionFactor.v1Buckets>0 ?
				(double)CorrectionFactor.v1Buckets/numBuckets : 1.0);
		return cardCF(rawEst, type, ks);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Getters: Counts-Only      ----------------*/
	/*--------------------------------------------------------------*/

	/** Basic linear counting estimate. */
	public double lcRaw(){return lcRawF;}
	/** Tier-compensated LC at the lowest non-full tier. */
	public double lcMin(){return lcMinF;}
	/** Primary DLC estimate (log-space exponential blend). CF-free. */
	public double dlcPure(){return dlcPureF;}
	public double dlcRaw(){return dlcRawF;}
	/** Single-best-tier DLC. CF-free. */
	public double dlcBest(){return dlcBestF;}
	/** 3-region DLC blend. CF-free. */
	public double dlcBlend3(){return dlc3bF;}
	/** DLC at lowest active tier = lcMin. */
	public double dlcLowest(){return dlcLowestF;}
	/** DLC with sbs fallback (instead of lcMin). Better at low cardinality with history. */
	public double dlcSbs(){return dlcSbsF;}
	/** HLL-style estimate (CF-corrected). */
	public double hll(){return hllCF;}
	/** HLL-style estimate, raw (no CF). */
	public double hllRaw(){return hllRawF;}
	/** HLL with history-corrected constant. */
	public double hllHistory(){return hllHistoryF;}
	/** Mean estimate (CF-corrected). Always plain counts-based, never history-corrected. */
	public double meanCF(){return meanCF;}
	/** Geometric mean estimate (CF-corrected). Always plain counts-based. */
	public double gmeanCF(){return gmeanCF;}
	/** Filled-bucket harmonic mean (CF-corrected). */
	public double hmeanCF(){return hmeanCF;}
	/** Mantissa-corrected harmonic mean (CF-corrected). Equals hmeanCF when no mantissa. */
	public double hmeanMCF(){return hmeanMCF;}
	/** Best available mean for standard output: mantissa-corrected when available, else plain. */
	public double meanMCF(){return hasMantissa ? meanMCF_ : meanCF;}
	/** Best available geometric mean: mantissa-corrected when available, else plain. */
	public double gmeanMCF(){return hasMantissa ? gmeanMCF_ : gmeanCF;}
	/** History-corrected mean (CF-corrected). Falls back to plain when no history. */
	public double meanHistCF(){return hasHistoryCorrection ? meanCorrCF : meanCF;}
	/** History-corrected geometric mean. Falls back to plain when no history. */
	public double gmeanHistCF(){return hasHistoryCorrection ? gmeanCorrCF : gmeanCF;}
	/** Hybrid DLL estimate (LC/Mean blend for mantissa-free classes). Plain, no history. */
	public double hybridDLL(){return hybridDLLF;}
	/** Hybrid DDL estimate (LC/Mean/HMeanM blend for mantissa classes). Plain. */
	public double hybridDDL(){return hybridDDLF;}
	/** History-corrected Hybrid DLL. Falls back to plain when no history. */
	public double hybridDLLHist(){return hybridDLLHistF;}
	/** History-corrected Hybrid DDL. Falls back to plain when no history. */
	public double hybridDDLHist(){return hybridDDLHistF;}
	/** Hybrid+2: SBS → Mean+H with DLC zone detection. */
	public double hybridPlus2(){return hybridPlus2F;}
	/** DLC-threshold hybrid. */
	public double dlcThreshHybrid(){return dThHybF;}
	/** MicroIndex-derived cardinality floor. 0 if micro disabled. */
	public long microCardinality(){return microEst;}

	/** DLC estimate at a given tier. */
	public double dlcTier(int tier){
		return dlcEstimate(tier, counts, V, numBuckets);
	}

	/*--------------------------------------------------------------*/
	/*----------------    SBS Computation        ----------------*/
	/*--------------------------------------------------------------*/

	/** Computes history-aware LC from per-bucket state indices.
	 *  Uses closed-form formula when USE_SBS_FORMULA is true, else SBS_CF_TABLE.
	 *  Falls back to lcRaw when table/formula and states are unavailable. */
	private static double computeSbs(final byte[] states, final int filled,
			final int numBuckets, final double lcRaw){
		if(states==null){return lcRaw;}

		// Formula mode: no table needed, B-independent
		if(CorrectionFactor.USE_SBS_FORMULA){
			double est=0;
			for(int i=0; i<states.length; i++){
				final int si=states[i];
				if(si>=0){est+=CorrectionFactor.sbsFormula(si, filled, numBuckets);}
			}
			return est;
		}

		// Table mode: original behavior
		final float[][] table=CorrectionFactor.SBS_CF_TABLE;
		if(table==null){return lcRaw;}
		final int tableBuckets=CorrectionFactor.sbsBuckets;
		final float[] row;
		if(numBuckets==tableBuckets){
			row=table[Math.min(filled, tableBuckets)];
		}else{
			final double pos=(double)filled*tableBuckets/numBuckets;
			final int lo=Math.max(1, Math.min((int)pos, tableBuckets-1));
			final int hi=Math.min(lo+1, tableBuckets);
			final double frac=pos-lo;
			row=new float[CorrectionFactor.SBS_STATES];
			final float[] rowLo=table[lo];
			final float[] rowHi=table[hi];
			for(int s=0; s<CorrectionFactor.SBS_STATES; s++){
				row[s]=(float)(rowLo[s]+(rowHi[s]-rowLo[s])*frac);
			}
		}
		double est=0;
		for(int i=0; i<states.length; i++){
			final int si=states[i];
			if(si>=0){est+=row[si];}
		}
		return est;
	}

	/** History-aware LC using multipliers: lcRaw × weightedAvg(CF per bucket state).
	 *  Falls back to lcRaw when table or states are unavailable. */
	private static double computeSbsMult(final byte[] states, final int filled,
			final int numBuckets, final double lcRaw){
		final float[][] table=CorrectionFactor.SBS_MULT_TABLE;
		if(table==null || states==null){return lcRaw;}
		final int tableBuckets=CorrectionFactor.sbsMultBuckets;
		final float[] row;
		if(numBuckets==tableBuckets){
			row=table[Math.min(filled, tableBuckets)];
		}else{
			final double pos=(double)filled*tableBuckets/numBuckets;
			final int lo=Math.max(1, Math.min((int)pos, tableBuckets-1));
			final int hi=Math.min(lo+1, tableBuckets);
			final double frac=pos-lo;
			row=new float[CorrectionFactor.SBS_STATES];
			final float[] rowLo=table[lo];
			final float[] rowHi=table[hi];
			for(int s=0; s<CorrectionFactor.SBS_STATES; s++){
				row[s]=(float)(rowLo[s]+(rowHi[s]-rowLo[s])*frac);
			}
		}
		double sumCF=0;
		int n=0;
		for(int i=0; i<states.length; i++){
			final int si=states[i];
			if(si>=0){sumCF+=row[si]; n++;}
		}
		if(n==0){return lcRaw;}
		return lcRaw*(sumCF/n);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Getters: History          ----------------*/
	/*--------------------------------------------------------------*/

	/** History-aware LC estimate. Falls back to lcRaw when history unavailable. */
	public double sbs(){return sbsF;}
	/** History-aware LC using multipliers. Falls back to lcRaw. */
	public double sbsMult(){return sbsMultF;}

	/*--------------------------------------------------------------*/
	/*----------------         Fields               ----------------*/
	/*--------------------------------------------------------------*/

	// --- Configuration (set once in constructor) ---
	final int numBuckets;       // total bucket count
	final int microBits;        // Long.bitCount(microIndex)
	final long added;           // total elements added
	final int nlzBits;          // NLZ bit width per bucket value
	final int histBits;         // history bit width (0=none)
	final int luckBits;         // luck bit width (0=none)
	final int mantissaBits;     // mantissa bit width (0=none)
	final int[] counts;         // NLZ histogram (int[64])
	final float[][] cfMatrix;   // per-class CF matrix
	final int cfBuckets;        // bucket count used to build cfMatrix
	final boolean hasMantissa;  // true when mantissa-corrected sums differ from plain

	// --- Counts-derived sums ---
	final double difSum;        // sum of dif values over filled buckets
	final double hllSumFilled;  // sum of 2^(-absNlz) over filled buckets
	final double gSum;          // sum of ln(dif) over filled buckets
	final int filled;           // number of filled buckets
	final int V;                // empty buckets (adjusted by microIndex)
	final double alpha_m;       // HLL bias correction: 0.7213/(1+1.079/B)
	final float correction;     // (filled+B) / (float)(2*B)

	// --- Micro estimate ---
	final int microEst;         // LC over 64 virtual micro-index buckets

	// --- Phase 1: counts-only estimates ---
	final double lcRawF;        // basic LC
	final double lcMinF;        // tier-compensated LC
	final double dlcPureF;      // pure tier blend (no lcMin transition)
	final double dlcRawF;       // primary DLC (dlcPure blended with lcMin), CF-free
	final double dlcBestF;      // best single-tier DLC, CF-free
	final double dlc3bF;        // 3-region DLC blend, CF-free
	final double dlcLowestF;    // DLC at lowest active tier = lcMin
	final double dlcSbsF;      // DLC with sbs fallback
	final double hllRawF;       // HLL raw (before CF)
	final double hllHistoryF;   // HLL with history-corrected constant
	final double meanRawF;      // Mean raw (before CF)
	final double gmeanRawF;     // GMean raw (before CF)
	final double hmeanRawF;     // HMean raw (filled-bucket, before CF)
	final double meanCF;        // Mean with CF correction
	final double gmeanCF;       // GMean with CF correction
	final double hmeanCF;       // HMean with CF correction
	final double hllCF;         // HLL with CF correction
	final double dlcCF;         // DLC with CF correction
	final double dlc3bCF;       // DLC3B with CF correction
	final double dlcBestCF;     // DLCBest with CF correction
	final double hybDLCcf;      // HybridDLC with CF correction
	final double hybDLC50CF;    // HybDLC50 with CF correction
	final double dThHybF;       // DLC-threshold hybrid

	// --- Phase 2a: history estimates (defaults to counts-only equivalents) ---
	final boolean hasHistoryCorrection; // true when per-state history corrections applied
	final double meanCorrRawF;  // history-corrected Mean raw (before CF)
	final double meanCorrCF;    // history-corrected Mean with CF
	final double gmeanCorrRawF; // history-corrected GMean raw
	final double gmeanCorrCF;   // history-corrected GMean with CF
	final double sbsF;       // history-aware LC
	final double sbsMultF;   // history-aware LC using multipliers
	final double lcForHybridF;  // LC value used for hybrid blending (sbs or lcMin)

	// --- Phase 2b: mantissa estimates (defaults to non-mantissa equivalents) ---
	final double hllSumFilledM; // mantissa-corrected harmonic indicator sum
	final double difSumM_;      // mantissa-corrected dif sum
	final double gSumM_;        // mantissa-corrected geometric sum
	final double hmeanMRawF;    // mantissa-corrected harmonic mean raw
	final double hmeanMCF;      // mantissa-corrected harmonic mean with CF
	final double meanMRawF;     // mantissa-corrected mean raw
	final double meanMCF_;      // mantissa-corrected mean with CF
	final double gmeanMRawF;    // mantissa-corrected geometric mean raw
	final double gmeanMCF_;     // mantissa-corrected geometric mean with CF

	// --- Final phase: hybrid estimates ---
	final double hybridDLLF;    // LC/Mean blend (for DLL types), plain (no history)
	final double hybridDDLF;    // LC/Mean/HMeanM blend (for DDL types), plain
	final double hybridDLLHistF; // History-corrected LC/Mean blend
	final double hybridDDLHistF; // History-corrected LC/Mean/HMeanM blend
	final double hybridPlus2F;  // SBS → Mean+H with DLC zone detection

}
