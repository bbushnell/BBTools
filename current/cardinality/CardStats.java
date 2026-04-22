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
		this(buckets_, counts_, nlzBits_, histBits_, luckBits_, mantissaBits_,
			numBuckets_, microIndex_, added_, cfMatrix_, cfBuckets_, mantissaOff_,
			Integer.MAX_VALUE, null, 1f, 1f, 1.0);
	}

	/**
	 * 12-base constructor plus per-class terminal Mean bias correction.
	 * Used by trackers that override terminalMeanCF()/terminalMeanPlusCF().
	 */
	CardStats(final char[] buckets_, final int[] counts_,
			final int nlzBits_, final int histBits_, final int luckBits_, final int mantissaBits_,
			final int numBuckets_, final long microIndex_, final long added_,
			final float[][] cfMatrix_, final int cfBuckets_,
			final double mantissaOff_,
			final float terminalMeanCF_, final float terminalMeanPlusCF_){
		this(buckets_, counts_, nlzBits_, histBits_, luckBits_, mantissaBits_,
			numBuckets_, microIndex_, added_, cfMatrix_, cfBuckets_, mantissaOff_,
			Integer.MAX_VALUE, null, terminalMeanCF_, terminalMeanPlusCF_, 1.0);
	}

	/**
	 * Constructor with overflow tier correction for DLC.
	 * @param nativeRelTiers_  number of native relative NLZ levels (3 for DLL2, 7 for DLL3); Integer.MAX_VALUE = no correction
	 * @param tierErrCoeffs_   polynomial coefficients for signed error vs log2(B/V_k), or null
	 */
	CardStats(final char[] buckets_, final int[] counts_,
			final int nlzBits_, final int histBits_, final int luckBits_, final int mantissaBits_,
			final int numBuckets_, final long microIndex_, final long added_,
			final float[][] cfMatrix_, final int cfBuckets_,
			final double mantissaOff_,
			final int nativeRelTiers_, final double[] tierErrCoeffs_){
		this(buckets_, counts_, nlzBits_, histBits_, luckBits_, mantissaBits_,
			numBuckets_, microIndex_, added_, cfMatrix_, cfBuckets_, mantissaOff_,
			nativeRelTiers_, tierErrCoeffs_, 1f, 1f, 1.0);
	}

	/**
	 * Full constructor with per-class terminal Mean bias correction.
	 * @param terminalMeanCF_      asymptotic meanRaw/trueCard ratio for plain Mean (1 = no correction)
	 * @param terminalMeanPlusCF_  asymptotic ratio for Mean+H (1 = no correction; only matters if histBits>0)
	 * @param tierScale_           DLC tier scale factor (1.0 standard, 1.5 compressed, 2.0 dual-hash)
	 */
	CardStats(final char[] buckets_, final int[] counts_,
			final int nlzBits_, final int histBits_, final int luckBits_, final int mantissaBits_,
			final int numBuckets_, final long microIndex_, final long added_,
			final float[][] cfMatrix_, final int cfBuckets_,
			final double mantissaOff_,
			final int nativeRelTiers_, final double[] tierErrCoeffs_,
			final float terminalMeanCF_, final float terminalMeanPlusCF_,
			final double tierScale_){

		assert(counts_!=null && counts_.length==66);
		assert((histBits_|luckBits_|mantissaBits_)==0 || buckets_!=null) :
			"Per-bucket array required when sub-NLZ bits are present";

		// Override resolution: static AbstractCardStats.OVERRIDE_* wins when > 0.
		// Set tmcf=1 tmpcf=1 on the command line during preliminary CF generation
		// so per-class bias corrections are disabled and raw ratios can be measured.
		final float tmCF=(OVERRIDE_TERMINAL_MEAN_CF>0 ? OVERRIDE_TERMINAL_MEAN_CF : terminalMeanCF_);
		final float tmPlusCF=(OVERRIDE_TERMINAL_MEAN_PLUS_CF>0 ? OVERRIDE_TERMINAL_MEAN_PLUS_CF : terminalMeanPlusCF_);

		/*--- Store configuration ---*/
		tierScale=tierScale_;
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
		nativeRelTiers=nativeRelTiers_;
		tierErrCoeffs=tierErrCoeffs_;

		/*--- Phase 1: counts-only estimates (all classes) ---*/

		// Derive fundamental sums from NLZ histogram (counts[0]=empties, counts[k+1]=absNlz k)
		final double[] sums=sumsFromCounts(counts);
		difSum=sums[0];
		hllSumFilled=sums[1];
		gSum=sums[2];
		filled=(int)sums[3];
		assert filled+counts[0]==numBuckets :
			"counts sum ("+filled+"+"+counts[0]+") != numBuckets ("+numBuckets+")";

		// Empty buckets from NLZ counts only (no microIndex adjustment to V)
		V=numBuckets-filled;

		// Overflow compensation for LC in io=t mode.
		// For DLL2 (nativeRelTiers=3): NLZ 0-2 are stored, NLZ >= 3 is dropped.
		// The count at absNlz=2 (the last stored tier) estimates the missing overflow:
		//   T3+T4+... ≈ T2 (geometric series sum).
		// Collision-adjusted: missing_fills = T2 * (V / B).
		final int lcV;
		if(nativeRelTiers<Integer.MAX_VALUE && V>0){
			// counts[k+1] = buckets at absNlz=k. The overflow boundary tier is nativeRelTiers-1.
			// Use the count at that tier as the overflow estimate.
			final int overflowBoundaryTier=nativeRelTiers-1; // absNlz=2 for DLL2
			final int tBoundary=(overflowBoundaryTier+1<counts.length) ? counts[overflowBoundaryTier+1] : 0;
			final double missingFills=tBoundary*((double)V/numBuckets);
			lcV=Math.max(0, (int)Math.round(V-missingFills));
		}else{
			lcV=V;
		}

		// Bias correction constants
		alpha_m=0.7213/(1.0+1.079/numBuckets);
		// NOTE: float cast matches CardinalityStats line 144 exactly — do not change to double
		correction=(filled+numBuckets)/(float)(numBuckets+numBuckets);

		// History-based cardinality floor: filled buckets + total set history bits.
		// Each filled bucket saw >=1 element; each set history bit saw >=1 additional.
		{
			int hbCount=0;
			if(histBits_>0 && buckets_!=null){
				final int hm=(1<<histBits_)-1;
				for(int i=0; i<numBuckets; i++){
					if(buckets_[i]!=0){hbCount+=Integer.bitCount(buckets_[i]&hm);}
				}
			}
			historyFloor=filled+hbCount;
		}

		// MicroIndex estimate: LC over 64 virtual buckets (calculated but not used for flooring)
		microEst=microEstimate(microBits);

		// 1. LC — pure linear counting, floored by microIndex and history bit count
		//    Uses lcV (overflow-compensated empties) when applicable.
		lcNoMicroF=lcRaw(lcV, numBuckets);
		lcRawF=Math.max(Math.max(lcNoMicroF, microBits), historyFloor);

		// 2. LCMin — tier-compensated LC, floored by microIndex and history bit count
		lcMinF=Math.max(Math.max(lcMin(counts, lcV, numBuckets, 0), microBits), historyFloor);

		// 3. DLC — primary CF-free estimator; info-power weighted blend (mode 2 = 0-param theory)
		//    dlcPure is the raw tier blend; dlcRaw blends it with lcMin at low occupancy.
		//    Must be set before any CF lookups since dlcRaw is used as the cardinality seed.
		dlcPureF=dlcPure(counts, V, numBuckets, DLC_INFO_POWER);
		dlcRawF=dlcBlendWithLcMin(dlcPureF, lcMinF, V, numBuckets, DLC_BLEND_LO, DLC_BLEND_HI);

		// 4. DLC variants
		dlcBestF=dlcBest(counts, V, numBuckets, lcMinF, nativeRelTiers, tierErrCoeffs);
		dlc3bF=dlcBlend3(counts, V, numBuckets, lcMinF);
		dlcLowestF=lcMinF; // tier-compensated LC at actual floor; no separate computation needed

		// 5. HLL — all-buckets harmonic mean with LC fallback
		hllRawF=hllEstimate(hllSumFilled, V, numBuckets, alpha_m);

		// 6. HLL with history-corrected constant
		hllHistoryF=hllHistory(hllSumFilled, V, numBuckets, alpha_m);

		// 7. Mean, GMean, HMean — raw (pre-CF) estimates from counts-derived sums
		//    Mean gets the class-specific terminal bias multiplied out so downstream
		//    CF tables can converge to 1.0 at high cardinality.
		meanRawF=meanEstimate(difSum, filled, numBuckets)*tmCF;
		gmeanRawF=gmeanEstimate(gSum, filled, numBuckets);
		hmeanRawF=hmeanEstimate(hllSumFilled, filled, alpha_m);

		// 8. CF-corrected estimates — use DLC raw as cardinality seed for table lookup
		//    Grab immutable CF snapshot once; all lookups use this local reference
		cfd=CorrectionFactor.cfData;
		final double keyScale=(cfd!=null && cfd.buckets>0 ?
				(double)cfd.buckets/numBuckets : 1.0);
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
			// Use per-class terminalMeanPlusCF_ (set by the tracker subclass) to
			// fold the Mean+H asymptotic bias into tierMult. The old global
			// StateTable.terminalCF(histBits_, 0) lookup is replaced by this
			// class-specific value so regenerated CF tables converge to 1.0.
			final double invTermCF=1.0/tmPlusCF;

			// Precompute tierMult for all (nlzBin, histPattern) combos to avoid per-bucket Math.pow
			final int maxNlzBin=histBits_+2; // nlzBin is clamped to histBits_+1
			final int numPatterns=1<<histBits_;
			final double[][] tierMultTable=new double[maxNlzBin][numPatterns];
			for(int nb=0; nb<maxNlzBin; nb++){
				for(int hp=0; hp<numPatterns; hp++){
					final double cf=StateTable.historyOffset(nb, histBits_, hp);
					tierMultTable[nb][hp]=Math.pow(2.0, -(cf+StateTable.CF_OFFSET))*invTermCF;
				}
			}

			double cDif=0, cGSum=0;
			sbsStatesArr=new byte[numBuckets];
			// HC accumulators: per-NLZ bucket count and per-history-bit set counts
			final int[] nlzBucketCount=new int[64];
			final int[][] nlzHbitSet=new int[histBits_][64];
			int bucketsWithHistory=0;
			for(int i=0; i<numBuckets; i++){
				final int val=buckets_[i];
				if(val==0){
					sbsStatesArr[i]=-1; // empty bucket
					continue;
				}
				final int absNlz=(val>>>histBits_)-1;
				final int histPattern=val&histMask;
				if(histPattern!=0){bucketsWithHistory++;}
				final int nlzBin=Math.min(absNlz, histBits_+1);

				// Per-bucket dif: scaled by tierScale for compressed-tier variants
				final double dif=(absNlz==0 ? (double)Long.MAX_VALUE
					: Math.pow(2.0, 63-absNlz*tierScale));

				// Per-state correction for Mean/GMean only (precomputed table lookup)
				final double tierMult=tierMultTable[nlzBin][histPattern];

				cDif+=dif*tierMult;
				cGSum+=Math.log(Math.max(1, dif*tierMult));

				// LC history state index for sbs computation
				sbsStatesArr[i]=(byte)CorrectionFactor.sbsStateIndex(nlzBin, histPattern, histBits_);

				// HC: accumulate per-NLZ bucket counts and per-bit history set counts
				if(absNlz>=0 && absNlz<64){
					nlzBucketCount[absNlz]++;
					for(int d=0; d<histBits_; d++){
						if((histPattern&(1<<(histBits_-1-d)))!=0){nlzHbitSet[d][absNlz]++;}
					}
				}
			}
			corrDifSum=cDif;
			corrGSum=cGSum;

			// HC: history-only per-tier exact LC with info-power weighting
			int maxNlzHC=0;
			for(int t=63; t>=0; t--){
				if(nlzBucketCount[t]>0){ maxNlzHC=t; break; }
			}
			final int maxTierHC=Math.min(maxNlzHC+2, 63);
			double hcSumW=0, hcSumWLogE=0;
			for(int t=0; t<=maxTierHC; t++){
				int hcBeff=0, hcUnseen=0;
				for(int d=0; d<histBits_; d++){
					final int sourceTier=t+d+1;
					if(sourceTier<64){
						hcBeff+=nlzBucketCount[sourceTier];
						hcUnseen+=(nlzBucketCount[sourceTier]-nlzHbitSet[d][sourceTier]);
					}
				}
				if(hcBeff>=8 && hcUnseen>=1 && hcUnseen<hcBeff){
					final double tierRatio=Math.pow(2.0, tierScale);
					final double est=Math.pow(2.0, (t+1)*tierScale)/(tierRatio-1.0)*(double)numBuckets*Math.log((double)hcBeff/hcUnseen);
					if(est>0 && !Double.isNaN(est)){
						final int hcOcc=hcBeff-hcUnseen;
						final double hcErr=SQRT_2_OVER_PI
							*Math.sqrt((double)hcOcc/((double)hcBeff*hcUnseen))
							/Math.log((double)hcBeff/hcUnseen);
						final double w=Math.pow(1.0/hcErr, HC_INFO_POWER);
						hcSumW+=w;
						hcSumWLogE+=w*Math.log(est);
					}
				}
			}
			final double hcRaw=(hcSumW>0 ? Math.exp(hcSumWLogE/hcSumW) : 0);
			hcF=(hcRaw>0 ? hcRaw*CorrectionFactor.hcCfFormula(hcRaw)*HC_SCALE : 0);
			historyCoverage=(filled>0 ? (double)bucketsWithHistory/filled : 0);
		}else{
			hcF=0;
			historyCoverage=0;
		}
		hasHistoryCorrection=hasHistCorr;

		// History-corrected Mean, GMean
		if(hasHistCorr){
			final int divH=Math.max(filled, 1);
			final double meanH=corrDifSum/divH;
			meanCorrRawF=2*(Long.MAX_VALUE/Tools.max(1.0, meanH))*divH*correction;
			meanCorrCF=meanCorrRawF*cardCF(meanCorrRawF, CorrectionFactor.MEANH, keyScale);
			final double gmeanH=Math.exp(corrGSum/divH);
			gmeanCorrRawF=2*(Long.MAX_VALUE/gmeanH)*divH*correction;
			gmeanCorrCF=gmeanCorrRawF*cardCF(gmeanCorrRawF, CorrectionFactor.GMEAN, keyScale);
		}else{
			meanCorrRawF=meanRawF;
			meanCorrCF=meanCF;
			gmeanCorrRawF=gmeanRawF;
			gmeanCorrCF=gmeanCF;
		}

		// SBS: history-aware LC from per-bucket state indices, floored by history floor
		final double sbsRaw=computeSbs(sbsStatesArr, filled, numBuckets, lcNoMicroF);
		sbsF=Math.max(sbsRaw, Math.max(microBits, historyFloor));
		sbsNoMicroF=sbsRaw;
		sbsMultF=computeSbsMult(sbsStatesArr, filled, numBuckets, lcNoMicroF);

		// If sbs table is available and enabled, use sbs for hybrid blending
		if(CardinalityTracker.USE_SBS_IN_HYBRID
				&& CorrectionFactor.SBS_CF_TABLE!=null
				&& histBits_>0 && sbsStatesArr!=null){
			lcForHybridTmp=sbsF;
		}
		lcForHybridF=lcForHybridTmp;

		// DSBS: DLC with sbs as fallback, blended in cardinality space using dlcRaw
		dlcSbsF=dlcBlendWithSbs(dlcPureF, sbsF, dlcRawF, numBuckets, DLCSBS_BLEND_LO, DLCSBS_BLEND_HI);

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
				if(val==0){continue;} // empty or floor-level bucket
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
			meanMRawF=2*(Long.MAX_VALUE/Tools.max(1.0, meanM))*divM*correctionM*tmCF;
			meanMCF_=meanMRawF*cardCF(meanMRawF, CorrectionFactor.MEANM, keyScale);
			final double gmeanM=Math.exp(gSumM/divM);
			gmeanMRawF=2*(Long.MAX_VALUE/gmeanM)*divM*correctionM;
			gmeanMCF_=gmeanMRawF*cardCF(gmeanMRawF, CorrectionFactor.GMEAN, keyScale);
		}else{
			// No mantissa: HMeanM = HMean with HMEAN CF (not HMEANM, which is trained on mantissa data).
			hllSumFilledM=hllSumFilled;
			difSumM_=difSum;
			gSumM_=gSum;
			hmeanMRawF=hmeanRawF;
			hmeanMCF=hmeanCF; // fall back to HMEAN-corrected, not raw
			meanMRawF=meanRawF;
			meanMCF_=meanCF;
			gmeanMRawF=gmeanRawF;
			gmeanMCF_=gmeanCF;
		}

		// hasMantissa: true when mantissa-corrected sums differ from plain sums
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
		hybridPlus2F=hybridDLL(lcForHybridF, dlcRawF, meanForHybridH, numBuckets,
				HYBRID2_BLEND_LO, HYBRID2_BLEND_HI);

		/*--- LDLC: DlcSbs blended with HC ---*/
		{
			final double B=(double)numBuckets;
			final double bLo=LDLC_B_LO*B, bHi=LDLC_B_HI*B;
			final double maxHcWeight=CardinalityTracker.LDLC_HC_WEIGHT;
			final boolean hcUsable=(hcF>0 && dlcSbsF>0);
			if(dlcSbsF<=bLo || !hcUsable){
				ldlcF=dlcSbsF;
			}else if(dlcSbsF<=bHi){
				final double t=(dlcSbsF-bLo)/(bHi-bLo);
				final double hcW=t*maxHcWeight;
				ldlcF=(1-hcW)*dlcSbsF+hcW*hcF;
			}else{
				ldlcF=(1-maxHcWeight)*dlcSbsF+maxHcWeight*hcF;
			}
		}
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
		if(!CorrectionFactor.USE_CORRECTION){return 1;}
		// Formula mode for Mean CF: bypasses table entirely
		if((CorrectionFactor.USE_MEAN_CF_FORMULA || CorrectionFactor.USE_FORMULAS) && type==CorrectionFactor.MEAN){
			return CorrectionFactor.meanCfFormula(dlcRawF, cfd!=null ? cfd.meanCoeffs : CorrectionFactor.meanCfCoeffs);
		}
		// Formula mode for MeanH CF (history-blended mean): bypasses table entirely
		if((CorrectionFactor.USE_MEAN_CF_FORMULA || CorrectionFactor.USE_FORMULAS)
				&& type==CorrectionFactor.MEANH){
			final double[] mc=(cfd!=null ? cfd.meanhCoeffs : CorrectionFactor.meanhCfCoeffs);
			if(mc!=null){return CorrectionFactor.meanCfFormula(dlcRawF, mc);}
		}
		// Formula mode for HMeanM CF: bypasses table entirely
		if((CorrectionFactor.USE_HMEANM_CF_FORMULA || CorrectionFactor.USE_FORMULAS)
				&& type==CorrectionFactor.HMEANM){
			final double[] hmc=(cfd!=null ? cfd.hmeanmCoeffs : CorrectionFactor.hmeanmCfCoeffs);
			if(hmc!=null){return CorrectionFactor.meanCfFormula(dlcRawF, hmc);}
		}
		// v5 table: use dlcRawF as seed (the primary DLC estimator)
		final float[][] mat=(cfd!=null ? cfd.matrix : CorrectionFactor.v1Matrix);
		final float[] keys=(cfd!=null ? cfd.keys : CorrectionFactor.v1Keys);
		final int cfBkt=(cfd!=null ? cfd.buckets : CorrectionFactor.v1Buckets);
		if(mat==null || type==CorrectionFactor.LINEAR || type>=mat.length){return 1;}
		final int iters=(dlcRawF*keyScale>MIN_SEED_CF_MULT*cfBkt ? DEFAULT_CF_ITERS : 1);
		return CorrectionFactor.getCF(dlcRawF, rawEst, mat[type], keys,
				iters, DEFAULT_CF_DIF, keyScale);
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
		final int cfBkt=(cfd!=null ? cfd.buckets : CorrectionFactor.v1Buckets);
		final double ks=(cfBkt>0 ? (double)cfBkt/numBuckets : 1.0);
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
	/** LDLC: DlcSbs blended with HC. */
	public double ldlc(){return ldlcF;}
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

		// Count actually contributing buckets (si>=0); floor-level entries have si=-1
		// and are skipped in the sum, so row selection must match.
		int contributing=0;
		for(int i=0; i<states.length; i++){
			if(states[i]>=0){contributing++;}
		}
		if(contributing==0){return lcRaw;}

		// Table mode (preferred when available; formula override skips table)
		final float[][] table=CorrectionFactor.USE_SBS_FORMULA ? null : CorrectionFactor.SBS_CF_TABLE;
		if(table!=null){
			final int tableBuckets=CorrectionFactor.sbsBuckets;
			final float[] row;
			if(numBuckets==tableBuckets){
				row=table[Math.min(contributing, tableBuckets)];
			}else{
				final double pos=(double)contributing*tableBuckets/numBuckets;
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
				if(si>=0 && si<row.length){est+=row[si];}
			}
			return est;
		}

		// Formula fallback: SBS must always use table or formula, never bare lcRaw
		{
			double est=0;
			for(int i=0; i<states.length; i++){
				final int si=states[i];
				if(si>=0 && si<CorrectionFactor.SBS_STATES){est+=CorrectionFactor.sbsFormula(si, contributing, numBuckets);}
			}
			return est;
		}
	}

	/** History-aware LC using multipliers: lcRaw × weightedAvg(CF per bucket state).
	 *  Falls back to lcRaw when table or states are unavailable. */
	private static double computeSbsMult(final byte[] states, final int filled,
			final int numBuckets, final double lcRaw){
		final float[][] table=CorrectionFactor.SBS_MULT_TABLE;
		if(table==null || states==null){return lcRaw;}

		// Count contributing buckets (same floor-level fix as computeSbs)
		int contributing=0;
		for(int i=0; i<states.length; i++){
			if(states[i]>=0){contributing++;}
		}
		if(contributing==0){return lcRaw;}

		final int tableBuckets=CorrectionFactor.sbsMultBuckets;
		final float[] row;
		if(numBuckets==tableBuckets){
			row=table[Math.min(contributing, tableBuckets)];
		}else{
			final double pos=(double)contributing*tableBuckets/numBuckets;
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
			if(si>=0 && si<row.length){sumCF+=row[si]; n++;}
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
	/** HC estimate (CF-corrected, history-only per-tier LC). 0 when no history. */
	public double hc(){return hcF;}
	/** Fraction of filled buckets with any history bit set. 0 when no history. */
	public double historyCoverage(){return historyCoverage;}

	/*--------------------------------------------------------------*/
	/*----------------         Fields               ----------------*/
	/*--------------------------------------------------------------*/

	// --- CF table snapshot (grabbed once in constructor for thread safety) ---
	private final CorrectionFactor.CFTableData cfd;

	// --- Configuration (set once in constructor) ---
	final double tierScale;     // DLC tier scale (1.0 standard, 1.5 compressed, 2.0 dual-hash)
	final int numBuckets;       // total bucket count
	final int microBits;        // Long.bitCount(microIndex)
	final long added;           // total elements added
	final int nlzBits;          // NLZ bit width per bucket value
	final int nativeRelTiers;   // native relative NLZ levels (3 for DLL2); Integer.MAX_VALUE = no correction
	final double[] tierErrCoeffs; // polynomial coefficients for tier signed error, or null
	final int histBits;         // history bit width (0=none)
	final int luckBits;         // luck bit width (0=none)
	final int mantissaBits;     // mantissa bit width (0=none)
	final int[] counts;         // NLZ histogram (int[64])
	final float[][] cfMatrix;   //TODO: Unused — remove if no downstream consumer needs post-construction access to the CF matrix
	final int cfBuckets;        //TODO: Unused — remove along with cfMatrix
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
	final double lcNoMicroF;    // pure LC without microIndex adjustment
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
	final int historyFloor;     // lower bound: filledBuckets + totalSetHistoryBits (0 when no history)
	final boolean hasHistoryCorrection; // true when per-state history corrections applied
	final double meanCorrRawF;  // history-corrected Mean raw (before CF)
	final double meanCorrCF;    // history-corrected Mean with CF
	final double gmeanCorrRawF; // history-corrected GMean raw
	final double gmeanCorrCF;   // history-corrected GMean with CF
	final double sbsF;       // history-aware LC
	final double sbsNoMicroF; // SBS without microIndex
	final double sbsMultF;   // history-aware LC using multipliers
	final double lcForHybridF;  // LC value used for hybrid blending (sbs or lcMin)
	final double hcF;        // HC estimate (CF-corrected, 0 when no history)
	final double historyCoverage; // fraction of filled buckets with any history bit set

	// --- Phase 2b: mantissa estimates (defaults to non-mantissa equivalents) ---
	final double hllSumFilledM; // mantissa-corrected harmonic indicator sum
	final double difSumM_;      //TODO: Unused — stored for diagnostics but never read; remove or expose via getter
	final double gSumM_;        //TODO: Unused — same as difSumM_
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
	final double ldlcF;         // DlcSbs blended with HC

}
