package cardinality;

import shared.Tools;

/**
 * Static computation methods for cardinality estimation.
 * <p>
 * No instance state.  Every method takes explicit parameters and returns a result.
 * CardStats extends this and holds the results as fields; subclasses never touch
 * these methods directly — they call new CardStats(...) and let the constructor
 * invoke these.
 * <p>
 * Two families of methods:
 * <ul>
 *   <li><b>Counts-only</b>: LC, LCMin, DLC, HLL, Mean — derived from int[64] NLZ histogram.</li>
 *   <li><b>CF lookup</b>: cardinality-indexed and occupancy-indexed correction factor retrieval.</li>
 * </ul>
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public abstract class AbstractCardStats {

	/*--------------------------------------------------------------*/
	/*----------------    Sums from NLZ Counts      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Derives difSum, hllSumFilled, gSum, and filled count from an NLZ histogram.
	 * These are the fundamental bucket sums needed by Mean, HMean, GMean, and HLL.
	 * <p>
	 * difSum: sum of dif values.  dif(k) = 2^(63-k) for k>=1, Long.MAX_VALUE for k=0.
	 * hllSumFilled: sum of 2^(-k) for filled buckets (HLL indicator sum).
	 * gSum: sum of ln(dif) for filled buckets (geometric mean input).
	 * filled: total non-empty buckets.
	 *
	 * @param counts  int[64] NLZ histogram: counts[k] = number of buckets with absNlz==k
	 * @return double[4]: {difSum, hllSumFilled, gSum, filled}
	 */
	/**
	 * Derive sums from NLZ histogram.  counts is int[66]: [0]=empties, [k+1]=buckets at absNlz=k.
	 * Returns {difSum, hllSumFilled, gSum, filled}.
	 */
	static double[] sumsFromCounts(final int[] counts){
		double difSum=0, hllSumFilled=0, gSum=0;
		int filled=0;
		// Skip counts[0] (empties); iterate counts[1..65] where absNlz = i-1
		for(int i=1; i<counts.length; i++){
			final int n=counts[i];
			if(n==0){continue;}
			filled+=n;
			final int absNlz=i-1;
			final long dif=(absNlz==0 ? Long.MAX_VALUE : (absNlz<WORDLEN ? 1L<<(WORDLEN-absNlz-1) : 1L));
			difSum+=n*(double)dif;
			hllSumFilled+=n*Math.pow(2.0, -absNlz);
			gSum+=n*Math.log(Math.max(1, dif));
		}
		return new double[]{difSum, hllSumFilled, gSum, filled};
	}

	/*--------------------------------------------------------------*/
	/*----------------       LC Estimators          ----------------*/
	/*--------------------------------------------------------------*/

	/** Basic linear counting: B * ln(B / max(V, 0.5)).  V = empty buckets. */
	static double lcRaw(final int V, final int B){
		return (double)B*Math.log((double)B/Math.max(V, 0.5));
	}

	/**
	 * Tier-compensated LC: finds the lowest tier with empty-equivalent buckets
	 * and returns DLC at that tier.
	 * <p>
	 * Walk the NLZ histogram accumulating V_k = V + counts[0] + ... + counts[k-1].
	 * The first tier k where V_k > 0 yields lcMin = 2^k * B * ln(B / V_k).
	 * If V > 0 (tier 0 has empties), lcMin = lcRaw.
	 *
	 * @param counts  NLZ histogram
	 * @param V       empty buckets (numBuckets - filled)
	 * @param B       total buckets
	 * @param lcFloor minimum LC estimate (from history bits or micro index)
	 * @return tier-compensated LC estimate
	 */
	static double lcMin(final int[] counts, final int V, final int B, final double lcFloor){
		if(V>0){
			return Math.max((double)B*Math.log((double)B/Math.max(V, 0.5)), lcFloor);
		}
		// All buckets filled at tier 0; walk upward. counts[k+1] = buckets at absNlz=k.
		int vk=V;
		for(int k=0; k<counts.length-1; k++){
			vk+=counts[k+1]; // counts[k+1] = buckets at tier k
			if(vk>0){
				return Math.max((1L<<(k+1))*(double)B*Math.log((double)B/Math.max(vk, 0.5)), lcFloor);
			}
		}
		return lcFloor;
	}

	/*--------------------------------------------------------------*/
	/*----------------       DLC Estimators         ----------------*/
	/*--------------------------------------------------------------*/

	/** DLC estimate at a single tier: 2^tier * B * ln(B / max(V_k, 0.5)). */
	static double dlcFromVk(final int tier, final int vk, final int B){
		return (1L<<tier)*(double)B*Math.log((double)B/Math.max(vk, 0.5));
	}

	/** Returns the lowest tier (absNlz) that has filled buckets. counts[k+1]=buckets at absNlz=k. */
	static int lowestActiveTier(final int[] counts){
		for(int k=0; k<counts.length-1; k++){
			if(counts[k+1]>0){return k;}
		}
		return 0;
	}

	/**
	 * DLC log-space exponential blend (the primary DLC estimator).
	 * <p>
	 * Blends all tier estimates using exponential weights centered on target occupancy,
	 * averaged in log-space for stability.  Smooth transition to lcMin at low occupancy.
	 * <p>
	 * This is the "one good DLC method" — the others (dlcOriginal, dlc45_55, dlc5tier,
	 * dlcExpLinear) are dead.
	 *
	 * @param counts   NLZ histogram
	 * @param V        empty buckets
	 * @param B        total buckets
	 * @param lcMin    tier-compensated LC (fallback at low occupancy)
	 * @return DLC estimate
	 */
	static double dlcLogSpaceBlend(final int[] counts, final int V, final int B,
			final double lcMin){
		if(V>=DLC_BLEND_HI*B){return lcMin;}
		final double target=B*DLC_TARGET_FRAC;
		final double alpha=DLC_ALPHA/B;
		final int minVK=Tools.max(1, DLC_MIN_VK, (int)(B*DLC_MIN_VK_FRACTION));
		final int maxVK=B-minVK;
		final int startTier=lowestActiveTier(counts);
		int vk=V;
		double sumW=0, sumWLogE=0;
		for(int tier=startTier; tier<NUM_DLC_TIERS; tier++){
			if(vk>=minVK && vk<=maxVK){
				final double est=dlcFromVk(tier, vk, B);
				final double w=Math.exp(-alpha*Math.abs(vk-target));
				sumW+=w;
				sumWLogE+=w*Math.log(est);
			}
			if(tier+1<counts.length && tier<NUM_DLC_TIERS-1){vk+=counts[tier+1];}
			if(vk>=B){break;}
		}
		if(sumW<=0){return lcMin;}
		final double blendEst=Math.exp(sumWLogE/sumW);
		if(V>DLC_BLEND_LO*B){
			final double t=(V-DLC_BLEND_LO*B)/((DLC_BLEND_HI-DLC_BLEND_LO)*B);
			return t*lcMin+(1-t)*blendEst;
		}
		return blendEst;
	}

	/**
	 * DLC 3-region blend with target=0.20 occupancy.
	 * Uses same log-space exponential weighting but with 20% free target instead of 25%.
	 */
	static double dlcBlend3(final int[] counts, final int V, final int B,
			final double lcMin){
		if(V>=0.5*B){return lcMin;}
		final double target=B*0.20;
		final double alpha=DLC_ALPHA/B;
		final int startTier=lowestActiveTier(counts);
		int vk=V;
		double sumW=0, sumWLogE=0;
		for(int tier=startTier; tier<NUM_DLC_TIERS; tier++){
			if(vk>=1 && vk<B){
				final double est=dlcFromVk(tier, vk, B);
				final double w=Math.exp(-alpha*Math.abs(vk-target));
				sumW+=w;
				sumWLogE+=w*Math.log(est);
			}
			if(tier+1<counts.length && tier<NUM_DLC_TIERS-1){vk+=counts[tier+1];}
			if(vk>=B){break;}
		}
		if(sumW<=0){return lcMin;}
		final double blendEst=Math.exp(sumWLogE/sumW);
		if(V>0.3*B){
			final double t=(V-0.3*B)/(0.2*B);
			return t*lcMin+(1-t)*blendEst;
		}
		return blendEst;
	}

	/**
	 * DLC best: single tier closest to 25% free (75% occupancy).
	 * Only averages when two adjacent tiers are equidistant from target.
	 */
	static double dlcBest(final int[] counts, final int V, final int B,
			final double lcMin){
		if(V>0.4*B){return lcMin;}
		final double target=B*0.25;
		final int startTier=lowestActiveTier(counts);
		int vk=V;
		int bestTier=startTier;
		double bestDist=Math.abs(vk-target);
		final int[] vks=new int[NUM_DLC_TIERS];
		vks[startTier]=vk;
		for(int tier=startTier+1; tier<NUM_DLC_TIERS; tier++){
			if(tier<counts.length){vk+=counts[tier];} // counts[tier] = buckets at absNlz=tier-1
			vks[tier]=vk;
			final double dist=Math.abs(vk-target);
			if(dist<bestDist && vk>=1 && vk<B){bestDist=dist; bestTier=tier;}
			if(vk>=B){break;}
		}
		final double estBest=dlcFromVk(bestTier, vks[bestTier], B);
		// Average with adjacent tier if equidistant
		if(bestTier>0){
			final double distLo=Math.abs(vks[bestTier-1]-target);
			if(Math.abs(distLo-bestDist)<0.5){
				return 0.5*estBest+0.5*dlcFromVk(bestTier-1, vks[bestTier-1], B);
			}
		}
		if(bestTier+1<NUM_DLC_TIERS && vks[bestTier+1]<B){
			final double distHi=Math.abs(vks[bestTier+1]-target);
			if(Math.abs(distHi-bestDist)<0.5){
				return 0.5*estBest+0.5*dlcFromVk(bestTier+1, vks[bestTier+1], B);
			}
		}
		return estBest;
	}

	/**
	 * Information-weighted DLC blend with selectable error model.
	 * <p>
	 * Per-tier expected error was measured empirically (LL6, B=2048, 200K instances)
	 * and fit with several models:
	 * <ul>
	 *   <li><b>Mode 0 — 2-param empirical (R²=0.999):</b>
	 *       {@code E[err] = a/√n + b/√V} where n=tier occupancy, V=B−n</li>
	 *   <li><b>Mode 1 — 1-param delta-method (R²=0.997):</b>
	 *       {@code E[err] = 0.736·√(n/(B·V)) / ln(B/V)}</li>
	 *   <li><b>Mode 2 — 0-param pure theory (R²=0.983, default):</b>
	 *       {@code E[err] = √(2/π)·√(n/(B·V)) / ln(B/V)}
	 *       <br>Delta-method propagation of binomial V variance through LC.
	 *       Coefficient √(2/π) = E[|Z|] for standard normal.</li>
	 *   <li><b>Mode 3 — legacy Gaussian:</b>
	 *       {@code w = exp(-α·|V_k - target|)}, the original dlcLogSpaceBlend weighting.</li>
	 * </ul>
	 * For modes 0-2: weight = (1/E[err])^power, blended in log-space.
	 * At low occupancy (V > DLC_BLEND_HI·B), falls back to lcMin.
	 *
	 * @param counts  int[66] NLZ histogram
	 * @param V       empty buckets
	 * @param B       total buckets
	 * @param lcMin   tier-compensated LC (fallback)
	 * @param power   exponent applied to inverse expected error
	 * @return information-weighted DLC estimate
	 */
	static double dlcInfoPow(final int[] counts, final int V, final int B,
			final double lcMin, final float power){
		if(V>=DLC_BLEND_HI*B){return lcMin;}
		final int minVK=Tools.max(1, DLC_MIN_VK, (int)(B*DLC_MIN_VK_FRACTION));
		final int maxVK=B-minVK;
		final int startTier=lowestActiveTier(counts);
		final double Bd=(double)B;
		final double target=B*DLC_TARGET_FRAC;
		final double alpha=DLC_ALPHA/B;
		int vk=V;
		double sumW=0, sumWLogE=0;
		for(int tier=startTier; tier<NUM_DLC_TIERS; tier++){
			if(vk>=minVK && vk<=maxVK){
				final double est=dlcFromVk(tier, vk, B);
				final int occ=B-vk;
				final double w;
				switch(DLC_INFO_MODE){
					case 0:{
						final double err=DLC_INFO_A/Math.sqrt(Math.max(occ, 1))
							+DLC_INFO_B/Math.sqrt(Math.max(vk, 1));
						w=Math.pow(1.0/err, power);
						break;
					}
					case 1:{
						final double err=0.736*Math.sqrt((double)occ/(Bd*vk))/Math.log(Bd/vk);
						w=Math.pow(1.0/err, power);
						break;
					}
					case 3:{
						w=Math.exp(-alpha*Math.abs(vk-target));
						break;
					}
					default:{
						final double err=SQRT_2_OVER_PI*Math.sqrt((double)occ/(Bd*vk))/Math.log(Bd/vk);
						w=Math.pow(1.0/err, power);
						break;
					}
				}
				sumW+=w; sumWLogE+=w*Math.log(est);
			}
			if(tier+1<counts.length && tier<NUM_DLC_TIERS-1){vk+=counts[tier+1];}
			if(vk>=B){break;}
		}
		if(sumW<=0){return lcMin;}
		final double blendEst=Math.exp(sumWLogE/sumW);
		if(V>DLC_BLEND_LO*B){
			final double t=(V-DLC_BLEND_LO*B)/((DLC_BLEND_HI-DLC_BLEND_LO)*B);
			return t*lcMin+(1-t)*blendEst;
		}
		return blendEst;
	}

	static final double SQRT_2_OVER_PI=Math.sqrt(2.0/Math.PI);

	/** DLC estimate at a given tier, using V_k computed from counts[1..tier]. counts is int[66]. */
	static double dlcEstimate(final int tier, final int[] counts, final int V, final int B){
		if(tier==0){return lcRaw(V, B);}
		int vk=V;
		for(int i=0; i<tier && i+1<counts.length; i++){
			vk+=counts[i+1]; // counts[i+1] = buckets at absNlz=i
		}
		return (1L<<tier)*(double)B*Math.log((double)B/Math.max(vk, 0.5));
	}

	/*--------------------------------------------------------------*/
	/*----------------     HLL / Mean Estimators    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * HLL-style all-buckets estimate with LC fallback.
	 * Uses the standard alpha_m correction and small-range LC correction.
	 *
	 * @param hllSumFilled  sum of 2^(-absNlz) over filled buckets
	 * @param V             empty buckets
	 * @param B             total buckets
	 * @param alpha_m       HLL bias correction: 0.7213/(1+1.079/B)
	 * @return HLL cardinality estimate
	 */
	static double hllEstimate(final double hllSumFilled, final int V, final int B,
			final double alpha_m){
		final double hllSum=hllSumFilled+V; // empty buckets contribute 2^0 = 1
		final double raw=2*alpha_m*(double)B*(double)B/hllSum;
		if(raw<2.5*B){
			// Small-range correction: plain LC (V=0 means no empties, use 0.5 to avoid log(inf))
			return (double)B*Math.log((double)B/Math.max(V, 0.5));
		}
		return raw;
	}

	/**
	 * HLL with history-corrected constant.
	 * Identical to standard HLL but multiplies alpha_m by the terminal CF
	 * to account for changed bias from per-state corrections.
	 */
	static double hllHistory(final double hllSumFilled, final int V, final int B,
			final double alpha_m){
		final double alpha_hist=alpha_m*HLL_HIST_TERMINAL_CF;
		final double hllSum=hllSumFilled+V;
		final double raw=2*alpha_hist*(double)B*(double)B/hllSum;
		if(raw<2.5*B){
			return (double)B*Math.log((double)B/Math.max(V, 0.5));
		}
		return raw;
	}

	/**
	 * Raw Mean estimate: 2 * (MAX_VALUE / mean_dif) * filled * correction.
	 * <p>
	 * The correction factor (filled+B)/(2*B) compensates for empty-bucket bias.
	 * NOTE: uses float cast on (B+B) to match existing CardinalityStats exactly.
	 *
	 * @param difSum  sum of dif values over filled buckets
	 * @param filled  number of filled buckets
	 * @param B       total buckets
	 * @return raw (pre-CF) mean cardinality estimate
	 */
	static double meanEstimate(final double difSum, final int filled, final int B){
		final int div=Math.max(filled, 1);
		final double mean=difSum/div;
		// Cast to float matches CardinalityStats line 144 exactly
		final float correction=(filled+B)/(float)(B+B);
		return 2*(Long.MAX_VALUE/Tools.max(1.0, mean))*div*correction;
	}

	/**
	 * Raw GMean (geometric mean) estimate.
	 * Same formula as Mean but uses geometric mean of dif values.
	 */
	static double gmeanEstimate(final double gSum, final int filled, final int B){
		final int div=Math.max(filled, 1);
		final double gmean=Math.exp(gSum/div);
		final float correction=(filled+B)/(float)(B+B);
		return 2*(Long.MAX_VALUE/gmean)*div*correction;
	}

	/**
	 * Filled-bucket harmonic mean (HMean).
	 * Unlike HLL, this uses only filled buckets (no empty-bucket contribution).
	 *
	 * @param hllSumFilled  sum of 2^(-absNlz) over filled buckets
	 * @param filled        number of filled buckets
	 * @param alpha_m       HLL bias correction constant
	 * @return filled-bucket harmonic mean estimate
	 */
	static double hmeanEstimate(final double hllSumFilled, final int filled,
			final double alpha_m){
		if(filled==0){return 0;}
		return 2*alpha_m*(double)filled*(double)filled/hllSumFilled;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Hybrid Estimators       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Hybrid estimator for mantissa-free classes (DLL3, DLL4).
	 * Blends lcForHybrid → meanCF using log interpolation over [0.20B, 5.0B].
	 * Zone detection uses lcMin (uncorrected); the blended VALUE uses lcForHybrid
	 * (which may be lcHist or lcMin depending on history availability).
	 *
	 * @param lcForHybrid  LC value for blending (lcHist when available, else lcMin)
	 * @param lcMin        tier-compensated LC for zone detection
	 * @param meanCF       CF-corrected mean estimate
	 * @param B            total buckets
	 * @return hybrid cardinality estimate
	 */
	static double hybridDLL(final double lcForHybrid, final double lcMin,
			final double meanCF, final int B){
		final double hb0=0.20*B, hb1=5.0*B;
		if(lcMin<=hb0){return lcForHybrid;}
		if(lcMin<hb1){
			final double t=Math.log(lcMin/hb0)/Math.log(hb1/hb0);
			return (1-t)*lcForHybrid+t*meanCF;
		}
		return meanCF;
	}

	/**
	 * Hybrid estimator for mantissa-having classes (DDL, DDL2, DDL8).
	 * Blends lcForHybrid → meanCF → hmeanMCF using log interpolation with smooth crossover.
	 */
	static double hybridDDL(final double lcForHybrid, final double lcMin,
			final double meanCF, final double hmeanMCF, final int B){
		final double hb0=0.20*B, hbMid1=1.0*B, hbMid2=2.5*B, hb1=5.0*B;
		if(lcMin<=hb0){return lcForHybrid;}
		final double t=Math.log(lcMin/hb0)/Math.log(hb1/hb0);
		if(lcMin<=hbMid1){
			return (1-t)*lcForHybrid+t*meanCF;
		}else if(lcMin<=hbMid2){
			final double mix=(hbMid2-lcMin)/(hbMid2-hbMid1);
			final double blended=meanCF*mix+hmeanMCF*(1-mix);
			return (1-t)*lcForHybrid+t*blended;
		}else if(lcMin<=hb1){
			return (1-t)*lcForHybrid+t*hmeanMCF;
		}
		return hmeanMCF;
	}

	/**
	 * HybridDLC: blends DLC → Mean using same formula as hybridDLL.
	 * DLC replaces LC as both the value and the zone detector.
	 */
	static double hybridDLC(final double dlcRaw, final double meanEst, final int B){
		final double hb0=0.20*B, hb1=5.0*B;
		if(dlcRaw<=hb0){return dlcRaw;}
		if(dlcRaw<hb1){
			final double t=Math.log(dlcRaw/hb0)/Math.log(hb1/hb0);
			return (1-t)*dlcRaw+t*meanEst;
		}
		return meanEst;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Micro Estimator        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * MicroIndex low-cardinality estimate: LC over 64 virtual buckets.
	 * Returns 0 when microIndex is 0 (disabled).
	 */
	static int microEstimate(final int microBits){
		if(microBits==0){return 0;}
		return (int)(64*Math.log((double)64/Math.max(Math.min(63, 64-microBits), 0.5)));
	}

	/*--------------------------------------------------------------*/
	/*----------------   Legacy Output Helpers      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds the rawEstimates() output array matching the old CardinalityStats.toArray() layout.
	 * <p>
	 * This exists ONLY for migration: subclass rawEstimates() calls this to produce
	 * output identical to the old code.  Will be deleted once calibration drivers
	 * are refactored to use named getters.
	 *
	 * @param s          the CardStats object with all computed estimates
	 * @param hybridEst  the hybrid estimate chosen by the subclass (hybridDLL or hybridDDL)
	 * @return double[] in the exact same column order as old CardinalityStats.toArray()
	 */
	static double[] buildLegacyArray(final CardStats s, final double hybridEst){
		final int total=17+NUM_DLC_TIERS;
		if(s.filled==0){
			final double[] z=new double[total];
			final double micro=s.microEst;
			for(int i=0; i<10; i++){z[i]=micro;}
			z[10]=micro;
			return z;
		}
		final double[] r=new double[total];
		// MicroIndex benefit flows through V (empty buckets) into LC/DLC/Hybrid.
		// No per-estimator Math.max(est, micro) floors — those distort calibration.
		r[0] =s.meanMCF();           // Mean
		r[1] =s.hmeanCF;             // HMean
		r[2] =s.hmeanMCF;            // HMeanM
		r[3] =s.gmeanMCF();          // GMean
		r[4] =s.hllCF;               // HLL
		r[5] =s.lcRawF;              // LC
		r[6] =hybridEst;             // Hybrid
		r[7] =s.hybDLC50CF;          // HybDLC50
		r[8] =s.dThHybF;             // DThHyb
		r[9] =s.dlcLowestF;          // LCmin
		r[10]=s.hllHistoryF;          // HLL with history-corrected constant
		r[11]=s.dlcCF;               // DLC
		r[12]=s.dlc3bCF;             // DLC3B
		r[13]=s.dlcBestCF;           // DLCBest
		r[14]=s.hybDLCcf;            // HybDLC
		r[15]=s.lcHistF;             // LCHist
		r[16]=s.lcHistMultF;          // LCHMult
		for(int t=0; t<NUM_DLC_TIERS; t++){
			r[17+t]=s.dlcTier(t);
		}
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Static Configuration    ----------------*/
	/*--------------------------------------------------------------*/

	/** Hash word length — 64 bits for long keys. */
	static final int WORDLEN=64;

	/** Number of DLC tiers in output. */
	public static final int NUM_DLC_TIERS=64;

	/** Exponential decay constant for DLC log-space blending. alpha = DLC_ALPHA / B. */
	public static float DLC_ALPHA=9.0f;
	/** Target occupancy fraction for DLC blending (0.25 = 75% full). */
	public static float DLC_TARGET_FRAC=0.25f;
	/** Below this V/B fraction, DLC blend is pure; above, transitions to lcMin. */
	public static float DLC_BLEND_LO=0.3f;
	/** Above this V/B fraction, DLC returns lcMin (too few filled buckets). */
	public static float DLC_BLEND_HI=0.5f;
	/** Minimum V_k fraction for tier to participate in DLC blend. */
	public static float DLC_MIN_VK_FRACTION=0.002f;
	/** Absolute minimum V_k for tier participation. */
	public static int DLC_MIN_VK=2;

	/** Coefficient for left arm of DLC information weight: a/sqrt(occ). */
	public static float DLC_INFO_A=0.719f;
	/** Coefficient for right arm of DLC information weight: b/sqrt(V). */
	public static float DLC_INFO_B=0.094f;
	/** Power to raise the information weight: w = (1/expectedErr)^power. */
	public static float DLC_INFO_POWER=4.5f;
	/** Power for HC (history counting) tier weighting. HC tiers have fewer effective
	 *  buckets than DLC tiers, so a lower power avoids over-concentrating on one tier. */
	public static float HC_INFO_POWER=2.0f;
	/** Error model mode: 0=2-param empirical, 1=1-param fitted, 2=0-param theory (default). */
	public static int DLC_INFO_MODE=2;

	/** LDLC double-blend zone multipliers (as fractions of B).
	 *  Blend A (LCHist→DLC): [LDLC_A_LO*B, LDLC_A_HI*B]
	 *  Blend B (HC ramps in): [LDLC_B_LO*B, LDLC_B_HI*B] */
	public static float LDLC_A_LO=1.0f, LDLC_A_HI=4.0f;
	public static float LDLC_B_LO=2.0f, LDLC_B_HI=5.0f;

	/** When false, microIndex does not adjust V for LC calculation. */
	public static boolean USE_MICRO_FOR_LC=true;

	/** When false, all microIndex values are treated as zero regardless of what callers pass.
	 *  Useful for comparing classes that populate microIndex against LL6 which does not. */
	public static boolean USE_MICRO_INDEX=true;

	/** Terminal HMean CF for history-corrected HLL. Default 1.0 = standard HLL. */
	public static double HLL_HIST_TERMINAL_CF=1.0;

	/** Max iterations for iterative CF refinement. */
	public static int DEFAULT_CF_ITERS=2;
	/** Convergence threshold for iterative CF. */
	public static double DEFAULT_CF_DIF=1e-6;
	/** Minimum DLC seed as multiple of B for iterative CF (below this, single lookup). */
	public static float MIN_SEED_CF_MULT=10.0f;

}
