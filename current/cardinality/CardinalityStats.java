package cardinality;

import java.util.Arrays;

import shared.Tools;
import structures.LongList;

/**
 * Messenger object holding pre-computed bucket statistics for cardinality estimation.
 * <p>
 * Created by a per-subclass {@code summarize()} call that scans the bucket array once,
 * accumulating the sums needed for all estimators. The subclass passes those sums here;
 * this class pre-computes all derived quantities in the constructor so that estimator
 * methods are simple field reads or trivial arithmetic.
 * <p>
 * Two hybrid formulae are provided:
 * <ul>
 *   <li>{@link #hybridDLL()} — mantissa-free (DLL3, DLL4): LC → Mean, log interpolation.</li>
 *   <li>{@link #hybridDDL()} — mantissa (DDL, DDL2, DDL8): LC → Mean → HMeanM blend.</li>
 * </ul>
 * The subclass chooses the appropriate one:
 * <pre>
 *   CardinalityStats s = summarize();
 *   return s.toArray(s.hybridDLL());   // DLL subclasses
 *   return s.toArray(s.hybridDDL());   // DDL subclasses
 * </pre>
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public class CardinalityStats {

	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a CardinalityStats from pre-accumulated bucket sums.
	 *
	 * @param difSum_         sum of dif values (mantissa-corrected) for filled buckets
	 * @param hllSumFilled_   sum of 2^(-absNlz) for filled buckets (integer NLZ)
	 * @param hllSumFilledM_  sum of 2^(-absNlz + fractional) for filled buckets (mantissa NLZ);
	 *                        pass hllSumFilled_ for mantissa-free classes
	 * @param gSum_           sum of log(dif) for filled buckets
	 * @param count_          number of filled buckets
	 * @param buckets_        total number of buckets
	 * @param sortBuf_        LongList of dif values; will be sorted in-place here
	 * @param cfMatrix_       per-class correction factor matrix
	 * @param cfBuckets_      bucket count used when building cfMatrix_
	 */
//	/** Backward-compatible constructor; uses occupancy-only CF, no microIndex. */
//	CardinalityStats(double difSum_, double hllSumFilled_, double hllSumFilledM_,
//	                 double gSum_, int count_, int buckets_,
//	                 LongList sortBuf_, float[][] cfMatrix_, int cfBuckets_){
//		this(difSum_, hllSumFilled_, hllSumFilledM_, gSum_, count_, buckets_,
//		     sortBuf_, cfMatrix_, cfBuckets_, null, null, 0);
//	}

//	/** Constructor with microIndex; uses occupancy-only CF (no cardinality table). */
//	CardinalityStats(double difSum_, double hllSumFilled_, double hllSumFilledM_,
//	                 double gSum_, int count_, int buckets_,
//	                 LongList sortBuf_, float[][] cfMatrix_, int cfBuckets_, long microIndex_){
//		this(difSum_, hllSumFilled_, hllSumFilledM_, gSum_, count_, buckets_,
//		     sortBuf_, cfMatrix_, cfBuckets_, null, null, microIndex_);
//	}

//	/** Backward-compatible constructor with cardinality table; microIndex defaults to 0. */
//	CardinalityStats(double difSum_, double hllSumFilled_, double hllSumFilledM_,
//	                 double gSum_, int count_, int buckets_,
//	                 LongList sortBuf_, float[][] cfMatrix_, int cfBuckets_,
//	                 float[][] matrixCard_, float[] cardKeys_){
//		this(difSum_, hllSumFilled_, hllSumFilledM_, gSum_, count_, buckets_,
//		     sortBuf_, cfMatrix_, cfBuckets_, matrixCard_, cardKeys_, 0);
//	}

	/**
	 * Full constructor with optional cardinality-indexed CF table, microIndex floor, and NLZ histogram.
	 * @param matrixCard_  cardinality CF matrix (null = use occupancy-only)
	 * @param cardKeys_    load-factor key array for matrixCard_ (null when matrixCard_ is null)
	 * @param microIndex_  64-bit micro-index for low-cardinality estimation; 0 = disabled
	 * @param nlzCounts_   absolute NLZ histogram (int[64]); null disables DLC output
	 */
	CardinalityStats(double difSum_, double hllSumFilled_, double hllSumFilledM_,
	                 double gSum_, int count_, int buckets_,
	                 LongList sortBuf_, float[][] cfMatrix_, int cfBuckets_,
	                 float[][] matrixCard_, float[] cardKeys_, long microIndex_,
	                 int[] nlzCounts_){
		// Set card data first so cf() calls below can use three-domain lookup.
		matrixCard=matrixCard_;
		cardKeys=cardKeys_;
		difSum=difSum_;
		hllSumFilled=hllSumFilled_;
		hllSumFilledM=hllSumFilledM_;
		gSum=gSum_;
		count=count_;
		buckets=buckets_;
		sortBuf=sortBuf_;
		cfMatrix=cfMatrix_;
		cfBuckets=cfBuckets_;

		V=buckets-count;
		div=Tools.max(count, 1);
		alpha_m=0.7213/(1.0+1.079/buckets);
		correction=(count+buckets)/(float)(buckets+buckets);

		mean   =difSum/div;
		gmean  =Math.exp(gSum/div);

		// Sort-dependent estimators: only computed when sortBuf is provided (USE_SORTBUF=true)
		if(sortBuf!=null){
			sortBuf.sort();
			median =Tools.max(1, sortBuf.median());
			mwa    =Tools.max(1.0, sortBuf.medianWeightedAverage());
			final int trim=count/256;
			final int trimLow=Math.max(0, trim-V);
			double mean99Sum=0;
			final int mean99N=count-trimLow-trim;
			if(mean99N>0){for(int i=trim; i<count-trimLow; i++){mean99Sum+=sortBuf.get(i);}}
			mean99=(mean99N>0 ? mean99Sum/mean99N : mean);
		}else{
			median=1;
			mwa=1;
			mean99=mean;
		}

		// HLL all-buckets sum: filled contribute 2^(-absNlz), empty contribute 1.0
		hllSum=hllSumFilled+V;

		// HLL-style estimate with LC fallback at low occupancy
		final double hmeanRaw=2*alpha_m*(double)buckets*(double)buckets/hllSum;
		double hmeanEstTmp=hmeanRaw;
		if(hmeanEstTmp<2.5*buckets && V>0){hmeanEstTmp=(double)buckets*Math.log((double)buckets/V);}
		hmeanEst=hmeanEstTmp;

		// Filled-bucket harmonic mean estimates
		hmeanPure =(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilled);
		hmeanPureM=(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilledM);

		// Raw cardinality estimates (before CF)
		lcPure    =buckets*Math.log((double)buckets/Math.max(V, 0.5));
		meanEst   =2*(Long.MAX_VALUE/Tools.max(1.0, mean))   *div*correction;
		gmeanEst  =2*(Long.MAX_VALUE/gmean)                  *div*correction;
		mwaEst    =2*(Long.MAX_VALUE/mwa)                    *div*correction;
		medianCorr=2*(Long.MAX_VALUE/(double)median)         *div*correction;
		mean99Est =2*(Long.MAX_VALUE/Tools.max(1.0, mean99)) *div*correction;

		// True for mantissa classes (DDL, DDL8); false for DLL3, DLL4 (hllSumFilledM == hllSumFilled)
		hasMantissa=(hllSumFilledM_!=hllSumFilled_);

		// MicroIndex low-cardinality estimate: LC over 64 virtual buckets
		final int microBits=(int)Long.bitCount(microIndex_);
		microEst=(microIndex_==0 ? 0 :
		          (int)(64*Math.log((double)64/Math.max(Math.min(63, 64-microBits), 0.5))));

		nlzCounts=nlzCounts_;

		// Cache DLC estimates for CF lookup (CF-free, so no circular dependency)
		// MUST be computed before cf() calls below
		dlc3bEst=dlcBlend3();
		dlcEst=dlcLogSpace025();

		// CF-corrected values used by hybrid
		meanEstCF   =meanEst   *cf(CorrectionFactor.MEAN);
		hmeanPureMCF=hasMantissa ? hmeanPureM*cf(CorrectionFactor.HMEANM) : hmeanPureM;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Estimators           ----------------*/
	/*--------------------------------------------------------------*/

	/** Mean cardinality estimate with CF correction. */
	public double mean()   {return meanEstCF;}

	/** HMean (filled-bucket harmonic mean) with optional CF. */
	public double hmean(boolean applyCF){
		return applyCF ? hmeanPure*cf(CorrectionFactor.HMEAN) : hmeanPure;
	}

	/** HMeanM (mantissa-corrected harmonic mean) with optional CF. */
	public double hmeanM(boolean applyCF){
		return applyCF ? hmeanPureM*cf(CorrectionFactor.HMEANM) : hmeanPureM;
	}

	/** GMean cardinality estimate with optional CF. */
	public double gmean(boolean applyCF){
		return applyCF ? gmeanEst*cf(CorrectionFactor.GMEAN) : gmeanEst;
	}

	/** HLL-style all-buckets estimate with optional CF. */
	public double hll(boolean applyCF){
		return applyCF ? hmeanEst*cf(CorrectionFactor.HLL) : hmeanEst;
	}

	/** LC (linear counting) estimate; no CF applied. */
	public double lc(){return lcPure;}

	/** MWA estimate with optional CF. */
	public double mwa(boolean applyCF){
		return applyCF ? mwaEst*cf(CorrectionFactor.MWA) : mwaEst;
	}

	/** MedianCorr estimate with optional CF. */
	public double medianCorr(boolean applyCF){
		return applyCF ? medianCorr*cf(CorrectionFactor.MEDCORR) : medianCorr;
	}

	/** Mean99 (trimmed mean) estimate with optional CF. */
	public double mean99(boolean applyCF){
		return applyCF ? mean99Est*cf(CorrectionFactor.MEAN99) : mean99Est;
	}

	/**
	 * Hybrid estimator for mantissa-free classes (DLL3, DLL4).
	 * Blends LC → Mean using log interpolation over [hb0, hb1].
	 */
	public double hybridDLL(){
		final double hb0=0.20*buckets, hb1=5.0*buckets;
		if(lcPure<=hb0){return lcPure;}
		if(lcPure<hb1){
			final double t=Math.log(lcPure/hb0)/Math.log(hb1/hb0);
			return (1-t)*lcPure+t*meanEstCF;
		}
		return meanEstCF;
	}

	/**
	 * Extended hybrid for DLL3 blending in HMeanM to suppress the high-cardinality hump.
	 * <p>
	 * NOTE: This method was developed to compensate for a hump in DLL3's Mean estimate at
	 * high cardinality (7x–15x buckets). The hump turned out to be caused by a bug:
	 * DLL3's summarize() was not passing matrixCard/cardKeys to CardinalityStats, so the
	 * cardinality-indexed CF table was never applied. Once that was fixed, the hump
	 * disappeared and this blend became unnecessary — HMeanM dives negative at high
	 * cardinality, making things worse. DLL3 now calls hybridDLL() instead.
	 *
	 * @deprecated Superseded by hybridDLL() once cardinality CF wiring was corrected.
	 */
	@Deprecated
	public double hybridDLL3(){
		// Phase 1: LC → Mean (same as hybridDLL)
		final double hb0=0.20*buckets, hb1=5.0*buckets;
		if(lcPure<=hb0){return lcPure;}
		if(lcPure<hb1){
			final double t=Math.log(lcPure/hb0)/Math.log(hb1/hb0);
			return (1-t)*lcPure+t*meanEstCF;
		}
		// Phase 2+: zone detection via meanEst*~CF (CF-independent, no circular dep)
		final double zone=meanEst*0.7;
		final double hb2=5.5*buckets, hb3=7.0*buckets, hb4mid=15.0*buckets, hb4=55.0*buckets;
		if(zone<=hb2){return meanEstCF;}
		final double blended=0.1*meanEstCF+0.9*hmeanPureMCF;
		if(zone<=hb3){
			final double u=Math.log(zone/hb2)/Math.log(hb3/hb2);
			return (1-u)*meanEstCF+u*blended;
		}
		if(zone<=hb4mid){return blended;}
		if(zone<=hb4){
			final double u=(zone-hb4mid)/(hb4-hb4mid);
			return (1-u)*blended+u*meanEstCF;
		}
		return meanEstCF;
	}

	/**
	 * Hybrid estimator blending DLC3B → Mean using the same formula as hybridDLL().
	 * DLC3B replaces LC as the low-cardinality anchor; Mean(CF) anchors the high end.
	 * Uses DLC3B estimate as both the value and the zone detector.
	 */
	public double hybridDLC(){
		final double dlc=dlcLogSpace025();
		final double hb0=0.20*buckets, hb1=5.0*buckets;
		if(dlc<=hb0){return dlc;}
		if(dlc<hb1){
			final double t=Math.log(dlc/hb0)/Math.log(hb1/hb0);
			return (1-t)*dlc+t*meanEst;
		}
		return meanEst;
	}

	/** Like hybridDLC but blends to 50/50 DLC+Mean instead of pure Mean at high cardinality. */
	private double hybDLC50(){
		final double dlc=dlcLogSpace025();
		final double high=0.5*dlc+0.5*meanEst;
		final double hb0=0.20*buckets, hb1=5.0*buckets;
		if(dlc<=hb0){return dlc;}
		if(dlc<hb1){
			final double t=Math.log(dlc/hb0)/Math.log(hb1/hb0);
			return (1-t)*dlc+t*high;
		}
		return high;
	}

	/**
	 * Hybrid estimator for mantissa-having classes (DDL, DDL2, DDL8).
	 * Blends LC → Mean → HMeanM using log interpolation with smooth Mean/HMeanM crossover.
	 */
	public double hybridDDL(){
		final double hb0=0.20*buckets, hbMid1=1.0*buckets, hbMid2=2.5*buckets, hb1=5.0*buckets;
		if(lcPure<=hb0){return lcPure;}
		final double t=Math.log(lcPure/hb0)/Math.log(hb1/hb0);
		if(lcPure<=hbMid1){
			return (1-t)*lcPure+t*meanEstCF;
		}else if(lcPure<=hbMid2){
			final double mix=(hbMid2-lcPure)/(hbMid2-hbMid1);
			final double blended=meanEstCF*mix+hmeanPureMCF*(1-mix);
			return (1-t)*lcPure+t*blended;
		}else if(lcPure<=hb1){
			return (1-t)*lcPure+t*hmeanPureMCF;
		}
		return hmeanPureMCF;
	}

	/**
	 * Returns microIndex-based cardinality floor, or 0 if USE_MICRO is false or microIndex was 0.
	 * Used by cardinality() as: card = Math.max(card, s.microCardinality()).
	 */
	public long microCardinality(){return CardinalityTracker.USE_MICRO ? microEst : 0;}

	/**
	 * Returns the raw estimates array: 11 standard estimators + 1 DLC combined + NUM_DLC_TIERS DLC tiers.
	 * The hybrid value is passed in by the caller (use hybridDLL() or hybridDDL()).
	 */
	public double[] toArray(double hybridEst){
		final int total=11+4+NUM_DLC_TIERS;
		final double micro=CardinalityTracker.USE_MICRO ? microEst : 0;
		if(count==0){final double[] z=new double[total]; for(int i=0;i<10;i++){z[i]=micro;} z[10]=microEst; return z;}
		final double[] r=new double[total];
		r[0] =Math.max(meanEstCF,                              micro);
		r[1] =Math.max(hmeanPure  *cf(CorrectionFactor.HMEAN), micro);
		r[2] =Math.max(hasMantissa ? hmeanPureM*cf(CorrectionFactor.HMEANM) : hmeanPureM, micro);
		r[3] =Math.max(gmeanEst   *cf(CorrectionFactor.GMEAN), micro);
		r[4] =Math.max(hmeanEst   *cf(CorrectionFactor.HLL),   micro);
		r[5] =Math.max(lcPure,                                 micro);
		r[6] =Math.max(hybridEst,                              micro);
		r[7] =hybDLC50()*cf(CorrectionFactor.HYBDLC50); // slot MWA → HybDLC50 with own CF
		r[8] =0; // disabled
		r[9] =0; // disabled
		r[10]=0; // disabled
		r[11]=dlcLogSpace025()*cf(CorrectionFactor.DLC); // DLC with CF
		r[12]=dlcBlend3();      // slot DLC3B → Chloe (log-space, target=0.20) ← NEW
		r[13]=dlcBest();
		r[14]=hybridDLC()*cf(CorrectionFactor.HYBDLC);
		for(int t=0; t<NUM_DLC_TIERS; t++){
			r[15+t]=dlcEstimate(t);
		}
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------       DLC Estimators         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * DynamicLC estimate for a given tier.
	 * <p>
	 * DLC(k) = 2^k * B * ln(B / max(V_k, 0.5))
	 * where V_k = (empty buckets) + sum(nlzCounts[0..k-1]).
	 * <p>
	 * DLC(0) = classic LC (uses empty buckets only).
	 * DLC(k) for k≥1 additionally counts buckets with absoluteNlz < k.
	 * Higher tiers remain accurate at higher cardinalities where lower tiers lose resolution.
	 *
	 * @param tier 0 = classic LC; k = extends to buckets with absoluteNlz < k
	 * @return cardinality estimate, or 0 if nlzCounts is null and tier > 0
	 */
	public double dlcEstimate(int tier){
		if(tier==0){return lcPure;}
		if(nlzCounts==null){return 0;}
		int vk=V; // empty buckets
		for(int i=0; i<tier && i<nlzCounts.length; i++){
			vk+=nlzCounts[i];
		}
		return (1L<<tier)*(double)buckets*Math.log((double)buckets/Math.max(vk, 0.5));
	}

	/**
	 * Combined DLC estimate: exponential-weighted blend of all nontrivial DLC tier estimates.
	 * <p>
	 * Each tier k is most accurate when its V_k (empty + buckets with absNlz &lt; k)
	 * is near B/2 (50% occupancy for that tier). Tiers are weighted by proximity to B/2
	 * with exponential decay: weight = 2^(-|V_k - B/2| / scale).
	 * <p>
	 * Below 25% DLC0 occupancy (V_0 &gt; 0.75*B), returns classic LC.
	 * Only tiers with V_k in [1, B-1] (nontrivially informative) participate.
	 */
	public double dlcCombined(){
		if(nlzCounts==null){return lcPure;}
		if(V>0.75*buckets){return lcPure;}
		final double half=buckets*0.5;
		final double scale=buckets*0.125;
		double sumW=0, sumWE=0;
		int vk=V;
		for(int tier=0; tier<NUM_DLC_TIERS; tier++){
			if(vk>=1 && vk<buckets){
				final double est=(1L<<tier)*(double)buckets*Math.log((double)buckets/Math.max(vk, 0.5));
				final double w=Math.pow(2, -Math.abs(vk-half)/scale);
				sumW+=w;
				sumWE+=w*est;
			}
			if(tier<nlzCounts.length){vk+=nlzCounts[tier];}
		}
		return sumW>0 ? sumWE/sumW : lcPure;
	}

	/**
	 * Combined DLC estimate: 3-region blend targeting 75% occupancy (25% free).
	 * <p>
	 * R0: Below 60% DLC0 occupancy — pure DLC0 (= LC).
	 * R1: DLC0 is still closest to target — blend 60% DLC0 + 40% DLC1.
	 * R2: DLC1+ becomes closest to target — 60% center + 40% split between flanks,
	 *     linearly interpolated by relative closeness to target occupancy.
	 */
	public double dlcBlend3(){
		if(nlzCounts==null){return lcPure;}
		// Smooth transition: above 50% free, pure LC; 30-50% free, blend with LC
		if(V>=0.5*buckets){return lcPure;}
		final double target=buckets*0.20; // 80% full = 20% free
		final double alpha=DLC_ALPHA/buckets; // alpha*B ≈ 9.2 (calibrated at 1024 buckets)
		int vk=V;
		double sumW=0, sumWLogE=0;
		for(int tier=0; tier<NUM_DLC_TIERS; tier++){
			if(vk>=1 && vk<buckets){
				final double est=dlcFromVk(tier, vk);
				final double w=Math.exp(-alpha*Math.abs(vk-target));
				sumW+=w;
				sumWLogE+=w*Math.log(est);
			}
			if(tier<nlzCounts.length && tier<NUM_DLC_TIERS-1){vk+=nlzCounts[tier];}
			if(vk>=buckets){break;}
		}
		if(sumW<=0){return lcPure;}
		final double blendEst=Math.exp(sumWLogE/sumW);
		// Smooth transition in the 30-50% free zone
		if(V>0.3*buckets){
			final double t=(V-0.3*buckets)/(0.2*buckets);
			return t*lcPure+(1-t)*blendEst;
		}
		return blendEst;
	}

	// ===== Comparison variants for plotting (temporary) =====

	/** Variant 0: Original 60/40 3-tier blend, target=0.25 */
	private double dlcOriginal(){
		if(nlzCounts==null){return lcPure;}
		if(V>0.4*buckets){return lcPure;}
		final double target=buckets*0.25;
		int vk=V; int bestTier=0; double bestDist=Math.abs(vk-target); int maxUseful=0;
		final int[] vks=new int[NUM_DLC_TIERS]; vks[0]=vk;
		for(int tier=1; tier<NUM_DLC_TIERS; tier++){
			if(tier-1<nlzCounts.length){vk+=nlzCounts[tier-1];}
			vks[tier]=vk; if(vk<buckets){maxUseful=tier;}
			double dist=Math.abs(vk-target);
			if(dist<bestDist && vk>=1 && vk<buckets){bestDist=dist; bestTier=tier;}
			if(vk>=buckets){break;}
		}
		final int mid=bestTier; final double estMid=dlcFromVk(mid, vks[mid]);
		if(mid==0){
			if(maxUseful>=1 && vks[1]<buckets){return 0.6*estMid+0.4*dlcFromVk(1, vks[1]);}
			return estMid;
		}
		boolean hasLo=(mid>0 && vks[mid-1]>=1), hasHi=(mid<maxUseful && vks[mid+1]<buckets);
		if(!hasLo && !hasHi){return estMid;}
		if(!hasLo){return 0.6*estMid+0.4*dlcFromVk(mid+1, vks[mid+1]);}
		if(!hasHi){return 0.6*estMid+0.4*dlcFromVk(mid-1, vks[mid-1]);}
		double distLo=Math.abs(vks[mid-1]-target), distHi=Math.abs(vks[mid+1]-target);
		double totalDist=distLo+distHi;
		return 0.6*estMid+(0.4*distHi/totalDist)*dlcFromVk(mid-1, vks[mid-1])+(0.4*distLo/totalDist)*dlcFromVk(mid+1, vks[mid+1]);
	}

	/** Variant 1: Agent1 — 45/55 center/flank weights, target=0.25 */
	private double dlc45_55(){
		if(nlzCounts==null){return lcPure;}
		if(V>0.4*buckets){return lcPure;}
		final double target=buckets*0.25;
		int vk=V; int bestTier=0; double bestDist=Math.abs(vk-target); int maxUseful=0;
		final int[] vks=new int[NUM_DLC_TIERS]; vks[0]=vk;
		for(int tier=1; tier<NUM_DLC_TIERS; tier++){
			if(tier-1<nlzCounts.length){vk+=nlzCounts[tier-1];}
			vks[tier]=vk; if(vk<buckets){maxUseful=tier;}
			double dist=Math.abs(vk-target);
			if(dist<bestDist && vk>=1 && vk<buckets){bestDist=dist; bestTier=tier;}
			if(vk>=buckets){break;}
		}
		final int mid=bestTier; final double estMid=dlcFromVk(mid, vks[mid]);
		if(mid==0){
			if(maxUseful>=1 && vks[1]<buckets){return 0.45*estMid+0.55*dlcFromVk(1, vks[1]);}
			return estMid;
		}
		boolean hasLo=(mid>0 && vks[mid-1]>=1), hasHi=(mid<maxUseful && vks[mid+1]<buckets);
		if(!hasLo && !hasHi){return estMid;}
		if(!hasLo){return 0.45*estMid+0.55*dlcFromVk(mid+1, vks[mid+1]);}
		if(!hasHi){return 0.45*estMid+0.55*dlcFromVk(mid-1, vks[mid-1]);}
		double distLo=Math.abs(vks[mid-1]-target), distHi=Math.abs(vks[mid+1]-target);
		double totalDist=distLo+distHi;
		return 0.45*estMid+(0.55*distHi/totalDist)*dlcFromVk(mid-1, vks[mid-1])+(0.55*distLo/totalDist)*dlcFromVk(mid+1, vks[mid+1]);
	}

	/** Variant 2: Agent4 — 5-tier inverse distance weighting, target=0.25 */
	private double dlc5tier(){
		if(nlzCounts==null){return lcPure;}
		if(V>0.4*buckets){return lcPure;}
		final double target=buckets*0.25;
		int vk=V; int bestTier=0; double bestDist=Math.abs(vk-target); int maxUseful=0;
		final int[] vks=new int[NUM_DLC_TIERS]; vks[0]=vk;
		for(int tier=1; tier<NUM_DLC_TIERS; tier++){
			if(tier-1<nlzCounts.length){vk+=nlzCounts[tier-1];}
			vks[tier]=vk; if(vk<buckets){maxUseful=tier;}
			double dist=Math.abs(vk-target);
			if(dist<bestDist && vk>=1 && vk<buckets){bestDist=dist; bestTier=tier;}
			if(vk>=buckets){break;}
		}
		// Use center + up to 2 on each side, weighted by 1/(dist+1)
		final int lo=Math.max(0, bestTier-2), hi=Math.min(maxUseful, bestTier+2);
		double sumW=0, sumWE=0;
		for(int t=lo; t<=hi; t++){
			if(vks[t]>=1 && vks[t]<buckets){
				double w=1.0/(Math.abs(vks[t]-target)+1.0);
				sumW+=w; sumWE+=w*dlcFromVk(t, vks[t]);
			}
		}
		return sumW>0 ? sumWE/sumW : lcPure;
	}

	/** Variant 3: Agent5 — all-tier exponential, LINEAR average, target=0.25 */
	private double dlcExpLinear(){
		if(nlzCounts==null){return lcPure;}
		if(V>0.4*buckets){return lcPure;}
		final double target=buckets*0.25;
		final double alpha=DLC_ALPHA/buckets;
		int vk=V; double sumW=0, sumWE=0;
		for(int tier=0; tier<NUM_DLC_TIERS; tier++){
			if(vk>=1 && vk<buckets){
				double est=dlcFromVk(tier, vk);
				double w=Math.exp(-alpha*Math.abs(vk-target));
				sumW+=w; sumWE+=w*est;
			}
			if(tier<nlzCounts.length && tier<NUM_DLC_TIERS-1){vk+=nlzCounts[tier];}
			if(vk>=buckets){break;}
		}
		return sumW>0 ? sumWE/sumW : lcPure;
	}

	/** Variant 4: Chloe — all-tier exponential, LOG-SPACE average, target=0.25, smooth transition */
	private double dlcLogSpace025(){
		if(nlzCounts==null){return lcPure;}
		if(V>=0.5*buckets){return lcPure;}
		final double target=buckets*0.25;
		final double alpha=DLC_ALPHA/buckets;
		int vk=V; double sumW=0, sumWLogE=0;
		for(int tier=0; tier<NUM_DLC_TIERS; tier++){
			if(vk>=1 && vk<buckets){
				double est=dlcFromVk(tier, vk);
				double w=Math.exp(-alpha*Math.abs(vk-target));
				sumW+=w; sumWLogE+=w*Math.log(est);
			}
			if(tier<nlzCounts.length && tier<NUM_DLC_TIERS-1){vk+=nlzCounts[tier];}
			if(vk>=buckets){break;}
		}
		if(sumW<=0){return lcPure;}
		double blendEst=Math.exp(sumWLogE/sumW);
		if(V>0.3*buckets){double t=(V-0.3*buckets)/(0.2*buckets); return t*lcPure+(1-t)*blendEst;}
		return blendEst;
	}

	// dlcBlend3() above is Variant 5: same as dlcLogSpace025 but target=0.20

	/**
	 * Combined DLC estimate: strictly uses the single tier closest to 75% occupancy (25% free).
	 * Only averages when two adjacent tiers are equidistant from target.
	 * Below 60% DLC0 occupancy, returns pure DLC0 (= LC).
	 */
	public double dlcBest(){
		if(nlzCounts==null){return lcPure;}
		if(V>0.4*buckets){return lcPure;}
		final double target=buckets*0.25; // 75% full = 25% free
		int vk=V;
		int bestTier=0;
		double bestDist=Math.abs(vk-target);
		final int[] vks=new int[NUM_DLC_TIERS];
		vks[0]=vk;
		for(int tier=1; tier<NUM_DLC_TIERS; tier++){
			if(tier-1<nlzCounts.length){vk+=nlzCounts[tier-1];}
			vks[tier]=vk;
			final double dist=Math.abs(vk-target);
			if(dist<bestDist && vk>=1 && vk<buckets){bestDist=dist; bestTier=tier;}
			if(vk>=buckets){break;}
		}
		final double estBest=dlcFromVk(bestTier, vks[bestTier]);
		if(bestTier>0){
			final double distLo=Math.abs(vks[bestTier-1]-target);
			if(Math.abs(distLo-bestDist)<0.5){return 0.5*estBest+0.5*dlcFromVk(bestTier-1, vks[bestTier-1]);}
		}
		if(bestTier+1<NUM_DLC_TIERS && vks[bestTier+1]<buckets){
			final double distHi=Math.abs(vks[bestTier+1]-target);
			if(Math.abs(distHi-bestDist)<0.5){return 0.5*estBest+0.5*dlcFromVk(bestTier+1, vks[bestTier+1]);}
		}
		return estBest;
	}

	/** DLC estimate from pre-computed V_k. */
	private double dlcFromVk(int tier, int vk){
		return (1L<<tier)*(double)buckets*Math.log((double)buckets/Math.max(vk, 0.5));
	}

	/** DLC with parameterized sine model correction.
	 *  @param ampScale  amplitude multiplier (1.0=full, 0.5=half). Applied as cf=(cf_full-1)*ampScale+1.
	 *  @param phaseShift  additional phase offset in radians (0=default, negative=left, positive=right). */
	private double dlcSineCorr(double ampScale, double phaseShift){
		final double dlc=dlcLogSpace025();
		if(dlc<=0){return dlc;}
		final double log2est=Math.log(dlc/buckets)/Math.log(2);
		final double amp=0.070234/Math.sqrt(buckets);
		final double sineErr=amp*Math.sin(2*Math.PI*log2est-0.176443+phaseShift)+0.000627;
		final double cfFull=1.0/(1.0+sineErr);
		final double cf=(cfFull-1.0)*ampScale+1.0; // scale deviation from 1: cf=(cfFull+1)/2 when ampScale=0.5
		return dlc*cf;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private double cf(int type){
		if(CorrectionFactor.tableVersion>=3){
			return CorrectionFactor.getV1CF(dlcEst, type);
		}
		if(CorrectionFactor.tableVersion>=1){
			return CorrectionFactor.getV1CF(dlc3bEst, type);
		}
		if(matrixCard!=null && CardinalityTracker.USE_CARD_CF && CorrectionFactor.USE_CORRECTION){
			return CorrectionFactor.getCF(cfMatrix, cfBuckets, matrixCard, cardKeys,
			                              count, buckets, meanEst, type);
		}
		return CorrectionFactor.getCF(cfMatrix, cfBuckets, count, buckets, type);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	// Inputs from bucket scan
	final double difSum;
	final double hllSumFilled;
	final double hllSumFilledM;
	final double gSum;
	final int count;
	final int buckets;
	final LongList sortBuf;
	final float[][] cfMatrix;
	final int cfBuckets;
	final float[][] matrixCard; // cardinality-indexed CF table; null = occupancy-only
	final float[] cardKeys;     // MeanEst keys for matrixCard binary search

	// Derived counts
	final int V;
	final int div;
	final double alpha_m;
	final double correction;
	final double hllSum;

	// Derived statistics
	final double mean;
	final double gmean;
	final double mean99;
	final double mwa;
	final long median;
	final double hmeanEst;    // HLL all-buckets estimate
	final double hmeanPure;   // filled-bucket HMean
	final double hmeanPureM;  // filled-bucket HMeanM (mantissa-corrected)
	final double lcPure;

	// Raw cardinality estimates (before CF)
	final double meanEst;
	final double gmeanEst;
	final double mwaEst;
	final double medianCorr;
	final double mean99Est;

	// True for mantissa classes (DDL, DDL8); false for mantissa-free classes (DLL3, DLL4)
	final boolean hasMantissa;

	// CF-corrected values pre-computed for hybrid blend
	final double meanEstCF;
	final double hmeanPureMCF;

	// MicroIndex-derived estimate; 0 if microIndex was 0 or not applicable
	final int microEst;

	// Absolute NLZ histogram from summarize(); null for legacy classes
	final int[] nlzCounts;

	// Cached DLC estimates (CF-free) for CF table lookup
	final double dlc3bEst;  // v1 key
	final double dlcEst;    // v3 key (dlcLogSpace025)

	/** Exponential decay constant for DLC log-space blending. alpha = DLC_ALPHA / buckets.
	 *  Controls how quickly tier weights decay away from the target occupancy.
	 *  Set via dlcalpha= parameter in DDLCalibrationDriver. */
	public static float DLC_ALPHA=9.216f;

	/** Number of DLC tiers included in toArray() output (DLC_0 through DLC_{N-1}). */
	public static final int NUM_DLC_TIERS=64;

}
