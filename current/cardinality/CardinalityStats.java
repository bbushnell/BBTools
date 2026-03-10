package cardinality;

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
	 * Full constructor with optional cardinality-indexed CF table and microIndex floor.
	 * @param matrixCard_  cardinality CF matrix (null = use occupancy-only)
	 * @param cardKeys_    load-factor key array for matrixCard_ (null when matrixCard_ is null)
	 * @param microIndex_  64-bit micro-index for low-cardinality estimation; 0 = disabled
	 */
	CardinalityStats(double difSum_, double hllSumFilled_, double hllSumFilledM_,
	                 double gSum_, int count_, int buckets_,
	                 LongList sortBuf_, float[][] cfMatrix_, int cfBuckets_,
	                 float[][] matrixCard_, float[] cardKeys_, long microIndex_){
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

		// Sort once here; all median/mwa/mean99 methods use the sorted buf
		sortBuf.sort();

		mean   =difSum/div;
		gmean  =Math.exp(gSum/div);
		median =Tools.max(1, sortBuf.median());
		mwa    =Tools.max(1.0, sortBuf.medianWeightedAverage());

		// Trimmed mean (mean99)
		final int trim=count/256;
		final int trimLow=Math.max(0, trim-V);
		double mean99Sum=0;
		final int mean99N=count-trimLow-trim;
		if(mean99N>0){for(int i=trim; i<count-trimLow; i++){mean99Sum+=sortBuf.get(i);}}
		mean99=(mean99N>0 ? mean99Sum/mean99N : mean);

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

		// CF-corrected values used by hybrid
		meanEstCF   =meanEst   *cf(CorrectionFactor.MEAN);
		hmeanPureMCF=hasMantissa ? hmeanPureM*cf(CorrectionFactor.HMEANM) : hmeanPureM;

		// MicroIndex low-cardinality estimate: LC over 64 virtual buckets
		final int microBits=(int)Long.bitCount(microIndex_);
		microEst=(microIndex_==0 ? 0 :
		          (int)(64*Math.log((double)64/Math.max(Math.min(63, 64-microBits), 0.5))));
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
	 * Returns the standard 10-element raw estimates array.
	 * The hybrid value is passed in by the caller (use hybridDLL() or hybridDDL()).
	 */
	public double[] toArray(double hybridEst){
		final double micro=CardinalityTracker.USE_MICRO ? microEst : 0;
		if(count==0){final double[] z=new double[11]; for(int i=0;i<10;i++){z[i]=micro;} z[10]=microEst; return z;}
		return new double[]{
			Math.max(meanEstCF,                              micro),
			Math.max(hmeanPure  *cf(CorrectionFactor.HMEAN), micro),
			Math.max(hasMantissa ? hmeanPureM*cf(CorrectionFactor.HMEANM) : hmeanPureM, micro),
			Math.max(gmeanEst   *cf(CorrectionFactor.GMEAN), micro),
			Math.max(hmeanEst   *cf(CorrectionFactor.HLL),   micro),
			Math.max(lcPure,                                 micro),
			Math.max(hybridEst,                              micro),
			Math.max(mwaEst     *cf(CorrectionFactor.MWA),   micro),
			Math.max(medianCorr *cf(CorrectionFactor.MEDCORR),micro),
			Math.max(mean99Est  *cf(CorrectionFactor.MEAN99), micro),
			microEst   // index 10: raw micro estimate, always computed regardless of USE_MICRO
		};
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private double cf(int type){
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

}
