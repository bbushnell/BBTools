package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DDL8: Dynamic DemiLog with 8-bit byte registers, parameterized by mantissa bits.
 * <p>
 * Each bucket: 1 byte, split as (exponentBits | mantissaBits) where exponentBits = 8 - mantissaBits.
 * Upper exponentBits hold (relNlz+1) clamped [0, maxRelNlzStored]; lower mantissaBits hold inverted mantissa.
 * Encoding: 0x00 = empty; non-empty = ((relNlz+1) << mantissaBits) | invMantissa.
 * <p>
 * mantissaBits=2: 6-bit exponent (64 tiers), 2-bit invMantissa.
 *   Used for LL8 (Mean estimator, ignores mantissa) and HLL8 (HLL estimator, ignores mantissa).
 *   Both get full 6-bit relative NLZ coverage; the 2 mantissa bits store precision for DDL8 use.
 * mantissaBits=3: 5-bit exponent (32 tiers), 3-bit invMantissa.
 * mantissaBits=4: 4-bit exponent (16 tiers), 4-bit invMantissa. Same range as DLL4.
 * <p>
 * All three variants produce the same 9-estimator output from rawEstimates():
 *   0=Mean(LL8), 1=HMean, 2=HMeanM(DDL8), 3=GMean, 4=HLL(HLL8), 5=LC, 6=Hybrid, 7=MWA, 8=MedianCorr
 * The LL8/HLL8/DDL8 distinction is just which column the caller examines.
 * <p>
 * Inverted mantissa convention (matches DLL2): higher stored invMantissa = higher fractional NLZ
 * = the hash is closer to the lower boundary of the NLZ tier = better sub-tier precision.
 * This makes the numeric byte comparison work directly for the update condition.
 * <p>
 * HMeanM fractional NLZ formula (DLL2-compatible):
 *   fractNlz = absNlz - 0.5 + invMantissa/mantissaScale
 *   hllSumFilledM += 2^(-(s>>>mantissaBits) + 1.5 - (s&mantissaMask)/mantissaScale - minZeros)
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicDemiLog8 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DynamicDemiLog8(){
		this(2048, 31, -1, 0);
	}

	DynamicDemiLog8(Parser p){
		super(p);
		maxArray=new byte[buckets];
		minZeroCount=buckets;
	}

	DynamicDemiLog8(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
		minZeroCount=buckets;
	}

	@Override
	public DynamicDemiLog8 copy(){return new DynamicDemiLog8(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns unsigned byte value for bucket i (0=empty, minNonEmpty-255=non-empty). */
	private int readBucket(final int i){return maxArray[i]&0xFF;}

	/** Writes byte value for bucket i. val must be in [0,255]. */
	private void writeBucket(final int i, final int val){maxArray[i]=(byte)val;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final double[] est=rawEstimates();
		final long card=Math.min(added, (long)est[6]); // Hybrid
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicDemiLog8)log);
	}

	/**
	 * Merges another DDL8 into this one (must have same mantissaBits).
	 * Converts each bucket to the new common minZeros frame, takes per-bucket max, recounts.
	 */
	public void add(DynamicDemiLog8 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			final int deltaA=(minZeros-newMinZeros)*minNonEmpty;
			final int deltaB=(log.minZeros-newMinZeros)*minNonEmpty;
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				writeBucket(i, Math.max(adjustStored(sA, deltaA), adjustStored(sB, deltaB)));
			}
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int s=readBucket(i);
				if(s>0){filledBuckets++;}
				if((s>>>mantissaBits)<2){minZeroCount++;} // empty (0) or tier-0 (relNlzStored=1)
			}
		}
	}

	/**
	 * Shifts a stored byte by delta (a multiple of minNonEmpty), clamping to valid range.
	 * delta is always <= 0 since newMinZeros >= oldMinZeros.
	 * Below-floor non-empty values clamp to tier-0 (minNonEmpty) rather than becoming empty.
	 */
	private int adjustStored(final int stored, final int delta){
		if(stored==0){return 0;}
		return Math.min(0xFF, Math.max(minNonEmpty, stored+delta));
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		if(Long.compareUnsigned(key, eeMask)>0){return;}
//		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);
		final int relNlz=nlz-minZeros;
//		if(relNlz<0){return;} // safety guard (eeMask should prevent this)

		// relNlzStored = relNlz+1, clamped to [1, maxRelNlzStored]
		final int newRelNlzStored=Math.min(relNlz+1, maxRelNlzStored);

		// Inverted mantissa: bits immediately after leading-1, inverted.
		// shift_offset=61 (precomputed); guard shift<0 for extreme NLZ (nlz>61).
		final int shift=shift_offset-nlz;
		final int invMantissa=(shift<0) ? 0 : (int)((~(key>>>shift))&mantissaMask);

		final int newStored=(newRelNlzStored<<mantissaBits)|invMantissa;
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
//		branch2++;
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		// Decrement minZeroCount when bucket leaves (empty ∪ tier-0) category.
		// Old bucket was in (empty ∪ tier-0) iff (oldStored >>> mantissaBits) < 2.
		// New bucket is tier-1+ iff newRelNlzStored >= 2.
		if((oldStored>>>mantissaBits)<2 && newRelNlzStored>=2){
			if(--minZeroCount<1){
				while(minZeroCount==0 && minZeros<wordlen){
					minZeros++;
					eeMask>>>=1;
					minZeroCount=countAndDecrement();
				}
			}
		}
	}

	/**
	 * Counts buckets that will become tier-0 after decrement (currently relNlzStored=2,
	 * i.e., stored in [2*minNonEmpty, 3*minNonEmpty-1]), then decrements all non-empty
	 * stored values by minNonEmpty (one relNlz tier).
	 * Only called when minZeroCount==0, guaranteeing all non-empty have relNlzStored>=2.
	 * Returns count of new tier-0 buckets (relNlzStored==1 after decrement).
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int i=0; i<buckets; i++){
			final int s=maxArray[i]&0xFF;
			if(s==0){continue;}
			final int decremented=s-minNonEmpty;
			maxArray[i]=(byte)decremented;
			if((decremented>>>mantissaBits)==1){newMinZeroCount++;}
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	/**
	 * Restores approximate hash magnitude from stored byte using the inverted mantissa,
	 * exactly mirroring DynamicDemiLog.restore() but for DDL8's byte encoding.
	 * Mean/GMean/MWA/Median all use this for sub-tier precision.
	 * HMean uses integer NLZ only (hllSumFilled); HMeanM uses fractional NLZ (hllSumFilledM).
	 */
	private long restoreDif(final int s){
		final int absNlz=(s>>>mantissaBits)-1+minZeros;
		if(absNlz==0){return Long.MAX_VALUE;}
		final long lowbits=(~s)&mantissaMask;        // uninvert: higher stored = smaller hash in tier
		final long mantissa=(1L<<mantissaBits)|lowbits;
		final int shift=wordlen-absNlz-mantissaBits-1;
		if(shift<0){return 1L;}
		return mantissa<<shift;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Returns cardinality estimates for calibration.
	 * Output order: 0=Mean(LL8), 1=HMean, 2=HMeanM(DDL8), 3=GMean, 4=HLL(HLL8),
	 *               5=LC, 6=Hybrid, 7=MWA, 8=MedianCorr
	 * <p>
	 * Mean/GMean/MWA/Median use integer absNlz (mantissa ignored) — LL8 mode.
	 * HLL uses integer absNlz across all buckets — HLL8 mode.
	 * HMeanM uses fractional absNlz from inverted mantissa — DDL8 mode.
	 */
	public double[] rawEstimates(){
		double difSum=0;
		double hllSumFilled=0;     // integer NLZ, filled buckets only (for HMean)
		double hllSumFilledM=0;    // fractional NLZ via mantissa, filled buckets (for HMeanM)
		double gSum=0;
		int count=0;
		sortBuf.clear();

		for(int i=0; i<buckets; i++){
			final int s=readBucket(i);
			if(s>0){
				// Mantissa-corrected dif for Mean/GMean/MWA/Median — mirrors DynamicDemiLog.restore().
				final long dif=restoreDif(s);
				difSum+=dif;
				hllSumFilled+=Math.pow(2.0, -(s>>>mantissaBits)+1-minZeros); // = 2^(-absNlz), integer NLZ for HMean
				// Fractional NLZ for DDL8/HMeanM: 2^(-absNlz + 0.5 - invMantissa/mantissaScale)
				// = 2^(-(s>>>mantissaBits) + 1.5 - (s&mantissaMask)/mantissaScale - minZeros)
				hllSumFilledM+=Math.pow(2.0, -(s>>>mantissaBits)+1.5-(s&mantissaMask)/mantissaScale-minZeros);
				gSum+=Math.log(Tools.max(1, dif));
				count++;
				sortBuf.add(dif);
			}
		}

		// All-buckets HLL sum: empty=1.0 (register=0), filled=2^(-absNlz)
		double hllSum=0;
		for(int i=0; i<buckets; i++){
			final int s=readBucket(i);
			if(s==0){
				hllSum+=1.0;
			}else{
				hllSum+=Math.pow(2.0, -(s>>>mantissaBits)+1-minZeros);
			}
		}

		final double alpha_m=0.7213/(1.0+1.079/buckets);
		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		final double gmean=Math.exp(gSum/div);
		sortBuf.sort();
		final long median=Tools.max(1, sortBuf.median());
		final double mwa=Tools.max(1.0, sortBuf.medianWeightedAverage());
		final int V=buckets-count;

		// HLL-style all-buckets estimate with LC fallback at low occupancy
		final double hmeanRaw=2*alpha_m*(double)buckets*(double)buckets/hllSum;
		double hmeanEst=hmeanRaw;
		if(hmeanEst<2.5*buckets && V>0){hmeanEst=(double)buckets*Math.log((double)buckets/V);}

		final double correction=(count+buckets)/(float)(buckets+buckets);
		final double hmeanPure =(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilled);
		final double hmeanPureM=(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilledM);

		final double meanEst   =2*(Long.MAX_VALUE/Tools.max(1.0, mean))*div*correction;
		final double gmeanEst  =2*(Long.MAX_VALUE/gmean)               *div*correction;
		final double mwaEst    =2*(Long.MAX_VALUE/mwa)                 *div*correction;
		final double medianCorr=2*(Long.MAX_VALUE/(double)median)      *div*correction;
		final double lcPure    =buckets*Math.log((double)buckets/Math.max(V, 0.5));

		final int trim=count/256;
		final int trimLow=Math.max(0, trim-V);
		double mean99Sum=0;
		final int mean99N=count-trimLow-trim;
		if(mean99N>0){for(int i=trim; i<count-trimLow; i++){mean99Sum+=sortBuf.get(i);}}
		final double mean99=(mean99N>0 ? mean99Sum/mean99N : mean);

		if(filledBuckets==0){return new double[10];}

		final double meanEstCF   =meanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEAN);
		final double hmeanPureMCF=hmeanPureM*CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEANM);

		final double hybridEst;
		final double hb0=0.5*buckets, hb1=1.5*buckets, hb2=3.0*buckets;
		if(lcPure<=hb0){
			hybridEst=lcPure;
		}else if(lcPure<=hb1){
			final double t=(lcPure-hb0)/(hb1-hb0);
			hybridEst=(1-t)*lcPure+t*meanEstCF;
		}else if(lcPure<=hb2){
			final double t=(lcPure-hb1)/(hb2-hb1);
			hybridEst=(1-t)*meanEstCF+t*hmeanPureMCF;
		}else{
			hybridEst=hmeanPureMCF;
		}

		final double mean99Est=2*(Long.MAX_VALUE/Tools.max(1.0, mean99))*div*correction;
		return new double[]{
			meanEstCF,
			hmeanPure  *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEAN),
			hmeanPureM *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEANM),
			gmeanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.GMEAN),
			hmeanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HLL),
			lcPure,
			hybridEst,
			mwaEst     *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MWA),
			medianCorr *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEDCORR),
			mean99Est  *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEAN99)
		};
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Byte register array: 1 byte per bucket, encoding ((relNlz+1)<<mantissaBits)|invMantissa. */
	private final byte[] maxArray;
	// mantissaBits, mantissaMask, minNonEmpty, maxRelNlzStored, mantissaScale: see static finals below

	private int minZeros=0;
	/** Count of (empty + tier-0) buckets; triggers minZeros floor advance when 0. */
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;
	private final LongList sortBuf=new LongList(buckets);
	// lastCardinality inherited from CardinalityTracker

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, added);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Hard-coded mantissa bits: 2-bit inverted mantissa (6-bit exponent = 64 tiers). */
	private static final int mantissaBits=2;
	private static final int mantissaMask=(1<<mantissaBits)-1;         // = 3
	/** Minimum stored byte value for a non-empty bucket (tier-0, invMantissa=0). */
	private static final int minNonEmpty=1<<mantissaBits;              // = 4
	/** Maximum relNlzStored before overflow clamp: 6-bit exponent = 63 tiers. */
	private static final int maxRelNlzStored=(1<<(8-mantissaBits))-1; // = 63
	/** Used in fractional NLZ formula: 1.5 - invMantissa/mantissaScale. */
	private static final double mantissaScale=1<<mantissaBits;        // = 4.0
	/** Precomputed shift offset: wordlen - mantissaBits - 1 = 61. */
	private static final int shift_offset=wordlen-mantissaBits-1;     // = 61

	public static double LC_CROSSOVER=0.75;
	public static double LC_SHARPNESS=20.0;
	public static boolean USE_LC=true;
	
	/** Default resource file for DDL correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDDL8.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the DDL correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

}
