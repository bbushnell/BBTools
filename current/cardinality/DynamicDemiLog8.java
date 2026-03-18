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

	/**
	 * Scans the bucket array once, accumulating all sums needed for estimation.
	 * This is the only per-subclass method required — all estimator logic
	 * lives in the returned CardinalityStats object.
	 */
	private CardinalityStats summarize(){
		double difSum=0;
		double hllSumFilled=0;
		double hllSumFilledM=0;
		double gSum=0;
		int count=0;
		if(USE_SORTBUF){
			if(sortBuf==null){sortBuf=new LongList(buckets);}
			sortBuf.clear();
		}
		if(nlzCounts==null){nlzCounts=new int[64];}
		else{java.util.Arrays.fill(nlzCounts, 0);}

		for(int i=0; i<buckets; i++){
			final int s=readBucket(i);
			if(s>0){
				final int absNlz=(s>>>mantissaBits)-1+minZeros;
				if(absNlz<64){nlzCounts[absNlz]++;}
				final long dif=restoreDif(s);
				difSum+=dif;
				hllSumFilled +=Math.pow(2.0, -absNlz);
				hllSumFilledM+=Math.pow(2.0, -(s>>>mantissaBits)+1.5-(s&mantissaMask)/mantissaScale-minZeros);
				gSum+=Math.log(Tools.max(1, dif));
				count++;
				if(USE_SORTBUF){sortBuf.add(dif);}
			}
		}
		return new CardinalityStats(difSum, hllSumFilled, hllSumFilledM,
		                            gSum, count, buckets, sortBuf, CF_MATRIX, CF_BUCKETS,
		                            CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys, 0, nlzCounts, minZeros);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardinalityStats s=summarize();
		final long card=Math.min(clampToAdded ? added : Long.MAX_VALUE, (long)s.hybridDDL());
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
		
		if(USE_MICRO){//Optional MicroIndex for low cardinality
			final long micro=(key>>bucketBits)&0x3FL;
			microIndex|=(1L<<micro);
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS) {return;}//Allows lazy array allocation
		}
		
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
	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(s.hybridDDL());
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
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

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
