package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicDemiLog cardinality estimator extending DemiLogLog (LogLog16) with an adaptive
 * early-exit mechanism. Maintains two extra integers for metadata:
 * minZeros: The minimum number of leading zeros seen across all buckets
 * minZeroCount: Number of buckets with minZeros
 * This drives a dynamic threshold that starts at zero and rises with cardinality, achieving
 * 99.9%+ early exit for large datasets without any loss of accuracy—skipped elements can
 * never improve a bucket that already holds a higher score.
 *
 * @author Brian Bushnell
 * @contributor Chloe
 * @date February 27, 2026
 */
public final class DynamicDemiLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a DynamicLogLog with default parameters.
	 * Uses 2048 buckets, k=31, random seed, and no minimum probability filtering. */
	DynamicDemiLog(){
		this(2048, 31, -1, 0);
	}

	/** Creates a DynamicLogLog with parameters parsed from command-line arguments.
	 * @param p Parser containing configuration from command-line flags */
	DynamicDemiLog(Parser p){
		super(p);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	/**
	 * Creates a DynamicLogLog with specified parameters.
	 * @param buckets_ Number of buckets (counters) for the hash table
	 * @param k_ K-mer length for sequence hashing
	 * @param seed Random number generator seed; -1 for random seed
	 * @param minProb_ Ignore k-mers with under this probability of being correct
	 */
	DynamicDemiLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	@Override
	public DynamicDemiLog copy(){return new DynamicDemiLog(buckets, k, -1, minProb);}

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
		sortBuf.clear();

		for(int i=0; i<maxArray.length; i++){
			final int max=maxArray[i];
			final long dif=restore(max);
			if(max>0 && dif>0){
				difSum+=dif;
				hllSumFilled +=Math.pow(2.0, -(max>>mantissabits));
				hllSumFilledM+=Math.pow(2.0, -(max>>mantissabits)+0.0-(max&mask)/1024.0);
				gSum+=Math.log(Tools.max(1, dif));
				count++;
				sortBuf.add(dif);
			}
		}
		return new CardinalityStats(difSum, hllSumFilled, hllSumFilledM,
		                            gSum, count, buckets, sortBuf, CF_MATRIX, CF_BUCKETS,
		                            CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys, microIndex);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardinalityStats s=summarize();
		long card=(long)s.hybridDDL();
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinalityStatic=lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicDemiLog)log);
	}

	public void add(DynamicDemiLog log){
		added+=log.added;
		branch1+=log.branch1;
		branch2+=log.branch2;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				final char maxA=maxArray[i], maxB=log.maxArray[i];
				final char countA=countArray[i], countB=log.countArray[i];
				maxArray[i]=Tools.max(maxA, maxB);
				if(maxA==maxB){countArray[i]=(char)Tools.min(countA+(int)countB, Character.MAX_VALUE);}
				else{countArray[i]=(maxA>maxB ? countA : countB);}
			}
			scanFrom(Math.max(minZeros, log.minZeros));
			//Recompute filledBuckets after merge; incremental tracking is not possible here
			filledBuckets=0;
			for(char c : maxArray){if(c>0){filledBuckets++;}}
		}
	}

	/** Compares two DynamicLogLog instances bucket by bucket.
	 * @return int[] {lower, equal, higher} bucket counts */
	public int[] compareTo(DynamicDemiLog blog){return compare(maxArray, blog.maxArray);}

	@Override
	public void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		//Earliest possible exit, fastest
		if(Long.compareUnsigned(key, eeMask)>0) {return;}
//		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);

		// Global early exit if not using eeMask
//		if(nlz<minZeros){return;}

		//Precalculate everything necessary
		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		final int score=(nlz<<mantissabits)+(int)((~(key>>>shift))&mask);//FP16 representation
		final int oldValue=maxArray[bucket];//Required memory read
		
		if(USE_MICRO){//Optional MicroIndex for low cardinality (causes speed decrease)
			final long micro=(key>>bucketBits)&0x3FL;
			microIndex|=(1L<<micro);
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS) {return;}//Allows lazy array allocation
		}
		
		//Optional early exit reduces writes, countArray access, and branches.
		//Expected to be usually taken, particularly when buckets is large.
		if(score<oldValue) {return;}
//		branch2++;
		lastCardinality=-1;
		final int newValue=Math.max(score, oldValue);
		final int nlzOld=(oldValue>>mantissabits);

//		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
		//Update bucket; required write to cached line only if score>oldValue
		maxArray[bucket]=(char)newValue;
		//Track filled bucket count: increment when bucket transitions from empty to non-empty
		if(oldValue==0 && newValue>0){filledBuckets++;}

		//Update count - totally optional, only for histograms
		//Can be debranched with clever math
		final char count=countArray[bucket];
		countArray[bucket]=(char)(oldValue>score ? count :
			oldValue==score ? Math.max(count, (char)(count+1)) : 1);

		//Update the dynamic early exit threshold
		if(nlz>nlzOld && nlzOld==minZeros && --minZeroCount<1) {//Promotion from bottom tier 
			/* NOTE - Due to a major Eclipse JDK 24 -> Java 8 target bytecode generation bug,
			 * enabling both of these assertions causes a 40% slowdown even with -da.
			 * Do not enable them in production.
			 */
			//			assert(minZeroCount>=0 && minZeros<=64) : minZeroCount+", "+minZeros;
			while(minZeroCount==0 && minZeros<wordlen) {//Scan for new tier
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countTermsInTier(minZeros, maxArray);
			}
			//			assert(minZeroCount>0 && minZeroCount<=buckets) : minZeroCount+", "+minZeros;
		}

	}

	private int scanFrom(int nlz) {
		minZeros=nlz-1;
		minZeroCount=0;
		eeMask=-1L;
		//Technically this could be 2-pass: find the lowest, THEN count that tier
		//But in practice that would be slower
		while(minZeroCount==0 && minZeros<wordlen) {
			minZeros++;
			eeMask>>>=1;
			minZeroCount=countTermsInTier(minZeros, maxArray);
		}
		return minZeros;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public char[] counts16(){return countArray;}

	/** Returns the number of buckets with any data (maxArray[i] > 0). O(1). */
	public int filledBuckets(){return filledBuckets;}

	/** Returns the fraction of buckets with any data (filled / total). O(1). */
	public double occupancy(){return (double)filledBuckets/buckets;}

	/**
	 * Returns cardinality estimates for calibration, without the added-count cap.
	 * Output order: 0=Mean, 1=HMean, 2=HMeanM, 3=GMean, 4=HLL, 5=LC, 6=Hybrid, 7=MWA, 8=MedianCorr, 9=Mean99
	 * Hybrid (index 6) is pre-corrected: CF already applied to its Mean and HMeanM components.
	 * LC (index 5) and Hybrid (index 6) receive no additional CF.
	 */
	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(Math.max(s.hybridDDL(), s.microCardinality()));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Restores floating-point compressed score to approximate original hash magnitude.
	 * @param score Compressed 16-bit value from maxArray
	 * @return Approximate original hash value before compression */
	private static long restore(int score){
		int leading=(int)(score>>>mantissabits);
		// nlz=0 case: true hash magnitude is in [2^63, 2^64-1] (unsigned), which
		// overflows signed long (shift of 53 bits causes Long.MIN_VALUE). Approximate
		// as Long.MAX_VALUE — this keeps the bucket contributing to value-based estimates
		// (giving meanEst≈1 for a single nlz=0 bucket) rather than being silently excluded.
		if(leading==0){return Long.MAX_VALUE;}
		long lowbits=(~score)&mask;
		long mantissa=(1L<<mantissabits)|lowbits;
		int shift=wordlen-leading-mantissabits-1;
		long original=mantissa<<shift;
		return original;
	}

	private static int countTermsInTier(final int nlz, final char[] array) {
		int sum=0;
		for(final char c : array){
			final int diff=(c>>mantissabits)^nlz;  //0 iff tier matches
			final int equalsBit=((~(diff|-diff))>>>31);  //1 if diff==0, else 0
			sum+=equalsBit;
		}
		return sum;
	}

	/** Compares two maxArrays bucket by bucket using branchless bit arithmetic.
	 * @return int[]{lower, equal, higher} bucket counts */
	public static int[] compare(char[] arrayA, char[] arrayB){
		int lower=0, equal=0, higher=0;
		for(int i=0; i<arrayA.length; i++){
			final int a=arrayA[i], b=arrayB[i];
			int dif=a-b;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher};
	}

	public static float ani(int lower, int equal, int higher) {
		//TODO
		return -1;
	}

	public static float completeness(int lower, int equal, int higher) {
		//TODO
		return -1;
	}

	public static float contam(int lower, int equal, int higher) {
		//TODO
		return -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compressed 16-bit bucket maxima; 6-bit leading-zero count + 10-bit mantissa. */
	private final char[] maxArray;
	/** Count of observations at the current maximum score for each bucket. */
	private final char[] countArray;
	/** Minimum leading-zero tier held by any bucket; dynamic early-exit threshold.
	 * Rises monotonically from 0 toward 63 as cardinality grows. */
	private int minZeros=0;
	/** Number of buckets with minZeros */
	private int minZeroCount;
	/** Number of buckets with any data (maxArray[i] > 0). Maintained incrementally
	 * in hashAndStore() for O(1) filledBuckets() and occupancy() calls. */
	private int filledBuckets=0;
	private long eeMask=-1L;
	/** Reusable sort buffer for rawEstimates(); avoids per-call allocation. */
	private final LongList sortBuf=new LongList(buckets);
	// lastCardinality inherited from CardinalityTracker

	/*--------------------------------------------------------------*/

	/** For tracking branch prediction rates; disable in production */
	public long branch1=0, branch2=0;
	public double branch1Rate() {return branch1/(double)Math.max(1, added);}
	public double branch2Rate() {return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Number of mantissa bits for floating-point compression; 10 is maximum. */
	private static final int mantissabits=10;
	private static final int mask=(1<<mantissabits)-1;
	private static final int offset=wordlen-mantissabits-1;

	//	/** Occupancy fraction (filled buckets / total buckets) at which LinearCounting and
	//	 * value-based estimates contribute equally in the sigmoid blend.  Below this point
	//	 * the blend leans toward LinearCounting; above it toward the value-based estimate. */

	/** Default resource file for DDL correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDDL.tsv.gz";
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
