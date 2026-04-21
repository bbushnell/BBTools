package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicDemiLog8v2: absolute-NLZ DemiLog with 16-bit char registers (6-bit NLZ + 10-bit mantissa).
 * <p>
 * Unlike DDL8 (relative storage), each bucket stores ABSOLUTE NLZ directly:
 *   score = (absNlz << mantissabits) | invMantissa, with score=0 meaning empty.
 * globalNLZ tracks the minimum absolute NLZ tier held by any bucket (dynamic early-exit threshold).
 * globalNLZ starts at 0 (not -1) because values are absolute.
 * <p>
 * Maintains two extra integers for metadata:
 * globalNLZ: The minimum absolute NLZ tier seen across all buckets.
 * minZeroCount: Number of buckets with NLZ == globalNLZ.
 * This drives a dynamic threshold that starts at zero and rises with cardinality, achieving
 * 99.9%+ early exit for large datasets without any loss of accuracy — skipped elements can
 * never improve a bucket that already holds a higher score.
 *
 * @author Brian Bushnell
 * @contributor Chloe
 * @date February 27, 2026
 */
public final class DynamicDemiLog8v2 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	DynamicDemiLog8v2(){
		this(2048, 31, -1, 0);
	}

	/** Construct from parsed command-line arguments. */
	DynamicDemiLog8v2(Parser p){
		super(p);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	/**
	 * Full constructor.
	 * @param buckets_ Number of buckets
	 * @param k_ Hash prefix length
	 * @param seed Random seed (-1 for default)
	 * @param minProb_ Minimum probability threshold
	 */
	DynamicDemiLog8v2(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	/** Create an independent copy with a fresh seed. */
	@Override
	public DynamicDemiLog8v2 copy(){return new DynamicDemiLog8v2(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Scans the bucket array once, accumulating all sums needed for estimation.
	 * This is the only per-subclass method required — all estimator logic
	 * lives in the returned CardinalityStats object.
	 */
	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		// DDL8v2 stores absolute NLZ directly: score = (absNlz << mantissabits) | invMantissa.
		// Pack into char[] as ((absNlz+1) << mantissabits) | invMant (+1 so val!=0 for filled).
		final char[] packedBuckets=new char[maxArray.length];
		int filledCount=0;
		for(int i=0; i<maxArray.length; i++){
			final int max=maxArray[i];
			if(max>0){
				final int absNlz=max>>mantissabits;
				if(absNlz<64){nlzCounts[absNlz+1]++;}
				filledCount++;
				packedBuckets[i]=(char)(((absNlz+1)<<mantissabits)|(max&mask));
			}
		}
		nlzCounts[0]=buckets-filledCount;
		// nlzBits=6, mantissaBits=2, mantissaOffset=0.5 (matches DDL8)
		return new CardStats(packedBuckets, nlzCounts, 6, 0, 0, mantissabits,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0.5);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		final double rawHyb=s.hybridDDL(); // CF already inside blend
		long card=(long)(rawHyb);
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinalityStatic=lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicDemiLog8v2)log);
	}

	/** Merges another DynamicDemiLog8v2 into this one by taking per-bucket max, then rescanning. */
	public void add(DynamicDemiLog8v2 log){
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
			scanFrom(Math.max(globalNLZ, log.globalNLZ));
			//Recompute filledBuckets after merge; incremental tracking is not possible here
			filledBuckets=0;
			for(char c : maxArray){if(c>0){filledBuckets++;}}
		}
	}

	/** Compares this instance's bucket array to another's bucket by bucket.
	 * @return int[] {lower, equal, higher} bucket counts */
	public int[] compareTo(DynamicDemiLog8v2 blog){return compare(maxArray, blog.maxArray);}

	@Override
	public void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		//Earliest possible exit, fastest
		if(Long.compareUnsigned(key, eeMask)>0){return;}
//		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);

		// Global early exit if not using eeMask
//		if(nlz<globalNLZ){return;}

		//Precalculate everything necessary
		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		final int score=Math.max(1, (nlz<<mantissabits)+(int)((~(key>>>shift))&mask));//FP16; min 1 so score=0 means empty
		final int oldValue=maxArray[bucket];//Required memory read

		if(LAZY_ALLOCATE){//Optional MicroIndex for low cardinality (causes speed decrease)
			final long micro=(key>>bucketBits)&0x3FL;
			microIndex|=(1L<<micro);
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}//Allows lazy array allocation
		}

		//Optional early exit reduces writes, countArray access, and branches.
		//Expected to be usually taken, particularly when buckets is large.
		if(score<oldValue){return;}
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
		if(nlz>nlzOld && nlzOld==globalNLZ && --minZeroCount<1){//Promotion from bottom tier
			/* NOTE - Due to a major Eclipse JDK 24 -> Java 8 target bytecode generation bug,
			 * enabling both of these assertions causes a 40% slowdown even with -da.
			 * Do not enable them in production.
			 */
//			assert(minZeroCount>=0 && globalNLZ<=64) : minZeroCount+", "+globalNLZ;
			while(minZeroCount==0 && globalNLZ<wordlen){//Scan for new tier
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countTermsInTier(globalNLZ, maxArray);
			}
//			assert(minZeroCount>0 && minZeroCount<=buckets) : minZeroCount+", "+globalNLZ;
		}
	}

	/** Rescans to find the correct globalNLZ starting from the given absolute NLZ value. */
	private int scanFrom(int nlz){
		globalNLZ=nlz-1;
		minZeroCount=0;
		eeMask=-1L;
		//Technically this could be 2-pass: find the lowest, THEN count that tier
		//But in practice that would be slower
		while(minZeroCount==0 && globalNLZ<wordlen){
			globalNLZ++;
			eeMask>>>=1;
			minZeroCount=countTermsInTier(globalNLZ, maxArray);
		}
		return globalNLZ;
	}

	/** Not used; CF correction handled via CF_MATRIX. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/** Returns the raw 16-bit countArray (count of observations per bucket). */
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
		final CardStats s=summarize();
		return AbstractCardStats.buildLegacyArray(s, s.hybridDDL());
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

	/** Counts buckets whose absolute NLZ tier equals nlz, using branchless bit arithmetic. */
	private static int countTermsInTier(final int nlz, final char[] array){
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

	/** Placeholder for ANI calculation; not yet implemented. */
	public static float ani(int lower, int equal, int higher){return -1;}
	/** Placeholder for completeness calculation; not yet implemented. */
	public static float completeness(int lower, int equal, int higher){return -1;}
	/** Placeholder for contamination calculation; not yet implemented. */
	public static float contam(int lower, int equal, int higher){return -1;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compressed 16-bit bucket maxima; 6-bit absolute NLZ + 10-bit inverted mantissa (score=0 means empty). */
	private final char[] maxArray;
	/** Count of observations at the current maximum score for each bucket. */
	private final char[] countArray;
	/** Minimum absolute NLZ tier held by any bucket; dynamic early-exit threshold.
	 * Starts at 0 (absolute values, not relative). Rises monotonically toward 63. */
	private int globalNLZ=0;
	/** Number of buckets with NLZ == globalNLZ. */
	private int minZeroCount;
	/** Number of buckets with any data (maxArray[i] > 0). Maintained incrementally
	 * in hashAndStore() for O(1) filledBuckets() and occupancy() calls. */
	private int filledBuckets=0;
	/** Early-exit mask: filters hashes whose NLZ is below globalNLZ. */
	private long eeMask=-1L;
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	/** For tracking branch prediction rates; disable in production. */
	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Number of mantissa bits for floating-point compression; 10 is maximum. */
	private static final int mantissabits=2;
	private static final int mask=(1<<mantissabits)-1;
	private static final int offset=wordlen-mantissabits-1;

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
