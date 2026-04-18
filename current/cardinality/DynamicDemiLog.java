package cardinality;

import parse.Parser;
import shared.Tools;
import structures.ByteBuilder;
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

	/** Public factory for creating a new empty DDL.
	 * @param buckets Number of buckets (power of 2)
	 * @param k K-mer length
	 * @param seed Hash seed
	 * @param minProb Minimum probability filter */
	public static DynamicDemiLog create(int buckets, int k, long seed, float minProb){
		return new DynamicDemiLog(buckets, k, seed, minProb);
	}

	/** Creates a DDL from a pre-filled maxArray (e.g., loaded from file).
	 * Reconstructs filledBuckets, minZeros, minZeroCount, and eeMask.
	 * @param loaded Pre-filled bucket values (will be copied)
	 * @param id_ Taxonomy or file ID
	 * @param name_ Optional name
	 * @param k_ K-mer length */
	public static DynamicDemiLog fromArray(char[] loaded, int id_, String name_, int k_){
		DynamicDemiLog ddl=new DynamicDemiLog(loaded.length, k_, defaultSeed, 0);
		System.arraycopy(loaded, 0, ddl.maxArray, 0, loaded.length);
		ddl.id=id_;
		ddl.name=name_;
		ddl.filledBuckets=0;
		for(char c : ddl.maxArray){if(c>0){ddl.filledBuckets++;}}
		ddl.scanFrom(0);
		return ddl;
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
	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		// DDL stores (absNlz << 10) | invMantissa directly — no tier promotion.
		// Pack into char[] as ((absNlz+1) << 10) | invMant (+1 so val!=0 for filled).
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
		// DDL: nlzBits=6, mantissaBits=10, mantissaOffset=0.0 (no +0.5 bias)
		return new CardStats(packedBuckets, nlzCounts, 6, 0, 0, mantissabits,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0.0);
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
		final int score=Math.max(1, (nlz<<mantissabits)+(int)((~(key>>>shift))&mask));//FP16; min 1 so score=0 means empty
		final int oldValue=maxArray[bucket];//Required memory read
		
		if(LAZY_ALLOCATE){//Optional MicroIndex for low cardinality (causes speed decrease)
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

	/** Returns the underlying maxArray for external comparison. */
	public char[] maxArray(){return maxArray;}

	/** Returns the number of buckets with any data (maxArray[i] > 0). O(1). */
	public int filledBuckets(){return filledBuckets;}

	/** Returns the fraction of buckets with any data (filled / total). O(1). */
	public double occupancy(){return (double)filledBuckets/buckets;}

	/** Appends this DDL as tab-separated A48-encoded bucket values.
	 * Each 16-bit value becomes 1-3 A48 characters. */
	public ByteBuilder toBytes(ByteBuilder bb){
		for(int i=0; i<maxArray.length; i++){
			if(i>0){bb.tab();}
			bb.appendA48((int)maxArray[i]);
		}
		return bb;
	}

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

	/** Compares two DDLs excluding empty-empty bucket pairs.
	 * @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailed(char[] arrayA, char[] arrayB){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<arrayA.length; i++){
			final int a=arrayA[i], b=arrayB[i];
			if(a==0 && b==0){bothEmpty++; continue;}
			int dif=a-b;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/** Instance version of compareDetailed. */
	public int[] compareToDetailed(DynamicDemiLog other){
		return compareDetailed(maxArray, other.maxArray);
	}

	/** Containment of the smaller set in the larger.
	 * equal / min(equal+lower, equal+higher), analogous to
	 * BBSketch's hits / min(queryDivisor, refDivisor).
	 * @param lower Buckets where A < B (from compareDetailed, excluding empty-empty)
	 * @param equal Matching buckets (excluding empty-empty)
	 * @param higher Buckets where A > B (excluding empty-empty) */
	public static float containment(int lower, int equal, int higher){
		int minDiv=Math.min(equal+lower, equal+higher);
		if(minDiv<=0 || equal<=0){return 0;}
		return Math.min(1f, (float)equal/minDiv);
	}

	/** Containment of A in B: fraction of B's buckets matched by A.
	 * If A has everything B has, this is 1.0 regardless of A's extra content. */
	public static float containmentAB(int lower, int equal, int higher){
		int denom=equal+lower;
		return denom>0 ? Math.min(1f, (float)equal/denom) : 0;
	}

	/** Containment of B in A: fraction of A's buckets matched by B. */
	public static float containmentBA(int lower, int equal, int higher){
		int denom=equal+higher;
		return denom>0 ? Math.min(1f, (float)equal/denom) : 0;
	}

	/** Estimates ANI from DDL bucket comparison.
	 * Uses max containment (containment of the smaller set in the larger),
	 * then ANI = containment^(1/k), as in BBSketch.
	 * @param lower Buckets where A < B (excluding empty-empty)
	 * @param equal Matching buckets (excluding empty-empty)
	 * @param higher Buckets where A > B (excluding empty-empty)
	 * @param k K-mer length used for hashing */
	public static float ani(int lower, int equal, int higher, int k){
		float c=containment(lower, equal, higher);
		return c>0 ? (float)Math.exp(Math.log(c)/k) : 0;
	}

	/** @deprecated Use ani(lower, equal, higher, k) instead */
	@Deprecated
	public static float ani(int lower, int equal, int higher){
		return ani(lower, equal, higher, 31);
	}

	/** Genomic completeness of A relative to B: does A's genome cover B's?
	 * If A is "higher" in a bucket, A had a k-mer B didn't have.
	 * The ratio (equal+higher)/(equal+lower) estimates |A|/|B|, clamped to [0,1].
	 * @return Estimated fraction of B's genome covered by A */
	public static float completeness(int lower, int equal, int higher){
		int denom=equal+lower;
		if(denom<=0){return 0;}
		return Math.min(1f, (float)(equal+higher)/denom);
	}

	/** Genomic completeness of B relative to A. */
	public static float completenessBA(int lower, int equal, int higher){
		return completeness(higher, equal, lower);
	}

	/** Contamination: requires multi-reference comparison; not supported pairwise.
	 * @return -1 (not applicable) */
	public static float contam(int lower, int equal, int higher){
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
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	/*--------------------------------------------------------------*/

	/** Taxonomy or file ID for this DDL; -1 if unassigned. */
	public int id=-1;
	/** Optional name (e.g., organism name). */
	public String name;

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
