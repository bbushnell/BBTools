package cardinality;

import parse.Parser;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * DynamicDemiLog cardinality estimator extending DemiLogLog (LogLog16) with an adaptive
 * early-exit mechanism. Each bucket stores a 16-bit absolute value:
 * (absNlz << 10) | invMantissa, where absNlz = number of leading zeros of the hash.
 *
 * globalNLZ tracks the minimum NLZ tier (floor) across all buckets, driving eeMask
 * for early exit. File I/O uses relative encoding via toBytesRelative()/fromArray(offset)
 * with an #offset header, but internal storage is always absolute.
 *
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

	/** Creates a DDL from a pre-filled array of absolute-encoded values (legacy format).
	 * Converts to relative encoding internally.
	 * @param loaded Pre-filled absolute bucket values
	 * @param id_ Taxonomy or file ID
	 * @param name_ Optional name
	 * @param k_ K-mer length */
	public static DynamicDemiLog fromArray(char[] loaded, int id_, String name_, int k_){
		return fromArray(loaded, id_, name_, k_, -1);
	}

	/** Creates a DDL from a pre-filled array.
	 * @param loaded Pre-filled bucket values
	 * @param id_ Taxonomy or file ID
	 * @param name_ Optional name
	 * @param k_ K-mer length
	 * @param offset If >=0, loaded values are relative to this globalNLZ.
	 *               If -1, loaded values are absolute (legacy format) and will be converted. */
	public static DynamicDemiLog fromArray(char[] loaded, int id_, String name_, int k_, int offset){
		DynamicDemiLog ddl=new DynamicDemiLog(loaded.length, k_, defaultSeed, 0);
		ddl.id=id_;
		ddl.name=name_;

		if(offset>=0){
			// New format: values are relative to offset. Convert to absolute.
			for(int i=0; i<loaded.length; i++){
				if(loaded[i]>0){
					final int relNlz=loaded[i]>>mantissabits;
					final int absNlz=relNlz+offset;
					ddl.maxArray[i]=(char)((absNlz<<mantissabits)|(loaded[i]&mask));
				}
			}
		}else{
			// Legacy format: values are already absolute.
			System.arraycopy(loaded, 0, ddl.maxArray, 0, loaded.length);
		}

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

	/** Merges another DDL into this one using absolute value max. */
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
			// Rescan for floor tier — matches old scanFrom(max(minZeros, log.minZeros))
			final int scanStart=Math.max(Math.max(0, globalNLZ+1), Math.max(0, log.globalNLZ+1));
			globalNLZ=scanStart-1;
			minZeroCount=0;
			eeMask=-1L;
			while(minZeroCount==0 && globalNLZ<wordlen){
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countTermsInTier(globalNLZ, maxArray);
			}
			filledBuckets=0;
			for(char c : maxArray){if(c>0){filledBuckets++;}}
		}
	}

	/** Compares two DynamicLogLog instances bucket by bucket.
	 * @return int[] {lower, equal, higher} bucket counts */
	public int[] compareTo(DynamicDemiLog blog){return compare(maxArray, blog.maxArray);}

	/** Instance version of compareDetailed.
	 * @return int[] {lower, equal, higher, bothEmpty} */
	public int[] compareToDetailed(DynamicDemiLog other){
		return compareDetailed(maxArray, other.maxArray);
	}

	@Override
	public void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		//Earliest possible exit, fastest
		if(Long.compareUnsigned(key, eeMask)>0) {return;}
//		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);

		//Precalculate everything necessary
		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		final int score=Math.max(1, (nlz<<mantissabits)+(int)((~(key>>>shift))&mask));
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
		if(nlz>nlzOld && nlzOld==globalNLZ && --minZeroCount<1){
			while(minZeroCount==0 && globalNLZ<wordlen){
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}

	}

	/** Counts buckets at the current tier (used by advancement loop). */
	private int countAndDecrement(){
		return countTermsInTier(globalNLZ, maxArray);
	}

	/** Scans from the given starting tier to find the first occupied tier and set globalNLZ. */
	private int scanFrom(int nlz) {
		globalNLZ=nlz-1;
		minZeroCount=0;
		eeMask=-1L;
		while(minZeroCount==0 && globalNLZ<wordlen) {
			globalNLZ++;
			eeMask>>>=1;
			minZeroCount=countTermsInTier(globalNLZ, maxArray);
		}
		return globalNLZ;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public char[] counts16(){return countArray;}

	/** Returns the underlying maxArray for internal comparison. */
	public char[] maxArray(){return maxArray;}

	/** Returns a copy of maxArray (absolute values). For DDLIndex and external consumers. */
	public char[] toAbsoluteArray(){
		return java.util.Arrays.copyOf(maxArray, maxArray.length);
	}

	/** Returns the number of buckets with any data (maxArray[i] > 0). O(1). */
	public int filledBuckets(){return filledBuckets;}

	/** Returns the fraction of buckets with any data (filled / total). O(1). */
	public double occupancy(){return (double)filledBuckets/buckets;}

	/** Returns the global NLZ floor (minimum absolute tier across all buckets). */
	public int getGlobalNLZ(){return globalNLZ;}

	/** Appends this DDL as tab-separated A48-encoded absolute bucket values. */
	public ByteBuilder toBytes(ByteBuilder bb){
		for(int i=0; i<maxArray.length; i++){
			if(i>0){bb.tab();}
			bb.appendA48((int)maxArray[i]);
		}
		return bb;
	}

	/** Appends relative-encoded bucket values as tab-separated A48.
	 * Converts from absolute internal storage to relative using globalNLZ as offset.
	 * Caller is responsible for writing the #offset header first. */
	public ByteBuilder toBytesRelative(ByteBuilder bb){
		for(int i=0; i<maxArray.length; i++){
			if(i>0){bb.tab();}
			final int abs=maxArray[i];
			if(abs==0){bb.appendA48(0);}
			else{
				final int absNlz=abs>>mantissabits;
				final int relNlz=absNlz-globalNLZ;
				bb.appendA48((relNlz<<mantissabits)|(abs&mask));
			}
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
	 * @param stored Compressed 16-bit value from maxArray (relative encoding)
	 * @return Approximate original hash value before compression */
	private long restore(int stored){
		final int absNlz=(stored>>mantissabits)+globalNLZ;
		if(absNlz<=0){return Long.MAX_VALUE;}
		long lowbits=(~stored)&mask;
		long mantissa=(1L<<mantissabits)|lowbits;
		int shift=wordlen-absNlz-mantissabits-1;
		if(shift<0){return 1L;}
		return mantissa<<shift;
	}

	/** Counts buckets whose absolute NLZ matches the given tier.
	 * Used during add() to find the new floor tier from absolute-encoded merged values. */
	private static int countTermsInTier(final int nlz, final char[] array) {
		int sum=0;
		for(final char c : array){
			final int diff=(c>>mantissabits)^nlz;  //0 iff tier matches
			final int equalsBit=((~(diff|-diff))>>>31);  //1 if diff==0, else 0
			sum+=equalsBit;
		}
		return sum;
	}

	/** Converts a relative stored value to its absolute equivalent. */
	private static int toAbsolute(int stored, int globalNLZ){
		if(stored==0){return 0;}
		return ((stored>>mantissabits)+globalNLZ)<<mantissabits | (stored&mask);
	}

	/** Compares two DDLs converting to absolute values.
	 * @return int[]{lower, equal, higher} bucket counts */
	private static int[] compareAbsolute(DynamicDemiLog a, DynamicDemiLog b){
		int lower=0, equal=0, higher=0;
		for(int i=0; i<a.maxArray.length; i++){
			final int va=toAbsolute(a.maxArray[i], a.globalNLZ);
			final int vb=toAbsolute(b.maxArray[i], b.globalNLZ);
			int dif=va-vb;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			lower+=nbit;
			higher+=hbit;
			equal+=1-nbit-hbit;
		}
		return new int[]{lower, equal, higher};
	}

	/** Compares two DDLs with absolute conversion, excluding empty-empty pairs.
	 * @return int[]{lower, equal, higher, bothEmpty} */
	private static int[] compareDetailedAbsolute(DynamicDemiLog a, DynamicDemiLog b){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<a.maxArray.length; i++){
			final int va=toAbsolute(a.maxArray[i], a.globalNLZ);
			final int vb=toAbsolute(b.maxArray[i], b.globalNLZ);
			if(va==0 && vb==0){bothEmpty++; continue;}
			int dif=va-vb;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			lower+=nbit;
			higher+=hbit;
			equal+=1-nbit-hbit;
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/** Compares two absolute-encoded maxArrays bucket by bucket using branchless bit arithmetic.
	 * Only valid when both arrays use the same absolute encoding (e.g., from toAbsoluteArray()).
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

	/** Compares two absolute-encoded arrays excluding empty-empty bucket pairs.
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

	/** Minimum reliable divisor for WKID calculation.
	 * When fewer than this many buckets contribute to the smaller genome's
	 * side of the comparison, results are unreliable. */
	private static final int MIN_DIVISOR=6;

	/** WKID: equal / (equal + min(lower, higher)).
	 * min(lower, higher) selects the mismatch count from the smaller genome —
	 * the genome whose "exceeds" count represents real mismatches rather than
	 * expected size differences. Analogous to BBSketch's hits / minDivisor.
	 * When the divisor < MIN_DIVISOR, pads to MIN_DIVISOR to prevent
	 * inflated scores from random bucket collisions.
	 * @param lower Buckets where A < B (from compareDetailed, excluding empty-empty)
	 * @param equal Matching buckets (excluding empty-empty)
	 * @param higher Buckets where A > B (excluding empty-empty) */
	public static float wkid(int lower, int equal, int higher){
		int div=Math.min(equal+lower, equal+higher);
		if(div<MIN_DIVISOR){div=MIN_DIVISOR;}
		if(equal<=0){return 0;}
		return Math.min(1f, (float)equal/div);
	}

	/** Containment of A in B: equal / (equal + higher).
	 * From A's perspective, higher = buckets where A exceeds B = A's mismatches.
	 * If A is smaller than B, this estimates what fraction of A's content
	 * is shared with B. */
	public static float containmentAB(int lower, int equal, int higher){
		int denom=equal+higher;
		return denom>0 ? Math.min(1f, (float)equal/denom) : 0;
	}

	/** Containment of B in A: equal / (equal + lower). */
	public static float containmentBA(int lower, int equal, int higher){
		int denom=equal+lower;
		return denom>0 ? Math.min(1f, (float)equal/denom) : 0;
	}

	/** Estimates ANI from DDL bucket comparison.
	 * WKID = equal / (equal + min(lower, higher)), then ANI = WKID^(1/k).
	 * min(lower, higher) selects the mismatch count from the smaller genome,
	 * analogous to BBSketch's minDivisor approach.
	 * @param lower Buckets where A < B (excluding empty-empty)
	 * @param equal Matching buckets (excluding empty-empty)
	 * @param higher Buckets where A > B (excluding empty-empty)
	 * @param k K-mer length used for hashing */
	public static float ani(int lower, int equal, int higher, int k){
		float w=wkid(lower, equal, higher);
		return w>0 ? (float)Math.exp(Math.log(w)/k) : 0;
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

	/** Compressed 16-bit bucket maxima; (relNlzStored << 10) | invMantissa
	 * where relNlzStored = nlz - globalNLZ. Stored value of 0 means empty. */
	private final char[] maxArray;
	/** Count of observations at the current maximum score for each bucket. */
	private final char[] countArray;
	/** Floor NLZ tier, matching old minZeros behavior. Starts at 0.
	 * Rises monotonically as cardinality increases and buckets advance past lower tiers. */
	private int globalNLZ=0;
	/** Number of buckets below floor level (stored < minNonEmpty), including empty;
	 * triggers globalNLZ floor advance when 0. */
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
	/** Minimum stored value for a non-empty bucket above floor (relNlzStored >= 1). */
	private static final int minNonEmpty=1<<mantissabits; // 1024
	/** Maximum relNlzStored before overflow clamp: 6-bit exponent = 63 tiers. */
	private static final int maxRelNlzStored=(1<<(16-mantissabits))-1; // 63

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
