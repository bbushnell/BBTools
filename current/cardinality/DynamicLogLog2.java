package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicLogLog2: 2-bit packed variant of DynamicLogLog using relative NLZ storage.
 * <p>
 * Packs 16 buckets into each int, 2 bits per bucket.
 * Encoding: 0 = empty; 1-3 = (relNlz + 1), where relNlz = absoluteNlz - globalNLZ - 1.
 * No mantissa — coarsest precision per bucket, but allows 2x bucket count vs DLL4
 * for the same memory footprint.
 * Overflow (relNlz >= 3) clamped to stored=3 and absorbed by CF matrix.
 * <p>
 * Memory: 2 bits/bucket × 2048 buckets = 4096 bits = 512 bytes.
 * 16 buckets per int; pure power-of-2 packing (no modulo, no wasted bits).
 * <p>
 * Key invariants:
 * - maxArray[i>>>4] holds 16 consecutive buckets, each in 2 bits.
 * - Bucket value 0 = empty; value 1-3 = (relNlz+1).
 * - absoluteNlz = stored + globalNLZ, always (no special cases).
 * - globalNLZ = -1 means nothing seen; >= 0 means all buckets have absNlz >= globalNLZ.
 * - minZeroCount tracks floor-level (stored=0 or stored=1) buckets; advances globalNLZ floor when 0.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicLogLog2 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	DynamicLogLog2(){
		this(2048, 31, -1, 0);
	}

	/** Construct from parsed command-line arguments. */
	DynamicLogLog2(Parser p){
		super(p);
		maxArray=new int[buckets>>>4];
		minZeroCount=buckets;
	}

	/**
	 * Full constructor.
	 * @param buckets_ Number of buckets (must be a power of 2 and multiple of 16)
	 * @param k_ Hash prefix length
	 * @param seed Random seed (-1 for default)
	 * @param minProb_ Minimum probability threshold
	 */
	DynamicLogLog2(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[buckets>>>4];
		minZeroCount=buckets;
	}

	/** Create an independent copy with a fresh seed. */
	@Override
	public DynamicLogLog2 copy(){return new DynamicLogLog2(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the 2-bit stored value for bucket i (0=empty, 1-3=relNlz+1). */
	private int readBucket(final int i){
		return (maxArray[i>>>4]>>>((i&15)<<1))&0x3;
	}

	/** Writes a 2-bit stored value for bucket i. val must be in [0,3]. */
	private void writeBucket(final int i, final int val){
		final int wordIdx=i>>>4;
		final int shift=(i&15)<<1;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x3<<shift))|((val&0x3)<<shift);
	}

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
		int filledCount=0;
		for(int i=0; i<buckets; i++){
			final int stored=readBucket(i);
			if(stored>0){
				final int absNlz=stored+globalNLZ;
				if(absNlz>=0 && absNlz<64){nlzCounts[absNlz+1]++; filledCount++;}
			}
		}
		nlzCounts[0]=buckets-filledCount;
		// DLL2 is counts-only: no history, luck, or mantissa bits. buckets=null.
		// nativeRelTiers=3: 2-bit DLL has 3 stored tiers (relNlz 0,1,2).
		// The top stored tier (startTier+2) loses overflows in io=t mode; dlcBest
		// applies bias correction only when it selects that specific tier.
		return new CardStats(null, nlzCounts, 0, 0, 0, 0,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				IGNORE_OVERFLOW ? 3 : Integer.MAX_VALUE,
				IGNORE_OVERFLOW ? CorrectionFactor.DLC_TIER_ERR_DLL2 : null);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		long card=Math.max(0, Math.round(s.dlcRaw()));
		card=LAZY_ALLOCATE ? Math.max(card, s.microCardinality()) : card;
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicLogLog2)log);
	}

	/**
	 * Merges another DynamicLogLog2 into this one.
	 * <p>
	 * Correct merge sequence:
	 * <ol>
	 *   <li>Re-frame both instances to newGlobalNLZ = max(this.globalNLZ, log.globalNLZ).
	 *       Buckets whose absNlz falls below the new floor correctly become empty (stored=0),
	 *       mirroring what countAndDecrement() would do one step at a time.
	 *       Note: Math.max(0,...) is correct here — NOT Math.max(1,...).
	 *       A below-floor bucket is genuinely empty in the new frame.</li>
	 *   <li>Take per-bucket max of the re-framed values.</li>
	 *   <li>Scan to recount filledBuckets and minZeroCount, then advance the floor
	 *       (possibly multiple times) if minZeroCount==0, exactly as hashAndStore does.</li>
	 * </ol>
	 * <p>
	 * <b>Multi-threaded accuracy warning:</b> With only 3 tiers, DLL2 is more sensitive
	 * to per-thread overflow than DLL4. Prefer a single synchronized instance for
	 * parallel streams.
	 */
	public void add(DynamicLogLog2 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to new relative frame: newStored = stored + (oldGlobalNLZ - newGlobalNLZ)
				final int nA=(sA==0 ? 0 : Math.max(0, Math.min(sA+(globalNLZ-newGlobalNLZ), 3)));
				final int nB=(sB==0 ? 0 : Math.max(0, Math.min(sB+(log.globalNLZ-newGlobalNLZ), 3)));
				writeBucket(i, Math.max(nA, nB));
			}
			globalNLZ=newGlobalNLZ;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int s=readBucket(i);
				if(s>0){filledBuckets++;}
				if(s==0||s==1){minZeroCount++;}
			}
			while(minZeroCount==0 && globalNLZ<wordlen){
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		if(Long.compareUnsigned(key, eeMask)>0){return;}
		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);
		final int relNlz=nlz-globalNLZ-1;

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){//Optional MicroIndex for low cardinality
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}//Allows lazy array allocation
		}

		// Stored = relNlz+1, clamped to [1,3] for overflow
		if(IGNORE_OVERFLOW && relNlz+1>3){return;} // silently ignore overflow
		final int newStored=Math.min(relNlz+1, 3);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		// oldRelNlz: 0 for both empty (stored=0) and tier-0 (stored=1).
		// minZeroCount decrements whenever a bucket leaves (empty or tier-0) for a higher tier.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		if(relNlz>oldRelNlz && oldRelNlz==0 && --minZeroCount<1){
			while(minZeroCount==0 && globalNLZ<wordlen){
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Counts buckets that will become tier-0 after decrement (currently stored=2),
	 * then decrements all non-empty stored values by 1.
	 * Only called when minZeroCount==0, guaranteeing all non-empty have stored>=2.
	 * Returns count of new tier-0 buckets (stored==1 after decrement).
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<maxArray.length; w++){
			int word=maxArray[w];
			if(word==0){continue;}
			int result=0;
			for(int b=0; b<16; b++){
				final int shift=b<<1;
				int stored=(word>>>shift)&0x3;
				if(stored>0){
					stored--; // decrement relative tier
					if(stored==1){newMinZeroCount++;} // new tier-0 after decrement
				}
				result|=(stored<<shift);
			}
			maxArray[w]=result;
		}
		return newMinZeroCount;
	}

	/** Number of buckets with stored > 0 (i.e., have seen at least one hash). */
	public int filledBuckets(){return filledBuckets;}
	/** Fraction of buckets that hold at least one observation. */
	public double occupancy(){return (double)filledBuckets/buckets;}

	/** Not used; CF correction handled via CF_MATRIX. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/** Compute all estimator values and return as a legacy-format array. */
	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 2-bit bucket array: 16 buckets per int, 0=empty, 1-3=relNlz+1. */
	private final int[] maxArray;
	/** Global floor tier: -1 means nothing seen; >= 0 means all buckets have absNlz >= globalNLZ. absNlz = stored + globalNLZ. */
	private int globalNLZ=-1;
	/** Count of floor-level (stored=0 or stored=1) buckets; triggers globalNLZ floor advance when 0. */
	private int minZeroCount;
	/** Count of buckets with stored > 0. */
	private int filledBuckets=0;
	/** Early-exit mask: filters hashes below the current global floor. */
	private long eeMask=-1L;
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	/** Diagnostic counter: hashes that pass the eeMask filter. */
	public long branch1=0, branch2=0;
	/** Rate of hashes passing eeMask relative to total added. */
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	/** Rate of hashes reaching branch2 relative to branch1. */
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum NLZ value (64-bit hash). */
	private static final int wordlen=64;
	/** Social promotion threshold (see DynamicLogLog3v2). 0=classic behavior. */
	public static int PROMOTE_THRESHOLD=0;

	/** When true, hashes that would overflow the 2-bit register are silently ignored.
	 *  Eliminates overflow corruption; the CF table absorbs resulting bias. */
	public static boolean IGNORE_OVERFLOW=false;

	/** Default resource file for DLL2 correction factors (using DLL3 table until DLL2-specific table is generated). */
	public static final String CF_FILE="?cardinalityCorrectionDLL2_iof.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the DLL2 correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

}
