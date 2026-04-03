package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicLogLog3: 3-bit packed variant of DynamicLogLog using relative NLZ storage.
 * <p>
 * Packs 10 buckets into each int, 3 bits per bucket.
 * Encoding: 0 = empty; 1-7 = (relNlz + 1), where relNlz = absoluteNlz - minZeros.
 * No mantissa — coarser precision per bucket, but allows more buckets for the same memory.
 * Overflow (relNlz >= 7) clamped to stored=7 and absorbed by CF matrix.
 * <p>
 * Memory: 3 bits/bucket × 8192 buckets = 24576 bits = 3072 bytes (~3KB).
 * 10 buckets per int; 2 bits per int wasted (bits 30-31 always zero).
 * Paper description: "3 bits per bucket."
 * <p>
 * Key invariants:
 * - maxArray[i/10] holds bucket i at bit position (i%10)*3, 3 bits wide.
 * - Bucket value 0 = empty; value 1-7 = (relNlz+1).
 * - absoluteNlz = (stored - 1) + minZeros for non-empty buckets.
 * - minZeroCount tracks empty + tier-0 (stored=1) buckets; advances minZeros floor when 0.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicLogLog3 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DynamicLogLog3(){
		this(2048, 31, -1, 0);
	}

	DynamicLogLog3(Parser p){
		super(p);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		xOverflow=buckets*Math.log(2.0*buckets)/256.0;
		overflowExpFactor=Math.exp(-xOverflow/buckets);
		storedOverflow=new int[64];
		if(FAST_COUNT){nlzCounts=new int[66];}
	}

	DynamicLogLog3(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		xOverflow=buckets*Math.log(2.0*buckets)/256.0;
		overflowExpFactor=Math.exp(-xOverflow/buckets);
		storedOverflow=new int[64];
		if(FAST_COUNT){nlzCounts=new int[66];}
	}

	@Override
	public DynamicLogLog3 copy(){return new DynamicLogLog3(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the 3-bit stored value for bucket i (0=empty, 1-7=relNlz+1). */
	private int readBucket(final int i){
		return (maxArray[i/10]>>>((i%10)*3))&0x7;
	}

	/** Writes a 3-bit stored value for bucket i. val must be in [0,7]. */
	private void writeBucket(final int i, final int val){
		final int wordIdx=i/10;
		final int shift=(i%10)*3;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x7<<shift))|((val&0x7)<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		final double rawHyb=s.hybridDLL();
		long card=(long)(rawHyb);
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicLogLog3)log);
	}

	/**
	 * Merges another DynamicLogLog3 into this one.
	 * <p>
	 * Correct merge sequence:
	 * <ol>
	 *   <li>Re-frame both instances to newMinZeros = max(this.minZeros, log.minZeros).
	 *       Buckets whose absNlz falls below the new floor correctly become empty (stored=0),
	 *       mirroring what countAndDecrement() would do one step at a time.
	 *       Note: Math.max(0,...) is correct here — NOT Math.max(1,...).
	 *       A below-floor bucket is genuinely empty in the new frame.</li>
	 *   <li>Take per-bucket max of the re-framed values.</li>
	 *   <li>Scan to recount filledBuckets and minZeroCount, then advance the floor
	 *       (possibly multiple times) if minZeroCount==0, exactly as hashAndStore does.</li>
	 * </ol>
	 * <p>
	 * <b>Multi-threaded accuracy warning:</b> Using clonal per-thread estimator copies
	 * (one DLL3 per thread, merged at the end) nondeterministically reduces accuracy
	 * and this cannot be compensated for.  The cause is architectural: each per-thread
	 * copy sees only a fraction of the data, so its sliding minZeros window promotes
	 * less aggressively than a single estimator seeing all data.  This increases net
	 * overflow events across the merged result.  With DLL3's narrow 7-tier range, the
	 * effect is severe: on a 12M-kmer dataset, t=1 gives ~10.8M, t=2 ~10.0M, t=4 ~8.2M,
	 * t=8 ~6.3M — roughly halving with each doubling of threads.
	 * For parallel use on the same stream, prefer a single synchronized instance.
	 */
	public void add(DynamicLogLog3 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to new relative frame: newStored = stored + (oldMinZeros - newMinZeros)
				final int nA=(sA==0 ? 0 : Math.max(0, Math.min(sA+(minZeros-newMinZeros), 7)));
				final int nB=(sB==0 ? 0 : Math.max(0, Math.min(sB+(log.minZeros-newMinZeros), 7)));
				writeBucket(i, Math.max(nA, nB));
			}
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
			if(FAST_COUNT){java.util.Arrays.fill(nlzCounts, 0);}
			final int phantomNlz=minZeros-1;
			for(int i=0; i<buckets; i++){
				final int s=readBucket(i);
				if(s>0){
					filledBuckets++;
					if(FAST_COUNT){final int absNlz=(s-1)+minZeros; if(absNlz<64){nlzCounts[absNlz+1]++;}}
				}else if(FAST_COUNT && minZeros>0 && phantomNlz>=0 && phantomNlz<64){
					nlzCounts[phantomNlz+1]++;
				}
				if(s==0 || (!EARLY_PROMOTE && s==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && minZeros<wordlen){
				minZeros++;
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
		final int relNlz=nlz-minZeros;
		
		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){//Optional MicroIndex for low cardinality
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS) {return;}//Allows lazy array allocation
		}

		// Stored = relNlz+1, clamped to [1,7] for overflow
		final int newStored=Math.min(relNlz+1, 7);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		if(FAST_COUNT){
			// Remove old slot from nlzCounts
			if(oldStored>0){
				final int oldAbsNlz=(oldStored-1)+minZeros;
				if(oldAbsNlz<64){nlzCounts[oldAbsNlz+1]--;}
			}else if(minZeros>0){
				final int pNlz=minZeros-1;
				if(pNlz<64){nlzCounts[pNlz+1]--;}
			}else{
				nlzCounts[0]--; // was truly empty
			}
			// Add new slot
			final int newAbsNlz=(newStored-1)+minZeros;
			if(newAbsNlz<64){nlzCounts[newAbsNlz+1]++;}
		}

		// minZeroCount decrements when a bucket leaves the tracked-zero category.
		// EARLY_PROMOTE=false (classic): tracks empty+tier-0; advances when all buckets >= 2.
		// EARLY_PROMOTE=true  (new):     tracks empty only;    advances when all buckets >= 1.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? oldStored==0 : (relNlz>oldRelNlz && oldRelNlz==0);
		if(shouldDecrement && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
				// Record overflow estimate before advancing: count buckets at stored=7
				// (the max tier, absNlz=6+minZeros). Half of those are estimated overflow victims.
				if(USE_STORED_OVERFLOW){
					final int nextCorrTier=7+minZeros;
					if(nextCorrTier<64){
						// Count buckets at stored=7 (current max tier = absNlz 6+minZeros).
						// Half are estimated overflow victims that should be in the next tier.
						// The raw topCount includes some overflow from previous eras; partially
						// adjust using half the previous correction to avoid over/undercorrection.
						int topCount=0;
						for(int i=0; i<buckets; i++){if(readBucket(i)==7){topCount++;}}
						// No cascading adjustment — raw topCount/2 only
						storedOverflow[nextCorrTier]=topCount/2;
					}
				}
				minZeros++;
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
			for(int b=0; b<10; b++){
				final int shift=b*3;
				int stored=(word>>>shift)&0x7;
				if(stored>0){
					stored--; // decrement relative tier
					if(EARLY_PROMOTE){
						if(stored==0){newMinZeroCount++; filledBuckets--;} // new empty after decrement
					}else{
						if(stored==1){newMinZeroCount++;} // new tier-0 after decrement (classic)
					}
				}
				result|=(stored<<shift);
			}
			maxArray[w]=result;
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Scans the bucket array once to populate the absolute NLZ histogram (nlzCounts),
	 * then delegates all sum computation to CardinalityStats.fromNlzCounts().
	 * <p>
	 * Phantom buckets (stored=0 when minZeros>0) are treated as absNlz = minZeros-1,
	 * one tier below the current floor. This ensures the nlzCounts distribution is
	 * identical regardless of EARLY_PROMOTE setting.
	 */
	private CardStats summarize(){
		if(!FAST_COUNT){
			if(nlzCounts==null){nlzCounts=new int[66];}
			else{java.util.Arrays.fill(nlzCounts, 0);}
			final int phantomNlz=minZeros-1;
			int filledCount=0;
			for(int i=0; i<buckets; i++){
				final int stored=readBucket(i);
				if(stored>0){
					final int absNlz=(stored-1)+minZeros;
					if(absNlz<64){nlzCounts[absNlz+1]++;}
					filledCount++;
				}else if(minZeros>0 && phantomNlz>=0 && phantomNlz<64){
					nlzCounts[phantomNlz+1]++;
					filledCount++;
				}
			}
			nlzCounts[0]=buckets-filledCount;
		}else{
			// FAST_COUNT=true: nlzCounts already maintained incrementally.
			// Recompute nlzCounts[0] to ensure empties are counted.
			int sum=0;
			for(int t=1; t<66; t++){sum+=nlzCounts[t];}
			nlzCounts[0]=buckets-sum;
		}
		final boolean doDebug=DEBUG_ONCE;
		if(doDebug){
			DEBUG_ONCE=false;
			System.err.println("DEBUG summarize(): minZeros="+minZeros+" filledBuckets="+filledBuckets+" buckets="+buckets);
			System.err.println("  CORRECT_OVERFLOW="+CORRECT_OVERFLOW+" xOverflow="+xOverflow+" expFactor="+overflowExpFactor);
			System.err.println("  nlzCounts (nonzero):");
			for(int t=0; t<66; t++){if(nlzCounts[t]>0){System.err.println("    tier "+(t-1)+": "+nlzCounts[t]);}}
		}
		final int[] counts;
		if(CORRECT_OVERFLOW && minZeros>=1){
			// Overflow for tier T accumulates BEFORE T is active; X_T is frozen when T activates.
			// Correction applies across all currently-active tiers in [max(7,minZeros), 6+minZeros].
			// Phantom tier (absNlz=minZeros-1) is excluded — it is below minZeros.
			final int lo=Math.max(7, minZeros);
			final int hi=Math.min(6+minZeros, 63);
			// Build reverse-cumulative array over nlz indices [1..65] (absNlz 0..64)
			// cumRaw[t] = sum of nlzCounts[t+1..65]  (i.e., buckets with absNlz >= t)
			final int[] cumRaw=new int[64];
			cumRaw[63]=nlzCounts[64]; // absNlz=63 stored at index 64
			for(int t=62; t>=0; t--){cumRaw[t]=cumRaw[t+1]+nlzCounts[t+1];}

			// Correct in cumulative space: each tier independently.
			// correctedCum[t] = cumRaw[t] + (B - cumRaw[t]) * (1 - exp(-X_t/B))
			final int[] corrCum=cumRaw.clone();
			for(int t=lo; t<=hi; t++){
				final double x=(USE_STORED_OVERFLOW && storedOverflow[t]>0)
					? storedOverflow[t]*OVERFLOW_SCALE : xOverflow*OVERFLOW_SCALE;
				final int addToCum=(int)Math.round(
					(buckets-cumRaw[t])*(1.0-Math.exp(-x/buckets)));
				corrCum[t]+=addToCum;
			}

			// Current-era ongoing overflow: the max tier (hi) always has ~50%
			// overflow from values with absNlz > hi clamped to stored=7.
			// Split the max tier: half stays, half goes to virtual tier hi+1.
			final int maxHi=Math.min(hi+1, 63);
			if(corrCum[hi]>0 && maxHi>hi){
				corrCum[maxHi]=corrCum[hi]/2;
			}

			// Derive corrected counts[] in the new int[66] format.
			counts=nlzCounts.clone();
			if(lo>0){counts[lo]=corrCum[lo-1]-corrCum[lo];} // absNlz=lo-1 at index lo
			for(int t=lo; t<maxHi; t++){counts[t+1]=corrCum[t]-corrCum[t+1];}
			counts[maxHi+1]=corrCum[maxHi];
			// Recompute empties
			int corrSum=0;
			for(int t=1; t<66; t++){corrSum+=counts[t];}
			counts[0]=buckets-corrSum;
			if(doDebug){
				System.err.println("  After correction [lo="+lo+" hi="+hi+"] (nonzero):");
				for(int t=0; t<66; t++){if(counts[t]>0){System.err.println("    tier "+(t-1)+": "+counts[t]);}}
			}
		}else{
			counts=nlzCounts;
		}
		lastRawNlz=nlzCounts.clone();
		lastCorrNlz=(counts==nlzCounts) ? lastRawNlz : counts;
		// DLL3 is counts-only: no history, luck, or mantissa bits. buckets=null.
		return new CardStats(null, counts, 0, 0, 0, 0,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0);
	}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 3-bit bucket array: 10 buckets per int, 0=empty, 1-7=relNlz+1. */
	private final int[] maxArray;
	/** Expected ghost overflow count per max tier: B*ln(2B)/256. Constant for given B. */
	private final double xOverflow;
	/** exp(-xOverflow/B): precomputed decay factor for overflow correction formula. */
	private final double overflowExpFactor;
	/** Per-tier stored overflow estimates. Recorded at each tier transition:
	 *  storedOverflow[t] = corrected count at tier t-1 / 2 at the moment t becomes active.
	 *  Used instead of constant xOverflow when USE_STORED_OVERFLOW=true. */
	private final int[] storedOverflow;
	private int minZeros=0;
	/** Count of (empty + tier-0) buckets; triggers minZeros floor advance when 0. */
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	public int getMinZeros(){return minZeros;}
	public int[] getStoredOverflow(){return storedOverflow;}
	/** Last raw and corrected nlzCounts from summarize(). Set after each rawEstimates() call. */
	int[] lastRawNlz, lastCorrNlz;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** When true, nlzCounts is maintained incrementally in hashAndStore() rather than rebuilt
	 *  in summarize(). Eliminates the O(buckets) scan per rawEstimates() call — ~32x speedup
	 *  for CF table generation. Disabled in production (false) to avoid the per-add overhead. */
	public static final boolean FAST_COUNT=false;
	/** Social promotion threshold — see DynamicLogLog3v2 for implementation.
	 * In DLL3, this field is parsed but has no effect (DLL3 uses classic promotion). */
	public static int PROMOTE_THRESHOLD=0;
	/** When true, advance tier as soon as all buckets are nonzero (stored>=1) rather than >=2.
	 * At advance, subtracts 1 from all buckets, which may reset some to empty (stored=0).
	 * This is safe because lcMin (tier-compensated LC) handles post-advance zero buckets correctly.
	 * Reduces tier-7+ overflow pollution in DLL3; requires new CF table generation when changed. */
	public static boolean EARLY_PROMOTE=true;
	/** When true, corrects the systematic underestimate caused by 3-bit relNlz overflow
	 *  clamping at tier 7+. Uses cumulative-space Poisson correction with per-DDL stored
	 *  overflow estimates. Each tier is corrected independently in reverse-cumulative space. */
	public static boolean CORRECT_OVERFLOW=true;
	/** Scale factor for overflow correction boost.  1.0 is theoretically correct and
	 *  may be best for gaussian weighting (mode 3), but for error-rate-based blending
	 *  (mode 2), 1.7 minimizes peak error and 1.8 minimizes terminal error.
	 *  Found empirically via agent investigation. */
	public static double OVERFLOW_SCALE=1.7;
	/** When true, use per-tier stored overflow estimates (recorded at tier transitions)
	 *  instead of the constant xOverflow formula. More accurate because it uses actual state. */
	public static boolean USE_STORED_OVERFLOW=true;

	/** Set to true externally to trigger one debug dump from summarize(), then auto-clears. */
	public static boolean DEBUG_ONCE=false;

	/** Default resource file for DLL3 correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDLL3.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=8192;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the DLL3 correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

}
