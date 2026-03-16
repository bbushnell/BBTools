package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicLogLog4: 4-bit packed variant of DynamicLogLog using relative NLZ storage.
 * <p>
 * Packs 8 buckets into each int, 4 bits per bucket (nibble packing).
 * Encoding: 0 = empty; 1-15 = (relNlz + 1), where relNlz = absoluteNlz - minZeros.
 * No mantissa — coarser precision per bucket, but allows 2x bucket count for the same memory.
 * Overflow (relNlz >= 15) clamped to stored=15 and absorbed by CF matrix.
 * <p>
 * Key question: does DLL4 with 4096 buckets beat DLL2 with 2048 buckets for the same
 * 4KB memory footprint? (4096 * 4 bits = 2048 * 8 bits = 2048 chars = 4KB)
 * <p>
 * Key invariants:
 * - maxArray[i>>>3] holds 8 consecutive buckets, each in 4 bits.
 * - Bucket value 0 = empty; value 1-15 = (relNlz+1).
 * - absoluteNlz = (stored - 1) + minZeros for non-empty buckets.
 * - minZeroCount tracks empty + tier-0 (stored=1) buckets; advances minZeros floor when 0.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DynamicLogLog4(){
		this(2048, 31, -1, 0);
	}

	DynamicLogLog4(Parser p){
		super(p);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
		if(FAST_COUNT){nlzCounts=new int[64];}
	}

	DynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
		if(FAST_COUNT){nlzCounts=new int[64];}
	}

	@Override
	public DynamicLogLog4 copy(){return new DynamicLogLog4(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the 4-bit stored value for bucket i (0=empty, 1-15=relNlz+1). */
	private int readBucket(final int i){
		return (maxArray[i>>>3]>>>((i&7)<<2))&0xF;
	}

	/** Writes a 4-bit stored value for bucket i. val must be in [0,15]. */
	private void writeBucket(final int i, final int val){
		final int wordIdx=i>>>3;
		final int shift=(i&7)<<2;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0xF<<shift))|((val&0xF)<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Scans the bucket array once to populate the absolute NLZ histogram (nlzCounts),
	 * then delegates all sum computation to CardinalityStats.fromNlzCounts().
	 * <p>
	 * Phantom buckets (stored=0 when minZeros>0) are treated as absNlz = minZeros-1,
	 * one tier below the current floor. This ensures the nlzCounts distribution is
	 * identical regardless of EARLY_PROMOTE setting.
	 */
	private CardinalityStats summarize(){
		if(!FAST_COUNT){
			if(nlzCounts==null){nlzCounts=new int[64];}
			else{java.util.Arrays.fill(nlzCounts, 0);}
			final int phantomNlz=minZeros-1;
			for(int i=0; i<buckets; i++){
				final int stored=readBucket(i);
				if(stored>0){
					final int absNlz=(stored-1)+minZeros;
					if(absNlz<64){nlzCounts[absNlz]++;}
				}else if(minZeros>0 && phantomNlz<64){
					nlzCounts[phantomNlz]++;
				}
			}
		}
		// FAST_COUNT=true: nlzCounts already maintained incrementally; use directly.
		lastRawNlz=nlzCounts.clone();
		lastCorrNlz=lastRawNlz; // DLL4 has no overflow correction
		return CardinalityStats.fromNlzCounts(nlzCounts, buckets, microIndex,
		                                      CF_MATRIX, CF_BUCKETS,
		                                      CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardinalityStats s=summarize();
		long card=(long)s.hybridDLL();
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicLogLog4)log);
	}

	/**
	 * Merges another DynamicLogLog4 into this one.
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
	 * (one DLL4 per thread, merged at the end) nondeterministically reduces accuracy
	 * and this cannot be compensated for.  The cause is architectural: each per-thread
	 * copy sees only a fraction of the data, so its sliding minZeros window promotes
	 * less aggressively than a single estimator seeing all data.  This increases net
	 * overflow events across the merged result.  With DLL4's 15-tier range the effect
	 * is minor in practice: on a 12M-kmer dataset, t=1 gives ~11.43M, t=2 ~11.43M,
	 * t=4 ~11.43M, t=8 ~11.42M — well under 0.1% degradation even at 8 threads.
	 * For parallel use on the same stream, prefer a single synchronized instance.
	 */
	public void add(DynamicLogLog4 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to new relative frame: newStored = stored + (oldMinZeros - newMinZeros)
				final int nA=(sA==0 ? 0 : Math.max(0, Math.min(sA+(minZeros-newMinZeros), 15)));
				final int nB=(sB==0 ? 0 : Math.max(0, Math.min(sB+(log.minZeros-newMinZeros), 15)));
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
					if(FAST_COUNT){final int absNlz=(s-1)+minZeros; if(absNlz<64){nlzCounts[absNlz]++;}}
				}else if(FAST_COUNT && minZeros>0 && phantomNlz<64){
					nlzCounts[phantomNlz]++;
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
		if(USE_MICRO){//Optional MicroIndex for low cardinality
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS) {return;}//Allows lazy array allocation
		}

		// Stored = relNlz+1, clamped to [1,15] for overflow
		final int newStored=Math.min(relNlz+1, 15);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		//		branch2++;
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		if(FAST_COUNT){
			// Remove old slot from nlzCounts
			if(oldStored>0){
				final int oldAbsNlz=(oldStored-1)+minZeros;
				if(oldAbsNlz<64){nlzCounts[oldAbsNlz]--;}
			}else if(minZeros>0){
				final int pNlz=minZeros-1;
				if(pNlz<64){nlzCounts[pNlz]--;}
			}
			// Add new slot
			final int newAbsNlz=(newStored-1)+minZeros;
			if(newAbsNlz<64){nlzCounts[newAbsNlz]++;}
		}

		// minZeroCount decrements when a bucket leaves the tracked-zero category.
		// EARLY_PROMOTE=false (classic): tracks empty+tier-0; advances when all buckets >= 2.
		// EARLY_PROMOTE=true  (new):     tracks empty only;    advances when all buckets >= 1.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? oldStored==0 : (relNlz>oldRelNlz && oldRelNlz==0);
		if(shouldDecrement && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
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
			for(int b=0; b<8; b++){
				final int shift=b<<2;
				int stored=(word>>>shift)&0xF;
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

	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(Math.max(s.hybridDLL(), s.microCardinality()));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 4-bit bucket array: 8 buckets per int, 0=empty, 1-15=relNlz+1. */
	private final int[] maxArray;
	private int minZeros=0;
	/** Count of (empty + tier-0) buckets; triggers minZeros floor advance when 0. */
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	public int getMinZeros(){return minZeros;}
	/** Last raw and corrected nlzCounts from summarize(). */
	int[] lastRawNlz, lastCorrNlz;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, added);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** When true, nlzCounts is maintained incrementally in hashAndStore() rather than rebuilt
	 *  in summarize(). Eliminates the O(buckets) scan per rawEstimates() call — ~32x speedup
	 *  for CF table generation. Disabled in production (false) to avoid the per-add overhead. */
	public static final boolean FAST_COUNT=true;
	/** Social promotion threshold (see DynamicLogLog3v2). 0=classic behavior. */
	public static int PROMOTE_THRESHOLD=0;
	/** When true, advance tier as soon as all buckets are nonzero (stored>=1) rather than >=2.
	 * At advance, subtracts 1 from all buckets, which may reset some to empty (stored=0).
	 * This is safe because lcMin (tier-compensated LC) handles post-advance zero buckets correctly.
	 * Reduces tier-15+ overflow pollution; requires new CF table generation when changed. */
	public static boolean EARLY_PROMOTE=true;


	/** Default resource file for DDL correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDLL4.tsv.gz";
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
