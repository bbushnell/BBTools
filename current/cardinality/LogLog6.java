package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * LogLog6: Classic 6-bit LogLog estimator with no tier promotion.
 * <p>
 * Stores one 6-bit absolute NLZ value per bucket in a byte array.
 * Encoding: 0 = empty; 1-63 = (absNlz + 1). Max representable absNlz = 62.
 * <p>
 * This is the "fair" HLL baseline for the DLL paper: same estimator math
 * (nlzCounts-based, same CardinalityStats pipeline) but without DLL's innovations
 * (no tier promotion, no eeMask early exit, 6 bits instead of 4).
 * Memory: 1 byte per bucket (2 bits wasted). 2048 buckets = 2KB.
 * <p>
 * Compared to DLL4 at the same memory budget:
 * <ul>
 *   <li>LL6 with 2048 buckets = 2KB = DLL4 with 4096 buckets (4 bits × 4096 = 2KB)</li>
 *   <li>DLL4 gets 2× the buckets for the same memory, reducing variance by √2</li>
 *   <li>DLL4 has eeMask early exit (~70% of hashes skipped at high cardinality)</li>
 *   <li>LL6 has no tier overflow (63 tiers vs DLL4's 15), but this rarely matters</li>
 * </ul>
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class LogLog6 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	LogLog6(){
		this(2048, 31, -1, 0);
	}

	LogLog6(Parser p){
		super(p);
		maxArray=new byte[buckets];
	}

	LogLog6(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
	}

	@Override
	public LogLog6 copy(){return new LogLog6(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Scans the bucket array to populate the absolute NLZ histogram (nlzCounts),
	 * then delegates to CardinalityStats.fromNlzCounts().
	 * No tier promotion: stored value IS the absolute encoding (0=empty, 1-63=absNlz+1).
	 */
	private CardStats summarize(){
		nlzCounts=new int[66];
		int filledCount=0;
		for(int i=0; i<buckets; i++){
			final int stored=maxArray[i]&0xFF;
			if(stored>0){
				final int absNlz=stored-1;
				if(absNlz<64){nlzCounts[absNlz+1]++;}
				filledCount++;
			}
		}
		nlzCounts[0]=buckets-filledCount;
		lastRawNlz=nlzCounts;
		lastCorrNlz=nlzCounts;
		// microIndex=0: LL6 doesn't use micro floor (keeps HLL comparison fair)
		return new CardStats(null, nlzCounts, 0, 0, 0, 0,
				buckets, 0, added, CF_MATRIX, CF_BUCKETS, 0);
	}

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
		add((LogLog6)log);
	}

	/** Merges another LogLog6 into this one via per-bucket max. */
	public void add(LogLog6 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				final int a=maxArray[i]&0xFF;
				final int b=log.maxArray[i]&0xFF;
				if(b>a){maxArray[i]=(byte)b;}
			}
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		// No eeMask early exit — every hash is processed (classic HLL behavior)

		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Stored = absNlz+1, clamped to [1,63]
		final int newStored=Math.min(nlz+1, 63);
		final int oldStored=maxArray[bucket]&0xFF;

		if(newStored<=oldStored){return;}
		lastCardinality=-1;
		if(oldStored==0){filledBuckets++;}
		maxArray[bucket]=(byte)newStored;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** One byte per bucket: 0=empty, 1-63=absNlz+1. 2 bits wasted per byte. */
	private final byte[] maxArray;
	private int filledBuckets=0;

	/** Last raw nlzCounts from summarize(). */
	int[] lastRawNlz, lastCorrNlz;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Default resource file for LL6 correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionLL6.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

}
