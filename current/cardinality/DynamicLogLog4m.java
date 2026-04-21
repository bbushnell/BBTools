package cardinality;

import shared.Tools;

/**
 * DynamicLogLog4m: DLL4 variant supporting arbitrary (non-power-of-2) bucket counts.
 * Uses modulo bucket selection (hash % buckets) instead of bit masking.
 * Paper-only: for equal-memory comparisons against UDLL6 at non-standard bucket counts.
 * <p>
 * Uses byte[] storage (1 register per byte, 4 bits used) for simplicity.
 * No merge support. No FAST_COUNT. No int packing.
 * <p>
 * Encoding: stored=0..14 = relNlz 0..14; stored=15 = overflow clamp.
 * globalNLZ = -1 means nothing seen; >= 0 means floor.
 * absoluteNlz = stored + globalNLZ, always.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicLogLog4m extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2560 buckets, k=31. */
	DynamicLogLog4m(){this(2560, 31, -1, 0);}

	/**
	 * Full constructor. Rounds super.buckets to next power of 2 for hash masking,
	 * but uses modBuckets (exact) for all register operations.
	 * @param buckets_ Actual number of buckets (may be non-power-of-2)
	 * @param k_ Hash prefix length
	 * @param seed Random seed (-1 for default)
	 * @param minProb_ Minimum probability threshold
	 */
	DynamicLogLog4m(int buckets_, int k_, long seed, float minProb_){
		super(Integer.highestOneBit(buckets_-1)<<1, k_, seed, minProb_); // next pow2
		modBuckets=buckets_;
		registers=new byte[modBuckets];
		minZeroCount=modBuckets;
	}

	/** Create an independent copy with a fresh seed. */
	@Override public DynamicLogLog4m copy(){return new DynamicLogLog4m(modBuckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);

		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets)); // modulo for non-power-of-2
		final int relNlz=nlz-globalNLZ;

		// MicroIndex from upper hash bits
		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int newStored=Math.min(relNlz, 15);
		final int oldStored=registers[bucket]&0xFF;

		if(newStored<=oldStored){return;}
		lastCardinality=-1;

		registers[bucket]=(byte)newStored;
		if(oldStored==0){filledBuckets++;}

		// EARLY_PROMOTE: advance when all buckets are non-empty
		if(oldStored==0 && --minZeroCount<1){
			while(minZeroCount==0 && globalNLZ<wordlen){
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Decrements all buckets with value at least 1 after a global floor advance.
	 * Returns the count of buckets that dropped to stored=0 (new floor-level).
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int i=0; i<modBuckets; i++){
			int stored=registers[i]&0xFF;
			if(stored>0){
				stored--;
				if(stored==0){newMinZeroCount++; filledBuckets--;}
			}
			registers[i]=(byte)stored;
		}
		return newMinZeroCount;
	}

	/**
	 * Build NLZ histogram for CardStats.
	 * absoluteNlz = stored + globalNLZ, always.
	 * Buckets with absNlz < 0 or >= 64 are excluded (counted as empty in nlzCounts[0]).
	 * @return CardStats with all estimator values computed
	 */
	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int absNlz=(registers[i]&0xFF)+globalNLZ;
			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;
		// DLL4m is counts-only: no history, luck, or mantissa bits. buckets=null.
		return new CardStats(null, nlzCounts, 0, 0, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		final double rawHyb=s.hybridDLL();
		long card=(long)(rawHyb*s.cf(rawHyb, CorrectionFactor.HYBRID));
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){throw new UnsupportedOperationException();}
	@Override public final float[] compensationFactorLogBucketsArray(){return null;}

	/** Number of buckets with stored > 0 (non-floor). */
	public int filledBuckets(){return filledBuckets;}
	/** Fraction of buckets with stored > 0. */
	public double occupancy(){return (double)filledBuckets/modBuckets;}
	/** Actual bucket count (may be non-power-of-2). */
	public int getModBuckets(){return modBuckets;}
	/** Compatibility accessor: returns globalNLZ+1 to match legacy minZeros convention. */
	public int getMinZeros(){return globalNLZ+1;}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Actual number of buckets (may be non-power-of-2). */
	private final int modBuckets;
	/** One byte per bucket; stored value is in [0,15] (upper nibble unused). */
	private final byte[] registers;
	/** Global NLZ floor. -1 = nothing seen; >= 0 = all buckets have absNlz >= globalNLZ. */
	private int globalNLZ=-1;
	/** Count of floor-level (stored=0) buckets; triggers globalNLZ floor advance when 0. */
	private int minZeroCount;
	/** Count of buckets with stored > 0. */
	private int filledBuckets=0;
	/** Early-exit mask: hashes above this value are skipped. */
	private long eeMask=-1L;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;

	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; shared with DynamicLogLog4. */
	private static float[][] CF_MATRIX=DynamicLogLog4.initializeCF(CF_BUCKETS);

}
