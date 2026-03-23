package cardinality;

import shared.Tools;

/**
 * DynamicLogLog4m: DLL4 variant supporting arbitrary (non-power-of-2) bucket counts.
 * Uses modulo bucket selection (hash % buckets) instead of bit masking.
 * Paper-only: for equal-memory comparisons against UDLL6 at non-standard bucket counts.
 * <p>
 * Uses byte[] storage (1 register per byte, 4 bits used) for simplicity.
 * No merge support. No FAST_COUNT. No int packing.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicLogLog4m extends CardinalityTracker {

	/** Actual number of buckets (may be non-power-of-2). */
	private final int modBuckets;
	private final byte[] registers;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	DynamicLogLog4m(){this(2560, 31, -1, 0);}
	DynamicLogLog4m(int buckets_, int k_, long seed, float minProb_){
		super(Integer.highestOneBit(buckets_-1)<<1, k_, seed, minProb_); // next pow2
		modBuckets=buckets_;
		registers=new byte[modBuckets];
		minZeroCount=modBuckets;
	}
	@Override public DynamicLogLog4m copy(){return new DynamicLogLog4m(modBuckets, k, -1, minProb);}

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);

		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
		final int relNlz=nlz-minZeros;

		// MicroIndex from upper hash bits
		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int newStored=Math.min(relNlz+1, 15);
		final int oldStored=registers[bucket]&0xFF;

		if(newStored<=oldStored){return;}
		lastCardinality=-1;

		registers[bucket]=(byte)newStored;
		if(oldStored==0){filledBuckets++;}

		// EARLY_PROMOTE: advance when all buckets are non-empty
		if(oldStored==0 && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

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

	private CardinalityStats summarize(){
		final int[] nlzCounts=new int[64];
		final int phantomNlz=minZeros-1;
		for(int i=0; i<modBuckets; i++){
			final int stored=registers[i]&0xFF;
			if(stored>0){
				final int absNlz=(stored-1)+minZeros;
				if(absNlz<64){nlzCounts[absNlz]++;}
			}else if(minZeros>0 && phantomNlz>=0 && phantomNlz<64){
				nlzCounts[phantomNlz]++;
			}
		}
		return CardinalityStats.fromNlzCounts(nlzCounts, modBuckets, microIndex,
		                                      CF_MATRIX, CF_BUCKETS,
		                                      CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardinalityStats s=summarize();
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

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/modBuckets;}
	public int getModBuckets(){return modBuckets;}
	public int getMinZeros(){return minZeros;}

	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(Math.max(s.hybridDLL(), s.microCardinality()));
	}

	private static final int wordlen=64;

	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=DynamicLogLog4.initializeCF(CF_BUCKETS);

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}
}
