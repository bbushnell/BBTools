package cardinality;

import shared.Tools;

/**
 * DualHashDynamicLogLog4: dual-hash variant of DynamicLogLog4.
 * <p>
 * Uses max(hash1, hash2) to suppress extreme NLZ outliers.
 * NLZ(max(h1,h2)) = min(NLZ(h1), NLZ(h2)), so P(NLZ >= k) = 4^{-k}
 * instead of 2^{-k}.  This compresses the NLZ tail, dramatically
 * reducing overflow with small exponents (3-bit or even 2-bit).
 * <p>
 * For DLC estimation, the tier multiplier changes from 2^tier to 4^tier.
 * This is handled by setting AbstractCardStats.TIER_SCALE=2 during
 * summarize/rawEstimates calls.
 * <p>
 * This is a proof-of-concept using 4-bit exponents (same as DLL4).
 * Once validated, the same technique can be applied to 3-bit and 2-bit
 * variants where overflow is a real problem.
 *
 * @author Brian Bushnell, Chloe
 * @date April 2026
 */
public final class DualHashDynamicLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DualHashDynamicLogLog4(){
		this(2048, 31, -1, 0);
	}

	DualHashDynamicLogLog4(parse.Parser p){
		super(p);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
	}

	DualHashDynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
	}

	@Override
	public DualHashDynamicLogLog4 copy(){return new DualHashDynamicLogLog4(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	private int readBucket(final int i){
		return (maxArray[i>>>3]>>>((i&7)<<2))&0xF;
	}

	private void writeBucket(final int i, final int val){
		final int wordIdx=i>>>3;
		final int shift=(i&7)<<2;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0xF<<shift))|((val&0xF)<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
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
		lastRawNlz=nlzCounts.clone();
		// Set TIER_SCALE=2 for dual-hash estimation (4^tier instead of 2^tier)
		final double savedScale=AbstractCardStats.TIER_SCALE;
		AbstractCardStats.TIER_SCALE=2;
		try{
			return new CardStats(null, nlzCounts, 0, 0, 0, 0,
					buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
					terminalMeanCF(), terminalMeanPlusCF());
		}finally{
			AbstractCardStats.TIER_SCALE=savedScale;
		}
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
		add((DualHashDynamicLogLog4)log);
	}

	public void add(DualHashDynamicLogLog4 log){
		added+=log.added;
		branch1+=log.branch1;
		branch2+=log.branch2;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				final int nA=(sA==0 ? 0 : Math.max(0, Math.min(sA+(minZeros-newMinZeros), 15)));
				final int nB=(sB==0 ? 0 : Math.max(0, Math.min(sB+(log.minZeros-newMinZeros), 15)));
				writeBucket(i, Math.max(nA, nB));
			}
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int s=readBucket(i);
				if(s>0){filledBuckets++;}
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
		// First hash
		final long rawKey1=number^hashXor;
		final long key1=Tools.hash64shift(rawKey1);

		// First early exit: if key1 has too few leading zeros, skip.
		// max(key1,key2) >= key1, so if key1 fails, max also fails.
		if(Long.compareUnsigned(key1, eeMask)>0){return;}
		branch1++;

		// Second hash with different seed
		final long rawKey2=number^(hashXor+HASH_OFFSET);
		final long key2=Tools.hash64shift(rawKey2);

		// Unsigned max of two hashes -> NLZ = min(NLZ1, NLZ2)
		final long key=(Long.compareUnsigned(key1, key2)>=0 ? key1 : key2);

		// Second early exit on the max
		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int nlz=Long.numberOfLeadingZeros(key);
		// Bucket from key1 (consistent element-to-bucket mapping)
		final int bucket=(int)(key1&bucketMask);
		final int relNlz=nlz-minZeros;

		// MicroIndex from key1 (element identity, not NLZ value)
		final long micro=(key1>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		if(IGNORE_OVERFLOW && relNlz+1>15){return;}
		final int newStored=Math.min(relNlz+1, 15);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		branch2++;
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

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

	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int wordIdx=0; wordIdx<maxArray.length; wordIdx++){
			final int oldWord=maxArray[wordIdx];
			if(oldWord==0){continue;}
			int newWord=0;
			for(int bucketIdx=0; bucketIdx<8; bucketIdx++){
				final int bucketOffset=bucketIdx<<2;
				int bucketVal=(oldWord>>>bucketOffset)&0xF;
				if(bucketVal>0){
					bucketVal--;
					if(EARLY_PROMOTE){
						if(bucketVal==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(bucketVal==1){newMinZeroCount++;}
					}
				}
				newWord|=(bucketVal<<bucketOffset);
			}
			maxArray[wordIdx]=newWord;
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		double[] legacy=AbstractCardStats.buildLegacyArray(s, hybridEst);
		// Pad to include WordEst/WordEstCV slots (unused, 0.0) for calibration driver compatibility
		double[] ext=new double[legacy.length+2];
		System.arraycopy(legacy, 0, ext, 0, legacy.length);
		return ext;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] maxArray;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return minZeros;}
	int[] lastRawNlz;
	private int[] nlzCounts;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Second hash seed offset: hash2 = hash(number ^ (hashXor + HASH_OFFSET)). */
	private static final long HASH_OFFSET=123456789L;

	public static boolean EARLY_PROMOTE=true;
	public static boolean IGNORE_OVERFLOW=false;

	public static final String CF_FILE="?cardinalityCorrectionDHDLL4.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null; // No CF table yet
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("DHDLL4: No CF table found (expected for new class); running without CF.");
			return CF_MATRIX=null;
		}
	}

	/** Stub: measure from preliminary CF table, then replace 1f with actual ratio. */
	@Override public float terminalMeanCF(){return 1f;}

}
