package cardinality;

import shared.Tools;

/**
 * CompressedDynamicLogLog3: 3-bit DLL variant with compressed NLZ tiers.
 * <p>
 * Packs 10 buckets per int at 3 bits each (7 usable levels + phantom).
 * Compresses the NLZ tail so that 3 storage bits cover a wider dynamic
 * range than a raw 3-bit single-hash register would.  Two compression
 * modes are supported, selected by the static final {@link #DUAL}:
 * <p>
 * <b>Mantissa mode (DUAL=false, current default):</b> One hash per
 * element.  A mantissa bit is derived from the top 16 bits after the
 * leading 1, thresholded at (2-sqrt(2)) so P(mantissa=1) = sqrt(2)-1.
 * halfNlz = 2*rawNlz + mantissa, tier = halfNlz/3.  Each tier divides
 * probability by 2*sqrt(2), giving effective dynamic range
 * (2*sqrt(2))^7 ~= 1448 — equivalent to ~10.5 bits of single-hash NLZ.
 * P(overflow per above-floor element) = (2*sqrt(2))^-7 ~= 7e-4.
 * TIER_SCALE=1.5 for DLC math.
 * <p>
 * <b>Dual-hash mode (DUAL=true):</b> Two independent hashes per element.
 * NLZ(max(h1,h2)) = min(NLZ(h1), NLZ(h2)), so P(NLZ >= k) = 4^-k
 * instead of 2^-k.  Equivalent in distribution to simply dividing a
 * single-hash NLZ by 2 — each stored level covers 2 NLZ bits.  7 levels
 * give dynamic range 4^7 = 16384.  P(overflow per above-floor element)
 * = 4^-7 ~= 6e-5, astronomically rare; no overflow correction needed.
 * TIER_SCALE=2 for DLC math.
 * <p>
 * Memory: 3 bits/bucket — 33% more buckets than DLL4 for the same memory.
 * Overflow-correction machinery ({@link #CORRECT_OVERFLOW},
 * storedOverflow[], xOverflow) is present but off by default; it is
 * primarily useful in mantissa mode at very high cardinality.
 *
 * @author Brian Bushnell, Chloe
 * @date April 2026
 */
public final class CompressedDynamicLogLog3 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	CompressedDynamicLogLog3(){
		this(2048, 31, -1, 0);
	}

	CompressedDynamicLogLog3(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	CompressedDynamicLogLog3(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*10<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*10-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*10;
		modBuckets=rounded>0 ? rounded : buckets;
		words=modBuckets/10;
		maxArray=new int[words];
		minZeroCount=modBuckets;
		xOverflow=modBuckets*Math.log(2.0*modBuckets)/256.0;
		storedOverflow=new int[64];
	}

	private static int roundToWords(int b){return Math.max(1, (b+5)/10);}

	@Override
	public CompressedDynamicLogLog3 copy(){return new CompressedDynamicLogLog3(modBuckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	private int readBucket(final int i){
		return (maxArray[i/10]>>>((i%10)*3))&0x7;
	}

	private void writeBucket(final int i, final int val){
		final int wordIdx=i/10;
		final int shift=(i%10)*3;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x7<<shift))|((val&0x7)<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		final int phantomNlz=minZeros-1;
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
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
		nlzCounts[0]=modBuckets-filledCount;

		final int[] counts;
		if(CORRECT_OVERFLOW && minZeros>=1){
			final int lo=Math.max(7, minZeros);
			final int hi=Math.min(6+minZeros, 63);
			final int[] cumRaw=new int[64];
			cumRaw[63]=nlzCounts[64];
			for(int t=62; t>=0; t--){cumRaw[t]=cumRaw[t+1]+nlzCounts[t+1];}
			final int[] corrCum=cumRaw.clone();
			for(int t=lo; t<=hi; t++){
				final double x=(storedOverflow[t]>0)
					? storedOverflow[t]*OVERFLOW_SCALE : xOverflow*OVERFLOW_SCALE;
				final int addToCum=(int)Math.round(
					(modBuckets-cumRaw[t])*(1.0-Math.exp(-x/modBuckets)));
				corrCum[t]+=addToCum;
			}
			final int maxHi=Math.min(hi+1, 63);
			if(corrCum[hi]>0 && maxHi>hi){
				corrCum[maxHi]=corrCum[hi]/2;
			}
			counts=nlzCounts.clone();
			if(lo>0){counts[lo]=corrCum[lo-1]-corrCum[lo];}
			for(int t=lo; t<maxHi; t++){counts[t+1]=corrCum[t]-corrCum[t+1];}
			counts[maxHi+1]=corrCum[maxHi];
			int corrSum=0;
			for(int t=1; t<66; t++){corrSum+=counts[t];}
			counts[0]=modBuckets-corrSum;
		}else{
			counts=nlzCounts;
		}

		return new CardStats(null, counts, 0, 0, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(),
				(DUAL ? 2.0 : 1.5));
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
		add((CompressedDynamicLogLog3)log);
	}

	public void add(CompressedDynamicLogLog3 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<modBuckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				final int nA=(sA==0 ? 0 : Math.max(0, Math.min(sA+(minZeros-newMinZeros), 7)));
				final int nB=(sB==0 ? 0 : Math.max(0, Math.min(sB+(log.minZeros-newMinZeros), 7)));
				writeBucket(i, Math.max(nA, nB));
			}
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
				final int s=readBucket(i);
				if(s>0){filledBuckets++;}
				if(s==0 || (!EARLY_PROMOTE && s==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && minZeros<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final int nlz;
		final int bucket;
		final long micro;

		if(DUAL){
			// --- Dual-hash mode: max(hash1, hash2), 4x per tier ---
			final long rawKey1=number^hashXor;
			final long key1=Tools.hash64shift(rawKey1);
			if(Long.compareUnsigned(key1, eeMask)>0){return;}

			final long rawKey2=number^(hashXor+HASH_OFFSET);
			final long key2=Tools.hash64shift(rawKey2);
			final long key=(Long.compareUnsigned(key1, key2)>=0 ? key1 : key2);
			if(Long.compareUnsigned(key, eeMask)>0){return;}

			nlz=Long.numberOfLeadingZeros(key);
			bucket=(int)(Long.remainderUnsigned(key1, modBuckets));
			micro=(key1>>>58)&0x3FL;
		}else{
			// --- Mantissa mode: single hash, half-NLZ tiers, 2*sqrt(2) per tier ---
			final long rawKey=number^hashXor;
			final long key=Tools.hash64shift(rawKey);
			if(Long.compareUnsigned(key, eeMask)>0){return;}

			final int rawNlz=Long.numberOfLeadingZeros(key);
			// Mantissa: top 20 bits after leading 1, thresholded at (2-sqrt(2))
			// to give each half-NLZ step exactly sqrt(2) probability ratio.
			// Threshold: (2-sqrt(2)) * 1048576 = 614242 (0.58579...)
			final int mBits;
			if(rawNlz>=43){
				mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
			}else{
				mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
			}
			final int mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
			// halfNlz encodes NLZ in half-steps; tier = halfNlz/3
			nlz=(2*rawNlz+mantissa)/3;  // absolute tier number
			bucket=(int)(Long.remainderUnsigned(key, modBuckets));
			micro=(key>>>58)&0x3FL;
		}

		final int relNlz=nlz-minZeros;
		if(relNlz<0){return;}

		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Stored = relNlz+1, clamped to [1,7]
		if(IGNORE_OVERFLOW && relNlz+1>7){return;}
		final int newStored=Math.min(relNlz+1, 7);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? oldStored==0 : (relNlz>oldRelNlz && oldRelNlz==0);
		if(shouldDecrement && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
				if(CORRECT_OVERFLOW){
					final int nextCorrTier=7+minZeros;
					if(nextCorrTier<64){
						int topCount=0;
						for(int i=0; i<modBuckets; i++){if(readBucket(i)==7){topCount++;}}
						storedOverflow[nextCorrTier]=topCount/2;
					}
				}
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	/** Advance the global floor by one tier and update eeMask. */
	private void advanceFloor(){
		if(DUAL){
			minZeros++;
			eeMask>>>=1;
		}else{
			// Mantissa mode: NLZ threshold per tier is (3*tier)/2 (integer div)
			// Shift alternates +1, +2 as tiers advance
			final int nlzShift=(minZeros%2==0) ? 1 : 2;
			minZeros++;
			eeMask>>>=nlzShift;
		}
	}

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
					stored--;
					if(EARLY_PROMOTE){
						if(stored==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(stored==1){newMinZeroCount++;}
					}
				}
				result|=(stored<<shift);
			}
			maxArray[w]=result;
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/modBuckets;}
	@Override public int actualBuckets(){return modBuckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		lastSummarized=s;
		final double hybridEst=s.hybridDLL();
		double[] legacy=AbstractCardStats.buildLegacyArray(s, hybridEst);
		double[] ext=new double[legacy.length+2];
		System.arraycopy(legacy, 0, ext, 0, legacy.length);
		return ext;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] maxArray;
	private final double xOverflow;
	private final int[] storedOverflow;
	private final int modBuckets;
	private final int words;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return minZeros;}
	int[] getNlzSnapshot(){summarize(); return nlzCounts.clone();}
	int readBucketAt(int i){return readBucket(i);}
	private int[] nlzCounts;
	private CardStats lastSummarized;

	public CardStats consumeLastSummarized(){
		final CardStats cs=lastSummarized;
		lastSummarized=null;
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	private static final long HASH_OFFSET=123456789L;
	/** Mantissa threshold for log-uniform half-NLZ steps: (2-sqrt(2)) * 1048576 */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	public static boolean EARLY_PROMOTE=true;
	public static boolean CORRECT_OVERFLOW=false;
	public static double OVERFLOW_SCALE=1.7;
	public static boolean IGNORE_OVERFLOW=false;

	/** When true, use dual-hash max (4x per tier).
	 *  When false, use single-hash with mantissa bit (2*sqrt(2) per tier). */
	public static final boolean DUAL=false;

	public static final String CF_FILE="?cardinalityCorrectionCDLL4.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null; // No CF table yet
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("CDLL3: No CF table found (expected for new class); running without CF.");
			return CF_MATRIX=null;
		}
	}

	@Override public float terminalMeanCF(){return 0.883441f;}
	@Override public float terminalMeanPlusCF(){return 1.0f;}

}
