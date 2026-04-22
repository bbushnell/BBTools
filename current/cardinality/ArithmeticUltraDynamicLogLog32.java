package cardinality;

import shared.Tools;

/**
 * ArithmeticUltraDynamicLogLog32: full arithmetic encoding of
 * standard-tier cardinality estimator with 3-state history.
 * <p>
 * Word layout (32-bit int, 6 buckets per word):
 *   No bit boundaries. Each bucket encodes one of 39 states
 *   (13 exponent × 3 history) packed via mixed-radix arithmetic:
 *   word = b0 + 39*b1 + 39^2*b2 + 39^3*b3 + 39^4*b4 + 39^5*b5
 *   39^6 = 3,518,743,761 < 2^32 = 4,294,967,296.
 * <p>
 * Exponent: 0 = phantom/empty, 1-12 = relNlz 0-11. Standard tiers
 * (tier = NLZ, TIER_SCALE = 1.0).
 * <p>
 * History: 3 states collapsed from UDLL6's 4-state 2-bit model.
 *   State 0 (00): no sub-tier observations.
 *   State 1 (01): observed element 2 tiers below only.
 *   State 2 (1x): observed element 1 tier below (±2 below).
 * State 10 (rarest in UDLL6) collapsed with 11 into state 2.
 * <p>
 * At 1KB: 256 words × 6 = 1536 buckets (vs UDLL6's 1280). 20% more
 * buckets at equal memory; error scales as 1/sqrt(B).
 * <p>
 * Modular bucket addressing: bucket = key % modBuckets (not power-of-2).
 *
 * @author Neptune
 * @date April 2026
 */
public final class ArithmeticUltraDynamicLogLog32 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ArithmeticUltraDynamicLogLog32(){
		this(2048, 31, -1, 0);
	}

	ArithmeticUltraDynamicLogLog32(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ArithmeticUltraDynamicLogLog32(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*6<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*6-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*6;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=modBuckets/6;
		packed=new int[packedLen];
		minZeroCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+5)/6);}

	@Override
	public ArithmeticUltraDynamicLogLog32 copy(){return new ArithmeticUltraDynamicLogLog32(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read combined bucket value (0-38) for bucket i. */
	int readCombined(final int i){
		final int wordIdx=i/6;
		final int pos=i%6;
		final long word=Integer.toUnsignedLong(packed[wordIdx]);
		return (int)((word/POW39[pos])%RADIX);
	}

	/** Read exponent for bucket i (0-12). */
	int readExponent(final int i){return readCombined(i)/3;}

	/** Read history state for bucket i (0-2). */
	int readHistory(final int i){return readCombined(i)%3;}

	/** Write combined bucket value for bucket i, preserving other buckets in word. */
	private void writeCombined(final int i, final int combined){
		assert(combined>=0 && combined<RADIX);
		final int wordIdx=i/6;
		final int pos=i%6;
		final long word=Integer.toUnsignedLong(packed[wordIdx]);
		final int old=(int)((word/POW39[pos])%RADIX);
		final long diff=(long)(combined-old)*POW39[pos];
		packed[wordIdx]=(int)(word+diff);
	}

	/** Write exponent and history for bucket i. */
	void writeBucket(final int i, final int exp, final int hist){
		assert(exp>=0 && exp<=MAX_EXP);
		assert(hist>=0 && hist<NUM_HIST);
		writeCombined(i, exp*NUM_HIST+hist);
	}

	/*--------------------------------------------------------------*/
	/*----------------     History Transitions      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Compute new history state after tier advance by delta.
	 * Derived from UDLL6's carry-shift with 10+11 collapsed to 1x:
	 *   delta=1: any → 2 (observed at -1)
	 *   delta=2: any → 1 (observed at -2 only)
	 *   delta>=3: any → 0 (shifted out)
	 */
	private static int advanceHistory(final int oldHist, final int delta){
		assert(delta>=1);
		if(delta==1){return 2;}
		if(delta==2){return 1;}
		return 0;
	}

	/**
	 * Set history for sub-tier observation at given distance below current tier.
	 *   diff=1 (-1 tier): any → 2
	 *   diff=2 (-2 tiers): 0→1, 1→1, 2→2
	 */
	private static int observeHistory(final int oldHist, final int diff){
		assert(diff>=1 && diff<=HISTORY_MARGIN);
		if(diff==1){return 2;}
		// diff==2: set low bit (01) without clearing high bit (1x)
		return Math.max(oldHist, 1);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		if(packedBuckets==null){packedBuckets=new char[modBuckets];}
		else{java.util.Arrays.fill(packedBuckets, (char)0);}

		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int combined=readCombined(i);
			final int exp=combined/3;
			final int hist=combined%3;
			final int absNlz=exp+globalNLZ;
			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				// Emit hist for CardStats: hbits=2, but only states 0-2 used.
				// Mask invalid history at low tiers.
				final int emitHist=emitHistMask(absNlz, hist);
				packedBuckets[i]=(char)(((absNlz+1)<<2)|emitHist);
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, 2/*hbits*/, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.0);
	}

	/**
	 * Map internal 3-state history (0,1,2) to 2-bit pattern for CardStats.
	 * State 2 ("1x" = observed at -1, unknown -2) is emitted as binary 11
	 * at tier 2+ (both bits valid), but as binary 10 at tier 1 (only high
	 * bit valid — sbsStateIndex rejects states with invalid low bits).
	 *   absNlz==0: no valid history
	 *   absNlz==1: only bit 1 valid → state 2→10, others→00
	 *   absNlz>=2: state 0→00, state 1→01, state 2→11
	 */
	private static int emitHistMask(final int absNlz, final int hist){
		if(absNlz==0){return 0;}
		if(absNlz==1){return (hist==2) ? 2 : 0;} // only high bit valid at tier 1
		if(hist==2){return 3;} // collapsed 1x → emit as 11 at tier 2+
		return hist; // state 0→00, state 1→01
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
		add((ArithmeticUltraDynamicLogLog32)log);
	}

	public void add(ArithmeticUltraDynamicLogLog32 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			for(int i=0; i<modBuckets; i++){
				final int expA=adjustExp(readExponent(i), globalNLZ, newGlobalNLZ);
				final int expB=adjustExp(log.readExponent(i), log.globalNLZ, newGlobalNLZ);
				final int histA=readHistory(i);
				final int histB=log.readHistory(i);
				if(expA>expB){
					writeBucket(i, expA, histA);
				}else if(expB>expA){
					writeBucket(i, expB, histB);
				}else{
					// Same exponent: merge history (take max — higher state = more info)
					writeBucket(i, expA, Math.max(histA, histB));
				}
			}
			globalNLZ=newGlobalNLZ;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
				final int exp=readExponent(i);
				if(exp>0){filledBuckets++;}
				if(exp==0 || (!EARLY_PROMOTE && exp==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && globalNLZ<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	private static int adjustExp(int exp, int oldFloor, int newFloor){
		if(exp==0){return 0;}
		return Math.max(0, Math.min(exp+(oldFloor-newFloor), MAX_EXP));
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
		if(!DISABLE_EEMASK && Long.compareUnsigned(key, eeMask)>0){return;}

		final int rawNlz=Long.numberOfLeadingZeros(key);
		// Standard tiers: tier = rawNlz directly (no mantissa compression)
		final int absTier=rawNlz;
		final int relTier=absTier-(globalNLZ+1);
		if(relTier<-HISTORY_MARGIN){return;}

		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int oldExp=readExponent(bucket);
		final int oldHist=readHistory(bucket);

		if(relTier>=0){
			final int newExp=Math.min(relTier+1, MAX_EXP);
			if(newExp>oldExp){
				// Tier advance
				final int delta=(oldExp>0) ? (newExp-oldExp) : (relTier+1);
				final int newHist;
				if(oldExp>0 || globalNLZ>=0){
					newHist=advanceHistory(oldHist, delta);
				}else{
					newHist=0; // First ever insertion, no meaningful history
				}
				lastCardinality=-1;
				writeBucket(bucket, newExp, newHist);
				if(oldExp==0){filledBuckets++;}

				final int oldRelTier=(oldExp==0 ? 0 : oldExp-1);
				final boolean shouldDecrement=EARLY_PROMOTE ? oldExp==0
					: (relTier>oldRelTier && oldRelTier==0);
				if(shouldDecrement && --minZeroCount<1){
					while(minZeroCount==0 && globalNLZ<wordlen){
						advanceFloor();
						minZeroCount=countAndDecrement();
					}
				}
				return;
			}
			// Sub-tier observation: element below current exponent
			if(newExp<oldExp){
				final int diff=oldExp-newExp;
				if(diff>=1 && diff<=HISTORY_MARGIN){
					final int newHist=observeHistory(oldHist, diff);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldExp, newHist);
					}
				}
			}
		}else if(relTier==-1){
			// One tier below floor: set history diff=1 relative to phantom,
			// or diff relative to current exponent
			if(oldExp>0){
				final int diff=oldExp+1; // distance from observed to current
				if(diff>=1 && diff<=HISTORY_MARGIN){
					final int newHist=observeHistory(oldHist, diff);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldExp, newHist);
					}
				}
			}else{
				// Phantom bucket: diff=1 from phantom at globalNLZ
				final int newHist=observeHistory(oldHist, 1);
				if(newHist!=oldHist){
					lastCardinality=-1;
					writeBucket(bucket, oldExp, newHist);
				}
			}
		}else{
			// relTier == -2: two tiers below floor
			final int diff=(oldExp>0) ? (oldExp+2) : 2;
			if(diff>=1 && diff<=HISTORY_MARGIN){
				final int newHist=observeHistory(oldHist, diff);
				if(newHist!=oldHist){
					lastCardinality=-1;
					writeBucket(bucket, oldExp, newHist);
				}
			}
		}
	}

	/** Advance the global floor by one tier. Update eeMask. */
	private void advanceFloor(){
		globalNLZ++;
		final int exitThreshold=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
	}

	/** Decrement all exponents by 1 across all words. Returns new minZeroCount.
	 *  Decodes all 6 buckets per word, decrements, re-encodes. */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			long word=Integer.toUnsignedLong(packed[w]);
			if(word==0){newMinZeroCount+=6; continue;} // all empty, skip

			// Decode 6 buckets (LSB-first via repeated mod/div)
			final int[] vals=new int[6];
			long rem=word;
			for(int b=0; b<6; b++){
				vals[b]=(int)(rem%RADIX);
				rem/=RADIX;
			}

			// Decrement exponents, preserve history
			for(int b=0; b<6; b++){
				final int exp=vals[b]/3;
				final int hist=vals[b]%3;
				if(exp>=1){
					final int newExp=exp-1;
					vals[b]=newExp*3+hist;
					if(EARLY_PROMOTE){
						if(newExp==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(newExp<=1){newMinZeroCount++;}
					}
				}else{
					newMinZeroCount++; // already phantom
				}
			}

			// Re-encode (MSB-first via Horner's method)
			long newWord=0;
			for(int b=5; b>=0; b--){
				newWord=newWord*RADIX+vals[b];
			}
			packed[w]=(int)newWord;
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/modBuckets;}

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

	/** Returns the cached CardStats from the last rawEstimates() call, then clears it. */
	public CardStats consumeLastSummarized(){
		final CardStats s=lastSummarized;
		lastSummarized=null;
		return s;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int modBuckets;
	private final int[] packed;
	private final int packedLen;
	private int globalNLZ=-1;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return globalNLZ+1;}
	private int[] nlzCounts;
	private char[] packedBuckets;
	private CardStats lastSummarized;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Number of history states: 3 (collapsed from UDLL6's 4). */
	static final int NUM_HIST=3;
	/** Number of exponent states: 13 (0=phantom, 1-12=real tiers). */
	static final int MAX_EXP=12;
	/** Combined states per bucket: 3×13=39. */
	static final int RADIX=NUM_HIST*(MAX_EXP+1);
	/** Tiers below floor that can still set history bits. */
	static final int HISTORY_MARGIN=2;

	/** Powers of 39 for positional decode. */
	static final long[] POW39=new long[6];
	static{
		POW39[0]=1;
		for(int i=1; i<6; i++){POW39[i]=POW39[i-1]*RADIX;}
		assert(POW39[5]*RADIX<=0xFFFFFFFFL) : "39^6 must fit in unsigned 32-bit";
	}

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionAUDLL32.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionAUDLL32_LC2BitHistSBS.tsv.gz";
	public static final String HSB_FILE="?cardinalityCorrectionAUDLL32_LC2BitHist.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("AUDLL32: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal CFs: 8k DDLs, 1536 buckets, maxmult=8192, tmcf=1 (Apr 21 2026).
	 *  Measured on cluster with production SBS (1M iters), HSB, collapsed-as-11 encoding. */
	@Override public float terminalMeanCF(){return 0.722158f;}
	@Override public float terminalMeanPlusCF(){return 0.895883f;}

}
