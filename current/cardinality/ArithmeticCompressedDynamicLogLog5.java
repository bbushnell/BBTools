package cardinality;

import shared.Tools;

/**
 * ArithmeticCompressedDynamicLogLog5: arithmetic-encoded exponents for
 * fractional bits per bucket.
 * <p>
 * Word layout (32-bit int, 6 buckets per word):
 *   - Bits [31:20]: 6 × 2-bit history (12 bits, fixed-width)
 *   - Bits [19:0]:  6 exponents packed radix-10 (10^6 = 1,000,000 < 2^20)
 * <p>
 * Each exponent: 0 = floor-level, 1-9 = relTier 0-8.
 * 10 exponent states (vs CDLL5's 8) from 20-bit arithmetic packing.
 * Zero wasted bits per word (vs CDLL5's 2 wasted).
 * <p>
 * Tier geometry: same as CDLL5 (halfNlz/3, TIER_SCALE=1.5) extended to 9 real tiers.
 * 2-bit history: same carry-shift model as CDLL5/UDLL6.
 *
 * @author UMP45
 * @date April 2026
 */
public final class ArithmeticCompressedDynamicLogLog5 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	ArithmeticCompressedDynamicLogLog5(){
		this(2048, 31, -1, 0);
	}

	/** Construct from parsed command-line arguments. */
	ArithmeticCompressedDynamicLogLog5(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	/**
	 * Full constructor.
	 * Bucket count is rounded up to the next multiple of 6 (complete words).
	 * @param buckets_ Number of buckets (rounded to next multiple of 6)
	 * @param k_ Hash prefix length
	 * @param seed Random seed (-1 for default)
	 * @param minProb_ Minimum probability threshold
	 */
	ArithmeticCompressedDynamicLogLog5(int buckets_, int k_, long seed, float minProb_){
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

	/** Create an independent copy with a fresh seed. */
	@Override
	public ArithmeticCompressedDynamicLogLog5 copy(){return new ArithmeticCompressedDynamicLogLog5(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read exponent for bucket i (0..buckets-1). Returns 0-9. */
	int readExponent(final int i){
		final int wordIdx=i/6;
		final int pos=i%6;
		int expBlock=packed[wordIdx]&EXP_MASK;
		return (expBlock/POW10[pos])%10;
	}

	/** Read 2-bit history for bucket i. */
	int readHistory(final int i){
		final int wordIdx=i/6;
		final int pos=i%6;
		return (packed[wordIdx]>>>(20+pos*2))&0x3;
	}

	/** Write exponent for bucket i (0-9), preserving history and other exponents. */
	void writeExponent(final int i, final int newExp){
		final int wordIdx=i/6;
		final int pos=i%6;
		int word=packed[wordIdx];
		int expBlock=word&EXP_MASK;
		int oldExp=(expBlock/POW10[pos])%10;
		expBlock+=(newExp-oldExp)*POW10[pos];
		packed[wordIdx]=(word&HIST_WORD_MASK)|(expBlock&EXP_MASK);
	}

	/** Write 2-bit history for bucket i, preserving exponents and other history. */
	void writeHistory(final int i, final int hist){
		final int wordIdx=i/6;
		final int pos=i%6;
		final int shift=20+pos*2;
		int word=packed[wordIdx];
		word=(word&~(0x3<<shift))|((hist&0x3)<<shift);
		packed[wordIdx]=word;
	}

	/** Write both exponent and history for bucket i. */
	void writeBucket(final int i, final int exp, final int hist){
		writeExponent(i, exp);
		writeHistory(i, hist);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Build NLZ histogram and per-bucket packed state for CardStats.
	 * Converts each bucket's (exp, hist) to absolute tier,
	 * emits into nlzCounts[absTier+1] and packedBuckets[i].
	 * @return CardStats with all estimator values computed
	 */
	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		if(packedBuckets==null){packedBuckets=new char[modBuckets];}
		else{java.util.Arrays.fill(packedBuckets, (char)0);}

		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int exp=readExponent(i);
			final int hist=readHistory(i);
			final int absTier=exp+globalNLZ;
			if(absTier>=0 && absTier<64){
				nlzCounts[absTier+1]++;
				final int emitHist=emitHistMask(absTier, hist);
				packedBuckets[i]=(char)(((absTier+1)<<2)|emitHist);
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.5);
	}

	/**
	 * Mask history bits to those valid at the given absolute tier.
	 * Tier 0: no sub-tiers exist, all hist invalid.
	 * Tier 1: only bit 1 (one tier below) is valid.
	 * Tier 2+: both bits valid.
	 */
	private static int emitHistMask(int absTier, int hist){
		if(absTier==0){return 0;}
		if(absTier==1){return hist&0x2;}
		return hist;
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
		add((ArithmeticCompressedDynamicLogLog5)log);
	}

	/**
	 * Merge another ACDLL5 into this one.
	 * Each bucket's exponent is reframed to the new common globalNLZ; ties merge history with OR.
	 */
	public void add(ArithmeticCompressedDynamicLogLog5 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newMinZeros=Math.max(globalNLZ, log.globalNLZ);
			for(int i=0; i<modBuckets; i++){
				int expA=adjustExp(readExponent(i), readHistory(i), globalNLZ, newMinZeros);
				int expB=adjustExp(log.readExponent(i), log.readHistory(i), log.globalNLZ, newMinZeros);
				int histA=readHistory(i);
				int histB=log.readHistory(i);
				// Adjust history during floor change (approximate: just keep)
				if(expA>expB){
					writeBucket(i, expA, histA);
				}else if(expB>expA){
					writeBucket(i, expB, histB);
				}else{
					writeBucket(i, expA, histA|histB);
				}
			}
			globalNLZ=newMinZeros;
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

	/** Reframe exponent from oldFloor to newFloor; clamped to [0, MAX_EXP]. */
	private static int adjustExp(int exp, int hist, int oldFloor, int newFloor){
		if(exp==0){return 0;}
		return Math.max(0, Math.min(exp+(oldFloor-newFloor), MAX_EXP));
	}

	/**
	 * Hash a value and store it in the appropriate bucket.
	 * Pipeline: hash → NLZ → mantissa → halfNlz → compressed tier →
	 * relTier (minus globalNLZ+1) → store or update history.
	 */
	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
		if(!DISABLE_EEMASK && Long.compareUnsigned(key, eeMask)>0){return;} // early-exit mask

		final int rawNlz=Long.numberOfLeadingZeros(key);
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		final int mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		final int halfNlz=2*rawNlz+mantissa;
		final int absTier=halfNlz/3;
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
				final int delta=(oldExp>0) ? (newExp-oldExp) : (relTier+1);
				final int carry=(oldExp>0 || globalNLZ>=0) ? HIST_CARRY : 0;
				final int newHist=((oldHist|carry)>>delta)&HIST_MASK;
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
			if(newExp<oldExp){
				final int diff=oldExp-newExp;
				if(diff>=1 && diff<=HBITS){
					final int bit=HBITS-diff;
					final int newHist=oldHist|(1<<bit);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeHistory(bucket, newHist);
					}
				}
			}
		}else if(relTier==-1){
			if(oldExp>0){
				final int diff=oldExp;
				if(diff>=1 && diff<=HBITS){
					final int bit=HBITS-diff;
					final int newHist=oldHist|(1<<bit);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeHistory(bucket, newHist);
					}
				}
			}
		}else{
			final int diff=(oldExp>0) ? (oldExp+1) : 1;
			if(diff>=1 && diff<=HBITS){
				final int bit=HBITS-diff;
				final int newHist=oldHist|(1<<bit);
				if(newHist!=oldHist){
					lastCardinality=-1;
					writeHistory(bucket, newHist);
				}
			}
		}
	}

	/**
	 * Advance global floor by one tier.
	 * Relaxes eeMask by HISTORY_MARGIN tiers so hashes at relTier in
	 * {-HISTORY_MARGIN, ..., -1} pass through for history updates.
	 */
	private void advanceFloor(){
		globalNLZ++;
		final int relaxedTier=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		final int minNlz=(3*relaxedTier)/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	/** Decrement all exponents by 1 across all words. Returns new minZeroCount. */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			int word=packed[w];
			int expBlock=word&EXP_MASK;
			int histBlock=word&HIST_WORD_MASK;
			if(expBlock==0){continue;} // all 6 exponents are 0, skip

			// Decode all 6 exponents, decrement, re-encode
			int newExpBlock=0;
			int mult=1;
			int remaining=expBlock;
			for(int b=0; b<6; b++){
				int exp=remaining%10;
				remaining/=10;
				if(exp>=1){
					exp--;
					if(EARLY_PROMOTE){
						if(exp==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(exp<=1){newMinZeroCount++;}
					}
				}
				newExpBlock+=exp*mult;
				mult*=10;
			}
			packed[w]=histBlock|(newExpBlock&EXP_MASK);
		}
		return newMinZeroCount;
	}

	/** Number of buckets with exponent > 0. */
	public int filledBuckets(){return filledBuckets;}

	/** Fraction of buckets that are occupied (exponent > 0). */
	public double occupancy(){return (double)filledBuckets/modBuckets;}

	/** Not used; CF correction handled via CF_MATRIX. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Compute all estimator values and return as a legacy-format array.
	 * Caches the CardStats for retrieval by consumeLastSummarized().
	 */
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

	/** Packed storage: 6 exponents (radix-10, bits [19:0]) + 6 × 2-bit history (bits [31:20]) per int. */
	private final int[] packed;
	/** Number of packed ints: ceil(modBuckets/6). */
	private final int packedLen;
	/** Actual bucket count (multiple of 6, may differ from super.buckets). */
	private final int modBuckets;
	/** Global floor tier: -1 means nothing seen; >=0 means all buckets have absNlz >= globalNLZ. absNlz = exp + globalNLZ. */
	private int globalNLZ=-1;
	/** Buckets at the global floor eligible for next advance. */
	private int minZeroCount;
	/** Count of buckets with exponent > 0. */
	private int filledBuckets=0;
	/** Early-exit mask: filters hashes below HISTORY_MARGIN tiers of floor. */
	private long eeMask=-1L;

	/** Compatibility accessor: returns globalNLZ+1 to match legacy minZeros convention. */
	public int getMinZeros(){return globalNLZ+1;}
	/** Reusable NLZ histogram buffer (index 0=empty, 1-65=absNlz 0-64). */
	private int[] nlzCounts;
	/** Reusable per-bucket packed state for CardStats: (absTier+1)<<2 | hist. */
	private char[] packedBuckets;
	/** Cached CardStats from last rawEstimates() call. */
	private CardStats lastSummarized;

	/** Return and clear the cached CardStats. Used by calibration drivers. */
	public CardStats consumeLastSummarized(){
		final CardStats cs=lastSummarized;
		lastSummarized=null;
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum NLZ value (64-bit hash). */
	private static final int wordlen=64;
	/** Number of history bits per bucket. */
	private static final int HBITS=2;
	/** Bitmask for history: (1<<HBITS)-1 = 0x3. */
	private static final int HIST_MASK=(1<<HBITS)-1;
	/** Carry bit injected on tier advance: 1<<HBITS = 4. */
	private static final int HIST_CARRY=1<<HBITS;
	/** eeMask relaxation: accept hashes this many tiers below floor for history. */
	static final int HISTORY_MARGIN=2;
	/** Maximum exponent value: 0=floor-level, 1-9=real tiers. */
	private static final int MAX_EXP=9;

	/** Radix-10 powers for single-bucket access: POW10[pos] = 10^pos. */
	private static final int[] POW10={1, 10, 100, 1000, 10000, 100000};

	/** Mask for exponent block (bottom 20 bits): holds 10^6 = 1,000,000 combinations. */
	private static final int EXP_MASK=0xFFFFF;
	/** Mask for history block (top 12 bits): 6 × 2-bit history fields. */
	private static final int HIST_WORD_MASK=0xFFF00000;

	/** Mantissa threshold for compressed tiers: (2-sqrt(2)) * 1048576 ≈ 614242. */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	/** If true, stored=0 triggers floor advance (vs stored<=1 when false). */
	public static boolean EARLY_PROMOTE=true;
	/** Debug flag: bypass eeMask filtering. Not for production. */
	public static boolean DISABLE_EEMASK=false;

	/** Auto-loaded v5 CF table resource path. */
	public static final String CF_FILE="?cardinalityCorrectionACDLL5.tsv.gz";
	/** SBS (per-fill LC correction) table for arithmetic-encoded 2-bit history. */
	public static final String SBS_FILE="?cardinalityCorrectionACDLL5_LC2BitHistSBS.tsv.gz";
	/** Bucket count the CF_MATRIX was generated for. */
	private static int CF_BUCKETS=2048;
	/** Per-cardinality correction factor table, set by initializeCF or setCFMatrix. */
	private static float[][] CF_MATRIX=null;

	/** Load the CF table from the resource file, scaled to the given bucket count. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("ACDLL5: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}

	/** Set the CF matrix directly (used by CardinalityParser after global load). */
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal CFs from cluster calibration: 512k DDLs, 1536 buckets, tmcf=1. */
	@Override public float terminalMeanCF(){return 0.878596f;}
	@Override public float terminalMeanPlusCF(){return 1.088939f;}

	/** HC weight for LDLC blend. Calibrated by Eru, 2026-04-22. */
	@Override public double ldlcHcWeight(){return 0.40;}
	/** HLDLC weight. Calibrated by UMP45, 2026-04-22: 8k DDLs, 1536 buckets.
	 *  Dead flat 0.70-0.78; Hybrid+2 contributes ~26%. */
	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.74f;}

}
