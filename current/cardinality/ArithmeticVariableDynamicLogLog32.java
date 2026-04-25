package cardinality;

import shared.Tools;

/**
 * ArithmeticVariableDynamicLogLog32: 32-bit arithmetic encoding with
 * variable-depth history per tier. Low tiers get standard 2-bit history
 * (where HC matters most for LDLC), high tiers get no history (pure DLC).
 * <p>
 * Word layout (32-bit int, 5 buckets per word):
 *   RADIX combined states per bucket, packed via mixed-radix arithmetic:
 *   word = s0 + R*s1 + R^2*s2 + R^3*s3 + R^4*s4
 *   RADIX^5 must fit in unsigned 32-bit.
 * <p>
 * State mapping (for HIST_TIER_LIMIT=H, MAX_EXP=M):
 *   State 0:              exp=0 (phantom/empty), no history
 *   States 1..4*H:        exp 1..H, 2-bit history (4 states per tier)
 *   States 4*H+1..RADIX-1: exp H+1..M, no history (1 state per tier)
 *   RADIX = 4*H + (M-H) + 1
 * <p>
 * When HIST_TIER_LIMIT=MAX_EXP, all tiers have 2-bit history (no variable part).
 * <p>
 * History: standard 4-state 2-bit model (same as UDLL6).
 * Inspired by Noire's VCDLL4 (variable 1-bit/0-bit history in 4-bit nibbles).
 * Modular bucket addressing: bucket = key % modBuckets.
 *
 * @author Neptune
 * @date April 2026
 */
public final class ArithmeticVariableDynamicLogLog32 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ArithmeticVariableDynamicLogLog32(){
		this(2048, 31, -1, 0);
	}

	ArithmeticVariableDynamicLogLog32(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ArithmeticVariableDynamicLogLog32(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*BPW<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*BPW-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*BPW;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=modBuckets/BPW;
		packed=new int[packedLen];
		minZeroCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+BPW-1)/BPW);}

	@Override
	public ArithmeticVariableDynamicLogLog32 copy(){return new ArithmeticVariableDynamicLogLog32(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------      State Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	/** Extract exponent from combined state. */
	static int stateToExp(int state){
		if(state==0){return 0;}
		final int histBound=4*HIST_TIER_LIMIT;
		if(state<=histBound){return (state-1)/4+1;}
		return state-3*HIST_TIER_LIMIT;
	}

	/** Extract 2-bit history from combined state. Returns 0 for phantom and 0-bit tiers. */
	static int stateToHist(int state){
		final int histBound=4*HIST_TIER_LIMIT;
		if(state<=0 || state>histBound){return 0;}
		return (state-1)%4;
	}

	/** Combine exponent and history into a state.
	 *  For exp > HIST_TIER_LIMIT, history is ignored (0-bit tier). */
	static int toState(int exp, int hist){
		if(exp==0){return 0;}
		if(exp<=HIST_TIER_LIMIT){return (exp-1)*4+(hist&3)+1;}
		return 4*HIST_TIER_LIMIT+(exp-HIST_TIER_LIMIT);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read combined state for bucket i. */
	int readCombined(final int i){
		final int wordIdx=i/BPW;
		final int pos=i%BPW;
		final long word=Integer.toUnsignedLong(packed[wordIdx]);
		return (int)((word/POW[pos])%RADIX);
	}

	/** Read exponent for bucket i. */
	int readExponent(final int i){return stateToExp(readCombined(i));}

	/** Read history for bucket i (0-3 for 2-bit tiers, 0 for 0-bit tiers). */
	int readHistory(final int i){return stateToHist(readCombined(i));}

	/** Write combined state for bucket i, preserving other buckets in word.
	 *  Uses long intermediate to avoid int overflow. */
	private void writeCombined(final int i, final int combined){
		assert(combined>=0 && combined<RADIX);
		final int wordIdx=i/BPW;
		final int pos=i%BPW;
		final long word=Integer.toUnsignedLong(packed[wordIdx]);
		final int old=(int)((word/POW[pos])%RADIX);
		final long diff=(long)(combined-old)*POW[pos];
		packed[wordIdx]=(int)(word+diff);
	}

	/** Write exponent and history for bucket i. */
	void writeBucket(final int i, final int exp, final int hist){
		assert(exp>=0 && exp<=MAX_EXP);
		writeCombined(i, toState(exp, hist));
	}

	/*--------------------------------------------------------------*/
	/*----------------     History Transitions      ----------------*/
	/*--------------------------------------------------------------*/

	/** Compute new 2-bit history after tier advance by delta.
	 *  Standard 4-state transitions (same as UDLL6):
	 *    delta=1: shift old -1 to -2, set new -1. hist = ((hist&gt;&gt;1)|2)
	 *    delta=2: old history drops, set -2 (old tier at -2). hist = 1
	 *    delta&gt;=3: all shifted out. hist = 0 */
	private static int advanceHistory(final int oldHist, final int delta){
		assert(delta>=1);
		if(delta==1){return ((oldHist>>1)|2);}
		if(delta==2){return 1;}
		return 0;
	}

	/** Set history for sub-tier observation at given distance below current tier.
	 *  diff=1 (-1 tier): set bit 1. hist |= 2
	 *  diff=2 (-2 tiers): set bit 0. hist |= 1 */
	private static int observeHistory(final int oldHist, final int diff){
		assert(diff>=1 && diff<=HISTORY_MARGIN);
		if(diff==1){return oldHist|2;}
		return oldHist|1;
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
			final int state=readCombined(i);
			final int exp=stateToExp(state);
			final int absNlz=exp+globalNLZ;
			if(absNlz<0 || absNlz>=64){continue;}
			filledCount++;
			nlzCounts[absNlz+1]++;
			if(exp<=HIST_TIER_LIMIT){
				// 2-bit history tiers: emit with actual history.
				final int hist=stateToHist(state);
				packedBuckets[i]=(char)(((absNlz+1)<<2)|hist);
			}else{
				// 0-bit tiers: emit with hist=0 (no history known).
				packedBuckets[i]=(char)((absNlz+1)<<2);
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, 2/*hbits*/, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.0);
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
		add((ArithmeticVariableDynamicLogLog32)log);
	}

	public void add(ArithmeticVariableDynamicLogLog32 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			for(int i=0; i<modBuckets; i++){
				final int expA=adjustExp(readExponent(i), globalNLZ, newGlobalNLZ);
				final int expB=adjustExp(log.readExponent(i), log.globalNLZ, newGlobalNLZ);
				int histA=readHistory(i);
				int histB=log.readHistory(i);
				if(expA>expB){
					writeBucket(i, expA, histA);
				}else if(expB>expA){
					writeBucket(i, expB, histB);
				}else{
					if(expA<=HIST_TIER_LIMIT){
						writeBucket(i, expA, histA|histB);
					}else{
						writeBucket(i, expA, 0);
					}
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
					newHist=0;
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
				if(oldExp>HIST_TIER_LIMIT){return;}
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
			if(oldExp>HIST_TIER_LIMIT){return;}
			if(oldExp>0){
				final int diff=oldExp;
				if(diff>=1 && diff<=HISTORY_MARGIN){
					final int newHist=observeHistory(oldHist, diff);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldExp, newHist);
					}
				}
			}else{
				final int newHist=observeHistory(oldHist, 1);
				if(newHist!=oldHist){
					lastCardinality=-1;
					writeBucket(bucket, oldExp, newHist);
				}
			}
		}else{
			// relTier == -2
			if(oldExp>HIST_TIER_LIMIT){return;}
			final int diff=(oldExp>0) ? (oldExp+1) : 1;
			if(diff>=1 && diff<=HISTORY_MARGIN){
				final int newHist=observeHistory(oldHist, diff);
				if(newHist!=oldHist){
					lastCardinality=-1;
					writeBucket(bucket, oldExp, newHist);
				}
			}
		}
	}

	private void advanceFloor(){
		globalNLZ++;
		final int exitThreshold=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
	}

	/** Decrement all exponents by 1 across all words. Returns new minZeroCount.
	 *  Full decode/re-encode path required.
	 *  Handles demotion at the 0-bit to 2-bit boundary. */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			long word=Integer.toUnsignedLong(packed[w]);
			if(word==0){newMinZeroCount+=BPW; continue;}

			final int[] states=new int[BPW];
			long rem=word;
			for(int b=0; b<BPW; b++){
				states[b]=(int)(rem%RADIX);
				rem/=RADIX;
			}

			for(int b=0; b<BPW; b++){
				final int exp=stateToExp(states[b]);
				final int hist=stateToHist(states[b]);
				if(exp>=1){
					final int newExp=exp-1;
					final int newHist;
					if(exp==HIST_TIER_LIMIT+1){
						// Demotion: 0-bit → 2-bit boundary
						if(DEMOTION_MODE==0){newHist=2;}
						else if(DEMOTION_MODE==1){newHist=0;}
						else{newHist=(w*BPW+b)&1;}
					}else{
						newHist=hist;
					}
					states[b]=toState(newExp, newHist);
					if(EARLY_PROMOTE){
						if(newExp==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(newExp<=1){newMinZeroCount++;}
					}
				}else{
					newMinZeroCount++;
				}
			}

			long newWord=0;
			for(int b=BPW-1; b>=0; b--){
				newWord=newWord*RADIX+states[b];
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
	/** Buckets per 32-bit word. */
	static final int BPW=5;
	/** Maximum exponent. */
	static int MAX_EXP=20;
	/** Last tier with 2-bit history. Set via histtiers= flag. */
	public static int HIST_TIER_LIMIT=20;
	/** Combined states per bucket: 1 + 4*HIST_TIER_LIMIT + (MAX_EXP-HIST_TIER_LIMIT). */
	static int RADIX=computeRadix();
	/** Tiers below floor that can still set history bits. */
	static final int HISTORY_MARGIN=2;

	/** Powers of RADIX for positional decode. */
	static long[] POW=computePow();

	static int computeRadix(){
		return 4*HIST_TIER_LIMIT+(MAX_EXP-HIST_TIER_LIMIT)+1;
	}

	static long[] computePow(){
		final int r=computeRadix();
		long[] pow=new long[BPW];
		pow[0]=1;
		for(int i=1; i<BPW; i++){pow[i]=pow[i-1]*r;}
		assert(pow[BPW-1]*r<=0xFFFFFFFFL) : r+"^"+BPW+" must fit in unsigned 32-bit";
		return pow;
	}

	/** Reconfigure radix after changing HIST_TIER_LIMIT or MAX_EXP. */
	public static void reconfigure(){
		RADIX=computeRadix();
		POW=computePow();
	}

	/** Demotion mode for 0-bit → 2-bit transition in countAndDecrement.
	 *  0=A (hist=2, assume -1 seen), 1=B (hist=0), 2=C (hist=bucketID&amp;1). */
	public static int DEMOTION_MODE=0;

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionAVDLL32.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionAVDLL32_LC2BitHistSBS.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("AVDLL32: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal CFs. Placeholder — recalibrate after config changes. */
	@Override public float terminalMeanCF(){return 0.720985f;}
	@Override public float terminalMeanPlusCF(){return 1.000266f;}

	@Override public double ldlcHcWeight(){return 0.30;}
	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.50f;}

}
