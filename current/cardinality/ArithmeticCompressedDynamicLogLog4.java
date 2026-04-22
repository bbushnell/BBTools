package cardinality;

import shared.Tools;

/**
 * ArithmeticCompressedDynamicLogLog4: arithmetic-encoded exponents with
 * 1-bit history, extending CDLL4's 8 exponent states to 10 via 3 extra
 * bits per word stored in a separate long[] array.
 * <p>
 * Word layout (32-bit int, 8 buckets per word):
 *   - Bits [31:24]: 8 × 1-bit history (fixed-width)
 *   - Bits [23:0]:  lower 24 bits of radix-10 exponent block
 * <p>
 * Extra bits (long[], 3 bits per word packed 21 per long):
 *   - Upper 3 bits of radix-10 exponent block per word
 *   - Combined: (extra3 << 24) | lower24 = 27-bit exponent block
 *   - 10^8 = 100,000,000 < 2^27 = 134,217,728
 * <p>
 * Each exponent: 0 = empty/phantom, 1-9 = relTier 0-8.
 * 10 exponent states (vs CDLL4's 8) from 27-bit arithmetic packing.
 * <p>
 * Tier geometry: same as CDLL4 (halfNlz/3, TIER_SCALE=1.5).
 * 1-bit history: same model as CDLL4 (set when sub-tier position
 * reaches top of tier).
 * <p>
 * Memory: 4.375 bits/bucket (35 bits per 8 buckets).
 * At 1KB: 232 words × 8 = 1856 buckets.
 *
 * @author Nahida
 * @date April 2026
 */
public final class ArithmeticCompressedDynamicLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ArithmeticCompressedDynamicLogLog4(){
		this(2048, 31, -1, 0);
	}

	ArithmeticCompressedDynamicLogLog4(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ArithmeticCompressedDynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*8<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*8-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*8;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=modBuckets/8;
		packed=new int[packedLen];
		extraBits=new long[(packedLen+20)/21];
		minZeroCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+7)/8);}

	@Override
	public ArithmeticCompressedDynamicLogLog4 copy(){return new ArithmeticCompressedDynamicLogLog4(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read exponent for bucket i (0..modBuckets-1). Returns 0-9. */
	int readExponent(final int i){
		final int wordIdx=i>>>3;
		final int pos=i&7;
		final int lower24=packed[wordIdx]&EXP_MASK;
		final int extra=readExtra(wordIdx);
		final int expBlock=lower24|(extra<<24);
		return (expBlock/POW10[pos])%10;
	}

	/** Read 1-bit history for bucket i. */
	int readHistory(final int i){
		final int wordIdx=i>>>3;
		final int pos=i&7;
		return (packed[wordIdx]>>>(24+pos))&1;
	}

	/** Write exponent for bucket i (0-9), preserving history and other exponents. */
	void writeExponent(final int i, final int newExp){
		final int wordIdx=i>>>3;
		final int pos=i&7;
		final int lower24=packed[wordIdx]&EXP_MASK;
		final int extra=readExtra(wordIdx);
		int expBlock=lower24|(extra<<24);
		final int oldExp=(expBlock/POW10[pos])%10;
		expBlock+=(newExp-oldExp)*POW10[pos];
		packed[wordIdx]=(packed[wordIdx]&HIST_WORD_MASK)|(expBlock&EXP_MASK);
		writeExtra(wordIdx, (expBlock>>>24)&0x7);
	}

	/** Write 1-bit history for bucket i, preserving exponents and other history. */
	void writeHistory(final int i, final int hist){
		final int wordIdx=i>>>3;
		final int pos=i&7;
		final int shift=24+pos;
		packed[wordIdx]=(packed[wordIdx]&~(1<<shift))|((hist&1)<<shift);
	}

	/** Write both exponent and history for bucket i. */
	void writeBucket(final int i, final int exp, final int hist){
		writeExponent(i, exp);
		writeHistory(i, hist);
	}

	/** Read 3 extra exponent bits for word w from long[] array. */
	private int readExtra(final int w){
		final int longIdx=w/21;
		final int shift=(w%21)*3;
		return (int)((extraBits[longIdx]>>>shift)&0x7);
	}

	/** Write 3 extra exponent bits for word w into long[] array. */
	private void writeExtra(final int w, final int val){
		final int longIdx=w/21;
		final int shift=(w%21)*3;
		extraBits[longIdx]=(extraBits[longIdx]&~(0x7L<<shift))|((long)(val&0x7)<<shift);
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
			final int exp=readExponent(i);
			final int hist=readHistory(i);
			final int absTier=exp+globalNLZ;
			if(absTier<0){continue;}
			filledCount++;
			if(absTier<64){
				nlzCounts[absTier+1]++;
				// Mask hist at absTier==0: no valid history below tier 0
				final int emitHist=(absTier==0) ? 0 : hist;
				packedBuckets[i]=(char)(((absTier+1)<<1)|emitHist);
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, 1, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.5);
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
		add((ArithmeticCompressedDynamicLogLog4)log);
	}

	public void add(ArithmeticCompressedDynamicLogLog4 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			for(int i=0; i<modBuckets; i++){
				int expA=adjustExp(readExponent(i), globalNLZ, newGlobalNLZ);
				int expB=adjustExp(log.readExponent(i), log.globalNLZ, newGlobalNLZ);
				int histA=readHistory(i);
				int histB=log.readHistory(i);
				if(expA>expB){
					writeBucket(i, expA, histA);
				}else if(expB>expA){
					writeBucket(i, expB, histB);
				}else{
					writeBucket(i, expA, histA|histB);
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
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF);
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		final int mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		final int halfNlz=2*rawNlz+mantissa;
		final int absTier=halfNlz/3;
		final int relTier=absTier-(globalNLZ+1);
		if(relTier<-1){return;}

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
				final int newHist=(delta==1) ? 1 : 0;
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
			// delta == -1: element from previous tier, set history
			if(newExp==oldExp-1){
				if(oldHist==0){lastCardinality=-1; writeHistory(bucket, 1);}
			}
		}else{
			// relTier == -1: element from one tier below floor
			if(oldExp<=1){
				if(oldHist==0){lastCardinality=-1; writeHistory(bucket, 1);}
			}
		}
	}

	private void advanceFloor(){
		globalNLZ++;
		final int relaxedTier=Math.max(0, globalNLZ);
		final int minNlz=(3*relaxedTier)/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	/** Decrement all exponents by 1 across all words. Returns new minZeroCount. */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			final int lower24=packed[w]&EXP_MASK;
			final int extra=readExtra(w);
			int expBlock=lower24|(extra<<24);
			final int histBits=packed[w]&HIST_WORD_MASK;
			if(expBlock==0){continue;}

			// Decode all 8 exponents, decrement, re-encode
			int newExpBlock=0;
			int mult=1;
			int remaining=expBlock;
			for(int b=0; b<8; b++){
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
			packed[w]=histBits|(newExpBlock&EXP_MASK);
			writeExtra(w, (newExpBlock>>>24)&0x7);
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

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int modBuckets;
	private final int[] packed;
	private final long[] extraBits;
	private final int packedLen;
	private int globalNLZ=-1;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return globalNLZ+1;}
	private int[] nlzCounts;
	private char[] packedBuckets;
	private CardStats lastSummarized;

	/** Return and clear the CardStats from the most recent summarize/rawEstimates call. */
	public CardStats consumeLastSummarized(){
		final CardStats cs=lastSummarized;
		lastSummarized=null;
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	private static final int MAX_EXP=9; // 0=phantom, 1-9=real tiers

	/** Radix-10 powers for 8-bucket access. */
	private static final int[] POW10={1, 10, 100, 1000, 10000, 100000, 1000000, 10000000};

	/** Mask for lower 24 bits (exponent block in packed word). */
	private static final int EXP_MASK=0xFFFFFF;
	/** Mask for upper 8 bits (history block in packed word). */
	private static final int HIST_WORD_MASK=0xFF000000;

	/** Mantissa threshold for log-uniform half-NLZ steps: (2-sqrt(2)) * 1048576 */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionACDLL4.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionACDLL4_LC1BitHistSBS.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal Mean CF — measured via ddlcalibrate2 cf=f tmcf=1 tmpcf=1 ddls=128k maxmult=8192, with SBS. */
	@Override public float terminalMeanCF(){return 0.878866f;}

	/** Terminal Mean+H CF — measured via ddlcalibrate2 cf=f tmcf=1 tmpcf=1 ddls=128k maxmult=8192, with SBS. */
	@Override public float terminalMeanPlusCF(){return 1.061534f;}

	/** HLDLC weight — placeholder until tuned. */
	@Override public float hldlcWeight(){return 0.325f;}

}
