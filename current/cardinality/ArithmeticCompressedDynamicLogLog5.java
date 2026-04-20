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
 * Each exponent: 0 = phantom, 1-9 = relTier 0-8.
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

	ArithmeticCompressedDynamicLogLog5(){
		this(2048, 31, -1, 0);
	}

	ArithmeticCompressedDynamicLogLog5(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

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

	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		if(packedBuckets==null){packedBuckets=new char[modBuckets];}
		else{java.util.Arrays.fill(packedBuckets, (char)0);}

		final int phantomTier=minZeros-1;
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int exp=readExponent(i);
			final int hist=readHistory(i);
			if(exp>0){
				final int absTier=(exp-1)+minZeros;
				if(absTier<64){
					nlzCounts[absTier+1]++;
					final int emitHist=emitHistMask(absTier, hist);
					packedBuckets[i]=(char)(((absTier+1)<<2)|emitHist);
				}
				filledCount++;
			}else if(minZeros>0 && phantomTier>=0 && phantomTier<64){
				nlzCounts[phantomTier+1]++;
				final int emitHist=emitHistMask(phantomTier, hist);
				packedBuckets[i]=(char)(((phantomTier+1)<<2)|emitHist);
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.5);
	}

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

	public void add(ArithmeticCompressedDynamicLogLog5 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<modBuckets; i++){
				int expA=adjustExp(readExponent(i), readHistory(i), minZeros, newMinZeros);
				int expB=adjustExp(log.readExponent(i), log.readHistory(i), log.minZeros, newMinZeros);
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
			minZeros=newMinZeros;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
				final int exp=readExponent(i);
				if(exp>0){filledBuckets++;}
				if(exp==0 || (!EARLY_PROMOTE && exp==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && minZeros<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	private static int adjustExp(int exp, int hist, int oldFloor, int newFloor){
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
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		final int mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		final int halfNlz=2*rawNlz+mantissa;
		final int absTier=halfNlz/3;
		final int relTier=absTier-minZeros;
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
				final int carry=(oldExp>0 || minZeros>0) ? HIST_CARRY : 0;
				final int newHist=((oldHist|carry)>>delta)&HIST_MASK;
				lastCardinality=-1;
				writeBucket(bucket, newExp, newHist);
				if(oldExp==0){filledBuckets++;}

				final int oldRelTier=(oldExp==0 ? 0 : oldExp-1);
				final boolean shouldDecrement=EARLY_PROMOTE ? oldExp==0
					: (relTier>oldRelTier && oldRelTier==0);
				if(shouldDecrement && --minZeroCount<1){
					while(minZeroCount==0 && minZeros<wordlen){
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

	private void advanceFloor(){
		minZeros++;
		final int relaxedTier=Math.max(0, minZeros-HISTORY_MARGIN);
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
	private final int packedLen;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return minZeros;}
	private int[] nlzCounts;
	private char[] packedBuckets;
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
	private static final int HBITS=2;
	private static final int HIST_MASK=(1<<HBITS)-1;
	private static final int HIST_CARRY=1<<HBITS;
	static final int HISTORY_MARGIN=2;
	private static final int MAX_EXP=9; // 0=phantom, 1-9=real tiers

	/** Radix-10 powers for single-bucket access. */
	private static final int[] POW10={1, 10, 100, 1000, 10000, 100000};

	/** Mask for exponent block (bottom 20 bits). */
	private static final int EXP_MASK=0xFFFFF;
	/** Mask for history block (top 12 bits). */
	private static final int HIST_WORD_MASK=0xFFF00000;

	/** Mantissa threshold for log-uniform half-NLZ steps: (2-sqrt(2)) * 65536 */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionACDLL5.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionACDLL5_LC2BitHistSBS.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("ACDLL5: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal CFs from cluster calibration: 512k DDLs, 1536 buckets, tmcf=1. */
	@Override public float terminalMeanCF(){return 0.878596f;}
	@Override public float terminalMeanPlusCF(){return 1.088939f;}

}
