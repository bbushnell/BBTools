package cardinality;

import shared.Tools;

/**
 * ArithmeticUltraDynamicLogLog33: 33-bit arithmetic encoding with
 * standard-tier cardinality estimator, 3-state history, 14 real tiers.
 * <p>
 * Extends AUDLL32's design by adding 1 extra bit per word in a parallel
 * int array, widening the packed value to 33 bits. This raises the radix
 * from 39 (13 exp × 3 hist) to 45 (15 exp × 3 hist), gaining 2 tiers
 * for minimal memory overhead.
 * <p>
 * Word layout (33-bit logical word: 32-bit int + 1 bit from extraBits[]):
 *   6 buckets per word, radix 45: word = b0 + 45*b1 + ... + 45^5*b5
 *   45^6 = 8,303,765,625 < 2^33 = 8,589,934,592.
 * <p>
 * Exponent: 0 = phantom/empty, 1-14 = relNlz 0-13. Standard tiers.
 * History: 3 collapsed states (same as AUDLL32).
 * <p>
 * At 1KB: 248 main words + 8 extra-bit ints = 1024 bytes.
 *   248 × 6 = 1488 buckets with 14 real tiers.
 *   (vs UDLL6's 1280 buckets with 15 tiers at 1KB.)
 * <p>
 * Modular bucket addressing: bucket = key % modBuckets.
 *
 * @author Neptune
 * @date April 2026
 */
public final class ArithmeticUltraDynamicLogLog33 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ArithmeticUltraDynamicLogLog33(){
		this(2048, 31, -1, 0);
	}

	ArithmeticUltraDynamicLogLog33(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ArithmeticUltraDynamicLogLog33(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*6<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*6-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*6;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=modBuckets/6;
		packed=new int[packedLen];
		extraBitsLen=(packedLen+31)>>>5;
		extraBits=new int[extraBitsLen];
		minZeroCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+5)/6);}

	@Override
	public ArithmeticUltraDynamicLogLog33 copy(){return new ArithmeticUltraDynamicLogLog33(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read the full 33-bit packed word for the given word index. */
	private long readWord33(final int wordIdx){
		final long extraBit=((long)(extraBits[wordIdx>>>5]>>>(wordIdx&31)))&1L;
		return (extraBit<<32)|Integer.toUnsignedLong(packed[wordIdx]);
	}

	/** Write a 33-bit packed word: low 32 bits to packed[], bit 32 to extraBits[]. */
	private void writeWord33(final int wordIdx, final long word33){
		packed[wordIdx]=(int)word33;
		final int ebWord=wordIdx>>>5;
		final int ebBit=wordIdx&31;
		final int newEB=(int)((word33>>>32)&1);
		extraBits[ebWord]=(extraBits[ebWord]&~(1<<ebBit))|(newEB<<ebBit);
	}

	/** Read combined bucket value (0-44) for bucket i. */
	int readCombined(final int i){
		final int wordIdx=i/6;
		final int pos=i%6;
		final long word=readWord33(wordIdx);
		return (int)((word/POW45[pos])%RADIX);
	}

	/** Read exponent for bucket i (0-14). */
	int readExponent(final int i){return readCombined(i)/NUM_HIST;}

	/** Read history state for bucket i (0-2). */
	int readHistory(final int i){return readCombined(i)%NUM_HIST;}

	/** Write combined bucket value for bucket i, preserving other buckets in word. */
	private void writeCombined(final int i, final int combined){
		assert(combined>=0 && combined<RADIX);
		final int wordIdx=i/6;
		final int pos=i%6;
		final long word=readWord33(wordIdx);
		final int old=(int)((word/POW45[pos])%RADIX);
		final long diff=(long)(combined-old)*POW45[pos];
		writeWord33(wordIdx, word+diff);
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

	private static int advanceHistory(final int oldHist, final int delta){
		assert(delta>=1);
		if(delta==1){return 2;}
		if(delta==2){return 1;}
		return 0;
	}

	private static int observeHistory(final int oldHist, final int diff){
		assert(diff>=1 && diff<=HISTORY_MARGIN);
		if(diff==1){return 2;}
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
			final int exp=combined/NUM_HIST;
			final int hist=combined%NUM_HIST;
			final int absNlz=exp+globalNLZ;
			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
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

	private static int emitHistMask(final int absNlz, final int hist){
		if(absNlz==0){return 0;}
		if(absNlz==1){return (hist==2) ? 2 : 0;}
		if(hist==2){return 3;}
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
		add((ArithmeticUltraDynamicLogLog33)log);
	}

	public void add(ArithmeticUltraDynamicLogLog33 log){
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
			if(oldExp>0){
				final int diff=oldExp+1;
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

	private void advanceFloor(){
		globalNLZ++;
		final int exitThreshold=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
	}

	/** Decrement all exponents by 1 across all words. Returns new minZeroCount.
	 *  Decodes all 6 buckets per 33-bit word, decrements, re-encodes. */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			long word=readWord33(w);
			if(word==0){newMinZeroCount+=6; continue;}

			final int[] vals=new int[6];
			long rem=word;
			for(int b=0; b<6; b++){
				vals[b]=(int)(rem%RADIX);
				rem/=RADIX;
			}

			for(int b=0; b<6; b++){
				final int exp=vals[b]/NUM_HIST;
				final int hist=vals[b]%NUM_HIST;
				if(exp>=1){
					final int newExp=exp-1;
					vals[b]=newExp*NUM_HIST+hist;
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
			for(int b=5; b>=0; b--){
				newWord=newWord*RADIX+vals[b];
			}
			writeWord33(w, newWord);
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
	private final int[] extraBits;
	private final int packedLen;
	private final int extraBitsLen;
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
	/** Number of exponent states: 15 (0=phantom, 1-14=real tiers). */
	static final int MAX_EXP=14;
	/** Combined states per bucket: 3×15=45. */
	static final int RADIX=NUM_HIST*(MAX_EXP+1);
	/** Tiers below floor that can still set history bits. */
	static final int HISTORY_MARGIN=2;

	/** Powers of 45 for positional decode. */
	static final long[] POW45=new long[6];
	static{
		POW45[0]=1;
		for(int i=1; i<6; i++){POW45[i]=POW45[i-1]*RADIX;}
		assert(POW45[5]*RADIX<=0x1FFFFFFFFL) : "45^6 must fit in unsigned 33-bit";
	}

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionAUDLL33.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionAUDLL33_LC2BitHistSBS.tsv.gz";
	public static final String HSB_FILE="?cardinalityCorrectionAUDLL33_LC2BitHist.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("AUDLL33: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal CFs measured from 128k DDLs, 1488 buckets, maxmult=8192. */
	@Override public float terminalMeanCF(){return 0.720997f;}
	@Override public float terminalMeanPlusCF(){return 0.894874f;}

}
