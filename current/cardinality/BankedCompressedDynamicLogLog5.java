package cardinality;

import shared.Tools;

/**
 * BankedCompressedDynamicLogLog5: CDLL5 with per-word banked exponents.
 * <p>
 * Layout: 6 buckets × 5-bit (3-bit tierPart + 2-bit hist) = 30 bits
 * + 2-bit bank exponent at bits 30-31 = 32 bits total (no waste).
 * <p>
 * Absolute tier = (tierPart - 1) + minZeros + bankExponent.
 * Bank promotion: when a register would overflow and bank < 3,
 * if all 6 registers are >= 1, subtract 1 from all and increment bank.
 * Global promotion absorption: words with bank > 0 just decrement bank.
 *
 * @author Eru, Brian Bushnell
 * @date April 2026
 */
public final class BankedCompressedDynamicLogLog5 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	BankedCompressedDynamicLogLog5(){
		this(2048, 31, -1, 0);
	}

	BankedCompressedDynamicLogLog5(parse.Parser p){
		super(p);
		packedLen=(buckets+5)/6;
		packed=new int[packedLen];
		minZeroCount=buckets;
	}

	BankedCompressedDynamicLogLog5(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		packedLen=(buckets+5)/6;
		packed=new int[packedLen];
		minZeroCount=buckets;
	}

	@Override
	public BankedCompressedDynamicLogLog5 copy(){return new BankedCompressedDynamicLogLog5(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	int readBucket(final int i){
		return (packed[i/6]>>>((i%6)*5))&0x1F;
	}

	private void writeBucket(final int i, final int tierPart, final int hist){
		final int idx=i/6;
		final int shift=(i%6)*5;
		final int val=((tierPart&0x7)<<2) | (hist&0x3);
		packed[idx]=(packed[idx]&~(0x1F<<shift))|(val<<shift);
	}

	private static int tierPart(int reg){return (reg>>>2)&0x7;}
	private static int histBits(int reg){return reg&0x3;}

	private int readBank(final int wordIdx){
		return (packed[wordIdx]>>>30)&0x3;
	}

	private void writeBank(final int wordIdx, final int val){
		packed[wordIdx]=(packed[wordIdx]&0x3FFFFFFF)|((val&0x3)<<30);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		if(packedBuckets==null){packedBuckets=new char[buckets];}
		else{java.util.Arrays.fill(packedBuckets, (char)0);}

		int filledCount=0;
		for(int i=0; i<buckets; i++){
			final int wordIdx=i/6;
			final int bank=readBank(wordIdx);
			final int reg=readBucket(i);
			final int tp=tierPart(reg);
			final int hist=histBits(reg);
			if(tp>0){
				final int absTier=(tp-1)+minZeros+bank;
				if(absTier<64){
					nlzCounts[absTier+1]++;
					final int emitHist=emitHistMask(absTier, hist);
					packedBuckets[i]=(char)(((absTier+1)<<2)|emitHist);
				}
				filledCount++;
			}else if(minZeros+bank>0){
				final int phantomTier=minZeros+bank-1;
				if(phantomTier>=0 && phantomTier<64){
					nlzCounts[phantomTier+1]++;
					final int emitHist=emitHistMask(phantomTier, hist);
					packedBuckets[i]=(char)(((phantomTier+1)<<2)|emitHist);
					filledCount++;
				}
			}
		}
		nlzCounts[0]=buckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
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
		add((BankedCompressedDynamicLogLog5)log);
	}

	public void add(BankedCompressedDynamicLogLog5 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int w=0; w<packedLen; w++){
				writeBank(w, 0);
			}
			for(int i=0; i<buckets; i++){
				final int wA=i/6, wB=i/6;
				final int bankA=readBank(wA);
				final int bankB=log.readBank(wB);
				final int regA=readBucket(i), regB=log.readBucket(i);
				final int tpA=tierPart(regA), tpB=tierPart(regB);
				final int absA=(tpA==0 ? -1 : (tpA-1)+minZeros+bankA);
				final int absB=(tpB==0 ? -1 : (tpB-1)+log.minZeros+bankB);
				final int absMax=Math.max(absA, absB);
				if(absMax<newMinZeros){
					writeBucket(i, 0, 0);
				}else{
					final int newTp=Math.min(absMax-newMinZeros+1, 7);
					final int hA=histBits(regA), hB=histBits(regB);
					final int mergedHist=(absA==absB) ? (hA|hB) : (absA>absB ? hA : hB);
					writeBucket(i, newTp, mergedHist);
				}
			}
			minZeros=newMinZeros;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int tp=tierPart(readBucket(i));
				if(tp>0){filledBuckets++;}
				if(tp==0 || (!EARLY_PROMOTE && tp==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && minZeros<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
		if(!DISABLE_EEMASK && Long.compareUnsigned(key, eeMask)>0){return;}

		final int rawNlz=Long.numberOfLeadingZeros(key);
		final int mantissa;
		if(rawNlz>=47){
			mantissa=0;
		}else{
			final int mBits=(int)((key>>>(46-rawNlz))&0xFFFF);
			mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		}
		final int halfNlz=2*rawNlz+mantissa;
		final int absTier=halfNlz/3;

		final int bucket=(int)(key&bucketMask);
		final int wordIdx=bucket/6;
		final int bank=readBank(wordIdx);
		final int relTier=absTier-minZeros-bank;

		if(relTier<-HISTORY_MARGIN){return;}

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int oldReg=readBucket(bucket);
		final int oldTierPart=tierPart(oldReg);
		final int oldHist=histBits(oldReg);

		if(relTier>=0){
			int newTierPart=Math.min(relTier+1, 7);

			// Try bank promotion if would overflow
			if(newTierPart>7 && bank<3 && canPromoteBank(wordIdx)){
				promoteBank(wordIdx);
				final int newBank=readBank(wordIdx);
				final int newRelTier=absTier-minZeros-newBank;
				newTierPart=Math.min(newRelTier+1, 7);
				// Re-read old bucket state after promotion
				final int promotedReg=readBucket(bucket);
				final int promotedTp=tierPart(promotedReg);
				final int promotedHist=histBits(promotedReg);
				if(newTierPart<=promotedTp){return;}

				final int delta=(promotedTp>0) ? (newTierPart-promotedTp) : (newRelTier+1);
				final int carry=(promotedTp>0 || minZeros+newBank>0) ? HIST_CARRY : 0;
				final int newHist=((promotedHist|carry)>>delta)&HIST_MASK;
				lastCardinality=-1;
				writeBucket(bucket, newTierPart, newHist);
				return;
			}

			if(newTierPart>oldTierPart){
				final int delta=(oldTierPart>0) ? (newTierPart-oldTierPart) : (relTier+1);
				final int carry=(oldTierPart>0 || minZeros+bank>0) ? HIST_CARRY : 0;
				final int newHist=((oldHist|carry)>>delta)&HIST_MASK;
				lastCardinality=-1;
				writeBucket(bucket, newTierPart, newHist);
				if(oldTierPart==0 && bank==0){filledBuckets++;}

				final int oldRelTier=(oldTierPart==0 ? 0 : oldTierPart-1);
				final boolean shouldDecrement=EARLY_PROMOTE
					? (oldTierPart==0 && bank==0)
					: (relTier>oldRelTier && oldRelTier==0 && bank==0);
				if(shouldDecrement && --minZeroCount<1){
					while(minZeroCount==0 && minZeros<wordlen){
						advanceFloor();
						minZeroCount=countAndDecrement();
					}
				}
				return;
			}
			if(newTierPart<oldTierPart){
				final int diff=oldTierPart-newTierPart;
				if(diff>=1 && diff<=HBITS){
					final int bit=HBITS-diff;
					final int newHist=oldHist|(1<<bit);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldTierPart, newHist);
					}
				}
			}
		}else if(relTier==-1){
			if(oldTierPart>0){
				final int diff=oldTierPart;
				if(diff>=1 && diff<=HBITS){
					final int bit=HBITS-diff;
					final int newHist=oldHist|(1<<bit);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldTierPart, newHist);
					}
				}
			}
		}else{
			// relTier == -2
			final int diff=(oldTierPart>0) ? (oldTierPart+1) : 1;
			if(diff>=1 && diff<=HBITS){
				final int bit=HBITS-diff;
				final int newHist=oldHist|(1<<bit);
				if(newHist!=oldHist){
					lastCardinality=-1;
					writeBucket(bucket, oldTierPart, newHist);
				}
			}
		}
	}

	private boolean canPromoteBank(final int wordIdx){
		final int word=packed[wordIdx]&0x3FFFFFFF;
		for(int b=0; b<6; b++){
			final int reg=(word>>>(b*5))&0x1F;
			if(tierPart(reg)==0){return false;}
		}
		return true;
	}

	private void promoteBank(final int wordIdx){
		final int oldBank=readBank(wordIdx);
		int word=packed[wordIdx]&0x3FFFFFFF;
		int result=0;
		for(int b=0; b<6; b++){
			final int shift=b*5;
			int reg=(word>>>shift)&0x1F;
			// Decrement tierPart by 1 (subtract 4), preserve hist bits
			reg=reg-4;
			result|=(reg<<shift);
		}
		packed[wordIdx]=result|((oldBank+1)<<30);
	}

	private void advanceFloor(){
		minZeros++;
		final int relaxedTier=Math.max(0, minZeros-HISTORY_MARGIN);
		final int minNlz=(3*relaxedTier)/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			final int bank=readBank(w);
			if(bank>0){
				writeBank(w, bank-1);
				if(bank==1){
					// Bank dropping to 0: count any stored=0 slots
					for(int b=0; b<6; b++){
						final int idx=w*6+b;
						if(idx>=buckets){break;}
						final int tp=tierPart(readBucket(idx));
						if(EARLY_PROMOTE){
							if(tp==0){newMinZeroCount++;}
						}else{
							if(tp==0 || tp==1){newMinZeroCount++;}
						}
					}
				}
			}else{
				// bank==0: decrement registers
				int word=packed[w];
				if((word&0x3FFFFFFF)==0){
					final int regsInWord=Math.min(6, buckets-w*6);
					newMinZeroCount+=regsInWord;
					continue;
				}
				int result=0;
				for(int b=0; b<6; b++){
					final int shift=b*5;
					int reg=(word>>>shift)&0x1F;
					final int tp=(reg>>>2)&0x7;
					if(tp>=1){
						reg=reg-4; // decrement tierPart, preserve hist
						final int newTp=(reg>>>2)&0x7;
						if(EARLY_PROMOTE){
							if(newTp==0){newMinZeroCount++; filledBuckets--;}
						}else{
							if(newTp<=1){newMinZeroCount++;}
						}
					}
					result|=(reg<<shift);
				}
				packed[w]=result;
			}
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){
		if(minZeros>0){return 1.0;}
		int empty=0;
		for(int w=0; w<packedLen; w++){
			if(readBank(w)>0){continue;}
			for(int b=0; b<6; b++){
				final int idx=w*6+b;
				if(idx>=buckets){break;}
				if(tierPart(readBucket(idx))==0){empty++;}
			}
		}
		return 1.0-(double)empty/buckets;
	}

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

	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*65536);

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionBCDLL5.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionCDLL5_LC2BitHistSBS.tsv.gz";
	public static final String HSB_FILE="?cardinalityCorrectionCDLL5_LC2BitHist.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("BCDLL5: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	@Override public float terminalMeanCF(){return 0.883501f;}
	@Override public float terminalMeanPlusCF(){return 1.094262f;}

}
