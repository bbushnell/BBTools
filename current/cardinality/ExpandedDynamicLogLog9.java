package cardinality;

import shared.Tools;

/**
 * ExpandedDynamicLogLog9: 9-bit packed cardinality estimator with
 * expanded tiers and 4-bit history.
 * <p>
 * Bits per bucket: 9 (7 per 64-bit long, 1 bit wasted per word)
 *   - Bits [8:4]: 5-bit tier part (0=phantom, 1-31=relTier 0-30)
 *   - Bits [3:0]: 4-bit history pattern
 *     bit 3 = observation at 1 tier below bucket
 *     bit 2 = observation at 2 tiers below bucket
 *     bit 1 = observation at 3 tiers below bucket
 *     bit 0 = observation at 4 tiers below bucket
 * <p>
 * Tier geometry: same as EDLL8. Each tier covers 0.5 NLZ (TIER_SCALE=0.5),
 * halfNlz = 2*NLZ + mantissa. Tier ratio = sqrt(2).
 * <p>
 * 4-bit history via carry-shift:
 * On advance by delta: newHist = ((oldHist | 16) >> delta) &amp; 0xF.
 * On sub-floor observation at diff tiers below: hist |= 1 &lt;&lt; (4-diff).
 * HISTORY_MARGIN=4.
 * <p>
 * At 1KB: 128 longs × 7 = 896 buckets (vs EDLL8's 1024).
 *
 * @author Nowi, Brian Bushnell
 * @date April 2026
 */
public final class ExpandedDynamicLogLog9 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ExpandedDynamicLogLog9(){
		this(2048, 31, -1, 0);
	}

	ExpandedDynamicLogLog9(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ExpandedDynamicLogLog9(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*BUCKETS_PER_WORD<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*BUCKETS_PER_WORD-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*BUCKETS_PER_WORD;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=modBuckets/BUCKETS_PER_WORD;
		packed=new long[packedLen];
		minZeroCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+BUCKETS_PER_WORD-1)/BUCKETS_PER_WORD);}

	@Override
	public ExpandedDynamicLogLog9 copy(){return new ExpandedDynamicLogLog9(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	int readBucket(final int i){
		return (int)((packed[i/BUCKETS_PER_WORD]>>>((i%BUCKETS_PER_WORD)*BITS_PER_BUCKET))&BUCKET_MASK);
	}

	private void writeBucket(final int i, final int tierPart, final int hist){
		final int idx=i/BUCKETS_PER_WORD;
		final int shift=(i%BUCKETS_PER_WORD)*BITS_PER_BUCKET;
		final long val=(long)(((tierPart&0x1F)<<HBITS) | (hist&HIST_MASK));
		packed[idx]=(packed[idx]&~(BUCKET_MASK<<shift))|(val<<shift);
	}

	private static int tierPart(int reg){return (reg>>>HBITS)&0x1F;}
	private static int histBits(int reg){return reg&HIST_MASK;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		if(packedBuckets==null){packedBuckets=new char[modBuckets];}
		else{java.util.Arrays.fill(packedBuckets, (char)0);}

		final int useHbits=EDLL_HBITS;
		final int shift=HBITS-useHbits;
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=readBucket(i);
			final int tp=tierPart(reg);
			final int hist=histBits(reg);
			final int absTier=tp+globalNLZ;
			if(absTier>=0 && absTier<64){
				nlzCounts[absTier+1]++;
				final int emitHist=emitHistMask(absTier, hist);
				final int usedHist=(useHbits>=HBITS) ? emitHist : (emitHist>>shift);
				packedBuckets[i]=(char)(((absTier+1)<<useHbits)|usedHist);
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, useHbits, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), TIER_SCALE);
	}

	/** Mask hist bits to ones valid at the given absolute tier.
	 *  At absTier 0: none valid. At absTier k (1-3): top k bits valid.
	 *  At absTier >= 4: all 4 bits valid. */
	private static int emitHistMask(int absTier, int hist){
		if(absTier<=0){return 0;}
		if(absTier>=HBITS){return hist;}
		return hist&(HIST_MASK<<(HBITS-absTier));
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
		add((ExpandedDynamicLogLog9)log);
	}

	public void add(ExpandedDynamicLogLog9 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			for(int i=0; i<modBuckets; i++){
				final int rA=adjustBucket(readBucket(i), globalNLZ, newGlobalNLZ);
				final int rB=adjustBucket(log.readBucket(i), log.globalNLZ, newGlobalNLZ);
				final int tpA=tierPart(rA), tpB=tierPart(rB);
				final int merged;
				if(tpA>tpB){merged=rA;}
				else if(tpB>tpA){merged=rB;}
				else{merged=(tpA<<HBITS)|(histBits(rA)|histBits(rB));}
				writeBucket(i, tierPart(merged), histBits(merged));
			}
			globalNLZ=newGlobalNLZ;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
				final int tp=tierPart(readBucket(i));
				if(tp>0){filledBuckets++;}
				if(tp==0 || (!EARLY_PROMOTE && tp==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && globalNLZ<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	private static int adjustBucket(int reg, int oldFloor, int newFloor){
		final int tp=tierPart(reg);
		if(tp==0){return reg;}
		final int hist=histBits(reg);
		final int newTp=Math.max(0, Math.min(tp+(oldFloor-newFloor), 31));
		return (newTp<<HBITS)|hist;
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
		final int absTier=2*rawNlz+mantissa;
		final int relTier=absTier-(globalNLZ+1);
		if(relTier<-HISTORY_MARGIN){return;}

		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int oldReg=readBucket(bucket);
		final int oldTierPart=tierPart(oldReg);
		final int oldHist=histBits(oldReg);

		if(relTier>=0){
			final int newTierPart=Math.min(relTier+1, 31);
			if(newTierPart>oldTierPart){
				final int delta=(oldTierPart>0) ? (newTierPart-oldTierPart) : (relTier+1);
				final int carry=(oldTierPart>0 || globalNLZ>=0) ? HIST_CARRY : 0;
				final int newHist=((oldHist|carry)>>delta)&HIST_MASK;
				lastCardinality=-1;
				writeBucket(bucket, newTierPart, newHist);
				if(oldTierPart==0){filledBuckets++;}

				final int oldRelTier=(oldTierPart==0 ? 0 : oldTierPart-1);
				final boolean shouldDecrement=EARLY_PROMOTE ? oldTierPart==0
					: (relTier>oldRelTier && oldRelTier==0);
				if(shouldDecrement && --minZeroCount<1){
					while(minZeroCount==0 && globalNLZ<wordlen){
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
		}else if(relTier>=-HISTORY_MARGIN){
			if(oldTierPart>0){
				final int diff=oldTierPart+(-relTier);
				if(diff>=1 && diff<=HBITS){
					final int bit=HBITS-diff;
					final int newHist=oldHist|(1<<bit);
					if(newHist!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldTierPart, newHist);
					}
				}
			}else{
				final int diff=-relTier;
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
	}

	private void advanceFloor(){
		globalNLZ++;
		final int relaxedTier=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		final int minNlz=relaxedTier/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			long word=packed[w];
			if(word==0){continue;}
			long result=0;
			for(int b=0; b<BUCKETS_PER_WORD; b++){
				final int shift=b*BITS_PER_BUCKET;
				int reg=(int)((word>>>shift)&BUCKET_MASK);
				final int tp=(reg>>>HBITS)&0x1F;
				if(tp>=1){
					reg=reg-(1<<HBITS);
					final int newTp=(reg>>>HBITS)&0x1F;
					if(EARLY_PROMOTE){
						if(newTp==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(newTp<=1){newMinZeroCount++;}
					}
				}
				result|=((long)reg<<shift);
			}
			packed[w]=result;
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

	private final long[] packed;
	private final int packedLen;
	private final int modBuckets;
	private int globalNLZ=-1;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return globalNLZ+1;}
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
	private static final int BITS_PER_BUCKET=9;
	private static final int BUCKETS_PER_WORD=7;
	private static final long BUCKET_MASK=(1L<<BITS_PER_BUCKET)-1; // 0x1FF
	private static final int HBITS=4;
	private static final int HIST_MASK=(1<<HBITS)-1; // 0xF
	private static final int HIST_CARRY=1<<HBITS;     // 16
	static final int HISTORY_MARGIN=4;

	private static final double TIER_SCALE=0.5;
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	public static int EDLL_HBITS=4;
	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionEDLL9.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionEDLL9_LC4BitHistSBS.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	@Override public float terminalMeanCF(){return 0.596929f;}
	@Override public float terminalMeanPlusCF(){return 0.751432f;}

	@Override public double ldlcHcWeight(){return 0.68;}
	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.68f;}


}
