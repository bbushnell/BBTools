package cardinality;

import shared.Tools;

/**
 * ExpandedDynamicLogLog8: 8-bit packed cardinality estimator with
 * expanded (wider) tiers and 3-bit history.
 * <p>
 * Bits per bucket: 8 (4 per 32-bit int, clean byte packing)
 *   - Bits [7:3]: 5-bit tier part (0=phantom, 1-31=relTier 0-30)
 *   - Bits [2:0]: 3-bit history pattern
 *     bit 2 = observation at 1 tier below bucket
 *     bit 1 = observation at 2 tiers below bucket
 *     bit 0 = observation at 3 tiers below bucket
 * <p>
 * Tier geometry: Expanded (narrower than standard). Each tier covers 0.5 NLZ
 * (TIER_SCALE=0.5), using a 1-bit mantissa for halfNlz subdivision:
 * halfNlz = 2*NLZ + mantissa, tier = halfNlz.
 * Tier ratio = sqrt(2). 31 usable tiers cover NLZ 0-15.5.
 * <p>
 * Narrower tiers mean buckets advance twice as often as standard, keeping
 * 3-bit history varied via frequent carry-shift transitions. This is the
 * opposite of compressed-tier designs (CDLL4) where wider tiers cause
 * history to saturate before the next advance.
 * <p>
 * 3-bit history extends CDLL5's 2-bit model via carry-shift:
 * On advance by delta: newHist = ((oldHist | 8) >> delta) &amp; 7.
 * On sub-floor observation at diff tiers below bucket: hist |= 1 &lt;&lt; (3-diff).
 * HISTORY_MARGIN=3: eeMask relaxed by 3 tiers so hashes at relTier in {-3,-2,-1}
 * reach hashAndStore for history updates.
 *
 * @author Nowi, Brian Bushnell
 * @date April 2026
 */
public final class ExpandedDynamicLogLog8 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ExpandedDynamicLogLog8(){
		this(2048, 31, -1, 0);
	}

	ExpandedDynamicLogLog8(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ExpandedDynamicLogLog8(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*4<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*4-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*4;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=modBuckets/4;
		packed=new int[packedLen];
		minZeroCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+3)/4);}

	@Override
	public ExpandedDynamicLogLog8 copy(){return new ExpandedDynamicLogLog8(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	int readBucket(final int i){
		return (packed[i/4]>>>((i&3)*8))&0xFF;
	}

	private void writeBucket(final int i, final int tierPart, final int hist){
		final int idx=i/4;
		final int shift=(i&3)*8;
		final int val=((tierPart&0x1F)<<3) | (hist&0x7);
		packed[idx]=(packed[idx]&~(0xFF<<shift))|(val<<shift);
	}

	private static int tierPart(int reg){return (reg>>>3)&0x1F;}
	private static int histBits(int reg){return reg&0x7;}

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
			final int reg=readBucket(i);
			final int tp=tierPart(reg);
			final int hist=histBits(reg);
			final int absTier=tp+globalNLZ;
			if(absTier>=0 && absTier<64){
				nlzCounts[absTier+1]++;
				final int emitHist=emitHistMask(absTier, hist);
				packedBuckets[i]=(char)(((absTier+1)<<3)|emitHist);
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, HBITS, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), TIER_SCALE);
	}

	/** Mask hist bits to ones valid at the given absolute tier.
	 *  At absTier=0: all hist invalid.
	 *  At absTier=1: only bit 2 (1 tier below) valid.
	 *  At absTier=2: bits 2 and 1 valid.
	 *  At absTier>=3: all 3 bits valid. */
	private static int emitHistMask(int absTier, int hist){
		if(absTier==0){return 0;}
		if(absTier==1){return hist&0x4;}
		if(absTier==2){return hist&0x6;}
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
		add((ExpandedDynamicLogLog8)log);
	}

	public void add(ExpandedDynamicLogLog8 log){
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
				else{merged=(tpA<<3)|(histBits(rA)|histBits(rB));}
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
		return (newTp<<3)|hist;
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
		final int absTier=2*rawNlz+mantissa; // Expanded tiers: tier = halfNlz
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
				// phantom: effective tier = minZeros-1, diff = 1+(-relTier-1) = -relTier
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
			int word=packed[w];
			if(word==0){continue;}
			int result=0;
			for(int b=0; b<4; b++){
				final int shift=b*8;
				int reg=(word>>>shift)&0xFF;
				final int tp=(reg>>>3)&0x1F;
				if(tp>=1){
					reg=reg-8; // decrement tierPart by 1 (subtract 1<<3), preserve hist
					final int newTp=(reg>>>3)&0x1F;
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

	private final int[] packed;
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
	private static final int HBITS=3;
	private static final int HIST_MASK=(1<<HBITS)-1; // 7
	private static final int HIST_CARRY=1<<HBITS;     // 8
	static final int HISTORY_MARGIN=3;

	/** TIER_SCALE for CardStats: each tier covers 0.5 NLZ (expanded tiers) */
	private static final double TIER_SCALE=0.5;

	/** Mantissa threshold for half-NLZ steps: same as CDLL5 (2-sqrt(2)) * 1048576 */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	// CF/SBS file references — HSB via hardcoded arrays in StateTable.CF_ETLL_3
	public static final String CF_FILE="?cardinalityCorrectionEDLL8.tsv.gz";
	public static final String SBS_FILE="?cardinalityCorrectionEDLL8_LC3BitHistSBS.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	@Override public float terminalMeanCF(){return 0.597305f;}
	@Override public float terminalMeanPlusCF(){return 0.582968f;}

}
