package cardinality;

import shared.Tools;

/**
 * CompressedDynamicLogLog5: 5-bit packed cardinality estimator with
 * mantissa-compressed tiers and 2-bit history.
 * <p>
 * Bits per bucket: 5 (6 per 32-bit int, 2 bits wasted)
 *   - Bits [4:2]: tier part (0=phantom, 1-7=relTier 0-6)
 *   - Bits [1:0]: 2-bit history pattern
 *     bit 1 = observation at 1 tier below bucket
 *     bit 0 = observation at 2 tiers below bucket
 * <p>
 * Tier compression: CDLL4's 2-sqrt(2) half-NLZ geometry.
 * halfNlz = 2*NLZ + mantissa, tier = halfNlz/3. TIER_SCALE = 1.5.
 * <p>
 * 2-bit history extends CDLL4's single-bit model via UDLL6's carry-shift:
 * On advance by delta: newHist = ((oldHist | 4) &gt;&gt; delta) &amp; 3.
 * On sub-floor observation at diff tiers below bucket: hist |= 1 &lt;&lt; (2-diff).
 * HISTORY_MARGIN=2: eeMask relaxed by 2 tiers so hashes at relTier in {-1,-2}
 * reach hashAndStore for history updates.
 * <p>
 * At same memory as CDLL4 (4-bit nibbles): CDLL5 has 20% fewer buckets
 * (6 per int vs 8 per int) in exchange for 2x history capacity. Targets
 * beating CDLL4 on low-complexity data where BDLL5 fails.
 *
 * @author Eru, Brian Bushnell
 * @date April 2026
 */
public final class CompressedDynamicLogLog5 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	CompressedDynamicLogLog5(){
		this(2048, 31, -1, 0);
	}

	CompressedDynamicLogLog5(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	CompressedDynamicLogLog5(int buckets_, int k_, long seed, float minProb_){
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
	public CompressedDynamicLogLog5 copy(){return new CompressedDynamicLogLog5(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read 5-bit bucket i (0..buckets-1). */
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
			final int reg=readBucket(i);
			final int tp=tierPart(reg);
			final int hist=histBits(reg);
			if(tp>0){
				final int absTier=(tp-1)+minZeros;
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

	/** Mask hist bits to ones valid at the given absolute tier.
	 *  At absTier=0: all hist invalid (-1, -2 don't exist).
	 *  At absTier=1: only bit 1 (-1 = 0, valid) remains.
	 *  At absTier>=2: both bits valid. */
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
		add((CompressedDynamicLogLog5)log);
	}

	public void add(CompressedDynamicLogLog5 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<modBuckets; i++){
				final int rA=adjustBucket(readBucket(i), minZeros, newMinZeros);
				final int rB=adjustBucket(log.readBucket(i), log.minZeros, newMinZeros);
				final int tpA=tierPart(rA), tpB=tierPart(rB);
				final int merged;
				if(tpA>tpB){merged=rA;}
				else if(tpB>tpA){merged=rB;}
				else{merged=(tpA<<2)|(histBits(rA)|histBits(rB));}
				writeBucket(i, tierPart(merged), histBits(merged));
			}
			minZeros=newMinZeros;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
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

	private static int adjustBucket(int reg, int oldFloor, int newFloor){
		final int tp=tierPart(reg);
		if(tp==0){return reg;}
		final int hist=histBits(reg);
		final int newTp=Math.max(0, Math.min(tp+(oldFloor-newFloor), 7));
		return (newTp<<2)|hist;
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
		// HISTORY_MARGIN=2: accept relTier in {-2, -1, 0, 1, ...}
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
			final int newTierPart=Math.min(relTier+1, 7);
			if(newTierPart>oldTierPart){
				// Tier advance. Unified delta: phantom's effective tier is minZeros-1,
				// so delta = absTier - (minZeros-1) = relTier+1 when oldTierPart==0.
				final int delta=(oldTierPart>0) ? (newTierPart-oldTierPart) : (relTier+1);
				// Carry-shift matches MantissaCompare2 MODE_CTLL (lines 213-215):
				// simulator uses carry=0 when stored==0 (no prior observation).
				// Equivalent here: suppress carry for a truly-fresh bucket
				// (oldTierPart==0 AND minZeros==0). Phantoms from floor advance
				// legitimately carry an observation at minZeros-1, so carry stays.
				final int carry=(oldTierPart>0 || minZeros>0) ? HIST_CARRY : 0;
				final int newHist=((oldHist|carry)>>delta)&HIST_MASK;
				lastCardinality=-1;
				writeBucket(bucket, newTierPart, newHist);
				if(oldTierPart==0){filledBuckets++;}

				final int oldRelTier=(oldTierPart==0 ? 0 : oldTierPart-1);
				final boolean shouldDecrement=EARLY_PROMOTE ? oldTierPart==0
					: (relTier>oldRelTier && oldRelTier==0);
				if(shouldDecrement && --minZeroCount<1){
					while(minZeroCount==0 && minZeros<wordlen){
						advanceFloor();
						minZeroCount=countAndDecrement();
					}
				}
				return;
			}
			if(newTierPart<oldTierPart){
				// Sub-tier observation: diff = oldTierPart - newTierPart.
				// newTierPart=0 shouldn't happen when relTier>=0, but guard anyway.
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
			// newTierPart == oldTierPart: same tier, no update
		}else if(relTier==-1){
			// Observation at minZeros-1. For oldTierPart>0 (real bucket at
			// absTier = oldTierPart-1+minZeros), diff = oldTierPart.
			// For oldTierPart=0 (phantom at effective minZeros-1), diff=0 — no-op.
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
			// relTier == -2: observation at minZeros-2.
			// For oldTierPart>0: diff = oldTierPart+1.
			// For oldTierPart=0 (phantom at minZeros-1): diff = 1.
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

	/** Advance floor by one tier; relax eeMask by HISTORY_MARGIN tiers so
	 *  hashes at relTier in [-HISTORY_MARGIN, -1] pass through for history. */
	private void advanceFloor(){
		minZeros++;
		final int relaxedTier=Math.max(0, minZeros-HISTORY_MARGIN);
		final int minNlz=(3*relaxedTier)/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			int word=packed[w];
			if(word==0){continue;}
			int result=0;
			for(int b=0; b<6; b++){
				final int shift=b*5;
				int reg=(word>>>shift)&0x1F;
				final int tp=(reg>>>2)&0x7;
				if(tp>=1){
					// Decrement tierPart by 1 (subtract 4 from 5-bit reg), preserve hist
					reg=reg-4;
					final int newTp=(reg>>>2)&0x7;
					if(EARLY_PROMOTE){
						if(newTp==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(newTp<=1){newMinZeroCount++;}
					}
				}
				// tp==0: already phantom, leave hist intact
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

	/** Mantissa threshold for log-uniform half-NLZ steps: (2-sqrt(2)) * 1048576 */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE="?cardinalityCorrectionCDLL5.tsv.gz";
	/** SBS (per-fill LC correction) table, 2-bit history with half-NLZ tier geometry. */
	public static final String SBS_FILE="?cardinalityCorrectionCDLL5_LC2BitHistSBS.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("CDLL5: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Measured 2026-04-20: 128k DDLs, 1536 modulo buckets, maxmult=8192, tmcf=1 tmpcf=1. */
	@Override public float terminalMeanCF(){return 0.883080f;}
	@Override public float terminalMeanPlusCF(){return 1.094268f;}

}
