package cardinality;

import shared.Tools;

/**
 * CompressedDynamicLogLog4: 4-bit nibble-packed cardinality estimator
 * with mantissa-compressed tiers and 1-bit history.
 * <p>
 * Bits per bucket: 4 (8 per int, nibble-packed)
 *   - Bits [3:1]: tier part (0=phantom, 1-7=relTier 0-6)
 *   - Bit [0]: history bit (sub-tier observation marker)
 * <p>
 * Tier compression: single hash with 16-bit mantissa threshold at
 * (2-sqrt(2)) fraction for log-uniform half-NLZ steps.
 * halfNlz = 2*NLZ + mantissa, tier = halfNlz/3.
 * Each tier = 2*sqrt(2) probability ratio.
 * TIER_SCALE = 1.5 for DLC formulas.
 * <p>
 * History bit: set when an element's sub-tier position (halfNlz%3)
 * reaches 2 (top of tier). Enables LDLC estimation via CardStats Phase 2a.
 * <p>
 * At same memory as UDLL6 (6-bit), CDLL4 has 50% more buckets.
 * Error scales as 1/sqrt(B), so the bucket advantage partially
 * compensates for fewer history bits (1 vs 2).
 *
 * @author Brian Bushnell, Chloe
 * @date April 2026
 */
public final class CompressedDynamicLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	CompressedDynamicLogLog4(){
		this(2048, 31, -1, 0);
	}

	CompressedDynamicLogLog4(parse.Parser p){
		super(p);
		maxArray=new int[(buckets+7)>>>3];
		minZeroCount=buckets;
	}

	CompressedDynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[(buckets+7)>>>3];
		minZeroCount=buckets;
	}

	@Override
	public CompressedDynamicLogLog4 copy(){return new CompressedDynamicLogLog4(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Nibble Access         ----------------*/
	/*--------------------------------------------------------------*/

	int readNibble(final int i){
		return (maxArray[i>>>3]>>>((i&7)<<2))&0xF;
	}

	private void writeNibble(final int i, final int val){
		final int wordIdx=i>>>3;
		final int shift=(i&7)<<2;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0xF<<shift))|((val&0xF)<<shift);
	}

	private static int tierPart(int nibble){return nibble>>>1;}
	private static int histBit(int nibble){return nibble&1;}
	private static int makeNibble(int tierPart, int hist){return (tierPart<<1)|hist;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		if(packedBuckets==null){packedBuckets=new char[buckets];}
		else{java.util.Arrays.fill(packedBuckets, (char)0);}

		final int phantomTier=minZeros-1;
		int filledCount=0;
		for(int i=0; i<buckets; i++){
			final int nib=readNibble(i);
			final int tp=tierPart(nib);
			final int hist=histBit(nib);
			if(tp>0){
				final int absTier=(tp-1)+minZeros;
				if(absTier<64){
					nlzCounts[absTier+1]++;
					packedBuckets[i]=(char)(((absTier+1)<<1)|hist);
				}
				filledCount++;
			}else if(minZeros>0 && phantomTier>=0 && phantomTier<64){
				nlzCounts[phantomTier+1]++;
				packedBuckets[i]=(char)(((phantomTier+1)<<1)|hist);
				filledCount++;
			}
		}
		nlzCounts[0]=buckets-filledCount;

		// Debug: count history bit distribution
		if(DEBUG_HIST && filledCount>0){
			int h0=0, h1=0;
			for(int i=0; i<buckets; i++){
				int nib=readNibble(i);
				if(tierPart(nib)>0){
					if(histBit(nib)==0) h0++; else h1++;
				}
			}
			System.err.println("CDLL4 hist: h0="+h0+" h1="+h1+
				" P(0)="+String.format("%.3f", (double)h0/(h0+h1))+
				" minZeros="+minZeros+" filled="+filledCount);
		}

		final double savedScale=AbstractCardStats.TIER_SCALE;
		AbstractCardStats.TIER_SCALE=1.5;
		try{
			return new CardStats(packedBuckets, nlzCounts, 0, 1, 0, 0,
					buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0);
		}finally{
			AbstractCardStats.TIER_SCALE=savedScale;
		}
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
		add((CompressedDynamicLogLog4)log);
	}

	public void add(CompressedDynamicLogLog4 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int nA=adjustNibble(readNibble(i), minZeros, newMinZeros);
				final int nB=adjustNibble(log.readNibble(i), log.minZeros, newMinZeros);
				final int tpA=tierPart(nA), tpB=tierPart(nB);
				final int merged;
				if(tpA>tpB){merged=nA;}
				else if(tpB>tpA){merged=nB;}
				else{merged=makeNibble(tpA, histBit(nA)|histBit(nB));}
				writeNibble(i, merged);
			}
			minZeros=newMinZeros;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int tp=tierPart(readNibble(i));
				if(tp>0){filledBuckets++;}
				if(tp==0 || (!EARLY_PROMOTE && tp==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && minZeros<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	private static int adjustNibble(int nib, int oldFloor, int newFloor){
		if(nib<=1){return nib;} // phantom: preserve as-is
		final int tp=tierPart(nib);
		final int hist=histBit(nib);
		final int newTp=Math.max(0, Math.min(tp+(oldFloor-newFloor), 7));
		return makeNibble(newTp, hist);
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
		if(!DISABLE_EEMASK && Long.compareUnsigned(key, eeMask)>0){cntEarlyExit++; return;}

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
		final int relTier=absTier-minZeros;
		// Allow relTier == -1 through: needed for history (delta == -1)
		if(relTier<-1){cntBelowFloor++; return;}

		final int bucket=(int)(key&bucketMask);
		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int oldNibble=readNibble(bucket);
		final int oldTierPart=tierPart(oldNibble);
		final int oldHist=histBit(oldNibble);

		if(relTier>=0){
			final int newTierPart=Math.min(relTier+1, 7);
			if(newTierPart>oldTierPart){
				// Tier advance: delta==1 → set hist, delta>=2 → clear hist.
				// For phantoms (tierPart=0): effective old tier = minZeros-1,
				// so delta = absTier - (minZeros-1) = relTier+1.
				final int delta=(oldTierPart>0) ? (newTierPart-oldTierPart) : (relTier+1);
				final int newHist=(delta==1) ? 1 : 0;
				if(newHist==1){cntAdvanceSet++;}else{cntAdvanceClear++;}
				lastCardinality=-1;
				writeNibble(bucket, makeNibble(newTierPart, newHist));
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
			// delta == -1: element from previous tier, set history
			if(newTierPart==oldTierPart-1){
				if(oldHist==0){cntDeltaNeg1Set++; lastCardinality=-1; writeNibble(bucket, oldNibble|1);}
				else{cntDeltaNeg1Skip++;}
			}
		}else{
			// relTier == -1: element from one tier below floor
			if(oldTierPart<=1){
				if(oldHist==0){cntDeltaNeg1Set++; lastCardinality=-1; writeNibble(bucket, oldNibble|1);}
				else{cntDeltaNeg1Skip++;}
			}
		}
	}

	/** Advance the global floor by one tier and update eeMask.
	 *  eeMask is relaxed by 1 tier to allow delta==-1 elements through for history. */
	private void advanceFloor(){
		minZeros++;
		// Relaxed: use minNlz for (minZeros-1) instead of minZeros
		final int relaxedTier=Math.max(0, minZeros-1);
		final int minNlz=(3*relaxedTier)/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<maxArray.length; w++){
			int word=maxArray[w];
			if(word==0){continue;}
			int result=0;
			for(int b=0; b<8; b++){
				final int shift=b<<2;
				int nib=(word>>>shift)&0xF;
				if(nib>=2){
					// Subtract 2: decrements tierPart by 1, preserves history bit
					nib-=2;
					final int tp=nib>>>1;
					if(EARLY_PROMOTE){
						if(tp==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(tp<=1){newMinZeroCount++;}
					}
				}
				// nib<2: already phantom, leave as-is
				result|=(nib<<shift);
			}
			maxArray[w]=result;
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

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

	private final int[] maxArray;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return minZeros;}
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
	/** Mantissa threshold for log-uniform half-NLZ steps: (2-sqrt(2)) * 65536 */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*65536);

	public static boolean EARLY_PROMOTE=true;
	public static boolean DEBUG_HIST=false;
	public static boolean DISABLE_EEMASK=false;

	// Debug counters (valid for ddls=1 only)
	public static long cntAdvanceSet=0;     // carry set history=1 during tier advance (delta==1)
	public static long cntAdvanceClear=0;   // carry cleared history during tier advance (delta>=2)
	public static long cntDeltaNeg1Set=0;   // delta==-1 set history (was 0)
	public static long cntDeltaNeg1Skip=0;  // delta==-1 skipped (already 1)
	public static long cntEarlyExit=0;      // filtered by eeMask
	public static long cntBelowFloor=0;     // filtered by relTier < -1

	public static void printDebugCounters(){
		System.err.println("=== CDLL4 History Debug ===");
		System.err.println("  advance set hist=1 (delta=1): "+cntAdvanceSet);
		System.err.println("  advance clear hist (delta>=2): "+cntAdvanceClear);
		System.err.println("  delta=-1 set hist: "+cntDeltaNeg1Set);
		System.err.println("  delta=-1 skip (already set): "+cntDeltaNeg1Skip);
		System.err.println("  early exit (eeMask): "+cntEarlyExit);
		System.err.println("  below floor (relTier<-1): "+cntBelowFloor);
	}

	public static final String CF_FILE="?cardinalityCorrectionCDLL4.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("CDLL4: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}

}
