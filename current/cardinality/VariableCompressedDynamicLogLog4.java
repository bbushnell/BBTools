package cardinality;

import shared.Tools;

/**
 * VariableCompressedDynamicLogLog4: 4-bit nibble-packed cardinality estimator
 * with mantissa-compressed tiers and variable history bits.
 * <p>
 * 16 states per nibble:
 *   - States 0-11 (nibbles 0x0-0xB): tierPart 0-5 with 1-bit history
 *     tierPart = nibble&gt;&gt;&gt;1, histBit = nibble&amp;1
 *   - States 12-15 (nibbles 0xC-0xF): tierPart 6-9 without history
 *     tierPart = nibble - 6
 * <p>
 * This gives 2 additional tiers vs CDLL4 (max tierPart 9 vs 7),
 * reducing overflow at the cost of losing history at high tiers.
 * <p>
 * Three demotion sub-variants control what happens when a bucket at
 * tierPart 6 (no history) demotes to tierPart 5 (has history) during
 * countAndDecrement:
 *   - Variant A (DEMOTION_MODE=0): hist = 1 (assume recently seen)
 *   - Variant B (DEMOTION_MODE=1): hist = 0 (assume not recently seen)
 *   - Variant C (DEMOTION_MODE=2): hist = bucketID &amp; 1 (deterministic noise)
 * <p>
 * Tier compression: single hash with 16-bit mantissa threshold at
 * (2-sqrt(2)) fraction for log-uniform half-NLZ steps.
 * halfNlz = 2*NLZ + mantissa, tier = halfNlz/3.
 * TIER_SCALE = 1.5 for DLC formulas.
 *
 * @author Noire
 * @date April 2026
 */
public final class VariableCompressedDynamicLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	VariableCompressedDynamicLogLog4(){
		this(2048, 31, -1, 0);
	}

	VariableCompressedDynamicLogLog4(parse.Parser p){
		super(p);
		maxArray=new int[(buckets+7)>>>3];
		minZeroCount=buckets;
	}

	VariableCompressedDynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[(buckets+7)>>>3];
		minZeroCount=buckets;
	}

	@Override
	public VariableCompressedDynamicLogLog4 copy(){return new VariableCompressedDynamicLogLog4(buckets, k, -1, minProb);}

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

	/** Extract tier part from nibble.
	 *  Standard range (0-11): tierPart = nibble>>>1 (yields 0-5).
	 *  Extended range (12-15): tierPart = nibble - 6 (yields 6-9). */
	static int tierPart(int nibble){
		return nibble<12 ? nibble>>>1 : nibble-6;
	}

	/** Extract history bit from nibble. Returns 0 or 1 for standard range (0-11).
	 *  Returns -1 for extended range (12-15, no history stored). */
	static int histBit(int nibble){
		return nibble<12 ? nibble&1 : -1;
	}

	/** Construct nibble from tierPart and history bit.
	 *  For tierPart 0-5: (tierPart&lt;&lt;1)|hist.
	 *  For tierPart 6-9: tierPart+6 (history not stored). */
	static int makeNibble(int tierPart, int hist){
		return tierPart<6 ? (tierPart<<1)|(hist&1) : tierPart+6;
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
			final int nib=readNibble(i);
			final int tp=tierPart(nib);
			final int absTier=tp+globalNLZ;
			if(absTier<0){continue;}
			filledCount++;
			if(absTier<64){
				nlzCounts[absTier+1]++;
				if(tp<6){
					// Standard tiers: emit into packedBuckets for HC/Mean+H/SBS.
					// Mask hist at absTier==0 (SBS validity).
					final int emitHist=(absTier==0) ? 0 : (nib&1);
					packedBuckets[i]=(char)(((absTier+1)<<1)|emitHist);
				}
				// Extended tiers (tp >= 6): contribute to nlzCounts (DLC/Mean)
				// but NOT packedBuckets — they have no history info and applying
				// HSB state-0 correction would corrupt HC/Mean+H/LDLC.
				// packedBuckets[i] stays 0 (treated as empty by history path).
			}
		}
		nlzCounts[0]=buckets-filledCount;

		return new CardStats(packedBuckets, nlzCounts, 0, 1, 0, 0,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
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
		add((VariableCompressedDynamicLogLog4)log);
	}

	public void add(VariableCompressedDynamicLogLog4 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			for(int i=0; i<buckets; i++){
				final int nA=adjustNibble(readNibble(i), globalNLZ, newGlobalNLZ);
				final int nB=adjustNibble(log.readNibble(i), log.globalNLZ, newGlobalNLZ);
				final int tpA=tierPart(nA), tpB=tierPart(nB);
				final int merged;
				if(tpA>tpB){merged=nA;}
				else if(tpB>tpA){merged=nB;}
				else{
					// Same tier: merge history bits if both are standard range
					if(nA<12 && nB<12){
						merged=makeNibble(tpA, (nA&1)|(nB&1));
					}else{
						merged=nA; // extended range, no history to merge
					}
				}
				writeNibble(i, merged);
			}
			globalNLZ=newGlobalNLZ;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int tp=tierPart(readNibble(i));
				if(tp>0){filledBuckets++;}
				if(tp==0 || (!EARLY_PROMOTE && tp==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && globalNLZ<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	/** Adjust a nibble when shifting to a new global floor (used by add/merge).
	 *  Extended nibbles that shift into standard range lose their history (assigned 0). */
	private static int adjustNibble(int nib, int oldFloor, int newFloor){
		if(nib<=1){return nib;}
		final int tp=tierPart(nib);
		final int hist=histBit(nib); // -1 for extended
		final int newTp=Math.max(0, Math.min(tp+(oldFloor-newFloor), 9));
		if(newTp<6){
			return makeNibble(newTp, hist<0 ? 0 : hist);
		}else{
			return newTp+6; // extended
		}
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

		final int bucket=(int)(key&bucketMask);
		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int oldNibble=readNibble(bucket);
		final int oldTierPart=tierPart(oldNibble);

		if(relTier>=0){
			final int newTierPart=Math.min(relTier+1, 9);
			if(newTierPart>oldTierPart){
				final int delta=(oldTierPart>0) ? (newTierPart-oldTierPart) : (relTier+1);
				final int newHist;
				if(newTierPart>=6){
					newHist=0; // extended: no history stored
				}else{
					newHist=(delta==1) ? 1 : 0;
				}
				lastCardinality=-1;
				writeNibble(bucket, makeNibble(newTierPart, newHist));
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
			// delta == -1: set history if bucket is in standard range
			if(newTierPart==oldTierPart-1 && oldNibble<12){
				if((oldNibble&1)==0){lastCardinality=-1; writeNibble(bucket, oldNibble|1);}
			}
			// Extended range (oldNibble >= 12): delta == -1 is a no-op
		}else{
			// relTier == -1: element one tier below floor
			if(oldTierPart<=1 && oldNibble<12){
				if((oldNibble&1)==0){lastCardinality=-1; writeNibble(bucket, oldNibble|1);}
			}
		}
	}

	private void advanceFloor(){
		globalNLZ++;
		final int relaxedTier=Math.max(0, globalNLZ);
		final int minNlz=(3*relaxedTier)/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	/** Decrement all non-phantom nibbles by one tier after floor advance.
	 *  Extended nibble 12 (tierPart 6) demotes to standard tier 5 with
	 *  variant-dependent history assignment. */
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
					if(nib>=13){
						// Extended tiers 7-9 → 6-8: subtract 1
						nib--;
					}else if(nib==12){
						// Extended tier 6 → standard tier 5: assign history per variant
						final int bucketIdx=w*8+b;
						final int demotionHist;
						if(DEMOTION_MODE==0){demotionHist=1;}         // A: set
						else if(DEMOTION_MODE==1){demotionHist=0;}    // B: clear
						else{demotionHist=bucketIdx&1;}               // C: bucket parity
						nib=10+demotionHist; // makeNibble(5, h) = (5<<1)|h
					}else{
						// Standard range (2-11): subtract 2 preserves history
						nib-=2;
					}
					final int tp=tierPart(nib);
					if(EARLY_PROMOTE){
						if(tp==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(tp<=1){newMinZeroCount++;}
					}
				}
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
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	/** Demotion mode for tier 6 → 5 transition in countAndDecrement.
	 *  0=A (hist=1), 1=B (hist=0), 2=C (hist=bucketID&amp;1). */
	public static int DEMOTION_MODE=0;

	public static boolean EARLY_PROMOTE=true;
	public static boolean DISABLE_EEMASK=false;

	public static final String CF_FILE=null; // placeholder until CF table exists
	public static final String SBS_FILE="?cardinalityCorrectionVCDLL4_LC1BitHistSBS.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal Mean CF. Measured via ddlcalibrate cf=f tmcf=1 tmpcf=1 ddls=8k maxmult=8192
	 *  on cluster with production SBS (1M iters). */
	@Override public float terminalMeanCF(){return 0.878809f;}

	/** Terminal Mean+H CF. Measured same run. */
	@Override public float terminalMeanPlusCF(){return 1.034382f;}
	@Override public float hldlcWeight(){return 0.325f;}

}
