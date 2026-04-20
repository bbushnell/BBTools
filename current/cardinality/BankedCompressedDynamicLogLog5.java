package cardinality;

import shared.Tools;

/**
 * BankedCompressedDynamicLogLog5: CDLL5 with per-word 2-bit bank exponents.
 * <p>
 * Packs 6 buckets into each 32-bit int with zero waste:
 * <ul>
 * <li>Bits [29:0]: 6 × 5-bit registers (bits [4:2]=tierPart, bits [1:0]=history)
 * <li>Bits [31:30]: 2-bit bank exponent (0-3) shared by all 6 registers
 * </ul>
 * Absolute tier = (tierPart - 1) + minZeros + bankExponent.
 * When minZeros+bank > 0, stored=0 means absNlz = minZeros+bank-1 (not empty).
 * <p>
 * Tier compression: halfNlz = 2*NLZ + mantissa, tier = halfNlz/3.
 * TIER_SCALE = 1.5 (each tier ≈ 2^1.5 ≈ 2.83× probability ratio).
 * <p>
 * 2-bit history via UDLL6-style carry-shift:
 * on advance by delta: newHist = ((oldHist | HIST_CARRY) >> delta) &amp; HIST_MASK.
 * On sub-floor observation at diff tiers below: hist |= 1 &lt;&lt; (HBITS-diff).
 * HISTORY_MARGIN=2: eeMask relaxed so hashes at relTier in {-2,-1} pass through.
 * <p>
 * Banking reduces overflow via two mechanisms:
 * <ul>
 * <li>Bank promotion: when a register would overflow (relTier+1 > 7) and
 *     bank &lt; 3, if all 6 registers are non-empty, subtract 1 from all
 *     tierParts and increment the bank. Local operation, 6 registers only.
 * <li>Global floor absorption: when minZeros advances, words with bank > 0
 *     decrement their bank instead of touching registers. Preserves history.
 * </ul>
 *
 * @author Eru, Brian Bushnell
 * @date April 2026
 */
public final class BankedCompressedDynamicLogLog5 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	BankedCompressedDynamicLogLog5(){
		this(2048, 31, -1, 0);
	}

	/** Construct from parsed command-line arguments. */
	BankedCompressedDynamicLogLog5(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	/**
	 * Full constructor.
	 * Bucket count is rounded up to the next multiple of 6 (complete words).
	 * @param buckets_ Number of buckets (rounded to next multiple of 6)
	 * @param k_ Hash prefix length
	 * @param seed Random seed (-1 for default)
	 * @param minProb_ Minimum probability threshold
	 */
	BankedCompressedDynamicLogLog5(int buckets_, int k_, long seed, float minProb_){
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

	/** Create an independent copy with a fresh seed. */
	@Override
	public BankedCompressedDynamicLogLog5 copy(){return new BankedCompressedDynamicLogLog5(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read 5-bit register for bucket i: bits [4:2]=tierPart, [1:0]=history. */
	int readBucket(final int i){
		return (packed[i/6]>>>((i%6)*5))&0x1F;
	}

	/** Write tierPart (0-7) and hist (0-3) into bucket i's 5-bit slot. */
	private void writeBucket(final int i, final int tierPart, final int hist){
		final int idx=i/6;
		final int shift=(i%6)*5;
		final int val=((tierPart&0x7)<<2)|(hist&0x3); // pack tp in [4:2], hist in [1:0]
		packed[idx]=(packed[idx]&~(0x1F<<shift))|(val<<shift);
	}

	/** Extract 3-bit tierPart (bits [4:2]) from a 5-bit register value. */
	private static int tierPart(int reg){return (reg>>>2)&0x7;}
	/** Extract 2-bit history (bits [1:0]) from a 5-bit register value. */
	private static int histBits(int reg){return reg&0x3;}

	/** Read 2-bit bank exponent from bits [31:30] of the packed word. */
	private int readBank(final int wordIdx){
		return (packed[wordIdx]>>>30)&0x3;
	}

	/** Write 2-bit bank exponent (0-3) into bits [31:30] of the packed word. */
	private void writeBank(final int wordIdx, final int val){
		packed[wordIdx]=(packed[wordIdx]&0x3FFFFFFF)|((val&0x3)<<30);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Build NLZ histogram and per-bucket packed state for CardStats.
	 * Converts each bucket's (tierPart, bank, hist) to absolute tier,
	 * emits into nlzCounts[absTier+1] and packedBuckets[i].
	 * @return CardStats with all estimator values computed
	 */
	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		if(packedBuckets==null){packedBuckets=new char[modBuckets];}
		else{java.util.Arrays.fill(packedBuckets, (char)0);}

		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
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
			}else if(minZeros+bank>0){ // stored=0 after floor advance: absNlz=minZeros+bank-1
				final int absTier=minZeros+bank-1;
				if(absTier>=0 && absTier<64){
					nlzCounts[absTier+1]++;
					final int emitHist=emitHistMask(absTier, hist);
					packedBuckets[i]=(char)(((absTier+1)<<2)|emitHist);
					filledCount++;
				}
			}
		}
		nlzCounts[0]=modBuckets-filledCount;
		assert(nlzCounts[0]>=0) : "Negative empties: "+nlzCounts[0]+" filled="+filledCount+" modBuckets="+modBuckets;

		return new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.5);
	}

	/**
	 * Mask history bits to those valid at the given absolute tier.
	 * Tier 0: no sub-tiers exist, all hist invalid.
	 * Tier 1: only bit 1 (one tier below) is valid.
	 * Tier 2+: both bits valid.
	 */
	private static int emitHistMask(int absTier, int hist){
		if(absTier==0){return 0;}
		if(absTier==1){return hist&0x2;} // only the -1 bit
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

	/**
	 * Merge another BCDLL5 into this one.
	 * Resets all bank exponents to 0 in the merged result;
	 * each bucket is converted to absolute tier then reframed.
	 */
	public void add(BankedCompressedDynamicLogLog5 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(packed!=log.packed){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			// Save banks before clearing — readBank after writeBank(w,0) returns 0
			final int[] savedBanks=new int[packedLen];
			for(int w=0; w<packedLen; w++){
				savedBanks[w]=readBank(w);
				writeBank(w, 0); // merged result uses bank=0
			}
			for(int i=0; i<modBuckets; i++){
				final int wordIdx=i/6;
				final int bankA=savedBanks[wordIdx]; // this's pre-merge bank
				final int bankB=log.readBank(wordIdx);
				final int regA=readBucket(i), regB=log.readBucket(i);
				final int tpA=tierPart(regA), tpB=tierPart(regB);
				final int absA=(tpA==0 ? -1 : (tpA-1)+minZeros+bankA);
				final int absB=(tpB==0 ? -1 : (tpB-1)+log.minZeros+bankB);
				final int absMax=Math.max(absA, absB);
				if(absMax<newMinZeros){
					writeBucket(i, 0, 0); // below new floor
				}else{
					final int newTp=Math.min(absMax-newMinZeros+1, 7);
					final int hA=histBits(regA), hB=histBits(regB);
					final int mergedHist=(absA==absB) ? (hA|hB) : (absA>absB ? hA : hB);
					writeBucket(i, newTp, mergedHist);
				}
			}
			minZeros=newMinZeros;
			filledBuckets=0; minZeroCount=0;
			for(int i=0; i<modBuckets; i++){// recount after merge (all banks are 0)
				final int tp=tierPart(readBucket(i));
				if(tp>0){filledBuckets++;}
				if(tp==0 || (!EARLY_PROMOTE && tp==1)){minZeroCount++;}
			}
			assert(filledBuckets>=0 && filledBuckets<=modBuckets) : "Bad filledBuckets="+filledBuckets;
			while(minZeroCount==0 && minZeros<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Hash a value and store it in the appropriate bucket.
	 * Pipeline: hash → NLZ → mantissa → halfNlz → compressed tier →
	 * relTier (minus minZeros and bank) → store or update history.
	 */
	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
		if(!DISABLE_EEMASK && Long.compareUnsigned(key, eeMask)>0){return;} // early-exit mask

		final int rawNlz=Long.numberOfLeadingZeros(key);
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		final int mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		final int halfNlz=2*rawNlz+mantissa;
		final int absTier=halfNlz/3; // compressed tier

		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
		final int wordIdx=bucket/6;
		final int bank=readBank(wordIdx);
		final int relTier=absTier-minZeros-bank; // tier relative to this word's floor

		if(relTier<-HISTORY_MARGIN){return;} // too far below floor even for history

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro); // 64-bit micro cardinality sketch
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int oldReg=readBucket(bucket);
		final int oldTierPart=tierPart(oldReg);
		final int oldHist=histBits(oldReg);

		if(relTier>=0){
			// Bank promotion: check overflow BEFORE clamping tierPart
			if(relTier+1>7 && bank<3 && canPromoteBank(wordIdx)){
				promoteBank(wordIdx); // decrement all 6 tierParts, increment bank
				final int newBank=readBank(wordIdx);
				final int newRelTier=absTier-minZeros-newBank; // recompute with new bank
				final int newTp=Math.min(newRelTier+1, 7);
				final int promotedReg=readBucket(bucket); // re-read after promotion
				final int promotedTp=tierPart(promotedReg);
				final int promotedHist=histBits(promotedReg);
				if(newTp<=promotedTp){return;} // no advancement after promotion
				final int delta=(promotedTp>0) ? (newTp-promotedTp) : (newRelTier+1);
				final int carry=(promotedTp>0 || minZeros+newBank>0) ? HIST_CARRY : 0;
				final int newHist=((promotedHist|carry)>>delta)&HIST_MASK; // carry-shift
				lastCardinality=-1;
				writeBucket(bucket, newTp, newHist);
				return;
			}

			final int newTierPart=Math.min(relTier+1, 7); // clamp after promotion check
			if(newTierPart>oldTierPart){ // tier advance
				final int delta=(oldTierPart>0) ? (newTierPart-oldTierPart) : (relTier+1);
				final int carry=(oldTierPart>0 || minZeros+bank>0) ? HIST_CARRY : 0;
				final int newHist=((oldHist|carry)>>delta)&HIST_MASK; // carry-shift
				lastCardinality=-1;
				writeBucket(bucket, newTierPart, newHist);
				if(oldTierPart==0 && bank==0){filledBuckets++;} // only count at true floor

				// Floor advancement check: only when bucket was at the global floor
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
			if(newTierPart<oldTierPart){ // sub-tier observation: set history bit
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
		}else if(relTier==-1){ // one below word's floor: history update only
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
		}else{ // relTier==-2: two below word's floor
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

	/** Returns true if all 6 registers in the word have tierPart >= 1. */
	private boolean canPromoteBank(final int wordIdx){
		final int word=packed[wordIdx]&0x3FFFFFFF; // mask off bank bits
		for(int b=0; b<6; b++){
			final int reg=(word>>>(b*5))&0x1F;
			if(tierPart(reg)==0){return false;}
		}
		return true;
	}

	/** Subtract 1 from all 6 tierParts in a word and increment its bank exponent.
	 *  Caller must verify canPromoteBank() and bank &lt; 3. */
	private void promoteBank(final int wordIdx){
		final int oldBank=readBank(wordIdx);
		assert(oldBank<3) : "bank overflow: bank="+oldBank;
		assert(canPromoteBank(wordIdx)) : "promoteBank called with empty register in word "+wordIdx;
		int word=packed[wordIdx]&0x3FFFFFFF; // mask off bank bits
		int result=0;
		for(int b=0; b<6; b++){
			final int shift=b*5;
			int reg=(word>>>shift)&0x1F;
			reg=reg-4; // subtract 4 = subtract 1 from tierPart (bits [4:2]), hist unchanged
			result|=(reg<<shift);
		}
		packed[wordIdx]=result|((oldBank+1)<<30); // write registers + new bank
	}

	/**
	 * Advance global floor by one tier.
	 * Relaxes eeMask by HISTORY_MARGIN tiers so hashes at relTier in
	 * {-HISTORY_MARGIN, ..., -1} pass through for history updates.
	 */
	private void advanceFloor(){
		minZeros++;
		final int relaxedTier=Math.max(0, minZeros-HISTORY_MARGIN);
		final int minNlz=(3*relaxedTier)/2; // convert compressed tier back to NLZ
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	/**
	 * Decrement all registers after a global floor advance.
	 * Words with bank > 0 absorb the advance by decrementing their bank,
	 * preserving register values and history bits. Only words at bank=0
	 * undergo destructive register decrement.
	 * @return New minZeroCount (buckets eligible for next floor advance)
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<packedLen; w++){
			final int bank=readBank(w);
			if(bank>0){
				writeBank(w, bank-1); // absorb: decrement bank, registers untouched
				if(bank==1){ // bank drops to 0: stored=0 slots become floor-level
					for(int b=0; b<6; b++){
						final int idx=w*6+b;
						if(idx>=modBuckets){break;}
						final int tp=tierPart(readBucket(idx));
						if(EARLY_PROMOTE){
							if(tp==0){newMinZeroCount++;}
						}else{
							if(tp==0 || tp==1){newMinZeroCount++;}
						}
					}
				}
			}else{
				int word=packed[w];
				if((word&0x3FFFFFFF)==0){ // all registers empty, fast path
					final int regsInWord=Math.min(6, modBuckets-w*6);
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
				packed[w]=result; // bank stays 0
			}
		}
		assert(newMinZeroCount>=0) : "Negative minZeroCount: "+newMinZeroCount;
		return newMinZeroCount;
	}

	/** Number of buckets with tierPart > 0 at bank=0 level. */
	public int filledBuckets(){return filledBuckets;}

	/**
	 * Fraction of buckets that are occupied.
	 * After minZeros > 0, all buckets hold valid data (occupancy=1.0).
	 * Before that, only buckets with tierPart > 0 in bank=0 words are occupied.
	 */
	public double occupancy(){
		if(minZeros>0){return 1.0;} // floor advanced: all buckets have data
		int empty=0;
		for(int w=0; w<packedLen; w++){
			if(readBank(w)>0){continue;} // banked words have no truly empty slots
			for(int b=0; b<6; b++){
				final int idx=w*6+b;
				if(idx>=modBuckets){break;}
				if(tierPart(readBucket(idx))==0){empty++;}
			}
		}
		return 1.0-(double)empty/modBuckets;
	}

	/** Not used; CF correction handled via CF_MATRIX. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Compute all estimator values and return as a legacy-format array.
	 * Caches the CardStats for retrieval by consumeLastSummarized().
	 */
	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		lastSummarized=s;
		final double hybridEst=s.hybridDLL();
		double[] legacy=AbstractCardStats.buildLegacyArray(s, hybridEst);
		double[] ext=new double[legacy.length+2]; // +2 reserved for future estimators
		System.arraycopy(legacy, 0, ext, 0, legacy.length);
		return ext;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed storage: 6 × 5-bit registers + 2-bit bank per 32-bit int. */
	private final int[] packed;
	/** Number of packed ints: ceil(modBuckets/6). */
	private final int packedLen;
	/** Actual bucket count (multiple of 6, may differ from super.buckets). */
	private final int modBuckets;
	/** Global floor tier. absNlz = stored-1 + minZeros + bank. */
	private int minZeros=0;
	/** Buckets at the global floor eligible for next advance. */
	private int minZeroCount;
	/** Count of buckets with tierPart > 0 at bank=0 level. */
	private int filledBuckets=0;
	/** Early-exit mask: filters hashes below HISTORY_MARGIN tiers of floor. */
	private long eeMask=-1L;

	/** Accessor for minZeros. */
	public int getMinZeros(){return minZeros;}
	/** Reusable NLZ histogram buffer (index 0=empty, 1-65=absNlz 0-64). */
	private int[] nlzCounts;
	/** Reusable per-bucket packed state for CardStats: (absTier+1)<<2 | hist. */
	private char[] packedBuckets;
	/** Cached CardStats from last rawEstimates() call. */
	private CardStats lastSummarized;

	/** Return and clear the cached CardStats. Used by calibration drivers. */
	public CardStats consumeLastSummarized(){
		final CardStats cs=lastSummarized;
		lastSummarized=null;
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum NLZ value (64-bit hash). */
	private static final int wordlen=64;
	/** Number of history bits per bucket. */
	private static final int HBITS=2;
	/** Bitmask for history: (1<<HBITS)-1 = 0x3. */
	private static final int HIST_MASK=(1<<HBITS)-1;
	/** Carry bit injected on tier advance: 1<<HBITS = 4. */
	private static final int HIST_CARRY=1<<HBITS;
	/** eeMask relaxation: accept hashes this many tiers below floor for history. */
	static final int HISTORY_MARGIN=2;

	/** Mantissa threshold for compressed tiers: (2-sqrt(2)) * 1048576 ≈ 614242. */
	private static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	/** If true, stored=0 triggers floor advance (vs stored<=1 when false). */
	public static boolean EARLY_PROMOTE=true;
	/** Debug flag: bypass eeMask filtering. Not for production. */
	public static boolean DISABLE_EEMASK=false;

	/** Auto-loaded v5 CF table resource path. */
	public static final String CF_FILE="?cardinalityCorrectionBCDLL5.tsv.gz";
	/** SBS (per-fill LC correction) table for banked compressed 2-bit history. */
	public static final String SBS_FILE="?cardinalityCorrectionBCDLL5_LC2BitHistSBS.tsv.gz";
	/** HSB (per-tier state-bias) table for banked compressed 2-bit history. */
	public static final String HSB_FILE="?cardinalityCorrectionBCDLL5_LC2BitHist.tsv.gz";
	/** Bucket count the CF_MATRIX was generated for. */
	private static int CF_BUCKETS=2048;
	/** Per-cardinality correction factor table, set by initializeCF or setCFMatrix. */
	private static float[][] CF_MATRIX=null;

	/** Load the CF table from the resource file, scaled to the given bucket count. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("BCDLL5: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}

	/** Set the CF matrix directly (used by CardinalityParser after global load). */
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Terminal Mean CF: 1/(raw Mean/true) at saturation.
	 *  Measured 2026-04-18: 512k DDLs, 1536 buckets, maxmult=8192, bug-fixed banks. */
	@Override public float terminalMeanCF(){return 0.880801f;}
	/** Terminal Mean+H CF: 1/(raw Mean+H/true) at saturation.
	 *  Measured 2026-04-18: same run as terminalMeanCF. */
	@Override public float terminalMeanPlusCF(){return 0.909594f;}

}
