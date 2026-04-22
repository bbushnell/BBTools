package cardinality;

import shared.Tools;

/**
 * CompressedDynamicLogLog5: 5-bit packed cardinality estimator with
 * mantissa-compressed tiers and 2-bit history.
 * <p>
 * Bits per bucket: 5 (6 per 32-bit int, 2 bits wasted)
 *   - Bits [4:2]: tier part (0=floor-level, 1-7=relTier 0-6)
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

	/** Default constructor: 2048 buckets, k=31. */
	CompressedDynamicLogLog5(){
		this(2048, 31, -1, 0);
	}

	/** Construct from parsed command-line arguments. */
	CompressedDynamicLogLog5(parse.Parser p){
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

	/** Create an independent copy with a fresh seed. */
	@Override
	public CompressedDynamicLogLog5 copy(){return new CompressedDynamicLogLog5(modBuckets, k, -1, minProb);}
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
		final int val=((tierPart&0x7)<<2)|(hist&0x3);
		packed[idx]=(packed[idx]&~(0x1F<<shift))|(val<<shift);
	}

	/** Extract 3-bit tierPart (bits [4:2]) from a 5-bit register value. */
	private static int tierPart(int reg){return (reg>>>2)&0x7;}
	/** Extract 2-bit history (bits [1:0]) from a 5-bit register value. */
	private static int histBits(int reg){return reg&0x3;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Build NLZ histogram and per-bucket packed state for CardStats.
	 * Converts each bucket's (tierPart, hist) to absolute tier,
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
			final int reg=readBucket(i);
			final int tp=tierPart(reg);
			final int hist=histBits(reg);
			final int absTier=tp+globalNLZ;
			if(absTier>=0 && absTier<64){
				nlzCounts[absTier+1]++;
				final int emitHist=emitHistMask(absTier, hist);
				packedBuckets[i]=(char)(((absTier+1)<<2)|emitHist);
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

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

	/**
	 * Merge another CDLL5 into this one.
	 * Each bucket is reframed to the new common globalNLZ; ties merge history with OR.
	 */
	public void add(CompressedDynamicLogLog5 log){
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
				else{merged=(tpA<<2)|(histBits(rA)|histBits(rB));}
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

	/** Reframe a 5-bit register from oldFloor to newFloor; returns adjusted register. */
	private static int adjustBucket(int reg, int oldFloor, int newFloor){
		final int tp=tierPart(reg);
		if(tp==0){return reg;}
		final int hist=histBits(reg);
		final int newTp=Math.max(0, Math.min(tp+(oldFloor-newFloor), 7));
		return (newTp<<2)|hist;
	}

	/**
	 * Hash a value and store it in the appropriate bucket.
	 * Pipeline: hash → NLZ → mantissa → halfNlz → compressed tier →
	 * relTier (minus globalNLZ+1) → store or update history.
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
		final int absTier=halfNlz/3;
		final int relTier=absTier-(globalNLZ+1);
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
				// Tier advance. Unified delta: floor-level bucket's effective tier is globalNLZ-1,
				// so delta = absTier - (globalNLZ-1) = relTier+1 when oldTierPart==0.
				final int delta=(oldTierPart>0) ? (newTierPart-oldTierPart) : (relTier+1);
				// Carry-shift matches MantissaCompare2 MODE_CTLL (lines 213-215):
				// simulator uses carry=0 when stored==0 (no prior observation).
				// Equivalent here: suppress carry for a truly-fresh bucket
				// (oldTierPart==0 AND globalNLZ<0). Floor-level buckets from floor advance
				// legitimately carry an observation at globalNLZ-1, so carry stays.
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
			// Observation at globalNLZ-1. For oldTierPart>0 (real bucket at
			// absTier = oldTierPart-1+globalNLZ), diff = oldTierPart.
			// For oldTierPart=0 (floor-level at effective globalNLZ-1), diff=0 — no-op.
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
			// relTier == -2: observation at globalNLZ-2.
			// For oldTierPart>0: diff = oldTierPart+1.
			// For oldTierPart=0 (floor-level at globalNLZ-1): diff = 1.
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

	/**
	 * Advance global floor by one tier.
	 * Relaxes eeMask by HISTORY_MARGIN tiers so hashes at relTier in
	 * {-HISTORY_MARGIN, ..., -1} pass through for history updates.
	 */
	private void advanceFloor(){
		globalNLZ++;
		final int relaxedTier=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		final int minNlz=(3*relaxedTier)/2;
		eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
	}

	/**
	 * Decrement all registers after a global floor advance.
	 * Floor-level buckets (tp==0) are left unchanged; history bits preserved.
	 * @return New minZeroCount (buckets eligible for next floor advance)
	 */
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
				// tp==0: already floor-level, leave hist intact
				result|=(reg<<shift);
			}
			packed[w]=result;
		}
		return newMinZeroCount;
	}

	/** Number of buckets with tierPart > 0. */
	public int filledBuckets(){return filledBuckets;}

	/** Fraction of buckets that are occupied (tierPart > 0). */
	public double occupancy(){return (double)filledBuckets/modBuckets;}

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
		double[] ext=new double[legacy.length+2];
		System.arraycopy(legacy, 0, ext, 0, legacy.length);
		return ext;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed storage: 6 × 5-bit registers per 32-bit int (2 bits wasted). */
	private final int[] packed;
	/** Number of packed ints: ceil(modBuckets/6). */
	private final int packedLen;
	/** Actual bucket count (multiple of 6, may differ from super.buckets). */
	private final int modBuckets;
	/** Global floor tier: -1 means nothing seen; >=0 means all buckets have absNlz >= globalNLZ. absNlz = tierPart + globalNLZ. */
	private int globalNLZ=-1;
	/** Buckets at the global floor eligible for next advance. */
	private int minZeroCount;
	/** Count of buckets with tierPart > 0. */
	private int filledBuckets=0;
	/** Early-exit mask: filters hashes below HISTORY_MARGIN tiers of floor. */
	private long eeMask=-1L;

	/** Compatibility accessor: returns globalNLZ+1 to match legacy minZeros convention. */
	public int getMinZeros(){return globalNLZ+1;}
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
	public static final String CF_FILE="?cardinalityCorrectionCDLL5.tsv.gz";
	/** SBS (per-fill LC correction) table, 2-bit history with half-NLZ tier geometry. */
	public static final String SBS_FILE="?cardinalityCorrectionCDLL5_LC2BitHistSBS.tsv.gz";
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
			System.err.println("CDLL5: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}

	/** Set the CF matrix directly (used by CardinalityParser after global load). */
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Measured 2026-04-20: 128k DDLs, 1536 modulo buckets, maxmult=8192, tmcf=1 tmpcf=1. */
	@Override public float terminalMeanCF(){return 0.883080f;}
	@Override public float terminalMeanPlusCF(){return 1.094268f;}

	/** HC weight for LDLC blend. Calibrated by Eru, 2026-04-22. */
	@Override public double ldlcHcWeight(){return 0.40;}

}
