package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * BankedDynamicLogLog5: BDLL3's banked 2-bit exponent system with UDLL6's 2-bit
 * history, packed 6 buckets per 32-bit word at 5 bits each plus a 2-bit bank.
 * <p>
 * Word layout (32 bits):
 * <ul>
 *   <li>Bits 0-4, 5-9, ..., 25-29: six 5-bit buckets
 *   <li>Bits 30-31: 2-bit bank exponent (shared across the word)
 * </ul>
 * Each 5-bit bucket:
 * <ul>
 *   <li>Bits 0-2 (low 3): stored value (0=empty/floor-level, 1-7 = localRelNlz+1)
 *   <li>Bits 3-4 (high 2): 2-bit history pattern (same semantics as UDLL6)
 * </ul>
 * <p>
 * Absolute NLZ = stored + globalNLZ + bank (always, no floor-level special cases).
 * globalNLZ==-1 means nothing seen; >=0 means all buckets have absNlz >= globalNLZ.
 * <p>
 * Bank promotion: when localRelNlz would overflow and all 6 buckets in the
 * word have stored>=1, subtract 1 from each stored (preserving hist bits)
 * and increment the bank. Effective per-word range expands from 7 to 10 tiers.
 * <p>
 * Global promotion absorption: for words with bank>0, decrement bank instead
 * of touching the buckets. This leaves all hist bits intact — they remain
 * valid because absNlz (= stored+globalNLZ+bank) is preserved across the promotion.
 * <p>
 * History update follows UDLL6's carry-shift rule for 2-bit hist:
 *   newHist = ((1&lt;&lt;2) | oldHist) &gt;&gt; delta &amp; 0x3, for delta in [1,3+].
 *
 * @author Brian Bushnell, Chloe
 * @date April 2026
 */
public final class BankedDynamicLogLog5 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	BankedDynamicLogLog5(){this(2048, 31, -1, 0);}

	/** Construct from parsed command-line arguments. */
	BankedDynamicLogLog5(Parser p){
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
	BankedDynamicLogLog5(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*6 <= 0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*6-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*6;
		modBuckets=rounded>0 ? rounded : buckets;
		words=modBuckets/6;
		maxArray=new int[words];
		minZeroCount=modBuckets;
		promoteThreshold=(PROMOTE_FRAC>0 ? (int)(modBuckets*PROMOTE_FRAC) : 0);
	}

	/** Rounds bucket count to the nearest multiple of 6 (minimum 6). */
	private static int roundToWords(int b){return Math.max(1, (b+3)/6);}

	/** Create an independent copy with a fresh seed. */
	@Override
	public BankedDynamicLogLog5 copy(){return new BankedDynamicLogLog5(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read the 3-bit stored value for bucket i (0=empty, 1-7=relNlz+1). */
	private int readStored(final int i){
		return (maxArray[i/6]>>>((i%6)*5))&0x7;
	}

	/** Read the 2-bit history pattern for bucket i. */
	private int readHist(final int i){
		return (maxArray[i/6]>>>((i%6)*5+3))&0x3;
	}

	/** Write stored (0-7) and hist (0-3) into bucket i's 5-bit slot. */
	private void writeBucket(final int i, final int stored, final int hist){
		final int wordIdx=i/6;
		final int shift=(i%6)*5;
		final int nibble=(stored&0x7)|((hist&0x3)<<3);
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x1F<<shift))|(nibble<<shift);
	}

	/** Read the 2-bit bank exponent from bits [31:30] of the given word. */
	private int readBank(final int wordIdx){
		return (maxArray[wordIdx]>>>30)&0x3;
	}

	/** Write the 2-bit bank exponent (0-3) into bits [31:30] of the given word. */
	private void writeBank(final int wordIdx, final int val){
		maxArray[wordIdx]=(maxArray[wordIdx]&0x3FFFFFFF)|((val&0x3)<<30);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Build absolute-NLZ histogram + packedBuckets (with 2-bit hist) for CardStats.
	 * @return CardStats with all estimator values computed
	 */
	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[modBuckets];
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int wordIdx=i/6;
			final int bank=readBank(wordIdx);
			final int stored=readStored(i);
			final int hist=readHist(i);
			final int absNlz=stored+globalNLZ+bank;
			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				final int emitHist=(absNlz==0) ? 0 : hist;
				packedBuckets[i]=(char)(((absNlz+1)<<2)|emitHist);
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;
		lastSummarized=new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
		return lastSummarized;
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
		add((BankedDynamicLogLog5)log);
	}

	/** Merge another BDLL5 into this one. Bank exponents reset to 0 in merged result. */
	public void add(BankedDynamicLogLog5 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			final int[] newHistA=new int[modBuckets];
			final int[] newHistB=new int[modBuckets];
			final int[] newStoredA=new int[modBuckets];
			final int[] newStoredB=new int[modBuckets];
			for(int i=0; i<modBuckets; i++){
				final int wordIdx=i/6;
				final int bankA=readBank(wordIdx);
				final int bankB=log.readBank(wordIdx);
				final int sA=readStored(i), hA=readHist(i);
				final int sB=log.readStored(i), hB=log.readHist(i);
				final int absA=sA+globalNLZ+bankA;
				final int absB=sB+log.globalNLZ+bankB;
				newStoredA[i]=(absA<newGlobalNLZ+1 ? 0 : Math.min(absA-newGlobalNLZ, 7));
				newStoredB[i]=(absB<newGlobalNLZ+1 ? 0 : Math.min(absB-newGlobalNLZ, 7));
				newHistA[i]=hA; newHistB[i]=hB;
			}
			for(int w=0; w<words; w++){writeBank(w, 0);}
			for(int i=0; i<modBuckets; i++){
				final int sA=newStoredA[i], sB=newStoredB[i];
				final int finalStored, finalHist;
				if(sA>sB){finalStored=sA; finalHist=newHistA[i];}
				else if(sB>sA){finalStored=sB; finalHist=newHistB[i];}
				else{finalStored=sA; finalHist=(newHistA[i]|newHistB[i])&0x3;}
				writeBucket(i, finalStored, finalHist);
			}
			globalNLZ=newGlobalNLZ;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
				final int s=readStored(i);
				if(s>0){filledBuckets++;}
				if(s==0 || (!EARLY_PROMOTE && s==1)){minZeroCount++;}
			}
			while(minZeroCount<=promoteThreshold && globalNLZ<wordlen){
				globalNLZ++;
				final int exitThreshold=Math.max(0, globalNLZ-HISTORY_MARGIN);
				eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Hash a value and store it in the appropriate bucket.
	 * Sub-floor observations update the 2-bit history when delta in {1, 2}.
	 */
	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		if(Long.compareUnsigned(key, eeMask)>0){return;}
		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
		final int wordIdx=bucket/6;
		int bank=readBank(wordIdx);
		int localRelNlz=nlz-globalNLZ-1-bank;

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Single-attempt bank promotion (matches BDLL3). Cascading breaks
		// minZeroCount tracking: shouldDecrement requires curBank==0, which
		// is never true after any cascade step (curBank>=1).
		if(localRelNlz>=7 && bank<3){
			if(canPromoteBank(wordIdx)){
				promoteBank(wordIdx);
				bank=readBank(wordIdx);
				localRelNlz=nlz-globalNLZ-1-bank;
			}
		}

		if(localRelNlz+1>7){totalOverflows++;}
		if(IGNORE_OVERFLOW && localRelNlz+1>7){return;}
		final int newStored=Math.min(localRelNlz+1, 7);

		final int oldStored=readStored(bucket);
		final int oldHist=readHist(bucket);

		if(newStored<=0){
			// Below current floor — only hist update possible.
			if(oldStored>0){
				// Real bucket: sub-max from current absNlz.
				final int delta=oldStored-newStored; // = oldStored - localRelNlz - 1; >=1
				if(delta==1 || delta==2){
					final int bit=1<<(2-delta); // delta=1→bit 1, delta=2→bit 0
					final int nh=oldHist|bit;
					if(nh!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldStored, nh);
					}
				}
			}else if(globalNLZ+bank>=0){
				// Floor-level: effective absNlz=(globalNLZ+bank). Sub-max delta
				// from floor's max to insert = -1-localRelNlz.
				//   localRelNlz=-2 → delta=1 → set bit 1.
				//   localRelNlz=-3 → delta=2 → set bit 0.
				//   localRelNlz=-1 → delta=0 (matches max), no change.
				final int delta=-1-localRelNlz;
				if(delta==1 || delta==2){
					final int bit=1<<(2-delta);
					final int nh=oldHist|bit;
					if(nh!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, 0, nh);
					}
				}
			}
			return;
		}

		if(newStored==oldStored){return;}
		if(newStored<oldStored){
			// hist update for sub-max observation: delta in {1, 2} → set bit.
			final int delta=oldStored-newStored;
			if(delta==1 || delta==2){
				final int bit=1<<(2-delta);
				final int nh=oldHist|bit;
				if(nh!=oldHist){
					lastCardinality=-1;
					writeBucket(bucket, oldStored, nh);
				}
			}
			return;
		}

		// newStored > oldStored: advance.
		branch2++;
		lastCardinality=-1;
		final int newHist;
		if(oldStored==0){
			if(globalNLZ+bank>=0){
				// Floor-level transitioning to real: carry-shift hist as if floor's
				// effective max was at absNlz=(globalNLZ+bank). Delta in absNlz
				// from floor's max to new max = newStored.
				newHist=(newStored>=3) ? 0 : (((1<<2)|oldHist)>>>newStored)&0x3;
			}else{
				// Truly first insert (no prior floor-level entry): no carry possible.
				newHist=0;
			}
		}else{
			final int delta=newStored-oldStored;
			// Carry-shift: ((1<<2) | oldHist) >> delta, masked to 2 bits.
			newHist=(delta>=3) ? 0 : (((1<<2)|oldHist)>>>delta)&0x3;
		}
		writeBucket(bucket, newStored, newHist);
		if(oldStored==0){filledBuckets++;}

		// minZeroCount tracking (same logic as BDLL3).
		final int curBank=bank;
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE
			? (oldStored==0 && curBank==0)
			: (localRelNlz>oldRelNlz && oldRelNlz==0 && curBank==0);
		if(shouldDecrement && --minZeroCount<=promoteThreshold){
			while(minZeroCount<=promoteThreshold && globalNLZ<wordlen){
				globalNLZ++;
				final int exitThreshold=Math.max(0, globalNLZ-HISTORY_MARGIN);
				eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/** Returns true if all 6 stored fields in the word are >= 1 (bank promotion eligible). */
	private boolean canPromoteBank(final int wordIdx){
		final int word=maxArray[wordIdx];
		for(int b=0; b<6; b++){
			if(((word>>>(b*5))&0x7)==0){return false;}
		}
		return true;
	}

	/** Subtract 1 from each stored field and increment bank. Preserves hist bits.
	 *  Caller must verify canPromoteBank() and bank < 3. */
	private void promoteBank(final int wordIdx){
		int word=maxArray[wordIdx];
		final int oldBank=(word>>>30)&0x3;
		int result=0;
		for(int b=0; b<6; b++){
			final int shift=b*5;
			final int nibble=(word>>>shift)&0x1F;
			int stored=nibble&0x7;
			final int hist=nibble&0x18; // already in position (bits 3-4)
			stored--; // guaranteed >= 1 by canPromoteBank check
			result|=((stored|hist)<<shift);
		}
		maxArray[wordIdx]=result|((oldBank+1)<<30);
	}

	/**
	 * Decrement all registers after a global floor advance.
	 * Words with bank > 0 absorb the advance by decrementing their bank.
	 * Words at bank == 0 undergo standard register decrement. Hist bits always preserved.
	 * @return New minZeroCount (buckets eligible for next floor advance)
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<words; w++){
			final int bank=readBank(w);
			if(bank>0){
				final int newBank=bank-1;
				writeBank(w, newBank);
				if(newBank==0){
					for(int b=0; b<6; b++){
						final int stored=(maxArray[w]>>>(b*5))&0x7;
						if(EARLY_PROMOTE){
							if(stored==0){newMinZeroCount++;}
						}else{
							if(stored==0 || stored==1){newMinZeroCount++;}
						}
					}
				}
			}else{
				int word=maxArray[w];
				if((word&0x3FFFFFFF)==0){
					final int regsInWord=Math.min(6, modBuckets-w*6);
					newMinZeroCount+=regsInWord;
					continue;
				}
				int result=0;
				for(int b=0; b<6; b++){
					final int shift=b*5;
					final int nibble=(word>>>shift)&0x1F;
					int stored=nibble&0x7;
					final int hist=nibble&0x18; // already in position (bits 3-4)
					if(stored>0){
						stored--;
						if(EARLY_PROMOTE){
							if(stored==0){newMinZeroCount++; filledBuckets--;}
						}else{
							if(stored==1){newMinZeroCount++;}
						}
					}
					result|=((stored|hist)<<shift);
				}
				maxArray[w]=result; // bank stays 0
			}
		}
		return newMinZeroCount;
	}

	/** Number of buckets with stored > 0 at bank=0 level. */
	public int filledBuckets(){return filledBuckets;}

	/**
	 * True occupancy: a bucket is "empty" only when stored==0 AND globalNLZ+bank&lt;0.
	 * Once globalNLZ>=0 (or any word has bank>0), floor-level buckets are informative and
	 * count as occupied. The raw filledBuckets field lags because EARLY_PROMOTE
	 * global promotion decrements it when stored drops 1→0, even though that bucket
	 * becomes a valid floor-level entry at the new floor.
	 */
	public double occupancy(){
		if(globalNLZ>=0){return 1.0;}
		int empty=0;
		for(int w=0; w<words; w++){
			if(readBank(w)>0){continue;}
			final int word=maxArray[w];
			for(int b=0; b<6; b++){
				if(((word>>>(b*5))&0x7)==0){empty++;}
			}
		}
		return 1.0-(double)empty/modBuckets;
	}

	/** Compatibility accessor: returns globalNLZ+1 to match legacy minZeros convention. */
	public int getMinZeros(){return globalNLZ+1;}

	/** Not used; CF correction handled via CF_MATRIX. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/** Compute all estimator values and return as a legacy-format array. */
	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

	/** Return and clear the cached CardStats. Used by calibration drivers. */
	public CardStats consumeLastSummarized(){
		final CardStats s=lastSummarized;
		lastSummarized=null;
		return s;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 5-bit bucket array with 2-bit bank exponent in bits 30-31. */
	private final int[] maxArray;
	/** Actual bucket count (multiple of 6, may be non-power-of-2). */
	private final int modBuckets;
	/** Number of 32-bit words (modBuckets / 6). */
	private final int words;
	/** Minimum fill threshold triggering global floor advance (fraction of modBuckets). */
	private final int promoteThreshold;
	/** Global floor tier: -1 means nothing seen; >=0 means all buckets have absNlz >= globalNLZ. */
	private int globalNLZ=-1;
	/** Buckets at the global floor eligible for next advance. */
	private int minZeroCount;
	/** Count of buckets with stored > 0 at bank=0 level. */
	private int filledBuckets=0;
	/** Early-exit mask: filters hashes below HISTORY_MARGIN tiers of floor. */
	private long eeMask=-1L;
	/** Cached CardStats from last summarize() call. */
	private CardStats lastSummarized;

	/** Total hashes that passed the eeMask check. */
	public long branch1=0;
	/** Total hashes that advanced a bucket (newStored > oldStored). */
	public long branch2=0;
	/** Total hashes that would overflow the 3-bit stored field. */
	public long totalOverflows=0;
	/** Rate of hashes passing eeMask per hash added. */
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	/** Rate of bucket advances per eeMask-passing hash. */
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum NLZ value (64-bit hash). */
	private static final int wordlen=64;
	/** Sub-NLZ window preserved below the floor. For 2-bit history we must keep
	 *  hashes up to 2 below globalNLZ so sub-max hist updates aren't filtered. */
	public static final int HISTORY_MARGIN=2;
	/** If true, stored=0 triggers floor advance (vs stored<=1 when false). */
	public static boolean EARLY_PROMOTE=true;
	/** Fraction of buckets at floor required before triggering global promotion. */
	public static float PROMOTE_FRAC=0.004f;
	/** When true, hashes that would overflow the register are silently ignored. */
	public static boolean IGNORE_OVERFLOW=false;
	/** When true, apply overflow correction to the NLZ histogram in summarize(). */
	public static boolean CORRECT_OVERFLOW=false;

	/** BDLL5-specific CF table, 16k ddls maxmult=4000 pf=0.004 (Apr 15 2026). */
	public static final String CF_FILE="?cardinalityCorrectionBDLL5.tsv.gz";
	/** BDLL5-specific SBS table, 4M iters 128t (Apr 15 2026). */
	public static final String SBS_FILE="?cardinalityCorrectionBDLL5_LC2BitHist.tsv.gz";
	/** Bucket count the CF_MATRIX was generated for. */
	private static int CF_BUCKETS=2048;
	/** Per-cardinality correction factor table, set by initializeCF or setCFMatrix. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);

	/** Load the CF table from the resource file, scaled to the given bucket count. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	/** Set the CF matrix directly (used by CardinalityParser after global load). */
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Measured 16k ddls maxmult=4000 buckets=2048, cf=f (Apr 15 2026). */
	@Override public float terminalMeanCF(){return 0.693755f;}
	/** Measured with bdll5_hsb_v2 HSB (HISTORY_MARGIN=2, cascade, uniform sampling),
	 *  32 ddls maxmult=512 cf=f (Apr 15 2026). 0.6644/(1+0.099) = 0.6046. */
	@Override public float terminalMeanPlusCF(){return 0.6046f;}

}
