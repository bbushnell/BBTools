package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * BankedDynamicLogLog5: BDLL3's banked 3-bit exponent system with UDLL6's 2-bit
 * history, packed 6 buckets per 32-bit word at 5 bits each plus a 2-bit bank.
 * <p>
 * Word layout (32 bits):
 * <ul>
 *   <li>Bits 0-4, 5-9, ..., 25-29: six 5-bit buckets
 *   <li>Bits 30-31: 2-bit bank exponent (shared across the word)
 * </ul>
 * Each 5-bit bucket:
 * <ul>
 *   <li>Bits 0-2 (low 3): stored value (0=empty/phantom, 1-7 = localRelNlz+1)
 *   <li>Bits 3-4 (high 2): 2-bit history pattern (same semantics as UDLL6)
 * </ul>
 * <p>
 * Absolute NLZ = (stored-1) + minZeros + bank when stored>0;
 * phantom absNlz = minZeros + bank - 1 when stored==0 and minZeros+bank>0.
 * <p>
 * Bank promotion: when localRelNlz would overflow and all 6 buckets in the
 * word have stored>=1, subtract 1 from each stored (preserving hist bits)
 * and increment the bank. Effective per-word range expands from 7 to 10 tiers.
 * <p>
 * Global promotion absorption: for words with bank>0, decrement bank instead
 * of touching the buckets. This leaves all hist bits intact — they remain
 * valid because absNlz is preserved across the promotion.
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

	BankedDynamicLogLog5(){this(2048, 31, -1, 0);}

	BankedDynamicLogLog5(Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

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

	private static int roundToWords(int b){return Math.max(1, (b+3)/6);}

	@Override
	public BankedDynamicLogLog5 copy(){return new BankedDynamicLogLog5(modBuckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the 3-bit stored value for bucket i (0=empty, 1-7=relNlz+1). */
	private int readStored(final int i){
		return (maxArray[i/6]>>>((i%6)*5))&0x7;
	}

	/** Reads the 2-bit hist pattern for bucket i. */
	private int readHist(final int i){
		return (maxArray[i/6]>>>((i%6)*5+3))&0x3;
	}

	/** Writes stored + hist together for bucket i. Both values masked. */
	private void writeBucket(final int i, final int stored, final int hist){
		final int wordIdx=i/6;
		final int shift=(i%6)*5;
		final int nibble=(stored&0x7)|((hist&0x3)<<3);
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x1F<<shift))|(nibble<<shift);
	}

	/** Reads the 2-bit bank exponent for the given word. */
	private int readBank(final int wordIdx){
		return (maxArray[wordIdx]>>>30)&0x3;
	}

	/** Writes the 2-bit bank exponent for the given word. */
	private void writeBank(final int wordIdx, final int val){
		maxArray[wordIdx]=(maxArray[wordIdx]&0x3FFFFFFF)|((val&0x3)<<30);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

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

	/** Merges another BDLL5 into this one. Bank exponents reset to 0. */
	public void add(BankedDynamicLogLog5 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
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
				final int absA=(sA==0 ? -1 : (sA-1)+minZeros+bankA);
				final int absB=(sB==0 ? -1 : (sB-1)+log.minZeros+bankB);
				newStoredA[i]=(absA<newMinZeros ? 0 : Math.min(absA-newMinZeros+1, 7));
				newStoredB[i]=(absB<newMinZeros ? 0 : Math.min(absB-newMinZeros+1, 7));
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
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
				final int s=readStored(i);
				if(s>0){filledBuckets++;}
				if(s==0 || (!EARLY_PROMOTE && s==1)){minZeroCount++;}
			}
			while(minZeroCount<=promoteThreshold && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

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
		int localRelNlz=nlz-minZeros-bank;

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Try bank promotion if localRelNlz would overflow
		if(localRelNlz>=7 && bank<3){
			if(canPromoteBank(wordIdx)){
				promoteBank(wordIdx);
				bank=readBank(wordIdx);
				localRelNlz=nlz-minZeros-bank;
			}
		}

		if(localRelNlz+1>7){totalOverflows++;}
		if(IGNORE_OVERFLOW && localRelNlz+1>7){return;}
		final int newStored=Math.min(localRelNlz+1, 7);

		final int oldStored=readStored(bucket);
		final int oldHist=readHist(bucket);

		if(newStored<=0){
			// Below current floor — only hist update possible.
			// localRelNlz is in {-1, -2, ...}. Update hist if oldStored>0
			// and newStored is within hist window of oldStored.
			if(oldStored>0){
				final int delta=oldStored-newStored; // = oldStored - localRelNlz - 1; >=1
				if(delta==1 || delta==2){
					final int bit=1<<(2-delta); // delta=1→bit 1, delta=2→bit 0
					final int nh=oldHist|bit;
					if(nh!=oldHist){
						lastCardinality=-1;
						writeBucket(bucket, oldStored, nh);
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
			// First insert: no prior bit exists — hist must be 0, not carry-shifted
			// from a fictional old top.
			newHist=0;
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
			while(minZeroCount<=promoteThreshold && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/** Returns true if all 6 stored fields in the word are >= 1. */
	private boolean canPromoteBank(final int wordIdx){
		final int word=maxArray[wordIdx];
		for(int b=0; b<6; b++){
			if(((word>>>(b*5))&0x7)==0){return false;}
		}
		return true;
	}

	/** Subtracts 1 from each stored field and increments bank. Preserves hist. */
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

	/** Global tier promotion. Bank>0 words decrement bank; bank==0 words
	 *  decrement each non-empty stored by 1. Hist bits always preserved. */
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
					newMinZeroCount+=6;
					continue;
				}
				int result=0;
				for(int b=0; b<6; b++){
					final int shift=b*5;
					final int nibble=(word>>>shift)&0x1F;
					int stored=nibble&0x7;
					final int hist=nibble&0x18;
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

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/modBuckets;}
	public int getMinZeros(){return minZeros;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/** Build absolute-NLZ histogram + packedBuckets (with 2-bit hist) for CardStats. */
	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[modBuckets];
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int wordIdx=i/6;
			final int bank=readBank(wordIdx);
			final int stored=readStored(i);
			final int hist=readHist(i);
			if(stored>0){
				final int absNlz=(stored-1)+minZeros+bank;
				if(absNlz<64){
					nlzCounts[absNlz+1]++;
					final int emitHist=(absNlz==0) ? 0 : hist;
					packedBuckets[i]=(char)(((absNlz+1)<<2)|emitHist);
					filledCount++;
				}
			}else if(minZeros+bank>0){
				final int absNlz=minZeros+bank-1;
				if(absNlz<64){
					nlzCounts[absNlz+1]++;
					final int emitHist=(absNlz==0) ? 0 : hist;
					packedBuckets[i]=(char)(((absNlz+1)<<2)|emitHist);
					filledCount++;
				}
			}
		}
		nlzCounts[0]=modBuckets-filledCount;
		lastSummarized=new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
		return lastSummarized;
	}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

	public CardStats consumeLastSummarized(){
		final CardStats s=lastSummarized;
		lastSummarized=null;
		return s;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] maxArray;
	private final int modBuckets;
	private final int words;
	private final int promoteThreshold;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;
	private CardStats lastSummarized;

	public long branch1=0, branch2=0;
	public long totalOverflows=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	public static boolean EARLY_PROMOTE=true;
	public static float PROMOTE_FRAC=0.004f;
	/** COF mode default: ignore overflow = false, correct overflow = false. */
	public static boolean IGNORE_OVERFLOW=false;
	public static boolean CORRECT_OVERFLOW=false;

	/** Placeholder: reuses DLL4's CF table until a BDLL5-specific one is generated. */
	public static final String CF_FILE="?cardinalityCorrectionDLL4.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Placeholder terminal Mean CF — will be measured in Phase 2. */
	@Override public float terminalMeanCF(){return 0.72f;}
	@Override public float terminalMeanPlusCF(){return 0.86f;}

}
