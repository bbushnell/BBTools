package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * BankedDynamicLogLog3: DLL3 variant with per-word banked shared exponents.
 * <p>
 * Like DLL3, packs 10 buckets into each int, 3 bits per bucket (30 bits).
 * The 2 unused bits (30-31) become a local bank exponent per word.
 * <p>
 * Absolute NLZ = globalMinZeros + bankExponent + (stored - 1).
 * <p>
 * Bank promotion: when a register would overflow (localRelNlz >= 7), if all
 * 10 registers in the word are >= 1, subtract 1 from all and increment bank.
 * This localizes overflow handling — only 10 registers shift, not the whole
 * structure. Extends effective per-word range from 7 to 10 tiers.
 * <p>
 * Global promotion absorption: when globalMinZeros advances, words with
 * bankExponent > 0 just decrement their bank instead of touching registers.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class BankedDynamicLogLog3 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	BankedDynamicLogLog3(){
		this(2048, 31, -1, 0);
	}

	BankedDynamicLogLog3(Parser p){
		super(p);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		xOverflow=buckets*Math.log(2.0*buckets)/256.0;
		overflowExpFactor=Math.exp(-xOverflow/buckets);
		storedOverflow=new int[64];
		if(FAST_COUNT){nlzCounts=new int[64];}
	}

	BankedDynamicLogLog3(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		xOverflow=buckets*Math.log(2.0*buckets)/256.0;
		overflowExpFactor=Math.exp(-xOverflow/buckets);
		storedOverflow=new int[64];
		if(FAST_COUNT){nlzCounts=new int[64];}
	}

	@Override
	public BankedDynamicLogLog3 copy(){return new BankedDynamicLogLog3(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the 3-bit stored value for bucket i (0=empty, 1-7=relNlz+1). */
	private int readBucket(final int i){
		return (maxArray[i/10]>>>((i%10)*3))&0x7;
	}

	/** Writes a 3-bit stored value for bucket i. val must be in [0,7]. */
	private void writeBucket(final int i, final int val){
		final int wordIdx=i/10;
		final int shift=(i%10)*3;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x7<<shift))|((val&0x7)<<shift);
	}

	/** Reads the 2-bit bank exponent for the word containing bucket i. */
	private int readBank(final int wordIdx){
		return (maxArray[wordIdx]>>>30)&0x3;
	}

	/** Writes the 2-bit bank exponent for a word. val must be in [0,3]. */
	private void writeBank(final int wordIdx, final int val){
		maxArray[wordIdx]=(maxArray[wordIdx]&0x3FFFFFFF)|((val&0x3)<<30);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardinalityStats s=summarize();
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
		add((BankedDynamicLogLog3)log);
	}

	/**
	 * Merges another BDLL3 into this one.
	 * Bank exponents reset to 0 in the merged result for simplicity.
	 * Each bucket is converted to absolute NLZ then reframed.
	 */
	public void add(BankedDynamicLogLog3 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			// Clear all bank exponents — merge resets to bank=0
			for(int w=0; w<maxArray.length; w++){
				writeBank(w, 0);
			}
			for(int i=0; i<buckets; i++){
				final int wordA=i/10, wordB=i/10;
				final int bankA=readBank(wordA);
				final int bankB=log.readBank(wordB);
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to absolute, then to new relative frame (bank=0)
				final int absA=(sA==0 ? -1 : (sA-1)+minZeros+bankA);
				final int absB=(sB==0 ? -1 : (sB-1)+log.minZeros+bankB);
				final int absMax=Math.max(absA, absB);
				final int newStored=(absMax<newMinZeros ? 0 : Math.min(absMax-newMinZeros+1, 7));
				writeBucket(i, newStored);
			}
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
			if(FAST_COUNT){java.util.Arrays.fill(nlzCounts, 0);}
			final int phantomNlz=minZeros-1;
			for(int i=0; i<buckets; i++){
				final int s=readBucket(i);
				if(s>0){
					filledBuckets++;
					if(FAST_COUNT){final int absNlz=(s-1)+minZeros; if(absNlz<64){nlzCounts[absNlz]++;}}
				}else if(FAST_COUNT && minZeros>0 && phantomNlz<64){
					nlzCounts[phantomNlz]++;
				}
				if(s==0 || (!EARLY_PROMOTE && s==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && minZeros<wordlen){
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
		final int bucket=(int)(key&bucketMask);
		final int wordIdx=bucket/10;
		int bank=readBank(wordIdx);
		int localRelNlz=nlz-minZeros-bank;

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Try bank promotion if overflow would occur
		while(localRelNlz>=7 && bank<3 && canPromoteBank(wordIdx)){
			promoteBank(wordIdx);
			localRelNlz--; // bank increased by 1, so localRelNlz decreases by 1
			bank++;
			// After promotion: bank is now bank+1, so localRelNlz = nlz - minZeros - (bank+1) = old - 1
			int localRelNlz2=nlz-minZeros-readBank(wordIdx); // recompute to be safe
			int bank2=readBank(wordIdx);
			assert(bank==bank2);
			assert(localRelNlz==localRelNlz2);
		}

		final int newStored=Math.min(localRelNlz+1, 7);
		if(newStored<1){return;} // below current floor
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		if(FAST_COUNT){
			// Remove old slot from nlzCounts
			if(oldStored>0){
				final int oldAbsNlz=(oldStored-1)+minZeros+bank;
				if(oldAbsNlz<64){nlzCounts[oldAbsNlz]--;}
			}else if(minZeros+bank>0){
				final int pNlz=minZeros+bank-1;
				if(pNlz<64){nlzCounts[pNlz]--;}
			}
			// Add new slot — use current bank (may have been promoted)
			final int curBank=readBank(wordIdx);
			final int newAbsNlz=(newStored-1)+minZeros+curBank;
			if(newAbsNlz<64){nlzCounts[newAbsNlz]++;}
		}

		// minZeroCount tracking
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? oldStored==0 : (localRelNlz>oldRelNlz && oldRelNlz==0);
		if(shouldDecrement && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
				if(USE_STORED_OVERFLOW){
					final int nextCorrTier=7+minZeros;
					if(nextCorrTier<64){
						int topCount=0;
						for(int i=0; i<buckets; i++){
							if(readBucket(i)==7){topCount++;}
						}
						storedOverflow[nextCorrTier]=topCount/2;
					}
				}
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/** Returns true if all 10 registers in the word are non-empty (stored >= 1). */
	private boolean canPromoteBank(final int wordIdx){
		final int word=maxArray[wordIdx]&0x3FFFFFFF; // mask off bank bits
		for(int b=0; b<10; b++){
			if(((word>>>(b*3))&0x7)==0){return false;}
		}
		return true;
	}

	/** Subtracts 1 from all 10 registers in a word and increments its bank exponent. */
	private void promoteBank(final int wordIdx){
		int word=maxArray[wordIdx];
		final int oldBank=(word>>>30)&0x3;
		int result=0;
		for(int b=0; b<10; b++){
			final int shift=b*3;
			int stored=(word>>>shift)&0x7;
			stored--; // guaranteed >= 1 by canPromoteBank check
			result|=(stored<<shift);
		}
		// Write registers + new bank
		maxArray[wordIdx]=result|((oldBank+1)<<30);
	}

	/**
	 * Global tier promotion. For words with bank > 0, absorb by decrementing bank.
	 * For words with bank == 0, decrement registers as in DLL3.
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<maxArray.length; w++){
			final int bank=readBank(w);
			if(bank>0){
				// Absorb global promotion: decrement bank, leave registers untouched
				writeBank(w, bank-1);
				// Count empties/tier-0 in this word (unchanged since registers didn't move)
				for(int b=0; b<10; b++){
					final int stored=(maxArray[w]>>>(b*3))&0x7;
					if(EARLY_PROMOTE){
						if(stored==0){newMinZeroCount++;}
					}else{
						if(stored==0 || stored==1){newMinZeroCount++;}
					}
				}
			}else{
				// bank==0: standard DLL3 decrement
				int word=maxArray[w];
				if((word&0x3FFFFFFF)==0){continue;} // all registers empty
				int result=0;
				for(int b=0; b<10; b++){
					final int shift=b*3;
					int stored=(word>>>shift)&0x7;
					if(stored>0){
						stored--;
						if(EARLY_PROMOTE){
							if(stored==0){newMinZeroCount++; filledBuckets--;}
						}else{
							if(stored==1){newMinZeroCount++;}
						}
					}
					result|=(stored<<shift);
				}
				maxArray[w]=result; // bank stays 0
			}
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Builds absolute NLZ histogram including bank exponents, delegates to CardinalityStats.
	 */
	private CardinalityStats summarize(){
		if(!FAST_COUNT){
			if(nlzCounts==null){nlzCounts=new int[64];}
			else{java.util.Arrays.fill(nlzCounts, 0);}
			for(int i=0; i<buckets; i++){
				final int wordIdx=i/10;
				final int bank=readBank(wordIdx);
				final int stored=readBucket(i);
				if(stored>0){
					final int absNlz=(stored-1)+minZeros+bank;
					if(absNlz<64){nlzCounts[absNlz]++;}
				}else{
					final int phantomNlz=minZeros+bank-1;
					if(phantomNlz>=0 && phantomNlz<64){nlzCounts[phantomNlz]++;}
				}
			}
		}
		final boolean doDebug=DEBUG_ONCE;
		if(doDebug){
			DEBUG_ONCE=false;
			System.err.println("DEBUG summarize(): minZeros="+minZeros+" filledBuckets="+filledBuckets+" buckets="+buckets);
			System.err.println("  nlzCounts (nonzero):");
			for(int t=0; t<64; t++){if(nlzCounts[t]>0){System.err.println("    tier "+t+": "+nlzCounts[t]);}}
		}
		final int[] counts;
		if(CORRECT_OVERFLOW && minZeros>=1){
			final int lo=Math.max(7, minZeros);
			final int hi=Math.min(6+minZeros, 63);
			final int[] cumRaw=new int[64];
			cumRaw[63]=nlzCounts[63];
			for(int t=62; t>=0; t--){cumRaw[t]=cumRaw[t+1]+nlzCounts[t];}
			final int[] corrCum=cumRaw.clone();
			for(int t=lo; t<=hi; t++){
				final double x=(USE_STORED_OVERFLOW && storedOverflow[t]>0)
					? storedOverflow[t]*OVERFLOW_SCALE : xOverflow*OVERFLOW_SCALE;
				final int addToCum=(int)Math.round(
					(buckets-cumRaw[t])*(1.0-Math.exp(-x/buckets)));
				corrCum[t]+=addToCum;
			}
			final int maxHi=Math.min(hi+1, 63);
			if(corrCum[hi]>0 && maxHi>hi){
				corrCum[maxHi]=corrCum[hi]/2;
			}
			counts=nlzCounts.clone();
			if(lo>0){counts[lo-1]=corrCum[lo-1]-corrCum[lo];}
			for(int t=lo; t<maxHi; t++){counts[t]=corrCum[t]-corrCum[t+1];}
			counts[maxHi]=corrCum[maxHi];
		}else{
			counts=nlzCounts;
		}
		lastRawNlz=nlzCounts.clone();
		lastCorrNlz=(counts==nlzCounts) ? lastRawNlz : counts;
		return CardinalityStats.fromNlzCounts(counts, buckets, microIndex,
		                                      CF_MATRIX, CF_BUCKETS,
		                                      CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys);
	}

	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(Math.max(s.hybridDLL(), s.microCardinality()));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 3-bit bucket array with 2-bit bank exponent in bits 30-31. */
	private final int[] maxArray;
	private final double xOverflow;
	private final double overflowExpFactor;
	private final int[] storedOverflow;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return minZeros;}
	public int[] getStoredOverflow(){return storedOverflow;}
	int[] lastRawNlz, lastCorrNlz;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	public static final boolean FAST_COUNT=false;
	public static boolean EARLY_PROMOTE=true;
	public static boolean CORRECT_OVERFLOW=true;
	public static double OVERFLOW_SCALE=1.0;
	public static boolean USE_STORED_OVERFLOW=true;
	public static boolean DEBUG_ONCE=false;

	/** Use DLL3's CF table initially — will need its own table eventually. */
	public static final String CF_FILE="?cardinalityCorrectionBDLL3_cof.tsv.gz";
	private static int CF_BUCKETS=8192;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

}
