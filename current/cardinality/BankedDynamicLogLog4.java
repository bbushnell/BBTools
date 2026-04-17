package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * BankedDynamicLogLog4: BDLL3's banked exponent system with 1-bit history,
 * packed 7 buckets per 32-bit word at 4 bits each plus a 4-bit bank.
 * <p>
 * Word layout (32 bits):
 * <ul>
 *   <li>Bits 0-3, 4-7, ..., 24-27: seven 4-bit buckets
 *   <li>Bits 28-31: 4-bit bank exponent (shared across the word)
 * </ul>
 * Each 4-bit bucket:
 * <ul>
 *   <li>Bits 0-2 (low 3): stored value (0=empty/phantom, 1-7 = localRelNlz+1)
 *   <li>Bit 3 (high 1): 1-bit history pattern
 * </ul>
 * <p>
 * Absolute NLZ = (stored-1) + minZeros + bank when stored>0;
 * phantom absNlz = minZeros + bank - 1 when stored==0 and minZeros+bank>0.
 * <p>
 * Bank promotion: when localRelNlz would overflow and all 7 buckets in the
 * word have stored>=1, subtract 1 from each stored (preserving hist bits)
 * and increment the bank. 4-bit bank gives range 0-15, extending effective
 * per-word range from 7 to 22 tiers.
 *
 * @author Brian Bushnell, Eru
 * @date April 2026
 */
public final class BankedDynamicLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	BankedDynamicLogLog4(){this(2048, 31, -1, 0);}

	BankedDynamicLogLog4(Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	BankedDynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*7 <= 0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*7-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*7;
		modBuckets=rounded>0 ? rounded : buckets;
		words=modBuckets/7;
		maxArray=new int[words];
		minZeroCount=modBuckets;
		promoteThreshold=(PROMOTE_FRAC>0 ? (int)(modBuckets*PROMOTE_FRAC) : 0);
	}

	private static int roundToWords(int b){return Math.max(1, (b+6)/7);}

	@Override
	public BankedDynamicLogLog4 copy(){return new BankedDynamicLogLog4(modBuckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	private int readStored(final int i){
		return (maxArray[i/7]>>>((i%7)*4))&0x7;
	}

	private int readHist(final int i){
		return (maxArray[i/7]>>>((i%7)*4+3))&0x1;
	}

	private void writeBucket(final int i, final int stored, final int hist){
		final int wordIdx=i/7;
		final int shift=(i%7)*4;
		final int nibble=(stored&0x7)|((hist&0x1)<<3);
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0xF<<shift))|(nibble<<shift);
	}

	private int readBank(final int wordIdx){
		return (maxArray[wordIdx]>>>28)&0xF;
	}

	private void writeBank(final int wordIdx, final int val){
		maxArray[wordIdx]=(maxArray[wordIdx]&0x0FFFFFFF)|((val&0xF)<<28);
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
		add((BankedDynamicLogLog4)log);
	}

	public void add(BankedDynamicLogLog4 log){
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
				final int wordIdx=i/7;
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
				else{finalStored=sA; finalHist=(newHistA[i]|newHistB[i])&0x1;}
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
				final int exitThreshold=Math.max(0, minZeros-HISTORY_MARGIN);
				eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
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
		final int wordIdx=bucket/7;
		int bank=readBank(wordIdx);
		int localRelNlz=nlz-minZeros-bank;

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		if(localRelNlz>=7 && bank<15){
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
			// Sub-floor observation. delta_abs from real bucket (oldStored-1)+minZeros+bank
			// is oldStored-newStored; from phantom at minZeros+bank-1 is -newStored.
			// 1-bit hist captures only delta_abs==1 events.
			if(oldStored>0){
				final int delta=oldStored-newStored;
				if(delta==1 && oldHist==0){
					lastCardinality=-1;
					writeBucket(bucket, oldStored, 1);
				}
			}else if(minZeros+bank>0 && newStored==-1 && oldHist==0){
				// Phantom sub-floor: delta_abs = -newStored = 1
				lastCardinality=-1;
				writeBucket(bucket, 0, 1);
			}
			return;
		}

		if(newStored==oldStored){return;}
		if(newStored<oldStored){
			final int delta=oldStored-newStored;
			if(delta==1){
				if(oldHist==0){
					lastCardinality=-1;
					writeBucket(bucket, oldStored, 1);
				}
			}
			return;
		}

		// newStored > oldStored: advance. Follows CDLL4 pattern:
		//   delta = (oldStored>0) ? (newStored-oldStored) : newStored
		// For oldStored==0, delta encodes observation vs phantom/implicit tier at
		// minZeros+bank-1. Preserves hist=1 on phantom→real delta_abs==1 transitions.
		// 1-bit hist: newHist=(delta==1)?1:0.
		branch2++;
		lastCardinality=-1;
		final int delta=(oldStored>0) ? (newStored-oldStored) : newStored;
		final int newHist=(delta==1) ? 1 : 0;
		writeBucket(bucket, newStored, newHist);
		if(oldStored==0){filledBuckets++;}

		final int curBank=bank;
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE
			? (oldStored==0 && curBank==0)
			: (localRelNlz>oldRelNlz && oldRelNlz==0 && curBank==0);
		if(shouldDecrement && --minZeroCount<=promoteThreshold){
			while(minZeroCount<=promoteThreshold && minZeros<wordlen){
				minZeros++;
				final int exitThreshold=Math.max(0, minZeros-HISTORY_MARGIN);
				eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
				minZeroCount=countAndDecrement();
			}
		}
	}

	private boolean canPromoteBank(final int wordIdx){
		final int word=maxArray[wordIdx];
		for(int b=0; b<7; b++){
			if(((word>>>(b*4))&0x7)==0){return false;}
		}
		return true;
	}

	private void promoteBank(final int wordIdx){
		int word=maxArray[wordIdx];
		final int oldBank=(word>>>28)&0xF;
		int result=0;
		for(int b=0; b<7; b++){
			final int shift=b*4;
			final int nibble=(word>>>shift)&0xF;
			int stored=nibble&0x7;
			final int hist=nibble&0x8;
			stored--;
			result|=((stored|hist)<<shift);
		}
		maxArray[wordIdx]=result|((oldBank+1)<<28);
	}

	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<words; w++){
			final int bank=readBank(w);
			if(bank>0){
				final int newBank=bank-1;
				writeBank(w, newBank);
				if(newBank==0){
					for(int b=0; b<7; b++){
						final int stored=(maxArray[w]>>>(b*4))&0x7;
						if(EARLY_PROMOTE){
							if(stored==0){newMinZeroCount++;}
						}else{
							if(stored==0 || stored==1){newMinZeroCount++;}
						}
					}
				}
			}else{
				int word=maxArray[w];
				if((word&0x0FFFFFFF)==0){
					newMinZeroCount+=7;
					continue;
				}
				int result=0;
				for(int b=0; b<7; b++){
					final int shift=b*4;
					final int nibble=(word>>>shift)&0xF;
					int stored=nibble&0x7;
					final int hist=nibble&0x8;
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
				maxArray[w]=result;
			}
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}

	/** True occupancy: a bucket is "empty" only when stored==0 AND minZeros+bank==0.
	 *  Once minZeros>0 (or any word has bank>0), phantom buckets are informative and
	 *  count as occupied. The raw filledBuckets field lags because EARLY_PROMOTE
	 *  global promotion decrements it when stored drops 1→0, even though that bucket
	 *  becomes a valid phantom at the new floor. */
	public double occupancy(){
		if(minZeros>0){return 1.0;}
		int empty=0;
		for(int w=0; w<words; w++){
			if(readBank(w)>0){continue;} // bank>0: all 7 buckets are phantoms, occupied
			final int word=maxArray[w];
			for(int b=0; b<7; b++){
				if(((word>>>(b*4))&0x7)==0){empty++;}
			}
		}
		return 1.0-(double)empty/modBuckets;
	}

	public int getMinZeros(){return minZeros;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[modBuckets];
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int wordIdx=i/7;
			final int bank=readBank(wordIdx);
			final int stored=readStored(i);
			final int hist=readHist(i);
			if(stored>0){
				final int absNlz=(stored-1)+minZeros+bank;
				if(absNlz<64){
					nlzCounts[absNlz+1]++;
					final int emitHist=(absNlz==0) ? 0 : hist;
					packedBuckets[i]=(char)(((absNlz+1)<<1)|emitHist);
				}
				filledCount++;
			}else if(minZeros+bank>0){
				final int absNlz=minZeros+bank-1;
				if(absNlz<64){
					nlzCounts[absNlz+1]++;
					final int emitHist=(absNlz==0) ? 0 : hist;
					packedBuckets[i]=(char)(((absNlz+1)<<1)|emitHist);
				}
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;
		lastSummarized=new CardStats(packedBuckets, nlzCounts, 0, 1, 0, 0,
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
	/** Conservative eeMask margin: lag by this many tiers so sub-floor hashes
	 *  reach the hot path to populate history bits. Equals HIST_BITS (1 here). */
	private static final int HISTORY_MARGIN=1;
	public static boolean EARLY_PROMOTE=true;
	public static float PROMOTE_FRAC=0.004f;
	public static boolean IGNORE_OVERFLOW=false;
	public static boolean CORRECT_OVERFLOW=false;

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

	/** Measured: raw Mean at saturation ≈ 1.4415 × trueCard (16k ddls, maxmult=4000, co=f). */
	@Override public float terminalMeanCF(){return 0.693722f;}
	@Override public float terminalMeanPlusCF(){return 0.86f;}

}
