package cardinality;

import shared.Tools;

/**
 * BankedCompressedDynamicLogLog3: CDLL3 with per-word banked shared exponents.
 * <p>
 * Combines CDLL3's compressed NLZ tiers (mantissa mode, TIER_SCALE=1.5)
 * with BDLL3's 2-bit bank exponent in bits 30-31 of each packed word.
 * 10 buckets × 3 bits = 30 bits + 2 bank bits = 32 bits per word.
 * <p>
 * Absolute NLZ = globalNLZ + bankExponent + stored.
 * <p>
 * Bank promotion: when a register would overflow (localRelNlz >= 7), if all
 * 10 registers in the word are >= 1, subtract 1 from all and increment bank.
 * Extends effective per-word range from 7 to 10 tiers.
 *
 * @author Brian Bushnell, Chloe
 * @date April 2026
 */
public final class BankedCompressedDynamicLogLog3 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	BankedCompressedDynamicLogLog3(){
		this(2048, 31, -1, 0);
	}

	BankedCompressedDynamicLogLog3(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	BankedCompressedDynamicLogLog3(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*10<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*10-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*10;
		modBuckets=rounded>0 ? rounded : buckets;
		words=modBuckets/10;
		maxArray=new int[words];
		minZeroCount=modBuckets;
		xOverflow=modBuckets*Math.log(2.0*modBuckets)/256.0;
		storedOverflow=new int[64];
	}

	private static int roundToWords(int b){return Math.max(1, (b+5)/10);}

	@Override
	public BankedCompressedDynamicLogLog3 copy(){return new BankedCompressedDynamicLogLog3(modBuckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	private int readBucket(final int i){
		return (maxArray[i/10]>>>((i%10)*3))&0x7;
	}

	private void writeBucket(final int i, final int val){
		final int wordIdx=i/10;
		final int shift=(i%10)*3;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x7<<shift))|((val&0x7)<<shift);
	}

	private int readBank(final int wordIdx){
		return (maxArray[wordIdx]>>>30)&0x3;
	}

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
		add((BankedCompressedDynamicLogLog3)log);
	}

	public void add(BankedCompressedDynamicLogLog3 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			// Save banks before clearing — readBank after writeBank(w,0) returns 0
			final int[] savedBanks=new int[maxArray.length];
			for(int w=0; w<maxArray.length; w++){
				savedBanks[w]=readBank(w);
				writeBank(w, 0);
			}
			for(int i=0; i<modBuckets; i++){
				final int wordIdx=i/10;
				final int bankA=savedBanks[wordIdx];
				final int bankB=log.readBank(wordIdx);
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				final int absA=sA+globalNLZ+bankA;
				final int absB=sB+log.globalNLZ+bankB;
				final int absMax=Math.max(absA, absB);
				final int newStored=(absMax<=newGlobalNLZ ? 0 : Math.min(absMax-newGlobalNLZ, 7));
				writeBucket(i, newStored);
			}
			globalNLZ=newGlobalNLZ;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<modBuckets; i++){
				final int s=readBucket(i);
				if(s>0){filledBuckets++;}
				if(s==0 || (!EARLY_PROMOTE && s==1)){minZeroCount++;}
			}
			while(minZeroCount==0 && globalNLZ<wordlen){
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int rawNlz=Long.numberOfLeadingZeros(key);
		final int mantissa;
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		final int nlz=(2*rawNlz+mantissa)/3;

		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
		final int wordIdx=bucket/10;

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		int bank=readBank(wordIdx);
		int localRelNlz=nlz-globalNLZ-1-bank;

		if(localRelNlz>=7 && bank<3){
			cntBankPromoteAttempt++;
			if(canPromoteBank(wordIdx)){
				cntBankPromoteSuccess++;
				promoteBank(wordIdx);
				bank=readBank(wordIdx);
				localRelNlz=nlz-globalNLZ-1-bank;
			}
		}

		if(localRelNlz+1>7){cntOverflows++;}
		if(IGNORE_OVERFLOW && localRelNlz+1>7){return;}
		final int newStored=Math.min(localRelNlz+1, 7);
		if(newStored<1){return;}
		final int curBank=readBank(wordIdx);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		lastCardinality=-1;
		cntStores++;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? (oldStored==0 && curBank==0) : (localRelNlz>oldRelNlz && oldRelNlz==0 && curBank==0);
		if(shouldDecrement && --minZeroCount<1){
			while(minZeroCount==0 && globalNLZ<wordlen){
				cntGlobalAdvance++;
				if(CORRECT_OVERFLOW){
					int[] topByBank=new int[4];
					for(int i=0; i<modBuckets; i++){
						if(readBucket(i)==7){topByBank[readBank(i/10)]++;}
					}
					for(int b=0; b<4; b++){
						int tier=8+globalNLZ+b;
						if(tier<64){storedOverflow[tier]+=topByBank[b]/2;}
					}
				}
				advanceFloor();
				minZeroCount=countAndDecrement();
			}
		}
	}

	private void advanceFloor(){
		final int nlzShift=((globalNLZ+1)%2==0) ? 1 : 2;
		globalNLZ++;
		eeMask>>>=nlzShift;
	}

	private boolean canPromoteBank(final int wordIdx){
		final int word=maxArray[wordIdx]&0x3FFFFFFF;
		for(int b=0; b<10; b++){
			if(((word>>>(b*3))&0x7)==0){return false;}
		}
		return true;
	}

	private void promoteBank(final int wordIdx){
		int word=maxArray[wordIdx];
		final int oldBank=(word>>>30)&0x3;
		int result=0;
		for(int b=0; b<10; b++){
			final int shift=b*3;
			int stored=(word>>>shift)&0x7;
			stored--;
			result|=(stored<<shift);
		}
		maxArray[wordIdx]=result|((oldBank+1)<<30);
	}

	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<maxArray.length; w++){
			final int bank=readBank(w);
			if(bank>0){
				cntBankAbsorb++;
				writeBank(w, bank-1);
				if(bank-1==0){
					for(int b=0; b<10; b++){
						final int stored=(maxArray[w]>>>(b*3))&0x7;
						if(EARLY_PROMOTE){
							if(stored==0){newMinZeroCount++;}
						}else{
							if(stored==0 || stored==1){newMinZeroCount++;}
						}
					}
				}
			}else{
				cntRegDecrement++;
				int word=maxArray[w];
				if((word&0x3FFFFFFF)==0){
					final int regsInWord=Math.min(10, modBuckets-w*10);
					newMinZeroCount+=regsInWord;
					continue;
				}
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
				maxArray[w]=result;
			}
		}
		return newMinZeroCount;
	}

	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int wordIdx=i/10;
			final int bank=readBank(wordIdx);
			final int stored=readBucket(i);
			final int absNlz=stored+globalNLZ+bank;
			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				filledCount++;
			}
		}
		nlzCounts[0]=modBuckets-filledCount;

		final int[] counts;
		if(CORRECT_OVERFLOW && globalNLZ>=0){
			final int lo=Math.max(7, globalNLZ+1);
			final int hi=Math.min(10+globalNLZ, 63);
			final int[] cumRaw=new int[64];
			cumRaw[63]=nlzCounts[64];
			for(int t=62; t>=0; t--){cumRaw[t]=cumRaw[t+1]+nlzCounts[t+1];}
			final int[] corrCum=cumRaw.clone();
			for(int t=lo; t<=hi; t++){
				final double x=(storedOverflow[t]>0)
					? storedOverflow[t]*OVERFLOW_SCALE : xOverflow*OVERFLOW_SCALE;
				final int addToCum=(int)Math.round(
					(modBuckets-cumRaw[t])*(1.0-Math.exp(-x/modBuckets)));
				corrCum[t]+=addToCum;
			}
			final int maxHi=Math.min(hi+1, 63);
			if(corrCum[hi]>0 && maxHi>hi){
				corrCum[maxHi]=corrCum[hi]/2;
			}
			counts=nlzCounts.clone();
			if(lo>0){counts[lo]=corrCum[lo-1]-corrCum[lo];}
			for(int t=lo; t<maxHi; t++){counts[t+1]=corrCum[t]-corrCum[t+1];}
			counts[maxHi+1]=corrCum[maxHi];
			int corrSum=0;
			for(int t=1; t<66; t++){corrSum+=counts[t];}
			counts[0]=modBuckets-corrSum;
		}else{
			counts=nlzCounts;
		}

		return new CardStats(null, counts, 0, 0, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.5);
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){
		if(globalNLZ>=0){return 1.0;}
		int empty=0;
		for(int w=0; w<maxArray.length; w++){
			if(readBank(w)>0){continue;}
			final int word=maxArray[w];
			for(int b=0; b<10; b++){
				if(((word>>>(b*3))&0x7)==0){empty++;}
			}
		}
		return 1.0-(double)empty/modBuckets;
	}

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
	private final double xOverflow;
	private final int[] storedOverflow;
	private final int modBuckets;
	private final int words;
	private int globalNLZ=-1;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public int getMinZeros(){return globalNLZ+1;}
	@Override public int actualBuckets(){return modBuckets;}
	int[] getNlzSnapshot(){summarize(); return nlzCounts.clone();}
	int readBucketAt(int i){return readBucket(i);}
	int readBankAt(int w){return readBank(w);}
	private int[] nlzCounts;
	private CardStats lastSummarized;

	public long cntStores=0;
	public long cntOverflows=0;
	public long cntBankPromoteAttempt=0;
	public long cntBankPromoteSuccess=0;
	public long cntBankAbsorb=0;
	public long cntRegDecrement=0;
	public long cntGlobalAdvance=0;

	public int[] bankDistribution(){
		int[] d=new int[4];
		for(int w=0; w<words; w++){d[readBank(w)]++;}
		return d;
	}

	public void printDiagnostics(){
		System.err.println("BCDLL3 diag: stores="+cntStores
			+" overflows="+cntOverflows
			+" bankAttempt="+cntBankPromoteAttempt
			+" bankSuccess="+cntBankPromoteSuccess
			+" bankAbsorb="+cntBankAbsorb
			+" regDecrement="+cntRegDecrement
			+" globalAdvance="+cntGlobalAdvance
			+" globalNLZ="+globalNLZ
			+" added="+added);
	}

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

	public static boolean EARLY_PROMOTE=true;
	public static boolean CORRECT_OVERFLOW=false;
	public static double OVERFLOW_SCALE=1.7;
	public static boolean IGNORE_OVERFLOW=false;

	public static final String CF_FILE="?cardinalityCorrectionBCDLL3.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=null;
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		try{
			return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		}catch(RuntimeException e){
			System.err.println("BCDLL3: No CF table found; running without CF.");
			return CF_MATRIX=null;
		}
	}

	/** Measured 2026-04-18 via ddlcalibrate2 tmcf=1 cf=f ddls=512k buckets=2560 maxmult=8192. */
	@Override public float terminalMeanCF(){return 0.880620f;}
	@Override public float terminalMeanPlusCF(){return 1.0f;}

}
