package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicLogLog3v3: 3-bit packed DLL with DLL3-compatible encoding and social promotion.
 * <p>
 * Based on DLL3v2, but uses DLL3's countAndDecrement logic (decrement all stored>0).
 * With EARLY_PROMOTE=false and PROMOTE_FRAC=0, produces identical results to DLL3 ep=f.
 * Social promotion (PROMOTE_FRAC>0) allows early floor advancement for reduced overflow.
 * <p>
 * Packs 10 buckets into each int, 3 bits per bucket.
 * Encoding: 0 = empty/phantom; 1-7 = (relNlz + 1).
 * After floor advancement, stored=0 means "phantom" (one tier below floor).
 *
 * @author Brian Bushnell, Chloe, Eru
 * @date April 2026
 */
public final class DynamicLogLog3v3 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DynamicLogLog3v3(){
		this(2048, 31, -1, 0);
	}

	DynamicLogLog3v3(Parser p){
		super(p);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		promoteThreshold=(PROMOTE_FRAC>0 ? (int)(buckets*PROMOTE_FRAC) : PROMOTE_THRESHOLD);
		xOverflow=buckets*Math.log(2.0*buckets)/256.0;
		overflowExpFactor=Math.exp(-xOverflow/buckets);
		storedOverflow=new int[64];
	}

	DynamicLogLog3v3(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		promoteThreshold=(PROMOTE_FRAC>0 ? (int)(buckets*PROMOTE_FRAC) : PROMOTE_THRESHOLD);
		xOverflow=buckets*Math.log(2.0*buckets)/256.0;
		overflowExpFactor=Math.exp(-xOverflow/buckets);
		storedOverflow=new int[64];
	}

	@Override
	public DynamicLogLog3v3 copy(){return new DynamicLogLog3v3(buckets, k, -1, minProb);}

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

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		double rawHyb=s.hybridDLL();
		if(IGNORE_OVERFLOW && USE_IO_BIAS){rawHyb=ioBiasCorrect(rawHyb, buckets);}
		long card=(long)(rawHyb);
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicLogLog3v3)log);
	}

	/**
	 * Merges another DynamicLogLog3v3 into this one.
	 */
	public void add(DynamicLogLog3v3 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				final int nA=(sA==0 ? 0 : Math.max(0, Math.min(sA+(minZeros-newMinZeros), 7)));
				final int nB=(sB==0 ? 0 : Math.max(0, Math.min(sB+(log.minZeros-newMinZeros), 7)));
				writeBucket(i, Math.max(nA, nB));
			}
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int s=readBucket(i);
				if(s>0){filledBuckets++;}
				if(EARLY_PROMOTE){
					if(s==0){minZeroCount++;}
				}else{
					if(s==0||s==1){minZeroCount++;}
				}
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
		final int bucket=(int)(key&bucketMask);
		final int relNlz=nlz-minZeros;

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Stored = relNlz+1, clamped to [1,7] for overflow
		if(IGNORE_OVERFLOW && relNlz+1>7){return;} // silently ignore overflow
		final int newStored=Math.min(relNlz+1, 7);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		// EARLY_PROMOTE=false (classic): tracks empty+tier-0; advances when all buckets >= 2.
		// EARLY_PROMOTE=true  (new):     tracks empty only;    advances when all buckets >= 1.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? oldStored==0 : (relNlz>oldRelNlz && oldRelNlz==0);
		if(shouldDecrement && --minZeroCount<=promoteThreshold){
			while(minZeroCount<=promoteThreshold && minZeros<wordlen){
				if(USE_STORED_OVERFLOW){
					final int nextCorrTier=7+minZeros;
					if(nextCorrTier<64){
						int topCount=0;
						for(int i=0; i<buckets; i++){if(readBucket(i)==7){topCount++;}}
						storedOverflow[nextCorrTier]=topCount/2;
					}
				}
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Decrements all non-empty stored values by 1.
	 * Uses DLL3-compatible logic: stored=1 becomes stored=0 (phantom).
	 * Returns count of tracked buckets after decrement:
	 *   EARLY_PROMOTE=true:  count of new empties (stored=1→0)
	 *   EARLY_PROMOTE=false: count of new tier-0  (stored=2→1)
	 * Pre-existing stored=0 (phantoms from social promotion) are correctly skipped.
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<maxArray.length; w++){
			int word=maxArray[w];
			if(word==0){continue;}
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
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public int minZeros(){return minZeros;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Scans the bucket array once to populate the absolute NLZ histogram (nlzCounts),
	 * then applies Poisson overflow correction if CORRECT_OVERFLOW is enabled.
	 * Phantom buckets (stored=0 when minZeros>0) are treated as absNlz = minZeros-1.
	 */
	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}

		final int phantomNlz=minZeros-1;
		int filledCount=0;
		for(int i=0; i<buckets; i++){
			final int stored=readBucket(i);
			if(stored>0){
				final int absNlz=(stored-1)+minZeros;
				if(absNlz<64){nlzCounts[absNlz+1]++;}
				filledCount++;
			}else if(minZeros>0 && phantomNlz>=0 && phantomNlz<64){
				nlzCounts[phantomNlz+1]++;
				filledCount++;
			}
		}
		nlzCounts[0]=buckets-filledCount;
		final int[] counts;
		if(CORRECT_OVERFLOW && minZeros>=1){
			final int lo=Math.max(7, minZeros);
			final int hi=Math.min(6+minZeros, 63);
			final int[] cumRaw=new int[64];
			cumRaw[63]=nlzCounts[64];
			for(int t=62; t>=0; t--){cumRaw[t]=cumRaw[t+1]+nlzCounts[t+1];}
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
			if(lo>0){counts[lo]=corrCum[lo-1]-corrCum[lo];}
			for(int t=lo; t<maxHi; t++){counts[t+1]=corrCum[t]-corrCum[t+1];}
			counts[maxHi+1]=corrCum[maxHi];
			int corrSum=0;
			for(int t=1; t<66; t++){corrSum+=counts[t];}
			counts[0]=buckets-corrSum;
		}else{
			counts=nlzCounts;
		}
		return new CardStats(null, counts, 0, 0, 0, 0, buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0);
	}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		double[] r=AbstractCardStats.buildLegacyArray(s, hybridEst);
		if(IGNORE_OVERFLOW && USE_IO_BIAS){
			final double mult=ioBiasMult(hybridEst, buckets);
			for(int i=0; i<r.length; i++){r[i]*=mult;}
		}
		return r;
	}

	public CardStats lastStats;

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 3-bit bucket array: 10 buckets per int, 0=empty, 1-7=relNlz+1. */
	private final int[] maxArray;
	/** Per-instance promotion threshold (computed from PROMOTE_FRAC*buckets or PROMOTE_THRESHOLD). */
	private final int promoteThreshold;
	private final double xOverflow;
	private final double overflowExpFactor;
	private final int[] storedOverflow;
	private int minZeros=0;
	/** Count of tracked buckets; triggers minZeros floor advance when <=promoteThreshold. */
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	public static int PROMOTE_THRESHOLD=0;
	public static float PROMOTE_FRAC=0.02f;
	public static boolean EARLY_PROMOTE=false;
	public static boolean CORRECT_OVERFLOW=true;
	public static double OVERFLOW_SCALE=1.7;
	public static boolean USE_STORED_OVERFLOW=true;
	public static boolean IGNORE_OVERFLOW=false;
	public static boolean USE_IO_BIAS=true;

	/** Per-tier bias constants for IO mode. Indexed by (tier-IO_TIER_MIN)/IO_TIER_STEP.
	 *  Can be overridden at runtime via iobfile= flag. */
	private static double[] IO_BIAS={
		-0.015172, -0.016103, -0.015661, -0.014902, -0.014328, // tiers -8..-4
		-0.013666, -0.016286, -0.017821, -0.020524, -0.028997, // tiers -3..1
		-0.045063, -0.070147, -0.095604, -0.115183, -0.128107, // tiers 2..6
		-0.132166, -0.129972, -0.126304, -0.126153, -0.126579, // tiers 7..11
		-0.126523                                                // tier 12
	};
	static double IO_TIER_STEP=1.0;
	static final int IO_TIER_MIN=-8, IO_TIER_MAX=12;

	public static void loadIOBias(String path){
		try{
			java.util.ArrayList<Double> vals=new java.util.ArrayList<>();
			java.io.BufferedReader br=new java.io.BufferedReader(new java.io.FileReader(path));
			for(String line=br.readLine(); line!=null; line=br.readLine()){
				line=line.trim();
				if(line.isEmpty() || line.startsWith("#")){continue;}
				String[] parts=line.split("\\s+");
				vals.add(Double.parseDouble(parts[0]));
			}
			br.close();
			final int fullCount=IO_TIER_MAX-IO_TIER_MIN+1; // 21
			final int halfCount=(IO_TIER_MAX-IO_TIER_MIN)*2+1; // 41
			if(vals.size()==fullCount){
				IO_TIER_STEP=1.0;
			}else if(vals.size()==halfCount){
				IO_TIER_STEP=0.5;
			}else{
				throw new RuntimeException("IO_BIAS file "+path+" has "+vals.size()+" values, expected "+fullCount+" or "+halfCount);
			}
			IO_BIAS=new double[vals.size()];
			for(int i=0; i<vals.size(); i++){IO_BIAS[i]=vals.get(i);}
		}catch(java.io.IOException e){throw new RuntimeException("Failed to load IO_BIAS from "+path, e);}
	}

	public static void setIOBias(double[] bias){IO_BIAS=bias;}

	static double ioBiasCorrect(double raw, int buckets){
		return raw*ioBiasMult(raw, buckets);
	}

	static double ioBiasMult(double raw, int buckets){
		if(raw<=0){return 1;}
		double ftier=Math.log(raw/buckets)/Math.log(2);
		double fidx=(ftier-IO_TIER_MIN)/IO_TIER_STEP;
		int lo=(int)Math.floor(fidx);
		lo=Math.max(0, Math.min(IO_BIAS.length-1, lo));
		int hi=Math.min(lo+1, IO_BIAS.length-1);
		double frac=fidx-Math.floor(fidx);
		double bias=IO_BIAS[lo]+(IO_BIAS[hi]-IO_BIAS[lo])*frac;
		return 1.0/(1.0+bias);
	}

	/** Default resource file for DLL3v3 correction factors (uses DLL3's table initially). */
	public static final String CF_FILE="?cardinalityCorrectionDLL3.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		CF_MATRIX_CARD=CorrectionFactor.lastCardMatrix;
		CF_CARD_KEYS=CorrectionFactor.lastCardKeys;
		return CF_MATRIX;
	}

	static float[][] CF_MATRIX_CARD;
	static float[] CF_CARD_KEYS;
	private int[] nlzCounts;
}
