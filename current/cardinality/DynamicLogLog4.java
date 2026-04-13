package cardinality;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicLogLog4: 4-bit packed variant of DynamicLogLog using relative NLZ storage.
 * <p>
 * Packs 8 buckets into each int, 4 bits per bucket (nibble packing).
 * Encoding: 0 = empty; 1-15 = (relNlz + 1), where relNlz = absoluteNlz - minZeros.
 * No mantissa — coarser precision per bucket, but allows 2x bucket count for the same memory.
 * Overflow (relNlz >= 15) clamped to stored=15 and absorbed by CF matrix.
 * <p>
 * Key question: does DLL4 with 4096 buckets beat DLL2 with 2048 buckets for the same
 * 4KB memory footprint? (4096 * 4 bits = 2048 * 8 bits = 2048 chars = 4KB)
 * <p>
 * Key invariants:
 * - maxArray[i>>>3] holds 8 consecutive buckets, each in 4 bits.
 * - Bucket value 0 = empty; value 1-15 = (relNlz+1).
 * - absoluteNlz = (stored - 1) + minZeros for non-empty buckets.
 * - minZeroCount tracks empty + tier-0 (stored=1) buckets; advances minZeros floor when 0.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DynamicLogLog4(){
		this(2048, 31, -1, 0);
	}

	DynamicLogLog4(Parser p){
		super(p);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
		if(FAST_COUNT){nlzCounts=new int[66];}
	}

	DynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
		if(FAST_COUNT){nlzCounts=new int[66];}
	}

	@Override
	public DynamicLogLog4 copy(){return new DynamicLogLog4(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the 4-bit stored value for bucket i (0=empty, 1-15=relNlz+1). */
	private int readBucket(final int i){
		return (maxArray[i>>>3]>>>((i&7)<<2))&0xF;
	}

	/** Writes a 4-bit stored value for bucket i. val must be in [0,15]. */
	private void writeBucket(final int i, final int val){
		final int wordIdx=i>>>3;
		final int shift=(i&7)<<2;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0xF<<shift))|((val&0xF)<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Populates the NLZ histogram (nlzCounts) and builds a CardStats.
	 * <p>
	 * Phantom buckets (stored=0 when minZeros>0) are treated as absNlz = minZeros-1,
	 * one tier below the current floor. This ensures the nlzCounts distribution is
	 * identical regardless of EARLY_PROMOTE setting.
	 */
	private CardStats summarize(){
		if(!FAST_COUNT){
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
		}else{
			// FAST_COUNT=true: nlzCounts already maintained incrementally.
			int sum=0; for(int t=1; t<66; t++){sum+=nlzCounts[t];} nlzCounts[0]=buckets-sum;
		}
		lastRawNlz=nlzCounts.clone();
		int effectiveBuckets=buckets;
		if(CLAMP_OUTLIERS){effectiveBuckets-=clampOutlierTiers(nlzCounts);}
		lastCorrNlz=nlzCounts.clone();
		// DLL4 is counts-only: no history, luck, or mantissa bits. buckets=null.
		return new CardStats(null, nlzCounts, 0, 0, 0, 0,
				effectiveBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
	}

	/**
	 * Clamp outlier tiers in the NLZ histogram.
	 * Computes the mean occupied tier, then sets a ceiling at meanTier + log2(filled).
	 * Buckets above the ceiling are moved down to the ceiling tier.
	 * Total bucket count is preserved.
	 */
	/** Move outlier tier buckets down to the mean tier. Returns 0 (no bucket count change). */
	private static int clampOutlierTiers(final int[] counts){
		// Compute filled count and mean tier (counts[t+1] = buckets at absNlz=t)
		int filled=0;
		double tierSum=0;
		int maxTier=0;
		for(int t=1; t<counts.length; t++){
			if(counts[t]>0){
				filled+=counts[t];
				tierSum+=counts[t]*(t-1); // absNlz = t-1
				maxTier=t;
			}
		}
		if(filled<2){return 0;} // nothing to clamp
		final double meanTier=tierSum/filled;
		final int ceiling=Math.max(1, (int)(meanTier+Math.log(filled)/Math.log(2))-3); // clamp floor to 1 to protect card=1
		if(ceiling>=maxTier){return 0;} // nothing above ceiling

		// Move outlier buckets above ceiling down to mean tier
		final int meanIdx=(int)Math.round(meanTier)+1; // +1 for index offset
		int excess=0;
		for(int t=ceiling+1; t<counts.length; t++){
			excess+=counts[t];
			counts[t]=0;
		}
		counts[meanIdx]+=excess;
		return 0;
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		final double rawHyb=s.hybridDLL(); // CF already inside blend (on meanEstCF)
		long card=(long)(rawHyb);
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicLogLog4)log);
	}

	/**
	 * Merges another DynamicLogLog4 into this one.
	 * <p>
	 * Correct merge sequence:
	 * <ol>
	 *   <li>Re-frame both instances to newMinZeros = max(this.minZeros, log.minZeros).
	 *       Buckets whose absNlz falls below the new floor correctly become empty (stored=0),
	 *       mirroring what countAndDecrement() would do one step at a time.
	 *       Note: Math.max(0,...) is correct here — NOT Math.max(1,...).
	 *       A below-floor bucket is genuinely empty in the new frame.</li>
	 *   <li>Take per-bucket max of the re-framed values.</li>
	 *   <li>Scan to recount filledBuckets and minZeroCount, then advance the floor
	 *       (possibly multiple times) if minZeroCount==0, exactly as hashAndStore does.</li>
	 * </ol>
	 * <p>
	 * <b>Multi-threaded accuracy warning:</b> Using clonal per-thread estimator copies
	 * (one DLL4 per thread, merged at the end) nondeterministically reduces accuracy
	 * and this cannot be compensated for.  The cause is architectural: each per-thread
	 * copy sees only a fraction of the data, so its sliding minZeros window promotes
	 * less aggressively than a single estimator seeing all data.  This increases net
	 * overflow events across the merged result.  With DLL4's 15-tier range the effect
	 * is minor in practice: on a 12M-kmer dataset, t=1 gives ~11.43M, t=2 ~11.43M,
	 * t=4 ~11.43M, t=8 ~11.42M — well under 0.1% degradation even at 8 threads.
	 * For parallel use on the same stream, prefer a single synchronized instance.
	 */
	public void add(DynamicLogLog4 log){
		added+=log.added;
		branch1+=log.branch1;
		branch2+=log.branch2;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to new relative frame: newStored = stored + (oldMinZeros - newMinZeros)
				final int nA=(sA==0 ? 0 : Math.max(0, Math.min(sA+(minZeros-newMinZeros), 15)));
				final int nB=(sB==0 ? 0 : Math.max(0, Math.min(sB+(log.minZeros-newMinZeros), 15)));
				writeBucket(i, Math.max(nA, nB));
			}
			minZeros=newMinZeros;
			filledBuckets=0;
			minZeroCount=0;
//			if(FAST_COUNT){java.util.Arrays.fill(nlzCounts, 0);}
//			final int phantomNlz=minZeros-1;
			for(int i=0; i<buckets; i++){
				final int s=readBucket(i);
				if(s>0){
					filledBuckets++;
//					if(FAST_COUNT){final int absNlz=(s-1)+minZeros; if(absNlz<64){nlzCounts[absNlz]++;}}
				}
//				else if(FAST_COUNT && minZeros>0 && phantomNlz<64){
//					nlzCounts[phantomNlz]++;
//				}
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
		final int relNlz=nlz-minZeros;

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){//Optional MicroIndex for low cardinality
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS) {return;}//Allows lazy array allocation
		}

		// Stored = relNlz+1, clamped to [1,15] for overflow
		if(IGNORE_OVERFLOW && relNlz+1>15){return;} // silently ignore overflow
		final int newStored=Math.min(relNlz+1, 15);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		branch2++;
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		if(FAST_COUNT){
			// Remove old slot from nlzCounts
			if(oldStored>0){
				final int oldAbsNlz=(oldStored-1)+minZeros;
				if(oldAbsNlz<64){nlzCounts[oldAbsNlz+1]--;}
			}else if(minZeros>0){
				final int pNlz=minZeros-1;
				if(pNlz<64){nlzCounts[pNlz+1]--;}
			}else{
				nlzCounts[0]--; // was truly empty
			}
			// Add new slot
			final int newAbsNlz=(newStored-1)+minZeros;
			if(newAbsNlz<64){nlzCounts[newAbsNlz+1]++;}
		}

		// minZeroCount decrements when a bucket leaves the tracked-zero category.
		// EARLY_PROMOTE=false (classic): tracks empty+tier-0; advances when all buckets >= 2.
		// EARLY_PROMOTE=true  (new):     tracks empty only;    advances when all buckets >= 1.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? oldStored==0 : (relNlz>oldRelNlz && oldRelNlz==0);
		if(shouldDecrement && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Decrements all buckets with value at least 1.
	 * Counts buckets that became min-tier (value 0 for EARLY_PROMOTE, value 1 for classic).
	 * Only called when minZeroCount==0, so all buckets will be >=1.
	 * Returns count of new minimum-tier buckets.
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int wordIdx=0; wordIdx<maxArray.length; wordIdx++){
			final int oldWord=maxArray[wordIdx];
			if(oldWord==0){continue;}
			int newWord=0;
			for(int bucketIdx=0; bucketIdx<8; bucketIdx++){
				final int bucketOffset=bucketIdx<<2;
				int bucketVal=(oldWord>>>bucketOffset)&0xF;
				if(bucketVal>0){
					bucketVal--;
					if(EARLY_PROMOTE){
						if(bucketVal==0){newMinZeroCount++; filledBuckets--;}
					}else{
						if(bucketVal==1){newMinZeroCount++;}
					}
				}
				newWord|=(bucketVal<<bucketOffset);
			}
			maxArray[wordIdx]=newWord;
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
		final double hybridEst=s.hybridDLL();
		double[] legacy=AbstractCardStats.buildLegacyArray(s, hybridEst);
		// Append word-based estimates if table is loaded
		if(wordTable!=null && WORD_TABLE_TIERS>0){
			double wordEst=rawWordEstimateShuffle()*WORD_TERMINAL_CORRECTION;
			double wordEstCV=rawWordEstimateCV()*WORD_TERMINAL_CORRECTION;
			// Apply CF correction to WordEst using DLC as seed estimate
			wordEst*=wordEstCF(s.dlcRawF, wordEst);
			// Extend array: add WordEst, WordEstCV at the end
			double[] ext=new double[legacy.length+2];
			System.arraycopy(legacy, 0, ext, 0, legacy.length);
			ext[legacy.length]=wordEst;
			ext[legacy.length+1]=wordEstCV;
			return ext;
		}
		// Pad to include WordEst/WordEstCV slots (unused, 0.0) for calibration driver compatibility
		double[] ext=new double[legacy.length+2];
		System.arraycopy(legacy, 0, ext, 0, legacy.length);
		return ext;
	}

	/** Returns the CF correction factor for WordEst.
	 *  Uses the v1/v5 CF table with DLC as the seed estimate for key lookup.
	 *  Returns 1.0 if CF is disabled or no table is loaded. */
	private static double wordEstCF(double dlcSeed, double rawWordEst){
		if(!CorrectionFactor.USE_CORRECTION){return 1;}
		final float[][] mat=CorrectionFactor.v1Matrix;
		final float[] keys=CorrectionFactor.v1Keys;
		if(mat==null || CorrectionFactor.WORDEST>=mat.length){return 1;}
		return CorrectionFactor.getCF(dlcSeed, rawWordEst, mat[CorrectionFactor.WORDEST], keys,
			AbstractCardStats.DEFAULT_CF_ITERS, AbstractCardStats.DEFAULT_CF_DIF);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Word-Based Estimation     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Estimates cardinality using canonical per-word state table lookup.
	 * For each 4-bucket word: pack raw registers into 16-bit key, remap
	 * to canonical index. Tier = minZeros + wordMin. Sums stateAvg.
	 */
	public double rawWordEstimate(){
		if(wordTable==null || WORD_TABLE_TIERS<=0){return 0;}
		final int numWords=buckets>>>2;
		double sum=0;
		for(int w=0; w<numWords; w++){
			final int base=w<<2;
			int r0=readBucket(base), r1=readBucket(base+1);
			int r2=readBucket(base+2), r3=readBucket(base+3);
			if(r0==0 && r1==0 && r2==0 && r3==0 && minZeros==0){continue;}
			int wMin=Math.min(Math.min(r0, r1), Math.min(r2, r3));
			int tier=Math.min(minZeros+wMin, WORD_TABLE_TIERS-1);
			int rawKey=r0|(r1<<4)|(r2<<8)|(r3<<12);
			int idx=WORD_REMAP[rawKey];
			sum+=wordTable[tier][idx];
		}
		return sum;
	}

	/**
	 * CV-weighted word-based estimate. Weights each word's contribution
	 * by inverse CV raised to CV_POWER.
	 */
	public double rawWordEstimateCV(){
		if(wordTable==null || WORD_TABLE_TIERS<=0 || wordCVTable==null){return rawWordEstimate();}
		final int numWords=buckets>>>2;
		double weightedSum=0, weightSum=0;
		int filledWords=0;
		for(int w=0; w<numWords; w++){
			final int base=w<<2;
			int r0=readBucket(base), r1=readBucket(base+1);
			int r2=readBucket(base+2), r3=readBucket(base+3);
			if(r0==0 && r1==0 && r2==0 && r3==0 && minZeros==0){continue;}
			filledWords++;
			int wMin=Math.min(Math.min(r0, r1), Math.min(r2, r3));
			int tier=Math.min(minZeros+wMin, WORD_TABLE_TIERS-1);
			int rawKey=r0|(r1<<4)|(r2<<8)|(r3<<12);
			int idx=WORD_REMAP[rawKey];
			double value=wordTable[tier][idx];
			double cv=wordCVTable[tier][idx];
			double weight=Math.pow(1.0/Math.max(cv, 0.001), CV_POWER);
			weightedSum+=value*weight;
			weightSum+=weight;
		}
		if(weightSum<=0){return 0;}
		return (weightedSum/weightSum)*filledWords;
	}

	/**
	 * Shuffled multi-context word estimate.  Each pass groups all buckets into
	 * words of 4 and sums per-word stateAvg.  Pass 0 uses the natural order;
	 * passes 1..SHUFFLE_PASSES-1 use a Fisher-Yates shuffle with a fixed seed
	 * so every bucket appears in SHUFFLE_PASSES different word contexts.
	 * Final result is the average across all passes.
	 * Reduces variance by ~sqrt(SHUFFLE_PASSES) compared to rawWordEstimate().
	 */
	public double rawWordEstimateShuffle(){
		if(wordTable==null || WORD_TABLE_TIERS<=0){return 0;}
		final int numWords=buckets>>>2;
		final int[] indices=new int[buckets];
		for(int i=0; i<buckets; i++){indices[i]=i;}

		double totalSum=0;
		final java.util.Random shuffleRng=new java.util.Random(SHUFFLE_SEED);
		for(int pass=0; pass<SHUFFLE_PASSES; pass++){
			if(pass>0){
				// Fisher-Yates shuffle
				for(int i=buckets-1; i>0; i--){
					final int j=shuffleRng.nextInt(i+1);
					final int tmp=indices[i]; indices[i]=indices[j]; indices[j]=tmp;
				}
			}
			double passSum=0;
			for(int w=0; w<numWords; w++){
				final int base=w<<2;
				final int r0=readBucket(indices[base]), r1=readBucket(indices[base+1]);
				final int r2=readBucket(indices[base+2]), r3=readBucket(indices[base+3]);
				if(r0==0 && r1==0 && r2==0 && r3==0 && minZeros==0){continue;}
				final int wMin=Math.min(Math.min(r0, r1), Math.min(r2, r3));
				final int tier=Math.min(minZeros+wMin, WORD_TABLE_TIERS-1);
				final int rawKey=r0|(r1<<4)|(r2<<8)|(r3<<12);
				final int idx=WORD_REMAP[rawKey];
				passSum+=wordTable[tier][idx];
			}
			totalSum+=passSum;
		}
		return totalSum/SHUFFLE_PASSES;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 4-bit bucket array: 8 buckets per int, 0=empty, 1-15=relNlz+1. */
	private final int[] maxArray;
	private int minZeros=0;
	/** Count of (empty + tier-0) buckets; triggers minZeros floor advance when 0. */
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	public int getMinZeros(){return minZeros;}
	/** Last raw and corrected nlzCounts from summarize(). */
	int[] lastRawNlz, lastCorrNlz;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** When true, nlzCounts is maintained incrementally in hashAndStore() rather than rebuilt
	 *  in summarize(). Eliminates the O(buckets) scan per rawEstimates() call — ~32x speedup
	 *  for CF table generation. Disabled in production (false) to avoid the per-add overhead. */
	public static final boolean FAST_COUNT=false;
	/** Social promotion threshold (see DynamicLogLog3v2). 0=classic behavior. */
	public static int PROMOTE_THRESHOLD=0;
	/** When true, advance tier as soon as all buckets are nonzero (stored>=1) rather than >=2.
	 * At advance, subtracts 1 from all buckets, which may reset some to empty (stored=0).
	 * This is safe because lcMin (tier-compensated LC) handles post-advance zero buckets correctly.
	 * Reduces tier-15+ overflow pollution; requires new CF table generation when changed. */
	public static boolean EARLY_PROMOTE=true;
	/** When true, silently drop hashes that would overflow (relNlz >= 15).
	 * Eliminates overflow pollution at the cost of a small systematic bias
	 * that the CF table corrects. */
	public static boolean IGNORE_OVERFLOW=false;
	/** When true, clamp outlier tiers in the NLZ histogram before passing to CardStats.
	 * Buckets above meanTier + log2(filledBuckets) are moved down to the ceiling tier.
	 * Reduces overflow pollution and early-lucky-hash distortion. */
	public static boolean CLAMP_OUTLIERS=false;


	/** Default resource file for DDL correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDLL4.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the DDL correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	/*--------------------------------------------------------------*/
	/*----------------   Word State Table Statics   ----------------*/
	/*--------------------------------------------------------------*/

	/** Remap: raw 16-bit register state -> compact canonical index. */
	static final char[] WORD_REMAP=new char[65536];
	/** Number of unique canonical states. */
	static final int WORD_NUM_CANONICAL;
	static {
		structures.IntHashMap map=new structures.IntHashMap(1024);
		for(int state=0; state<65536; state++){
			int r0=state&0xF, r1=(state>>4)&0xF, r2=(state>>8)&0xF, r3=(state>>12)&0xF;
			int wMin=Math.min(Math.min(r0, r1), Math.min(r2, r3));
			int d0=r0-wMin, d1=r1-wMin, d2=r2-wMin, d3=r3-wMin;
			if(d0>d1){int t=d0; d0=d1; d1=t;}
			if(d2>d3){int t=d2; d2=d3; d3=t;}
			if(d0>d2){int t=d0; d0=d2; d2=t;}
			if(d1>d3){int t=d1; d1=d3; d3=t;}
			if(d1>d2){int t=d1; d1=d2; d2=t;}
			int canonKey=d0|(d1<<4)|(d2<<8)|(d3<<12);
			if(!map.containsKey(canonKey)){map.put(canonKey, map.size());}
			WORD_REMAP[state]=(char)map.get(canonKey);
		}
		WORD_NUM_CANONICAL=map.size();
	}

	/** Per-tier word state averages [tier][numCanonical]. null until loaded. */
	static double[][] wordTable;
	/** Per-tier word state CVs [tier][numCanonical]. null until loaded. */
	static double[][] wordCVTable;
	/** Number of tiers in the word table. */
	static int WORD_TABLE_TIERS=0;
	/** Terminal correction multiplier applied to word estimates. */
	public static double WORD_TERMINAL_CORRECTION=1.0;
	/** CV weighting power: weight = (1/CV)^CV_POWER. */
	public static double CV_POWER=2.0;
	/** Number of shuffle passes for rawWordEstimateShuffle(). Pass 0=natural order. */
	public static int SHUFFLE_PASSES=1;
	/** Fixed RNG seed for deterministic shuffle order. */
	public static final long SHUFFLE_SEED=48577L;

	static final String WORD_TABLE_FILE="stateTableDLL4Word.tsv.gz";

	/**
	 * Loads the DLL4 word state table from resources.
	 * Sparse format: #tier \t idx \t canonKey \t stateAvg \t count \t cv
	 * Populates wordTable[tier][idx] = stateAvg.
	 */
	public static synchronized void loadWordTable(){
		if(wordTable!=null){return;}
		String path=dna.Data.findPath("?"+WORD_TABLE_FILE);
		if(path==null){
			System.err.println("WARNING: DLL4 word state table not found: "+WORD_TABLE_FILE);
			return;
		}
		FileFormat ff=FileFormat.testInput(path, null, false);
		if(ff==null){return;}
		ByteFile bf=ByteFile.makeByteFile(ff, 1);
		if(bf==null){return;}

		int maxTier=0;
		java.util.ArrayList<byte[]> allLines=new java.util.ArrayList<>();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			allLines.add(line);
			String s=new String(line);
			int tier=Integer.parseInt(s.substring(0, s.indexOf('\t')));
			if(tier>maxTier){maxTier=tier;}
		}
		bf.close();

		WORD_TABLE_TIERS=maxTier+1;
		wordTable=new double[WORD_TABLE_TIERS][WORD_NUM_CANONICAL];
		wordCVTable=new double[WORD_TABLE_TIERS][WORD_NUM_CANONICAL];
		for(int t=0; t<WORD_TABLE_TIERS; t++){
			java.util.Arrays.fill(wordCVTable[t], 1.0);
		}

		int loaded=0;
		for(byte[] line : allLines){
			String s=new String(line);
			String[] parts=s.split("\t");
			if(parts.length<6){continue;}
			int tier=Integer.parseInt(parts[0]);
			int idx=Integer.parseInt(parts[1]);
			// parts[2] = canonKey (for reference)
			double avg=Double.parseDouble(parts[3]);
			// parts[4] = count
			double cv=Double.parseDouble(parts[5]);
			if(tier<WORD_TABLE_TIERS && idx>=0 && idx<WORD_NUM_CANONICAL){
				wordTable[tier][idx]=avg;
				wordCVTable[tier][idx]=cv;
				loaded++;
			}
		}

		System.err.println("Loaded DLL4 word table: "+WORD_TABLE_TIERS+" tiers, "
			+loaded+" entries, "+WORD_NUM_CANONICAL+" canonical states");
	}

	/** Stub: measure from preliminary CF table, then replace 1f with actual ratio. */
	@Override public float terminalMeanCF(){return 1f;}

}
