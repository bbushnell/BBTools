package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicLogLog3v2: 3-bit packed variant of DynamicLogLog using relative NLZ storage.
 * <p>
 * Packs 10 buckets into each int, 3 bits per bucket.
 * Encoding: 0 = empty; 1-7 = (relNlz + 1), where relNlz = absoluteNlz - minZeros.
 * No mantissa — coarser precision per bucket, but allows more buckets for the same memory.
 * Overflow (relNlz >= 7) clamped to stored=7 and absorbed by CF matrix.
 * <p>
 * Memory: 3 bits/bucket × 8192 buckets = 24576 bits = 3072 bytes (~3KB).
 * 10 buckets per int; 2 bits per int wasted (bits 30-31 always zero).
 * Paper description: "3 bits per bucket."
 * <p>
 * Key invariants:
 * - maxArray[i/10] holds bucket i at bit position (i%10)*3, 3 bits wide.
 * - Bucket value 0 = empty; value 1-7 = (relNlz+1).
 * - absoluteNlz = (stored - 1) + minZeros for non-empty buckets.
 * - minZeroCount tracks empty + tier-0 (stored=1) buckets; advances minZeros floor when 0.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicLogLog3v2 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DynamicLogLog3v2(){
		this(2048, 31, -1, 0);
	}

	DynamicLogLog3v2(Parser p){
		super(p);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		promoteThreshold=(PROMOTE_FRAC>0 ? (int)(buckets*PROMOTE_FRAC) : PROMOTE_THRESHOLD);
	}

	DynamicLogLog3v2(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[(buckets+9)/10];
		minZeroCount=buckets;
		promoteThreshold=(PROMOTE_FRAC>0 ? (int)(buckets*PROMOTE_FRAC) : PROMOTE_THRESHOLD);
	}

	@Override
	public DynamicLogLog3v2 copy(){return new DynamicLogLog3v2(buckets, k, -1, minProb);}

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
		add((DynamicLogLog3v2)log);
	}

	/**
	 * Merges another DynamicLogLog3v2 into this one.
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
	 * (one DLL3 per thread, merged at the end) nondeterministically reduces accuracy
	 * and this cannot be compensated for.  The cause is architectural: each per-thread
	 * copy sees only a fraction of the data, so its sliding minZeros window promotes
	 * less aggressively than a single estimator seeing all data.  This increases net
	 * overflow events across the merged result.  With DLL3's narrow 7-tier range, the
	 * effect is severe: on a 12M-kmer dataset, t=1 gives ~10.8M, t=2 ~10.0M, t=4 ~8.2M,
	 * t=8 ~6.3M — roughly halving with each doubling of threads.
	 * For parallel use on the same stream, prefer a single synchronized instance.
	 */
	public void add(DynamicLogLog3v2 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to new relative frame: newStored = stored + (oldMinZeros - newMinZeros)
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
				if(s==0||s==1){minZeroCount++;}
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

		// oldRelNlz: 0 for both empty (stored=0) and tier-0 (stored=1).
		// minZeroCount decrements whenever a bucket leaves (empty or tier-0) for a higher tier.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		if(relNlz>oldRelNlz && oldRelNlz==0 && --minZeroCount<=promoteThreshold){
			while(minZeroCount<=promoteThreshold && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
		minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Counts buckets that will become tier-0 after decrement (currently stored=2),
	 * then decrements all non-empty stored values by 1.
	 * Only called when minZeroCount==0, guaranteeing all non-empty have stored>=2.
	 * Returns count of new tier-0 buckets (stored==1 after decrement).
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
				if(stored>1){
					stored--; // decrement relative tier
					if(stored==1){newMinZeroCount++;} // new tier-0 after decrement
				}else if(stored==1){
					if(RESET_ON_PROMOTE){
						stored=0; // reset to empty instead of retaining
						filledBuckets--; // track that this bucket is now empty
					}
					newMinZeroCount++; // empty or tier-0, either way counts toward minZeroCount
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
	 * then delegates all sum computation to CardStats.fromNlzCounts().
	 * <p>
	 * Phantom buckets (stored=0 when minZeros>0) are treated as absNlz = minZeros-1,
	 * one tier below the current floor. This ensures the nlzCounts distribution is
	 * identical regardless of promotion setting.
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
		return new CardStats(null, nlzCounts, 0, 0, 0, 0, buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
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

	/** Last CardStats from rawEstimates(); for debugging. */
	public CardStats lastStats;

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 3-bit bucket array: 10 buckets per int, 0=empty, 1-7=relNlz+1. */
	private final int[] maxArray;
	/** Per-instance promotion threshold (computed from PROMOTE_FRAC*buckets or PROMOTE_THRESHOLD). */
	private final int promoteThreshold;
	private int minZeros=0;
	/** Count of (empty + tier-0) buckets; triggers minZeros floor advance when 0. */
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;
	// sortBuf inherited from CardinalityTracker (lazy, gated by USE_SORTBUF)
	// lastCardinality inherited from CardinalityTracker

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Fallback absolute social promotion threshold when PROMOTE_FRAC==0.
	 * 0 = classic behavior. Shared with DLL3/DLL4 for parser compatibility. */
	public static int PROMOTE_THRESHOLD=0;
	/** Fraction of buckets to use as social promotion threshold.
	 * When >0, overrides PROMOTE_THRESHOLD: promoteThreshold = (int)(buckets * PROMOTE_FRAC).
	 * E.g. 0.05 = promote when last 5% of tier-0 remain (default). */
	public static float PROMOTE_FRAC=0.02f;
	/** When true, social promotion resets tier-0 (stored=1) buckets to empty (stored=0)
	 *  instead of retaining them. Tests whether CF tables converge better without
	 *  retained stale tier-0 values biasing the estimator after floor advancement. */
	public static boolean RESET_ON_PROMOTE=false;
	/** When true, silently drop hashes that would overflow (relNlz >= 7). */
	public static boolean IGNORE_OVERFLOW=false;
	/** When true, apply per-tier IO_BIAS correction in IO mode. */
	public static boolean USE_IO_BIAS=true;

	/** Per-tier bias constants for IO mode. Indexed by (tier-IO_TIER_MIN)/IO_TIER_STEP.
	 *  Can be overridden at runtime via iobfile= flag. */
	private static double[] IO_BIAS={
		-0.007270, -0.006966, -0.007194, -0.007460, -0.007832, // tiers -8..-4
		-0.008336, -0.009006, -0.010032, -0.013286, -0.023293, // tiers -3..1
		-0.045825, -0.069857, -0.090856, -0.111244, -0.128496, // tiers 2..6
		-0.138177, -0.139975, -0.140810, -0.141716, -0.141567, // tiers 7..11
		-0.141567                                                // tier 12
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

	/** Returns IO-bias-corrected estimate: raw / (1 + bias[tier]). */
	static double ioBiasCorrect(double raw, int buckets){
		return raw*ioBiasMult(raw, buckets);
	}

	/** Returns the multiplicative correction factor for IO mode. */
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

	/** Default resource file for DLL3 correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDLL3.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the DLL3v2 correction factor matrix from CF_FILE. Also captures cardinality table. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
		CF_MATRIX_CARD=CorrectionFactor.lastCardMatrix;
		CF_CARD_KEYS=CorrectionFactor.lastCardKeys;
		return CF_MATRIX;
	}

	/** Cardinality-indexed CF table from bipartite CF file; null for legacy files. */
	private static float[][] CF_MATRIX_CARD=null;
	/** MeanEst key array for CF_MATRIX_CARD binary search. */
	private static float[] CF_CARD_KEYS=null;

	/** Stub: measure from preliminary CF table, then replace 1f with actual ratio. */
	@Override public float terminalMeanCF(){return 1f;}

}
