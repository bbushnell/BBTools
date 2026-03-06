package cardinality;

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
		this(4096, 31, -1, 0);
	}

	DynamicLogLog4(Parser p){
		super(p);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
	}

	DynamicLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[buckets>>>3];
		minZeroCount=buckets;
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

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final double[] est=rawEstimates();
		final long card=Math.min(added, (long)est[6]); // Hybrid
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
	 * Converts relative stored values to absolute, takes per-bucket max, re-relativizes.
	 */
	public void add(DynamicLogLog4 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			for(int i=0; i<buckets; i++){
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to new relative frame: newStored = stored + (oldMinZeros - newMinZeros)
				final int nA=(sA==0 ? 0 : Math.max(1, Math.min(sA+(minZeros-newMinZeros), 15)));
				final int nB=(sB==0 ? 0 : Math.max(1, Math.min(sB+(log.minZeros-newMinZeros), 15)));
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
		if(relNlz<0){return;} // safety guard (eeMask should prevent this)

		// Stored = relNlz+1, clamped to [1,15] for overflow
		final int newStored=Math.min(relNlz+1, 15);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		branch2++;
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		// oldRelNlz: 0 for both empty (stored=0) and tier-0 (stored=1).
		// minZeroCount decrements whenever a bucket leaves (empty or tier-0) for a higher tier.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		if(relNlz>oldRelNlz && oldRelNlz==0 && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
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
			for(int b=0; b<8; b++){
				final int shift=b<<2;
				int stored=(word>>>shift)&0xF;
				if(stored>0){
					stored--; // decrement relative tier
					if(stored==1){newMinZeroCount++;} // new tier-0 after decrement
				}
				result|=(stored<<shift);
			}
			maxArray[w]=result;
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Returns cardinality estimates for calibration.
	 * Output order: 0=Mean, 1=HMean, 2=HMeanM(=HMean; no mantissa in DLL4),
	 *               3=GMean, 4=HLL, 5=LC, 6=Hybrid, 7=MWA, 8=MedianCorr
	 */
	public double[] rawEstimates(){
		double difSum=0;
		double hllSumFilled=0;
		double gSum=0;
		int count=0;
		sortBuf.clear();

		for(int i=0; i<buckets; i++){
			final int stored=readBucket(i);
			if(stored>0){
				final int absNlz=(stored-1)+minZeros;
				// Approximate magnitude matching DLL2 restore() convention:
				// absNlz=0 → Long.MAX_VALUE (edge case), else 1L<<(wordlen-absNlz-1)
				final long dif;
				if(absNlz==0){dif=Long.MAX_VALUE;}
				else if(absNlz<wordlen){dif=1L<<(wordlen-absNlz-1);}
				else{dif=1L;}
				difSum+=dif;
				hllSumFilled+=Math.pow(2.0, -absNlz);
				gSum+=Math.log(Tools.max(1, dif));
				count++;
				sortBuf.add(dif);
			}
		}

		// All-buckets HLL: empty=1.0 (register=0 convention), filled=2^(-absNlz)
		double hllSum=0;
		for(int i=0; i<buckets; i++){
			final int stored=readBucket(i);
			if(stored==0){
				hllSum+=1.0;
			}else{
				hllSum+=Math.pow(2.0, -(stored-1)-minZeros);
			}
		}

		final double alpha_m=0.7213/(1.0+1.079/buckets);
		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		final double gmean=Math.exp(gSum/div);
		sortBuf.sort();
		final long median=Tools.max(1, sortBuf.median());
		final double mwa=Tools.max(1.0, sortBuf.medianWeightedAverage());
		final int V=buckets-count;

		// HLL-style all-buckets estimate with LC fallback at low occupancy
		final double hmeanRaw=2*alpha_m*(double)buckets*(double)buckets/hllSum;
		double hmeanEst=hmeanRaw;
		if(hmeanEst<2.5*buckets && V>0){hmeanEst=(double)buckets*Math.log((double)buckets/V);}

		final double correction=(count+buckets)/(float)(buckets+buckets);
		// Filled-bucket HMean (no fractional mantissa in DLL4, so HMeanM = HMean)
		final double hmeanPure=(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilled);
		final double hmeanPureM=hmeanPure; // no mantissa correction possible

		final double meanEst   =2*(Long.MAX_VALUE/Tools.max(1.0, mean))*div*correction;
		final double gmeanEst  =2*(Long.MAX_VALUE/gmean)               *div*correction;
		final double mwaEst    =2*(Long.MAX_VALUE/mwa)                 *div*correction;
		final double medianCorr=2*(Long.MAX_VALUE/(double)median)      *div*correction;
		final double lcPure    =buckets*Math.log((double)buckets/Math.max(V, 0.5));

		if(filledBuckets==0){return new double[9];}

		final double meanEstCF   =meanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEAN);
		final double hmeanPureMCF=hmeanPureM*CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEANM);

		final double hybridEst;
		final double hb0=0.5*buckets, hb1=1.5*buckets, hb2=3.0*buckets;
		if(lcPure<=hb0){
			hybridEst=lcPure;
		}else if(lcPure<=hb1){
			final double t=(lcPure-hb0)/(hb1-hb0);
			hybridEst=(1-t)*lcPure+t*meanEstCF;
		}else if(lcPure<=hb2){
			final double t=(lcPure-hb1)/(hb2-hb1);
			hybridEst=(1-t)*meanEstCF+t*hmeanPureMCF;
		}else{
			hybridEst=hmeanPureMCF;
		}

		return new double[]{
			meanEstCF,
			hmeanPure  *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEAN),
			hmeanPureM, // = hmeanPure; no extra CF since it equals HMean
			gmeanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.GMEAN),
			hmeanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HLL),
			lcPure,
			hybridEst,
			mwaEst     *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MWA),
			medianCorr *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEDCORR)
		};
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
	private final LongList sortBuf=new LongList(buckets);
	// lastCardinality inherited from CardinalityTracker

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, added);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;

	public static double LC_CROSSOVER=0.75;
	public static double LC_SHARPNESS=20.0;
	public static boolean USE_LC=true;
	
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

}
