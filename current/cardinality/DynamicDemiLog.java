package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicDemiLog cardinality estimator extending DemiLogLog (LogLog16) with an adaptive
 * early-exit mechanism. Maintains two extra integers for metadata:
 * minZeros: The minimum number of leading zeros seen across all buckets
 * minZeroCount: Number of buckets with minZeros
 * This drives a dynamic threshold that starts at zero and rises with cardinality, achieving
 * 99.9%+ early exit for large datasets without any loss of accuracy—skipped elements can
 * never improve a bucket that already holds a higher score.
 *
 * @author Brian Bushnell
 * @contributor Chloe
 * @date February 27, 2026
 */
public final class DynamicDemiLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a DynamicLogLog with default parameters.
	 * Uses 2048 buckets, k=31, random seed, and no minimum probability filtering. */
	DynamicDemiLog(){
		this(2048, 31, -1, 0);
	}

	/** Creates a DynamicLogLog with parameters parsed from command-line arguments.
	 * @param p Parser containing configuration from command-line flags */
	DynamicDemiLog(Parser p){
		super(p);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	/**
	 * Creates a DynamicLogLog with specified parameters.
	 * @param buckets_ Number of buckets (counters) for the hash table
	 * @param k_ K-mer length for sequence hashing
	 * @param seed Random number generator seed; -1 for random seed
	 * @param minProb_ Ignore k-mers with under this probability of being correct
	 */
	DynamicDemiLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	@Override
	public DynamicDemiLog copy(){return new DynamicDemiLog(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0) {return lastCardinality;}
		double difSum=0;
		double hllSumFilled=0;
		double gSum=0;
		double rSum=0;
		double estLogSum=0;
		int count=0;
		LongList list=new LongList(buckets);

		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
			if(max>0 && val>0){
				long dif=val;
				difSum+=dif;
				hllSumFilled+=Math.pow(2.0, -(max>>mantissabits));
				gSum+=Math.log(Tools.max(1, dif));
				rSum+=Math.sqrt(dif);
				count++;
				double est=2*(Long.MAX_VALUE/(double)dif);
				estLogSum+=Math.log(est);
				list.add(dif);
			}
		}
		final double alpha_m=0.7213/(1.0+1.079/buckets);
		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		double gmean=Math.exp(gSum/div);
		list.sort();
		final long median=list.median();
		final double mwa=list.medianWeightedAverage();
		// HLL-style filled-bucket HMean — already a full cardinality estimate.
		final double hmean=(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilled);

		final double correction=(count+buckets)/(float)(buckets+buckets);
		final int V=buckets-count;

		// HMean is already a full estimate; all other estimators work via a dif proxy.
		double total;
		int cfType;
		if(USE_HMEAN && count>0){
			total=hmean;
			cfType=CorrectionFactor.HMEAN;
		}else{
			final double proxy=(USE_MEDIAN ? median : USE_MWA ? mwa : USE_GMEAN ? gmean : mean);
			cfType=(USE_MEDIAN ? CorrectionFactor.MEDCORR :
				USE_MWA ? CorrectionFactor.MWA :
					USE_GMEAN ? CorrectionFactor.GMEAN : CorrectionFactor.MEAN);
			total=2*(Long.MAX_VALUE/proxy)*div*correction;
		}
		total*=CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets, cfType);

		// LinearCounting correction for sparse regimes.
		// Factor of 2: canonical k-mers use max(kmer, revcomp), halving the effective hash space.
		final double lcEstimate=(V>0 ? 2.0*buckets*Math.log((double)buckets/V) : total);

		// Sigmoid blend: smoothly transitions from LinearCounting (w=0) to value-based (w=1).
		double estimate=total;
		if(USE_LC){
			final double occ=(double)count/buckets;
			final double w=1.0/(1.0+Math.exp(-LC_SHARPNESS*(occ-LC_CROSSOVER)));
			estimate=(1-w)*lcEstimate+w*total;
		}

		long cardinality=Math.min(added, (long)estimate);
		lastCardinalityStatic=lastCardinality=cardinality;
		return cardinality;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicDemiLog)log);
	}

	public void add(DynamicDemiLog log){
		added+=log.added;
		branch1+=log.branch1;
		branch2+=log.branch2;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				final char maxA=maxArray[i], maxB=log.maxArray[i];
				final char countA=countArray[i], countB=log.countArray[i];
				maxArray[i]=Tools.max(maxA, maxB);
				if(maxA==maxB){countArray[i]=(char)Tools.min(countA+(int)countB, Character.MAX_VALUE);}
				else{countArray[i]=(maxA>maxB ? countA : countB);}
			}
			scanFrom(Math.max(minZeros, log.minZeros));
			//Recompute filledBuckets after merge; incremental tracking is not possible here
			filledBuckets=0;
			for(char c : maxArray){if(c>0){filledBuckets++;}}
		}
	}

	/** Compares two DynamicLogLog instances bucket by bucket.
	 * @return int[] {lower, equal, higher} bucket counts */
	public int[] compareTo(DynamicDemiLog blog){return compare(maxArray, blog.maxArray);}

	@Override
	public void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		//Earliest possible exit, fastest
		if(Long.compareUnsigned(key, eeMask)>0) {return;}
//		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);

		// Global early exit if not using eeMask
//		if(nlz<minZeros){return;}

		//Precalculate everything necessary
		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		final int score=(nlz<<mantissabits)+(int)((~(key>>>shift))&mask);//FP16 representation
		final int oldValue=maxArray[bucket];//Required memory read

		//Optional early exit reduces writes, countArray access, and branches.
		//Expected to be usually taken, particularly when buckets is large.
		if(score<oldValue) {return;}
//		branch2++;
		lastCardinality=-1;
		final int newValue=Math.max(score, oldValue);
		final int nlzOld=(oldValue>>mantissabits);

//		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
		//Update bucket; required write to cached line only if score>oldValue
		maxArray[bucket]=(char)newValue;
		//Track filled bucket count: increment when bucket transitions from empty to non-empty
		if(oldValue==0 && newValue>0){filledBuckets++;}

		//Update count - totally optional, only for histograms
		//Can be debranched with clever math
		final char count=countArray[bucket];
		countArray[bucket]=(char)(oldValue>score ? count :
			oldValue==score ? Math.max(count, (char)(count+1)) : 1);

		//Update the dynamic early exit threshold
		if(nlz>nlzOld && nlzOld==minZeros && --minZeroCount<1) {//Promotion from bottom tier 
			/* NOTE - Due to a major Eclipse JDK 24 -> Java 8 target bytecode generation bug,
			 * enabling both of these assertions causes a 40% slowdown even with -da.
			 * Do not enable them in production.
			 */
			//			assert(minZeroCount>=0 && minZeros<=64) : minZeroCount+", "+minZeros;
			while(minZeroCount==0 && minZeros<wordlen) {//Scan for new tier
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countTermsInTier(minZeros, maxArray);
			}
			//			assert(minZeroCount>0 && minZeroCount<=buckets) : minZeroCount+", "+minZeros;
		}

	}

	private int scanFrom(int nlz) {
		minZeros=nlz-1;
		minZeroCount=0;
		eeMask=-1L;
		//Technically this could be 2-pass: find the lowest, THEN count that tier
		//But in practice that would be slower
		while(minZeroCount==0 && minZeros<wordlen) {
			minZeros++;
			eeMask>>>=1;
			minZeroCount=countTermsInTier(minZeros, maxArray);
		}
		return minZeros;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public char[] counts16(){return countArray;}

	/** Returns the number of buckets with any data (maxArray[i] > 0). O(1). */
	public int filledBuckets(){return filledBuckets;}

	/** Returns the fraction of buckets with any data (filled / total). O(1). */
	public double occupancy(){return (double)filledBuckets/buckets;}

	/**
	 * Returns cardinality estimates for calibration, without the added-count cap.
	 * Output order: 0=Mean, 1=HMean, 2=HMeanM, 3=GMean, 4=HLL, 5=LC, 6=Hybrid, 7=MWA, 8=MedianCorr, 9=Mean99
	 * Hybrid (index 6) is pre-corrected: CF already applied to its Mean and HMeanM components.
	 * LC (index 5) and Hybrid (index 6) receive no additional CF.
	 */
	public double[] rawEstimates(){
		double difSum=0;
		double hllSumFilled=0;
		double hllSumFilledM=0;
		double gSum=0;
		int count=0;
		sortBuf.clear();
		final LongList list=sortBuf;

		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
			if(max>0 && val>0){
				long dif=val;
				difSum+=dif;
				hllSumFilled+=Math.pow(2.0, -(max>>mantissabits));
				hllSumFilledM+=Math.pow(2.0, -(max>>mantissabits)+0.5-(max&mask)/1024.0);
				gSum+=Math.log(Tools.max(1, dif));
				count++;
				list.add(dif);
			}
		}

		// HLL-style sum over ALL buckets (including empty); alpha_m bias correction from Flajolet et al.
		double hllSum=0;
		for(int i=0; i<maxArray.length; i++){
			hllSum+=Math.pow(2.0, -(maxArray[i]>>mantissabits));
		}
		final double alpha_m=0.7213/(1.0+1.079/buckets);

		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		double gmean=Math.exp(gSum/div);
		list.sort();
		final long median=Tools.max(1, list.median());
		final double mwa=Tools.max(1.0, list.medianWeightedAverage());

		final int V=buckets-count;

		// HLL: all-buckets formula with LC small-range fallback.
		final double hmeanRaw=2*alpha_m*(double)buckets*(double)buckets/hllSum;
		double hmeanEst=hmeanRaw;
		if(hmeanEst<2.5*buckets && V>0){hmeanEst=(double)buckets*Math.log((double)buckets/V);}

		final double correction=(count+buckets)/(float)(buckets+buckets);
		final double hmeanPure =(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilled);
		final double hmeanPureM=(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilledM);
		final double meanEst   =2*(Long.MAX_VALUE/Tools.max(1.0, mean))*div*correction;
		final double gmeanEst  =2*(Long.MAX_VALUE/gmean)               *div*correction;
		final double mwaEst    =2*(Long.MAX_VALUE/mwa)                 *div*correction;
		final double medianCorr=2*(Long.MAX_VALUE/(double)median)      *div*correction;
		// Pure LC: clamped denominator prevents division by zero; no CF ever applied.
		final double lcPure    =buckets*Math.log((double)buckets/Math.max(V, 0.5));

		// Mean99: trimmed mean ignoring upper and lower 1/256 of values.
		// "Upper" = highest estimates = smallest dif (left of sorted list): trim count/256.
		// "Lower" = lowest estimates = largest dif (right of sorted list): empty buckets (V)
		//   already represent this tail, so only trim max(0, count/256 - V) from the right.
		final int trim=count/256;
		final int trimLow=Math.max(0, trim-V);
		final int mean99Start=trim;
		final int mean99End=count-trimLow; // exclusive
		double mean99Sum=0;
		final int mean99N=mean99End-mean99Start;
		if(mean99N>0){
			for(int i=mean99Start; i<mean99End; i++){mean99Sum+=list.get(i);}
		}
		final double mean99=(mean99N>0 ? mean99Sum/mean99N : mean);

		if(filledBuckets==0){return new double[10];}

		// Apply per-occupancy correction factors; getCF returns 1 when USE_CORRECTION=false.
		final double meanEstCF    =meanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEAN);
		final double hmeanPureMCF =hmeanPureM*CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEANM);

		// Hybrid: LC → Mean → HMeanM blend, driven by lcPure as a bucket-count-independent cardinality proxy.
		// [0, 0.5b]=LC, [0.5b, 1.5b]=LC→Mean, [1.5b, 3b]=Mean→HMeanM, [3b+]=HMeanM.
		// CF already applied to meanEstCF and hmeanPureMCF; hybridEst needs no additional CF.
		final double hybridEst;
		final double hb0=0.25*buckets, hb1=1.7*buckets, hb2=3.2*buckets;
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

		final double mean99Est=2*(Long.MAX_VALUE/Tools.max(1.0, mean99))*div*correction;
		return new double[]{
			meanEstCF,
			hmeanPure *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEAN),
			hmeanPureMCF,
			gmeanEst  *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.GMEAN),
			hmeanEst  *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HLL),
			lcPure,    // LC: no CF applied
			hybridEst, // Hybrid: CF already applied to components
			mwaEst    *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MWA),
			medianCorr*CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEDCORR),
			mean99Est *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEAN99)
		};
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Restores floating-point compressed score to approximate original hash magnitude.
	 * @param score Compressed 16-bit value from maxArray
	 * @return Approximate original hash value before compression */
	private static long restore(int score){
		int leading=(int)(score>>>mantissabits);
		// nlz=0 case: true hash magnitude is in [2^63, 2^64-1] (unsigned), which
		// overflows signed long (shift of 53 bits causes Long.MIN_VALUE). Approximate
		// as Long.MAX_VALUE — this keeps the bucket contributing to value-based estimates
		// (giving meanEst≈1 for a single nlz=0 bucket) rather than being silently excluded.
		if(leading==0){return Long.MAX_VALUE;}
		long lowbits=(~score)&mask;
		long mantissa=(1L<<mantissabits)|lowbits;
		int shift=wordlen-leading-mantissabits-1;
		long original=mantissa<<shift;
		return original;
	}

	private static int countTermsInTier(final int nlz, final char[] array) {
		int sum=0;
		for(final char c : array){
			final int diff=(c>>mantissabits)^nlz;  //0 iff tier matches
			final int equalsBit=((~(diff|-diff))>>>31);  //1 if diff==0, else 0
			sum+=equalsBit;
		}
		return sum;
	}

	/** Compares two maxArrays bucket by bucket using branchless bit arithmetic.
	 * @return int[]{lower, equal, higher} bucket counts */
	public static int[] compare(char[] arrayA, char[] arrayB){
		int lower=0, equal=0, higher=0;
		for(int i=0; i<arrayA.length; i++){
			final int a=arrayA[i], b=arrayB[i];
			int dif=a-b;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher};
	}

	public static float ani(int lower, int equal, int higher) {
		//TODO
		return -1;
	}

	public static float completeness(int lower, int equal, int higher) {
		//TODO
		return -1;
	}

	public static float contam(int lower, int equal, int higher) {
		//TODO
		return -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compressed 16-bit bucket maxima; 6-bit leading-zero count + 10-bit mantissa. */
	private final char[] maxArray;
	/** Count of observations at the current maximum score for each bucket. */
	private final char[] countArray;
	/** Minimum leading-zero tier held by any bucket; dynamic early-exit threshold.
	 * Rises monotonically from 0 toward 63 as cardinality grows. */
	private int minZeros=0;
	/** Number of buckets with minZeros */
	private int minZeroCount;
	/** Number of buckets with any data (maxArray[i] > 0). Maintained incrementally
	 * in hashAndStore() for O(1) filledBuckets() and occupancy() calls. */
	private int filledBuckets=0;
	private long eeMask=-1L;
	/** Reusable sort buffer for rawEstimates(); avoids per-call allocation. */
	private final LongList sortBuf=new LongList(buckets);
	// lastCardinality inherited from CardinalityTracker

	/*--------------------------------------------------------------*/

	/** For tracking branch prediction rates; disable in production */
	public long branch1=0, branch2=0;
	public double branch1Rate() {return branch1/(double)Math.max(1, added);}
	public double branch2Rate() {return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Number of mantissa bits for floating-point compression; 10 is maximum. */
	private static final int mantissabits=10;
	private static final int mask=(1<<mantissabits)-1;
	private static final int offset=wordlen-mantissabits-1;

	//	/** Occupancy fraction (filled buckets / total buckets) at which LinearCounting and
	//	 * value-based estimates contribute equally in the sigmoid blend.  Below this point
	//	 * the blend leans toward LinearCounting; above it toward the value-based estimate. */

	/** Default resource file for DDL correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDDL.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the DDL correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	public static double LC_CROSSOVER=0.75;
	/** Steepness of the sigmoid transition between LinearCounting and value-based estimation.
	 * Higher values create a sharper switch; lower values create a more gradual blend.
	 * At 20, the blend moves from 10% to 90% value-based weight over roughly ±11 occupancy
	 * percentage points around LC_CROSSOVER.  At 5 it spans roughly ±44 points. */
	public static double LC_SHARPNESS=20.0;
	public static boolean USE_LC=true;

}
