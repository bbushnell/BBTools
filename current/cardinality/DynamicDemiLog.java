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
		double difSum=0;
		double hSum=0;
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
				hSum+=1.0/Tools.max(1, dif);
				gSum+=Math.log(Tools.max(1, dif));
				rSum+=Math.sqrt(dif);
				count++;
				double est=2*(Long.MAX_VALUE/(double)dif);
				estLogSum+=Math.log(est);
				list.add(dif);
			}
		}
		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		double hmean=hSum/div;
		double gmean=gSum/div;
		double rmean=rSum/div;
		hmean=1.0/hmean;
		gmean=Math.exp(gmean);
		rmean=rmean*rmean;
		list.sort();
		final long median=list.median();
		final double mwa=list.medianWeightedAverage();

		final double proxy=(USE_MEAN ? mean : USE_MEDIAN ? median : USE_MWA ? mwa :
			USE_HMEAN ? hmean : USE_GMEAN ? gmean : mean);

		final double estimatePerSet=2*(Long.MAX_VALUE/proxy);
		final double total=estimatePerSet*div*((count+buckets)/(float)(buckets+buckets));

		final double estSum=div*Math.exp(estLogSum/(Tools.max(div, 1)));
		final double medianEst=2*(Long.MAX_VALUE/(double)median)*div;

		// LinearCounting correction for sparse regimes: uses empty-bucket occupancy,
		// which carries more information than bucket values when cardinality << buckets.
		// Mathematically equivalent to Bloom filter occupancy estimation with hashes=1:
		//   n = -buckets * ln(1 - occupancy)
		final int V=buckets-count;  // empty buckets
		// Factor of 2: canonical k-mers use max(kmer, revcomp), halving the effective hash
		// space. Value-based estimate already corrects via 2*(MAX/proxy); LC needs the same.
		final double lcEstimate=(V>0 ? 2.0*buckets*Math.log((double)buckets/V) : total);

		// Sigmoid blend: smoothly transitions from LinearCounting (w=0) to value-based (w=1).
		// LC_CROSSOVER: occupancy fraction at which the two estimators contribute equally (w=0.5).
		//   Lower = trust LC longer; higher = switch to value-based sooner.
		// LC_SHARPNESS: steepness of the sigmoid. Higher = sharper transition (more like a hard
		//   switch); lower = more gradual blend. At sharpness 20, the transition from 10% to 90%
		//   value-based weight spans roughly ±11 occupancy percentage points around the crossover.
		// TODO: derive optimal values mathematically from relative estimator variances.
		double estimate=total;
		if(USE_LC) {
			final double occ=(double)count/buckets;
			final double w=1.0/(1.0+Math.exp(-LC_SHARPNESS*(occ-LC_CROSSOVER)));
			final double blended=(1-w)*lcEstimate+w*total;
			estimate=blended;
		}

		long cardinality=Math.min(added, (long)estimate);
		lastCardinality=cardinality;
		return cardinality;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicDemiLog)log);
	}

	public void add(DynamicDemiLog log){
		added+=log.added;
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
		//It's possible to exit early before calculating nlz but inefficient
		final int nlz=Long.numberOfLeadingZeros(key);

		// Global early exit
		if(nlz<minZeros){return;}

		//Precalculate everything necessary
		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		final int score=(nlz<<mantissabits)+(int)((~(key>>>shift))&mask); //FP16 representation
		final int oldValue=maxArray[bucket];//Required memory read

		//Optional early exit reduces writes, countArray access, and branches.
		//Expected to be usually taken, particularly when buckets is large.
		if(score<oldValue) {return;}

		final int newValue=Math.max(score, oldValue);
		final int nlzOld=(oldValue>>mantissabits);

		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
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
			/* 
			 * NOTE - Due to a major Eclipse JDK 24 -> Java 8 target bytecode generation bug,
			 * enabling both of these assertions causes a 40% slowdown even with -da.
			 * Do not enable them in production.
			 */
			//			assert(minZeroCount>=0 && minZeros<=64) : minZeroCount+", "+minZeros;
			while(minZeroCount==0 && minZeros<wordlen) {//Scan for new tier
				minZeros++;
				minZeroCount=countTermsInTier(minZeros, maxArray);
			}
			//			assert(minZeroCount>0 && minZeroCount<=buckets) : minZeroCount+", "+minZeros;
		}

	}

	private int scanFrom(int nlz) {
		minZeros=nlz-1;
		minZeroCount=0;
		//Technically this could be 2-pass: find the lowest, THEN count that tier
		//But in practice that would be slower
		while(minZeroCount==0 && minZeros<wordlen) {
			minZeros++;
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
	 * Returns all cardinality estimates for calibration, without the added-count cap.
	 * All value-based estimates use the same bucket-correction factor as cardinality().
	 * Blended uses mean as the value-based component (matching USE_MEAN default).
	 * <p>
	 * Array indices:
	 * 0=mean, 1=hmean, 2=gmean, 3=rmean, 4=mwa, 5=medianCorr,
	 * 6=medianLegacy (no correction), 7=estSum, 8=lc, 9=blended
	 */
	/**
	 * Array indices: 0=Mean, 1=HMean(HLL-style), 2=GMean, 3=RMean, 4=MWA,
	 * 5=MedianCorr, 6=MedianLeg, 7=EstSum, 8=LCHybrid, 9=LCTrue, 10=Blended
	 */
	public double[] rawEstimates(){
		double difSum=0;
		double gSum=0;
		double rSum=0;
		double estLogSum=0;
		int count=0;
		sortBuf.clear();
		final LongList list=sortBuf;

		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
			if(max>0 && val>0){
				long dif=val;
				difSum+=dif;
				gSum+=Math.log(Tools.max(1, dif));
				rSum+=Math.sqrt(dif);
				count++;
				double est=2*(Long.MAX_VALUE/(double)dif);
				estLogSum+=Math.log(est);
				list.add(dif);
			}
		}

		// HLL-style harmonic mean: sum 2^(-nlz) over ALL buckets (including empty).
		// Empty bucket score=0 -> nlz=0 -> 2^0=1, same as HyperLogLog register M_j=0.
		// alpha_m bias correction from Flajolet et al.
		double hllSum=0;
		for(int i=0; i<maxArray.length; i++){
			hllSum+=Math.pow(2.0, -(maxArray[i]>>mantissabits));
		}
		final double alpha_m=0.7213/(1.0+1.079/buckets);

		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		double gmean=gSum/div;
		double rmean=rSum/div;
		gmean=Math.exp(gmean);
		rmean=rmean*rmean;
		list.sort();
		final long median=Tools.max(1, list.median());
		final double mwa=Tools.max(1.0, list.medianWeightedAverage());

		final int V=buckets-count;

		// HLL estimate; fall back to LC for small estimates (standard HLL small-range correction).
		// Main formula: 2x retained — the HLL sum over all m buckets needs the factor to account
		// for the full signed 64-bit hash range including nlz=0 values.
		// LC fallback: 1x — restore() now returns Long.MAX_VALUE for nlz=0 (rather than overflowing
		// to negative and being excluded), so V correctly counts actual empty buckets; no 2x needed.
		double hmeanEst=2*alpha_m*(double)buckets*(double)buckets/hllSum;
		if(hmeanEst<2.5*buckets && V>0){
			hmeanEst=(double)buckets*Math.log((double)buckets/V);
		}

		final double correction=(count+buckets)/(float)(buckets+buckets);
		final double meanEst     =2*(Long.MAX_VALUE/Tools.max(1.0, mean)) *div*correction *MEAN_FACTOR;
		final double gmeanEst    =2*(Long.MAX_VALUE/gmean)                *div*correction *GMEAN_FACTOR;
		final double rmeanEst    =2*(Long.MAX_VALUE/Tools.max(1.0, rmean))*div*correction *RMEAN_FACTOR;
		final double mwaEst      =2*(Long.MAX_VALUE/mwa)                  *div*correction;  // MWA: skip factor; odd-list bug pending fix
		final double medianCorr  =2*(Long.MAX_VALUE/(double)median)       *div*correction *MEDIAN_FACTOR;
		final double medianLeg   =2*(Long.MAX_VALUE/(double)median)       *div            *MEDIAN_FACTOR;
		hmeanEst*=HMEAN_FACTOR;
		final double estSum      =div*Math.exp(estLogSum/Tools.max(div, 1))*ESTSUM_FACTOR;

		// LCHybrid: original behavior — falls back to meanEst when all buckets full.
		// Factor of 2 removed: restore() now returns Long.MAX_VALUE for nlz=0 (fixing the exclusion
		// bug), so V correctly reflects actual empty buckets. For production canonical k-mer use,
		// the factor-of-2 question is handled separately in cardinality().
		final double lcHybrid=(V>0 ? buckets*Math.log((double)buckets/V) : meanEst);

		// LCTrue: honest LC — uses max(1,V) so it doesn't blow up at full occupancy.
		// No factor of 2: restore() now returns Long.MAX_VALUE for nlz=0 scores (rather than
		// overflowing to negative and being excluded), so V correctly reflects actual empty buckets.
		final double lcTrue=buckets*Math.log((double)buckets/Math.max(1, V));

		double blended=meanEst;
		if(USE_LC){
			final double occ=(double)count/buckets;
			final double w=1.0/(1.0+Math.exp(-LC_SHARPNESS*(occ-LC_CROSSOVER)));
			blended=(1-w)*lcTrue+w*meanEst;
		}

		return new double[]{meanEst, hmeanEst, gmeanEst, rmeanEst, mwaEst,
			medianCorr, medianLeg, estSum, lcHybrid, lcTrue, blended};
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
	/** Reusable sort buffer for rawEstimates(); avoids per-call allocation. */
	private final LongList sortBuf=new LongList(buckets);

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Number of mantissa bits for floating-point compression; 10 is maximum. */
	private static final int mantissabits=10;
	private static final int mask=(1<<mantissabits)-1;
	private static final int offset=wordlen-mantissabits-1;

	/** Occupancy fraction (filled buckets / total buckets) at which LinearCounting and
	 * value-based estimates contribute equally in the sigmoid blend.  Below this point
	 * the blend leans toward LinearCounting; above it toward the value-based estimate. */
	public static double LC_CROSSOVER=0.75;
	/** Steepness of the sigmoid transition between LinearCounting and value-based estimation.
	 * Higher values create a sharper switch; lower values create a more gradual blend.
	 * At 20, the blend moves from 10% to 90% value-based weight over roughly ±11 occupancy
	 * percentage points around LC_CROSSOVER.  At 5 it spans roughly ±44 points. */
	public static double LC_SHARPNESS=20.0;
	public static boolean USE_LC=true;

	/**
	 * High-cardinality correction factors: flat multipliers applied to each estimator
	 * to remove systematic bias at maximum occupancy (load=100x, measured over 4000 DDLs,
	 * 2048 buckets, 100x cardinality via DDLCalibrationDriver, March 2026).
	 * <p>
	 * Each factor = 1 / (1 + mean_signed_error_at_max_cardinality).
	 * Applying these makes estimators unbiased at high load; some residual bias
	 * remains at intermediate load and will be addressed with occupancy-dependent curves.
	 * <p>
	 * MWA omitted: medianWeightedAverage() has an odd-list bug (count+=2*mult for
	 * center element; should be count+=mult), making it unreliable until fixed.
	 * EstSum omitted: fundamentally flawed (confirmed), pending removal.
	 */
	public static final double MEAN_FACTOR   = 0.99901637;
	public static final double HMEAN_FACTOR  = 0.99967171;
	public static final double GMEAN_FACTOR  = 0.56072095;
	public static final double RMEAN_FACTOR  = 0.78455962;
	public static final double MEDIAN_FACTOR = 1.67743074;
	/** EstSum: same high-load correction as GMean (both are geometric-mean-like). */
	public static final double ESTSUM_FACTOR = 0.56072095;

}
