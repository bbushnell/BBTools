package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * Inverted LogLog cardinality estimator.
 * <p>
 * Where standard LogLog uses leading zeros (NLZ) as the bucket selector and stores
 * the maximum NLZ score per bucket, this estimator uses NLZ as the bucket selector
 * but stores the maximum lower-32-bits value seen in each bucket.  The two fields
 * (high bits for bucket selection, low bits for the stored value) are orthogonal and
 * carry independent information.
 * <p>
 * At high cardinality, low-NLZ buckets saturate (max low-bits ≈ MAX_VALUE) while
 * high-NLZ buckets remain empty.  In between, a band of informative buckets exists
 * where diff[i] = MAX_VALUE - maxArray[i] follows the geometric progression
 * E[diff[i]] ≈ MAX_VALUE * 2^(i+1) / cardinality.  On a log2 scale this is a
 * straight line with slope exactly 1 (known from theory) — only the intercept
 * changes with cardinality.  Cardinality is estimated from the geometric mean of
 * per-bucket estimates across all informative buckets.
 * <p>
 * Default size: 32 buckets × 32 bits = 128 bytes (same footprint as 128-bucket 8-bit LogLog).
 *
 * @author Brian Bushnell
 * @contributor Chloe
 * @date March 2026
 */
public final class InvertedLogLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an InvertedLogLog with default parameters. */
	InvertedLogLog(){
		this(2048, 31, -1, 0);
	}

	/** Creates an InvertedLogLog with parameters parsed from command line arguments. */
	InvertedLogLog(Parser p){
		super(p);
		maxArray=new int[BUCKETS];
	}

	/**
	 * Creates an InvertedLogLog with specified parameters.
	 * Bucket count and k-mer length are accepted for API compatibility but the
	 * number of buckets is fixed at BUCKETS=32.
	 */
	InvertedLogLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[BUCKETS];
	}

	@Override
	public InvertedLogLog copy(){return new InvertedLogLog(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Hashes the input and updates the bucket for the NLZ of the hashed value.
	 * Stores the maximum lower-32-bits value seen in that bucket.
	 */
	@Override
	public void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);
		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(nlz<BUCKETS ? nlz : BUCKETS-1);
		final int lowBits=(int)(key);  // lower 32 bits
		if(Integer.compareUnsigned(lowBits, maxArray[bucket])>0){
			maxArray[bucket]=lowBits;
		}
		added++;
	}

	/**
	 * Estimates cardinality from the median of per-bucket log estimates.
	 * <p>
	 * For each informative bucket i (neither saturated nor empty):
	 *   diff[i] = UINT_MAX - maxArray[i]  (unsigned)
	 *   log estimate: logEst_i = log(UINT_MAX+1) + (i+1)*log(2) - log(diff[i])
	 * <p>
	 * Median is used instead of mean to avoid the upward bias from Jensen's inequality
	 * (E[log(diff)] < log(E[diff])) and to suppress outliers from near-saturated or
	 * near-empty buckets.  Also computes MWA (triangular-weighted average around median)
	 * for comparison; currently returns median.
	 */
	@Override
	public long cardinality(){
		final double[] logEsts=new double[BUCKETS];
		int count=0;
		for(int i=0; i<BUCKETS; i++){
			if(maxArray[i]==0){continue;}  // empty: never updated
			final long diff=UINT_MAX-(maxArray[i]&UINT_MASK);
			if(diff==0){continue;}  // saturated
			logEsts[count++]=LOG_UINT_MAX_PLUS1+(i+1)*LOG2-Math.log(diff);
		}
		if(count==0){return 0;}

		// Sort informative estimates; median is the L1-optimal robust center
		java.util.Arrays.sort(logEsts, 0, count);
//		final double medianLogEst=(count%2==1)
//			? logEsts[count/2]
//			: (logEsts[count/2-1]+logEsts[count/2])*0.5;

//		 MWA: triangular weights peaking at median, for future comparison
		 double mwaLogEst=medianWeightedAverage(logEsts, count);

		final long cardinality=(long)Math.exp(mwaLogEst);
		lastCardinality=cardinality;
		return cardinality;
	}

	/** Triangular-weighted average of sorted values (weight peaks at center). */
	private static double medianWeightedAverage(final double[] sorted, final int n){
		if(n==1){return sorted[0];}
		final int half=n/2;
		double sum=0;
		long weight=0;
		for(int i=0, j=n-1; i<half; i++, j--){
			final int mult=i+1;
			sum+=(sorted[i]+sorted[j])*mult;
			weight+=2*mult;
		}
		if((n&1)==1){
			final int mult=half+1;
			sum+=sorted[half]*mult;
			weight+=mult;
		}
		return sum/weight;
	}

	/** Merges another InvertedLogLog into this one by taking element-wise unsigned maxima. */
	@Override
	public void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((InvertedLogLog)log);
	}

	public void add(InvertedLogLog log){
		added+=log.added;
		if(maxArray!=log.maxArray){
			for(int i=0; i<BUCKETS; i++){
				if(Integer.compareUnsigned(log.maxArray[i], maxArray[i])>0){
					maxArray[i]=log.maxArray[i];
				}
			}
		}
	}

	@Override
	public float[] compensationFactorLogBucketsArray(){return null;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum lower-32-bits value seen per NLZ bucket.  Stored as signed int, compared unsigned. */
	private final int[] maxArray;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Fixed number of buckets.  Covers NLZ values 0..31; NLZ>=32 is clamped to bucket 31. */
	private static final int BUCKETS=32;
	/** Unsigned max of a 32-bit value. */
	private static final long UINT_MAX=0xFFFFFFFFL;
	/** Mask to convert signed int to unsigned long. */
	private static final long UINT_MASK=0xFFFFFFFFL;
	/** log((UINT_MAX+1)) = log(2^32) = 32*log(2), precomputed. */
	private static final double LOG_UINT_MAX_PLUS1=32*Math.log(2);
	/** log(2), precomputed. */
	private static final double LOG2=Math.log(2);

}
