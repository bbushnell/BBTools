package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * Simplified LogLog cardinality tracker using per-bucket max hash values.
 * Stripped-down variant of BBLog without geometric-mean or median estimates.
 *
 * @author Brian Bushnell
 * @date Feb 20, 2020
 */
public final class BBLog_simple extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Default constructor: 2048 buckets, k=31. */
	BBLog_simple(){this(2048, 31, -1, 0f);}

	/** Construct from parsed command-line arguments. */
	BBLog_simple(Parser p){
		super(p);
		maxArray=new long[buckets];
		counts=(trackCounts ? new int[buckets] : null);
	}

	/**
	 * Full constructor.
	 * @param buckets_ Number of buckets for hash partitioning
	 * @param k_ K-mer length
	 * @param seed Random seed (-1 for default)
	 * @param minProb_ Minimum probability threshold
	 */
	BBLog_simple(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new long[buckets];
		counts=(trackCounts ? new int[buckets] : null);
	}

	/** Create an independent copy with a fresh seed. */
	@Override
	public BBLog_simple copy(){return new BBLog_simple(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Estimates cardinality from bucket maxima using arithmetic-mean LogLog.
	 * @return Estimated number of unique elements
	 */
	@Override
	public final long cardinality(){
		double difSum=0;
		int count=0;
		for(int i=0; i<maxArray.length; i++){
			final long val=maxArray[i];
			if(val>0){
				final long dif=Long.MAX_VALUE-val;
				difSum+=dif;
				count++;
			}
		}
		final double mean=difSum/Tools.max(count, 1);
		final double estimatePerSet=2*(Long.MAX_VALUE/mean);
		final double total=estimatePerSet*count*((count+buckets)/(float)(buckets+buckets));

		final long cardinality=Math.min(added, (long)(total));
		lastCardinalityStatic=cardinality;
		return cardinality;
	}

	/** Returns per-bucket count array, or null if count tracking is disabled. */
	@Override
	public int[] getCounts(){return counts;}
	
	/** Merge another tracker into this one (must be BBLog_simple). */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((BBLog_simple)log);
	}

	/** Merge another BBLog_simple by taking per-bucket max hash values. */
	public void add(BBLog_simple log){
		added+=log.added;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
			}
		}
	}
	
	/** Hash a value and store the maximum hash per bucket. */
	@Override
	public void hashAndStore(final long number){
		final long key=Tools.hash64shift(number);
		final int bucket=(int)(key&bucketMask);
		maxArray[bucket]=Tools.max(key, maxArray[bucket]);
	}
	
	/** No compensation factors used by this estimator. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum hash value observed per bucket. */
	private final long[] maxArray;
	/** Per-bucket occurrence counts; null if count tracking is disabled. */
	private final int[] counts;


}
