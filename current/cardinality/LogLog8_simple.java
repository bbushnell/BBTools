package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * Simple LogLog cardinality estimator using 8-bit buckets for memory efficiency.
 * Implements basic LogLog algorithm without mantissa tracking, storing only
 * the number of leading zeros (exponent) in each bucket. This reduces memory
 * usage compared to full LogLog implementations but may sacrifice some accuracy.
 *
 * @author Brian Bushnell
 * @date Mar 10, 2020
 */
public final class LogLog8_simple extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a LogLog with default parameters: 2048 buckets, k=31 */
	LogLog8_simple(){
		this(2048, 31, -1, 0);
	}
	
	/** Creates a LogLog with parameters parsed from command line arguments */
	LogLog8_simple(Parser p){
		super(p);
		maxArray=new byte[buckets];
	}
	
	/**
	 * Creates a LogLog with specified parameters.
	 *
	 * @param buckets_ Number of buckets (counters)
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	LogLog8_simple(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
	}
	
	@Override
	public LogLog8_simple copy(){return new LogLog8_simple(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Restores an approximate original value from a leading-zero count.
	 * Assumes mantissa=1 (no fractional bits stored).
	 * @param value Number of leading zeros stored in bucket
	 * @return Approximate original value
	 */
	private long restore(final int value){
		final int leading=value; //Number of leading zeros
		final long mantissa=1; //1.xxxx but in this case the X's are all zero
		final int shift=64-leading-1; //Amount to left shift the mantissa
		final long original=mantissa<<shift; //Restored original number
		return original;
	}
	
	/**
	 * Estimates cardinality from leading-zero statistics across all buckets.
	 * Uses arithmetic mean of non-zero buckets with empirically-derived
	 * mantissa correction and empty-bucket compensation.
	 */
	@Override
	public final long cardinality(){
		double sum=0;
		int count=0;

		for(int i=0; i<maxArray.length; i++){
			final int max=maxArray[i];
			final long val=restore(max);
			if(max>0 && val>0){
				sum+=val;
				count++;
			}
		}

		final int subsets=count;//Could be set to count or buckets
		final double mean=sum/Tools.max(subsets, 1);

		//What to use as the value from the counters
		final double proxy=mean;

		final double estimatePerSet=2*(Long.MAX_VALUE/proxy);
		final double mantissaFactor=0.7213428177;//Empirically derived
		final double emptyBucketModifier=((count+buckets)/(float)(buckets+buckets));//Approximate; overestimate
		final double total=estimatePerSet*subsets*mantissaFactor*emptyBucketModifier;

		final long cardinality=Math.min(added, (long)(total));
		lastCardinalityStatic=cardinality;
		return cardinality;
	}
	
	/** Merges another tracker into this one by taking maximum of each bucket.
	 * @param log Tracker to merge (must be LogLog8_simple) */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog8_simple)log);
	}
	
	/** Merges another LogLog8_simple by element-wise maximum of bucket arrays. */
	public void add(LogLog8_simple log){
		added+=log.added;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
			}
		}
	}
	
	/** {@inheritDoc} */
	@Override
	public void hashAndStore(final long number){
		final long key=Tools.hash64shift(number);
		final byte leading=(byte)Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);
		maxArray[bucket]=Tools.max(leading, maxArray[bucket]);
	}
	
	/** Returns null; this tracker does not use per-bucket compensation factors. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Bucket array storing maximum leading zeros count for each hash bucket */
	private final byte[] maxArray;

}
