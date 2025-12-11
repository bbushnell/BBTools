package cardinality;

import shared.Parser;
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
	public LogLog8_simple copy() {return new LogLog8_simple(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Restores the approximate original value from the stored leading zeros count.
	 * Since this implementation has no mantissa, assumes mantissa = 1.0000...
	 * and reconstructs the original number by left-shifting 1 by (64-leadingZeros-1).
	 *
	 * @param value Number of leading zeros stored in bucket
	 * @return Approximate original value that produced this leading zeros count
	 */
	private long restore(int value){
		final int leading=value; //Number of leading zeros
		long mantissa=1; //1.xxxx but in this case the X's are all zero
		int shift=64-leading-1; //Amount to left shift the mantissa
		long original=mantissa<<shift; //Restored original number
		return original;
	}
	
	/**
	 * Calculates cardinality estimate using harmonic mean of restored bucket values.
	 * Applies empirically-derived correction factors for mantissa approximation
	 * and empty bucket compensation. Uses mean of non-zero buckets rather than
	 * all buckets to reduce bias from uninitialized counters.
	 *
	 * @return Estimated number of unique elements added to this tracker
	 */
	@Override
	public final long cardinality(){
		double sum=0;
		int count=0;
		
		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
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
		
		long cardinality=(long)(total);
		lastCardinality=cardinality;
		return cardinality;
	}
	
	/** Merges another tracker into this one by taking maximum of each bucket.
	 * @param log Tracker to merge (must be LogLog8_simple) */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog8_simple)log);
	}
	
	/**
	 * Merges another LogLog8_simple tracker by taking element-wise maximum
	 * of bucket arrays. This preserves the maximum leading zeros count
	 * seen for each bucket across both trackers.
	 * @param log LogLog8_simple tracker to merge
	 */
	public void add(LogLog8_simple log){
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
			}
		}
	}
	
	/**
	 * Hashes a number and updates the corresponding bucket with leading zeros count.
	 * Uses Tools.hash64shift for hashing, extracts bucket from low bits, and
	 * stores the maximum leading zeros count seen for this bucket.
	 * @param number Value to hash and track
	 */
	@Override
	public void hashAndStore(final long number){
		final long key=Tools.hash64shift(number);
		final byte leading=(byte)Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);
		maxArray[bucket]=Tools.max(leading, maxArray[bucket]);
	}
	
	/**
	 * Returns compensation factors for bucket count bias correction.
	 * This implementation returns null, indicating no special compensation.
	 * @return null (no compensation array)
	 */
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Bucket array storing maximum leading zeros count for each hash bucket */
	private final byte[] maxArray;

}
