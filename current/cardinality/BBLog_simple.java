package cardinality;

import shared.Parser;
import shared.Tools;

/**
 * Simplified LogLog-based cardinality tracker that estimates unique element counts
 * by storing the maximum hash observed in each bucket. Optimized for low memory
 * overhead while providing approximate distinct counts.
 *
 * @author Brian Bushnell
 * @date Feb 20, 2020
 */
public final class BBLog_simple extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	BBLog_simple(){
		this(2048, 31, -1, 0f);
	}
	
	BBLog_simple(Parser p){
		super(p);
		maxArray=new long[buckets];
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	BBLog_simple(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new long[buckets];
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	@Override
	public BBLog_simple copy() {return new BBLog_simple(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Estimates cardinality from bucket maxima using a LogLog-style mean difference
	 * calculation, caches the result in lastCardinality, and returns the estimate.
	 * @return Estimated number of unique elements
	 */
	@Override
	public final long cardinality(){
		double difSum=0;
		int count=0;
		
		for(int i=0; i<maxArray.length; i++){
			long val=maxArray[i];
			if(val>0){
				long dif=Long.MAX_VALUE-val;
				difSum+=dif;
				count++;
			}
		}
			
		final double mean=difSum/Tools.max(count, 1);
		final double estimatePerSet=2*(Long.MAX_VALUE/mean);
		final double total=estimatePerSet*count*((count+buckets)/(float)(buckets+buckets));
		
		long cardinality=(long)(total);
		lastCardinality=cardinality;
		return cardinality;
	}

	/** Returns the optional per-bucket count array if count tracking is enabled.
	 * @return Count array, or null when counts are not tracked */
	@Override
	public int[] getCounts(){
		return counts;
	}
	
	/**
	 * Merges another tracker of the same type into this one by delegating to the
	 * type-specific add implementation.
	 * @param log CardinalityTracker to merge (must be BBLog_simple)
	 */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((BBLog_simple)log);
	}
	
	public void add(BBLog_simple log){
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
			}
		}
	}
	
	/**
	 * Hashes the given number with Tools.hash64shift, selects a bucket, and stores the
	 * maximum hash observed for that bucket.
	 * @param number Value to hash and record
	 */
	@Override
	public void hashAndStore(final long number){
//		if(number%SKIPMOD!=0){return;}
//		final long key=hash(number, tables[((int)number)&numTablesMask]);
		final long key=Tools.hash64shift(number);
		
//		if(key<minKey){return;}
		final int bucket=(int)(key&bucketMask);
		
		{
			maxArray[bucket]=Tools.max(key, maxArray[bucket]);
		}
	}
	
	/**
	 * Returns compensation factors for log buckets; this implementation performs no
	 * compensation and returns null.
	 * @return null (no compensation factors)
	 */
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Buckets storing the maximum hash values used for cardinality estimation. */
	private final long[] maxArray;
	private final int[] counts;
	
//	private static long minKey=(long)(0.75f*Long.MAX_VALUE); //non-atomic 15% faster without this
	
}
