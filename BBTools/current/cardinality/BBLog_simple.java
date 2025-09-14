package cardinality;

import shared.Parser;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Feb 20, 2020
 *
 */
public final class BBLog_simple extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create a LogLog with default parameters */
	BBLog_simple(){
		this(2048, 31, -1, 0f);
	}
	
	/** Create a LogLog with parsed parameters */
	BBLog_simple(Parser p){
		super(p);
		maxArray=new long[buckets];
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	/**
	 * Create a LogLog with specified parameters
	 * @param buckets_ Number of buckets (counters)
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	BBLog_simple(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new long[buckets];
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
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

	@Override
	public int[] getCounts(){
		return counts;
	}
	
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
	
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maintains state.  These are the actual buckets. */
	private final long[] maxArray;
	/** Counts associated with values in the buckets. */
	private final int[] counts;
	
//	private static long minKey=(long)(0.75f*Long.MAX_VALUE); //non-atomic 15% faster without this
	
}
