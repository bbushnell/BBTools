package cardinality;

import shared.Parser;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Mar 10, 2020
 *
 */
public final class LogLog8_simple extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create a LogLog with default parameters */
	LogLog8_simple(){
		this(2048, 31, -1, 0);
	}
	
	/** Create a LogLog with parsed parameters */
	LogLog8_simple(Parser p){
		super(p);
		maxArray=new byte[buckets];
	}
	
	/**
	 * Create a LogLog with specified parameters
	 * @param buckets_ Number of buckets (counters)
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	LogLog8_simple(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Restores floating point to integer.
	 * This subclass has no mantissa so only the exponent is restored. */
	private long restore(int value){
		final int leading=value; //Number of leading zeros
		long mantissa=1; //1.xxxx but in this case the X's are all zero
		int shift=64-leading-1; //Amount to left shift the mantissa
		long original=mantissa<<shift; //Restored original number
		return original;
	}
	
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
	
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog8_simple)log);
	}
	
	/** @See add(CardinalityTracker) */
	public void add(LogLog8_simple log){
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
			}
		}
	}
	
	@Override
	public void hashAndStore(final long number){
		final long key=Tools.hash64shift(number);
		final byte leading=(byte)Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);
		maxArray[bucket]=Tools.max(leading, maxArray[bucket]);
	}
	
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Maintains state.  These are the actual buckets. */
	private final byte[] maxArray;

}
