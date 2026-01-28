package cardinality;

import shared.Parser;
import shared.Tools;
import structures.LongList;

/**
 * An 8-bit implementation of the LogLog cardinality estimation algorithm for
 * memory-efficient probabilistic counting of unique elements.
 * Uses byte arrays to track maximum leading zero counts across buckets,
 * providing space-efficient cardinality estimation with O(log log n) complexity.
 *
 * @author Brian Bushnell
 * @date Mar 6, 2020
 */
public final class LogLog8 extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a LogLog8 with default parameters (2048 buckets, k=31, random seed)
	 */
	LogLog8(){
		this(2048, 31, -1, 0);
	}
	
	/** Creates a LogLog8 with parameters parsed from command line arguments.
	 * @param p Parser containing configuration parameters */
	LogLog8(Parser p){
		super(p);
		maxArray=new byte[buckets];
	}
	
	/**
	 * Creates a LogLog8 with specified cardinality tracking parameters.
	 *
	 * @param buckets_ Number of buckets for cardinality estimation
	 * @param k_ K-mer length for hashing
	 * @param seed Random number generator seed; -1 for random seed
	 * @param minProb_ Minimum probability threshold for element inclusion
	 */
	LogLog8(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
	}
	
	@Override
	public LogLog8 copy() {return new LogLog8(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
//	/** Restores floating point to integer */
//	private long restore2(int score, int bucket){
//		final int bucketBits=Integer.numberOfTrailingZeros(buckets);
//		int leading=score;
//		long mantissa=(1L<<bucketBits)|bucket;
//		int shift=wordlen-leading-bucketBits-1;
//		long original=mantissa<<shift;
//		return original;
//	}
	
	/**
	 * Restores the original hash value from a leading zero count.
	 * Converts the stored leading zero score back to an approximate original value
	 * using bit manipulation to reconstruct the mantissa and shift operations.
	 *
	 * @param score Number of leading zeros in the original hash
	 * @return Reconstructed approximate original hash value
	 */
	private long restore(int score){
		int leading=score;
		long mantissa=1;
		int shift=64-leading-1;
		long original=mantissa<<shift;
		return original;
	}
	
	/**
	 * Estimates the total number of unique elements seen by this tracker.
	 * Uses statistical analysis of leading zero counts across all buckets
	 * to calculate cardinality. Applies conversion factor and bucket occupancy
	 * adjustments for improved accuracy.
	 *
	 * @return Estimated cardinality of unique elements
	 */
	@Override
	public final long cardinality(){
		double difSum=0;
		double estLogSum=0;
		int count=0;
		LongList list=new LongList(buckets);
		
		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
			if(max>0 && val>0){
				final long dif=val;
				difSum+=dif;
				count++;
				double est=2*(Long.MAX_VALUE/(double)dif)*SKIPMOD;
				estLogSum+=Math.log(est);
				list.add(dif);
			}
		}
			
		final int div=count;//Could also be count be that causes problems
		final double mean=difSum/Tools.max(div, 1);
		list.sort();
		final long median=list.median();
		final double mwa=list.medianWeightedAverage();
		
		//What to use as the value from the counters 
		final double proxy=mean;
		
//		assert(false) : mean+", "+median+", "+difSum+", "+list;
		
		final double estimatePerSet=2*(Long.MAX_VALUE/proxy)*SKIPMOD;
		
		//12000000        16635460.58 //8k sims, 100k reads, 128k buckets
		//16635789.16  //16k sims, 100k reads, 128k buckets
		//16635901.26  //64k sims
		//16635476.90  //128k
		//16635631.18  //256k
		//16635645.59  //512k
		
		//0.72134379048167576520498945661274
		//0.72134281774212774758647730006006
		final double conversionFactor=0.7213428177;
		final double total=conversionFactor*estimatePerSet*div*((count+buckets)/(float)(buckets+buckets));

//		final double estSum=div*Math.exp(estLogSum/(Tools.max(div, 1)));
//		double medianEst=2*(Long.MAX_VALUE/(double)median)*SKIPMOD*div;
		
//		new Exception().printStackTrace();
		
//		System.err.println(maxArray);
////		Overall, it looks like "total" is the best, then "estSum", then "medianEst" is the worst, in terms of variance.
//		System.err.println("difSum="+difSum+", count="+count+", mean="+mean+", est="+estimatePerSet+", total="+(long)total);
//		System.err.println("estSum="+(long)estSum+", median="+median+", medianEst="+(long)medianEst);
		
		long cardinality=Math.min(added, (long)(total));
		lastCardinality=cardinality;
		return cardinality;
	}
	
	/** Merges another cardinality tracker into this one by taking maximum values.
	 * @param log The cardinality tracker to merge (must be LogLog8) */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog8)log);
	}
	
	public void add(LogLog8 log){
		added+=log.added;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
			}
		}
	}
	
	@Override
	public void hashAndStore(final long number){
//		if(number%SKIPMOD!=0){return;} //Slows down moderately
		long key=number;
		
//		key=hash(key, tables[((int)number)&numTablesMask]);
		
		key=Tools.hash64shift(key);
//		if(key<0 || key>maxHashedValue){return;}//Slows things down by 50% lot, mysteriously
		byte leading=(byte)(Long.numberOfLeadingZeros(key)&63);//mask is used to keep number in 6 bits 
		
//		counts[leading]++;
		
//		if(leading<3){return;}//Slows things down slightly
//		final int bucket=(int)((number&Integer.MAX_VALUE)%buckets);
		final int bucket=(int)(key&bucketMask);
		
		maxArray[bucket]=Tools.max(leading, maxArray[bucket]);
	}
	
	/**
	 * Returns null as LogLog8 does not use compensation factors.
	 * This implementation relies on the conversion factor in cardinality()
	 * rather than bucket-specific compensation.
	 * @return null (no compensation factors used)
	 */
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Array storing maximum leading zero counts for each bucket */
	private final byte[] maxArray;

}
