package cardinality;

import java.util.concurrent.atomic.AtomicLongArray;

import shared.Parser;
import shared.Tools;
import structures.LongList;

/**
 * LogLog-based cardinality estimator that tracks the maximum hash value per bucket.
 * Uses a probabilistic approach to estimate the number of unique k-mers in large datasets
 * with constant memory usage proportional to the number of buckets.
 *
 * @author Brian Bushnell
 * @date Feb 20, 2020
 */
public final class BBLog extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a LogLog estimator with default parameters.
	 * Uses 2048 buckets, k-mer length 31, random seed, and no minimum probability filter. */
	BBLog(){
		this(2048, 31, -1, 0);
	}
	
	/**
	 * Creates a LogLog estimator with parameters parsed from command-line arguments.
	 * Initializes appropriate storage arrays based on atomic and count tracking settings.
	 * @param p Parser containing configuration parameters
	 */
	BBLog(Parser p){
		super(p);
		maxArrayA=(atomic ? new AtomicLongArray(buckets) : null);
		maxArray=(atomic ? null : new long[buckets]);
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	/**
	 * Creates a LogLog estimator with explicitly specified parameters.
	 *
	 * @param buckets_ Number of buckets for hash partitioning
	 * @param k_ K-mer length for sequence processing
	 * @param seed Random number generator seed; -1 for random seed
	 * @param minProb_ Minimum probability threshold for k-mer inclusion
	 */
	BBLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArrayA=(atomic ? new AtomicLongArray(buckets) : null);
		maxArray=(atomic ? null : new long[buckets]);
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	@Override
	public BBLog copy() {return new BBLog(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Estimates the total cardinality using LogLog algorithm.
	 * Computes multiple estimates from bucket max values: arithmetic mean-based,
	 * geometric mean-based, and median-based. Returns the arithmetic mean estimate
	 * which typically provides the best variance characteristics.
	 *
	 * @return Estimated number of unique elements
	 */
	@Override
	public final long cardinality(){
		double difSum=0;
		double estLogSum=0;
		int count=0;
		LongList list=new LongList(buckets);
		//assert(atomic);
		if(atomic){
			for(int i=0; i<maxArrayA.length(); i++){
				long val=maxArrayA.get(i);
				if(val>0){
//					System.err.println("val="+val);
					long dif=Long.MAX_VALUE-val;
					difSum+=dif;
					count++;
					double est=2*(Long.MAX_VALUE/(double)dif)*SKIPMOD;
					estLogSum+=Math.log(est);
					list.add(dif);
				}
			}
		}else{
			for(int i=0; i<maxArray.length; i++){
				long val=maxArray[i];
				if(val>0){
					long dif=Long.MAX_VALUE-val;
					difSum+=dif;
					count++;
					double est=2*(Long.MAX_VALUE/(double)dif)*SKIPMOD;
					estLogSum+=Math.log(est);
					list.add(dif);
				}
			}
		}
		int div=count;//Could also be count be that causes problems
		final double mean=difSum/Tools.max(div, 1);
		final double estimatePerSet=2*(Long.MAX_VALUE/mean)*SKIPMOD;
		final double total=estimatePerSet*div*((count+buckets)/(float)(buckets+buckets));

		final double estSum=div*Math.exp(estLogSum/(Tools.max(div, 1)));
		list.sort();
		long median=list.median();
		double medianEst=2*(Long.MAX_VALUE/(double)median)*SKIPMOD*div;
		
//		new Exception().printStackTrace();
		
//		System.err.println(maxArray);
//		//Overall, it looks like "total" is the best, then "estSum", then "medianEst" is the worst, in terms of variance.
//		System.err.println("difSum="+difSum+", count="+count+", mean="+mean+", est="+estimatePerSet+", total="+(long)total);
//		System.err.println("estSum="+(long)estSum+", median="+median+", medianEst="+(long)medianEst);
		
		long cardinality=Math.min(added, (long)(total));
		lastCardinality=cardinality;
		return cardinality;
	}

	/** Returns the count array tracking occurrences of maximum values per bucket.
	 * @return Array of counts, or null if count tracking is disabled */
	@Override
	public int[] getCounts(){
		return counts;
	}
	
	/**
	 * Merges another cardinality tracker into this one.
	 * Casts to BBLog and delegates to the typed add method.
	 * @param log The cardinality tracker to merge
	 */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((BBLog)log);
	}
	
	public void add(BBLog log){
		added+=log.added;
		if(atomic && maxArrayA!=log.maxArrayA){
			for(int i=0; i<buckets; i++){
				maxArrayA.set(i, Tools.max(maxArrayA.get(i), log.maxArrayA.get(i)));
			}
		}else 
		if(maxArray!=log.maxArray){
			if(counts==null){
				for(int i=0; i<buckets; i++){
					maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
				}
			}else{
				for(int i=0; i<buckets; i++){
					final long a=maxArray[i], b=log.maxArray[i];
					if(a==b){
						counts[i]+=log.counts[i];
					}else if(b>a){
						maxArray[i]=b;
						counts[i]=log.counts[i];
					}
				}
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
		
		if(atomic){
			long x=maxArrayA.get(bucket);
			while(key>x){
				boolean b=maxArrayA.compareAndSet(bucket, x, key);
				if(b){x=key;}
				else{x=maxArrayA.get(bucket);}
			}
		}else{
			if(trackCounts){
				if(key>maxArray[bucket]){
					maxArray[bucket]=key;
					counts[bucket]=1;
				}else if(key==maxArray[bucket]){
					counts[bucket]++;
				}
			}else{
				maxArray[bucket]=Tools.max(key, maxArray[bucket]);
			}
		}
	}
	
	/**
	 * Returns compensation factors for LogLog estimation accuracy.
	 * This implementation returns null as compensation is not used.
	 * @return null (no compensation factors)
	 */
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Non-atomic array storing maximum hash values per bucket */
	private final long[] maxArray;
	/** Atomic array storing maximum hash values per bucket for concurrent access */
	private final AtomicLongArray maxArrayA;
	/** Array tracking occurrence counts of maximum values per bucket */
	private final int[] counts;
	
//	private static long minKey=(long)(0.75f*Long.MAX_VALUE); //non-atomic 15% faster without this
	
}
