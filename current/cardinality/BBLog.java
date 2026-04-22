package cardinality;

import parse.Parser;
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
	
	/** Default constructor: 2048 buckets, k=31. */
	BBLog(){this(2048, 31, -1, 0);}
	
	/** Construct from parsed command-line arguments. */
	BBLog(Parser p){
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
	BBLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new long[buckets];
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	/** Create an independent copy with a fresh seed. */
	@Override
	public BBLog copy(){return new BBLog(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Estimates cardinality from bucket maxima using arithmetic-mean LogLog.
	 * Also computes geometric-mean and median estimates (unused but retained
	 * for comparison); returns the arithmetic-mean estimate.
	 * @return Estimated number of unique elements
	 */
	@Override
	public final long cardinality(){
		double difSum=0;
		double estLogSum=0;
		int count=0;
		final LongList list=new LongList(buckets);
		for(int i=0; i<maxArray.length; i++){
			final long val=maxArray[i];
			if(val>0){
				final long dif=Long.MAX_VALUE-val;
				difSum+=dif;
				count++;
				final double est=2*(Long.MAX_VALUE/(double)dif);
				estLogSum+=Math.log(est);
				list.add(dif);
			}
		}
		final int div=count;
		final double mean=difSum/Tools.max(div, 1);
		final double estimatePerSet=2*(Long.MAX_VALUE/mean);
		final double total=estimatePerSet*div*((count+buckets)/(float)(buckets+buckets));

		final double estSum=div*Math.exp(estLogSum/(Tools.max(div, 1)));
		list.sort();
		final long median=list.median();
		final double medianEst=2*(Long.MAX_VALUE/(double)median)*div;

		final long cardinality=Math.min(added, (long)(total));
		lastCardinalityStatic=cardinality;
		return cardinality;
	}

	/** Returns per-bucket count array, or null if count tracking is disabled. */
	@Override
	public int[] getCounts(){return counts;}
	
	/** Merge another tracker into this one (must be BBLog). */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((BBLog)log);
	}

	/** Merge another BBLog by taking per-bucket max hash values. */
	public void add(BBLog log){
		added+=log.added;
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
	
	/** Hash a value and store the maximum hash per bucket. */
	@Override
	public void hashAndStore(final long number){
		final long key=Tools.hash64shift(number);
		final int bucket=(int)(key&bucketMask);
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
