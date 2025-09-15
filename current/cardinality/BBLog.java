package cardinality;

import java.util.concurrent.atomic.AtomicLongArray;

import shared.Parser;
import shared.Tools;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Feb 20, 2020
 *
 */
public final class BBLog extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create a LogLog with default parameters */
	BBLog(){
		this(2048, 31, -1, 0);
	}
	
	/** Create a LogLog with parsed parameters */
	BBLog(Parser p){
		super(p);
		maxArrayA=(atomic ? new AtomicLongArray(buckets) : null);
		maxArray=(atomic ? null : new long[buckets]);
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	/**
	 * Create a LogLog with specified parameters
	 * @param buckets_ Number of buckets (counters)
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	BBLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArrayA=(atomic ? new AtomicLongArray(buckets) : null);
		maxArray=(atomic ? null : new long[buckets]);
		counts=(trackCounts ? new int[buckets] : null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
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
		add((BBLog)log);
	}
	
	public void add(BBLog log){
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
	
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maintains state.  These are the actual buckets. */
	private final long[] maxArray;
	/** Atomic version of maxArray. */
	private final AtomicLongArray maxArrayA;
	/** Counts associated with values in the buckets. */
	private final int[] counts;
	
//	private static long minKey=(long)(0.75f*Long.MAX_VALUE); //non-atomic 15% faster without this
	
}
