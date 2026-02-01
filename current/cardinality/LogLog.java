package cardinality;

import java.util.concurrent.atomic.AtomicIntegerArray;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * Basic LogLog cardinality estimator for counting distinct k-mers or sequences.
 * Uses leading zero counts from hashed values to estimate cardinality.
 * Less accurate than other variants but simple and fast.
 *
 * @author Brian Bushnell
 * @date Sep 30, 2015
 */
public final class LogLog extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a LogLog with default parameters: 2048 buckets, k=31, random seed, minProb=0
	 */
	LogLog(){
		this(2048, 31, -1, 0);
	}
	
	/** Creates a LogLog with parameters from command-line parser.
	 * @param p Parser containing loglog configuration parameters */
	LogLog(Parser p){
		super(p);
		//assert(atomic);
		maxArrayA=(atomic ? new AtomicIntegerArray(buckets) : null);
		maxArray=(atomic ? null : new int[buckets]);
	}
	
	/**
	 * Creates a LogLog with specified parameters.
	 *
	 * @param buckets_ Number of buckets (counters) for tracking
	 * @param k_ K-mer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore k-mers with under this probability of being correct
	 */
	LogLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		//assert(atomic);
		maxArrayA=(atomic ? new AtomicIntegerArray(buckets) : null);
		maxArray=(atomic ? null : new int[buckets]);
	}
	
	@Override
	public LogLog copy() {return new LogLog(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	//Restores floating point to integer
	private long restore(int score){
		int leading=score;
		long mantissa=1;
		int shift=64-leading-1;
		long original=mantissa<<shift;
		return original;
	}
	
	@Override
	public final long cardinality(){
		double difSum=0;
		double estLogSum=0;
		int count=0;
		LongList list=new LongList(buckets);
		//assert(atomic);
		if(atomic){
			for(int i=0; i<maxArrayA.length(); i++){
				int max=maxArrayA.get(i);
				long val=restore(max);
				if(max>0 && val>0){
//					long val=restore(max);
//					System.err.println("val="+val);
					final long dif=val;
					difSum+=dif;
					count++;
					double est=2*(Long.MAX_VALUE/(double)dif)*SKIPMOD;
					estLogSum+=Math.log(est);
					list.add(dif);
				}
			}
		}else{
			for(int i=0; i<maxArray.length; i++){
				int max=maxArray[i];
				long val=restore(max);
				if(max>0 && val>0){
//					long val=restore(max);
					final long dif=val;
					difSum+=dif;
					count++;
					double est=2*(Long.MAX_VALUE/(double)dif)*SKIPMOD;
					estLogSum+=Math.log(est);
					list.add(dif);
				}
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
	
//	@Override
//	public final long cardinality(){
//		long sum=0;
//		//assert(atomic);
//		if(atomic){
//			for(int i=0; i<maxArray.length(); i++){
//				sum+=maxArray.get(i);
//			}
//		}else{
//			for(int i=0; i<maxArray2.length; i++){
//				sum+=maxArray2[i];
//			}
//		}
//		double mean=sum/(double)buckets;
//		long cardinality=(long)((((Math.pow(2, mean)-1)*buckets*SKIPMOD))/1.258275);
//		lastCardinality=cardinality;
//		return cardinality;
//	}
	
	public final long cardinalityH(){
		double sum=0;
		for(int i=0; i<maxArrayA.length(); i++){
			int x=Tools.max(1, maxArrayA.get(i));
			sum+=1.0/x;
		}
		double mean=buckets/sum;
		return (long)((Math.pow(2, mean)*buckets*SKIPMOD));
	}
	
	/**
	 * Adds another CardinalityTracker to this LogLog by merging bucket values.
	 * Takes the maximum leading zero count for each bucket position.
	 * @param log CardinalityTracker to merge (must be LogLog instance)
	 */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog)log);
	}
	
	public void add(LogLog log){
		added+=log.added;
		if(atomic && maxArrayA!=log.maxArrayA){
			for(int i=0; i<buckets; i++){
				maxArrayA.set(i, Tools.max(maxArrayA.get(i), log.maxArrayA.get(i)));
			}
		}else if(maxArray!=log.maxArray){
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
		int leading=Long.numberOfLeadingZeros(key);
		
//		counts[leading]++;
		
//		if(leading<3){return;}//Slows things down slightly
//		final int bucket=(int)((number&Integer.MAX_VALUE)%buckets);
		final int bucket=(int)(key&bucketMask);
		
		if(atomic){
			int x=maxArrayA.get(bucket);
			while(leading>x){
				boolean b=maxArrayA.compareAndSet(bucket, x, leading);
				if(b){x=leading;}
				else{x=maxArrayA.get(bucket);}
			}
		}else{
			maxArray[bucket]=Tools.max(leading, maxArray[bucket]);
		}
	}
	
	/** Returns the compensation factor array for LogLog bucket count adjustments */
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return compensationFactorLogBucketsArray;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] maxArray;
	/**
	 * Atomic array storing maximum leading zero counts for each bucket (atomic mode)
	 */
	private final AtomicIntegerArray maxArrayA;
	
	private static final float[] compensationFactorLogBucketsArray={
			0.053699781f, 0.49556874f, 0.742263622f, 0.861204899f, 0.926038294f,
			0.967001269f, 0.982949748f, 0.992495155f, 0.996185775f, 0.998077246f
		};

}
