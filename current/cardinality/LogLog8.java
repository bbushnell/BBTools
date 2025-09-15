package cardinality;

import shared.Parser;
import shared.Tools;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Mar 6, 2020
 *
 */
public final class LogLog8 extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create a LogLog with default parameters */
	LogLog8(){
		this(2048, 31, -1, 0);
	}
	
	/** Create a LogLog with parsed parameters */
	LogLog8(Parser p){
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
	LogLog8(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
	}
	
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
	
	/** Restores floating point to integer */
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
		
		long cardinality=(long)(total);
		lastCardinality=cardinality;
		return cardinality;
	}
	
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog8)log);
	}
	
	public void add(LogLog8 log){
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
