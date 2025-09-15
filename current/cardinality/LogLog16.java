package cardinality;

import shared.Parser;
import shared.Tools;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Mar 6, 2020
 *
 */
public final class LogLog16 extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create a LogLog with default parameters */
	LogLog16(){
		this(2048, 31, -1, 0);
	}
	
	/** Create a LogLog with parsed parameters */
	LogLog16(Parser p){
		super(p);
		maxArray=new char[buckets];
	}
	
	/**
	 * Create a LogLog with specified parameters
	 * @param buckets_ Number of buckets (counters)
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	LogLog16(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	//Restores floating point to integer
	private long restore(int score){
		long lowbits=(~score)&mask;
		int leading=(int)(score>>>mantissabits);
		long mantissa=(1L<<mantissabits)|lowbits;
		int shift=wordlen-leading-mantissabits-1;
		long original=mantissa<<shift;
		return original;
	}
	
	@Override
	public final long cardinality(){
		double difSum=0;
		double hSum=0;
		double gSum=0;
		double rSum=0;
		double estLogSum=0;
		int count=0;
		LongList list=new LongList(buckets);
		
		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
			if(max>0 && val>0){
				long dif=val;
				difSum+=dif;
				hSum+=1.0/Tools.max(1, dif);
				gSum+=Math.log(Tools.max(1, dif));
				rSum+=Math.sqrt(dif);
				count++;
				double est=2*(Long.MAX_VALUE/(double)dif)*SKIPMOD;
				estLogSum+=Math.log(est);
				list.add(dif);
			}
		}
		final int div=Tools.max(count, 1);//Could be count or buckets but one causes problems
		final double mean=difSum/div;
		double hmean=hSum/div;
		double gmean=gSum/div;
		double rmean=rSum/div;
		hmean=1.0/hmean;
		gmean=Math.exp(gmean);
		rmean=rmean*rmean;
		list.sort();
		final long median=list.median();
		final double mwa=list.medianWeightedAverage();
		
		//What to use as the value from the counters 
		final double proxy=(USE_MEAN ? mean : USE_MEDIAN ? median : USE_MWA ? mwa : USE_HMEAN ? hmean : USE_GMEAN ? gmean : mean);
		
		final double estimatePerSet=2*(Long.MAX_VALUE/proxy)*SKIPMOD;
		final double total=estimatePerSet*div*((count+buckets)/(float)(buckets+buckets));

		final double estSum=div*Math.exp(estLogSum/(Tools.max(div, 1)));
		double medianEst=2*(Long.MAX_VALUE/(double)median)*SKIPMOD*div;
		
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
		add((LogLog16)log);
	}
	
	public void add(LogLog16 log){
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
		int leading=Long.numberOfLeadingZeros(key)&63;//mask is used to keep number in 6 bits 
		
//		counts[leading]++;
		
//		if(leading<3){return;}//Speeds up by 20%, even more at 4.  Slows at 2.
		
		int shift=wordlen-leading-mantissabits-1;
		
		int score=(leading<<mantissabits)+(int)((~(key>>>shift))&mask);
//		assert(false) : "\n"+Long.toBinaryString(key)+", leading="+leading+", shift="+shift+"\n"+Long.toBinaryString(score);
		
		//+"\n"+score+"\n"+restore(score);
		
//		final int bucket=(int)((number&Integer.MAX_VALUE)%buckets);
		final int bucket=(int)(key&bucketMask);
		
		int newValue=Tools.max(score, maxArray[bucket]);
		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
		maxArray[bucket]=(char)newValue;
	}
	
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final char[] maxArray;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void setMantissaBits(int x){
		assert(x>=0 && x<25);
		assert(x+6<32);
		mantissabits=x;
		mask=(1<<mantissabits)-1;
	}

	private static final int wordlen=64;
	
	/** Precision or mantissa bits.
	 * This should not be changed.  As long as it is >10 the result will be accurate.
	 * At low values like 2 the cardinality estimate becomes too high due to a loss of precision,
	 * and would need a fixed multiplier.  
	 */
	private static int mantissabits=10;//10 is the max possible
	private static int mask=(1<<mantissabits)-1;
	
}
