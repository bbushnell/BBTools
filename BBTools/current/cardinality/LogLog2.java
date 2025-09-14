package cardinality;

import java.util.concurrent.atomic.AtomicIntegerArray;

import shared.Parser;
import shared.Tools;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Feb 20, 2020
 *
 */
public final class LogLog2 extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create a LogLog with default parameters */
	LogLog2(){
		this(2048, 31, -1, 0);
	}
	
	/** Create a LogLog with parsed parameters */
	LogLog2(Parser p){
		super(p);
		//assert(atomic);
		maxArrayA=(atomic ? new AtomicIntegerArray(buckets) : null);
		maxArray=(atomic ? null : new int[buckets]);
	}
	
	/**
	 * Create a LogLog with specified parameters
	 * @param buckets_ Number of buckets (counters)
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	LogLog2(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		//assert(atomic);
		maxArrayA=(atomic ? new AtomicIntegerArray(buckets) : null);
		maxArray=(atomic ? null : new int[buckets]);
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
		//assert(atomic);
		if(atomic){
			for(int i=0; i<maxArrayA.length(); i++){
				int max=maxArrayA.get(i);
				long val=restore(max);
				if(max>0 && val>0){
//					long val=restore(max);
//					System.err.println("val="+val);
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
		}else{
			for(int i=0; i<maxArray.length; i++){
				int max=maxArray[i];
				long val=restore(max);
				if(max>0 && val>0){
//					long val=restore(max);
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
//		double mean=sum/((1<<mantissabits)*(double)buckets);
//		double correction=0.56326183361037098678934414274035;//0.56403894240204307426602541326855;
//		//Better: //0.56326183361037098678934414274035
//		long cardinality=(long)((((Math.pow(2, mean)-1)*buckets*SKIPMOD))*correction);
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
	
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog2)log);
	}
	
	public void add(LogLog2 log){
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
		
//		if(leading<3){return;}//Speeds up by 20%, even more at 4.  Slows at 2.
		
		int shift=wordlen-leading-mantissabits-1;
		
		int score=(leading<<mantissabits)+(int)((~(key>>>shift))&mask);
//		assert(false) : "\n"+Long.toBinaryString(key)+", leading="+leading+", shift="+shift+"\n"+Long.toBinaryString(score);
		
		//+"\n"+score+"\n"+restore(score);
		
//		final int bucket=(int)((number&Integer.MAX_VALUE)%buckets);
		final int bucket=(int)(key&bucketMask);
		
		if(atomic){
			int x=maxArrayA.get(bucket);
			while(score>x){
//				System.err.println("\n"+Long.toBinaryString(key)+", leading="+leading+", score="+score+", x="+x+"\n"+Long.toBinaryString(score));
//				System.err.println("\n"+Long.toBinaryString(restore(score)));
				boolean b=maxArrayA.compareAndSet(bucket, x, score);
				if(b){x=score;}
				else{x=maxArrayA.get(bucket);}
//				assert(leading<9);
			}
		}else{
			maxArray[bucket]=Tools.max(score, maxArray[bucket]);
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
	private final int[] maxArray;
	/** Atomic version of maxArray. */
	private final AtomicIntegerArray maxArrayA;
	
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
	private static int mantissabits=20;
	private static int mask=(1<<mantissabits)-1;
//	private static final int shift=wordlen-leading-mantissabits-1;
	
}
