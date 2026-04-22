package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * LogLog cardinality estimation with floating-point compression.
 * Stores leading-zero counts plus mantissa bits for improved precision
 * over basic LogLog. Configurable mantissa width (default 20 bits).
 *
 * @author Brian Bushnell
 * @date Feb 20, 2020
 */
public final class LogLog2 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	LogLog2(){
		this(2048, 31, -1, 0);
	}

	/** Construct from parsed command-line arguments. */
	LogLog2(Parser p){
		super(p);
		maxArray=new int[buckets];
	}

	/**
	 * Full constructor.
	 * @param buckets_ Number of buckets (counters) for hash distribution
	 * @param k_ K-mer length
	 * @param seed Random number generator seed; -1 for random seed
	 * @param minProb_ Minimum probability threshold for k-mer inclusion
	 */
	LogLog2(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new int[buckets];
	}

	@Override
	public LogLog2 copy(){return new LogLog2(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Restores a compressed score back to the original hash magnitude. */
	private long restore(int score){
		final long lowbits=(~score)&mask;
		final int leading=(int)(score>>>mantissabits);
		final long mantissa=(1L<<mantissabits)|lowbits;
		final int shift=wordlen-leading-mantissabits-1;
		final long original=mantissa<<shift;
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
		final LongList list=new LongList(buckets);
		{
			for(int i=0; i<maxArray.length; i++){
				final int max=maxArray[i];
				final long val=restore(max);
				if(max>0 && val>0){
//					long val=restore(max);
					final long dif=val;
					difSum+=dif;
					hSum+=1.0/Tools.max(1, dif);
					gSum+=Math.log(Tools.max(1, dif));
					rSum+=Math.sqrt(dif);
					count++;
					final double est=2*(Long.MAX_VALUE/(double)dif);
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

		final double estimatePerSet=2*(Long.MAX_VALUE/proxy);
		final double total=estimatePerSet*div*((count+buckets)/(float)(buckets+buckets));

		final double estSum=div*Math.exp(estLogSum/(Tools.max(div, 1)));
		final double medianEst=2*(Long.MAX_VALUE/(double)median)*div;

//		new Exception().printStackTrace();

//		System.err.println(maxArray);
////		Overall, it looks like "total" is the best, then "estSum", then "medianEst" is the worst, in terms of variance.
//		System.err.println("difSum="+difSum+", count="+count+", mean="+mean+", est="+estimatePerSet+", total="+(long)total);
//		System.err.println("estSum="+(long)estSum+", median="+median+", medianEst="+(long)medianEst);

		final long cardinality=Math.min(added, (long)(total));
		lastCardinalityStatic=cardinality;
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

	/** Harmonic-mean cardinality estimate (alternative estimator). */
	public final long cardinalityH(){
		double sum=0;
		for(int i=0; i<maxArray.length; i++){
			final int x=Tools.max(1, maxArray[i]);
			sum+=1.0/x;
		}
		final double mean=buckets/sum;
		return (long)((Math.pow(2, mean)*buckets));
	}

	/**
	 * Merges another CardinalityTracker into this LogLog2 via per-bucket max.
	 * @param log CardinalityTracker to merge (must be LogLog2 instance)
	 */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog2)log);
	}

	/** Merges another LogLog2 into this one via per-bucket max. */
	public void add(LogLog2 log){
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
		final int leading=Long.numberOfLeadingZeros(key);

//		counts[leading]++;

//		if(leading<3){return;}//Speeds up by 20%, even more at 4.  Slows at 2.

		final int shift=wordlen-leading-mantissabits-1;

		final int score=(leading<<mantissabits)+(int)((~(key>>>shift))&mask);
//		assert(false) : "\n"+Long.toBinaryString(key)+", leading="+leading+", shift="+shift+"\n"+Long.toBinaryString(score);

		//+"\n"+score+"\n"+restore(score);

//		final int bucket=(int)((number&Integer.MAX_VALUE)%buckets);
		final int bucket=(int)(key&bucketMask);

		maxArray[bucket]=Tools.max(score, maxArray[bucket]);
	}

	/** Returns null; LogLog2 does not use compensation factors. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){
		return null;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] maxArray;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Sets the mantissa bit width and updates the derived mask. */
	public static void setMantissaBits(int x){
		assert(x>=0 && x<25);
		assert(x+6<32);
		mantissabits=x;
		mask=(1<<mantissabits)-1;
	}

	private static final int wordlen=64;
	/** Number of mantissa bits for floating-point compression (default 20). */
	private static int mantissabits=20;
	private static int mask=(1<<mantissabits)-1;
//	private static final int shift=wordlen-leading-mantissabits-1;

}
