package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * LogLog cardinality estimator using 16-bit counters with 10-bit mantissa compression.
 * Implements a probabilistic data structure for approximating set cardinality with
 * controlled memory usage and precision tradeoffs. Uses floating-point compression
 * to store larger values in 16-bit counters while maintaining reasonable accuracy.
 *
 * @author Brian Bushnell
 * @date Mar 6, 2020
 */
public final class LogLog16 extends CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a LogLog16 with default parameters.
	 * Uses 2048 buckets, k=31, random seed, and no minimum probability filtering. */
	LogLog16(){
		this(2048, 31, -1, 0);
	}
	
	/** Creates a LogLog16 with parameters parsed from command-line arguments.
	 * @param p Parser containing configuration from command-line flags */
	LogLog16(Parser p){
		super(p);
		maxArray=new char[buckets];
		countArray=new char[buckets];
	}
	
	/**
	 * Creates a LogLog16 with specified parameters.
	 *
	 * @param buckets_ Number of buckets (counters) for the hash table
	 * @param k_ K-mer length for sequence hashing
	 * @param seed Random number generator seed; -1 for random seed
	 * @param minProb_ Ignore k-mers with under this probability of being correct
	 */
	LogLog16(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
		countArray=new char[buckets];
	}
	
	@Override
	public LogLog16 copy(){return new LogLog16(buckets, k, -1, minProb);}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Restores an approximate original value from a compressed 16-bit score.
	 * Extracts leading-zero count and mantissa bits, then reconstructs the value. */
	private long restore(final int score){
		final long lowbits=(~score)&mask;
		final int leading=(int)(score>>>mantissabits);
		final long mantissa=(1L<<mantissabits)|lowbits;
		final int shift=wordlen-leading-mantissabits-1;
		final long original=mantissa<<shift;
		return original;
	}
	
	/**
	 * Estimates cardinality using restored 16-bit bucket values.
	 * Supports multiple proxy statistics (mean, harmonic, geometric, median, MWA)
	 * selected via static flags in CardinalityTracker.
	 */
	@Override
	public final long cardinality(){
		double difSum=0;
		double hSum=0;
		double gSum=0;
		double rSum=0;
		double estLogSum=0;
		int count=0;
		final LongList list=new LongList(buckets);

		for(int i=0; i<maxArray.length; i++){
			final int max=maxArray[i];
			final long val=restore(max);
			if(max>0 && val>0){
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
		double medianEst=2*(Long.MAX_VALUE/(double)median)*div;

//		new Exception().printStackTrace();

//		System.err.println(maxArray);
////		Overall, it looks like "total" is the best, then "estSum", then "medianEst" is the worst, in terms of variance.
//		System.err.println("difSum="+difSum+", count="+count+", mean="+mean+", est="+estimatePerSet+", total="+(long)total);
//		System.err.println("estSum="+(long)estSum+", median="+median+", medianEst="+(long)medianEst);

		final long cardinality=Math.min(added, (long)(total));
		lastCardinalityStatic=cardinality;
		return cardinality;
	}
	
	/** Merges another tracker into this one by taking per-bucket maximums.
	 * @param log Tracker to merge (must be LogLog16) */
	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog16)log);
	}

	/** Merges another LogLog16 by element-wise maximum of bucket arrays.
	 * When buckets tie, counts are summed; otherwise the winner's count is kept. */
	public void add(LogLog16 log){
		added+=log.added;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				final char maxA=maxArray[i], maxB=log.maxArray[i];
				final char countA=countArray[i], countB=log.countArray[i];
				maxArray[i]=Tools.max(maxA, maxB);
				if(maxA==maxB){countArray[i]=(char)Tools.min(countA+(int)countB, Character.MAX_VALUE);}
				else{countArray[i]=(maxA>maxB ? countA : countB);}
			}
		}
	}
	
	/** Compares bucket arrays, returning {lower, equal, higher} counts. */
	public int[] compareTo(LogLog16 blog){
		return compare(maxArray, blog.maxArray);
	}

	/** Counts how many buckets in arrayA are lower, equal, or higher than arrayB. */
	public static int[] compare(char[] arrayA, char[] arrayB){
		int lower=0, equal=0, higher=0;
		for(int i=0; i<arrayA.length; i++){
			final int a=arrayA[i], b=arrayB[i];
			final int dif=a-b;
			final int nbit=(dif>>>31);
			final int hbit=((-dif)>>>31);
			final int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher};
	}
	
	/** {@inheritDoc} */
	@Override
	public void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
//		if(key<0 || key>0x07FFFFFFFFFFFFFFL){return;}//Super slow!
		final int nlz=Long.numberOfLeadingZeros(key)&63;//mask is used to keep number in 6 bits
		if(nlz<4){return;}

//		counts[leading]++;

//		int shift=wordlen-leading-mantissabits-1;
		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		final int score=(nlz<<mantissabits)+(int)((~(key>>>shift))&mask);
		final int oldValue=maxArray[bucket];
		final int newValue=Math.max(score, oldValue);

		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
		maxArray[bucket]=(char)newValue;
		final char count=countArray[bucket];
		countArray[bucket]=(char)(oldValue>score ? count :
			oldValue==score ? Math.max(count, (char)(count+1)) : 1);
	}
	
	/** Returns null; this tracker does not use per-bucket compensation factors. */
	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/** Returns the per-bucket hit count array. */
	public char[] counts16(){return countArray;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final char[] maxArray;
	private final char[] countArray;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

//	public static void setMantissaBits(int x){
//		assert(x>=0 && x<=10);
//		assert(x+6<=16);
//		mantissabits=x;
//		mask=(1<<mantissabits)-1;
//	}

	private static final int wordlen=64;

	/** Number of mantissa bits for floating-point compression; 10 is maximum. */
	private static final int mantissabits=10;
	private static final int mask=(1<<mantissabits)-1;
	private static final int offset=wordlen-mantissabits-1;
	
}
