package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicLogLog cardinality estimator extending DemiLogLog (LogLog16) with an adaptive
 * early-exit mechanism. Maintains a LogLogLog meta-structure: an int[64] histogram tracking
 * how many buckets hold each leading-zero tier. This drives a dynamic minimum threshold
 * (minZeros) that starts at zero and rises with cardinality, achieving 99.9%+ early exit
 * for large datasets without any loss of accuracyâ€”skipped elements can never improve a
 * bucket that already holds a higher score.
 *
 * @author Brian Bushnell
 * @contributor Chloe
 * @date February 27, 2026
 */
public final class DynamicLogLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a DynamicLogLog with default parameters.
	 * Uses 2048 buckets, k=31, random seed, and no minimum probability filtering. */
	DynamicLogLog(){
		this(2048, 31, -1, 0);
	}

	/** Creates a DynamicLogLog with parameters parsed from command-line arguments.
	 * @param p Parser containing configuration from command-line flags */
	DynamicLogLog(Parser p){
		super(p);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		logloglog=new int[wordlen];
		logloglog[0]=buckets;
		minZeroCount=buckets;
	}

	/**
	 * Creates a DynamicLogLog with specified parameters.
	 * @param buckets_ Number of buckets (counters) for the hash table
	 * @param k_ K-mer length for sequence hashing
	 * @param seed Random number generator seed; -1 for random seed
	 * @param minProb_ Ignore k-mers with under this probability of being correct
	 */
	DynamicLogLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		logloglog=new int[wordlen];
		logloglog[0]=buckets;
		minZeroCount=buckets;
	}

	@Override
	public DynamicLogLog copy(){return new DynamicLogLog(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Restores floating-point compressed score to approximate original hash magnitude.
	 * @param score Compressed 16-bit value from maxArray
	 * @return Approximate original hash value before compression */
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
		final int div=Tools.max(count, 1);
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

		final double proxy=(USE_MEAN ? mean : USE_MEDIAN ? median : USE_MWA ? mwa : USE_HMEAN ? hmean : USE_GMEAN ? gmean : mean);

		final double estimatePerSet=2*(Long.MAX_VALUE/proxy)*SKIPMOD;
		final double total=estimatePerSet*div*((count+buckets)/(float)(buckets+buckets));

		final double estSum=div*Math.exp(estLogSum/(Tools.max(div, 1)));
		double medianEst=2*(Long.MAX_VALUE/(double)median)*SKIPMOD*div;

		long cardinality=Math.min(added, (long)(total));
		lastCardinality=cardinality;
		return cardinality;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicLogLog)log);
	}

	public void add(DynamicLogLog log){
		added+=log.added;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				final char maxA=maxArray[i], maxB=log.maxArray[i];
				final char countA=countArray[i], countB=log.countArray[i];
				maxArray[i]=Tools.max(maxA, maxB);
				if(maxA==maxB){countArray[i]=(char)Tools.min(countA+(int)countB, Character.MAX_VALUE);}
				else{countArray[i]=(maxA>maxB ? countA : countB);}
			}
			rebuildLLL();
		}
	}

	/** Rebuilds the lll histogram and minZeros from current maxArray contents.
	 * Required after a merge since bucket values may change non-monotonically. */
	private void rebuildLLL(){
		java.util.Arrays.fill(logloglog, 0);
		for(int i=0; i<buckets; i++){logloglog[maxArray[i]>>mantissabits]++;}
		minZeros=0;
		while(logloglog[minZeros]==0){minZeros++;}
	}

	/** Compares two DynamicLogLog instances bucket by bucket.
	 * @return int[]{lower, equal, higher} bucket counts */
	public int[] compareTo(DynamicLogLog blog){return compare(maxArray, blog.maxArray);}

	/** Compares two maxArrays bucket by bucket using branchless bit arithmetic.
	 * @return int[]{lower, equal, higher} bucket counts */
	public static int[] compare(char[] arrayA, char[] arrayB){
		int lower=0, equal=0, higher=0;
		for(int i=0; i<arrayA.length; i++){
			final int a=arrayA[i], b=arrayB[i];
			int dif=a-b;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher};
	}

	/**
	 * Hashes a number and updates the appropriate bucket counter.
	 * Dynamic early exit: if leading-zero count is below minZeros, the element
	 * cannot improve any bucket and is discarded immediatelyâ€”pure register arithmetic,
	 * always-taken branch, no memory access. When a bucket's leading-zero tier increases,
	 * lll is updated in O(1); if the current minimum tier empties, minZeros advances.
	 *
	 * @param number The value to hash and track
	 */
	public void hashAndStore0(final long number){
		final long key=Tools.hash64shift(number);
		final int nlz=Long.numberOfLeadingZeros(key);
		if(nlz<minZeros){return;}

		final int shift=offset-nlz;
		final int score=(nlz<<mantissabits)+(int)((~(key>>>shift))&mask);
		final int bucket=(int)(key&bucketMask);

		final int oldValue=maxArray[bucket];
		final int newValue=Math.max(score, oldValue);
		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
		maxArray[bucket]=(char)newValue;

		final int nlzOld;
		if(newValue>oldValue && (nlzOld=(oldValue>>mantissabits))<nlz){
			logloglog[nlzOld]--;
			logloglog[nlz]++;
			while(logloglog[minZeros]==0){minZeros++;}
		}

		final char count=countArray[bucket];
		countArray[bucket]=(char)(oldValue>score ? count :
			oldValue==score ? Math.max(count, (char)(count+1)) : 1);
	}

	@Override
	public void hashAndStore(final long number){
		final long key=Tools.hash64shift(number);
		final int nlz=Long.numberOfLeadingZeros(key);

		// Global early exit
		if(nlz<minZeros){return;}

		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		final int score=(nlz<<mantissabits)+(int)((~(key>>>shift))&mask);
		final int oldValue=maxArray[bucket];
		if(score<oldValue) {return;}//Optional early exit.
		
		final int newValue=Math.max(score, oldValue);
		final int nlzOld=(oldValue>>mantissabits);
		
		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
		maxArray[bucket]=(char)newValue;
		final char count=countArray[bucket];
		countArray[bucket]=(char)(oldValue>score ? count : 
			oldValue==score ? Math.max(count, (char)(count+1)) : 1);
		
		if(nlz>nlzOld && nlzOld==minZeros) {//Promotion
//			minZeroCount--;
//			assert(minZeroCount>=0 && minZeros<=64) : minZeroCount+", "+minZeros;
//			while(minZeroCount<1) {
//				minZeros++;
//				minZeroCount=countTermsInTier(minZeros, maxArray);
//				assert(minZeroCount>=0 && minZeros<=64) : minZeroCount+", "+minZeros;
//			}
//			assert(minZeroCount>0 && minZeroCount<=buckets) : minZeroCount+", "+minZeros;
			
			minZeroCount--;
			while(minZeroCount==0 && minZeros<wordlen) {
				minZeros++;
				minZeroCount=countTermsInTier(minZeros, maxArray);
			}
			
//			if(--minZeroCount<1){
//				do{
//					minZeros++;
//					minZeroCount=countTermsInTier(minZeros, maxArray);
//				}while(minZeroCount==0 && minZeros<wordlen);
//			}
		}
		
//		// Update the histogram and global threshold
//		if(nlzOld<nlz){
//			logloglog[nlzOld]--;
//			logloglog[nlz]++;
//			while(logloglog[minZeros]==0){minZeros++;}
//		}
		
	}
	
	private static int countTermsInTier(final int nlz, final char[] array) {
		int sum=0;
		for(final char c : array){
			final int diff=(c>>mantissabits)^nlz;  // 0 iff tier matches
			final int equalsBit=(~(diff|-diff))>>>31;  // 1 if diff==0, else 0
			sum+=equalsBit;
		}
		return sum;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public char[] counts16(){return countArray;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compressed 16-bit bucket maxima; 6-bit leading-zero count + 10-bit mantissa. */
	private final char[] maxArray;
	/** Count of observations at the current maximum score for each bucket. */
	private final char[] countArray;
	/**
	 * LogLogLog meta-structure: logloglog[k] = number of buckets whose current maximum
	 * has k leading zeros. Initialized with logloglog[0]=buckets; sum always equals buckets.
	 * Updated in O(1) amortized time per hashAndStore; minZeros advances at most 64
	 * times over the entire lifetime of the structure.
	 */
	private final int[] logloglog;
	/** Minimum leading-zero tier held by any bucket; dynamic early-exit threshold.
	 * Rises monotonically from 0 toward 63 as cardinality grows. */
	private int minZeros=0;
	/** Number of buckets with minZeros */
	private int minZeroCount;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Number of mantissa bits for floating-point compression; 10 is maximum. */
	private static final int mantissabits=10;
	private static final int mask=(1<<mantissabits)-1;
	private static final int offset=wordlen-mantissabits-1;

}
