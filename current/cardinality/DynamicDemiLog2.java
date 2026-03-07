package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * DynamicDemiLog2: variant of DynamicDemiLog using RELATIVE NLZ storage.
 * <p>
 * Each bucket stores (absoluteNLZ - minZeros) instead of absoluteNLZ.
 * When minZeros advances, all bucket values are decremented by 1 NLZ (= 2048 in FP16
 * with 10 mantissa bits), keeping the stored range small and enabling future 4-bit packing.
 * <p>
 * Overflow: if relNlz would exceed the representable range, it is clamped.
 * For the current 6-bit NLZ (max 63 relative), overflow is essentially impossible.
 * For future DLL4 (4-bit, max 15 relative), overflow is absorbed by the CF matrix.
 * <p>
 * Key invariants:
 * - maxArray[i] == 0 means empty bucket (same sentinel as DDL)
 * - maxArray[i] > 0 stores (relNlz << mantissabits | mantissa), relNlz = absoluteNlz - minZeros
 * - restore() adds minZeros back to recover the absolute NLZ
 * - hllSum/hllSumFilled calculations apply minZeros correction in exponent
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class DynamicDemiLog2 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DynamicDemiLog2(){
		this(2048, 31, -1, 0);
	}

	DynamicDemiLog2(Parser p){
		super(p);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	DynamicDemiLog2(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
		countArray=new char[buckets];
		minZeroCount=buckets;
	}

	@Override
	public DynamicDemiLog2 copy(){return new DynamicDemiLog2(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		double difSum=0;
		double hllSumFilled=0;
		double gSum=0;
		int count=0;
		LongList list=new LongList(buckets);

		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
			if(max>0 && val>0){
				long dif=val;
				difSum+=dif;
				// Absolute NLZ = relative NLZ + minZeros
				hllSumFilled+=Math.pow(2.0, -(max>>mantissabits)-minZeros);
				gSum+=Math.log(Tools.max(1, dif));
				count++;
				list.add(dif);
			}
		}
		final double alpha_m=0.7213/(1.0+1.079/buckets);
		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		double gmean=Math.exp(gSum/div);
		list.sort();
		final long median=list.median();
		final double mwa=list.medianWeightedAverage();
		// HLL-style filled-bucket HMean
		final double hmean=(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilled);

		final double correction=(count+buckets)/(float)(buckets+buckets);
		final int V=buckets-count;

		double total;
		int cfType;
		if(USE_HMEAN && count>0){
			total=hmean;
			cfType=CorrectionFactor.HMEAN;
		}else{
			final double proxy=(USE_MEDIAN ? median : USE_MWA ? mwa : USE_GMEAN ? gmean : mean);
			cfType=(USE_MEDIAN ? CorrectionFactor.MEDCORR :
				USE_MWA ? CorrectionFactor.MWA :
				USE_GMEAN ? CorrectionFactor.GMEAN : CorrectionFactor.MEAN);
			total=2*(Long.MAX_VALUE/proxy)*div*correction;
		}
		total*=CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,cfType);

		final double lcEstimate=(V>0 ? 2.0*buckets*Math.log((double)buckets/V) : total);
		double estimate=total;
		if(USE_LC){
			final double occ=(double)count/buckets;
			final double w=1.0/(1.0+Math.exp(-LC_SHARPNESS*(occ-LC_CROSSOVER)));
			estimate=(1-w)*lcEstimate+w*total;
		}

		long cardinality=Math.min(added, (long)estimate);
		lastCardinalityStatic=lastCardinality=cardinality;
		return cardinality;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((DynamicDemiLog2)log);
	}

	/**
	 * Promotes this object's minZeros by 1: subtract 1 NLZ (1<<mantissabits) from every
	 * non-empty bucket. Buckets that would underflow become empty (0), count cleared.
	 */
	private void promote(){
		minZeros++;
		eeMask>>>=1;
		minZeroCount=0;
		filledBuckets=0;
		for(int i=0; i<buckets; i++){
			final int v=maxArray[i];
			if(v==0){continue;}
			final int promoted=v-(1<<mantissabits);
			if(promoted<=0){
				maxArray[i]=0;
				countArray[i]=0;
			}else{
				maxArray[i]=(char)promoted;
				filledBuckets++;
				if((promoted>>mantissabits)==0){minZeroCount++;}
			}
		}
	}

	/**
	 * Merges another DynamicDemiLog2 into this one.
	 * Promotes whichever has lower minZeros until they match, then combines normally.
	 */
	public void add(DynamicDemiLog2 log){
		added+=log.added;
		branch1+=log.branch1;
		branch2+=log.branch2;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			// Promote this up to log's minZeros
			while(minZeros<log.minZeros){promote();}

			// If log is lower, compute promoted log values inline (same math as promote())
			final int logPromoSteps=minZeros-log.minZeros;
			final int logShift=logPromoSteps*(1<<mantissabits);

			for(int i=0; i<buckets; i++){
				final int a=maxArray[i];

				final int rawB=log.maxArray[i];
				final int b=(rawB==0 || rawB<=logShift) ? 0 : rawB-logShift;

				final int newMax=Math.max(a, b);
				maxArray[i]=(char)newMax;

				if(newMax==0){
					countArray[i]=0;
				}else if(a==b){
					countArray[i]=(char)Tools.min(countArray[i]+(int)log.countArray[i], Character.MAX_VALUE);
				}else if(b>a){
					countArray[i]=log.countArray[i];
				}
				// else a>b: keep countArray[i]
			}

			filledBuckets=0;
			minZeroCount=0;
			for(char c : maxArray){
				if(c>0){
					filledBuckets++;
					if((c>>mantissabits)==0){minZeroCount++;}
				}
			}
		}
	}

	public int[] compareTo(DynamicDemiLog2 blog){return compare(maxArray, blog.maxArray);}

	@Override
	public void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		//Earliest possible exit, fastest
		if(Long.compareUnsigned(key, eeMask)>0) {return;}
//		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);

		final int bucket=(int)(key&bucketMask);
		final int shift=offset-nlz;
		// Store RELATIVE NLZ = absoluteNlz - minZeros
		final int relNlz=nlz-minZeros;
		final int score=(relNlz<<mantissabits)+(int)((~(key>>>shift))&mask);//FP16 representation
		final int oldValue=maxArray[bucket];//Required memory read

		//Optional early exit reduces writes, countArray access, and branches.
		//Expected to be usually taken, particularly when buckets is large.
		if(score<oldValue) {return;}
//		branch2++;
		lastCardinality=-1;
		final int newValue=Math.max(score, oldValue);
		final int nlzOld=(oldValue>>mantissabits); // relative NLZ of old value

//		assert(newValue>=0 && newValue<=Character.MAX_VALUE) : newValue;
		//Update bucket; required write to cached line only if score>oldValue
		maxArray[bucket]=(char)newValue;
		//Track filled bucket count: increment when bucket transitions from empty to non-empty
		if(oldValue==0 && newValue>0){filledBuckets++;}

		final char count=countArray[bucket];
		countArray[bucket]=(char)(oldValue>score ? count :
			oldValue==score ? Math.max(count, (char)(count+1)) : 1);

		// nlzOld==0 means the old bucket was at the relative minimum tier.
		if(relNlz>nlzOld && nlzOld==0 && --minZeroCount<1){
			// No more buckets at the old minimum; advance floor and decrement all values.
			while(minZeroCount==0 && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
				// Count buckets at new minimum (relNlz==1 before decrement → 0 after)
				// and decrement all non-empty buckets by 1 NLZ.
				minZeroCount=countAndDecrementTerms(1, maxArray, -(1<<mantissabits));
			}
		}
	}

	private int scanFrom(int nlz){
		minZeros=nlz-1;
		minZeroCount=0;
		eeMask=-1L;
		// NOTE: maxArray values are relative to the minZeros being scanned for.
		// After this call, values must be re-relativized if minZeros changed significantly.
		// For now, treat as if re-scanning from scratch (values already converted by add()).
		while(minZeroCount==0 && minZeros<wordlen){
			minZeros++;
			eeMask>>>=1;
			// Count relative-NLZ==0 buckets (those at the new minimum).
			int count=0;
			for(final char c : maxArray){if(c>0 && (c>>mantissabits)==0){count++;}}
			minZeroCount=count;
		}
		return minZeros;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public char[] counts16(){return countArray;}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	/**
	 * Returns cardinality estimates for calibration.
	 * Output order: 0=Mean, 1=HMean, 2=HMeanM, 3=GMean, 4=HLL, 5=LC, 6=Hybrid, 7=MWA, 8=MedianCorr
	 */
	public double[] rawEstimates(){
		double difSum=0;
		double hllSumFilled=0;
		double hllSumFilledM=0;
		double gSum=0;
		int count=0;
		sortBuf.clear();
		final LongList list=sortBuf;

		for(int i=0; i<maxArray.length; i++){
			int max=maxArray[i];
			long val=restore(max);
			if(max>0 && val>0){
				long dif=val;
				difSum+=dif;
				// Absolute NLZ = relNlz + minZeros
				hllSumFilled +=Math.pow(2.0, -(max>>mantissabits)-minZeros);
				hllSumFilledM+=Math.pow(2.0, -(max>>mantissabits)+0.5-(max&mask)/1024.0-minZeros);
				gSum+=Math.log(Tools.max(1, dif));
				count++;
				list.add(dif);
			}
		}

		// All-buckets HLL sum: empty buckets contribute 2^0=1 (absolute NLZ=0 convention).
		double hllSum=0;
		for(int i=0; i<maxArray.length; i++){
			final int v=maxArray[i];
			if(v==0){
				hllSum+=1.0; // empty: standard HLL register=0
			}else{
				hllSum+=Math.pow(2.0, -(v>>mantissabits)-minZeros);
			}
		}
		final double alpha_m=0.7213/(1.0+1.079/buckets);

		final int div=Tools.max(count, 1);
		final double mean=difSum/div;
		double gmean=Math.exp(gSum/div);
		list.sort();
		final long median=Tools.max(1, list.median());
		final double mwa=Tools.max(1.0, list.medianWeightedAverage());

		final int V=buckets-count;

		final double hmeanRaw=2*alpha_m*(double)buckets*(double)buckets/hllSum;
		double hmeanEst=hmeanRaw;
		if(hmeanEst<2.5*buckets && V>0){hmeanEst=(double)buckets*Math.log((double)buckets/V);}

		final double correction=(count+buckets)/(float)(buckets+buckets);
		final double hmeanPure =(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilled);
		final double hmeanPureM=(count==0 ? 0 : 2*alpha_m*(double)count*(double)count/hllSumFilledM);
		final double meanEst   =2*(Long.MAX_VALUE/Tools.max(1.0, mean))*div*correction;
		final double gmeanEst  =2*(Long.MAX_VALUE/gmean)               *div*correction;
		final double mwaEst    =2*(Long.MAX_VALUE/mwa)                 *div*correction;
		final double medianCorr=2*(Long.MAX_VALUE/(double)median)      *div*correction;
		final double lcPure    =buckets*Math.log((double)buckets/Math.max(V, 0.5));

		final int trim=count/256;
		final int trimLow=Math.max(0, trim-V);
		double mean99Sum=0;
		final int mean99N=count-trimLow-trim;
		if(mean99N>0){for(int i=trim; i<count-trimLow; i++){mean99Sum+=list.get(i);}}
		final double mean99=(mean99N>0 ? mean99Sum/mean99N : mean);

		if(filledBuckets==0){return new double[10];}

		final double meanEstCF    =meanEst   *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEAN);
		final double hmeanPureMCF =hmeanPureM*CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEANM);

		final double hybridEst;
		final double hb0=0.5*buckets, hb1=1.5*buckets, hb2=3.0*buckets;
		if(lcPure<=hb0){
			hybridEst=lcPure;
		}else if(lcPure<=hb1){
			final double t=(lcPure-hb0)/(hb1-hb0);
			hybridEst=(1-t)*lcPure+t*meanEstCF;
		}else if(lcPure<=hb2){
			final double t=(lcPure-hb1)/(hb2-hb1);
			hybridEst=(1-t)*meanEstCF+t*hmeanPureMCF;
		}else{
			hybridEst=hmeanPureMCF;
		}

		final double mean99Est=2*(Long.MAX_VALUE/Tools.max(1.0, mean99))*div*correction;
		return new double[]{
			meanEstCF,
			hmeanPure *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HMEAN),
			hmeanPureMCF,
			gmeanEst  *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.GMEAN),
			hmeanEst  *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.HLL),
			lcPure,
			hybridEst,
			mwaEst    *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MWA),
			medianCorr*CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEDCORR),
			mean99Est *CorrectionFactor.getCF(CF_MATRIX, CF_BUCKETS, count, buckets,CorrectionFactor.MEAN99)
		};
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Restores relative compressed score to approximate original hash magnitude.
	 * Adds minZeros to recover absolute NLZ before computing the value.
	 */
	private long restore(int score){
		final int relLeading=(int)(score>>>mantissabits);
		final int leading=relLeading+minZeros; // absolute NLZ
		if(relLeading==0 && minZeros==0){return Long.MAX_VALUE;} // absolute NLZ=0 edge case
		long lowbits=(~score)&mask;
		long mantissa=(1L<<mantissabits)|lowbits;
		int shift=wordlen-leading-mantissabits-1;
		if(shift<0){return Long.MAX_VALUE;} // guards against nlz>53 overflow
		long original=mantissa<<shift;
		return original;
	}

	/**
	 * Counts buckets at target relative NLZ tier (BEFORE decrement), then decrements
	 * all non-empty buckets by incr.  Pass nlz=1, incr=-(1<<mantissabits) when advancing
	 * minZeros: counts the new minimum tier and shifts all values down by 1 NLZ.
	 */
	private static int countAndDecrementTerms(final int nlz, final char[] array, final int incr){
		int sum=0;
		for(int i=0; i<array.length; i++){
			final char c=array[i];
			if(c==0){continue;}
			final int diff=(c>>mantissabits)^nlz;
			final int equalsBit=((~(diff|-diff))>>>31);
			sum+=equalsBit;
			array[i]=(char)(c+incr);
		}
		return sum;
	}

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

	public static float ani(int lower, int equal, int higher){return -1;}
	public static float completeness(int lower, int equal, int higher){return -1;}
	public static float contam(int lower, int equal, int higher){return -1;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final char[] maxArray;
	private final char[] countArray;
	private int minZeros=0;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;
	private final LongList sortBuf=new LongList(buckets);
	// lastCardinality inherited from CardinalityTracker

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, added);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	private static final int mantissabits=10;
	private static final int mask=(1<<mantissabits)-1;
	private static final int offset=wordlen-mantissabits-1;
	
	/** Default resource file for DDL correction factors. */
	public static final String CF_FILE="?cardinalityCorrectionDDL.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the DDL correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	public static double LC_CROSSOVER=0.75;
	public static double LC_SHARPNESS=20.0;
	public static boolean USE_LC=true;

}
