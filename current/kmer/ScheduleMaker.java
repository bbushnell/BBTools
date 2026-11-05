package kmer;

import java.util.Arrays;

import shared.Primes;
import shared.Shared;
import shared.Tools;
import structures.IntList;

/**
 * Plans prime-sized primary arrays for sharded k-mer tables without allocating them.
 * Uses a caller-supplied bytes-per-entry estimate and a JVM-heap heuristic; this is
 * not a measurement or guarantee of the complete tables' live memory requirement.
 * Optional prefilter reservations and a caller-enforced bound on concurrent resizes
 * reduce the primary-array budget. Returned schedules are independent mutable arrays.
 * @author Brian Bushnell
 */
public class ScheduleMaker{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a ScheduleMaker with basic memory allocation parameters.
	 * Uses default values for advanced filtering options.
	 *
	 * @param ways_ Number of hash table shards for parallel processing (auto-calculated if {@code <1})
	 * @param bytesPerKmer_ Caller-estimated bytes per primary entry; not measured heap usage
	 * @param prealloc_ Whether to pre-allocate maximum size or allow dynamic resizing
	 * @param memRatio_ Multiplier on the calculated table budget; nonpositive selects 1
	 */
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_){
		this(ways_, bytesPerKmer_, prealloc_, memRatio_, 0, 0, 0, 0);
	}

	/**
	 * Creates a ScheduleMaker with custom initial table size.
	 * Extends basic constructor with control over starting hash table capacity.
	 *
	 * @param ways_ Number of hash table shards for parallel processing
	 * @param bytesPerKmer_ Caller-estimated bytes per primary entry
	 * @param prealloc_ Whether to pre-allocate maximum size or allow dynamic resizing
	 * @param memRatio_ Multiplier on the calculated table budget; nonpositive selects 1
	 * @param initialSize_ Positive lower ladder bound, rounded down to a prime; nonpositive uses the default,
	 * and preallocation ignores this bound
	 */
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_, int initialSize_){
		this(ways_, bytesPerKmer_, prealloc_, memRatio_, initialSize_, 0, 0, 0);
	}

	/**
	 * Full constructor with prefiltering and memory override capabilities.
	 * Calculates memory allocation considering prefiltering passes and manual memory limits.
	 *
	 * @param ways_ Number of hash table shards for parallel processing
	 * @param bytesPerKmer_ Caller-estimated bytes per primary entry
	 * @param prealloc_ Whether to pre-allocate maximum size or allow dynamic resizing
	 * @param memRatio_ Multiplier on the calculated table budget; nonpositive selects 1
	 * @param initialSize_ Positive lower ladder bound, rounded down to a prime; nonpositive uses the default,
	 * and preallocation ignores this bound
	 * @param prepasses_ Number of prefiltering passes; callers use -1 for automatic passes
	 * @param prefilterFraction_ Fraction defining alternating filter budgets with its complement;
	 * the smaller budget is subtracted when planning the main tables
	 * @param filterMemoryOverride_ Requested bytes for each filter allocation; nonpositive selects automatic budgeting
	 */
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_, int initialSize_,
			int prepasses_, double prefilterFraction_, long filterMemoryOverride_){
		this(ways_, bytesPerKmer_, prealloc_, memRatio_, initialSize_, prepasses_,
				prefilterFraction_, filterMemoryOverride_, -1);
	}

	/**
	 * Budgets old resize arrays for a caller-proven bound on concurrent resizers.
	 * A negative bound retains the legacy all-shards assumption. Positive bounds
	 * require the caller to enforce its writer fanout; preallocation needs no old arrays.
	 * Other arguments follow the eight-argument constructor's contract.
	 * @param maxConcurrentResizes_ Maximum simultaneously resizing shards, or -1
	 * @see #ScheduleMaker(int, int, boolean, double, int, int, double, long)
	 */
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_, int initialSize_,
			int prepasses_, double prefilterFraction_, long filterMemoryOverride_, int maxConcurrentResizes_){
		assert(maxConcurrentResizes_==-1 || maxConcurrentResizes_>0) :
			"Resize headroom needs a positive writer bound or -1 for all shards: "+maxConcurrentResizes_;
		bytesPerKmer=bytesPerKmer_;
		prealloc=prealloc_;
		memRatio=(float)(memRatio_<=0 ? 1 : memRatio_);
		prepasses=prepasses_;
		prefilter=(prepasses!=0);
		prefilterFraction=prefilter ? prefilterFraction_ : 0;
		filterMemoryOverride=filterMemoryOverride_;
		assert(prefilter==prefilterFraction>0) : prefilter+", "+prefilterFraction_+", "+prefilterFraction;
//		assert(false && prefilter==prefilterFraction>0) : prefilter+", "+prefilterFraction_+", "+prefilterFraction;
		if(!prefilter){
			filterMemory0=filterMemory1=0;
		}else if(filterMemoryOverride>0){
			filterMemory0=filterMemory1=filterMemoryOverride;
		}else{
			double low=Tools.min(prefilterFraction, 1-prefilterFraction);
			double high=1-low;
			if(prepasses<0 || (prepasses&1)==1){//odd passes
				filterMemory0=(long)(usableMemory*low);
				filterMemory1=(long)(usableMemory*high);
			}else{//even passes
				filterMemory0=(long)(usableMemory*high);
				filterMemory1=(long)(usableMemory*low);
			}
		}
		tableMemory=(long)(usableMemory*.95-Tools.min(filterMemory0, filterMemory1));

		if(ways_<1){
			long maxKmers=(2*tableMemory)/bytesPerKmer;
			long minWays=Tools.min(10000, maxKmers/Integer.MAX_VALUE);
			ways_=(int)Tools.max(31, (int)(Tools.min(96, Shared.threads())*2.5), minWays);
			ways_=(int)Primes.primeAtLeast(ways_);
			assert(ways_>0);
			//		System.err.println("ways="+ways_);
		}
		ways=ways_;

		maxConcurrentResizes=(maxConcurrentResizes_<0 ? ways : Tools.min(ways, maxConcurrentResizes_));
		//Probe tails are fixed allocations, so reserve them after scaling the byte budget.
		final double tailBytes=(ways+(prealloc ? 0L : maxConcurrentResizes))*(double)HashArray.extra*bytesPerKmer;
		final double maxSize0=(tableMemory*0.95*memRatio-tailBytes)/((double)bytesPerKmer*ways);
		assert(maxPrime>1 && maxSize0>2) :
			"\nmaxPrime="+maxPrime+", maxSize0="+maxSize0+", tableMemory="+tableMemory+", usableMemory="+usableMemory+
			", \nprepasses="+prepasses+", filterMemory0="+filterMemory0+", filterMemory1="+filterMemory1+", prefilterFraction="+prefilterFraction+
			", \nmemRatio="+memRatio+", bytesPerKmer="+bytesPerKmer+", ways="+ways+
			", \ninitialSize="+initialSize_+", initialSizeDefault="+initialSizeDefault+", prealloc="+prealloc;

		//The last resize retains old arrays at approximately finalSize/5.
		//Unknown callers reserve this for every shard. A bounded loader reserves
		//it only for the shards its writers can resize simultaneously.
		lastSizeFraction=prealloc ? 1.0 : (maxConcurrentResizes==ways ?
				resizeMult/(1.0+resizeMult) : ways/(ways+maxConcurrentResizes/resizeMult));
		if(!Double.isFinite(maxSize0) || maxSize0*lastSizeFraction<2){
			throw new IllegalArgumentException("No feasible primary-array capacity after filter and probe-tail reservations: "+
				"tableMemory="+tableMemory+", memRatio="+memRatio+", tailBytes="+tailBytes+
				", ways="+ways+", bytesPerKmer="+bytesPerKmer);
		}
		maxSize=Primes.primeAtMost((int)Tools.min(maxPrime, maxSize0*lastSizeFraction));
		initialSize=(prealloc ? maxSize : Primes.primeAtMost(initialSize_>0 ? initialSize_ : initialSizeDefault));

		estimatedKmerCapacity=(long)(maxSize*HashArray.maxLoadFactorFinal*0.97*ways);

//		System.err.println(Arrays.toString(makeSchedule()));

//		System.err.println("ways="+ways+", maxSize="+maxSize+", estimatedKmerCapacity="+estimatedKmerCapacity+", "+Arrays.toString(makeSchedule()));
//
//		assert(false) :
//			"\nmaxPrime="+maxPrime+", maxSize0="+maxSize0+", tableMemory="+tableMemory+", usableMemory="+usableMemory+
//			", \nprepasses="+prepasses+", filterMemory0="+filterMemory0+", filterMemory1="+filterMemory1+", prefilterFraction="+prefilterFraction+
//			", \nmemRatio="+memRatio+", bytesPerKmer="+bytesPerKmer+", ways="+ways+
//			", \ninitialSize="+initialSize_+", initialSizeDefault="+initialSizeDefault+", prealloc="+prealloc+
//			", \nmaxSize="+maxSize+", initialSize="+initialSize+", estimatedKmerCapacity="+estimatedKmerCapacity+
//			", \n"+Arrays.toString(makeSchedule());
//		assert(false) : Arrays.toString(makeSchedule());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds a nondecreasing capacity ladder ending at the budgeted maximum.
	 * Descends by five for the last resize and by three for earlier resizes, then
	 * reverses and rounds upward to primes (except capacity one). The requested
	 * prime-rounded initial size is a lower stopping bound and need not appear in the result.
	 * @return A new capacity array, or one maximum entry when preallocated or too small to resize
	 */
	public int[] makeSchedule(){
		//Prime rounding preserves order; duplicate capacities are harmless to table resizers.
		if(prealloc || maxSize<2L*initialSize){return validateSchedule(new int[] {maxSize});}
		IntList list=new IntList(10);
		list.add(maxSize);
		for(double x=maxSize*invResizeMult; x>=initialSize; x=x*invResizeMult2){
			list.add((int)x);
		}
		if(list.size()>1 && list.lastElement()>=2*initialSize){
			list.add(initialSize);
		}
		list.reverse();
//		if(list.lastElement()*2L<maxPrime){list.add(2*maxSize);}//This ensures that the program will crash rather than garbage-collecting for a long time
		int[] array=list.toArray();
//		if(initialSize>2 && array.length>2 && array[0]>initialSize){array[0]=initialSize;}
		assert(Tools.isSorted(array)) : Arrays.toString(array);
		for(int i=0; i<array.length; i++){
			array[i]=array[i]==1 ? 1 : (int)Tools.min(maxPrime, Primes.primeAtLeast(array[i]));
		}
		return validateSchedule(array);
	}

	/** Checks the actual prime-rounded resize ladder against its reserved budget. */
	private int[] validateSchedule(final int[] array){
		assert(array.length>0 && Tools.isSorted(array) && array[array.length-1]==maxSize) :
			"HashArray.nextScheduleSize requires an ascending ladder ending at its budgeted maximum: "+Arrays.toString(array);
		final long oldSize=(array.length>1 ? array[array.length-2]+(long)HashArray.extra : 0);
		final double peakBytes=((maxSize+(double)HashArray.extra)*ways+oldSize*maxConcurrentResizes)*bytesPerKmer;
		//Probe tails were reserved explicitly; the inner 5% margin covers prime rounding.
		assert(peakBytes<=tableMemory*(double)memRatio) :
			"Prime-rounded table arrays plus bounded old arrays exceed the reserved table budget: peak="+
			peakBytes+", budget="+(tableMemory*(double)memRatio)+", "+this;
		return array;
	}

	/** Explains the scheduled capacity without allocating any kmer table arrays. */
	@Override
	public String toString(){
		assert(ways>0 && maxConcurrentResizes>0) : "Schedule diagnostics require a valid shard/writer budget.";
		return "Kmer table schedule: ways="+ways+", maxConcurrentResizes="+maxConcurrentResizes+
			", resizeMult="+resizeMult+", lastSizeFraction="+lastSizeFraction+
			", bytesPerKmer="+bytesPerKmer+", tableMemory="+tableMemory+
			", maxSize="+maxSize+", estimatedKmerCapacity="+estimatedKmerCapacity;
	}

//	public static int[] makeScheduleStatic(int initialSize, int maxSize, boolean autoResize){
//		if(!autoResize || initialSize>=maxSize){return null;}
//		IntList list=new IntList(10);
//		list.add((int)(maxSize*0.8));
//		for(long x=maxSize; x>=initialSize; x=x/5){
//			list.add((int)x);
//		}
//		if(list.size()>1 && list.lastElement()>=2*initialSize){
//			list.add(initialSize);
//		}
//		list.reverse();
//		int[] array=list.toArray();
//		for(int i=0; i<array.length; i++){
//			array[i]=array[i]==1 ? 1 : (int)Tools.min(maxPrime, Primes.primeAtLeast(array[i]));
//		}
//		return array;
//	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary multiplier for hash table resizing (5.0x growth factor) */
	final double resizeMult=5.0;
	/** Secondary multiplier for intermediate resizing steps (3.0x growth factor) */
	final double resizeMult2=3.0;
	/** Inverse of primary resize multiplier for calculating smaller sizes */
	final double invResizeMult=1/resizeMult;
	/** Inverse of secondary resize multiplier for intermediate size calculations */
	final double invResizeMult2=1/resizeMult2;
	/** Fraction of tail-free capacity left for final arrays after allowing concurrent old arrays. */
	final double lastSizeFraction;
	/** Caller-enforced maximum number of shards retaining old arrays at once. */
	public final int maxConcurrentResizes;

	/** Configured maximum JVM heap, not currently free memory. */
	final long memory=Runtime.getRuntime().maxMemory();
	/** Ratio of initial heap size to maximum heap size from JVM configuration */
	final double xmsRatio=Shared.xmsRatio();

	//TODO: Add term for JDK (Oracle/Open) and version.
	/** Heap heuristic after allowing space for other work and garbage collection. */
	final long usableMemory=(long)Tools.max(((memory-96000000)*
			(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);

	/**
	 * Requested bytes per filter allocation; nonpositive values select automatic budgeting.
	 */
	final long filterMemoryOverride;

	/** Alternating filter allocation budgets; these describe planned bytes, not allocations. */
	final long filterMemory0, filterMemory1;
	/**
	 * Memory available for main k-mer hash tables after reserving prefilter memory
	 */
	public final long tableMemory;
	/** Fraction defining alternating filter budgets with its complement, or zero if disabled. */
	final double prefilterFraction;
	/** Requested pass count, zero for disabled or negative for automatic; this planner does not run passes. */
	final int prepasses;
	/** Whether prefiltering is enabled; zero disables it and negative passes mean automatic. */
	final boolean prefilter;
	/** Caller estimate of bytes per primary entry, excluding separately budgeted probe tails. */
	public final int bytesPerKmer;
	/** Estimated total k-mer capacity across all hash table shards */
	public final long estimatedKmerCapacity;

	/** Number of hash table shards for parallel processing */
	public final int ways;

	/** Lower stopping bound for geometric descent; not necessarily the first returned size. */
	final int initialSize;
	/** Maximum size each hash table can grow to */
	public final int maxSize;

	/**
	 * Whether hash tables should pre-allocate maximum size or resize dynamically
	 */
	final boolean prealloc;
	/**
	 * Memory utilization ratio controlling how much calculated memory to actually use
	 */
	public final float memRatio;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Default starting size for hash tables when no custom size is specified */
	static final int initialSizeDefault=128000;
	/**
	 * Conservative prime ceiling below the int limit, leaving room for probe tails.
	 */
	static final int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE-100-20);

}
