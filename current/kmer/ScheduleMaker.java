package kmer;

import java.util.Arrays;

import shared.Primes;
import shared.Shared;
import shared.Tools;
import structures.IntList;

/**
 * Calculates optimal memory allocation and resizing strategy for k-mer hash tables.
 * Generates prime-sized hash table schedule with intelligent memory management
 * based on available system memory, prefiltering requirements, and parallel processing needs.
 * @author Brian Bushnell
 */
public class ScheduleMaker {
	
	/**
	 * Creates a ScheduleMaker with basic memory allocation parameters.
	 * Uses default values for advanced filtering options.
	 *
	 * @param ways_ Number of hash table shards for parallel processing (auto-calculated if <1)
	 * @param bytesPerKmer_ Memory footprint per k-mer entry in bytes
	 * @param prealloc_ Whether to pre-allocate maximum size or allow dynamic resizing
	 * @param memRatio_ Memory utilization ratio (1.0 = use all available calculated memory)
	 */
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_){
		this(ways_, bytesPerKmer_, prealloc_, memRatio_, 0, 0, 0, 0);
	}
	
	/**
	 * Creates a ScheduleMaker with custom initial table size.
	 * Extends basic constructor with control over starting hash table capacity.
	 *
	 * @param ways_ Number of hash table shards for parallel processing
	 * @param bytesPerKmer_ Memory footprint per k-mer entry in bytes
	 * @param prealloc_ Whether to pre-allocate maximum size or allow dynamic resizing
	 * @param memRatio_ Memory utilization ratio
	 * @param initialSize_ Starting size for hash tables before resizing
	 */
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_, int initialSize_){
		this(ways_, bytesPerKmer_, prealloc_, memRatio_, initialSize_, 0, 0, 0);
	}
	
	/**
	 * Full constructor with prefiltering and memory override capabilities.
	 * Calculates memory allocation considering prefiltering passes and manual memory limits.
	 *
	 * @param ways_ Number of hash table shards for parallel processing
	 * @param bytesPerKmer_ Memory footprint per k-mer entry in bytes
	 * @param prealloc_ Whether to pre-allocate maximum size or allow dynamic resizing
	 * @param memRatio_ Memory utilization ratio
	 * @param initialSize_ Starting size for hash tables before resizing
	 * @param prepasses_ Number of prefiltering passes (enables prefiltering if >0)
	 * @param prefilterFraction_ Fraction of memory allocated to prefiltering vs main tables
	 * @param filterMemoryOverride_ Manual override for filter memory allocation (0 = auto-calculate)
	 */
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_, int initialSize_, 
			int prepasses_, double prefilterFraction_, long filterMemoryOverride_){
		bytesPerKmer=bytesPerKmer_;
		prealloc=prealloc_;
		memRatio=(float)(memRatio_<=0 ? 1 : memRatio_);
		prepasses=prepasses_;
		prefilter=(prepasses>0);
		prefilterFraction=prefilter ? prefilterFraction_ : 0;
		assert(prefilter==prefilterFraction>0) : prefilter+", "+prefilterFraction_+", "+prefilterFraction;
//		assert(false && prefilter==prefilterFraction>0) : prefilter+", "+prefilterFraction_+", "+prefilterFraction;
		if(prepasses<1){
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

		final double maxSize0=(tableMemory*0.95*memRatio)/(bytesPerKmer*ways);
		assert(maxPrime>1 && maxSize0>2) : 
			"\nmaxPrime="+maxPrime+", maxSize0="+maxSize0+", tableMemory="+tableMemory+", usableMemory="+usableMemory+
			", \nprepasses="+prepasses+", filterMemory0="+filterMemory0+", filterMemory1="+filterMemory1+", prefilterFraction="+prefilterFraction+
			", \nmemRatio="+memRatio+", bytesPerKmer="+bytesPerKmer+", ways="+ways+
			", \ninitialSize="+initialSize_+", initialSizeDefault="+initialSizeDefault+", prealloc="+prealloc;
		
		lastSizeFraction=prealloc ? 1.0 : resizeMult/(1.0+resizeMult);
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
	
	/**
	 * Generates an array of prime-sized hash table capacities for progressive resizing.
	 * Creates a schedule that starts at initialSize and grows to maxSize in geometric steps.
	 * All sizes are converted to prime numbers for optimal hash distribution.
	 * @return Array of prime table sizes in ascending order, or single maxSize if prealloc=true
	 */
	public int[] makeSchedule(){
		if(prealloc || maxSize<2L*initialSize){return new int[] {maxSize};}
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
		return array;
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

	/** Primary multiplier for hash table resizing (5.0x growth factor) */
	final double resizeMult=5.0;
	/** Secondary multiplier for intermediate resizing steps (3.0x growth factor) */
	final double resizeMult2=3.0;
	/** Inverse of primary resize multiplier for calculating smaller sizes */
	final double invResizeMult=1/resizeMult;
	/** Inverse of secondary resize multiplier for intermediate size calculations */
	final double invResizeMult2=1/resizeMult2;
	/** Fraction of maximum size to use for final table capacity calculation */
	final double lastSizeFraction;
	
	/** Total JVM heap memory available at runtime */
	final long memory=Runtime.getRuntime().maxMemory();
	/** Ratio of initial heap size to maximum heap size from JVM configuration */
	final double xmsRatio=Shared.xmsRatio();
	
	//TODO: Add term for JDK (Oracle/Open) and version.
	/** Calculated usable memory after reserving space for JVM overhead and GC */
	final long usableMemory=(long)Tools.max(((memory-96000000)*
			(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
	
	/**
	 * Manual override for prefilter memory allocation (0 = use automatic calculation)
	 */
	final long filterMemoryOverride=0;
	
	/** Memory allocated for second prefilter pass */
	/** Memory allocated for first prefilter pass */
	final long filterMemory0, filterMemory1;
	/**
	 * Memory available for main k-mer hash tables after reserving prefilter memory
	 */
	final long tableMemory;
	/** Fraction of total memory allocated to prefiltering operations */
	final double prefilterFraction;
	/** Number of prefiltering passes to perform before main table construction */
	final int prepasses;
	/** Whether prefiltering is enabled (true when prepasses > 0) */
	final boolean prefilter;
	/** Memory footprint per k-mer entry in the hash tables */
	final int bytesPerKmer;
	/** Estimated total k-mer capacity across all hash table shards */
	public final long estimatedKmerCapacity;
	
	/** Number of hash table shards for parallel processing */
	public final int ways;

	/** Starting size for hash tables before any resizing operations */
	final int initialSize;
	/** Maximum size each hash table can grow to */
	final int maxSize;
	
	/**
	 * Whether hash tables should pre-allocate maximum size or resize dynamically
	 */
	final boolean prealloc;
	/**
	 * Memory utilization ratio controlling how much calculated memory to actually use
	 */
	final float memRatio;

	/** Default starting size for hash tables when no custom size is specified */
	static final int initialSizeDefault=128000;
	/**
	 * Largest prime number that fits in a 32-bit integer, used as upper bound for table sizes
	 */
	static final int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE-100-20);
	
}
