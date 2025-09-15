package kmer;

import java.util.Arrays;

import shared.Primes;
import shared.Shared;
import shared.Tools;
import structures.IntList;

public class ScheduleMaker {
	
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_){
		this(ways_, bytesPerKmer_, prealloc_, memRatio_, 0, 0, 0, 0);
	}
	
	public ScheduleMaker(int ways_, int bytesPerKmer_, boolean prealloc_, double memRatio_, int initialSize_){
		this(ways_, bytesPerKmer_, prealloc_, memRatio_, initialSize_, 0, 0, 0);
	}
	
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

	final double resizeMult=5.0;
	final double resizeMult2=3.0;
	final double invResizeMult=1/resizeMult;
	final double invResizeMult2=1/resizeMult2;
	final double lastSizeFraction;
	
	final long memory=Runtime.getRuntime().maxMemory();
	final double xmsRatio=Shared.xmsRatio();
	
	//TODO: Add term for JDK (Oracle/Open) and version.
	final long usableMemory=(long)Tools.max(((memory-96000000)*
			(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
	
	final long filterMemoryOverride=0;
	
	final long filterMemory0, filterMemory1;
	final long tableMemory;
	final double prefilterFraction;
	final int prepasses;
	final boolean prefilter;
	final int bytesPerKmer;
	public final long estimatedKmerCapacity;
	
	public final int ways;

	final int initialSize;
	final int maxSize;
	
	final boolean prealloc;
	final float memRatio;

	static final int initialSizeDefault=128000;
	static final int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE-100-20);
	
}
