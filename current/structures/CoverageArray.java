package structures;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser;
import parse.LineParser1;
import shared.Shared;
import shared.Tools;


/**
 * Base class for per-position coverage arrays.
 * Implementations store nonnegative coverage depths for a single chromosome.
 * Concrete classes differ in storage growth and aggregate bounds: some cached
 * aggregate operations scan the touched minIndex..maxIndex span, while others
 * scan the full logical length.
 */
public abstract class CoverageArray implements Serializable{

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Serialization version identifier. */
	private static final long serialVersionUID=-7175422489330746676L;

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reads a serialized CoverageArray from a .ca file.
	 * @param fname Input path ending in or resolving to a .ca extension.
	 * @return Deserialized coverage array.
	 */
	public static final CoverageArray read(String fname){

		if(!fname.contains(".ca")){
			throw new RuntimeException();
//			ca=new CoverageArray2();
//			ca.load(new TsvCoverageFile(fname));
//			return ca;
		}

		fname=ReadWrite.findFileExtension(fname);
//		System.err.println("Found "+fname);

		return ReadWrite.read(CoverageArray.class, fname, true);

//		if(fname.endsWith(".ca2") || fname.contains(".ca2.")){return ReadWrite.read(CoverageArray2.class, fname);}
//		else if(fname.endsWith(".ca") || fname.contains(".ca.")){return ReadWrite.read(CoverageArray1.class, fname);}
//		else{return ReadWrite.read(CoverageArray.class, fname);}
	}

	public static final Class<? extends CoverageArray> getType(boolean atomic, boolean bits32){
		if(atomic && bits32){return CoverageArray3A.class;}
		else if(!atomic && bits32){return CoverageArray3.class;}
		else if(atomic && !bits32){return CoverageArray2A.class;}
		else if(!atomic && !bits32){return CoverageArray2.class;}
		throw new RuntimeException("No type.");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a coverage array shell for one chromosome.
	 * @param chrom Chromosome identifier.
	 * @param len Logical number of bases represented by this array.
	 */
	public CoverageArray(int chrom, int len){
		chromosome=chrom;
		length=len;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments coverage at specified location by given amount.
	 * @param loc Position to increment
	 * @param amt Amount to add to coverage
	 */
	public abstract void increment(int loc, int amt);

	/** Increments coverage at specified location by 1.
	 * @param loc Position to increment */
	public abstract void increment(int loc);

	/**
	 * Increments an inclusive location range by 1.
	 * @param min First position to increment.
	 * @param max Last position to increment.
	 */
	public final void incrementRange(int min, int max){incrementRange(min, max, 1);}

	/**
	 * Increments an inclusive location range.
	 * @param min First position to increment.
	 * @param max Last position to increment.
	 * @param amt Amount to add at each position.
	 */
	public abstract void incrementRange(int min, int max, int amt);

	/**
	 * Increments an inclusive location range while applying any synchronization
	 * required by the implementation.
	 * @param min First position to increment.
	 * @param max Last position to increment.
	 * @param amt Amount to add at each position.
	 */
	public abstract void incrementRangeSynchronized(int min, int max, int amt);

	/**
	 * Increments each adjacent half-open range encoded as start, stop pairs.
	 * @param ranges Flat list of [start, stop) pairs.
	 * @param amt Amount to add at each position.
	 */
	public void incrementRanges(IntList ranges, int amt){
		for(int i=0; i<ranges.size; i+=2){
			int a=ranges.get(i), b=ranges.get(i+1);
			incrementRange(a, b-1, amt);
		}
	}

	/**
	 * Sets coverage at a location.
	 * @param loc Position to set.
	 * @param val Coverage depth to store.
	 */
	public abstract void set(int loc, int val);

	/**
	 * Returns coverage at a location.
	 * @param loc Position to query.
	 * @return Stored coverage depth.
	 */
	public abstract int get(int loc);

	/**
	 * Resizes the backing storage when supported by the implementation.
	 * @param newlen New logical or backing length.
	 */
	public abstract void resize(int newlen);

	/**
	 * Converts coverage into averaged graph blocks over a bounded interval.
	 * @param blocksize Requested block size, or nonpositive to choose an
	 * Excel-sized automatic block.
	 * @param min Minimum requested position.
	 * @param max Maximum requested position.
	 * @return Rows of rounded mean position and mean coverage.
	 */
	public final double[][] toGraph(int blocksize, int min, int max){

		min=max(min, minIndex);
		max=min(max, maxIndex);
		int length=max-min;

		ArrayList<double[]> list=new ArrayList<double[]>();

		int block;

		if(blocksize<=0){
//			block=((array.length+62999)/63000);//For Excel
//			block=((length+62999)/63000);//For Excel
			block=((length+31499)/31500);//For Excel
		}else{
			block=blocksize;
		}
		block=max(block, 1);

		int current=0;
		double[] sum=new double[2];
		for(int loc=min; loc<=max; loc++){
			if(current==block){
				for(int i=0; i<sum.length; i++){
					sum[i]=sum[i]/current;
				}
				sum[0]=Math.round(sum[0]);
				list.add(sum);
				sum=new double[2];
				current=0;
			}

			sum[0]+=loc;
			sum[1]+=get(loc);

			current++;
		}

		if(current>0){
			for(int i=0; i<sum.length; i++){
				sum[i]=sum[i]/current;
			}
			sum[0]=Math.round(sum[0]);
			list.add(sum);
		}
		return list.toArray(new double[0][]);

	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static final void print(double[][] data){

//		data=stats.Smoother.weightedAveragePlank(data, 24);
		assert(false) : "Smoother disabled in this code purely to reduce dependancies.";
		StringBuilder sb=new StringBuilder(data.length*20);
		for(double[] d : data){
			sb.append(Tools.format("%d\t%.2f\n",(int)d[0],d[1]));
		}
		System.out.print(sb);
	}

	/**
	 * Creates a concrete coverage array of the requested implementation class.
	 * @param num Chromosome identifier.
	 * @param size Logical size.
	 * @param c Implementation class.
	 * @return New coverage array.
	 */
	public static CoverageArray makeArray(int num, int size, Class<? extends CoverageArray> c){
		if(c==CoverageArray2.class){
			return new CoverageArray2(num, size);
		}else if(c==CoverageArray3.class){
			return new CoverageArray3(num, size);
		}else if(c==CoverageArray2A.class){
			return new CoverageArray2A(num, size);
		}else if(c==CoverageArray3A.class){
			return new CoverageArray3A(num, size);
		}
		throw new RuntimeException("Unhandled class: "+c);
	}

	//TODO: Was extremely slow due to string processing. Now, should be fast, but needs verification that LineParser change is correct.
	//This has been modified to allow concise cov files missing fields in lines when they are expected to be the same (or +1).
	/**
	 * Loads a tabular depth file into per-sequence coverage arrays.
	 * Missing sequence, position, or depth fields reuse the previous sequence,
	 * next position, or previous depth respectively.
	 * @param ffdepth Depth input file.
	 * @param c Concrete array implementation to instantiate.
	 * @return Map from sequence name to coverage array.
	 */
	public static HashMap<String, CoverageArray> loadDepth(FileFormat ffdepth, Class<? extends CoverageArray> c){
		ByteFile bf=ByteFile.makeByteFile(ffdepth);
		HashMap<String, CoverageArray> map=new HashMap<String, CoverageArray>();

		String prevName=null;
		int prevPos=-1;
		int prevDepth=0;//TODO: I don't really need the 'prev' variables.

		CoverageArray prevArray=null;
		LineParser lp=new LineParser1('\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line[0]!='#'){
//				String[] split=Tools.tabPattern.split(new String(line));
//				String name=split[0];
//				int pos=Integer.parseInt(split[1]);
//				int depth=Integer.parseInt(split[2]);
				lp.set(line);
				final String name=(Tools.startsWith(line, prevName) || lp.length(0)<1 ? prevName : lp.parseString(0));
				final int pos=(lp.length(1)<1 ? prevPos+1 : lp.parseInt(1));
				final int depth=(lp.length(2)<1 ? prevDepth : lp.parseInt(2));
				CoverageArray current;
				if(name.equals(prevName)){
					current=prevArray;
				}else{
					assert(!map.containsKey(name)) : name; //Could do a lookup but should not be needed
					current=makeArray(map.size()+1, 64, c);
					map.put(name, current);
					prevName=name;
					prevArray=current;
				}
				if(depth>0){current.set(pos, depth);}
				prevPos=pos;
				prevDepth=depth;
			}
		}
		return map;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	public final int basesUnderAverageCoverage(final double avg, final int window){
		if(underWindowAverage>-1){return underWindowAverage;}

		//Possibly this should be special-cased to give a real result or return length
		if(length<window){return underWindowAverage=length;}
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0){return underWindowAverage=length;}

		final long limit=(long)Math.ceil(window*avg);
		long covSum=0;
		int baseCount=0;
		for(int i=0; i<window; i++){
			covSum+=get(i);
		}

		boolean below=false;
		int lastStop=-1, lastStart=0;
		for(int a=0, b=window; b<length; a++, b++){
			if(covSum>=limit){
				if(below){//end range
					baseCount=b-Tools.max(lastStop+1, lastStart);
					lastStop=b-1;
					below=false;
				}
			}else{
				if(!below){//start range
					lastStart=a;
					below=true;
				}
			}
			covSum-=get(a);
			assert(covSum>=0);
			covSum+=get(b);
		}

		if(below){//end range
			baseCount=length()-Tools.max(lastStop, lastStart);
		}

		assert(baseCount>=0) : baseCount+", "+avg+", "+window;
		return underWindowAverage=baseCount;
	}

	public final long sum(){
		if(sum>=0){return sum;}
		sum=0;
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0){return sum;}
		for(int i=minIndex; i<=maxIndex; i++){sum+=get(i);}
		return sum;
	}

	public final int covered(int minDepth){
		if(covered>=0){return covered;}
		calcSumAndCovered(minDepth);
		return covered;
	}

	public final long calcSumAndCovered(int minDepth){
		if(sum>=0 && covered>=0){return sum;}
		sum=0;
		covered=0;
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0){return sum;}
		for(int i=minIndex; i<=maxIndex; i++){
			int d=get(i);
			sum+=d;
			covered+=(d>=minDepth ? 1 : 0);
		}
		return sum;
	}

	//Don't use devSum() here!
	public final double standardDeviation(){
		if(stdev>=0){return stdev;}
		int length=length();
		if(length<2){return 0;}
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0){return stdev=0;}
		long sum=sum();
		double avg=sum/(double)length;
		double sumdev2=0;
		for(int i=0; i<length; i++){
			int x=get(i);
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return stdev=Math.sqrt(sumdev2/length);
	}

	/**
	 * Returns the sum of squared deviations from the supplied mean.
	 * Reuses the cached value only when it was computed for this exact mean;
	 * callers must still reset the cache with unkeyed setDevSum(-1) after
	 * changing coverage because element updates do not invalidate it. Legacy
	 * serialized arrays have no keyed mean and must recompute once.
	 * @param globalMean Mean to compare each coverage value against.
	 * @return Sum of squared deviations for this array.
	 */
	public final double devSum(double globalMean){
		if(devSum>=0 && devSumHasMean && Double.compare(devSumMean, globalMean)==0){
			return devSum;
		}
		final int length=length();
		int usedDif=maxIndex-minIndex;
		if(usedDif<0){return cacheDevSum(globalMean*globalMean*length, globalMean);}
		double sum=0;
		if(length<1){return cacheDevSum(sum, globalMean);}
		for(int i=0; i<length; i++){
			int x=get(i);
			double dev=globalMean-x;
			sum+=(dev*dev);
		}
		return cacheDevSum(sum, globalMean);
	}

	/** Sorts independent storage because toArray() can expose live coordinates.
	 * @return Median coverage depth while preserving coordinate order and touched bounds.
	 */
	public final int median(){
		if(median>=0){return median;}
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0){return 0;}
		Object o=toMedianArray();
		if(o.getClass()==int[].class){
			int[] array=(int[])o;
			Shared.sort(array);
			Tools.reverseInPlace(array);
			median=array[array.length/2];
		}else{
			char[] array=(char[])o;
			Arrays.sort(array);
			Tools.reverseInPlace(array);
			median=array[array.length/2];
		}
		return median;
	}

	/** Returns fresh writable storage for median() sorting.
	 * median() may mutate the returned array; the default implementation
	 * defensively copies public toArray().
	 */
	protected Object toMedianArray(){
		Object o=toArray();
		if(o.getClass()==int[].class){return Arrays.copyOf((int[])o, ((int[])o).length);}
		return Arrays.copyOf((char[])o, ((char[])o).length);
	}

	public abstract Object toArray();

	/** Returns string representation of the coverage array.
	 * @return String describing the coverage data */
	@Override
	public abstract String toString();

	static final long min(long x, long y){return x<y ? x : y;}
	static final long max(long x, long y){return x>y ? x : y;}
	static final int min(int x, int y){return x<y ? x : y;}
	static final int max(int x, int y){return x>y ? x : y;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public int maxIndex=-1;
	public int minIndex=Integer.MAX_VALUE;
	private final int length;//Note: This is arbitrary if resizing is allowed, as used in old classes...
	public int length(){return length;}
	public abstract int arrayLength();

	public void setUnderWindowAverage(int x){underWindowAverage=x;}
	public void setMedian(int x){median=x;}
	public void setCovered(int x){covered=x;}
	public void setSum(long x){sum=x;}
	public void setStdev(double x){stdev=x;}
	/**
	 * Stores an externally computed, unkeyed deviation sum.
	 * The next devSum(mean) call must recompute because this value has no mean.
	 */
	public void setDevSum(double x){
		devSum=x;
		devSumHasMean=false;
	}

	/**
	 * Stores an externally computed deviation sum for the supplied mean.
	 * @param x Sum of squared deviations computed for globalMean.
	 * @param globalMean Mean that was used to compute x.
	 */
	public void setDevSum(double x, double globalMean){
		devSum=x;
		devSumMean=globalMean;
		devSumHasMean=true;
	}
	private double cacheDevSum(double x, double globalMean){
		setDevSum(x, globalMean);
		return devSum;
	}

	//TODO: Mutating concrete classes do not invalidate these cached aggregates,
	//and covered(minDepth) and basesUnderAverageCoverage(avg, window) reuse one
	//cached value across parameter changes. Values can be stale unless a caller
	//resets them.
	private int underWindowAverage=-1;
	private int median=-1;
	private int covered=-1;
	private long sum=-1;
	private double stdev=-1;
	private double devSum=-1;
	private double devSumMean=0;
	// False in old serialized objects, so stale unkeyed devSum values are not reused.
	private boolean devSumHasMean=false;

	/** Optional chromosome identifier for this coverage array */
	public int chromosome;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean OVERFLOWED=false;

}
