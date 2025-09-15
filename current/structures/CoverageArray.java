package structures;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.LineParser;
import shared.LineParser1;
import shared.Shared;
import shared.Tools;


public abstract class CoverageArray implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7175422489330746676L;
	
	
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
		if(atomic && bits32) {return CoverageArray3A.class;}
		else if(!atomic && bits32) {return CoverageArray3.class;}
		else if(atomic && !bits32) {return CoverageArray2A.class;}
		else if(!atomic && !bits32) {return CoverageArray2.class;}
		throw new RuntimeException("No type.");
	}
	
	public CoverageArray(int chrom, int len){
		chromosome=chrom;
		length=len;
	}
	
	/**
	 * @param loc
	 * @param amt
	 */
	public abstract void increment(int loc, int amt);
	
	/**
	 * @param loc
	 */
	public abstract void increment(int loc);

	public final void incrementRange(int min, int max){incrementRange(min, max, 1);}
	public abstract void incrementRange(int min, int max, int amt);
	public abstract void incrementRangeSynchronized(int min, int max, int amt);
	
	public void incrementRanges(IntList ranges, int amt){
		for(int i=0; i<ranges.size; i+=2){
			int a=ranges.get(i), b=ranges.get(i+1);
			incrementRange(a, b-1, 1);
		}
	}
	
	public abstract void set(int loc, int val);
	
	public abstract int get(int loc);
	
	public abstract void resize(int newlen);
	
	
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
		
		return list.toArray(new double[0][]);
		
	}
	
	
	public static final void print(double[][] data){
		
//		data=stats.Smoother.weightedAveragePlank(data, 24);
		assert(false) : "Smoother disabled in this code purely to reduce dependancies.";
		StringBuilder sb=new StringBuilder(data.length*20);
		for(double[] d : data){
			sb.append(Tools.format("%d\t%.2f\n",(int)d[0],d[1]));
		}
		System.out.print(sb);
	}
	
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
	public static HashMap<String, CoverageArray> loadDepth(FileFormat ffdepth, Class<? extends CoverageArray> c) {
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
	
	public final int basesUnderAverageCoverage(final double avg, final int window){
		if(underWindowAverage>-1) {return underWindowAverage;}
		
		//Possibly this should be special-cased to give a real result or return length
		if(length<window){return underWindowAverage=length;}
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0) {return underWindowAverage=length;}
		
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
	
	public final long sum() {
		if(sum>=0) {return sum;}
		sum=0;
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0) {return sum;}
		for(int i=minIndex; i<=maxIndex; i++) {sum+=get(i);}
		return sum;
	}
	
	public final int covered(int minDepth) {
		if(covered>=0) {return covered;}
		calcSumAndCovered(minDepth);
		return covered;
	}
	
	public final long calcSumAndCovered(int minDepth) {
		if(sum>=0 && covered>=0) {return sum;}
		sum=0;
		covered=0;
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0) {return sum;}
		for(int i=minIndex; i<=maxIndex; i++) {
			int d=get(i);
			sum+=d;
			covered+=(d>=minDepth ? 1 : 0);
		}
		return sum;
	}
	
	//Don't use devSum() here!
	public final double standardDeviation(){
		if(stdev>=0) {return stdev;}
		int length=length();
		if(length<2){return 0;}
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0) {return stdev=0;}
		long sum=sum();
		double avg=sum/(double)length;
		double sumdev2=0;
		for(int i=minIndex; i<length; i++){
			int x=get(i);
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return stdev=Math.sqrt(sumdev2/length);
	}
	
	public final double devSum(double globalMean){
		if(devSum>=0) {return devSum;}
		final int length=length();
		int usedDif=maxIndex-minIndex;
		if(usedDif<0) {return devSum=globalMean*globalMean*length;}
		devSum=0;
		if(length<1){return 0;}
		for(int i=minIndex; i<length; i++){
			int x=get(i);
			double dev=globalMean-x;
			devSum+=(dev*dev);
		}
		return devSum;
	}
	
	//Note that empty arrays have a median of zero and do not need sorting
	//Also note that sorting will mess up minIndex and maxIndex
	public final int median(){
		if(median>=0) {return median;}
		final int usedDif=maxIndex-minIndex;
		if(usedDif<0) {return 0;}
		Object o=toArray();
		if(o.getClass()==int[].class) {
			int[] array=(int[])o;
			Shared.sort(array);
			Tools.reverseInPlace(array);
			median=array[array.length/2];
		}else {
			char[] array=(char[])o;
			Arrays.sort(array);
			Tools.reverseInPlace(array);
			median=array[array.length/2];
		}
		//Change the range after sorting
		minIndex=0;
		maxIndex=minIndex+usedDif;
		return median;
	}
	
	public abstract Object toArray();
	
	@Override
	public abstract String toString();
	
	static final long min(long x, long y){return x<y ? x : y;}
	static final long max(long x, long y){return x>y ? x : y;}
	static final int min(int x, int y){return x<y ? x : y;}
	static final int max(int x, int y){return x>y ? x : y;}
	
	public int maxIndex=-1;
	public int minIndex=Integer.MAX_VALUE;
	private final int length;//Note: This is arbitrary if resizing is allowed, as used in old classes...
	public int length(){return length;}
	public abstract int arrayLength();

	public void setUnderWindowAverage(int x) {underWindowAverage=x;}
	public void setMedian(int x) {median=x;}
	public void setCovered(int x) {covered=x;}
	public void setSum(long x) {sum=x;}
	public void setStdev(double x) {stdev=x;}
	public void setDevSum(double x) {devSum=x;}
	
	private int underWindowAverage=-1;
	private int median=-1;
	private int covered=-1;
	private long sum=-1;
	private double stdev=-1;
	private double devSum=-1;
	
	/** Optional */
	public int chromosome;
	
	private static boolean OVERFLOWED=false;
	
}
