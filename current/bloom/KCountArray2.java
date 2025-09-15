package bloom;

import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 */
public class KCountArray2 {
	
	public static void main(String[] args){
		KCountArray2 kca=new KCountArray2(1024, 16);
	}
	
	public KCountArray2(long cells_, int bits_){
		this(cells_, bits_, 0);
	}
		
	public KCountArray2(long cells_, int bits_, int gap_){
		gap=gap_;
		assert(bits_<=32);
		assert(Integer.bitCount(bits_)==1);
		assert(Long.bitCount(cells_)==1);
		
		while(bits_*cells_<32*numArrays){
			assert(false);
			bits_*=2;
		} //Increases bits per cell so that at minimum each array is size 1
		
		assert(bits_!=32);
		
		cells=cells_;
		cellBits=bits_;
		valueMask=~((-1)<<cellBits);
		maxValue=min(Integer.MAX_VALUE, ~((-1)<<min(cellBits,31)));
		cellsPerWord=32/cellBits;
		indexShift=Integer.numberOfTrailingZeros(cellsPerWord);
		long words=cells/cellsPerWord;
		int wordsPerArray=(int)(words/numArrays);
		matrix=new int[numArrays][wordsPerArray];
		
		if(verbose){
			System.out.println("cells:   \t"+cells);
			System.out.println("cellBits:\t"+cellBits);
			System.out.println("valueMask:\t"+Long.toHexString(valueMask));
			System.out.println("maxValue:\t"+maxValue);
			System.out.println("cellsPerWord:\t"+cellsPerWord);
			System.out.println("indexShift:\t"+indexShift);
			System.out.println("words:   \t"+words);
			System.out.println("wordsPerArray:\t"+wordsPerArray);
			System.out.println("numArrays:\t"+numArrays);


			long mem=words*4;
			if(mem<(1<<30)){
				System.out.println("memory:   \t"+Tools.format("%.2f MB", mem*1d/(1<<20)));
			}else{
				System.out.println("memory:   \t"+Tools.format("%.2f GB", mem*1d/(1<<30)));
			}
		}
	}
	
	public int read(long key){
//		System.out.println("key="+key);
		int arrayNum=(int)(key&arrayMask);
//		System.out.println("array="+arrayNum);
		key>>>=arrayBits;
//		System.out.println("key2="+key);
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
//		System.out.println("index="+index);
		int word=array[index];
//		System.out.println("word="+Integer.toHexString(word));
		int cellShift=(int)(cellBits*key);
//		System.out.println("cellShift="+cellShift);
		return (int)((word>>>cellShift)&valueMask);
	}
	
	public void write(long key, int value){
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
	}
	
	public int increment(long key, int incr){
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		int value=((word>>>cellShift)&valueMask);
		if(value==0 && incr>0){cellsUsed++;}
		else if(incr<0 && value+incr==0){cellsUsed--;}
		value=min(value+incr, maxValue);
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		return (int)value;
	}
	
	/** Returns unincremented value */
	public int increment2(long key, int incr){
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		final int value=((word>>>cellShift)&valueMask);
		final int value2=min(value+incr, maxValue);
		word=(value2<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		return value;
	}
	
	public long[] transformToFrequency(){
		long[] freq=new long[100000];
		int maxFreq=freq.length-1;

		if(cellBits!=32){
			assert(cellBits>0);
			for(int[] array : matrix){
				for(int i=0; i<array.length; i++){
					int word=array[i];
					int j=cellsPerWord;
					//				System.out.println("initial: word = "+word+", j = "+Integer.toHexString(j)+", cellbits="+cellBits);
					for(; word!=0; j--){
						int x=word&valueMask;
						int x2=(int)min(x, maxFreq);
						freq[x2]++;
						word=(word>>>cellBits);
						//					System.out.println("word = "+word+", j = "+Integer.toHexString(j)+", cellbits="+cellBits);
					}
					freq[0]+=j;
				}
			}
		}else{
			for(int[] array : matrix){
				for(int i=0; i<array.length; i++){
					int word=array[i];
					int x2=(int)min(word, maxFreq);
					freq[x2]++;
				}
			}
		}
		return freq;
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append("[");
		String comma="";
		for(int[] array : matrix){
			for(int i=0; i<array.length; i++){
				int word=array[i];
				for(int j=0; j<cellsPerWord; j++){
					int x=word&valueMask;
					sb.append(comma);
					sb.append(x);
					word>>>=cellBits;
					comma=", ";
				}
			}
		}
		sb.append("]");
		return sb.toString();
	}
	
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	public long cellsUsed(int mindepth){
		long count=0;
		for(int[] array : matrix){
			if(array!=null){
				for(int word : array){
					while(word>0){
						int x=word&valueMask;
						if(x>=mindepth){count++;}
						word>>>=cellBits;
					}
				}
			}
		}
		return count;
	}
	
	public String mem(){
		long mem=(cells*cellBits)/8;
		if(mem<(1<<20)){
			return (Tools.format("%.2f KB", mem*1d/(1<<10)));
		}else if(mem<(1<<30)){
			return (Tools.format("%.2f MB", mem*1d/(1<<20)));
		}else{
			return (Tools.format("%.2f GB", mem*1d/(1<<30)));
		}
	}
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	public static final long min(long x, long y){return x<y ? x : y;}
	public static final long max(long x, long y){return x>y ? x : y;}
	
	private long cellsUsed;
	
	public final long cells;
	public final int cellBits;
	public final int maxValue;
	public final int gap; //Set this for convenience on gapped tables to make sure you're using the right table.
	
	private final int cellsPerWord;
	private final int indexShift;
	private final int valueMask;
	private final int[][] matrix;
	
	private static final int arrayBits=2;
	private static final int numArrays=1<<arrayBits;
	private static final int arrayMask=numArrays-1;
	
	public static boolean verbose=false;
	
	
}
