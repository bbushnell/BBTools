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
	
	/**
	 * Creates a count array with specified cell count and bit width.
	 * @param cells_ Number of cells in the array (must be power of 2)
	 * @param bits_ Bits per cell for count storage (must be power of 2, max 32)
	 */
	public KCountArray2(long cells_, int bits_){
		this(cells_, bits_, 0);
	}
		
	/**
	 * Creates a count array with specified parameters including gap setting.
	 * Distributes cells across multiple arrays for load balancing and initializes
	 * bit manipulation constants for efficient packing/unpacking.
	 *
	 * @param cells_ Number of cells in the array (must be power of 2)
	 * @param bits_ Bits per cell for count storage (must be power of 2, max 32)
	 * @param gap_ Gap parameter for convenience tracking on gapped tables
	 */
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
	
	/**
	 * Reads the count value stored at the specified key.
	 * Uses array distribution and bit manipulation to extract the packed value.
	 * @param key Hash key identifying the cell location
	 * @return Count value stored at the key location
	 */
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
	
	/**
	 * Writes a count value to the specified key location.
	 * Uses bit masking to update only the target cell bits within the packed word.
	 * @param key Hash key identifying the cell location
	 * @param value Count value to store (will be clamped to maxValue)
	 */
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
	
	/**
	 * Increments the count at the specified key by the given amount.
	 * Tracks cell usage statistics and clamps result to maxValue.
	 *
	 * @param key Hash key identifying the cell location
	 * @param incr Amount to add to current count (can be negative)
	 * @return New count value after increment
	 */
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
	
	/**
	 * Generates frequency histogram of count values across all cells.
	 * Iterates through all packed cells and counts occurrences of each count value.
	 * @return Array where index is count value and value is frequency of that count
	 */
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
	
	/** Calculates fraction of cells containing non-zero counts.
	 * @return Proportion of cells that have been used (0.0 to 1.0) */
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	/**
	 * Calculates fraction of cells with counts at or above minimum depth.
	 * @param mindepth Minimum count threshold
	 * @return Proportion of cells meeting the depth requirement
	 */
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	/**
	 * Counts number of cells with counts at or above minimum depth.
	 * Unpacks all cells and counts those meeting the threshold.
	 * @param mindepth Minimum count threshold
	 * @return Number of cells with counts >= mindepth
	 */
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
	
	/** Returns formatted string describing memory usage of this array.
	 * @return Memory usage string with appropriate units (KB/MB/GB) */
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
	
	/** Returns the smaller of two integers */
	public static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger of two integers */
	public static final int max(int x, int y){return x>y ? x : y;}
	/** Returns the smaller of two longs */
	public static final long min(long x, long y){return x<y ? x : y;}
	/** Returns the larger of two longs */
	public static final long max(long x, long y){return x>y ? x : y;}
	
	/** Number of cells currently containing non-zero counts */
	private long cellsUsed;
	
	/** Total number of cells in the array */
	public final long cells;
	/** Number of bits allocated per cell for count storage */
	public final int cellBits;
	/** Maximum count value that can be stored per cell */
	public final int maxValue;
	/** Gap parameter for convenience tracking on gapped k-mer tables */
	public final int gap; //Set this for convenience on gapped tables to make sure you're using the right table.
	
	/** Number of cells that fit in a 32-bit word given current cellBits */
	private final int cellsPerWord;
	/** Bit shift amount for extracting word index from key */
	private final int indexShift;
	/** Bit mask for extracting cell values from packed words */
	private final int valueMask;
	/**
	 * 2D array matrix storing bit-packed count values across multiple sub-arrays
	 */
	private final int[][] matrix;
	
	/** Number of bits used for array selection in key distribution */
	private static final int arrayBits=2;
	/** Number of parallel arrays for load distribution */
	private static final int numArrays=1<<arrayBits;
	/** Bit mask for extracting array selection bits from key */
	private static final int arrayMask=numArrays-1;
	
	/** Whether to print detailed debugging information during construction */
	public static boolean verbose=false;
	
	
}
