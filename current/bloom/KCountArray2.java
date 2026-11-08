package bloom;

import shared.Tools;

/**
 * Legacy direct-address bit-packed count array for k-mer frequencies.
 * The low two key bits choose one of four rows and the remaining bits address
 * a packed cell in that row. Keys must be in [0,cells). The table is not
 * thread-safe; the cached nonzero-cell count is maintained only by
 * {@link #increment(long, int)}.
 *
 * @author Brian Bushnell
 * @date Jul 5, 2012
 */
public class KCountArray2{
	
	/** Tiny constructor smoke test. */
	public static void main(String[] args){
		KCountArray2 kca=new KCountArray2(1024, 16);
	}
	
	/** Uses the geometry contract of {@link #KCountArray2(long, int, int)} with gap zero. */
	public KCountArray2(long cells_, int bits_){
		this(cells_, bits_, 0);
	}

	/** Constructs a bit-packed count array with a remembered split-kmer gap.
	 * Each row must contain at least one packed word and its word count must fit
	 * a Java array index. Assertions check part of this contract; disabling them
	 * does not make unsupported shapes valid. The legacy assertion-disabled path
	 * can widen the requested cell width for undersized tables.
	 * @param cells_ Positive power-of-two direct-address cell count
	 * @param bits_ Requested bits per cell: 1, 2, 4, 8 or 16, with {@code bits_*cells_>=128}
	 * @param gap_ Split-kmer gap cached for callers that check table geometry */
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
	
	/** Returns the count stored at direct-address {@code key}; valid keys are in [0,{@link #cells}). */
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
	
	/** Stores {@code value} at direct-address {@code key} without updating the
	 * cached nonzero-cell estimate.  Caller must keep {@code key} in
	 * [0,{@link #cells}) and {@code value} in [0,{@link #maxValue}]; out-of-range
	 * values are not masked to this cell before shifting. */
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
	
	/** Adds {@code incr} to direct-address {@code key}, saturates at {@link #maxValue},
	 * and returns the new count.  Caller must keep {@code key} in [0,{@link #cells})
	 * and avoid negative underflow; positive overflow saturates, but negative
	 * values are not floored at zero. */
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
		value=(int)min(value+(long)incr, maxValue);//long add: int value+incr can overflow before saturation.
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		return (int)value;
	}
	
	/** Adds {@code incr} to direct-address {@code key}, saturates at {@link #maxValue},
	 * and returns the old count.  Caller must keep {@code key} in [0,{@link #cells})
	 * and avoid negative underflow; positive overflow saturates, but negative
	 * values are not floored at zero.
	 * This old-value fast path intentionally leaves {@link #cellsUsed} unchanged;
	 * use {@link #cellsUsed(int)} to scan exact occupancy after old-value writes. */
	public int increment2(long key, int incr){
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		final int value=((word>>>cellShift)&valueMask);
		final int value2=(int)min(value+(long)incr, maxValue);//long add: int value+incr can overflow before saturation.
		word=(value2<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		return value;
	}
	
	/** Returns a 100000-bin histogram of all cells, including zero-count cells,
	 * clamping counts above 99999 to the last bucket. */
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
	
	/** Lists cells in row-storage order, not increasing-key order. */
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
	
	/** Returns the cached approximate nonzero fraction tracked by {@link #increment(long, int)} only. */
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	/** Scans for the exact fraction of cells meeting positive threshold {@code mindepth}. */
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	/** Scans for cells meeting positive threshold {@code mindepth}. */
	public long cellsUsed(int mindepth){
		long count=0;
		for(int[] array : matrix){
			if(array!=null){
				for(int word : array){
					long word2=word&0xFFFFFFFFL;
					while(word2!=0){
						int x=(int)(word2&valueMask);
						if(x>=mindepth){count++;}
						word2>>>=cellBits;
					}
				}
			}
		}
		return count;
	}
	
	/** Returns a human-readable estimate of packed counter storage. */
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
	
	/** Returns the smaller integer. */
	public static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger integer. */
	public static final int max(int x, int y){return x>y ? x : y;}
	/** Returns the smaller long. */
	public static final long min(long x, long y){return x<y ? x : y;}
	/** Returns the larger long. */
	public static final long max(long x, long y){return x>y ? x : y;}
	
	/** Cached occupancy updated by increment's zero/nonzero transitions, not by write or increment2. */
	private long cellsUsed;
	
	/** Direct-address cell count. */
	public final long cells;
	/** Effective bits per packed cell. */
	public final int cellBits;
	/** Largest supported value in one cell. */
	public final int maxValue;
	/** Split-kmer gap cached for callers that check they are using the matching table. */
	public final int gap;
	
	/** Packed cells per 32-bit word. */
	private final int cellsPerWord;
	/** Bit shift from per-row cell number to word index. */
	private final int indexShift;
	/** Mask for one packed cell. */
	private final int valueMask;
	/** Four packed counter rows selected by low key bits. */
	private final int[][] matrix;
	
	/** Number of low key bits used to choose a row. */
	private static final int arrayBits=2;
	/** Number of independent rows. */
	private static final int numArrays=1<<arrayBits;
	/** Mask for the low row-selection bits. */
	private static final int arrayMask=numArrays-1;
	
	/** Print constructor geometry diagnostics. */
	public static boolean verbose=false;
	
	
}
