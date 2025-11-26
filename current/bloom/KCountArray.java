package bloom;

import java.io.Serializable;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;

import dna.AminoAcid;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 */
public abstract class KCountArray implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1590374813059942002L;

	/**
	 * Creates a new KCountArray with single hash function.
	 * @param cells_ Number of storage cells
	 * @param cbits_ Bits per cell (must be power of 2)
	 * @return New KCountArray instance
	 */
	public static KCountArray makeNew(long cells_, int cbits_){
		return makeNew(cells_, cbits_, 1);
	}
	
	/**
	 * Creates a new KCountArray with specified hash functions.
	 *
	 * @param cells_ Number of storage cells
	 * @param cbits_ Bits per cell (must be power of 2)
	 * @param hashes_ Number of hash functions to use
	 * @return New KCountArray instance
	 */
	public static KCountArray makeNew(long cells_, int cbits_, int hashes_){
		return makeNew(cells_, cbits_, hashes_, null, 0);
	}
	
	//TODO: Get rid of keys_ arg.
	/**
	 * Creates a new KCountArray with prefilter support.
	 * Uses KCountArray7MTA implementation with optional prefiltering for memory efficiency.
	 *
	 * @param cells_ Number of storage cells
	 * @param cbits_ Bits per cell (must be power of 2)
	 * @param hashes_ Number of hash functions to use
	 * @param prefilter Optional prefilter array for reducing memory usage
	 * @param prefilterLimit_ Minimum count threshold for prefilter
	 * @return New initialized KCountArray instance
	 */
	public static KCountArray makeNew(long cells_, int cbits_, int hashes_, KCountArray prefilter, int prefilterLimit_){
		KCountArray kca=new KCountArray7MTA(cells_, cbits_, hashes_, prefilter, prefilterLimit_);
		kca.initialize();
		return kca;
	}
		
	/**
	 * Constructor with default 64 arrays configuration.
	 * Sets up bit manipulation constants and validates parameters.
	 * @param cells_ Number of storage cells (must be power of 2)
	 * @param cbits_ Bits per cell (must be power of 2, ≤32)
	 */
	protected KCountArray(final long cells_, int cbits_){
		assert(cbits_<=32);
		assert(Integer.bitCount(cbits_)==1);
		assert(Long.bitCount(cells_)==1) || this.getClass()==KCountArray7MT.class : this.getClass();
		
		numArrays=64;
		arrayBits=31-Integer.numberOfLeadingZeros(numArrays);
		arrayMask=numArrays-1;
		
		while(cbits_*cells_<32*numArrays){
			assert(false) : cells_+", "+cbits_+", "+numArrays+", "+(cbits_*cells_)+"<"+(32*numArrays);
			cbits_*=2;
		} //Increases bits per cell so that at minimum each array is size 1
		
		assert(cbits_<=32);
		
		cells=cells_;
		cellBits=cbits_;
		valueMask=(cellBits==32 ? Integer.MAX_VALUE : ~((-1)<<cellBits));
		maxValue=min(Integer.MAX_VALUE, ~((-1)<<min(cellBits,31)));
		cellsPerWord=32/cellBits;
		indexShift=Integer.numberOfTrailingZeros(cellsPerWord);
		cellMask=cellsPerWord-1;
		
		if(verbose){
			System.out.println(description());
		}
	}

	/**
	 * Constructor with configurable array count.
	 * Validates parameters and computes bit manipulation constants for efficient indexing.
	 *
	 * @param cells_ Number of storage cells (must be power of 2)
	 * @param cbits_ Bits per cell (must be power of 2, ≤32)
	 * @param arrays_ Number of parallel arrays (must be power of 2)
	 */
	protected KCountArray(final long cells_, int cbits_, int arrays_){
		assert(cbits_<=32);
		assert(Integer.bitCount(cbits_)==1);
		assert(Long.bitCount(cells_)==1) || this.getClass()==KCountArray7MT.class || this.getClass()==KCountArray7MTA.class || this.getClass()==KCountArray8MT.class;

		numArrays=arrays_;
		assert(Integer.bitCount(numArrays)==1) : numArrays+", "+cells_+", "+cbits_;
		arrayBits=31-Integer.numberOfLeadingZeros(numArrays);
		arrayMask=numArrays-1;

		while(cbits_*cells_<32*numArrays){
			assert(false) : cells_+", "+cbits_+", "+numArrays+", "+(cbits_*cells_)+"<"+(32*numArrays);
			cbits_*=2;
		} //Increases bits per cell so that at minimum each array is size 1

		assert(cbits_<=32) : "Why?";

		cells=cells_;
		cellBits=cbits_;
		valueMask=(cellBits==32 ? Integer.MAX_VALUE : ~((-1)<<cellBits));
		maxValue=min(Integer.MAX_VALUE, ~((-1)<<min(cellBits,31)));
		cellsPerWord=32/cellBits;
		indexShift=Integer.numberOfTrailingZeros(cellsPerWord);
		cellMask=cellsPerWord-1;

		if(verbose){
			System.out.println(description());
		}
	}

	/**
	 * Reads the count value for a k-mer key.
	 * @param key K-mer encoded as long integer
	 * @return Count value stored for this k-mer
	 */
	public abstract int read(long key);
	/**
	 * Reads count for long k-mer represented as array.
	 * Default implementation throws exception - override in subclasses supporting long k-mers.
	 *
	 * @param keys K-mer encoded as long array
	 * @return Count value for this k-mer
	 * @throws RuntimeException Always thrown in base implementation
	 */
	public int read(long keys[]){throw new RuntimeException("Unimplemented.");}
	/**
	 * Reads count for k-mer with optional canonical conversion.
	 *
	 * @param key K-mer encoded as long integer
	 * @param k K-mer length in bases
	 * @param makeCanonical Whether to convert to canonical form before lookup
	 * @return Count value for this k-mer
	 */
	public final int read(long key, int k, boolean makeCanonical){return read(makeCanonical ? makeCanonical2(key, k) : key);}

	/**
	 * Sets the count value for a k-mer key.
	 * @param key K-mer encoded as long integer
	 * @param value Count value to store
	 */
	public abstract void write(long key, int value);

	//TODO:  Consider adding a boolean for return old value.
	/** Increments the count for a k-mer by 1.
	 * @param key K-mer encoded as long integer */
	public final void increment(long key){increment(key, 1);}
	/** Decrements the count for a k-mer by 1.
	 * @param key K-mer encoded as long integer */
	public final void decrement(long key){decrement(key, 1);}

	/** Returns nothing for simplicity. */
	public abstract void increment(long key, int incr);
	
	/** Returns unincremented value */
	public abstract int incrementAndReturnUnincremented(long key, int incr);
	
//	/** Returns unincremented value */
//	public final int incrementAndReturnUnincremented(Kmer kmer, int incr){
//		return incrementAndReturnUnincremented(kmer.xor(), incr);
//	}
	
	//For long kmers.
	/**
	 * Atomically increments count for long k-mer and returns previous value.
	 * Default implementation throws exception - override in subclasses supporting long k-mers.
	 *
	 * @param keys K-mer encoded as long array
	 * @param incr Amount to increment
	 * @return Count value before increment
	 * @throws RuntimeException Always thrown in base implementation
	 */
	public int incrementAndReturnUnincremented(long[] keys, int incr){
		throw new RuntimeException("Unimplemented.");
	}
	
	/** Optional method. */
	public void decrement(long key, int incr){
		throw new RuntimeException("This class "+getClass().getName()+" does not support decrement.");
	}
	
	/**
	 * Reads k-mer count with improved precision using adjacent k-mer information.
	 * Averages counts of left and right adjacent k-mers to reduce hash collision effects.
	 *
	 * @param key K-mer encoded as long integer
	 * @param k K-mer length in bases (must be ≤32)
	 * @param makeCanonical Whether to convert to canonical form before lookup
	 * @return More accurate count estimate using adjacent k-mers
	 */
	public final int readPrecise(long key, int k, boolean makeCanonical){
		assert(k<=32);
		int b=read(makeCanonical ? makeCanonical2(key, k) : key);
		if(b<1){return b;}
		int a=readLeft(key, k, makeCanonical);
		if(a>=b){return b;}
		int c=readRight(key, k, makeCanonical);
		if(c>=b){return b;}
		return (int)(((long)a+(long)c)/2);
//		return max(a, c);
//		int mid=Tools.min(a, b, c);
//		System.out.println("a="+a+", b="+b+", c="+c+" -> "+mid);
//		return mid;
	}
	
	/**
	 * Reads k-mer count using minimum of center and adjacent k-mer counts.
	 * Conservative estimate that returns the lowest count among adjacent k-mers.
	 *
	 * @param key K-mer encoded as long integer
	 * @param k K-mer length in bases (must be ≤32)
	 * @param makeCanonical Whether to convert to canonical form before lookup
	 * @return Minimum count among this k-mer and its neighbors
	 */
	public final int readPreciseMin(long key, int k, boolean makeCanonical){
		assert(k<=32);
		int b=read(makeCanonical ? makeCanonical2(key, k) : key);
		if(b<1){return b;}
		int a=readLeft(key, k, makeCanonical);
		if(a<1){return a;}
		int c=readRight(key, k, makeCanonical);
		return Tools.min(a, b, c);
	}
	
	/**
	 * @param key Kmer to evaluate
	 * @return Sum of counts of all 4 possible left-adjacent kmers
	 */
	public int readLeft(long key, int k, boolean makeCanonical){throw new RuntimeException("Unsupported.");}
	/**
	 * @param key Kmer to evaluate
	 * @return Sum of counts of all 4 possible right-adjacent kmers
	 */
	public int readRight(long key, int k, boolean makeCanonical){throw new RuntimeException("Unsupported.");}
	/**
	 * @param key Kmer to evaluate
	 * @return Array of counts of all 4 possible left-adjacent kmers
	 */
	public int[] readAllLeft(final long key, final int k, boolean makeCanonical, int[] rvec){throw new RuntimeException("Unsupported.");}
	/**
	 * @param key Kmer to evaluate
	 * @return Array of counts of all 4 possible right-adjacent kmers
	 */
	public int[] readAllRight(final long key, final int k, boolean makeCanonical, int[] rvec){throw new RuntimeException("Unsupported.");}

	//Appears to never be used?  Should be a LongList now anyway
	/**
	 * Increments counts for multiple k-mers in synchronized block.
	 * Deprecated - use LongList version instead for better performance.
	 * @param keys Array of k-mer keys to increment
	 * @deprecated Use increment(LongList) instead
	 */
	@Deprecated
	public void increment(long[] keys){
		synchronized(this){
			for(long key : keys){
				increment(key);
			}
		}
	}
	
	//TODO: Optionally, add flag to eliminate duplicates or bulk-add them here
	/**
	 * Increments counts for multiple k-mers from LongList.
	 * Thread-safe bulk increment operation.
	 * @param keys LongList containing k-mer keys to increment
	 */
	public void increment(LongList keys){
		synchronized(this){
			final long[] array=keys.array;
			for(int i=0; i<keys.size; i++){
				increment(array[i]);
			}
		}
	}
	
	/**
	 * Generates frequency histogram of count values.
	 * Returns array where index represents count value and value represents frequency.
	 * @return Frequency distribution of counts across all cells
	 */
	public abstract long[] transformToFrequency();
	/**
	 * Generates frequency histogram from integer matrix storage.
	 * Unpacks bit-packed cells and counts frequency of each value.
	 * @param matrix Integer arrays containing bit-packed count data
	 * @return Frequency distribution of counts (index=count, value=frequency)
	 */
	public final long[] transformToFrequency(int[][] matrix){
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
	
	/**
	 * Generates frequency histogram from AtomicIntegerArray storage.
	 * Thread-safe version that unpacks bit-packed cells and counts frequencies.
	 * @param matrix AtomicIntegerArrays containing bit-packed count data
	 * @return Frequency distribution of counts (index=count, value=frequency)
	 */
	public final long[] transformToFrequency(AtomicIntegerArray[] matrix){
		long[] freq=new long[100000];
		int maxFreq=freq.length-1;

		if(cellBits!=32){
			assert(cellBits>0);
			for(AtomicIntegerArray array : matrix){
				for(int i=0; i<array.length(); i++){
					int word=array.get(i);
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
			for(AtomicIntegerArray array : matrix){
				for(int i=0; i<array.length(); i++){
					int word=array.get(i);
					int x2=(int)min(word, maxFreq);
					freq[x2]++;
				}
			}
		}
		return freq;
	}
	
	/**
	 * Creates detailed description of array configuration and statistics.
	 * Includes cell count, bit width, memory usage, and utilization metrics.
	 * @return ByteBuilder containing formatted configuration details
	 */
	public final ByteBuilder description(){
		ByteBuilder sb=new ByteBuilder();
		long words=cells/cellsPerWord;
		int wordsPerArray=(int)(words/numArrays);
		sb.append("cells:   \t"+cells).append('\n');
		sb.append("cellBits:\t"+cellBits).append('\n');
		sb.append("valueMask:\t"+Long.toHexString(valueMask)).append('\n');
		sb.append("maxValue:\t"+maxValue).append('\n');
		sb.append("cellsPerWord:\t"+cellsPerWord).append('\n');
		sb.append("indexShift:\t"+indexShift).append('\n');
		sb.append("words:   \t"+words).append('\n');
		sb.append("wordsPerArray:\t"+wordsPerArray).append('\n');
		sb.append("numArrays:\t"+numArrays).append('\n');
		sb.append("Memory:   \t"+mem()).append('\n');
		sb.append("Usage:    \t"+Tools.format("%.3f%%",usedFraction()*100));
		return sb;
	}
	
	/** Creates compact summary of memory usage and utilization.
	 * @return String with memory, cell count, and usage percentage */
	public final String toShortString(){
		return "mem = "+mem()+"   \tcells = "+toKMG(cells)+"   \tused = "+Tools.format("%.3f%%",usedFraction()*100);
	}
	
	/**
	 * Creates compact summary including hash function count.
	 * @param hashes Number of hash functions used
	 * @return String with hash count, memory, cells, and usage percentage
	 */
	public final String toShortString(int hashes){
		return ("hashes = "+hashes+"   \t ")+
				"mem = "+mem()+"   \tcells = "+toKMG(cells)+"   \tused = "+Tools.format("%.3f%%",usedFraction()*100);
	}

	@Override
	public final String toString(){
		return description().toString();
	}
	
	/**
	 * Returns string representation of array contents.
	 * Implementation varies by subclass storage format.
	 * @return String showing stored count values
	 */
	public abstract CharSequence toContentsString();
	
	/** Calculates fraction of cells with non-zero counts.
	 * @return Proportion of occupied cells (0.0 to 1.0) */
	public abstract double usedFraction();
	
	/**
	 * Calculates fraction of cells with counts at or above threshold.
	 * @param mindepth Minimum count threshold
	 * @return Proportion of cells with count ≥ mindepth (0.0 to 1.0)
	 */
	public abstract double usedFraction(int mindepth);
	
	/**
	 * Counts cells with counts at or above threshold.
	 * @param mindepth Minimum count threshold
	 * @return Number of cells with count ≥ mindepth
	 */
	public abstract long cellsUsed(int mindepth);
	
	/**
	 * Estimates number of unique k-mers based on utilization and hash count.
	 * Uses probabilistic formula accounting for hash collisions.
	 * @param hashes Number of hash functions used
	 * @return Estimated count of unique k-mers stored
	 */
	public final double estimateUniqueKmers(int hashes){
		double f=usedFraction();
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	/**
	 * Estimates unique k-mers with minimum depth threshold.
	 * Uses probabilistic formula on cells meeting depth requirement.
	 *
	 * @param hashes Number of hash functions used
	 * @param mindepth Minimum count threshold
	 * @return Estimated count of unique k-mers with count ≥ mindepth
	 */
	public final double estimateUniqueKmers(int hashes, int mindepth){
//		assert(false) : this.getClass().getName();
		double f=usedFraction(mindepth);
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	/**
	 * Estimates unique k-mers from precomputed utilization fraction.
	 * Avoids recalculating utilization for repeated estimations.
	 *
	 * @param hashes Number of hash functions used
	 * @param usedFraction Fraction of cells with non-zero counts
	 * @return Estimated count of unique k-mers
	 */
	public final double estimateUniqueKmersFromUsedFraction(int hashes, double usedFraction){
		double f=usedFraction;
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	/**
	 * Calculates and formats memory usage with appropriate units.
	 * Returns KB, MB, or GB based on size.
	 * @return Formatted memory usage string with units
	 */
	public final String mem(){
		long mem=(cells*cellBits)/8;
		if(mem<(1<<20)){
			return (Tools.format("%.2f KB", mem*1d/(1<<10)));
		}else if(mem<(1<<30)){
			return (Tools.format("%.2f MB", mem*1d/(1<<20)));
		}else{
			return (Tools.format("%.2f GB", mem*1d/(1<<30)));
		}
	}
	
	/**
	 * Formats large numbers with K/M/B suffixes.
	 * Uses 1000-based units for readability.
	 * @param x Number to format
	 * @return Formatted string with appropriate suffix
	 */
	public static String toKMG(long x){
		double div=1;
		String ext="";
		if(x>10000000000L){
			div=1000000000L;
			ext="B";
		}else if(x>10000000){
			div=1000000;
			ext="M";
		}else if(x>100000){
			div=1000;
			ext="K";
		}
		return Tools.format("%.2f", x/div)+ext;
	}
	
	/**
	 * Allocates matrix of AtomicIntegerArrays using multiple threads.
	 * Parallelizes memory allocation to reduce initialization time for large arrays.
	 *
	 * @param numArrays Number of parallel arrays to create
	 * @param wordsPerArray Size of each array in integers
	 * @return Matrix of initialized AtomicIntegerArrays
	 */
	static final AtomicIntegerArray[] allocMatrix(final int numArrays, final int wordsPerArray){
		final AtomicIntegerArray[] matrix=new AtomicIntegerArray[numArrays];
		final AllocThread[] array=new AllocThread[Tools.min(Tools.max(Shared.threads()/2, 1), numArrays)];
		final AtomicInteger next=new AtomicInteger(0);
		for(int i=0; i<array.length; i++){
			array[i]=new AllocThread(matrix, next, wordsPerArray);
		}
		for(int i=0; i<array.length; i++){array[i].start();}
		for(AllocThread at : array){
			while(at.getState()!=Thread.State.TERMINATED){
				try {
					at.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		return matrix;
	}
	
	/** Worker thread for parallel AtomicIntegerArray allocation.
	 * Reduces memory allocation time by distributing work across threads. */
	private static class AllocThread extends Thread{
		
		/**
		 * Creates allocation thread with shared work parameters.
		 * @param matrix_ Shared matrix to populate
		 * @param next_ Atomic counter for work distribution
		 * @param wordsPerArray_ Size of each array to allocate
		 */
		AllocThread(AtomicIntegerArray[] matrix_, AtomicInteger next_, int wordsPerArray_){
			matrix=matrix_;
			next=next_;
			wordsPerArray=wordsPerArray_;
		}
		
		@Override
		public void run(){
			int x=next.getAndIncrement();
			while(x<matrix.length){
				matrix[x]=new AtomicIntegerArray(wordsPerArray);
				x=next.getAndIncrement();
			}
		}
		
		/** Shared matrix being populated by allocation threads */
		private final AtomicIntegerArray[] matrix;
		/** Atomic counter for distributing work among threads */
		private final AtomicInteger next;
		/** Size of each AtomicIntegerArray to allocate */
		private final int wordsPerArray;
		
	}
	
	
//	long hash(long x, int y){throw new RuntimeException("Not supported.");}
	/**
	 * Computes hash value for k-mer key with hash function index.
	 * Different subclasses implement different hash strategies.
	 *
	 * @param x K-mer encoded as long integer
	 * @param y Hash function index
	 * @return Hash value for array indexing
	 */
	abstract long hash(long x, int y);
	
	/**
	 * Returns minimum of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Smaller of x and y
	 */
	public static final int min(int x, int y){return x<y ? x : y;}
	/**
	 * Returns maximum of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Larger of x and y
	 */
	public static final int max(int x, int y){return x>y ? x : y;}
	/**
	 * Returns minimum of two longs.
	 * @param x First long
	 * @param y Second long
	 * @return Smaller of x and y
	 */
	public static final long min(long x, long y){return x<y ? x : y;}
	/**
	 * Returns maximum of two longs.
	 * @param x First long
	 * @param y Second long
	 * @return Larger of x and y
	 */
	public static final long max(long x, long y){return x>y ? x : y;}
	
	/** Any necessary initialization. */
	public void initialize(){}
	
	/** Any necessary shutdown steps. */
	public void shutdown(){}
	
	/** Total number of storage cells in the array */
	public final long cells;
	/** Number of bits used per storage cell */
	public final int cellBits;
	/** Originally this was different than valueMask in the case that valueMask was negative, but now they are the same. */
	public final int maxValue;
	
	/** Number of cells that fit in one 32-bit integer word */
	protected final int cellsPerWord;
	/** Bit shift amount for converting cell index to word index */
	protected final int indexShift;
	/** Bit mask for extracting cell position within a word */
	protected final int cellMask;
	/** Bit mask for extracting cell value from packed word */
	protected final int valueMask;
	
	/** Minimum number of parallel arrays based on thread count */
	protected static int minArrays=calcMinArrays();
	/** Number of bits needed to index into array selection */
	protected final int arrayBits;
	/** Total number of parallel arrays used for storage */
	protected final int numArrays;
	/** Bit mask for selecting which parallel array to use */
	protected final int arrayMask;
	
	/** Whether to print configuration details during construction */
	public static boolean verbose=false;
	
	/**
	 * Calculates minimum array count based on available threads.
	 * Ensures result is power of 2 and at least 2.
	 * @return Minimum number of arrays (power of 2, ≥2)
	 */
	private static final int calcMinArrays(){
		int x=Tools.max(Shared.threads(), 2);
		while(Integer.bitCount(x)!=1){x++;}
		return x;
	}
	
	/**
	 * Tests if k-mer is in canonical form (lexicographically greater than reverse complement).
	 * @param key K-mer encoded as long integer
	 * @param k K-mer length in bases (4 < k ≤ 32)
	 * @return true if k-mer is canonical, false otherwise
	 */
	public static final boolean isCanonical(long key, int k){
		assert(k>3 && k<=32);
		long b=AminoAcid.reverseComplementBinaryFast(key, k);
		return key>=b;
	}
	
	/** Assumes that the key is not canonical */
	public static final long makeCanonical(final long key, final int k){
		assert(k>3 && k<=32);
//		assert(!isCanonical(key, k));
		final long r=AminoAcid.reverseComplementBinaryFast(key, k);
		assert(r>=key);
//		assert(isCanonical(r, k));
//		assert(AminoAcid.reverseComplementBinaryFast(r, k)==key);
		return r;
	}
	
	
	/**
	 * Converts k-mer to canonical form with canonicality test.
	 * Returns input if already canonical, otherwise returns reverse complement.
	 *
	 * @param key K-mer encoded as long integer
	 * @param k K-mer length in bases (4 < k ≤ 32)
	 * @return Canonical form of the k-mer
	 */
	public static final long makeCanonical2(final long key, final int k){
		assert(k>3 && k<=32);
		if(isCanonical(key, k)){return key;}
		long r=AminoAcid.reverseComplementBinaryFast(key, k);
//		assert(isCanonical(r, k)) : k+"\n"+Long.toBinaryString(key)+"\n"+Long.toBinaryString(r)+"\n"+Long.toBinaryString(AminoAcid.reverseComplementBinaryFast(r, k));
//		assert(AminoAcid.reverseComplementBinaryFast(r, k)==key) : k+"\n"+Long.toBinaryString(key)+"\n"+Long.toBinaryString(r)+"\n"+Long.toBinaryString(AminoAcid.reverseComplementBinaryFast(r, k));
		return r;
	}
	
	/**
	 * Returns associated prefilter array if supported.
	 * Default implementation throws exception - override in subclasses.
	 * @return Prefilter KCountArray instance
	 * @throws RuntimeException If not supported by subclass
	 */
	public KCountArray prefilter(){
		throw new RuntimeException("TODO: Override");
	}
	
	/**
	 * Clears or removes prefilter to free memory.
	 * Default implementation throws exception - override in subclasses.
	 * @throws RuntimeException If not supported by subclass
	 */
	public void purgeFilter(){
		throw new RuntimeException("TODO: Override");
	}
	
	/** Increases accuracy of overloaded multi-bit tables */
	public static boolean LOCKED_INCREMENT=false;
	/** Flag indicating if LOCKED_INCREMENT has been explicitly configured */
	public static boolean SET_LOCKED_INCREMENT=false;
	
}
