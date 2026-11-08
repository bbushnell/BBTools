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
 * Abstract packed-count Bloom array with optional multi-lane hashing and
 * prefilter support.
 *
 * @author Brian Bushnell
 * @date Jul 5, 2012
 */
public abstract class KCountArray implements Serializable{

	/** Independent secondary lane seed for an exact primitive kmer key; cell selection still mixes every lane. */
	public static long secondaryHash(final long key){
		return Tools.splitMix64(key^0x9E3779B97F4A7C15L);
	}
	
	/** Serialization identifier. */
	private static final long serialVersionUID=1590374813059942002L;

	/** Returns an initialized array with one Bloom lane and no prefilter. */
	public static KCountArray makeNew(long cells_, int cbits_){
		return makeNew(cells_, cbits_, 1);
	}
	
	/** Returns an initialized array with {@code hashes_} Bloom lanes and no prefilter. */
	public static KCountArray makeNew(long cells_, int cbits_, int hashes_){
		return makeNew(cells_, cbits_, hashes_, null, 0);
	}
	
	/** Returns the canonical initialized implementation, optionally with prefilter gating. */
	public static KCountArray makeNew(long cells_, int cbits_, int hashes_, KCountArray prefilter, int prefilterLimit_){
		KCountArray kca=new KCountArray7MTA(cells_, cbits_, hashes_, prefilter, prefilterLimit_);
		kca.initialize();
		return kca;
	}
		
	/** Initializes packed-cell indexing for the legacy fixed-64-partition layout. */
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

	/** Initializes packed-cell indexing for a caller-selected partition count. */
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

	/** Returns the visible count for a primitive key. */
	public abstract int read(long key);
	/** Returns the count for a primitive key with an optional independent secondary hash. */
	public int read(long key, long key2){return read(key);}
	/** Throws unless a subclass supports multiword long-kmer keys. */
	public int read(long keys[]){throw new RuntimeException("Unimplemented.");}
	/** Returns a primitive-key count after optional reverse-complement canonicalization. */
	public final int read(long key, int k, boolean makeCanonical){return read(makeCanonical ? makeCanonical2(key, k) : key);}

	/** Stores an exact value for a primitive key when the implementation supports it. */
	public abstract void write(long key, int value);

	/** Increments a primitive key once. */
	public final void increment(long key){increment(key, 1);}
	/** Increments a dual-hash key by 1. Named to avoid ambiguity with increment(key, amount). */
	public final void incrementDual(long key, long key2){increment(key, key2, 1);}
	/** Decrements a primitive key once if supported. */
	public final void decrement(long key){decrement(key, 1);}

	/** Adds {@code incr} to a primitive key without returning the prior count. */
	public abstract void increment(long key, int incr);
	/** Increments a dual-hash key. Implementations without dual-hash support use the primary key. */
	public void increment(long key, long key2, int incr){increment(key, incr);}
	
	/** Adds {@code incr} to a primitive key and returns the prior count. */
	public abstract int incrementAndReturnUnincremented(long key, int incr);
	/** Returns the unincremented value for a dual-hash key. */
	public int incrementAndReturnUnincremented(long key, long key2, int incr){
		return incrementAndReturnUnincremented(key, incr);
	}
	
//	/** Returns unincremented value */
//	public final int incrementAndReturnUnincremented(Kmer kmer, int incr){
//		return incrementAndReturnUnincremented(kmer.xor(), incr);
//	}
	
	/** Adds {@code incr} to a multiword long-kmer key and returns the prior count. */
	public int incrementAndReturnUnincremented(long[] keys, int incr){
		throw new RuntimeException("Unimplemented.");
	}
	
	/** Throws unless the implementation supports decrementing primitive keys. */
	public void decrement(long key, int incr){
		throw new RuntimeException("This class "+getClass().getName()+" does not support decrement.");
	}
	
	/** Smooths a center key by its left and right single-base neighborhoods. */
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
	
	/** Returns the minimum of a center key and its single-base neighborhoods. */
	public final int readPreciseMin(long key, int k, boolean makeCanonical){
		assert(k<=32);
		int b=read(makeCanonical ? makeCanonical2(key, k) : key);
		if(b<1){return b;}
		int a=readLeft(key, k, makeCanonical);
		if(a<1){return a;}
		int c=readRight(key, k, makeCanonical);
		return Tools.min(a, b, c);
	}
	
	/** Throws unless the implementation can sum the four left neighbors. */
	public int readLeft(long key, int k, boolean makeCanonical){throw new RuntimeException("Unsupported.");}
	/** Throws unless the implementation can sum the four right neighbors. */
	public int readRight(long key, int k, boolean makeCanonical){throw new RuntimeException("Unsupported.");}
	/** Throws unless the implementation can list the four left-neighbor counts. */
	public int[] readAllLeft(final long key, final int k, boolean makeCanonical, int[] rvec){throw new RuntimeException("Unsupported.");}
	/** Throws unless the implementation can list the four right-neighbor counts. */
	public int[] readAllRight(final long key, final int k, boolean makeCanonical, int[] rvec){throw new RuntimeException("Unsupported.");}

	/** Increments all array keys while holding this object's monitor.
	 * @deprecated Prefer {@link #increment(LongList)} to avoid array copying. */
	@Deprecated
	public void increment(long[] keys){
		synchronized(this){
			for(long key : keys){
				increment(key);
			}
		}
	}
	
	/** Increments every key currently present in a LongList under one monitor. */
	public void increment(LongList keys){
		synchronized(this){
			final long[] array=keys.array;
			for(int i=0; i<keys.size; i++){
				increment(array[i]);
			}
		}
	}
	
	/** Returns a histogram over the implementation's stored cell counts. */
	public abstract long[] transformToFrequency();
	/** Unpacks an int matrix into a clamped count histogram. */
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
	
	/** Unpacks an AtomicIntegerArray matrix into a clamped count histogram. */
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
	
	/** Returns cell geometry, memory and occupancy diagnostics. */
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
	
	/** Returns a compact memory/cell/occupancy summary. */
	public final String toShortString(){
		return "mem = "+mem()+"   \tcells = "+toKMG(cells)+"   \tused = "+Tools.format("%.3f%%",usedFraction()*100);
	}
	
	/** Returns a compact memory/cell/occupancy summary with the lane count. */
	public final String toShortString(int hashes){
		return ("hashes = "+hashes+"   \t ")+
				"mem = "+mem()+"   \tcells = "+toKMG(cells)+"   \tused = "+Tools.format("%.3f%%",usedFraction()*100);
	}

	@Override
	public final String toString(){
		return description().toString();
	}
	
	/** Renders stored cell values in implementation-defined order. */
	public abstract CharSequence toContentsString();
	
	/** Returns the fraction of cells reported as nonzero. */
	public abstract double usedFraction();
	
	/** Returns the fraction of cells reported at or above a positive minimum depth. */
	public abstract double usedFraction(int mindepth);
	
	/** Returns the number of cells reported at or above a positive minimum depth. */
	public abstract long cellsUsed(int mindepth);
	
	/** Estimates unique kmers from nonzero-cell occupancy and Bloom lane count. */
	public final double estimateUniqueKmers(int hashes){
		double f=usedFraction();
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	/** Estimates unique kmers from cells meeting a positive minimum depth. */
	public final double estimateUniqueKmers(int hashes, int mindepth){
//		assert(false) : this.getClass().getName();
		double f=usedFraction(mindepth);
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	/** Estimates unique kmers from a caller-supplied occupancy fraction. */
	public final double estimateUniqueKmersFromUsedFraction(int hashes, double usedFraction){
		double f=usedFraction;
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	/** Returns packed cell storage size in binary units. */
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
	
	/** Formats a count with decimal K/M/B suffixes. */
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
	
	/** Allocates a matrix of AtomicIntegerArray rows in parallel. */
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
				try{
					at.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
		return matrix;
	}
	
	/** Worker for parallel AtomicIntegerArray allocation. */
	private static class AllocThread extends Thread{
		
		/** Binds this worker to a matrix and row counter. */
		AllocThread(AtomicIntegerArray[] matrix_, AtomicInteger next_, int wordsPerArray_){
			matrix=matrix_;
			next=next_;
			wordsPerArray=wordsPerArray_;
		}
		
		/** Allocates every matrix row claimed from the shared counter. */
		@Override
		public void run(){
			int x=next.getAndIncrement();
			while(x<matrix.length){
				matrix[x]=new AtomicIntegerArray(wordsPerArray);
				x=next.getAndIncrement();
			}
		}
		
		/** Shared matrix being populated by allocation threads. */
		private final AtomicIntegerArray[] matrix;
		/** Atomic counter for distributing work among threads. */
		private final AtomicInteger next;
		/** Size of each AtomicIntegerArray to allocate. */
		private final int wordsPerArray;
		
	}
	
	
	/** Applies one implementation-specific Bloom lane hash. */
	abstract long hash(long x, int y);
	
	/** Returns the smaller int. */
	public static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger int. */
	public static final int max(int x, int y){return x>y ? x : y;}
	/** Returns the smaller long. */
	public static final long min(long x, long y){return x<y ? x : y;}
	/** Returns the larger long. */
	public static final long max(long x, long y){return x>y ? x : y;}
	
	/** Performs any implementation-specific initialization. */
	public void initialize(){}
	
	/** Performs any implementation-specific finalization. */
	public void shutdown(){}
	
	/** Total number of storage cells in the array. */
	public final long cells;
	/** Number of bits used per storage cell. */
	public final int cellBits;
	/** Originally this was different than valueMask in the case that valueMask was negative, but now they are the same. */
	public final int maxValue;
	
	/** Number of cells that fit in one 32-bit integer word. */
	protected final int cellsPerWord;
	/** Bit shift amount for converting cell index to word index. */
	protected final int indexShift;
	/** Bit mask for extracting cell position within a word. */
	protected final int cellMask;
	/** Bit mask for extracting cell value from packed word. */
	protected final int valueMask;
	
	/** Minimum number of parallel arrays based on thread count. */
	protected static int minArrays=calcMinArrays();
	/** Number of bits needed to index into array selection. */
	protected final int arrayBits;
	/** Total number of parallel arrays used for storage. */
	protected final int numArrays;
	/** Bit mask for selecting which parallel array to use. */
	protected final int arrayMask;
	
	/** Whether to print configuration details during construction. */
	public static boolean verbose=false;
	
	/** Returns the smallest power of two at least as large as the thread count and 2. */
	private static final int calcMinArrays(){
		int x=Tools.max(Shared.threads(), 2);
		while(Integer.bitCount(x)!=1){x++;}
		return x;
	}
	
	/** Returns whether a short key is greater than or equal to its reverse complement. */
	public static final boolean isCanonical(long key, int k){
		assert(k>3 && k<=32);
		long b=AminoAcid.reverseComplementBinaryFast(key, k);
		return key>=b;
	}
	
	/** Returns the reverse complement of a known noncanonical short key. */
	public static final long makeCanonical(final long key, final int k){
		assert(k>3 && k<=32);
//		assert(!isCanonical(key, k));
		final long r=AminoAcid.reverseComplementBinaryFast(key, k);
		assert(r>=key);
//		assert(isCanonical(r, k));
//		assert(AminoAcid.reverseComplementBinaryFast(r, k)==key);
		return r;
	}
	
	
	/** Returns a short key or its reverse complement, whichever is canonical. */
	public static final long makeCanonical2(final long key, final int k){
		assert(k>3 && k<=32);
		if(isCanonical(key, k)){return key;}
		long r=AminoAcid.reverseComplementBinaryFast(key, k);
//		assert(isCanonical(r, k)) : k+"\n"+Long.toBinaryString(key)+"\n"+Long.toBinaryString(r)+"\n"+Long.toBinaryString(AminoAcid.reverseComplementBinaryFast(r, k));
//		assert(AminoAcid.reverseComplementBinaryFast(r, k)==key) : k+"\n"+Long.toBinaryString(key)+"\n"+Long.toBinaryString(r)+"\n"+Long.toBinaryString(AminoAcid.reverseComplementBinaryFast(r, k));
		return r;
	}
	
	/** Returns the associated prefilter when the implementation supports one. */
	public KCountArray prefilter(){
		throw new RuntimeException("TODO: Override");
	}
	
	/** Clears the associated prefilter when the implementation supports one. */
	public void purgeFilter(){
		throw new RuntimeException("TODO: Override");
	}
	
	/** Enables synchronized increments to improve overloaded multi-bit table accuracy. */
	public static boolean LOCKED_INCREMENT=false;
	/** Tracks whether LOCKED_INCREMENT has been explicitly configured. */
	public static boolean SET_LOCKED_INCREMENT=false;
	
}
