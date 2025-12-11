package bloom;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;


/**
 * Multi-threaded k-mer counting array with buffered writes and prime-sized arrays.
 * Uses dedicated WriteThread instances for each array to enable parallel increments
 * while maintaining thread safety. Employs prime numbers for array lengths to reduce
 * hash collisions and includes prefilter support for two-stage counting.
 *
 * @author Brian Bushnell
 * @date Aug 17, 2012
 */
public class KCountArray8MT extends KCountArray {
	
	/** Serial version ID for serialization compatibility */
	private static final long serialVersionUID = -3146298383509476887L;

	public static void main(String[] args){
		long cells=Long.parseLong(args[0]);
		int bits=Integer.parseInt(args[1]);
		int hashes=Integer.parseInt(args[3]);
		
		verbose=false;
		
		KCountArray8MT kca=new KCountArray8MT(cells, bits, hashes, null);
		
		System.out.println(kca.read(0));
		kca.increment(0);
		System.out.println(kca.read(0));
		kca.increment(0);
		System.out.println(kca.read(0));
		System.out.println();
		
		System.out.println(kca.read(1));
		kca.increment(1);
		System.out.println(kca.read(1));
		kca.increment(1);
		System.out.println(kca.read(1));
		System.out.println();
		
		System.out.println(kca.read(100));
		kca.increment(100);
		System.out.println(kca.read(100));
		kca.increment(100);
		System.out.println(kca.read(100));
		kca.increment(100);
		System.out.println(kca.read(100));
		System.out.println();
		

		System.out.println(kca.read(150));
		kca.increment(150);
		System.out.println(kca.read(150));
		System.out.println();
		
	}
		
	public KCountArray8MT(long cells_, int bits_, int hashes_, KCountArray prefilter_){
		super(getPrimeCells(cells_, bits_), bits_, getDesiredArrays(cells_, bits_));
//		verbose=false;
//		assert(false);
		cellsPerArray=cells/numArrays;
		wordsPerArray=(int)((cellsPerArray%cellsPerWord)==0 ? (cellsPerArray/cellsPerWord) : (cellsPerArray/cellsPerWord+1));
		cellMod=cellsPerArray;
		hashes=hashes_;
//		System.out.println("cells="+cells+", words="+words+", wordsPerArray="+wordsPerArray+", numArrays="+numArrays+", hashes="+hashes);
//		assert(false);
		matrix=new int[numArrays][];
		prefilter=prefilter_;
		assert(prefilter!=null);
		assert(hashes>0 && hashes<=hashMasks.length);
	}
	
	private static int getDesiredArrays(long desiredCells, int bits){
		
		long words=Tools.max((desiredCells*bits+31)/32, minArrays);
		int arrays=minArrays;
		while(words/arrays>=Integer.MAX_VALUE){
			arrays*=2;
		}
		return arrays;
	}
	
	private static long getPrimeCells(long desiredCells, int bits){
		
		int arrays=getDesiredArrays(desiredCells, bits);
		
		long x=(desiredCells+arrays-1)/arrays;
		long x2=Primes.primeAtMost(x);
		return x2*arrays;
	}
	
	/**
	 * Reads the count for a k-mer using multiple hash functions.
	 * Returns the minimum count across all hash positions after checking prefilter.
	 * @param rawKey The k-mer key to read
	 * @return Minimum count across all hash functions
	 */
	@Override
	public int read(final long rawKey){
		assert(finished);
		if(verbose){System.err.println("Reading raw key "+rawKey);}
		int pre=0;
		if(prefilter!=null){
			pre=prefilter.read(rawKey);
			if(pre<prefilter.maxValue){return pre;}
		}
		long key2=hash(rawKey, 0);
		int min=readHashed(key2);
		for(int i=1; i<hashes && min>0; i++){
			if(verbose){System.err.println("Reading. i="+i+", key2="+key2);}
			key2=Long.rotateRight(key2, hashBits);
			key2=hash(key2, i);
			if(verbose){System.err.println("Rot/hash. i="+i+", key2="+key2);}
			min=min(min, readHashed(key2));
		}
		return min;
	}
	
	private int readHashed(long key){
		if(verbose){System.err.print("Reading hashed key "+key);}
//		System.out.println("key="+key);
		int arrayNum=(int)(key&arrayMask);
		key=(key>>>arrayBits)%(cellMod);
//		key=(key>>>(arrayBits+1))%(cellMod);
//		System.out.println("array="+arrayNum);
//		System.out.println("key2="+key);
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
//		assert(false) : indexShift;
//		System.out.println("index="+index);
		int word=array[index];
//		System.out.println("word="+Integer.toHexString(word));
		assert(word>>>(cellBits*key) == word>>>(cellBits*(key&cellMask)));
//		int cellShift=(int)(cellBits*(key&cellMask));
		int cellShift=(int)(cellBits*key);
		if(verbose){System.err.println(", array="+arrayNum+", index="+index+", cellShift="+(cellShift%32)+", value="+((int)((word>>>cellShift)&valueMask)));}
//		System.out.println("cellShift="+cellShift);
		return (int)((word>>>cellShift)&valueMask);
	}
	
	/**
	 * Direct write operation is not supported for this thread-safe implementation.
	 * Use increment operations instead for proper thread safety.
	 *
	 * @param key The k-mer key
	 * @param value The value to write
	 * @throws RuntimeException Always thrown as operation is not allowed
	 */
	@Override
	public void write(final long key, int value){
		throw new RuntimeException("Not allowed for this class.");
	}
	
	/**
	 * Increments counts for multiple k-mers with optimized batching.
	 * Pre-hashes keys outside critical section and filters through prefilter
	 * before passing to buffered increment system.
	 * @param keys Array of k-mer keys to increment
	 */
	@Override
	/** This should increase speed by doing the first hash outside the critical section, but it does not seem to help. */
	public void increment(long[] keys){
		if(prefilter==null){
			for(int i=0; i<keys.length; i++){
				keys[i]=hash(keys[i], 0);
			}
			synchronized(buffers){
				for(long key : keys){
					incrementPartiallyHashed(key);
				}
			}
		}else{
			int j=0;
			for(int i=0; i<keys.length; i++){
				long key=keys[i];
				int x=prefilter.read(key);
				if(x==prefilter.maxValue){
					keys[j]=hash(key, 0);
					j++;
				}
			}
			synchronized(buffers){
				for(int i=0; i<j; i++){
					incrementPartiallyHashed(keys[i]);
				}
			}
		}
	}
	
	//Slow
	/**
	 * Increments a k-mer count by a specified amount using repeated single increments.
	 * This is a slower implementation that calls increment0 multiple times.
	 * @param rawKey The k-mer key to increment
	 * @param amt Amount to increment by
	 */
	@Override
	public void increment(final long rawKey, int amt){
		for(int i=0; i<amt; i++){increment0(rawKey);}
	}
	
	public void increment0(final long rawKey){
		if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}
		if(prefilter!=null){
			int pre=prefilter.read(rawKey);
			if(pre<prefilter.maxValue){return;}
		}
		
		long key2=rawKey;
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
//			assert(readHashed(key2)==0);
			
			int bnum=(int)(key2&arrayMask);
			long[] array=buffers[bnum];
			int loc=bufferlen[bnum];
			array[loc]=key2;
			bufferlen[bnum]++;
			if(verbose){System.err.println("bufferlen["+bnum+"] = "+bufferlen[bnum]);}
			if(bufferlen[bnum]>=array.length){

				if(verbose){System.err.println("Moving array.");}
				bufferlen[bnum]=0;
				buffers[bnum]=new long[array.length];

				writers[bnum].add(array);
				if(verbose){System.err.println("Moved.");}
			}
//			assert(read(rawKey)<=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
//			assert(readHashed(key2)>=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
			key2=Long.rotateRight(key2, hashBits);
		}
	}
	
	private void incrementPartiallyHashed(final long pKey){
		if(verbose){System.err.println("\n*** Incrementing key "+pKey+" ***");}
		
		long key2=pKey;
		
		{
			int bnum=(int)(key2&arrayMask);
			long[] array=buffers[bnum];
			int loc=bufferlen[bnum];
			array[loc]=key2;
			bufferlen[bnum]++;
			if(verbose){System.err.println("bufferlen["+bnum+"] = "+bufferlen[bnum]);}
			if(bufferlen[bnum]>=array.length){

				if(verbose){System.err.println("Moving array.");}
				bufferlen[bnum]=0;
				buffers[bnum]=new long[array.length];

				writers[bnum].add(array);
				if(verbose){System.err.println("Moved.");}
			}
		}
		
		for(int i=1; i<hashes; i++){
			key2=Long.rotateRight(key2, hashBits);
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
//			assert(readHashed(key2)==0);
			
			int bnum=(int)(key2&arrayMask);
			long[] array=buffers[bnum];
			int loc=bufferlen[bnum];
			array[loc]=key2;
			bufferlen[bnum]++;
			if(verbose){System.err.println("bufferlen["+bnum+"] = "+bufferlen[bnum]);}
			if(bufferlen[bnum]>=array.length){

				if(verbose){System.err.println("Moving array.");}
				bufferlen[bnum]=0;
				buffers[bnum]=new long[array.length];

				writers[bnum].add(array);
				if(verbose){System.err.println("Moved.");}
			}
		}
	}
	
	/**
	 * Increment and return operation is not supported for this implementation.
	 * The buffered multi-threaded design makes atomic read-increment-return unfeasible.
	 *
	 * @param key The k-mer key
	 * @param incr Increment amount
	 * @return Never returns as operation throws exception
	 * @throws RuntimeException Always thrown as operation is not supported
	 */
	@Override
	public int incrementAndReturnUnincremented(long key, int incr){
		throw new RuntimeException("Operation not supported.");
	}
	
	/** Converts count data to frequency histogram.
	 * @return Array where index represents count value and value represents frequency */
	@Override
	public long[] transformToFrequency(){
		return transformToFrequency(matrix);
	}
	
	/**
	 * Creates string representation of all count values in the arrays.
	 * Iterates through each array extracting individual cell values for debugging.
	 * @return ByteBuilder containing comma-separated list of all count values
	 */
	@Override
	public ByteBuilder toContentsString(){
		ByteBuilder sb=new ByteBuilder();
		sb.append('[');
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
		sb.append(']');
		return sb;
	}
	
	/** Calculates fraction of cells that contain non-zero counts.
	 * @return Ratio of used cells to total cells */
	@Override
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	/**
	 * Calculates fraction of cells with counts at or above threshold.
	 * @param mindepth Minimum count threshold
	 * @return Ratio of cells meeting threshold to total cells
	 */
	@Override
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	/**
	 * Counts cells with values at or above the specified threshold.
	 * Iterates through all arrays examining each packed cell value.
	 * @param mindepth Minimum count threshold
	 * @return Number of cells meeting the threshold
	 */
	@Override
	public long cellsUsed(int mindepth){
		long count=0;
//		System.out.println("A: "+cellBits+", "+Integer.toBinaryString(valueMask));
		for(int[] array : matrix){
//			System.out.println("B");
			if(array!=null){
//				System.out.println("C");
				for(int word : array){
//					System.out.println("D: "+Integer.toBinaryString(word));
					while(word>0){
						int x=word&valueMask;
//						System.out.println("E: "+x+", "+mindepth);
						if(x>=mindepth){count++;}
						word>>>=cellBits;
					}
				}
			}
		}
		return count;
	}
	
	/**
	 * Hash function that applies row-specific transformations using precomputed masks.
	 * Uses double hashing for the first row to improve distribution quality.
	 *
	 * @param key Input key to hash
	 * @param row Hash function row (0-based)
	 * @return Hashed key value
	 */
	@Override
	final long hash(long key, int row){
		int cell=(int)((Long.MAX_VALUE&key)%(hashArrayLength-1));
//		int cell=(int)(hashCellMask&(key));
		
		if(row==0){//Doublehash only first time
			key=key^hashMasks[(row+4)%hashMasks.length][cell];
			cell=(int)(hashCellMask&(key>>5));
//			cell=(int)(hashCellMask&(key>>hashBits));
//			cell=(int)((Long.MAX_VALUE&key)%(hashArrayLength-1));
		}
		
		return key^hashMasks[row][cell];
	}
	
	/**
	 * Creates randomized hash masks for the hash functions.
	 * Each mask has exactly 16 bits set in both upper and lower 32-bit halves
	 * to ensure good bit distribution. Uses collision detection to prevent
	 * duplicate patterns that could reduce hash quality.
	 *
	 * @param rows Number of hash function rows
	 * @param cols Number of mask values per row
	 * @return 2D array of hash masks with balanced bit patterns
	 */
	private static long[][] makeMasks(int rows, int cols) {
		
		long seed;
		synchronized(KCountArray8MT.class){
			seed=counter;
			counter++;
		}
		
		Timer t=new Timer();
		long[][] r=new long[rows][cols];
		Random randy=Shared.threadLocalRandom(seed);
		for(int i=0; i<r.length; i++){
			fillMasks(r[i], randy);
		}
		t.stop();
		if(t.elapsed>200000000L){System.out.println("Mask-creation time: "+t);}
		return r;
	}
	
	private static void fillMasks(long[] r, Random randy) {
//		for(int i=0; i<r.length; i++){
//			long x=0;
//			while(Long.bitCount(x&0xFFFFFFFF)!=16){
//				x=randy.nextLong();
//			}
//			r[i]=(x&Long.MAX_VALUE);
//		}
		
		final int hlen=(1<<hashBits);
		assert(r.length==hlen);
		int[] count1=new int[hlen];
		int[] count2=new int[hlen];
		final long mask=hlen-1;

		for(int i=0; i<r.length; i++){
			long x=0;
			int y=0;
			int z=0;
			while(Long.bitCount(x&0xFFFFFFFFL)!=16){
				x=randy.nextLong();
				while(Long.bitCount(x&0xFFFFFFFFL)<16){
					x|=(1L<<randy.nextInt(32));
				}
				while(Long.bitCount(x&0xFFFFFFFFL)>16){
					x&=(~(1L<<randy.nextInt(32)));
				}
				while(Long.bitCount(x&0xFFFFFFFF00000000L)<16){
					x|=(1L<<(randy.nextInt(32)+32));
				}
				while(Long.bitCount(x&0xFFFFFFFF00000000L)>16){
					x&=(~(1L<<(randy.nextInt(32)+32)));
				}
				
//				System.out.print(".");
//				y=(((int)(x&mask))^i);
				y=(((int)(x&mask)));
				z=(int)((x>>hashBits)&mask);
				if(count1[y]>0 || count2[z]>0){
					x=0;
				}
			}
//			System.out.println(Long.toBinaryString(x));
			r[i]=(x&Long.MAX_VALUE);
			count1[y]++;
			count2[z]++;
		}
		
	}
	
	
	/** Initializes and starts all WriteThread instances.
	 * Each thread will allocate its portion of the matrix and begin processing. */
	@Override
	public void initialize(){
		for(int i=0; i<writers.length; i++){
			writers[i]=new WriteThread(i);
			writers[i].start();

//			while(!writers[i].isAlive()){
//				System.out.print(".");
//			}
		}
	}
	
	/**
	 * Gracefully shuts down all WriteThread instances and finalizes counts.
	 * Flushes remaining buffers, sends poison messages, and waits for completion.
	 * Aggregates per-thread cell usage statistics into final count.
	 */
	@Override
	public void shutdown(){
		if(finished){return;}
		synchronized(this){
			if(finished){return;}
			
			//Clear buffers
			for(int i=0; i<numArrays; i++){
				long[] array=buffers[i];
				int len=bufferlen[i];
				buffers[i]=null;
				bufferlen[i]=0;
				
				if(len<array.length){
					array=Arrays.copyOf(array, len);
				}
				
				if(array.length>0){
					writers[i].add(array);
				}
			}
			
			//Add poison
			for(WriteThread wt : writers){
				wt.add(poison);
			}
			
			//Wait for termination
			for(WriteThread wt : writers){
//				System.out.println("wt"+wt.num+" is alive: "+wt.isAlive());
				while(wt.isAlive()){
//					System.out.println("wt"+wt.num+" is alive: "+wt.isAlive());
					try {
						wt.join(10000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if(wt.isAlive()){System.err.println(wt.getClass().getCanonicalName()+" is taking a long time to die.");}
				}
				cellsUsed+=wt.cellsUsedPersonal;
//				System.out.println("cellsUsed="+cellsUsed);
			}
			
			assert(!finished);
			finished=true;
		}
	}
	
	private class WriteThread extends Thread{
		
		public WriteThread(int tnum){
			num=tnum;
		}
		
		/**
		 * Main thread execution loop that processes key batches until shutdown.
		 * Allocates local array for NUMA optimization and processes incoming
		 * key arrays from the blocking queue until poison message received.
		 */
		@Override
		public void run(){
			assert(matrix[num]==null);
			array=new int[wordsPerArray]; //Makes NUMA systems use local memory.
			
			matrix[num]=array;
			
			long[] keys=null;
			while(!shutdown){

				if(verbose){System.err.println(" - Reading keys for wt"+num+".");}
				while(keys==null){
					try {
						keys=writeQueue.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if(keys==poison){
//					assert(false);
					shutdown=true;
				}else{
					for(long key : keys){
						incrementHashedLocal(key);
					}
				}
//				System.out.println(" -- Read keys for   wt"+num+". poison="+(keys==poison)+", len="+keys.length);
				if(verbose){System.err.println(" -- Read keys for   wt"+num+". (success)");}
				keys=null;
				if(verbose){System.err.println("shutdown="+shutdown);}
			}

//			System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> I died: "+shutdown+", "+(keys==null)+".");
//			assert(false) : ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> I died: "+shutdown+", "+(keys==null)+".";
			
			array=null;
		}
		
		void add(long[] keys){
//			assert(isAlive());
			assert(!shutdown);
			if(shutdown){return;}
//			assert(keys!=poison);
			if(verbose){System.err.println(" + Adding keys to wt"+num+".");}
			boolean success=false;
			while(!success){
				try {
					writeQueue.put(keys);
					success=true;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(verbose){System.err.println(" ++ Added keys to wt"+num+". (success)");}
		}
		
		private int incrementHashedLocal(long key){
			assert((key&arrayMask)==num);
			key=(key>>>arrayBits)%(cellMod);
//			key=(key>>>(arrayBits+1))%(cellMod);
			int index=(int)(key>>>indexShift);
			int word=array[index];
			int cellShift=(int)(cellBits*key);
			int value=((word>>>cellShift)&valueMask);
			if(value==0){cellsUsedPersonal++;}
			value=min(value+1, maxValue);
			word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
			array[index]=word;
			return value;
		}
		
		private int[] array;
		private final int num;
		public long cellsUsedPersonal=0;
		
		public ArrayBlockingQueue<long[]> writeQueue=new ArrayBlockingQueue<long[]>(16);
		public boolean shutdown=false;
		
	}
	
	
	public long cellsUsed(){return cellsUsed;}
	
	private boolean finished=false;
	
	private long cellsUsed;
	final int[][] matrix;
	private final WriteThread[] writers=new WriteThread[numArrays];
	private final int hashes;
	final int wordsPerArray;
	private final long cellsPerArray;
	final long cellMod;
	private final long[][] hashMasks=makeMasks(8, hashArrayLength);
	
	private final long[][] buffers=new long[numArrays][500];
	private final int[] bufferlen=new int[numArrays];
	
	public final KCountArray prefilter;
	
	private static final int hashBits=6;
	private static final int hashArrayLength=1<<hashBits;
	private static final int hashCellMask=hashArrayLength-1;
	static final long[] poison=new long[0];
	
	private static long counter=0;
	
}
