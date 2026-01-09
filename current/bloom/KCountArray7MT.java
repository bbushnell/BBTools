package bloom;

import java.util.Arrays;
import shared.Random;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;


/**
 * Multi-threaded k-mer count array using prime-sized array dimensions
 * for improved hash distribution.
 * Employs separate write threads for each array partition to enable
 * concurrent increment operations.
 * Extends KCountArray with multi-threading support and uses prime number
 * array lengths to reduce hash collisions.
 *
 * @author Brian Bushnell
 * @date Aug 17, 2012
 */
public class KCountArray7MT extends KCountArray {
	
	/** Serial version UID for serialization compatibility */
	private static final long serialVersionUID = -8767643111803866913L;

	public static void main(String[] args){
		long cells=Long.parseLong(args[0]);
		int bits=Integer.parseInt(args[1]);
		int hashes=Integer.parseInt(args[2]);
		
		verbose=false;
		
		KCountArray7MT kca=new KCountArray7MT(cells, bits, hashes);
		
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
		
	public KCountArray7MT(long cells_, int bits_, int hashes_){
		super(getPrimeCells(cells_, bits_), bits_, getDesiredArrays(cells_, bits_));
//		verbose=false;
//		assert(false);
//		System.out.println(cells);
		cellsPerArray=cells/numArrays;
		wordsPerArray=(int)((cellsPerArray%cellsPerWord)==0 ? (cellsPerArray/cellsPerWord) : (cellsPerArray/cellsPerWord+1));
		cellMod=cellsPerArray;
		hashes=hashes_;
//		System.out.println("cells="+cells+", words="+words+", wordsPerArray="+wordsPerArray+", numArrays="+numArrays+", hashes="+hashes);
//		assert(false);
		matrix=new int[numArrays][];
		assert(hashes>0 && hashes<=hashMasks.length);
	}
	
	private static int getDesiredArrays(long desiredCells, int bits){
		
		long words=Tools.max((desiredCells*bits+31)/32, minArrays);
		int arrays=minArrays;
		while(words/arrays>=Integer.MAX_VALUE){
			arrays*=2;
		}
//		assert(false) : arrays;
		return arrays;
	}
	
	private static long getPrimeCells(long desiredCells, int bits){
		
		int arrays=getDesiredArrays(desiredCells, bits);
		
		long x=(desiredCells+arrays-1)/arrays;
		long x2=Primes.primeAtMost(x);
		return x2*arrays;
	}
	
	/**
	 * Reads the minimum count value across all hash functions for a key.
	 * Applies multiple hash functions and returns the minimum count found,
	 * implementing bloom filter semantics where the minimum represents the true count.
	 *
	 * @param rawKey The raw key to look up
	 * @return Minimum count across all hash functions
	 */
	@Override
	public int read(final long rawKey){
		assert(finished);
		if(verbose){System.err.println("Reading raw key "+rawKey);}
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
	 * Write operation not supported for this multi-threaded implementation.
	 * Use increment operations instead for thread-safe count modifications.
	 *
	 * @param key The key to write to
	 * @param value The value to write
	 * @throws RuntimeException Always thrown as operation is not allowed
	 */
	@Override
	public void write(final long key, int value){
		throw new RuntimeException("Not allowed for this class.");
	}
	
	/**
	 * Increments counts for multiple keys with pre-hashing optimization.
	 * Pre-computes first hash outside the critical section for better performance,
	 * then processes all keys within a synchronized block.
	 * @param keys Array of keys to increment
	 */
	@Override
	/** This should increase speed by doing the first hash outside the critical section, but it does not seem to help. */
	public void increment(long[] keys){
		for(int i=0; i<keys.length; i++){
			keys[i]=hash(keys[i], 0);
		}
		synchronized(buffers){
			for(long key : keys){
				incrementPartiallyHashed(key);
			}
		}
	}
	
	//Slow
	/**
	 * Increments a key's count by the specified amount.
	 * Implemented by calling increment0 multiple times, which is inefficient
	 * but maintains correctness for the multi-threaded implementation.
	 *
	 * @param rawKey The key to increment
	 * @param amt Number of times to increment
	 */
	@Override
	public void increment(final long rawKey, int amt){
		for(int i=0; i<amt; i++){increment0(rawKey);}
	}
	
	public void increment0(final long rawKey){
		if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}
		
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
	 * Increment and return operation not supported for this implementation.
	 * The multi-threaded nature makes atomic increment-and-return operations complex.
	 *
	 * @param key The key to increment
	 * @param incr Increment amount
	 * @return Never returns
	 * @throws RuntimeException Always thrown as operation is not supported
	 */
	@Override
	public int incrementAndReturnUnincremented(long key, int incr){
		throw new RuntimeException("Operation not supported.");
	}
	
	/**
	 * Transforms count data to frequency histogram.
	 * Delegates to the parent class implementation using this instance's matrix.
	 * @return Array containing frequency counts for each possible count value
	 */
	@Override
	public long[] transformToFrequency(){
		return transformToFrequency(matrix);
	}
	
	/**
	 * Creates a string representation of all array contents.
	 * Iterates through all arrays and cells, extracting individual count values
	 * and formatting them as a comma-separated list enclosed in brackets.
	 * @return ByteBuilder containing formatted array contents
	 */
	@Override
	public ByteBuilder toContentsString(){
		ByteBuilder sb=new ByteBuilder();
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
		return sb;
	}
	
	/** Calculates the fraction of cells that contain non-zero counts.
	 * @return Proportion of cells in use (0.0 to 1.0) */
	@Override
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	/**
	 * Calculates the fraction of cells with counts at or above the specified depth.
	 * @param mindepth Minimum count threshold
	 * @return Proportion of cells meeting the depth requirement
	 */
	@Override
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	/**
	 * Counts cells with count values at or above the specified minimum depth.
	 * Iterates through all arrays and extracts individual cell counts from packed words,
	 * counting those that meet the depth threshold.
	 *
	 * @param mindepth Minimum count value to include
	 * @return Number of cells with counts >= mindepth
	 */
	@Override
	public long cellsUsed(int mindepth){
		long count=0;
//		System.out.println("A");
		for(int[] array : matrix){
//			System.out.println("B");
			if(array!=null){
//				System.out.println("C");
				for(int word : array){
//					System.out.println("D: "+Integer.toHexString(word));
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
	 * Applies hash function for the specified row to distribute keys across arrays.
	 * Uses double hashing on the first row (row 0) for improved distribution,
	 * then XORs with pre-computed hash masks for subsequent rows.
	 *
	 * @param key The key to hash
	 * @param row Hash function row/iteration number
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
	 * Creates hash mask arrays with balanced bit patterns for hash distribution.
	 * Generates random masks where each 32-bit half has exactly 16 set bits,
	 * ensuring no collision in the lower hash bits across all masks.
	 *
	 * @param rows Number of hash function rows
	 * @param cols Number of hash mask columns
	 * @return 2D array of hash masks with balanced bit patterns
	 */
	private static long[][] makeMasks(int rows, int cols) {
		
		long seed;
		synchronized(KCountArray7MT.class){
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
	
	
	/** Starts all write threads for concurrent array processing.
	 * Creates and launches one write thread per array partition. */
	@Override
	public void initialize(){
		for(int i=0; i<writers.length; i++){
			writers[i]=new WriteThread(i);
			writers[i].start();

//			while(!writers[i].isAlive()){
//				System.out.print(".");
//			}
		}
//		assert(false) : writers.length;
	}
	
	/**
	 * Shuts down all write threads and finalizes count data.
	 * Flushes remaining buffered operations, sends poison pills to threads,
	 * waits for thread termination, and accumulates final cell usage counts.
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
//				System.err.println("wt.cellsUsedPersonal="+wt.cellsUsedPersonal);
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
		 * Main write thread execution loop.
		 * Allocates local array memory, processes key batches from the queue,
		 * and performs local increment operations until shutdown signal received.
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
	
	private static final int hashBits=6;
	private static final int hashArrayLength=1<<hashBits;
	private static final int hashCellMask=hashArrayLength-1;
	static final long[] poison=new long[0];
	
	private static long counter=0;
	
}
