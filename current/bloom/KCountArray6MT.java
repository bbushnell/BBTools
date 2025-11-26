package bloom;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;


/**
 * 
 * Uses hashing rather than direct-mapping to support longer kmers.
 * 
 * @author Brian Bushnell
 * @date Aug 17, 2012
 *
 */
public class KCountArray6MT extends KCountArray {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -1524266549200637631L;

	/**
	 * Test program demonstrating basic counting operations.
	 * Creates a KCountArray6MT instance and tests increment/read operations on sample keys.
	 * @param args Command-line arguments: cells, bits, gap, hashes
	 */
	public static void main(String[] args){
		long cells=Long.parseLong(args[0]);
		int bits=Integer.parseInt(args[1]);
		int gap=Integer.parseInt(args[2]);
		int hashes=Integer.parseInt(args[3]);
		
		verbose=false;
		
		KCountArray6MT kca=new KCountArray6MT(cells, bits, gap, hashes);
		
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
		
	/**
	 * Constructs a multithreaded k-mer counting array with hash-based storage.
	 * Initializes multiple arrays for parallel processing and sets up hash parameters.
	 *
	 * @param cells_ Total number of cells in the counting array
	 * @param bits_ Bits per cell for storing count values
	 * @param gap_ Gap parameter (currently unused in this implementation)
	 * @param hashes_ Number of hash functions to use for bloom filter behavior
	 */
	public KCountArray6MT(long cells_, int bits_, int gap_, int hashes_){
		super(cells_, bits_);
//		verbose=false;
		long words=cells/cellsPerWord;
		assert(words/numArrays<=Integer.MAX_VALUE);
		wordsPerArray=(int)(words/numArrays);
		cellsPerArray=cells/numArrays;
		cellMod=cellsPerArray-1;
		hashes=hashes_;
//		System.out.println("cells="+cells+", words="+words+", wordsPerArray="+wordsPerArray+", numArrays="+numArrays+", hashes="+hashes);
//		assert(false);
		matrix=new int[numArrays][];
		assert(hashes>0 && hashes<=hashMasks.length);
	}
	
	@Override
	public int read(final long rawKey){
		assert(finished);
		if(verbose){System.err.println("Reading raw key "+rawKey);}

		long key1=hash(rawKey, 3);
		int arrayNum=(int)(key1&arrayMask);
		long key2=hash(rawKey, 0);
		
		int min=readHashed(key2, arrayNum);
		for(int i=1; i<hashes && min>0; i++){
			if(verbose){System.err.println("Reading. i="+i+", key2="+key2);}
			key2=Long.rotateRight(key2, hashBits);
			key2=hash(key2, i);
			if(verbose){System.err.println("Rot/hash. i="+i+", key2="+key2);}
			min=min(min, readHashed(key2, arrayNum));
		}
		return min;
	}
	
	/**
	 * Reads count from a specific hash position in the specified array.
	 * Applies modulo operation and bit shifting to locate the cell value.
	 *
	 * @param key Hashed key value
	 * @param arrayNum Index of the array to read from
	 * @return Count value stored at the hashed position
	 */
	int readHashed(long key, int arrayNum){
		if(verbose){System.err.print("Reading hashed key "+key);}
//		System.out.println("key="+key);
//		int arrayNum=(int)(key&arrayMask);
		key=(key&Long.MAX_VALUE)%(cellMod);
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
	
	@Override
	public void write(final long key, int value){
		throw new RuntimeException("Not allowed for this class.");
	}
	
	//Slow
	@Override
	public void increment(final long rawKey, int amt){
		for(int i=0; i<amt; i++){increment0(rawKey);}
	}
	
	/**
	 * Increments the count for a k-mer by 1 using buffered writes.
	 * Adds the key to a buffer for the appropriate writer thread to process.
	 * Creates new buffer and queues work when buffer becomes full.
	 * @param rawKey The k-mer key to increment
	 */
	public void increment0(final long rawKey){
		if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}

		long key1=hash(rawKey, 3);

		if(verbose){System.err.println("key2="+key1+", value="+read(rawKey));}

		int bnum=(int)(key1&arrayMask);
		long[] array=buffers[bnum];
		int loc=bufferlen[bnum];
		
//		key2=Long.rotateRight(key2, hashBits);
//		array[loc]=key2;
		
		array[loc]=rawKey;
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
	
	/** Returns unincremented value */
	@Override
	public int incrementAndReturnUnincremented(long key, int incr){
		throw new RuntimeException("Operation not supported.");
	}
	
	@Override
	public long[] transformToFrequency(){
		return transformToFrequency(matrix);
	}
	
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
	
	@Override
	public double usedFraction(){return cellsUsed/(double)cells;}
	
	@Override
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	@Override
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
	
	
	@Override
	final long hash(long key, int row){
		int cell=(int)((Long.MAX_VALUE&key)%(hashArrayLength-1));
//		int cell=(int)(hashCellMask&(key));
		
		if(row==0){//Doublehash only first time
			key=key^hashMasks[(row+4)%hashMasks.length][cell];
			cell=(int)(hashCellMask&(key>>4));
//			cell=(int)(hashCellMask&(key>>hashBits));
//			cell=(int)((Long.MAX_VALUE&key)%(hashArrayLength-1));
		}
		
		return key^hashMasks[row][cell];
	}
	
	/**
	 * @param i
	 * @param j
	 * @return
	 */
	private static long[][] makeMasks(int rows, int cols) {
		
		long seed;
		synchronized(KCountArray6MT.class){
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
	
	/**
	 * Fills an array with hash masks having exactly 16 bits set in each 32-bit half.
	 * Ensures uniform bit distribution and uniqueness constraints for hash quality.
	 * @param r Array to fill with hash masks
	 * @param randy Random number generator for mask generation
	 */
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
	
	/**
	 * Worker thread that processes increment operations for a specific array partition.
	 * Receives batched key arrays through a blocking queue and applies increments locally.
	 * Maintains its own array partition to avoid synchronization overhead.
	 */
	private class WriteThread extends Thread{
		
		/** Constructs a WriteThread for processing a specific array partition.
		 * @param tnum Thread number identifying which array partition to handle */
		public WriteThread(int tnum){
			num=tnum;
		}
		
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
						incrementRawLocal(key);
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
		
		/**
		 * Adds a batch of keys to this thread's processing queue.
		 * Blocks if queue is full until space becomes available.
		 * @param keys Array of keys to process for increment operations
		 */
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
		
		/**
		 * Increments count for a raw key using all hash functions in this thread's array.
		 * Applies each hash function and increments the corresponding cell.
		 * Returns the expected minimum count after all increments complete.
		 *
		 * @param rawKey The raw k-mer key to increment
		 * @return Expected minimum count value after increment
		 */
		private int incrementRawLocal(long rawKey){
//			verbose=(rawKey==32662670693L);
			if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}
//			verbose=true;
			assert(1>0);
			
			long key2=rawKey;
			if(hashes==1){
				key2=hash(key2, 0);
//				int x=incrementHashedIfAtMost(key2, 1, maxValue-1);
				int x=incrementHashedLocal(key2);
				assert(x>=1) : "original=?, new should be >="+(1)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
				return x;
			}
			
			int min=0;
//			final int min=read(rawKey);
//			if(min>=maxValue){return maxValue;}
			
			assert(key2==rawKey);
			for(int i=0; i<hashes; i++){
				key2=hash(key2, i);
				if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2, num));}
//				int x=incrementHashedIfAtMost(key2, 1, min);
				int x=incrementHashedLocal(key2);
				assert(x>=min+1) : "i="+i+", original="+min+", new should be <="+(min+1)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
				if(verbose){System.err.println("postIncr value="+readHashed(key2, num));}
//				assert(read(rawKey)<=min+1) : "i="+i+", original="+min+", new should be <="+(min+1)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
//				assert(readHashed(key2)>=min+1) : "i="+i+", original="+min+", new should be <="+(min+1)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
				key2=Long.rotateRight(key2, hashBits);
			}
//			assert(read(rawKey)==min+1) : "original="+min+", new should be "+(min+1)+", new="+read(rawKey)+", max="+maxValue;
//			assert(false);
			return min(min+1, maxValue);
		}
		
		/**
		 * Increments a single hashed position in this thread's local array.
		 * Extracts the current value, increments it, and stores back with bit manipulation.
		 * Tracks cell usage statistics when cells transition from zero to non-zero.
		 *
		 * @param key Hashed key identifying the cell position
		 * @return New count value after increment (capped at maxValue)
		 */
		private int incrementHashedLocal(long key){
//			assert((key&arrayMask)==num);
			key=(key&Long.MAX_VALUE)%(cellMod);
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
		
		/** Local array partition managed by this thread */
		private int[] array;
		/** Thread number identifying which array partition this thread manages */
		private final int num;
		/** Count of non-zero cells in this thread's array partition */
		public long cellsUsedPersonal=0;
		
		/** Queue for receiving key batches to process */
		public ArrayBlockingQueue<long[]> writeQueue=new ArrayBlockingQueue<long[]>(16);
		/** Flag indicating this thread should terminate processing */
		public boolean shutdown=false;
		
	}
	
	
	/** Returns the total number of non-zero cells across all arrays.
	 * @return Total count of cells containing non-zero values */
	public long cellsUsed(){return cellsUsed;}
	
	/** Flag indicating the counting array has been finalized */
	private boolean finished=false;
	
	/** Total number of non-zero cells across all arrays */
	private long cellsUsed;
	/** 2D matrix storing the actual count data partitioned across arrays */
	final int[][] matrix;
	/** Array of worker threads, one per array partition */
	private final WriteThread[] writers=new WriteThread[numArrays];
	/** Number of hash functions used for bloom filter behavior */
	final int hashes;
	/** Number of 32-bit words in each array partition */
	final int wordsPerArray;
	/** Number of cells in each array partition */
	private final long cellsPerArray;
	/** Modulo value for cell addressing within each array */
	final long cellMod;
	/** Precomputed hash masks for all hash functions */
	private final long[][] hashMasks=makeMasks(8, hashArrayLength);
	
	/** Key buffers for batching writes to each thread */
	private final long[][] buffers=new long[numArrays][1000];
	/** Current length of each key buffer */
	private final int[] bufferlen=new int[numArrays];
	
	/** Number of bits used for hash array indexing */
	private static final int hashBits=6;
	/** Length of each hash mask array */
	private static final int hashArrayLength=1<<hashBits;
	/** Bit mask for hash array cell addressing */
	private static final int hashCellMask=hashArrayLength-1;
	/** Sentinel value signaling threads to terminate */
	static final long[] poison=new long[0];
	
	/** Counter for generating unique seeds for hash mask creation */
	private static long counter=0;
	
}
