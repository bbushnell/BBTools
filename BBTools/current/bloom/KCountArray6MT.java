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
	
	private class WriteThread extends Thread{
		
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
	final int hashes;
	final int wordsPerArray;
	private final long cellsPerArray;
	final long cellMod;
	private final long[][] hashMasks=makeMasks(8, hashArrayLength);
	
	private final long[][] buffers=new long[numArrays][1000];
	private final int[] bufferlen=new int[numArrays];
	
	private static final int hashBits=6;
	private static final int hashArrayLength=1<<hashBits;
	private static final int hashCellMask=hashArrayLength-1;
	static final long[] poison=new long[0];
	
	private static long counter=0;
	
}
