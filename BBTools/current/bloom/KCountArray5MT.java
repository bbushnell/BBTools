package bloom;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Shared;
import shared.Timer;


/**
 * 
 * Uses hashing rather than direct-mapping to support longer kmers.
 * 
 * @author Brian Bushnell
 * @date Aug 17, 2012
 *
 */
public class KCountArray5MT extends KCountArray {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -5456926900022701212L;

	public static void main(String[] args){
		long cells=Long.parseLong(args[0]);
		int bits=Integer.parseInt(args[1]);
		int gap=Integer.parseInt(args[2]);
		int hashes=Integer.parseInt(args[3]);
		
		verbose=false;
		
		KCountArray5MT kca=new KCountArray5MT(cells, bits, gap, hashes);
		
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
		
	public KCountArray5MT(long cells_, int bits_, int gap_, int hashes_){
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
	
	int readHashed(long key){
		assert(finished);
		if(verbose){System.err.print("Reading hashed key "+key);}
//		System.out.println("key="+key);
		int arrayNum=(int)(key&arrayMask);
		key=(key>>>arrayBits)%(cellMod);
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

		buffer[bufferlen]=hash(rawKey, 0);
		bufferlen++;

		if(bufferlen>=buffer.length){

			if(verbose){System.err.println("Moving array.");}

			for(int w=0; w<writers.length; w++){
				writers[w].add(buffer);
			}
			bufferlen=0;
			buffer=new long[buffer.length];
			if(verbose){System.err.println("Moved.");}
		}

	}
	

	@Override
	public synchronized void increment(long[] keys){
		for(int i=0; i<keys.length; i++){
			keys[i]=hash(keys[i],0);
		}
		for(int w=0; w<writers.length; w++){
			writers[w].add(keys);
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
	public String toContentsString(){
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
		synchronized(KCountArray5MT.class){
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

			while(!writers[i].isAlive()){
				System.out.print(".");
			}
		}
	}
	
	@Override
	public void shutdown(){
		if(finished){return;}
		synchronized(this){
			if(finished){return;}
			
			//Clear buffer
			if(bufferlen<buffer.length){
				buffer=Arrays.copyOf(buffer, bufferlen);
			}
			
			if(buffer.length>0){
				for(int i=0; i<writers.length; i++)
				writers[i].add(buffer);
			}
			buffer=null;
			bufferlen=0;
			
			
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
//			assert(false);
			matrix[num]=array;
			
//			assert(num==1);
			
			long[] keys=null;
			while(!shutdown){
//				assert(false);

				if(verbose){System.err.println(" - Reading keys for wt"+num+".");}
				while(keys==null){
//					System.out.println("Searching for keys.");
					try {
						keys=writeQueue.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
//					System.out.println("*******************************************Found keys: "+keys.length);
//					assert(false);
				}
				if(keys==poison){
//					assert(false) : "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~0  ";
					shutdown=true;
				}else{
//					assert(false) : "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1  ";
					for(long rawKey : keys){
//						assert(false) : "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2  ";
						if(verbose){System.err.println("Writer "+num+" considering raw key "+rawKey);}
						long key2=rawKey;
//						int y=read(rawKey);
						if((key2&arrayMask)==num){
							int x=incrementHashedLocal(key2);
							assert(x>=0) : "i="+0+", original=?, new should be >=0, new="+readHashed(key2)+", max="+maxValue+", key="+rawKey;
							if(verbose){System.err.println("postIncr value="+readHashed(key2));}

//							assert(false) : "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~4  ";
						}
						
						for(int i=1; i<hashes; i++){
							key2=Long.rotateRight(key2, hashBits);
							key2=hash(key2, i);
//							assert(false) : "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3  ";
							if(verbose){System.err.println("rawKey="+rawKey+", i="+i+", key2="+key2+", value="+readHashed(key2));}
							if((key2&arrayMask)==num){
								int x=incrementHashedLocal(key2);
								assert(x>=0) : "i="+i+", original=?, new should be >=0, new="+readHashed(key2)+", max="+maxValue+", key="+rawKey;
								if(verbose){System.err.println("postIncr value="+readHashed(key2));}

//								assert(false) : "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~4  ";
							}

//							assert(false) : "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5  ";
						}
//						int z=read(rawKey);
//						assert(hashes!=1 || !b || z==maxValue || z==y+1) : "b="+b+", y="+y+", z="+z+", rawKey="+rawKey+", num="+num;
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
			
//			assert(!shutdown);
//			if(shutdown){return;}
			
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
		
		private int incrementHashedLocal(final long key_){
			if(verbose){System.err.println("\n*** wt"+num+" incrementing hashed key "+key_+" ***");}
			assert((key_&arrayMask)==num);
			long key=(key_>>>arrayBits)%(cellMod);
			int index=(int)(key>>>indexShift);
			int word=array[index];
			int cellShift=(int)(cellBits*key);
			int value=((word>>>cellShift)&valueMask);
			if(value==0){cellsUsedPersonal++;}
			value=min(value+1, maxValue);
			word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
			array[index]=word;
			if(verbose){System.err.println("\n*** wt"+num+" Incremented hashed key "+key_+".  Value = "+readHashed(key_)+" ***");}
			return value;
		}
		
		private int[] array;
		private final int num;
		public long cellsUsedPersonal=0;
		
		public ArrayBlockingQueue<long[]> writeQueue=new ArrayBlockingQueue<long[]>(8);
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
	
	private long[] buffer=new long[2000];
	private int bufferlen=0;
	
	private static final int hashBits=6;
	private static final int hashArrayLength=1<<hashBits;
	private static final int hashCellMask=hashArrayLength-1;
	static final long[] poison=new long[0];
	
	private static long counter=0;
	
}
