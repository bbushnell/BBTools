package bloom;

import java.util.Random;

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
public class KCountArray4 extends KCountArray {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -1418539960644885681L;

	public static void main(String[] args){
		long cells=Long.parseLong(args[0]);
		int bits=Integer.parseInt(args[1]);
		int gap=Integer.parseInt(args[2]);
		int hashes=Integer.parseInt(args[3]);
		
		verbose=false;
		
		KCountArray4 kca=new KCountArray4(cells, bits, gap, hashes);
		
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
		
	public KCountArray4(long cells_, int bits_, int gap_, int hashes_){
		super(cells_, bits_);
		long words=cells/cellsPerWord;
		assert(words/numArrays<=Integer.MAX_VALUE);
		int wordsPerArray=(int)(words/numArrays);
		hashes=hashes_;
//		System.out.println("cells="+cells+", words="+words+", wordsPerArray="+wordsPerArray+", numArrays="+numArrays+", hashes="+hashes);
//		assert(false);
		matrix=new int[numArrays][wordsPerArray];
		assert(hashes>0 && hashes<=hashMasks.length);
	}
	
	@Override
	public int read(final long rawKey){
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
		key=((key&Long.MAX_VALUE)%(cells-1));
//		System.out.println("key="+key);
		int arrayNum=(int)(key&arrayMask);
//		System.out.println("array="+arrayNum);
		key>>>=arrayBits;
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
	
	@Override
	public void increment(final long rawKey, int incr){
//		verbose=(rawKey==32662670693L);
		if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}
//		verbose=true;
		assert(incr>0);
		
		long key2=rawKey;
		if(hashes==1){
			key2=hash(key2, 0);
			int x=incrementHashedIfAtMost(key2, incr, maxValue-1);
			assert(x>=incr) : "original=?, new should be >="+(incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
			//return x;
		}
		
		final int min=read(rawKey);
		if(min>=maxValue){return /*maxValue*/;}
		
		assert(key2==rawKey);
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
			int x=incrementHashedIfAtMost(key2, incr, min);
			assert(x>=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
			if(verbose){System.err.println("postIncr value="+readHashed(key2));}
//			assert(read(rawKey)<=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
//			assert(readHashed(key2)>=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
			key2=Long.rotateRight(key2, hashBits);
		}
//		assert(read(rawKey)==min+incr) : "original="+min+", new should be "+(min+incr)+", new="+read(rawKey)+", max="+maxValue;
//		assert(false);
		//return min(min+incr, maxValue);
	}
	
	/** Returns unincremented value */
	@Override
	public int incrementAndReturnUnincremented(long rawKey, int incr){
//		verbose=(rawKey==32662670693L);
		if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}
//		verbose=true;
		assert(incr>0);
		
		long key2=rawKey;
		if(hashes==1){
			key2=hash(key2, 0);
			int x=incrementHashedIfAtMost(key2, incr, maxValue-1);
			assert(x>=incr) : "original=?, new should be >="+(incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
			return x;
		}
		
		final int min=read(rawKey);
		if(min>=maxValue){return maxValue;}
		
		assert(key2==rawKey);
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
			int x=incrementHashedIfAtMost(key2, incr, min);
			assert(x>=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
			if(verbose){System.err.println("postIncr value="+readHashed(key2));}
//			assert(read(rawKey)<=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
//			assert(readHashed(key2)>=min+incr) : "i="+i+", original="+min+", new should be <="+(min+incr)+", new="+read(rawKey)+", max="+maxValue+", key="+rawKey;
			key2=Long.rotateRight(key2, hashBits);
		}
//		assert(read(rawKey)==min+incr) : "original="+min+", new should be "+(min+incr)+", new="+read(rawKey)+", max="+maxValue;
//		assert(false);
		return min;
	}
	
	private int incrementHashedIfAtMost(long key, int incr, int lim){
		if(verbose){System.err.print("incrementing hashed key "+key);}
		key=((key&Long.MAX_VALUE)%(cells-1));
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		int value=((word>>>cellShift)&valueMask);
		if(verbose){System.err.println(", array="+arrayNum+", index="+index+", cellShift="+(cellShift%32)+", value="+value+", limit="+lim);}
		if(value>lim){return value;}
		if(value==0 && incr>0){cellsUsed++;}
		value=min(value+incr, maxValue);
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		return value;
	}
	
	private int incrementHashed(long key, int incr){
		assert(incr>0);
		int arrayNum=(int)(key&arrayMask);
		key>>>=arrayBits;
		int[] array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
		int word=array[index];
		int cellShift=(int)(cellBits*key);
		int value=((word>>>cellShift)&valueMask);
		if(value==0 && incr>0){cellsUsed++;}
		value=min(value+incr, maxValue);
		word=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		array[index]=word;
		return value;
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
		synchronized(KCountArray4.class){
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
	
	public long cellsUsed(){return cellsUsed;}

	private long cellsUsed;
	private final int[][] matrix;
	private final int hashes;
	
	
	private static final int hashBits=6;
	private static final int hashArrayLength=1<<hashBits;
	private static final int hashCellMask=hashArrayLength-1;
	private final long[][] hashMasks=makeMasks(8, hashArrayLength);
	
	private static long counter=0;
	
}
