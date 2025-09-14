package kmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Primes;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Stores kmers in a long[] and values in an int[][], with a victim cache.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
 */
public abstract class HashArray extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	HashArray(int[] schedule_, long coreMask_, boolean twod_){
		schedule=schedule_;
		autoResize=schedule.length>1;
		prime=schedule[0];
		
		sizeLimit=(long)((schedule.length==1 ? maxLoadFactorFinal : maxLoadFactor)*prime);
		array=allocLong1D(prime+extra);
		victims=new HashForest(Tools.max(10, prime/victimRatio), autoResize, twod_);
		Arrays.fill(array, NOT_PRESENT);
		twoD=twod_;
		coreMask=coreMask_;
//		coreMask2=coreMask_|3;
		coreMask2=coreMask_; //Simplifies fast fill
	}
	
	HashArray(int initialSize, long coreMask_, boolean autoResize_, boolean twod_){
		schedule=null;
		prime=initialSize;
		sizeLimit=(long)(maxLoadFactor*prime);
		array=allocLong1D(prime+extra);
		victims=new HashForest(Tools.max(10, initialSize/victimRatio), autoResize_, twod_);
		Arrays.fill(array, NOT_PRESENT);
		autoResize=autoResize_;
		twoD=twod_;
		coreMask=coreMask_;
//		coreMask2=coreMask_|3;
		coreMask2=coreMask_; //Simplifies fast fill
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Test Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			int[] old=getValues(kmer, new int[1]);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Fetched "+Arrays.toString(old));}
//			x=set0(kmer, v);
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+v+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int set_Test(final long kmer, final int v[]){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			final int[] singleton=new int[1];
//			int[] old=getValues(kmer, singleton);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Before: old="+Arrays.toString(old)+", v="+Arrays.toString(v));}
//			x=set0(kmer, v);
//			if(verbose){System.err.println("After:  old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", get()="+Arrays.toString(getValues(kmer, singleton)));}
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v[0]==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
////			int[] vals=getValues(kmer, null);
////			assert(vals==null || contains(kmer, vals));
////			x=setIfNotPresent(kmer, v);
////			assert(contains(kmer, vals));
////			assert(contains(kmer, v));
//			x=0;
//			assert(false);
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=setIfNotPresent0(kmer, v);
//			assert((old<1 && contains(kmer, v)) || (old>0 && contains(kmer, old))) : kmer+", "+old+", "+v;
//		}
//		return x;
//	}
	
	/*--------------------------------------------------------------*/
	
	public final int kmerToCell(long kmer){
		final int cell=(int)((kmer&coreMask)%prime);
		return cell;
	}
	
	@Override
	public final int set(final long kmer, final int[] v, final int vlen){
		int cell=kmerToCell(kmer);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				if(verbose){System.err.println("A2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
				insertValue(kmer, v, cell, vlen);
				if(verbose){System.err.println("A2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				return 0;
			}else if(n==NOT_PRESENT){
				if(verbose){System.err.println("B2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
				array[cell]=kmer;
				insertValue(kmer, v, cell, vlen);
				if(verbose){System.err.println("B2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				size++;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
		if(verbose){System.err.println("C2: Adding "+kmer+", "+v+", "+cell);}
		final int x=victims.set(kmer, v, vlen);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		if(verbose){System.err.println("C2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
		return x;
	}
	
	@Override
	public final int set(final long kmer, final int v){
		int cell=kmerToCell(kmer);
		
//		assert(TESTMODE);
//		ll.add(kmer);
//		il.add(v);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				if(verbose){System.err.println("A1: Adding "+kmer+", "+v+", "+cell);}
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("A1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				return 0;
			}else if(n==NOT_PRESENT){
				if(verbose){System.err.println("B1: Adding "+kmer+", "+v+", "+cell);}
				array[cell]=kmer;
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("B1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				size++;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
		if(verbose){System.err.println("C1: Adding "+kmer+", "+v+", "+cell+
				"; victims.get(kmer)="+Arrays.toString(victims.getValues(kmer, new int[1])));}
		final int x=victims.set(kmer, v);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		if(verbose){System.err.println("C1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1]))+
				"; victims.get(kmer)="+Arrays.toString(victims.getValues(kmer, new int[1])));}
		return x;
	}
	
	@Override
	public final int setIfNotPresent(long kmer, int value){
		int cell=kmerToCell(kmer);
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				return 0;
			}else if(n==NOT_PRESENT){
				array[cell]=kmer;
				insertValue(kmer, value, cell);
				size++;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
//		System.err.println("size="+size+", prime="+prime+", limit="+sizeLimit);
		int x=victims.setIfNotPresent(kmer, value);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public final int getValue(long kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}
	
	public final int getValue(long kmer, int startCell){
		int cell=findKmer(kmer, startCell);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}
	
	@Override
	public final int[] getValues(long kmer, int[] singleton){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){
			singleton[0]=NOT_PRESENT;
			return singleton;
		}
		if(cell==HASH_COLLISION){return victims.getValues(kmer, singleton);}
		return readCellValues(cell, singleton);
	}
	
	@Override
	public final boolean contains(long kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return false;}
		if(cell==HASH_COLLISION){return victims.contains(kmer);}
		return true;
	}
	
	public final long getKmer(int cell) {
		return array[cell];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		assert(owners==null);
		owners=allocAtomicInt(array.length);
		for(int i=0; i<array.length; i++){
			owners.set(i, NO_OWNER);
		}
		victims.initializeOwnership();
	}
	
	@Override
	public final void clearOwnership(){
		owners=null;
		victims.clearOwnership();
	}
	
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.setOwner(kmer, newOwner);}
		return setOwner(kmer, newOwner, cell);
	}
	
	public final int setOwner(final long kmer, final int newOwner, final int cell){
		assert(array[cell]==kmer);
		final int original=owners.get(cell);
		int current=original;
		while(current<newOwner){
			boolean success=owners.compareAndSet(cell, current, newOwner);
			if(!success){current=owners.get(cell);}
			else{current=newOwner;}
		}
		assert(current>=original) : "original="+original+", current="+current+", newOwner="+newOwner+", re-read="+owners.get(cell);
		return current;
	}
	
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.clearOwner(kmer, owner);}
		return clearOwner(kmer, owner, cell);
	}
	
	public final boolean clearOwner(final long kmer, final int owner, final int cell){
		assert(array[cell]==kmer);
		boolean success=owners.compareAndSet(cell, owner, NO_OWNER);
		return success;
	}
	
	@Override
	public final int getOwner(final long kmer){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.getOwner(kmer);}
		return getCellOwner(cell);
	}
	
	public final int getCellOwner(final int cell){
		return owners.get(cell);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	protected abstract void insertValue(final long kmer, final int v, final int cell);

//	protected abstract void insertValue(final long kmer, final int[] vals, final int cell);
	
	/** This is for IntList3 support with HashArrayHybridFast */
	protected abstract void insertValue(final long kmer, final int[] vals, final int cell, final int vlen);
	
	protected abstract int readCellValue(int cell);
	protected abstract int[] readCellValues(int cell, int[] singleton);
	
	@Override
	final Object get(long kmer){
		throw new RuntimeException("Unimplemented.");
	}
	
	final int findKmer(long kmer){
		return findKmer(kmer, kmerToCell(kmer));
	}
	
	final int findKmer(final long kmer, final int startCell){
		int cell=startCell;
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array[cell];
			if(n==kmer){return cell;}
			else if(n==NOT_PRESENT){return NOT_PRESENT;}
		}
		return HASH_COLLISION;
	}
	
	final int findKmerOrEmpty(long kmer){
		int cell=kmerToCell(kmer);
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array[cell];
			if(n==kmer || n==NOT_PRESENT){return cell;}
		}
		return HASH_COLLISION;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	final boolean canResize() {return true;}
	
	@Override
	final public long size() {return size;}
	
	@Override
	final public int arrayLength() {return array.length;}
	
	@Override
	protected abstract void resize();
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		if(twoD){
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT){
					tsw.print(toText(kmer, readCellValues(i, singleton), k).append('\n'));
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT && (mincount<2 || readCellValue(i)>=mincount)){
					tsw.print(toText(kmer, readCellValue(i), k).append('\n'));
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}
	
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(twoD){
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					bsw.printlnKmer(kmer, readCellValues(i, singleton), k);
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT && (mincount<2 || readCellValue(i)>=mincount)){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					bsw.printlnKmer(kmer, readCellValue(i), k);
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}
	
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		if(twoD){
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					toBytes(kmer, readCellValues(i, singleton), k, bb);
					bb.nl();
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT && (mincount<2 || readCellValue(i)>=mincount)){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					toBytes(kmer, readCellValue(i), k, bb);
					bb.nl();
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
		}
		return true;
	}
	
	@Override
	public final void fillHistogram(long[] ca, int max){
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=NOT_PRESENT){
				int count=Tools.min(readCellValue(i), max);
				ca[count]++;
			}
		}
		if(victims!=null){
			victims.fillHistogram(ca, max);
		}
	}
	
	@Override
	public void fillHistogram(SuperLongList sll){
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=NOT_PRESENT){
				int count=readCellValue(i);
				sll.add(count);
			}
		}
		if(victims!=null){
			victims.fillHistogram(sll);
		}
	}
	
	@Override
	public final void countGC(long[] gcCounts, int max){
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=NOT_PRESENT){
				int count=readCellValue(i);
				int index=Tools.min(count, max);
				gcCounts[index]+=gc(kmer);
			}
		}
		if(victims!=null){
			victims.countGC(gcCounts, max);
		}
	}
	
	public HashForest victims(){
		return victims;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	AtomicIntegerArray owners;
	long[] array;
	int prime;
	long size=0;
	long sizeLimit;
	final HashForest victims;
	final boolean autoResize;
	public final boolean twoD;
	private final Lock lock=new ReentrantLock();
	private final long coreMask;//for ways
	private final long coreMask2;//for cells
	
	protected int nextScheduleSize(){
		if(schedulePos<schedule.length-1){schedulePos++;}
		return schedule[schedulePos];
	}
	
	protected boolean atMaxSize(){
		return schedulePos>=schedule.length-1;
	}
	
	protected final int[] schedule;
	private int schedulePos=0;
	
	public long[] array(){return array;}
	
	public AtomicIntegerArray owners() {return owners;}
	@Override
	final Lock getLock(){return lock;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final static int victimRatio=16; //Initial divisor for victim cache size; it self-resizes.
	final static int extra=60; //Amazingly, increasing this gave increasing returns past 300.  Old default was 21.  Could allow higher maxLoadFactorFinal and smaller victim cache.
	final static int maxPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra-20);
	final static float resizeMult=2f; //Resize by a minimum of this much; not needed for schedule
	final static float minLoadFactor=0.58f; //Resize by enough to get the load above this factor; not needed for schedule
	final static float maxLoadFactor=0.88f; //Reaching this load triggers resizing
	final static float maxLoadFactorFinal=0.95f; //Reaching this load triggers killing
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	
}
