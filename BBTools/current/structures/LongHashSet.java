package structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public final class LongHashSet{
	
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		LongHashSet set=new LongHashSet(20, 0.7f);
		HashSet<Long> set2=new HashSet<Long>(20, 0.7f);
		ArrayList<Long> list=new ArrayList<Long>();
		ArrayList<Long> list2=new ArrayList<Long>();
		for(long i=0; i<1000; i++){
			assert(!set.contains(i));
			assert(!set2.contains(i));
			list.add(Long.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			long r=randy2.nextLong();
			list2.add(r);
		}
		
		for(long x : list){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		
		for(long x : list){
			assert(set.contains(x));
			assert(set2.contains(x));
			set.remove(x);
			set2.remove(x);
			assert(!set.contains(x));
			assert(!set2.contains(x));
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		
		for(long x : list2){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		
		for(long x : list2){
			assert(set.contains(x));
			assert(set2.contains(x));
			set.remove(x);
			set2.remove(x);
			assert(!set.contains(x));
			assert(!set2.contains(x));
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		
		int count=4000000;
		int runs=32;
		LongList ll=new LongList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextLong());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("LongHashSet:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					set.add(x);
//				}
				final long[] y=ll.array;
				for(int z=0; z<count; z++){
					final long value=y[z];
					set.add(value);
					set.contains(value);
					set.remove(value);
					set.add(value);
				}
//				for(long x : ll.array){
//					set.remove(x);
//				}
//				set.clear();
//				assert(set.isEmpty());
//				System.err.println("Finished run "+i);
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();
			
//			System.err.println("HashSet:");
//			t.start();
//			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					set2.add(x);
//				}
//				for(long x : ll.array){
//					set2.remove(x);
//				}
//				assert(set2.isEmpty());
////				System.err.println("Finished run "+i);
//			}
//			t.stop();
//			System.err.println(t);
//			Shared.printMemory();
		}
		t.stop();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public LongHashSet(){
		this(256);
	}
	
	public LongHashSet(int initialSize){
		this(initialSize, 0.7f);
	}
	
	public LongHashSet(int initialSize, float loadFactor_){
		invalid=randy.nextLong()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0) : "Attempting to initialize a "+getClass().getSimpleName()+" of size<1.";
		assert(loadFactor_>0 && loadFactor_<1) : "Attempting to initialize a "+getClass().getSimpleName()+" with invalid load factor: "+loadFactor_;
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void clear(){
		if(size<1){return;}
		Arrays.fill(array, invalid);
		size=0;
	}
	
	public boolean contains(long value){
		return value==invalid ? false : findCell(value)>=0;
	}
	
	/**
	 * Add this value to the set.
	 * @param value
	 * @return true if the value was added, false if it was already contained.
	 */
	public boolean add(long value){
		if(value==invalid){resetInvalid();}
		int cell=findCellOrEmpty(value);
		if(array[cell]==invalid){
			array[cell]=value;
			size++;
			if(size>sizeLimit){resize();}
			return true;
		}
		assert(array[cell]==value);
		return false;
	}
	
	/**
	 * Remove this value from the set.
	 * @param value
	 * @return true if the value was removed, false if it was not present.
	 */
	public boolean remove(long value){
		if(value==invalid){return false;}
		final int cell=findCell(value);
		if(cell<0){return false;}
		assert(array[cell]==value);
		array[cell]=invalid;
		size--;
		
		rehashFrom(cell);
		return true;
	}
	
	public int size(){return size;}
	
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<array.length; i++){
			if(array[i]!=invalid){
				sb.append(comma+"("+i+", "+array[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<array.length; i++){
			if(array[i]!=invalid){
				sb.append(comma+array[i]);
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public long[] toArray(){
		long[] x=new long[array.length];
		int i=0;
		for(long v : array){
			x[i]=v;
			i++;
		}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean verify(){
		int numValues=0;
		int numFound=0;
		for(int i=0; i<array.length; i++){
			final long value=array[i];
			if(value!=invalid){
				numValues++;
				final int cell=findCell(value);
				if(i==cell){
					numFound++;
				}else{
					return false;
				}
			}
		}
		return numValues==numFound && numValues==size;
	}
	
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=array.length;
		for(int cell=initial+1; cell<limit; cell++){
			final long x=array[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final long x=array[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
	}
	
	private boolean rehashCell(final int cell){
		final long value=array[cell];
		assert(value!=invalid);
		if(value==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(value);
		if(cell==dest){return false;}
		assert(array[dest]==invalid);
		array[cell]=invalid;
		array[dest]=value;
		return true;
	}
	
	private void resetInvalid(){
		final long old=invalid;
		long x=invalid;
		while(x==old || contains(x)){x=randy.nextLong()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<array.length; i++){
			if(array[i]==old){array[i]=invalid;}
		}
	}
	
	private int findCell(final long value){
		if(value==invalid){return -1;}
		
		final int limit=array.length, initial=(int)((value&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=array[cell];
			if(x==value){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=array[cell];
			if(x==value){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}
	
	private int findCellOrEmpty(final long value){
		assert(value!=invalid) : "Collision - this should have been intercepted.";
		
		final int limit=array.length, initial=(int)((value&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=array[cell];
			if(x==value || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=array[cell];
			if(x==value || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}
	
	private final void resize(){
		assert(size>=sizeLimit);
		resize(array.length*2L+1);
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2);
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra);
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;
		
		final int size3=(int)(newPrime+extra);
		sizeLimit=(int)(modulus*loadFactor);
		final long[] old=array;
		array=KillSwitch.allocLong1D(size3);
		Arrays.fill(array, invalid);
		
//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);
		
		if(size<1){return;}
		
		size=0;
		for(long value : old){
			if(value!=invalid){
				add(value);
			}
		}
	}

	public final void resizeDestructive(int newSize){
		size=0;
		sizeLimit=0;
		array=null;
		resize(newSize);
	}
	
	/** Number of values that can be held without resizing */
	public int capacity(){
		return sizeLimit;
	}

	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private long[] array;
	private int size=0;
	/** Value for empty cells */
	private long invalid;
	private int modulus;
	private int sizeLimit;
	private final float loadFactor;
	
	private static final Random randy=new Random(1);
	private static final long MASK=Long.MAX_VALUE;
	private static final long MINMASK=Long.MIN_VALUE;
	
	private static final int extra=10;

	
}
