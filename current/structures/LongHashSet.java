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
	
	/**
	 * Creates a LongHashSet with default initial size of 256 and load factor 0.7
	 */
	public LongHashSet(){
		this(256);
	}
	
	/** Creates a LongHashSet with specified initial size and default load factor 0.7.
	 * @param initialSize Initial capacity of the hash table */
	public LongHashSet(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Creates a LongHashSet with specified initial size and load factor.
	 * Load factor is constrained to range [0.25, 0.90] for performance.
	 * @param initialSize Initial capacity of the hash table
	 * @param loadFactor_ Target load factor before resizing (0.0-1.0)
	 */
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
	
	/** Removes all values from the set by filling the array with invalid markers */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(array, invalid);
		size=0;
	}
	
	/**
	 * Tests whether the set contains the specified value.
	 * @param value The value to search for
	 * @return true if the value is present in the set, false otherwise
	 */
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
	
	/** Returns the number of values currently stored in the set */
	public int size(){return size;}
	
	/** Returns true if the set contains no values */
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns a detailed string showing array indices and their values.
	 * Useful for debugging hash table structure and collision patterns.
	 * @return String in format [(index, value), (index, value), ...]
	 */
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
	
	/** Returns a string representation showing only the stored values.
	 * @return String in format [value1, value2, value3, ...] */
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
	
	/**
	 * Returns a copy of the internal array including both valid values and invalid markers.
	 * The returned array has the same length as the internal hash table.
	 * @return Array copy containing all stored values and invalid markers
	 */
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
	
	/**
	 * Verifies the internal consistency of the hash table structure.
	 * Checks that all values can be found at their expected positions.
	 * @return true if the hash table structure is valid, false otherwise
	 */
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
	
	/**
	 * Rehashes all elements starting from the specified position to maintain
	 * linear probing sequence integrity after a removal operation.
	 * @param initial The starting position for rehashing
	 */
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
	
	/**
	 * Rehashes a single cell to its correct position.
	 * Used during removal operations to fill gaps in the probe sequence.
	 * @param cell The cell index to rehash
	 * @return true if the cell was moved to a different position
	 */
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
	
	/** Generates a new invalid marker value when collision occurs with current marker.
	 * Updates all existing invalid markers in the array to the new value. */
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
	
	/**
	 * Locates the cell containing the specified value using linear probing.
	 * @param value The value to locate
	 * @return Cell index if found, -1 if not found
	 */
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
	
	/**
	 * Finds the cell containing the value or the first empty cell for insertion.
	 * Uses linear probing starting from the hash position.
	 *
	 * @param value The value to locate or position for insertion
	 * @return Cell index containing the value or first empty cell
	 * @throws RuntimeException if no empty cells are found
	 */
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
	
	/** Doubles the hash table size plus one when size limit is reached */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(array.length*2L+1);
	}
	
	/**
	 * Resizes the hash table to accommodate the specified minimum size.
	 * Uses prime numbers for optimal hash distribution and rehashes all existing values.
	 * @param size2 Minimum required size for the new hash table
	 */
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

	/**
	 * Clears the set and resizes to the specified capacity.
	 * All existing values are lost in this operation.
	 * @param newSize New capacity for the hash table
	 */
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
	
	/** Internal hash table array storing the long values */
	private long[] array;
	/** Current number of values stored in the set */
	private int size=0;
	/** Value for empty cells */
	private long invalid;
	/** Prime number used as modulus for hash function calculations */
	private int modulus;
	/** Maximum number of values before triggering resize operation */
	private int sizeLimit;
	/** Load factor threshold for determining when to resize the hash table */
	private final float loadFactor;
	
	/** Random number generator used for generating invalid marker values */
	private static final Random randy=new Random(1);
	/** Bit mask for ensuring positive hash values (Long.MAX_VALUE) */
	private static final long MASK=Long.MAX_VALUE;
	/** Bit mask for ensuring invalid markers are negative (Long.MIN_VALUE) */
	private static final long MINMASK=Long.MIN_VALUE;
	
	/** Additional buffer size added to hash table capacity */
	private static final int extra=10;

	
}
