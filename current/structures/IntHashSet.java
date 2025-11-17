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
 * @date June 7, 2017
 *
 */
public class IntHashSet{
	
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		IntHashSet set=new IntHashSet(20, 0.7f);
		HashSet<Integer> set2=new HashSet<Integer>(20, 0.7f);
		ArrayList<Integer> list=new ArrayList<Integer>();
		ArrayList<Integer> list2=new ArrayList<Integer>();
		for(int i=0; i<1000; i++){
			assert(!set.contains(i));
			assert(!set2.contains(i));
			list.add(Integer.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			int r=randy2.nextInt();
			list2.add(r);
		}
		
		for(int x : list){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		
		for(int x : list){
			assert(set.contains(x));
			assert(set2.contains(x));
			set.remove(x);
			set2.remove(x);
			assert(!set.contains(x));
			assert(!set2.contains(x));
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		
		for(int x : list2){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		
		for(int x : list2){
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
		IntList ll=new IntList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextInt());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("IntHashSet:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(int x : ll.array){
//					set.add(x);
//				}
				final int[] y=ll.array;
				for(int z=0; z<count; z++){
					final int value=y[z];
					set.add(value);
					set.contains(value);
					set.remove(value);
					set.add(value);
				}
//				for(int x : ll.array){
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
//				for(int x : ll.array){
//					set2.add(x);
//				}
//				for(int x : ll.array){
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
	
	/** Creates an IntHashSet with default capacity of 256 and load factor of 0.7 */
	public IntHashSet(){
		this(256);
	}
	
	/** Creates an IntHashSet with specified initial capacity and default load factor of 0.7.
	 * @param initialSize Initial capacity for the hash table */
	public IntHashSet(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Creates an IntHashSet with specified initial capacity and load factor.
	 * Load factor is clamped between 0.25 and 0.90 for optimal performance.
	 * @param initialSize Initial capacity for the hash table
	 * @param loadFactor_ Target load factor before resizing (0.25-0.90)
	 */
	public IntHashSet(int initialSize, float loadFactor_){
		invalid=randy.nextInt()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Removes all elements from the set by filling array with invalid sentinel value
	 */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(array, invalid);
		size=0;
	}
	
	/**
	 * Tests whether the set contains the specified value.
	 * @param value The integer value to search for
	 * @return true if the value is present in the set, false otherwise
	 */
	public boolean contains(int value){
		return value==invalid ? false : findCell(value)>=0;
	}
	
	/**
	 * Add this value to the set.
	 * @param value
	 * @return true if the value was added, false if it was already contained.
	 */
	public boolean add(int value){
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
	 * Add this list to the set.
	 * @param list
	 * @return Number of elements added.
	 */
	public int addAll(IntList list){
		int added=0;
		for(int i=0; i<list.size; i++){
			boolean b=add(list.get(i));
			if(b){added++;}
		}
		return added;
	}
	
	//Does not get added to the list in an IntHashSetList.  For rehashing.  Cannot be overriden.
	/**
	 * Internal add method used during rehashing operations.
	 * Identical to add() but marked final to prevent override in subclasses.
	 * Used specifically for rehashing to avoid triggering list updates in subclasses.
	 *
	 * @param value The integer value to add
	 * @return true if the value was added, false if already contained
	 */
	protected final boolean addSpecial(int value){
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
	public boolean remove(int value){
		if(value==invalid){return false;}
		final int cell=findCell(value);
		if(cell<0){return false;}
		assert(array[cell]==value);
		array[cell]=invalid;
		size--;
		
		rehashFrom(cell);
		return true;
	}

	/** Returns the number of elements currently in the set */
	public int size(){return size;}
	/** Returns the maximum size before next resize based on current load factor */
	public int sizeLimit(){return sizeLimit;}
	
	/** Returns true if the set contains no elements */
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns a string representation showing array positions and values.
	 * Format: [(index, value), ...] for occupied cells only.
	 * @return String showing internal array structure
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
	
	/**
	 * Returns a string representation showing only the stored values.
	 * Format: [value1, value2, ...] similar to standard collections.
	 * @return String showing stored values only
	 */
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
	 * Returns an array containing all elements in the set.
	 * Order is not guaranteed and depends on internal hash table structure.
	 * @return Array containing all set elements
	 */
	public int[] toArray(){
		int[] x=new int[size];
		int i=0;
		for(int v : array){
			if(v!=invalid){
				x[i]=v;
				i++;
			}
		}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Validates internal consistency of the hash table structure.
	 * Checks that all stored values can be found at their expected positions.
	 * Used primarily for debugging and testing.
	 * @return true if the internal structure is consistent, false otherwise
	 */
	public boolean verify(){
		int numValues=0;
		int numFound=0;
		for(int i=0; i<array.length; i++){ //Possible bug: should this be findCell(value) not findCell(i)?
			final int value=array[i];
			if(value!=invalid){
				numValues++;
				final int cell=findCell(i); //Possible bug: searching for index i, should be searching for value?
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
	 * Rehashes elements starting from the specified position after a removal.
	 * Uses wrap-around search to maintain proper hash table invariants.
	 * @param initial Starting position for rehashing operation
	 */
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=array.length;
		for(int cell=initial+1; cell<limit; cell++){
			final int x=array[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final int x=array[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
	}
	
	private boolean rehashCell(final int cell){
		final int value=array[cell];
		assert(value!=invalid);
		if(value==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(value);
		if(cell==dest){return false;}
		assert(array[dest]==invalid);
		array[cell]=invalid;
		array[dest]=value;
		return true;
	}
	
	/**
	 * Generates a new invalid sentinel value and updates all empty cells.
	 * Called when a user attempts to add the current invalid value to avoid conflicts.
	 * Ensures the invalid value is always negative and not already in the set.
	 */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<array.length; i++){
			if(array[i]==old){array[i]=invalid;}
		}
	}
	
	private int findCell(final int value){
		if(value==invalid){return -1;}
		
		final int limit=array.length, initial=(int)((value&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final int x=array[cell];
			if(x==value){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=array[cell];
			if(x==value){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}
	
	private int findCellOrEmpty(final int value){
		assert(value!=invalid) : "Collision - this should have been intercepted.";
		
		final int limit=array.length, initial=(int)((value&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final int x=array[cell];
			if(x==value || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=array[cell];
			if(x==value || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}
	
	/** Doubles the hash table size plus one and rehashes all elements */
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
		final int[] old=array;
		array=KillSwitch.allocInt1D(size3);
		Arrays.fill(array, invalid);
		
//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);
		
		if(size<1){return;}
		
		size=0;
		for(int value : old){
			if(value!=invalid){
				addSpecial(value);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Hash table array storing integer values and invalid sentinels for empty cells
	 */
	private int[] array;
	/** Current number of elements stored in the set */
	private int size=0;
	/** Value for empty cells */
	private int invalid;
	/**
	 * Prime number used as modulus for hash function to distribute values evenly
	 */
	private int modulus;
	/** Maximum number of elements before triggering resize based on load factor */
	private int sizeLimit;
	/** Target ratio of occupied cells to total capacity before resizing */
	private final float loadFactor;
	
	/**
	 * Random number generator for creating invalid sentinel values, seeded for reproducibility
	 */
	private static final Random randy=new Random(1);
	/** Bit mask (Integer.MAX_VALUE) used to ensure hash values are positive */
	private static final int MASK=Integer.MAX_VALUE;
	/**
	 * Bit mask (Integer.MIN_VALUE) used to ensure invalid sentinels are negative
	 */
	private static final int MINMASK=Integer.MIN_VALUE;
	
	/**
	 * Additional capacity beyond prime modulus to reduce collisions at array boundaries
	 */
	private static final int extra=10;
	
}
