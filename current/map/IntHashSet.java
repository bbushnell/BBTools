package map;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import shared.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntList;

/**
 * Hash set implementation optimized for integer values using open addressing and linear probing.
 * Maintains a randomly chosen negative sentinel value for empty cells and dynamically resizes based on load factor.
 * Every int is storable: inserting the current sentinel replaces that marker.
 * Instances are unsynchronized; mutations require exclusive access. Resizing
 * preserves contents but can change the backing-table iteration order.
 * @author Brian Bushnell
 * @date June 7, 2017
 */
public class IntHashSet{

	/** Runs legacy assertion checks followed by a large, fixed-size stress benchmark.
	 * @param args Unused */
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
				final int[] y=ll.array;
				for(int z=0; z<count; z++){
					final int value=y[z];
					set.add(value);
					set.contains(value);
					set.remove(value);
					set.add(value);
				}
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();

		}
		t.stop();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty set with default capacity. */
	public IntHashSet(){
		this(256);
	}

	/** Creates an empty set with the specified initial capacity.
	 * @param initialSize Positive minimum prime modulus before adding slack cells */
	public IntHashSet(int initialSize){
		this(initialSize, 0.7f);
	}

	/** Creates an empty set with the specified capacity and load factor.
	 * @param initialSize Positive minimum prime modulus before adding slack cells
	 * @param loadFactor_ Requested fraction in (0,1), clamped to [0.25,0.90] */
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

	/** Removes all values without shrinking the backing array or changing its sentinel. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(array, invalid);
		size=0;
	}

	/** Tests whether an integer is present.
	 * @param value Integer to query
	 * @return True when value is stored */
	public boolean contains(int value){
		return value==invalid ? false : findCell(value)>=0;
	}

	/**
	 * Adds the specified value to the set if not already present.
	 * Triggers resize if adding would exceed the load factor threshold.
	 * @param value The integer value to add
	 * @return true if the value was added, false if already contained
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
	 * Adds all values from the specified IntList to this set.
	 * @param list IntList containing values to add
	 * @return Number of elements actually added (excludes duplicates)
	 */
	public int addAll(IntList list){
		int added=0;
		for(int i=0; i<list.size; i++){
			boolean b=add(list.get(i));
			if(b){added++;}
		}
		return added;
	}

	/** Re-adds a value while rebuilding a probe cluster without dispatching to IntHashSetList.add().
	 * @param value Value to restore
	 * @return true if the value was absent and added */
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
	 * Removes the specified value from the set if present.
	 * Triggers rehashing of subsequent elements to maintain lookup integrity.
	 * @param value The integer value to remove
	 * @return true if the value was removed, false if not present
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

	/** Returns the number of distinct stored integers.
	 * @return Occupied cell count */
	public int size(){return size;}
	/** Returns the insertion threshold for the current allocation.
	 * @return Maximum entry count before an insertion requests growth */
	public int sizeLimit(){return sizeLimit;}

	/** Tests whether no integers are stored.
	 * @return True when size() is zero */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Formats stored integers in backing-array order.
	 * @return The same bracketed representation as toStringListView() */
	@Override
	public String toString(){
		return toStringListView();
	}

	/** Formats occupied cells as (cell index, value) pairs in backing-array order.
	 * @return New bracketed diagnostic string */
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

	/** Formats stored values in backing-array order, which is not sorted.
	 * @return New bracketed, comma-separated value string */
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

	/** Copies stored integers into a dense array in backing-table order.
	 * @return New independent array of size() integers, without empty-cell markers */
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

	/** Checks every occupied cell's probe reachability and agreement with size.
	 * Scans physical capacity; intended for explicit diagnostics rather than hot paths.
	 * @return True if the observed entry count and lookup locations are consistent */
	public boolean verify(){
		int numValues=0;
		int numFound=0;
		for(int i=0; i<array.length; i++){
			final int value=array[i];
			if(value!=invalid){
				numValues++;
				// findCell takes a value, not a backing-array index; matches LongHashSet.verify(). [map/IntHashSet#001 FIXED]
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

	/** Rehashes the contiguous probe cluster after a deletion.
	 * @param initial Cell just emptied */
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

	/** Moves one live value to the first reachable hole after a deletion, preserving size.
	 * @param cell Occupied cell index
	 * @return True if the value moved and its old cell was cleared */
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

	/** Chooses a new empty-cell sentinel and rewrites old sentinels in the backing array. */
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

	/** Finds a value by linear probing, stopping at the first empty cell.
	 * @param value Integer to query
	 * @return Matching cell index, or -1 when absent */
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

	/** Finds the value or first empty cell using the same probe sequence as lookup.
	 * @param value Integer unequal to the current empty-cell sentinel
	 * @return Matching or empty cell index
	 * @throws RuntimeException If no matching or empty cell exists */
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

	/** Doubles the backing table around the next prime capacity. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(array.length*2L+1);
	}

	/** Rebuilds the backing array, preserving values without invoking subclass add().
	 * Uses a prime modulus plus slack; caps the prime below Integer.MAX_VALUE-extra
	 * and asserts that it exceeds the old modulus.
	 * @param size2 Positive requested modulus, greater than the current entry count */
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

	/** Stored integers; empty cells contain invalid. */
	private int[] array;
	/** Number of occupied cells. */
	private int size=0;
	/** Empty-cell sentinel, always negative and regenerated if a user attempts to insert it. */
	private int invalid;
	/** Prime modulus for initial probe positions, excluding physical slack. */
	private int modulus;
	/** Entry threshold above which insertion requests growth. */
	private int sizeLimit;
	/** Requested load fraction clamped to [0.25,0.90]. */
	private final float loadFactor;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Seeded generator for replacement empty-cell markers. */
	private static final Random randy=shared.Shared.random(1);
	/** Mask clearing a value's sign bit before taking the prime modulus. */
	private static final int MASK=Integer.MAX_VALUE;
	/** Sign-bit mask used when choosing a sentinel. */
	private static final int MINMASK=Integer.MIN_VALUE;

	/** Extra physical probe cells beyond the prime modulus. */
	private static final int extra=10;

}
