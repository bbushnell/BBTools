package map;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntList;

/**
 * Primitive int-to-long map with a prime modulus and linear probing.
 * Every int key and long value is supported; the negative empty-cell marker is
 * replaced when inserted as a key. put() inserts only when absent, whereas
 * increment() changes existing values with wrapping long arithmetic.
 * Missing lookups return -1, which can also be stored; contains() distinguishes it.
 * Unsynchronized: mutations require exclusive access. Diagnostic positive-count
 * helpers impose stronger restrictions than the general storage operations.
 *
 * @author Brian Bushnell
 * @date December 12, 2017
 */
public final class IntLongHashMap{

	/** Runs legacy assertion checks and a large fixed-size counting stress benchmark.
	 * @param args Unused */
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		IntLongHashMap map=new IntLongHashMap(20, 0.7f);
		HashMap<Integer, Long> map2=new HashMap<Integer, Long>(20, 0.7f);
		ArrayList<Integer> list=new ArrayList<Integer>();
		ArrayList<Integer> list2=new ArrayList<Integer>();
//		ArrayList<Integer> vals=new ArrayList<Integer>();
		for(int i=0; i<1000; i++){
			assert(!map.contains(i));
			assert(!map2.containsKey(i));
			list.add(Integer.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			int r=randy2.nextInt();
			list2.add(r);
		}

		for(int x : list){
			map.put(x, 2L*x);
			map2.put(x, (2L*x));
			assert(map.get(x)==(2L*x));
			assert(map2.get(x)==(2L*x));
		}

		for(int x : list){
			assert(map.get(x)==(2L*x));
			assert(map2.get(x)==(2L*x));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());

		for(int x : list2){
			map.put(x, (2L*x));
			map2.put(x, (2L*x));
			assert(map.get(x)==((2L*x)));
			assert(map2.get(x)==((2L*x)));
		}

		for(int x : list2){
			assert(map.get(x)==((2L*x)));
			assert(map2.get(x)==((2L*x)));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());

		int count=4000000;
		int runs=32;
		IntList ll=new IntList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextInt());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("LongHashMap:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map.add(x);
//				}
				final int[] y=ll.array;
				for(int z=0; z<count; z++){
					final int key=y[z];
					map.add(key);
					map.contains(key);
					map.remove(key);
					map.add(key);
				}
//				for(long x : ll.array){
//					map.remove(x);
//				}
//				map.clear();
//				assert(map.isEmpty());
//				System.err.println("Finished run "+i);
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();

//			System.err.println("HashMap:");
//			t.start();
//			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map2.add(x);
//				}
//				for(long x : ll.array){
//					map2.remove(x);
//				}
//				assert(map2.isEmpty());
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

	/** Creates an empty map with default capacity. */
	public IntLongHashMap(){
		this(256);
	}

	/** Creates an empty map with the specified initial capacity.
	 * @param initialSize Positive minimum prime modulus before adding slack cells */
	public IntLongHashMap(int initialSize){
		this(initialSize, 0.7f);
	}

	/** Creates an empty map with the specified capacity and load factor.
	 * @param initialSize Positive minimum prime modulus before adding slack cells
	 * @param loadFactor_ Requested fraction in (0,1), clamped to [0.25,0.90] */
	public IntLongHashMap(int initialSize, float loadFactor_){
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

	/** Removes all mappings and zeroes values without shrinking the arrays. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
//		assert(verify()); //123
	}

	/** Tests key occupancy independently of the stored long value.
	 * @param key Integer key to query
	 * @return True when key is present */
	public boolean contains(int key){
//		assert(verify()); //123
		return key==invalid ? false : findCell(key)>=0;
	}

	/** Alias for contains(int).
	 * @param key Integer key to query
	 * @return True when key is present */
	public boolean containsKey(int key){
		return contains(key);
	}

	/** Looks up a stored long without modifying the map.
	 * @param key Integer key to query
	 * @return Stored value or -1 when absent; contains() disambiguates a stored -1 */
	public long get(int key){
//		assert(verify()); //123
		long value=-1;
		if(key!=invalid){
			int cell=findCell(key);
			if(cell>=0){value=values[cell];}
		}
		return value;
	}

	/**
	 * Increments the value for the specified key by 1.
	 * Creates new mapping with value 1 if key doesn't exist.
	 * Long.MAX_VALUE wraps to Long.MIN_VALUE; this variant does not saturate.
	 * @param key The key to increment
	 * @return The new value after incrementing
	 */
	public long add(int key){
		return increment(key, 1);
	}

	/**
	 * Increments the value for the specified key by the given amount.
	 * Creates new mapping with the increment value if key doesn't exist.
	 * Automatically resizes the hash table when load threshold is exceeded.
	 * Addition uses wrapping long arithmetic, including both overflow directions.
	 * @param key The key to increment
	 * @param incr The amount to increment by
	 * @return The new value after incrementing
	 */
	public long increment(int key, long incr){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			values[cell]=incr;
			size++;
//			assert(verify()); //123
			if(size>sizeLimit){resize();}
//			assert(verify()); //123
			return incr;
		}else{
			values[cell]+=incr;
//			assert(verify()); //123
			return values[cell];
		}
	}

	/** Adds source values by key using wrapping arithmetic, inserting absent keys.
	 * A distinct source is preserved; self-merging doubles each value with wrapping.
	 * @param map Non-null source, stable during this call unless it is this map itself */
	public void incrementAll(IntLongHashMap map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				increment(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * Associates the specified value with the specified key.
	 * Does not overwrite existing mappings - returns false if key already exists.
	 * @param key The key to map
	 * @param value The value to associate with the key
	 * @return true if the key was newly added, false if key already existed
	 */
	public boolean put(int key, long value){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			values[cell]=value;
			size++;
			if(size>sizeLimit){resize();}
//			assert(verify()); //123
			return true;
		}
		assert(keys[cell]==key);
//		assert(verify()); //123
		return false;
	}

	/**
	 * Removes the mapping for the specified key from the hash map.
	 * Performs rehashing of subsequent entries to maintain probe sequence integrity.
	 * @param key The key to remove
	 * @return The previous value associated with the key, or -1 if not found
	 */
	public long remove(int key){
//		assert(verify()); //123
		if(key==invalid){return -1;}
		final int cell=findCell(key);
		if(cell<0){return -1;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		final long value=values[cell];
		values[cell]=0;
		size--;

		rehashFrom(cell);
//		assert(verify()); //123
		return value;
	}

	/** Returns occupied key count.
	 * @return Number of mappings */
	public int size(){return size;}

	/** Tests whether no mappings are stored.
	 * @return True when size() is zero */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Formats stored keys, not values, in backing-array order.
	 * @return The same bracketed representation as toStringListView() */
	@Override
	public String toString(){
		return toStringListView();
	}

	/** Formats occupied cells as (cell index, key, value) tuples.
	 * @return New bracketed diagnostic string in backing-array order */
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+"("+i+", "+keys[i]+", "+values[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}

	/** Formats keys in backing-array order, which is not sorted.
	 * @return New bracketed key string */
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+keys[i]);
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}

	/** Copies all stored keys into a dense array.
	 * @return New independent int array in backing-table order */
	public int[] toArray(){
		int[] x=KillSwitch.allocInt1D(size);
		int i=0;
		for(int key : keys){
			if(key!=invalid){
				x[i]=key;
				i++;
			}
		}
		return x;
	}

	/** Exports qualifying keys as longs, using a positive-counting helper contract.
	 * Asserts every occupied key is nonnegative and every occupied value is nonzero.
	 * The threshold must be positive to exclude empty zero-valued cells.
	 * @param thresh Positive minimum stored value to include
	 * @return Newly allocated long array of qualifying int keys in backing-table order */
	public long[] toArray(long thresh){
		// Active asserts encode the positive-counting contract for this exporter:
		// Stored keys must be nonnegative and values nonzero; put()/increment() do not enforce it.
		int len=0;
//		assert(verify());
		for(int i=0; i<values.length; i++){
			assert((values[i]==0)==(keys[i]==invalid)) : i+", "+values[i]+", "+keys[i]+", "+invalid+"\n"+toStringSetView();
			assert((keys[i]<0)==((keys[i]==invalid))) : toStringSetView();
			if(values[i]>=thresh){
				assert(keys[i]>=0) : "\nNegative key ("+keys[i]+", "+values[i]+", "+i+") for thresh "+thresh+":\n"+toStringSetView();
				len++;
			}
		}
		long[] x=KillSwitch.allocLong1D(len);
		for(int i=0, j=0; j<len; i++){
			if(values[i]>=thresh){
				x[j]=keys[i];
				assert(keys[i]>=0) : "\nNegative key ("+keys[i]+", "+values[i]+", "+i+") for thresh "+thresh+":\n"+toStringSetView();
				j++;
			}
		}
		return x;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks probe integrity and the stronger positive-counting invariant.
	 * Requires occupied values to be at least one, although ordinary put permits any long.
	 * @return True if counts, key reachability and positive-value invariants hold
	 * @throws AssertionError On a violated invariant when assertions are enabled */
	public boolean verify(){
		// repOK-style dev checker for the positive-counting invariant: present key
		// => value >= 1, empty cell => value == 0, and each key is reachable.
		if(keys==null){return true;}
		int numValues=0;
		int numFound=0;
		for(int i=0; i<keys.length; i++){
			final int key=keys[i];
			final long value=values[i];

			if(key==invalid){
				if(value!=0){
					assert(false) : i+", "+key+", "+value;
					return false;
				}
			}else{
				numValues++;
				if(value<1){
					assert(false) : i+", "+key+", "+value;
					return false;
				}
				final int cell=findCell(key);
				if(i==cell){
					numFound++;
				}else{
					assert(false) : i+", "+key+", "+value+", "+cell+"\n"+((cell>=0) ? keys[cell]+", "+values[cell]+"\n" : "");
					return false;
				}
			}
		}
		boolean pass=(numValues==numFound && numValues==size);
		assert(pass) : numValues+", "+numFound+", "+size;
		return pass;
	}

	/** Rehashes the contiguous probe cluster after a deletion.
	 * @param initial Cell just emptied */
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=keys.length;
		for(int cell=initial+1; cell<limit; cell++){
			final long x=keys[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
	}

	/** Moves a live key/value pair to the first reachable hole, preserving size.
	 * @param cell Occupied cell index
	 * @return True if the entry moved and its old cell was cleared */
	private boolean rehashCell(final int cell){
		final int key=keys[cell];
		final long value=values[cell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(key);
		if(cell==dest){return false;}
		assert(keys[dest]==invalid);
		keys[cell]=invalid;
		values[cell]=0;
		keys[dest]=key;
		values[dest]=value;
		return true;
	}

	/** Chooses a new empty-cell sentinel and rewrites old sentinels in the key array. */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				assert(values[i]==0);
				keys[i]=invalid;
			}
		}
	}

	/** Searches from the key's home cell, stopping at the first empty cell.
	 * @param key Integer key to query
	 * @return Matching cell index, or -1 when absent */
	private int findCell(final int key){
		if(key==invalid){return -1;}

		final int limit=keys.length, initial=(key&MASK)%modulus;
		for(int cell=initial; cell<limit; cell++){
			final long x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}

	/** Searches for a key or the first empty cell using the lookup probe sequence.
	 * @param key Integer key unequal to the current empty-cell marker
	 * @return Matching or empty cell index
	 * @throws RuntimeException If no matching or empty cell exists */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";

		final int limit=keys.length, initial=(key&MASK)%modulus;
		for(int cell=initial; cell<limit; cell++){
			final long x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}

	/** Doubles the backing table around the next prime capacity. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}

	/** Rebuilds the arrays with a larger prime modulus, preserving mappings.
	 * Caps the prime below Integer.MAX_VALUE-extra and asserts that it grows.
	 * Borrowed array references are replaced; iteration order may change.
	 * @param size2 Positive requested modulus, greater than the entry count */
	private final void resize(final long size2){
//		assert(verify()); //123
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2);
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra);
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;

		final int size3=(int)(newPrime+extra);
		sizeLimit=(int)(modulus*loadFactor);
		final int[] oldKeys=keys;
		final long[] oldValues=values;
		keys=KillSwitch.allocInt1D(size3);
		values=KillSwitch.allocLong1D(size3);
		Arrays.fill(keys, invalid);

//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldKeys.length; i++){
			int key=oldKeys[i];
			if(key!=invalid){
				put(key, oldValues[i]);
			}
		}
//		assert(verify()); //123
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows raw key storage; external mutation can invalidate probe invariants.
	 * Resizing replaces it; clear() mutates the same array in place.
	 * @return Map-owned keys with invalid() in empty cells */
	public int[] keys(){return keys;}

	/** Borrows raw values parallel to keys(); empty cells contain zero.
	 * Resizing replaces it; clear() mutates the same array in place.
	 * @return Map-owned value array; use keys() to determine occupancy */
	public long[] values(){return values;}

	/** Returns the current empty-cell marker, which changes if inserted as a real key.
	 * @return Negative sentinel absent from the stored keys */
	public int invalid(){return invalid;}

	/** Adds values from occupied cells using wrapping long arithmetic.
	 * @return Sum modulo 2^64, interpreted as a signed long */
	public long sum(){
		long sum=0;
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				sum+=values[i];
			}
		}
		return sum;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Key storage with invalid in empty cells. */
	private int[] keys;
	/** Values parallel to keys; empty cells contain zero. */
	private long[] values;
	/** Number of occupied cells. */
	private int size=0;
	/** Empty-cell sentinel for the keys array. */
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
	/** Mask clearing the sign bit before taking the modulus. */
	private static final int MASK=Integer.MAX_VALUE;
	/** Sign-bit mask used when choosing a marker. */
	private static final int MINMASK=Integer.MIN_VALUE;

	/** Extra physical probe cells beyond the logical modulus. */
	private static final int extra=10;

}
