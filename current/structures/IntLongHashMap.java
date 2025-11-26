package structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date December 12, 2017
 *
 */
public final class IntLongHashMap{
	
	/**
	 * Test harness comparing performance against Java HashMap.
	 * Benchmarks insertion, lookup, and removal operations on large datasets.
	 * @param args Command line arguments (unused)
	 */
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
	
	/**
	 * Constructs a hash map with default initial capacity of 256 and load factor 0.7.
	 */
	public IntLongHashMap(){
		this(256);
	}
	
	/** Constructs a hash map with specified initial capacity and default load factor 0.7.
	 * @param initialSize Initial capacity of the hash table */
	public IntLongHashMap(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Constructs a hash map with specified initial capacity and load factor.
	 * Load factor is clamped to range [0.25, 0.90] for optimal performance.
	 * @param initialSize Initial capacity of the hash table
	 * @param loadFactor_ Target load factor for resize threshold
	 */
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
	
	/** Removes all key-value mappings from the hash map.
	 * Resets all cells to invalid state and clears values. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
//		assert(verify()); //123
	}
	
	/**
	 * Tests if the specified key is present in the hash map.
	 * @param key The key to test for presence
	 * @return true if the key exists in the map, false otherwise
	 */
	public boolean contains(int key){
//		assert(verify()); //123
		return key==invalid ? false : findCell(key)>=0;
	}
	
	/**
	 * Tests if the specified key is present in the hash map.
	 * Alias for contains() method for HashMap compatibility.
	 * @param key The key to test for presence
	 * @return true if the key exists in the map, false otherwise
	 */
	public boolean containsKey(int key){
		return contains(key);
	}
	
	/**
	 * Retrieves the value associated with the specified key.
	 * @param key The key to look up
	 * @return The associated value, or -1 if key not found
	 */
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
	 * Increment this key's value by 1.
	 * @param key
	 * @return New value
	 */
	public long add(int key){
		return increment(key, 1);
	}
	
	/**
	 * Increment this key's value by incr.
	 * @param key
	 * @param incr
	 * @return New value
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
	
	/**
	 * Adds all key-value pairs from another map to this map.
	 * For each key in the source map, increments this map's value by the source value.
	 * @param map The source map to add from
	 */
	public void incrementAll(IntLongHashMap map) {
		for(int i=0; i<map.keys.length; i++) {
			if(map.keys[i]!=map.invalid) {
				increment(map.keys[i], map.values[i]);
			}
		}
	}
	
	/**
	 * Map this key to value.
	 * @param key
	 * @param value
	 * @return true if the key was added, false if it was already contained.
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
	 * Remove this key from the map.
	 * @param key
	 * @return Old value.
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
	
	/** Returns the number of key-value mappings in the hash map */
	public int size(){return size;}
	
	/** Returns true if the hash map contains no key-value mappings */
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns a string representation showing internal structure.
	 * Displays cell index, key, and value for each occupied cell.
	 * @return String showing internal hash table state
	 */
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
	
	/**
	 * Returns a string representation as a list of keys.
	 * Shows only the keys present in the map, not internal structure.
	 * @return Comma-separated list of keys in brackets
	 */
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
	
	/**
	 * Returns an array containing all keys in the hash map.
	 * The order of keys in the array is not guaranteed.
	 * @return Array of all keys present in the map
	 */
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
	
	/**
	 * Returns an array of keys whose values meet or exceed the threshold.
	 * Only keys with associated values >= thresh are included.
	 * @param thresh Minimum value threshold for inclusion
	 * @return Array of keys with values at or above threshold
	 */
	public long[] toArray(long thresh){
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
	
	/**
	 * Verifies internal consistency of the hash map data structure.
	 * Checks that all keys can be found at their correct hash positions
	 * and that size tracking is accurate. Used for debugging and testing.
	 * @return true if structure is consistent, false if corrupted
	 */
	public boolean verify(){
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
	
	/**
	 * Rehashes all entries after a removal to maintain probe sequence integrity.
	 * Starting from the removal position, moves entries backward to fill gaps
	 * created by the removed entry. Essential for open addressing correctness.
	 * @param initial Position where entry was removed
	 */
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
	
	/**
	 * Rehashes a single cell to its optimal position.
	 * Moves the key-value pair to the earliest available cell in its probe sequence.
	 * @param cell The cell to rehash
	 * @return true if the entry was moved, false if already in optimal position
	 */
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
	
	/**
	 * Generates a new invalid key value when the current one conflicts.
	 * Updates all empty cells to use the new invalid marker.
	 * Required when a key matching the invalid value is inserted.
	 */
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
	
	/**
	 * Finds the cell containing the specified key using linear probing.
	 * Searches from the hash position forward, wrapping around if necessary.
	 * @param key The key to locate
	 * @return Cell index containing the key, or -1 if not found
	 */
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
	
	/**
	 * Finds the cell containing the key or the first empty cell in its probe sequence.
	 * Used for insertions to locate where a key should be placed.
	 * @param key The key to locate or place
	 * @return Cell index for the key or first available empty cell
	 */
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
	
	/** Resizes the hash table when load factor is exceeded.
	 * Doubles capacity plus 1 and rehashes all existing entries. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
	/**
	 * Resizes hash table to specified capacity using prime number sizing.
	 * Finds the smallest prime >= size2 and rehashes all entries to new table.
	 * Handles integer overflow by capping at maximum safe prime.
	 * @param size2 Target minimum capacity
	 */
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

	/** Returns direct reference to the internal keys array */
	public int[] keys() {return keys;}

	/** Returns direct reference to the internal values array */
	public long[] values() {return values;}

	/** Returns the current invalid key marker value */
	public int invalid() {return invalid;}

	/**
	 * Computes the sum of all values in the hash map.
	 * Iterates through all entries and accumulates their values.
	 * @return Sum of all values currently stored in the map
	 */
	public long sum() {
		long sum=0;
		for(int i=0; i<keys.length; i++) {
			if(keys[i]!=invalid) {
				sum+=values[i];
			}
		}
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Hash table array storing keys, with invalid values marking empty cells */
	private int[] keys;
	/** Parallel array storing values corresponding to keys array */
	private long[] values;
	/** Current number of key-value pairs stored in the hash map */
	private int size=0;
	/** Value for empty cells */
	private int invalid;
	/** Prime number used as modulus for hash function */
	private int modulus;
	/** Maximum size before triggering resize, based on load factor */
	private int sizeLimit;
	/** Target ratio of occupied cells to total capacity */
	private final float loadFactor;
	
	/** Random number generator for generating invalid key values */
	private static final Random randy=new Random(1);
	/** Bit mask for ensuring non-negative hash values (Integer.MAX_VALUE) */
	private static final int MASK=Integer.MAX_VALUE;
	/** Bit mask ensuring invalid keys are negative (Integer.MIN_VALUE) */
	private static final int MINMASK=Integer.MIN_VALUE;
	
	/** Additional capacity beyond prime modulus to reduce collisions */
	private static final int extra=10;
	
}
