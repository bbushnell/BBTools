package map;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Hash map with Object keys and primitive int values.
 * Uses open addressing with linear probing for collision resolution.
 * Caches hash codes to avoid expensive equals() calls during probing.
 * Uses power-of-2 sizing for fast modulo via bitwise AND.
 * 
 * Significantly more memory-efficient than HashMap<K, Integer> by:
 * - Storing primitive int values instead of Integer objects
 * - Using arrays instead of Entry objects
 * - Avoiding per-entry object overhead
 * 
 * Thread-safety: Not thread-safe. External synchronization required for concurrent access.
 * 
 * @author Isla
 * @date November 2, 2025
 * 
 * @param <K> Key type - must properly implement hashCode() and equals()
 */
public final class ObjectIntMap<K> implements Serializable {

	private static final long serialVersionUID = 1L;

	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;
		ArrayList<String> list=randomStrings(size, 1, 50);
		test(list);
		bench(list, repeats);
	}

	private static ArrayList<String> randomStrings(int size, int minLen, int maxLen){
		// Generate random strings
		Shared.printMemory();
		ArrayList<String> strings=new ArrayList<String>(size);
		Random randy=new Random(12345);
		int range=maxLen-minLen+1;
		for(int i=0; i<size; i++){
			int len=randy.nextInt(range)+minLen;
			StringBuilder sb=new StringBuilder(len);
			for(int j=0; j<len; j++){
				sb.append((char)('0'+randy.nextInt(76)));
			}
			strings.add(sb.toString());
		}
		Shared.printMemory();
		System.err.println("Generated "+size+" random strings");
		return strings;
	}

	/**
	 * Tests correctness by comparing ObjectIntMap against HashMap.
	 * Creates random list, inserts them in both maps, and validates all values match.
	 */
	private static void test(ArrayList<? extends Object> list){
		System.err.println("\n*** Testing Correctness ***");

		// Build HashMap
		HashMap<Object, Integer> hashMap=new HashMap<Object, Integer>();
		for(int i=0; i<list.size(); i++){
			hashMap.put(list.get(i), i);
		}

		// Build ObjectIntMap
		ObjectIntMap<Object> objectMap=new ObjectIntMap<Object>();
		for(int i=0; i<list.size(); i++){
			objectMap.put(list.get(i), i);
		}

		System.err.println("HashMap size: "+hashMap.size());
		System.err.println("ObjectIntMap size: "+objectMap.size());
		assert(hashMap.size()==objectMap.size()) : "Size mismatch!";

		// Validate all values match
		int errors=0;
		for(int i=0; i<list.size(); i++){
			Object key=list.get(i);
			Integer hashVal=hashMap.get(key);
			int objVal=objectMap.get(key);
			if(hashVal==null || hashVal!=objVal){
				System.err.println("ERROR at index "+i+": key='"+key+"', HashMap="+hashVal+", ObjectIntMap="+objVal);
				errors++;
				if(errors>10){break;} // Don't spam too much
			}
		}

		if(errors==0){
			System.err.println("*** PASS: All values match! ***");
		}else{
			System.err.println("*** FAIL: "+errors+" mismatches found! ***");
			System.exit(1);
		}
	}

	/**
	 * Benchmarks performance against HashMap<Object, Integer>.
	 */
	private static void bench(ArrayList<? extends Object> list, int repeats){
		System.gc();
		Timer t=new Timer();

		{
			System.err.println("\n*** ObjectIntMap<Object> ***");
			Shared.printMemory();
			t.start();
			long sum=0;
			ObjectIntMap<Object> map=null;
			for(int r=0; r<repeats; r++){
				map=new ObjectIntMap<Object>();
				for(int i=0; i<list.size(); i++){
					map.put(list.get(i), i);
				}
				for(int i=0; i<list.size(); i++){
					sum+=map.get(list.get(i));
				}
				//				for(int i=0; i<list.size(); i++){
				//					map.remove(list.get(i));
				//				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+map.size()+", sum="+sum);
			Shared.printMemory();
			map=null;
			System.gc();
		}

		{
			System.err.println("\n*** HashMap<Object, Integer> ***");
			Shared.printMemory();
			t.start();
			long sum=0;
			HashMap<Object, Integer> map=null;
			for(int r=0; r<repeats; r++){
				map=new HashMap<Object, Integer>();
				for(int i=0; i<list.size(); i++){
					map.put(list.get(i), i);
				}
				for(int i=0; i<list.size(); i++){
					sum+=map.get(list.get(i));
				}
				//				for(int i=0; i<list.size(); i++){
				//					map.remove(list.get(i));
				//				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+map.size()+", sum="+sum);
			Shared.printMemory();
			map=null;
			System.gc();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a new map with default initial capacity (256) and load factor (0.7).
	 */
	public ObjectIntMap(){
		this(256);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 */
	public ObjectIntMap(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param loadFactor Load factor (0.25-0.90) - map resizes when size exceeds capacity*loadFactor
	 */
	public ObjectIntMap(int initialSize, float loadFactor_){
		assert(initialSize>0) : "Initial size must be positive";
		assert(loadFactor_>0 && loadFactor_<1) : "Load factor must be between 0 and 1";
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Removes all entries from the map.
	 */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, null);
		Arrays.fill(values, 0);
		Arrays.fill(hashes, 0);
		size=0;
	}

	/**
	 * Gets the value associated with the given key.
	 * @param key Key to look up
	 * @return Value associated with key, or -1 if key not found
	 */
	public int get(K key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
	}

	/**
	 * Checks if the map contains the given key.
	 * @param key Key to check
	 * @return true if key is present
	 */
	public boolean contains(K key){
		return findCell(key)>=0;
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Key to insert/update
	 * @param value Value to associate with key
	 * @return Previous value associated with key, or -1 if key was not present
	 */
	public int put(K key, int value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * @param map Source map to copy from
	 */
	public void putAll(ObjectIntMap<K> map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=null){
				put(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Key to insert/update (must not be null)
	 * @param value Value to associate with key
	 * @return Previous value associated with key, or -1 if key was not present
	 */
	public int set(K key, int value){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		final int oldV=values[cell];
		values[cell]=value;
		if(keys[cell]==null){
			keys[cell]=key;
			hashes[cell]=hash;
			size++;
			if(size>sizeLimit){resize();}
		}
		return oldV;
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Key to insert/update (must not be null)
	 * @param value Value to associate with key
	 * @param hash Hashcode of key
	 * @return Previous value associated with key, or -1 if key was not present
	 */
	private int set(final K key, final int value, final int hash){
		assert(key!=null) : "Null keys not supported";
		final int cell=findCellOrEmpty(key, hash);
		final int oldV=values[cell];
		values[cell]=value;
		if(keys[cell]==null){
			keys[cell]=key;
			hashes[cell]=hash;
			size++;
			if(size>sizeLimit){resize();}
		}
		return oldV;
	}

	/**
	 * Increments the value associated with the key by 1.
	 * If key is not present, inserts it with value 1.
	 * @param key Key to increment
	 * @return New value after increment
	 */
	public int increment(K key){
		return increment(key, 1);
	}

	/**
	 * Increments the value associated with the key by the specified amount.
	 * If key is not present, inserts it with the increment value.
	 * Values are clamped to Integer.MAX_VALUE on overflow.
	 * @param key Key to increment
	 * @param incr Amount to add (can be negative)
	 * @return New value after increment
	 */
	public int increment(K key, int incr){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		final int oldV=values[cell];
		final int value=oldV+incr;
		values[cell]=Tools.min(Integer.MAX_VALUE, value);
		if(keys[cell]==null){
			keys[cell]=key;
			hashes[cell]=hash;
			size++;
			if(size>sizeLimit){resize();}
		}
		return values[cell];
	}

	/**
	 * Increments all entries in this map by corresponding values from another map.
	 * Keys not present in this map are inserted with their values from the source map.
	 * @param map Source map containing increment values
	 */
	public void incrementAll(ObjectIntMap<K> map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=null){
				increment(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * For each key in the source map, sets this map's value to the maximum
	 * of the current value and the source value.
	 * @param map Source map to compare against
	 */
	public void setToMax(ObjectIntMap<K> map){
		for(int i=0; i<map.keys.length; i++){
			final K key=map.keys[i];
			if(key!=null){
				put(key, Tools.max(map.values[i], get(key)));
			}
		}
	}

	/**
	 * Removes the entry with the specified key.
	 * @param key Key to remove
	 * @return true if key was present and removed, false if key was not found
	 */
	public boolean remove(K key){
		if(key==null){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(keys[cell].equals(key));
		keys[cell]=null;
		values[cell]=0;
		hashes[cell]=0;
		size--;

		rehashFrom(cell);
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Rehashes entries after a removal to maintain probe sequence integrity.
	 * @param initial Starting position of removed entry
	 */
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=keys.length;
		for(int cell=initial+1; cell<limit; cell++){
			final K key=keys[cell];
			if(key==null){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final K key=keys[cell];
			if(key==null){return;}
			rehashCell(cell);
		}
	}

	/**
	 * Attempts to move an entry to its ideal position.
	 * @param cell Position of entry to rehash
	 * @return true if entry was moved
	 */
	private boolean rehashCell(final int cell){
		final K key=keys[cell];
		final int value=values[cell];
		final int hash=hashes[cell];
		assert(key!=null);
		final int dest=findCellOrEmpty(key, hash);
		if(cell==dest){return false;}
		assert(keys[dest]==null);
		keys[cell]=null;
		values[cell]=0;
		hashes[cell]=0;
		keys[dest]=key;
		values[dest]=value;
		hashes[dest]=hash;

		return true;
	}

	/**
	 * Finds the cell containing the given key.
	 * Uses linear probing to resolve collisions.
	 * Checks hash code equality before calling expensive equals().
	 * @param key Key to search for
	 * @return Cell index if found, -1 if not found
	 */
	private int findCell(final K key){
		if(key==null){return -1;}

		final int limit=keys.length;
		final int hash=Tools.hash32plus(key.hashCode());
		final int initial=hash & mask;

		for(int cell=initial; cell<limit; cell++){
			final K x=keys[cell];
			if(x==null){return -1;}
			if(hashes[cell]==hash && x.equals(key)){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final K x=keys[cell];
			if(x==null){return -1;}
			if(hashes[cell]==hash && x.equals(key)){return cell;}
		}
		return -1;
	}

	/**
	 * Finds the cell containing the key, or an empty cell where it can be inserted.
	 * Checks hash code equality before calling expensive equals().
	 * @param key Key to search for
	 * @param hash Pre-computed hash code of key
	 * @return Cell index (either containing key or empty)
	 */
	private int findCellOrEmpty(final K key, final int hash){
		assert(key!=null) : "Null keys not supported";

		final int limit=keys.length;
		final int initial=hash & mask;

		for(int cell=initial; cell<limit; cell++){
			final K x=keys[cell];
			if(x==null || (hashes[cell]==hash && x.equals(key))){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final K x=keys[cell];
			if(x==null || (hashes[cell]==hash && x.equals(key))){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}

	/**
	 * Doubles the map capacity and rehashes all entries.
	 */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L);
	}

	/**
	 * Resizes the map to at least the specified size and rehashes all entries.
	 * Actual size will be rounded up to the next power of 2.
	 * @param size2 Minimum new size
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;

		long size3=Long.highestOneBit(size2);
		if(size3<size2){size3<<=1;}
		mask=(int)(size3-1);
		size3=Math.min(size3+extra, Shared.SAFE_ARRAY_LEN);
		if((keys!=null && size3<=keys.length) || size3>Shared.SAFE_ARRAY_LEN){
			throw new RuntimeException("Map hit capacity at "+size);
		}
		final float loadFactor2=(size3<Shared.SAFE_ARRAY_LEN ? loadFactor : 0.85f);
		sizeLimit=(int)((size3-extra)*loadFactor2);

		@SuppressWarnings("unchecked")
		final K[] oldK=keys;
		final int[] oldV=values;
		final int[] oldH=hashes;
		@SuppressWarnings("unchecked")
		K[] temp=(K[])new Object[(int)size3];
		keys=temp;
		values=KillSwitch.allocInt1D((int)size3);
		hashes=KillSwitch.allocInt1D((int)size3);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final K k=oldK[i];
			if(k!=null){
				final int v=oldV[i];
				final int h=oldH[i];
				set(k, v, h);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the internal key array.
	 * WARNING: Contains null entries for empty cells. Use with caution.
	 * @return Internal key array
	 */
	public K[] keys(){return keys;}

	/**
	 * Returns the internal value array.
	 * WARNING: Contains 0 for empty cells corresponding to null keys.
	 * @return Internal value array
	 */
	public int[] values(){return values;}

	/**
	 * Returns the number of key-value pairs in the map.
	 * @return Number of entries
	 */
	public int size(){return size;}

	/**
	 * Checks if the map is empty.
	 * @return true if map contains no entries
	 */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Array of keys (null for empty cells) */
	private K[] keys;
	/** Array of values (parallel to keys array) */
	private int[] values;
	/** Array of cached hash codes (parallel to keys array, 0 for empty cells) */
	private int[] hashes;
	/** Number of entries in the map */
	private int size=0;
	/** Bit mask for fast modulo (always power of 2 minus 1) */
	private int mask;
	/** Size threshold for triggering resize (capacity * loadFactor) */
	private int sizeLimit;
	/** Load factor (fraction of capacity before resize) */
	private final float loadFactor;

	/** Extra space beyond power-of-2 size to reduce wrap-around collisions */
	private static final int extra=10;

}