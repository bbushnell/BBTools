package map;

import java.io.Serializable;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Hash map with Object keys and Object values.
 * Uses open addressing with linear probing for collision resolution.
 * Caches hash codes to avoid expensive equals() calls during probing.
 * Uses power-of-2 sizing for fast modulo via bitwise AND.
 *
 * Parallel arrays retain keys and values without allocating per-entry objects.
 *
 * Null keys cannot be stored; null lookups/removals behave as absent keys.
 * Null values are allowed, so use contains() to
 * distinguish absence from a stored null. Keys must retain stable hashes and
 * equality while stored. Not thread-safe; synchronize accesses that overlap a mutation.
 *
 * @author Isla
 * @date November 2, 2025
 *
 * @param <K> Key type - must properly implement hashCode() and equals()
 * @param <V> Value type
 */
public final class ObjectMap<K, V> implements Serializable{

	/** Serialization version for this map layout. */
	private static final long serialVersionUID=1L;

	/** Runs a legacy correctness comparison and allocation/timing benchmark on random strings.
	 * @param args Optional positional nonnegative key count and positive repetition count */
	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;
		ArrayList<String> list=randomStrings(size, 1, 50);
		test(list);
		bench(list, repeats);
	}

	/** Generates reproducible uppercase A-Z strings without ensuring uniqueness.
	 * @param size Nonnegative number of strings
	 * @param minLen Nonnegative minimum string length
	 * @param maxLen Inclusive maximum length, at least minLen; range must fit an int
	 * @return Newly allocated list and strings */
	private static ArrayList<String> randomStrings(int size, int minLen, int maxLen){
		// Generate random strings
		Shared.printMemory();
		ArrayList<String> strings=new ArrayList<String>(size);
		Random randy=shared.Shared.random(12345);
		int range=maxLen-minLen+1;
		for(int i=0; i<size; i++){
			int len=randy.nextInt(range)+minLen;
			StringBuilder sb=new StringBuilder(len);
			for(int j=0; j<len; j++){
				sb.append((char)('A'+randy.nextInt(26)));
			}
			strings.add(sb.toString());
		}
		Shared.printMemory();
		System.err.println("Generated "+size+" random strings");
		return strings;
	}

	/**
	 * Compares insertion/lookup with HashMap using each list element as both key and value.
	 * Reports mismatches to stderr and exits with status one; the size check also
	 * requires assertions. The supplied list is not modified.
	 * @param list Non-null list of non-null keys with stable equality and hash codes
	 */
	private static void test(ArrayList<? extends Object> list){
		System.err.println("\n*** Testing Correctness ***");

		// Build HashMap
		HashMap<Object, Object> hashMap=new HashMap<Object, Object>();
		for(int i=0; i<list.size(); i++){
			hashMap.put(list.get(i), list.get(i));
		}

		// Build ObjectObjectMap
		ObjectMap<Object, Object> objectMap=new ObjectMap<Object, Object>(Object.class, Object.class);
		for(int i=0; i<list.size(); i++){
			objectMap.put(list.get(i), list.get(i));
		}

		System.err.println("HashMap size: "+hashMap.size());
		System.err.println("ObjectObjectMap size: "+objectMap.size());
		assert(hashMap.size()==objectMap.size()) : "Size mismatch!";

		// Validate all values match
		int errors=0;
		for(int i=0; i<list.size(); i++){
			Object key=list.get(i);
			Object hashVal=hashMap.get(key);
			Object objVal=objectMap.get(key);
			if(hashVal==null || !hashVal.equals(objVal)){
				System.err.println("ERROR at index "+i+": key='"+key+"', HashMap="+hashVal+", ObjectObjectMap="+objVal);
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
	 * Benchmarks performance against {@code HashMap<Object, Object>}.
	 * Reconstructs each map per repetition, inserts/looks up/removes the supplied keys,
	 * prints timings and requests GC. This is an ad hoc benchmark, not a throughput oracle.
	 * @param list Non-null list of non-null stable keys; not modified
	 * @param repeats Positive number of cycles per map implementation
	 */
	private static void bench(ArrayList<? extends Object> list, int repeats){
		System.gc();
		Timer t=new Timer();

		{
			System.err.println("\n*** ObjectObjectMap<Object, Object> ***");
			Shared.printMemory();
			t.start();
			ObjectMap<Object, Object> map=null;
			for(int r=0; r<repeats; r++){
				map=new ObjectMap<Object, Object>(Object.class, Object.class);
				for(int i=0; i<list.size(); i++){
					map.put(list.get(i), list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					Object x=map.get(list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					map.remove(list.get(i));
				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+map.size());
			Shared.printMemory();
			map=null;
			System.gc();
		}

		{
			System.err.println("\n*** HashMap<Object, Object> ***");
			Shared.printMemory();
			t.start();
			HashMap<Object, Object> map=null;
			for(int r=0; r<repeats; r++){
				map=new HashMap<Object, Object>();
				for(int i=0; i<list.size(); i++){
					map.put(list.get(i), list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					Object x=map.get(list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					map.remove(list.get(i));
				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+map.size());
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
	 * Uses the component classes to create typed arrays returned by keys()/values().
	 * @param typeK Non-null reference component class assignable from every stored key
	 * @param typeV Non-null reference component class assignable from every stored value
	 */
	public ObjectMap(Class<K> typeK, Class<V> typeV){
		this(256, typeK, typeV);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Positive minimum logical table size, rounded up to a power of two
	 * @param typeK Component class of the keys (see {@link #ObjectMap(Class, Class)}).
	 * @param typeV Component class of the values (see {@link #ObjectMap(Class, Class)}).
	 */
	public ObjectMap(int initialSize, Class<K> typeK, Class<V> typeV){
		this(initialSize, 0.7f, typeK, typeV);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Positive minimum logical table size, rounded up to a power of two
	 * @param loadFactor_ Requested load factor in (0,1), clamped to [0.25,0.90]
	 * @param typeK Component class of the keys (see {@link #ObjectMap(Class, Class)}).
	 * @param typeV Component class of the values (see {@link #ObjectMap(Class, Class)}).
	 */
	public ObjectMap(int initialSize, float loadFactor_, Class<K> typeK, Class<V> typeV){
		assert(initialSize>0) : "Initial size must be positive";
		assert(loadFactor_>0 && loadFactor_<1) : "Load factor must be between 0 and 1";
		assert(typeK!=null && typeV!=null) : "Component types must be provided (e.g. String.class) to allocate typed backing arrays";
		this.typeK=typeK;
		this.typeV=typeV;
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Removes all entries and releases key/value references without shrinking the arrays.
	 */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, null);
		Arrays.fill(values, null);
		Arrays.fill(hashes, 0);
		size=0;
	}

	/**
	 * Gets the value associated with the given key.
	 * @param key Lookup key; null is treated as absent
	 * @return Stored value, or null when absent; contains() distinguishes a stored null
	 */
	public V get(K key){
		int cell=findCell(key);
		return cell<0 ? null : values[cell];
	}

	/**
	 * Checks if the map contains the given key.
	 * @param key Lookup key; null is treated as absent
	 * @return true if key is present
	 */
	public boolean contains(K key){
		return findCell(key)>=0;
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Non-null key compatible with the constructor's key component type
	 * @param value Value compatible with the value component type, or null
	 * @return Previous value, or null for either an absent key or a stored null
	 */
	public V put(K key, V value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * Shared keys take the source value; destination-only keys remain. Key/value objects
	 * are shared, not cloned. A distinct source is preserved; self-copy is a no-op.
	 * @param map Non-null source stable during this call, with compatible key/value types
	 */
	public void putAll(ObjectMap<K, V> map){
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
	 * @param value Value to associate with key (may be null)
	 * @return Previous value, or null for either an absent key or a stored null
	 */
	public V set(K key, V value){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		final V oldV=values[cell];
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
	 * Internal set method that takes pre-computed hash.
	 * Used during resize to avoid recomputing hashes.
	 * @param key Non-null key compatible with the backing array's component type
	 * @param value Compatible value, or null
	 * @param hash Precomputed Tools.hash32plus(key.hashCode()), not the raw hashCode()
	 * @return Previous value, or null for either an absent key or a stored null
	 */
	private V set(K key, V value, int hash){
		assert(key!=null) : "Null keys not supported";
		final int cell=findCellOrEmpty(key, hash);
		final V oldV=values[cell];
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
	 * Removes the entry with the specified key.
	 * @param key Key to remove; null is treated as absent
	 * @return Previous value, or null for either an absent key or a stored null
	 */
	public V remove(K key){
		if(key==null){return null;}
		final int cell=findCell(key);
		if(cell<0){return null;}
		assert(keys[cell].equals(key));
		final V oldV=values[cell];
		keys[cell]=null;
		values[cell]=null;
		hashes[cell]=0;
		size--;

		rehashFrom(cell);
		return oldV;
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
	 * Moves an entry to the first reachable hole after deletion, preserving its cached
	 * hash and value. Clears the vacated cell without changing the number of entries.
	 * @param cell Occupied cell index
	 * @return true if entry was moved
	 */
	private boolean rehashCell(final int cell){
		final K key=keys[cell];
		final V value=values[cell];
		final int hash=hashes[cell];
		assert(key!=null);
		final int dest=findCellOrEmpty(key, hash);
		if(cell==dest){return false;}
		assert(keys[dest]==null);
		keys[cell]=null;
		values[cell]=null;
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
	 * @param key Lookup key; null is treated as absent
	 * @return Cell index if found, -1 if not found
	 */
	private int findCell(final K key){
		if(key==null){return -1;}

		final int limit=keys.length;
		final int hash=Tools.hash32plus(key.hashCode());
		final int initial=hash&mask;

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
	 * @param key Non-null key to search for
	 * @param hash Precomputed Tools.hash32plus(key.hashCode()), not the raw hashCode()
	 * @return Cell index (either containing key or empty)
	 * @throws RuntimeException If no matching or empty cell exists
	 */
	private int findCellOrEmpty(final K key, final int hash){
		assert(key!=null) : "Null keys not supported";

		//Cached hashes avoid equals() calls on mismatches. Tools.hash32plus caps the initial
		//index below SAFE_ARRAY_LEN even at the terminal, physically capped table size.
		final int limit=keys.length;
		final int initial=hash&mask;

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
		//Double logical capacity, excluding slack; doubling the full array would round
		//ordinary growth up to four times the logical capacity. [map/IntLongHashMap2#002]
		resize((keys.length-extra)*2L);
	}

	/**
	 * Rebuilds all three arrays, preserving entries and cached hashes.
	 * Rounds the logical request up to a power of two, adds slack and caps physical
	 * length at Shared.SAFE_ARRAY_LEN. The final capped allocation uses load factor 0.85.
	 * Borrowed array identities are replaced; referenced key/value objects are retained.
	 * @param size2 Positive requested logical size, strictly greater than the entry count
	 * @throws RuntimeException If the capped allocation cannot grow the existing table
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
		@SuppressWarnings("unchecked")
		final V[] oldV=values;
		final int[] oldH=hashes;
		//Option A [map/ObjectMap#002]: allocate REAL K[]/V[] via reflection from the stored component types, so the
		//backing arrays' runtime types are genuinely K[]/V[] (e.g. String[], byte[][]) and keys()/values() return usable
		//typed arrays instead of an Object[] that CCEs at the call site. typeK/typeV are cheap cached class literals.
		@SuppressWarnings("unchecked")
		K[] tempK=(K[])Array.newInstance(typeK, (int)size3);
		@SuppressWarnings("unchecked")
		V[] tempV=(V[])Array.newInstance(typeV, (int)size3);
		keys=tempK;
		values=tempV;
		hashes=KillSwitch.allocInt1D((int)size3);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final K k=oldK[i];
			if(k!=null){
				final V v=oldV[i];
				final int h=oldH[i];
				set(k, v, h);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows the live typed key array; null entries mark empty cells.
	 * Do not modify its entries or retain it across resizing. Stored keys must retain
	 * stable hashes/equality even when callers never write directly to this array.
	 * @return Map-owned array with the constructor-supplied key component type
	 */
	public K[] keys(){return keys;}

	/**
	 * Borrows the live typed value array. Null can mean an empty slot or a stored null;
	 * test the parallel key for occupancy. Resizing replaces the backing array.
	 * @return Map-owned array with the constructor-supplied value component type
	 */
	public V[] values(){return values;}

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
	/** Array of values (parallel to keys array, may contain null values) */
	private V[] values;
	/** Array of cached hash codes (parallel to keys array, 0 for empty cells) */
	private int[] hashes;
	/** Number of entries in the map */
	private int size=0;
	/** Bit mask for fast modulo (always power of 2 minus 1) */
	private int mask;
	/** Size threshold for triggering resize (capacity * loadFactor) */
	private int sizeLimit;
	/** Requested load clamped to [0.25,0.90]; terminal allocation uses 0.85 instead. */
	private final float loadFactor;
	/** Component class of K, for allocating a real K[] backing via reflection (Array.newInstance). [map/ObjectMap#002] */
	private final Class<K> typeK;
	/** Component class of V, for allocating a real V[] backing via reflection (Array.newInstance). [map/ObjectMap#002] */
	private final Class<V> typeV;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Extra physical probe cells beyond the power-of-two logical capacity. */
	private static final int extra=10;

}
