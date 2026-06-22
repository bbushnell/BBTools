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
 * More memory-efficient than HashMap<K, V> by:
 * - Using arrays instead of Entry objects
 * - Avoiding per-entry object overhead
 * - Better cache locality with linear probing
 * 
 * Thread-safety: Not thread-safe. External synchronization required for concurrent access.
 * 
 * @author Isla
 * @date November 2, 2025
 * 
 * @param <K> Key type - must properly implement hashCode() and equals()
 * @param <V> Value type
 */
public final class ObjectMap<K, V> implements Serializable {

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
	 * Tests correctness by comparing ObjectObjectMap against HashMap.
	 * Creates random list, inserts them in both maps, and validates all values match.
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
	 * Benchmarks performance against HashMap<Object, Object>.
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
	 * @param typeK Component class of the KEYS (e.g. String.class, byte[].class) — for a real K[] backing.
	 * @param typeV Component class of the VALUES — for a real V[] backing. Both let keys()/values() return genuine
	 *        typed arrays instead of an Object[] that CCEs at the call site. Pass class literals; they are cheap
	 *        cached constants. [map/ObjectMap#002]
	 */
	public ObjectMap(Class<K> typeK, Class<V> typeV){
		this(256, typeK, typeV);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param typeK Component class of the keys (see {@link #ObjectMap(Class, Class)}).
	 * @param typeV Component class of the values (see {@link #ObjectMap(Class, Class)}).
	 */
	public ObjectMap(int initialSize, Class<K> typeK, Class<V> typeV){
		this(initialSize, 0.7f, typeK, typeV);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param loadFactor Load factor (0.25-0.90) - map resizes when size exceeds capacity*loadFactor
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
	 * Removes all entries from the map.
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
	 * @param key Key to look up
	 * @return Value associated with key, or null if key not found
	 */
	public V get(K key){
		int cell=findCell(key);
		return cell<0 ? null : values[cell];
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
	 * @return Previous value associated with key, or null if key was not present
	 */
	public V put(K key, V value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * @param map Source map to copy from
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
	 * @return Previous value associated with key, or null if key was not present
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
	 * @param key Key to remove
	 * @return Previous value associated with key, or null if key was not found
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
	 * Attempts to move an entry to its ideal position.
	 * @param cell Position of entry to rehash
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

		//CLEVER [verified]: hashes[] caches each key's hash so the probe compares cheap ints first
		//(hashes[cell]==hash) and only calls equals() on a hash match. initial=hash & mask is in range for ANY mask:
		//hash=Tools.hash32plus(...) is CAPPED in [0, 0x7FFFF800) < SAFE_ARRAY_LEN, so even at the top tier
		//(mask=Integer.MAX_VALUE after the 2x-growth fix) hash & mask == hash is a valid index (see hash32plus javadoc).
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
		//grow 2x: double the LOGICAL pow2 capacity (keys.length-extra==pow2), not keys.length (=pow2+extra).
		//keys.length*2L overshot the power-of-2 boundary -> rounded up to 4*pow2 (same family bug as the pow2 maps,
		//anchor [map/IntLongHashMap2#002]); fixed to true 2x growth. Output unchanged; ~2x less memory.
		resize((keys.length-extra)*2L);
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
	 * Returns the internal key array (the live backing array; null entries mark empty cells).
	 * This is a GENUINE K[] (e.g. String[], byte[][]) because the backing is allocated via Array.newInstance(typeK,...)
	 * from the key component class passed to the constructor [map/ObjectMap#002] — usable directly (assigned to K[],
	 * iterated, .length) without the ClassCastException a plain Object[] backing would cause at the call site.
	 * WARNING: live array, contains nulls for empty cells. Use with caution.
	 * @return Internal key array (a real K[])
	 */
	public K[] keys(){return keys;}

	/**
	 * Returns the internal value array (the live backing array; null entries mark empty cells).
	 * This is a GENUINE V[] (e.g. String[], byte[][]) because the backing is allocated via Array.newInstance(typeV,...)
	 * from the value component class passed to the constructor [map/ObjectMap#002] — usable directly (assigned to V[],
	 * iterated, .length) without the ClassCastException a plain Object[] backing would cause at the call site.
	 * WARNING: live array, contains nulls for empty cells. Use with caution.
	 * @return Internal value array (a real V[])
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
	/** Load factor (fraction of capacity before resize) */
	private final float loadFactor;
	/** Component class of K, for allocating a real K[] backing via reflection (Array.newInstance). [map/ObjectMap#002] */
	private final Class<K> typeK;
	/** Component class of V, for allocating a real V[] backing via reflection (Array.newInstance). [map/ObjectMap#002] */
	private final Class<V> typeV;

	/** Extra space beyond power-of-2 size to reduce wrap-around collisions */
	private static final int extra=10;

}