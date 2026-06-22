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
import structures.ByteBuilder;

/**
 * Hash map with Object keys and primitive double values.
 * Uses open addressing with linear probing for collision resolution.
 * Caches hash codes to avoid expensive equals() calls during probing.
 * Uses power-of-2 sizing for fast modulo via bitwise AND.
 * 
 * Significantly more memory-efficient than HashMap<K, Double> by:
 * - Storing primitive double values instead of Double objects
 * - Using arrays instead of Entry objects
 * - Avoiding per-entry object overhead
 * 
 * Thread-safety: Not thread-safe. External synchronization required for concurrent access.
 * 
 * @author Brian
 * @date Feb 14, 2026
 * 
 * @param <K> Key type - must properly implement hashCode() and equals()
 */
public final class ObjectDoubleMap<K> implements Serializable {

	private static final long serialVersionUID = 1L;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a new map with default initial capacity (256) and load factor (0.7).
	 * @param type Component class of the keys (e.g. String.class, byte[].class) — used to allocate a real K[] backing
	 *        array via reflection so keys() returns a genuine typed array (a plain Object[] backing CCEs at the call
	 *        site under a concrete K). Pass a class literal; it is a cheap cached constant. [map/ObjectDoubleMap#002]
	 */
	public ObjectDoubleMap(Class<K> type){
		this(256, type);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param type Component class of the keys (see {@link #ObjectDoubleMap(Class)}).
	 */
	public ObjectDoubleMap(int initialSize, Class<K> type){
		this(initialSize, 0.7f, type);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param loadFactor Load factor (0.25-0.90) - map resizes when size exceeds capacity*loadFactor
	 * @param type Component class of the keys (see {@link #ObjectDoubleMap(Class)}).
	 */
	public ObjectDoubleMap(int initialSize, float loadFactor_, Class<K> type){
		assert(initialSize>0) : "Initial size must be positive";
		assert(loadFactor_>0 && loadFactor_<1) : "Load factor must be between 0 and 1";
		assert(type!=null) : "Component type must be provided (e.g. String.class) to allocate a typed backing array";
		this.type=type;
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
	public double get(K key){
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
	public double put(K key, double value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * @param map Source map to copy from
	 */
	public void putAll(ObjectDoubleMap<K> map){
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
	public double set(K key, double value){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		final double oldV=values[cell];
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
	private double set(final K key, final double value, final int hash){
		assert(key!=null) : "Null keys not supported";
		final int cell=findCellOrEmpty(key, hash);
		final double oldV=values[cell];
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
	public double increment(K key){
		return increment(key, 1);
	}

	/**
	 * Increments the value associated with the key by the specified amount.
	 * If key is not present, inserts it with the increment value.
	 * @param key Key to increment
	 * @param incr Amount to add (can be negative)
	 * @return New value after increment
	 */
	public double increment(K key, double incr){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		final double oldV=values[cell];
		final double value=oldV+incr;
		values[cell]=value;
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
	public void incrementAll(ObjectDoubleMap<K> map){
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
	public void setToMax(ObjectDoubleMap<K> map){
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
		final double value=values[cell];
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
		final double[] oldV=values;
		final int[] oldH=hashes;
		//Option A [map/ObjectDoubleMap#002]: allocate a REAL K[] via reflection from the stored component type, so the
		//backing array's runtime type is genuinely K[] (e.g. String[], byte[][]) and keys() returns a usable typed
		//array instead of an Object[] that CCEs at the call site. `type` is a cheap cached class literal.
		@SuppressWarnings("unchecked")
		K[] temp=(K[])Array.newInstance(type, (int)size3);
		keys=temp;
		values=KillSwitch.allocDouble1D((int)size3);
		Arrays.fill(values, -1);
		hashes=KillSwitch.allocInt1D((int)size3);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final K k=oldK[i];
			if(k!=null){
				final double v=oldV[i];
				final int h=oldH[i];
				set(k, v, h);
			}
		}
	}
	
	@Override
	public String toString(){
		return toBytes().toString();
	}
	
	public ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		bb.append('{');
		for(int i=0; i<keys.length; i++) {
			K key=keys[i];
			if(key!=null) {
				if(bb.length>1) {bb.comma().space();}
				bb.append(key.toString()).equals().append(values[i], 3);
			}
		}
		return bb;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the internal key array (the live backing array; null entries mark empty cells).
	 * This is a GENUINE K[] (e.g. String[], byte[][]) because the backing is allocated via Array.newInstance(type,...)
	 * from the component class passed to the constructor [map/ObjectDoubleMap#002] — so it can be used directly (assigned
	 * to K[], iterated, .length) without the ClassCastException a plain Object[] backing would cause at the call site.
	 * WARNING: live array, contains nulls for empty cells. Use with caution.
	 * @return Internal key array (a real K[])
	 */
	public K[] keys(){return keys;}

	/**
	 * Returns the internal value array.
	 * WARNING: Contains 0 for empty cells corresponding to null keys.
	 * @return Internal value array
	 */
	public double[] values(){return values;}

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
	private double[] values;
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
	/** Component class of K, for allocating a real K[] backing via reflection (Array.newInstance). [map/ObjectDoubleMap#002] */
	private final Class<K> type;

	/** Extra space beyond power-of-2 size to reduce wrap-around collisions */
	private static final int extra=10;

}