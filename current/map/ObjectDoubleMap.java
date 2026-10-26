package map;

import java.io.Serializable;
import java.lang.reflect.Array;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Hash map with Object keys and primitive double values.
 * Uses open addressing with linear probing for collision resolution.
 * Caches hash codes to avoid expensive equals() calls during probing.
 * Uses power-of-2 sizing for fast modulo via bitwise AND.
 *
 * Parallel arrays avoid boxed values and per-entry objects; keys themselves are retained.
 *
 * Null keys cannot be stored; null lookups and removals behave as absent keys.
 * Stored keys must retain stable equality and hash codes. All double values may
 * be stored, including -1 and NaN; arithmetic follows Java floating-point rules.
 * Not thread-safe. External synchronization is required when mutations overlap
 * another access; safely published, immutable instances may be read concurrently.
 *
 * @author Brian
 * @date Feb 14, 2026
 *
 * @param <K> Key type - must properly implement hashCode() and equals()
 */
public final class ObjectDoubleMap<K> implements Serializable{

	/** Serialization version for this map layout. */
	private static final long serialVersionUID=1L;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a new map with default initial capacity (256) and load factor (0.7).
	 * @param type Non-null reference component class assignable from every stored key;
	 * used to allocate correctly typed arrays returned by keys() [map/ObjectDoubleMap#002]
	 */
	public ObjectDoubleMap(Class<K> type){
		this(256, type);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Positive minimum logical table size, rounded up to a power of two
	 * @param type Component class of the keys (see {@link #ObjectDoubleMap(Class)}).
	 */
	public ObjectDoubleMap(int initialSize, Class<K> type){
		this(initialSize, 0.7f, type);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Positive minimum logical table size, rounded up to a power of two
	 * @param loadFactor_ Requested load factor in (0,1), clamped to [0.25,0.90]
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
	 * Removes all entries and releases key references without shrinking the arrays.
	 * Empty values become zero, unlike newly allocated arrays whose empty values are -1.
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
	 * @param key Lookup key; null is treated as absent
	 * @return Stored value, or -1 when absent; contains() distinguishes a stored -1
	 */
	public double get(K key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
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
	 * @param key Non-null key compatible with the constructor's component type
	 * @param value Value to associate with key
	 * @return Previous slot value: for absent keys this is -1 after allocation or
	 * resizing, and 0 in cleared slots. Use contains() before insertion when
	 * presence matters; the historical slot-dependent return is preserved.
	 */
	public double put(K key, double value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * Shared keys take the source value; destination-only keys remain. Key references
	 * are shared, not cloned. A distinct source is preserved; self-copy is a no-op.
	 * @param map Non-null source stable during this operation, with keys compatible with this map's type
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
	 * @return Prior slot value, with the history-dependent absent return documented by put()
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
	 * @param hash Precomputed Tools.hash32plus(key.hashCode())
	 * @return Prior slot value, with the history-dependent absent return documented by put()
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
	 * An absent key starts from positive zero. Uses ordinary double addition, including
	 * infinity, NaN and signed-zero behavior; no overflow clamping is performed.
	 * @param key Non-null key compatible with the constructor's component type
	 * @param incr Signed amount to add, including non-finite values
	 * @return New value after increment
	 */
	public double increment(K key, double incr){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		final double oldV=(keys[cell]==null ? 0 : values[cell]);
		final double value=oldV+incr;
		values[cell]=value;
		if(keys[cell]==null){
			keys[cell]=key;
			hashes[cell]=hash;
			size++;
			if(size>sizeLimit){resize();}
		}
		return value;//A resize can move the key away from the original cell.
	}

	/**
	 * Increments all entries in this map by corresponding values from another map.
	 * Absent keys start from positive zero. Key references are shared, not cloned.
	 * A distinct source is preserved; a self-merge doubles each value using double addition.
	 * @param map Non-null source, stable during this call unless it is this map itself
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
	 * Keys absent from this map copy the source value exactly, including negative values.
	 * Shared keys use Math.max semantics: NaN propagates and positive zero wins over
	 * negative zero. Destination-only keys remain. Key references are shared and a
	 * distinct source is preserved; self-merging leaves values unchanged.
	 * @param map Non-null source stable during this operation, with compatible key types
	 */
	public void setToMax(ObjectDoubleMap<K> map){
		for(int i=0; i<map.keys.length; i++){
			final K key=map.keys[i];
			if(key!=null){
				final int hash=Tools.hash32plus(key.hashCode());
				final int cell=findCellOrEmpty(key, hash);
				final double value=map.values[i];
				if(keys[cell]==null){
					keys[cell]=key;
					hashes[cell]=hash;
					values[cell]=value;
					size++;
					if(size>sizeLimit){resize();}
				}else{values[cell]=Tools.max(value, values[cell]);}
			}
		}
	}

	/**
	 * Removes the entry with the specified key.
	 * @param key Key to remove; null is treated as absent
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
	 * Moves an entry to the first reachable hole after deletion, preserving its cached
	 * hash and value. Clears the vacated cell without changing the number of entries.
	 * @param cell Occupied cell index
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
	 * New empty values are -1; borrowed array references are replaced.
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

	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	/** Formats occupied entries through toBytes(); key text is not escaped.
	 * @return Brace-delimited key=value entries in backing-array order */
	@Override
	public String toString(){
		return toBytes().toString();
	}

	/** Appends occupied entries with ByteBuilder's three-fractional-digit formatting.
	 * Key text is not escaped; the result is diagnostic text, not a round-trip format.
	 * @return Newly allocated builder containing brace-delimited key=value entries */
	public ByteBuilder toBytes(){
		ByteBuilder bb=new ByteBuilder();
		bb.append('{');
		for(int i=0; i<keys.length; i++){
			K key=keys[i];
			if(key!=null){
				if(bb.length>1){bb.comma().space();}
				bb.append(key.toString()).equals().append(values[i], 3);
			}
		}
		return bb.append('}');
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows the live typed key array; null entries mark empty cells.
	 * Do not modify its entries or retain it across resizing. Mutating a stored key's
	 * equality/hash state invalidates the cached hashes even if the array is untouched.
	 * @return Map-owned array with the constructor-supplied component type
	 */
	public K[] keys(){return keys;}

	/**
	 * Returns the internal value array.
	 * Empty slots can contain -1 or 0; test the parallel key for occupancy.
	 * @return Borrowed value array, which may be replaced during resizing
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
	/** Requested load clamped to [0.25,0.90]; terminal allocation uses 0.85 instead. */
	private final float loadFactor;
	/** Component class of K, for allocating a real K[] backing via reflection (Array.newInstance). [map/ObjectDoubleMap#002] */
	private final Class<K> type;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Extra physical probe cells beyond the power-of-two logical capacity. */
	private static final int extra=10;

}
