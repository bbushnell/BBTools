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
 * Hash map with primitive long keys and Object values.
 * Uses open addressing with linear probing for collision resolution.
 * Uses power-of-2 sizing and hash mixing for fast, well-distributed hashing.
 * 
 * Ideal for kmer to data structure mappings (e.g., kmer → IntList of positions).
 * Hash mixing ensures adjacent kmers (like AAAAAA vs AAAAAC) distribute evenly.
 * 
 * More memory-efficient than HashMap<Long, V> by:
 * - Storing primitive long keys instead of Long objects
 * - Using arrays instead of Entry objects
 * - Avoiding per-entry object overhead
 * - Better cache locality with linear probing
 * 
 * Thread-safety: Not thread-safe. External synchronization required for concurrent access.
 * 
 * @author Isla
 * @date November 2, 2025
 * 
 * @param <V> Value type
 */
public final class LongObjectMap<V> implements Serializable {

	private static final long serialVersionUID = 1L;

	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;

		// Generate random keys with collisions
		System.err.println("Generating "+size+" random keys...");
		Random randy=shared.Shared.random(12345);
		long[] keys=new long[size];
		for(int i=0; i<size; i++){
			long x=randy.nextLong() & Long.MAX_VALUE;
			keys[i]=(long)Math.sqrt(x);
		}

		// Generate random string values (unrelated to keys)
		System.err.println("Generating "+size+" random values...");
		ArrayList<String> values=new ArrayList<String>(size);
		for(int i=0; i<size; i++){
			int len=randy.nextInt(50)+1;
			StringBuilder sb=new StringBuilder(len);
			for(int j=0; j<len; j++){
				sb.append((char)('A'+randy.nextInt(26)));
			}
			values.add(sb.toString());
		}
		System.err.println("Keys and values generated.");

		bench(keys, values, repeats);
	}

	private static void bench(long[] keys, ArrayList<String> values, int repeats){
		final int size=keys.length;
		System.gc();
		Timer t=new Timer();

		{
			System.err.println("\n*** LongObjectMap<String> ***");
			Shared.printMemory();
			t.start();
			LongObjectMap<String> map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new LongObjectMap<String>(String.class);
				for(int i=0; i<size; i++){
					map.put(keys[i], values.get(i));
				}
				for(int i=0; i<size; i++){
					String s=map.get(keys[i]);
					if(s!=null){sum+=s.length();}
				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+map.size()+", sum="+sum);
			Shared.printMemory();
			map=null;
			System.gc();
		}

		{
			System.err.println("\n*** HashMap<Long, String> ***");
			Shared.printMemory();
			t.start();
			HashMap<Long, String> map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new HashMap<Long, String>();
				for(int i=0; i<size; i++){
					map.put(keys[i], values.get(i));
				}
				for(int i=0; i<size; i++){
					String s=map.get(keys[i]);
					if(s!=null){sum+=s.length();}
				}
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
	 * @param type Component class of the VALUES (e.g. String.class, byte[].class, long[].class) — used to allocate a
	 *        real V[] backing array via reflection so values() returns a genuine typed array (a plain Object[] backing
	 *        CCEs at the call site under a concrete V). Pass a class literal; it is a cheap cached constant. [map/LongObjectMap#002]
	 */
	public LongObjectMap(Class<V> type){
		this(256, type);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param type Component class of the values (see {@link #LongObjectMap(Class)}).
	 */
	public LongObjectMap(int initialSize, Class<V> type){
		this(initialSize, 0.7f, type);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param loadFactor Load factor (0.25-0.90) - map resizes when size exceeds capacity*loadFactor
	 * @param type Component class of the values (see {@link #LongObjectMap(Class)}).
	 */
	public LongObjectMap(int initialSize, float loadFactor_, Class<V> type){
		invalid=randy.nextLong()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
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
		Arrays.fill(keys, invalid);
		Arrays.fill(values, null);
		size=0;
	}

	/**
	 * Gets the value associated with the given key.
	 * @param key Key to look up
	 * @return Value associated with key, or null if key not found
	 */
	public V get(long key){
		int cell=findCell(key);
		return cell<0 ? null : values[cell];
	}

	/**
	 * Checks if the map contains the given key.
	 * @param key Key to check
	 * @return true if key is present
	 */
	public boolean contains(long key){
		return findCell(key)>=0;
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Key to insert/update
	 * @param value Value to associate with key (may be null)
	 * @return Previous value associated with key, or null if key was not present
	 */
	public V put(long key, V value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * @param map Source map to copy from
	 */
	public void putAll(LongObjectMap<V> map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				put(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Key to insert/update
	 * @param value Value to associate with key (may be null)
	 * @return Previous value associated with key, or null if key was not present
	 */
	public V set(long key, V value){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final V oldV=values[cell];
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
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
	public V remove(long key){
		if(key==invalid){return null;}
		final int cell=findCell(key);
		if(cell<0){return null;}
		assert(keys[cell]==key);
		final V oldV=values[cell];
		keys[cell]=invalid;
		values[cell]=null;
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
			final long key=keys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final long key=keys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
	}

	/**
	 * Attempts to move an entry to its ideal position.
	 * @param cell Position of entry to rehash
	 * @return true if entry was moved
	 */
	private boolean rehashCell(final int cell){
		final long key=keys[cell];
		final V value=values[cell];
		assert(key!=invalid);
		final int dest=findCellOrEmpty(key);
		if(cell==dest){return false;}
		assert(keys[dest]==invalid);
		keys[cell]=invalid;
		values[cell]=null;
		keys[dest]=key;
		values[dest]=value;

		return true;
	}

	/**
	 * Resets the invalid sentinel when it collides with a real key.
	 */
	private void resetInvalid(){
		final long old=invalid;
		long x=invalid;
		while(x==old || contains(x)){x=randy.nextLong()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				keys[i]=invalid;
			}
		}
	}

	/**
	 * Finds the cell containing the given key.
	 * Uses linear probing with hash mixing for collision resolution.
	 * Hash mixing ensures adjacent kmers distribute evenly.
	 * @param key Key to search for
	 * @return Cell index if found, -1 if not found
	 */
	private int findCell(final long key){
		if(key==invalid){return -1;}

		final int limit=keys.length;
		final long hash=Tools.hash64plus(key);
		final int initial=(int)(hash & mask);

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
	 * Finds the cell containing the key, or an empty cell where it can be inserted.
	 * @param key Key to search for
	 * @return Cell index (either containing key or empty)
	 */
	private int findCellOrEmpty(final long key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";

		//CLEVER [verified]: no separate occupancy array - empty cells hold the negative `invalid` sentinel and keys are
		//compared directly (x==key || x==invalid), one long compare per cell. initial=(int)(hash & mask) is in range for
		//ANY mask: mask is a non-negative int (pow2-1, <=Integer.MAX_VALUE), so hash64plus(key) & mask <= mask < 2^31 and
		//the (int) cast can't go negative. set() guards key==invalid via resetInvalid() (re-rolls the sentinel + remaps
		//old cells), so a real key can never alias an empty cell.
		final int limit=keys.length;
		final long hash=Tools.hash64plus(key);
		final int initial=(int)(hash & mask);

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

		final long[] oldK=keys;
		@SuppressWarnings("unchecked")
		final V[] oldV=values;
		keys=KillSwitch.allocLong1D((int)size3);
		//Option A [map/LongObjectMap#002]: allocate a REAL V[] via reflection from the stored component type, so the
		//backing array's runtime type is genuinely V[] (e.g. String[], byte[][]) and values() returns a usable typed
		//array instead of an Object[] that CCEs at the call site. `type` is a cheap cached class literal.
		@SuppressWarnings("unchecked")
		V[] temp=(V[])Array.newInstance(type, (int)size3);
		values=temp;
		Arrays.fill(keys, invalid);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final long k=oldK[i];
			if(k!=invalid){
				final V v=oldV[i];
				set(k, v);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns array of all keys currently in the map.
	 * @return Array containing all keys (no invalid sentinels)
	 */
	public long[] toArray(){
		long[] x=KillSwitch.allocLong1D(size);
		int i=0;
		for(long key : keys){
			if(key!=invalid){
				x[i]=key;
				i++;
			}
		}
		return x;
	}

	/**
	 * Returns the internal key array.
	 * WARNING: Contains invalid sentinel for empty cells. Use with caution.
	 * @return Internal key array
	 */
	public long[] keys(){return keys;}

	/**
	 * Returns the internal value array (the live backing array; null entries mark empty cells).
	 * This is a GENUINE V[] (e.g. String[], byte[][]) because the backing is allocated via Array.newInstance(type,...)
	 * from the component class passed to the constructor [map/LongObjectMap#002] — so it can be used directly (assigned
	 * to V[], iterated, .length) without the ClassCastException a plain Object[] backing would cause at the call site.
	 * WARNING: live array, contains nulls for empty cells. Use with caution.
	 * @return Internal value array (a real V[])
	 */
	public V[] values(){return values;}

	/**
	 * Returns the invalid sentinel value used for empty cells.
	 * @return Invalid sentinel (always negative)
	 */
	public long invalid(){return invalid;}

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

	/** Array of long keys (invalid sentinel for empty cells) */
	private long[] keys;
	/** Array of Object values (parallel to keys array, may contain nulls) */
	private V[] values;
	/** Number of entries in the map */
	private int size=0;
	/** Invalid sentinel value for empty cells (always negative) */
	private long invalid;
	/** Bit mask for fast modulo (always power of 2 minus 1) */
	private int mask;
	/** Size threshold for triggering resize (capacity * loadFactor) */
	private int sizeLimit;
	/** Load factor (fraction of capacity before resize) */
	private final float loadFactor;
	/** Component class of V, for allocating a real V[] backing via reflection (Array.newInstance). [map/LongObjectMap#002] */
	private final Class<V> type;

	/** Extra space beyond power-of-2 size to reduce wrap-around collisions */
	private static final int extra=10;
	/** Mask to ensure invalid sentinel is negative */
	private static final long MINMASK=Long.MIN_VALUE;

	/** Random number generator for invalid sentinel */
	private static final Random randy=shared.Shared.random(1);

}