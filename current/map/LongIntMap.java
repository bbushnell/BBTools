package map;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Hash map with primitive long keys and int values.
 * Uses open addressing with linear probing for collision resolution.
 * Uses power-of-2 sizing and hash mixing for fast, well-distributed hashing.
 *
 * Ideal for kmer counting where keys are kmers (long) and values are counts (int).
 *
 * Parallel primitive arrays avoid boxed keys, boxed values and per-entry objects.
 *
 * Supports every long key and int value; an absent lookup returns -1, so use
 * contains() when that sentinel can also be a stored value.
 * Not thread-safe; publish an immutable map for concurrent reads or synchronize
 * access externally. Backing arrays are borrowed and may change during resizing.
 *
 * @author Isla
 * @date November 2, 2025
 */
public final class LongIntMap implements Serializable{

	/** Serialization version for this map layout. */
	private static final long serialVersionUID=1L;

	/** Runs the legacy allocation/put/get comparison and prints timings to stderr.
	 * @param args Optional positional nonnegative key count and positive repetition count */
	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;

		// Generate random keys with collisions
		System.err.println("Generating "+size+" random keys...");
		Random randy=shared.Shared.random(12345);
		long[] keys=new long[size];
		for(int i=0; i<size; i++){
			long x=randy.nextLong()&Long.MAX_VALUE;
			keys[i]=(long)Math.sqrt(x);
		}
		System.err.println("Keys generated.");

		bench(keys, repeats);
	}

	/** Compares repeated construction and lookup with a boxed HashMap.
	 * Duplicate keys overwrite earlier indices in both maps; keys are not modified.
	 * @param keys Non-null input keys, retained only for the duration of this call
	 * @param repeats Positive number of construction/lookup cycles for each map */
	private static void bench(long[] keys, int repeats){
		final int size=keys.length;
		System.gc();
		Timer t=new Timer();

		{
			System.err.println("\n*** LongIntMap ***");
			Shared.printMemory();
			t.start();
			LongIntMap map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new LongIntMap();
				for(int i=0; i<size; i++){
					map.put(keys[i], i);
				}
				for(int i=0; i<size; i++){
					sum+=map.get(keys[i]);
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
			System.err.println("\n*** HashMap<Long, Integer> ***");
			Shared.printMemory();
			t.start();
			HashMap<Long, Integer> map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new HashMap<Long, Integer>();
				for(int i=0; i<size; i++){
					map.put(keys[i], i);
				}
				for(int i=0; i<size; i++){
					sum+=map.get(keys[i]);
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
	 */
	public LongIntMap(){
		this(256);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Positive minimum logical table size, rounded up to a power of two
	 */
	public LongIntMap(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Positive minimum logical table size, rounded up to a power of two
	 * @param loadFactor_ Requested load factor in (0,1), clamped to [0.25,0.90]
	 */
	public LongIntMap(int initialSize, float loadFactor_){
		invalid=randy.nextLong()|MINMASK;
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
	 * Removes all entries without changing allocated capacity or the sentinel.
	 */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
	}

	/**
	 * Gets the value associated with the given key.
	 * @param key Key to look up
	 * @return Value associated with key, or -1 if key not found
	 */
	public int get(long key){
		int cell=findCell(key);
		//-1 also can be stored; contains() distinguishes absence when values are unrestricted.
		return cell<0 ? -1 : values[cell];
	}

	/**
	 * Gets the value associated with {@code key}, returning the caller-supplied
	 * default when the key is absent. All int values, including -1, may be stored;
	 * this lets bit-vector callers choose zero as their absent value. A stored value
	 * can still equal the default; use contains() when presence itself matters.
	 * @param key Key to look up
	 * @param defaultValue Value returned when the key is absent
	 * @return Stored value, or {@code defaultValue} when absent
	 */
	public int getOrDefault(long key, int defaultValue){
		int cell=findCell(key);
		return cell<0 ? defaultValue : values[cell];
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
	 * @param value Value to associate with key
	 * @return Previous value associated with key, or 0 if key was not present.
	 *         Stored zero is indistinguishable from absence in this return value.
	 */
	public int put(long key, int value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * Shared keys take the source value; destination-only entries remain. The source
	 * is neither retained nor modified, and copying this map onto itself is a no-op.
	 * @param map Non-null source, stable for the duration of this operation
	 */
	public void putAll(LongIntMap map){
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
	 * @param value Value to associate with key
	 * @return Previous value associated with key, or 0 if key was not present.
	 *         Stored zero is indistinguishable from absence in this return value.
	 */
	public int set(long key, int value){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		//Empty slots contain zero after allocation, removal or clearing; key presence is separate from value.
		final int oldV=values[cell];
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
		return oldV;
	}

	/**
	 * Increments the value associated with the key by 1.
	 * If key is not present, inserts it with value 1.
	 * An existing Integer.MAX_VALUE count remains saturated.
	 * @param key Key to increment
	 * @return New value after increment
	 */
	public int increment(long key){
		return increment(key, 1);
	}

	/**
	 * Increments the value associated with the key by the specified amount.
	 * If key is not present, inserts it with the increment value.
	 * Positive overflow is clamped to Integer.MAX_VALUE. Negative underflow retains
	 * the historical int narrowing behavior; this is not a lower-saturating counter.
	 * @param key Key to increment
	 * @param incr Amount to add (can be negative)
	 * @return New value after increment
	 */
	public int increment(long key, int incr){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final long oldV=values[cell];
		//Widen before addition so upper overflow is detected before narrowing; absent slots start at zero.
		final int value=(int)Tools.min(Integer.MAX_VALUE, oldV+incr);
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
		return value;//A resize can move the key away from the original cell.
	}

	/**
	 * Increments all entries in this map by corresponding values from another map.
	 * Keys not present in this map are inserted with their values from the source map.
	 * Uses increment(long,int)'s upper saturation and lower wrapping behavior.
	 * A distinct source is preserved; using this map doubles each value with those rules.
	 * @param map Non-null source, externally stable unless it is this map itself
	 */
	public void incrementAll(LongIntMap map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				increment(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * For each key in the source map, sets this map's value to the maximum
	 * of the current value and the source value.
	 * Keys absent from this map copy the source value exactly, including negative values.
	 * Destination-only keys remain; a distinct source is preserved and self-merging is a no-op.
	 * @param map Non-null source, stable for the duration of this operation
	 */
	public void setToMax(LongIntMap map){
		for(int i=0; i<map.keys.length; i++){
			final long key=map.keys[i];
			if(key!=map.invalid){
				final int value=map.values[i];
				if(key==invalid){put(key, value); continue;}
				final int cell=findCellOrEmpty(key);
				if(keys[cell]==invalid){
					keys[cell]=key;
					values[cell]=value;
					size++;
					if(size>sizeLimit){resize();}
				}else{values[cell]=Tools.max(value, values[cell]);}
			}
		}
	}

	/**
	 * Removes the entry with the specified key.
	 * @param key Key to remove
	 * @return true if key was present and removed, false if key was not found
	 */
	public boolean remove(long key){
		if(key==invalid){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		values[cell]=0;
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
	 * Moves an entry to the first reachable hole in its probe sequence after deletion.
	 * Preserves its value and the entry count; the vacated slot becomes empty.
	 * @param cell Occupied cell index
	 * @return true if entry was moved
	 */
	private boolean rehashCell(final int cell){
		final long key=keys[cell];
		final int value=values[cell];
		assert(key!=invalid);
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
	 * Chooses a different negative sentinel absent from the map and rewrites empty keys.
	 * Existing entries, values and size are unchanged; borrowed key-array contents change.
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
	 * Stops at the first empty cell and wraps at the physical array end.
	 * @param key Key to search for
	 * @return Cell index if found, -1 if not found
	 */
	private int findCell(final long key){
		if(key==invalid){return -1;}

		final int limit=keys.length;
		final long hash=Tools.hash64plus(key);
		//two-segment wrap-around scan: [initial,limit) then [0,initial); returns -1 at the first invalid cell.
		final int initial=(int)(hash&mask);

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
	 * Requires at least one empty cell and preserves all map state.
	 * @param key Key to search for, unequal to the current empty-cell sentinel
	 * @return Cell index (either containing key or empty)
	 * @throws RuntimeException If a full table contains neither the key nor an empty cell
	 */
	private int findCellOrEmpty(final long key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";

		final int limit=keys.length;
		final long hash=Tools.hash64plus(key);
		//Tools.hash64plus caps the initial index below SAFE_ARRAY_LEN, including at the terminal tier.
		//The load threshold leaves an empty cell, so an intact table always terminates its probe.
		final int initial=(int)(hash&mask);

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
		//Exclude slack before doubling; otherwise power-of-two rounding can quadruple capacity.
		//This preserves ordinary twofold growth. [map/IntLongHashMap2#002]
		resize((keys.length-extra)*2L);
	}

	/**
	 * Rebuilds the arrays, preserving entries and replacing borrowed-array identities.
	 * Rounds the requested logical size up to a power of two, adds slack, then caps
	 * the physical length at Shared.SAFE_ARRAY_LEN. At the cap, the load factor is 0.85.
	 * @param size2 Positive requested logical size, strictly greater than the current entry count
	 * @throws RuntimeException If the capped allocation cannot grow the existing table
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;

		long size3=Long.highestOneBit(size2);
		if(size3<size2){size3<<=1;}
		//Keep the logical mask when physical length is capped: hash64plus already limits
		//the initial index below SAFE_ARRAY_LEN. The mask need not be capped. [map/LongIntMap#002]
		mask=(int)(size3-1);
		size3=Math.min(size3+extra, Shared.SAFE_ARRAY_LEN);
		if((keys!=null && size3<=keys.length) || size3>Shared.SAFE_ARRAY_LEN){
			throw new RuntimeException("Map hit capacity at "+size);
		}
		final float loadFactor2=(size3<Shared.SAFE_ARRAY_LEN ? loadFactor : 0.85f);
		sizeLimit=(int)((size3-extra)*loadFactor2);

		final long[] oldK=keys;
		final int[] oldV=values;
		keys=KillSwitch.allocLong1D((int)size3);
		values=KillSwitch.allocInt1D((int)size3);
		Arrays.fill(keys, invalid);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final long k=oldK[i];
			if(k!=invalid){
				final int v=oldV[i];
				set(k, v);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies keys into a new dense array in backing-table order, which is not sorted.
	 * @return Independently owned array of size() keys, with no empty-cell sentinels
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
	 * @return Borrowed internal key array; do not mutate it or retain it across resizing
	 */
	public long[] keys(){return keys;}

	/**
	 * Returns the internal value array.
	 * WARNING: Contains 0 for empty cells. Use with caution.
	 * @return Borrowed values parallel to keys(); only occupied key slots are entries
	 */
	public int[] values(){return values;}

	/**
	 * Returns the invalid sentinel value used for empty cells.
	 * @return Current negative sentinel; it may change when that long is inserted as a key
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
	/** Array of int values (parallel to keys array) */
	private int[] values;
	/** Number of entries in the map */
	private int size=0;
	/** Invalid sentinel value for empty cells (always negative) */
	private long invalid;
	/** Bit mask for fast modulo (always power of 2 minus 1) */
	private int mask;
	/** Maximum entry count before growth; derived from logical capacity and effective load. */
	private int sizeLimit;
	/** Requested load fraction clamped to [0.25, 0.90]; the terminal allocation uses 0.85. */
	private final float loadFactor;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Extra physical probe cells beyond the power-of-two logical capacity. */
	private static final int extra=10;
	/** Mask to ensure invalid sentinel is negative */
	private static final long MINMASK=Long.MIN_VALUE;

	/** Random number generator for invalid sentinel */
	private static final Random randy=shared.Shared.random(1);

}
