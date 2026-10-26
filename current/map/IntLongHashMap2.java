package map;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;

/**
 * Hash map with primitive int keys and long values.
 * Uses open addressing with linear probing for collision resolution.
 * Uses power-of-two logical sizing and the unrestricted Tools.hash32shift mixer.
 * Parallel primitive arrays avoid boxed entries. All int keys and long values are
 * supported; a negative empty-cell marker is replaced if inserted as a key.
 * Increment saturates upward but wraps downward; getAndIncrement wraps in both
 * directions. Unsynchronized: mutations require exclusive access.
 * @author Isla
 * @contributor Amber
 * @date February 4, 2026
 */
public final class IntLongHashMap2 implements Serializable{

	/** Serialization version for this map layout. */
	private static final long serialVersionUID=1L;

	/** Runs a legacy fill/lookup/removal benchmark with generated keys and reports to stderr.
	 * @param args Optional positional nonnegative key count and positive repetition count */
	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;

		// Generate random keys with collisions
		System.err.println("Generating "+size+" random keys...");
		Random randy=shared.Shared.random(12345);
		int[] keys=new int[size];
		for(int i=0; i<size; i++){
			keys[i]=(int)Math.sqrt(randy.nextLong()&Long.MAX_VALUE);
		}
		System.err.println("Keys generated.");

		bench(keys, repeats);
	}

	/** Reconstructs the map and runs fill/lookup/removal cycles, requesting GC between reports.
	 * Duplicate input keys overwrite equal-key values; the key array is not changed.
	 * @param keys Non-null input keys
	 * @param repeats Positive number of measured cycles */
	private static void bench(int[] keys, int repeats){
		final int size=keys.length;
		System.gc();
		Timer t=new Timer();
		boolean remove=true;
		{
			System.err.println("\n*** IntLongHashMap2 ***");
			Shared.printMemory();
			t.start();
			IntLongHashMap2 map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new IntLongHashMap2();
				for(int i=0; i<size; i++){map.put(keys[i], 2*keys[i]);}
				for(int i=0; i<size; i++){sum+=map.get(keys[i]);}
				if(remove){
					for(int i=0; i<size; i++){map.remove(keys[i]);}
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
	public IntLongHashMap2(){
		this(256);
	}

	/**
	 * Creates a new map with specified initial capacity and default load factor (0.7).
	 * @param initialSize Positive logical request, rounded up to a power of two; requests above 2^30 fail
	 */
	public IntLongHashMap2(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Positive logical request, rounded up to a power of two; requests above 2^30 fail
	 * @param loadFactor_ Requested fraction in (0,1), clamped to [0.25,0.90] and stored as float
	 */
	public IntLongHashMap2(int initialSize, double loadFactor_){
		invalid=randy.nextInt()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=(float)Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Removes all mappings without shrinking capacity or changing the sentinel.
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
	 * @return Stored value or -1 when absent; contains() disambiguates a stored -1
	 */
	public long get(int key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
	}

	/**
	 * Checks if the map contains the given key.
	 * @param key Key to check
	 * @return true if key is present
	 */
	public boolean contains(int key){
		return findCell(key)>=0;
	}

	/** Alias for contains(int).
	 * @param key Integer key to query
	 * @return True when key is present */
	public boolean containsKey(int key){
		return findCell(key)>=0;
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Key to insert/update
	 * @param value Value to associate with key
	 * @return Previous value, or 0 if key was not present (empty cells hold 0;
	 *         use contains() to disambiguate a stored 0). NOT -1. [map/IntLongHashMap2#001]
	 */
	public long put(int key, long value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * Shared keys take the source value; destination-only mappings remain.
	 * A distinct source is preserved; self-copy is a no-op.
	 * @param map Non-null source stable during this operation
	 */
	public void putAll(IntLongHashMap2 map){
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
	 * @return Previous value, or 0 if key was not present (empty cells hold 0;
	 *         use contains() to disambiguate a stored 0). NOT -1. [map/IntLongHashMap2#001]
	 */
	public long set(int key, long value){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final long oldV=values[cell];
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
	 * An existing Long.MAX_VALUE count remains saturated.
	 * @param key Key to increment
	 * @return New value after increment
	 */
	public long increment(int key){
		return increment(key, 1);
	}

	/**
	 * Increments the value associated with the key by the specified amount.
	 * If key is not present, inserts it with the increment value.
	 * Positive overflows clamp to Long.MAX_VALUE; negative underflows wrap unchanged.
	 * @param key Key to increment
	 * @param incr Amount to add (can be negative)
	 * @return New value after increment
	 */
	public long increment(int key, long incr){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final long oldV=values[cell];
		final long value=(incr>0 && oldV>Long.MAX_VALUE-incr) ? Long.MAX_VALUE : oldV+incr;
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
		return value;
	}

	/**
	 * Returns the old value and stores old+1 with a single lookup.
	 * Unlike increment(), addition wraps at Long.MAX_VALUE rather than saturating.
	 * An absent key starts from zero and is inserted with value one.
	 * @param key Key to find
	 * @return Old value associated with key, or 0 if key was not present (and is now inserted with value 1)
	 */
	public long getAndIncrement(int key){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final long oldV=values[cell];

		// Note: We don't clamp here because this is mainly for packed indices
		// where we expect controlled incrementing.
		values[cell]=oldV+1;

		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
		return oldV;
	}

	/**
	 * Increments all entries in this map by corresponding values from another map.
	 * Keys not present in this map are inserted with their values from the source map.
	 * Uses upper-saturating/lower-wrapping increments. A distinct source is preserved;
	 * self-merging doubles each value according to those arithmetic rules.
	 * @param map Non-null source stable during this call unless it is this map itself
	 */
	public void incrementAll(IntLongHashMap2 map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				increment(map.keys[i], map.values[i]);
			}
		}
	}

	/** Merges the source map by taking the maximum value for each shared key.
	 * Keys absent here are inserted with their exact source values.
	 * Destination-only mappings remain. A distinct source is preserved; self-merging is a no-op.
	 * @param map Non-null source stable during this operation
	 */
	public void setToMax(IntLongHashMap2 map){
		for(int i=0; i<map.keys.length; i++){
			final int key=map.keys[i];
			if(key!=map.invalid){
				final long value=map.values[i];
				if(key==invalid){put(key, value); continue;}
				final int cell=findCellOrEmpty(key);
				if(keys[cell]==invalid){
					values[cell]=value;
					keys[cell]=key;
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
	public boolean remove(int key){
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
			final int key=keys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final int key=keys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
	}

	/**
	 * Moves an entry to the first reachable hole after deletion, preserving its value and size.
	 * @param cell Occupied cell index
	 * @return true if entry was moved
	 */
	private boolean rehashCell(final int cell){
		final int key=keys[cell];
		final long value=values[cell];
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

	/** Chooses an unused negative sentinel and rewrites empty keys, preserving live entries. */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		Vector.changeAll(keys, old, x);
	}

	/**
	 * Finds the cell containing the given key.
	 * Starts from Tools.hash32shift(key) masked to the logical table, then probes linearly.
	 * @param key Key to search for
	 * @return Cell index if found, -1 if not found
	 */
	private int findCell(final int key){
		if(key==invalid){return -1;}
		// UPDATED: Using unrestricted hash32shift for better speed
		final int hash=Tools.hash32shift(key);
		final int initial=hash&mask;
		return Vector.findKeyScalar(keys, key, initial, invalid);
	}

	/**
	 * Finds the cell containing the key, or an empty cell where it can be inserted.
	 * @param key Integer key unequal to the current empty-cell sentinel
	 * @return Cell index (either containing key or empty)
	 */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		// UPDATED: Using unrestricted hash32shift for better speed
		final int hash=Tools.hash32shift(key);
		final int initial=hash&mask;
		return Vector.findKeyOrInvalidScalar(keys, key, initial, invalid);
	}

	/**
	 * Doubles the map capacity and rehashes all entries.
	 */
	private final void resize(){
		assert(size>=sizeLimit);
		//grow 2x: double the LOGICAL pow2 capacity (keys.length-extra==pow2), not keys.length (=pow2+extra).
		//keys.length*2L overshot the power-of-2 boundary and rounded up to 4*pow2 -> fixed [map/IntLongHashMap2#002].
		resize((keys.length-extra)*2L);
	}

	/**
	 * Rebuilds both arrays at a rounded-up power-of-two logical size plus slack.
	 * Preserves mappings but replaces borrowed arrays and can change iteration order.
	 * The largest positive int power-of-two tier is 2^30, with load at least 0.85;
	 * the next doubling fails the positive-size assertion rather than capping storage.
	 * @param size2 Positive requested logical size, greater than the current entry count
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;

		//No physical cap is applied: the positive-size assertion bounds the logical tier at 2^30.
		int newSize=Integer.highestOneBit((int)size2);
		if(newSize<size2){newSize<<=1;}
		assert(newSize>0) : "Overflow: "+size2;

		final int size3=newSize+extra;
		final float lf=(newSize<0x40000000 ? loadFactor : Math.max(loadFactor, 0.85f));
		sizeLimit=(int)(newSize*lf);
		mask=newSize-1;

		final int[] oldK=keys;
		final long[] oldV=values;
		keys=KillSwitch.allocInt1D(size3);
		values=KillSwitch.allocLong1D(size3);
		Arrays.fill(keys, invalid);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i];
			if(k!=invalid){
				final long v=oldV[i];
				set(k, v);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies all stored keys into a dense array in backing-table order, which is not sorted.
	 * @return New independent array of size() keys without empty-cell markers
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
	 * Borrows the mutable key array; external mutation can invalidate probe invariants.
	 * Resizing replaces the array; clear() mutates the same array in place.
	 * @return Map-owned keys with invalid() in empty cells
	 */
	public int[] keys(){return keys;}

	/**
	 * Borrows the mutable value array parallel to keys(); empty cells contain zero.
	 * Resizing replaces the array; clear() mutates the same array in place.
	 * @return Map-owned values; use the parallel keys to determine occupancy
	 */
	public long[] values(){return values;}

	/**
	 * Returns the invalid sentinel value used for empty cells.
	 * @return Current negative marker, which changes if inserted as a real key
	 */
	public int invalid(){return invalid;}

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

	/** Array of keys (invalid sentinel for empty cells) */
	private int[] keys;
	/** Array of values (parallel to keys array) */
	private long[] values;
	/** Number of entries in the map */
	private int size=0;
	/** Invalid sentinel value for empty cells (always negative) */
	private int invalid;
	/** Bit mask for fast modulo (always power of 2 minus 1) */
	private int mask;
	/** Size threshold for triggering resize (capacity * loadFactor) */
	private int sizeLimit;
	/** Requested load fraction clamped to [0.25,0.90]; the 2^30 tier uses at least 0.85. */
	private final float loadFactor;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Extra physical probe cells beyond the power-of-two logical size. */
	private static final int extra=10;
	/** Mask to ensure invalid sentinel is negative */
	private static final int MINMASK=Integer.MIN_VALUE;

	/** Random number generator for invalid sentinel */
	private static final Random randy=shared.Shared.random(1);

}
