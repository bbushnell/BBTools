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
 * Uses power-of-2 sizing and hash mixing for fast, well-distributed hashing.
 * * Performance optimized:
 * - Uses unrestricted Tools.hash32shift() for faster hashing.
 * - Stores primitive long values.
 * - Avoids object overhead.
 * * Thread-safety: Not thread-safe. External synchronization required for concurrent access.
 * * @author Isla
 * @contributor Amber
 * @date February 4, 2026
 */
public final class IntLongHashMap2 implements Serializable {

	private static final long serialVersionUID = 1L;

	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;
		
		// Generate random keys with collisions
		System.err.println("Generating "+size+" random keys...");
		Random randy=shared.Shared.random(12345);
		int[] keys=new int[size];
		for(int i=0; i<size; i++){
			keys[i]=(int)Math.sqrt(randy.nextLong() & Long.MAX_VALUE);
		}
		System.err.println("Keys generated.");
		
		bench(keys, repeats);
	}

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
				if(remove) {
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
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 */
	public IntLongHashMap2(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Creates a new map with specified initial capacity and load factor.
	 * @param initialSize Initial capacity (will be rounded up to next power of 2)
	 * @param loadFactor Load factor (0.25-0.90) - map resizes when size exceeds capacity*loadFactor
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
	 * Removes all entries from the map.
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
	
	public boolean containsKey(int key){
		return findCell(key)>=0;
	}

	/**
	 * Associates the specified value with the specified key.
	 * If the key already exists, updates its value.
	 * @param key Key to insert/update
	 * @param value Value to associate with key
	 * @return Previous value associated with key, or -1 if key was not present
	 */
	public long put(int key, long value){
		return set(key, value);
	}

	/**
	 * Copies all entries from another map into this map.
	 * @param map Source map to copy from
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
	 * @return Previous value associated with key, or -1 if key was not present
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
	 * @param key Key to increment
	 * @return New value after increment
	 */
	public long increment(int key){
		return increment(key, 1);
	}

	/**
	 * Increments the value associated with the key by the specified amount.
	 * If key is not present, inserts it with the increment value.
	 * Values are clamped to Long.MAX_VALUE on overflow.
	 * @param key Key to increment
	 * @param incr Amount to add (can be negative)
	 * @return New value after increment
	 */
	public long increment(int key, long incr){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final long oldV=values[cell];
		final long value=oldV+incr;
		values[cell]=Tools.min(Long.MAX_VALUE, value);
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
		return values[cell];
	}
	
	/**
	 * Optimization: Gets the current value and increments it by 1.
	 * Avoids a double lookup (get + set).
	 * Useful for filling PackedIndex.
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
	 * @param map Source map containing increment values
	 */
	public void incrementAll(IntLongHashMap2 map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				increment(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * For each key in the source map, sets this map's value to the maximum
	 * of the current value and the source value.
	 * @param map Source map to compare against
	 */
	public void setToMax(IntLongHashMap2 map){
		for(int i=0; i<map.keys.length; i++){
			final int key=map.keys[i];
			if(key!=map.invalid){
				put(key, Tools.max(map.values[i], get(key)));
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
	 * Attempts to move an entry to its ideal position.
	 * @param cell Position of entry to rehash
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
	
	/** Resets the invalid sentinel when it collides with a real key. */
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
	 * Uses linear probing with hash mixing for collision resolution.
	 * @param key Key to search for
	 * @return Cell index if found, -1 if not found
	 */
	private int findCell(final int key){
		if(key==invalid){return -1;}
		// UPDATED: Using unrestricted hash32shift for better speed
		final int hash=Tools.hash32shift(key);
		final int initial=hash & mask;
		return Vector.findKeyScalar(keys, key, initial, invalid);
	}

	/**
	 * Finds the cell containing the key, or an empty cell where it can be inserted.
	 * @param key Key to search for
	 * @return Cell index (either containing key or empty)
	 */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		// UPDATED: Using unrestricted hash32shift for better speed
		final int hash=Tools.hash32shift(key);
		final int initial=hash & mask;
		return Vector.findKeyOrInvalidScalar(keys, key, initial, invalid);
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

		// Round up to next power of 2
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
	 * Returns array of all keys currently in the map.
	 * @return Array containing all keys (no invalid sentinels)
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
	 * Returns the internal key array.
	 * WARNING: Contains invalid sentinel for empty cells. Use with caution.
	 * @return Internal key array
	 */
	public int[] keys(){return keys;}

	/**
	 * Returns the internal value array.
	 * WARNING: Contains 0 for empty cells. Use with caution.
	 * @return Internal value array
	 */
	public long[] values(){return values;}

	/**
	 * Returns the invalid sentinel value used for empty cells.
	 * @return Invalid sentinel (always negative)
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
	/** Load factor (fraction of capacity before resize) */
	private final float loadFactor;

	/** Extra space beyond power-of-2 size to reduce wrap-around collisions */
	private static final int extra=10;
	/** Mask to ensure invalid sentinel is negative */
	private static final int MINMASK=Integer.MIN_VALUE;

	/** Random number generator for invalid sentinel */
	private static final Random randy=shared.Shared.random(1);

}