package map;

import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;

/**
 * Int-to-int map that mixes int keys into a power-of-two open-addressed table.
 * Empty cells hold a rerollable negative sentinel, values in empty cells are zero,
 * get() returns -1 on misses, and put()/set() return zero on fresh insertions.
 * Positive increment overflow is capped at {@link Integer#MAX_VALUE}; negative
 * underflow wraps through signed int truncation. Instances are unsynchronized.
 *
 * @author Isla
 * @date November 2, 2025
 */
public final class IntHashMap2 implements IntHashMapInterface{

	/** Serialization version. */
	private static final long serialVersionUID=1L;

	/**
	 * Runs a destructive timing comparison against {@link IntHashMap} and {@link HashMap}.
	 * @param args Optional element count and repeat count
	 */
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

	/**
	 * Times repeated load/lookup/remove cycles on the supplied keys.
	 * @param keys Randomized keys to insert and query
	 * @param repeats Positive number of map construction cycles per implementation
	 */
	private static void bench(int[] keys, int repeats){
		final int size=keys.length;
		System.gc();
		Timer t=new Timer();
		boolean remove=true;
		{
			System.err.println("\n*** IntHashMap2 ***");
			Shared.printMemory();
			t.start();
			IntHashMap2 map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new IntHashMap2();
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

		{
			System.err.println("\n*** IntHashMap ***");
			Shared.printMemory();
			t.start();
			IntHashMap map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new IntHashMap();
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

		{
			System.err.println("\n*** HashMap<Integer, Integer> ***");
			Shared.printMemory();
			t.start();
			HashMap<Integer, Integer> map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new HashMap<Integer, Integer>();
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

	/** Creates a map with the default initial size and load factor. */
	public IntHashMap2(){
		this(256);
	}

	/**
	 * Creates a map with default load factor 0.7.
	 * @param initialSize Positive requested logical minimum size, rounded up to the next power of two
	 */
	public IntHashMap2(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Creates a map with caller-selected initial size and load factor.
	 * The load factor is clamped to [0.25,0.90].
	 * @param initialSize Positive requested logical minimum size, rounded up to the next power of two
	 * @param loadFactor_ Requested load in (0,1), clamped to [0.25,0.90] after validation
	 */
	public IntHashMap2(int initialSize, float loadFactor_){
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

	/** Removes all mappings, resets values to zero, and preserves the allocated raw arrays. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
	}

	/**
	 * Returns the value mapped from a key.
	 * @param key Key to look up
	 * @return Mapped value, or -1 if absent
	 */
	public int get(int key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
	}

	/**
	 * Tests whether a key is currently mapped.
	 * @param key Key to check
	 * @return true if {@code key} is live in this map
	 */
	public boolean contains(int key){
		return findCell(key)>=0;
	}

	/**
	 * Tests whether a key is currently mapped.
	 * @param key Key to check
	 * @return true if {@code key} is live in this map
	 */
	public boolean containsKey(int key){
		return findCell(key)>=0;
	}

	/**
	 * Returns 1 if the key is found, else 0.
	 * @param key Key to check
	 * @return 1 if the key is present, 0 if not
	 */
	public int containsKeyBinary(int key){
		//CLEVER [verified in-file]: branchless. findCell returns cell>=0 (found) or -1 (absent).
		//~cell>>>31 == 1 for any cell>=0 (sign bit of ~cell set), == 0 for cell==-1 (~-1==0).
		return (~(findCell(key))>>>31);
	}

	/**
	 * Stores a key/value pair.
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous value, or 0 if key was not present (empty cells hold 0;
	 *         use contains() to disambiguate a stored 0). NOT -1.
	 */
	public int put(int key, int value){
		return set(key, value);
	}

	/**
	 * Copies all live mappings from another map, overwriting matching keys here.
	 * @param map Source map to read without modifying
	 */
	public void putAll(IntHashMap2 map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				put(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * Stores a key/value pair, rerolling the sentinel if it collides with {@code key}.
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous value, or 0 if key was not present (empty cells hold 0;
	 *         use contains() to disambiguate a stored 0). NOT -1.
	 */
	public int set(int key, int value){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final int oldV=values[cell];//0 for a fresh (invalid) cell -> set/put return 0 when absent
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
		return oldV;
	}

	/**
	 * Adds one to a key's value, inserting 1 if absent.
	 * @param key Key to insert or update
	 * @return New value, capped at {@link Integer#MAX_VALUE} on positive overflow
	 */
	public int increment(int key){
		return increment(key, 1);
	}

	/**
	 * Adds a signed delta to a key's value, inserting the delta if absent.
	 * Positive overflow caps at {@link Integer#MAX_VALUE}; negative underflow wraps.
	 * @param key Key to insert or update
	 * @param incr Signed increment to add
	 * @return New value after overflow handling
	 */
	public int increment(int key, int incr){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final long oldV=values[cell];
		final int value=(int)Tools.min(Integer.MAX_VALUE, oldV+incr);
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
		return value;
	}

	/**
	 * Adds all live values from another map into this map.
	 * If {@code map} is this map, each value is doubled as its occupied cell is visited;
	 * distinct sources must remain structurally stable while their raw arrays are scanned.
	 * @param map Source map to scan for increments; may be this map
	 */
	public void incrementAll(IntHashMap2 map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				increment(map.keys[i], map.values[i]);
			}
		}
	}

	/** Merges the source map by taking the maximum value for each shared key.
	 * Keys absent here are inserted with their exact source values.
	 * @param map Source map to read without modifying
	 */
	public void setToMax(IntHashMap2 map){
		for(int i=0; i<map.keys.length; i++){
			final int key=map.keys[i];
			if(key!=map.invalid){
				final int value=map.values[i];
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
	 * Removes a key and closes the following probe cluster.
	 * @param key Key to remove
	 * @return true if a live mapping was removed
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
	 * Rehashes occupied cells after a deletion until the next empty cell.
	 * @param initial Cell that was just emptied
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
	 * Reinserts one occupied cell if its key can move earlier in its probe sequence.
	 * @param cell Occupied raw cell to test
	 * @return true if the mapping moved to a different cell
	 */
	private boolean rehashCell(final int cell){
		final int key=keys[cell];
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

	/** Rerolls the empty-cell sentinel and rewrites all cells holding the old sentinel. */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		Vector.changeAll(keys, old, x);
	}

	/**
	 * Locates a key by probing from {@code Tools.hash32plus(key) & mask}.
	 * @param key Key to search for
	 * @return Raw cell index if found, or -1 if absent
	 */
	public int findCell(final int key){
		if(key==invalid){return -1;}
		final int hash=Tools.hash32plus(key);
		final int initial=hash&mask;
		return Vector.findKeyScalar(keys, key, initial, invalid);
	}

	/**
	 * Locates a key or the empty cell that would receive it.
	 * @param key Non-sentinel key to search for
	 * @return Raw cell containing {@code key}, or a raw cell holding {@link #invalid}
	 */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		final int hash=Tools.hash32plus(key);
		final int initial=hash&mask;
		//probe (scalar, Java-8-safe) wraps [initial,len)->[0,initial); the extra=10 slots merely let
		//most probes skip the wrap branch -> correctness rests on loadFactor<1, not on `extra`.
		return Vector.findKeyOrInvalidScalar(keys, key, initial, invalid);
	}

	/** Doubles the logical power-of-two table size and rehashes live mappings. */
	private final void resize(){
		assert(size>=sizeLimit);
		//grow 2x: double the LOGICAL pow2 capacity (keys.length-extra==pow2), not keys.length (=pow2+extra).
		//keys.length*2L overshot the power-of-2 boundary and rounded up to 4*pow2 -> fixed to true 2x growth
		//(shared pow2-family over-allocation, anchor [map/IntLongHashMap2#002]). Output unchanged; ~2x less memory.
		resize((keys.length-extra)*2L);
	}

	/**
	 * Replaces raw storage with at least the specified logical size and rehashes live mappings.
	 * Requests near the int-array ceiling crash loudly rather than silently truncating.
	 * @param size2 Requested logical minimum size before power-of-two rounding
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;

		// Round up to next power of 2
		int newSize=Integer.highestOneBit((int)size2);
		if(newSize<size2){newSize<<=1;}
		//crash-loud ceiling: past ~2^30 entries the next doubling overflows the (int)size2 cast
		//-> newSize<=0 -> this assert fires (requires -ea, always on for BBTools).
		assert(newSize>0) : "Overflow: "+size2;

		final int size3=newSize+extra;
		final float lf=(newSize<0x40000000 ? loadFactor : Math.max(loadFactor, 0.85f));
		sizeLimit=(int)(newSize*lf);
		mask=newSize-1;

		final int[] oldK=keys;
		final int[] oldV=values;
		keys=KillSwitch.allocInt1D(size3);
		values=KillSwitch.allocInt1D(size3);
		Arrays.fill(keys, invalid);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i];
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
	 * Copies live keys into a dense new array.
	 * @return Owned key snapshot in raw table iteration order
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
	 * Borrows the raw key array.
	 * @return Borrowed raw key array with {@link #invalid()} in empty cells
	 */
	public int[] keys(){return keys;}

	/**
	 * Borrows the raw value array.
	 * @return Borrowed raw value array, index-aligned with {@link #keys()}
	 */
	public int[] values(){return values;}

	/**
	 * Returns the current empty-cell sentinel.
	 * @return Current negative empty-cell sentinel
	 */
	public int invalid(){return invalid;}

	/**
	 * Returns the number of live mappings.
	 * @return Number of live mappings
	 */
	public int size(){return size;}

	/**
	 * Tests whether the map has no live mappings.
	 * @return true when the map has no live mappings
	 */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Raw key table; empty cells contain {@link #invalid}. */
	private int[] keys;
	/** Raw int values, index-aligned with {@link #keys}. */
	private int[] values;
	/** Number of live mappings. */
	private int size=0;
	/** Negative empty-cell sentinel absent from the live key set. */
	private int invalid;
	/** Power-of-two mask applied to mixed key hashes. */
	private int mask;
	/** Resize threshold derived from the logical capacity and load factor. */
	private int sizeLimit;
	/** Clamped requested load factor. */
	private final float loadFactor;

	/** Physical overflow cells above the logical power-of-two mask range. */
	private static final int extra=10;
	/** Mask that forces sentinel candidates negative. */
	private static final int MINMASK=Integer.MIN_VALUE;

	/** Reproducible source for negative empty-cell sentinel candidates. */
	private static final Random randy=shared.Shared.random(1);

}
