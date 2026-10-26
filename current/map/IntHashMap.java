package map;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;

/**
 * Mutable int-to-int hash map backed by parallel primitive arrays.
 * <p>
 * Keys are stored in an open-addressed prime-sized table and resolved by linear
 * probing. One randomly chosen negative int is reserved as the empty-cell
 * sentinel; inserting that value is supported by replacing the sentinel across
 * the table before the insert. Missing-key {@link #get(int)} returns -1, so
 * callers that store -1 must use {@link #contains(int)} or
 * {@link #containsKey(int)} to distinguish absence from a present -1 value.
 * </p>
 * <p>
 * Instances are unsynchronized and mutate internal arrays in place. Raw arrays
 * returned by {@link #keys()} and {@link #values()} are borrowed, remain owned by
 * this map, are replaced by resize, and are retained but filled with empty cells by clear.
 * </p>
 *
 * @author Brian Bushnell
 * @date June 7, 2017
 */
public final class IntHashMap extends AbstractIntHashMap implements Serializable{

	/** Serialization ID for this concrete map layout. */
	private static final long serialVersionUID=5525726007591609843L;

	/**
	 * Runs a local allocation/time comparison against {@code HashMap<Integer,Integer>}.
	 * @param args optional first element overrides the default entry count
	 */
	public static void main(String[] args){
//		IntHashMap set=new IntHashMap(20, 0.7f);
//		test(set);
		int size=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		bench(size);
	}

	/** Fills this map and {@code HashMap<Integer,Integer>} with {@code size} sequential mappings.
	 * @param size number of sequential keys to insert into each map */
	private static void bench(int size){
		System.gc();
		Timer t=new Timer();

		{
			System.err.println("\nIntHashMap:");
			Shared.printMemory();
			t.start();
			IntHashMap map=new IntHashMap();
			for(int i=0; i<size; i++){
				map.put(i, 2*i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			map=null;
			System.gc();
		}

		{
			System.err.println("\nHashMap<Integer, Integer>:");
			Shared.printMemory();
			t.start();
			HashMap<Integer, Integer> map=new HashMap<Integer, Integer>();
			for(int i=0; i<size; i++){
				map.put(i, 2*i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			map=null;
			System.gc();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty map with default capacity. */
	public IntHashMap(){
		this(256);
	}

	/**
	 * Creates an empty map whose backing table can hold at least the requested cells.
	 * @param initialSize minimum prime table size; must be positive
	 */
	public IntHashMap(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Creates an empty map with a chosen minimum table size and resize threshold.
	 * @param initialSize minimum prime table size; must be positive
	 * @param loadFactor_ requested maximum load before resize; must be in (0,1), and is clamped to [0.25,0.90]
	 */
	public IntHashMap(int initialSize, float loadFactor_){
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

	/** Removes all mappings and clears occupied values to 0 without shrinking the backing arrays. */
	@Override
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
	}

	/**
	 * Returns the value for a key.
	 * @return the stored value, or -1 if absent or equal to the current empty-cell sentinel
	 */
	@Override
	public int get(int key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
	}

	/**
	 * Equivalent to {@link #set(int, int)}.
	 * @param key key to insert or replace
	 * @param value value to store
	 * @return the previous value at the resolved cell; 0 for a newly inserted key
	 */
	@Override
	public int put(int key, int value){return set(key, value);}

	/**
	 * Copies all live mappings from another map into this map.
	 * @param map source map; its arrays are read but not retained
	 */
	public void putAll(IntHashMap map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				put(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * Associates a key with a value.
	 * @param key key to insert or replace
	 * @param value value to store
	 * @return the previous value at the resolved cell; 0 for a newly inserted key
	 */
	@Override
	public int set(int key, int value){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final int oldV=values[cell];
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
//		assert(get(key)==value);//TODO: slow
		return oldV;
	}

	/**
	 * Increments a key by 1, inserting it with value 1 when absent.
	 * @param key key to increment
	 * @return the new stored value
	 */
	@Override
	public int increment(int key){
		return increment(key, 1);
	}

	/**
	 * Increments a key by an arbitrary signed delta.
	 * Positive overflow clamps to {@link Integer#MAX_VALUE}; negative underflow
	 * wraps through the final int cast. An absent key starts from 0.
	 * @param key key to increment
	 * @param incr signed value to add
	 * @return the new stored value
	 */
	@Override
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
//		assert(get(key)==value);//TODO: slow
		return value;
	}

	/**
	 * Adds every live source value into this map by key.
	 * @param map source map; its arrays are read but not retained
	 */
	public void incrementAll(IntHashMap map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				increment(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * Merges another map by retaining the maximum value for each key.
	 * Keys absent from this map are inserted with their exact source values,
	 * including negative values. Self-merges are safe and leave the map unchanged.
	 * @param map source map; its arrays are read but not retained
	 */
	public void setToMax(IntHashMap map){
		for(int i=0; i<map.keys.length; i++){
			final int key=map.keys[i];
			if(key!=map.invalid){
				final int value=map.values[i];
				// [map/IntHashMap#001 FIXED] Probe by key, not by the source array index.
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
	 * Removes one mapping, then rehashes the following probe cluster to preserve reachability.
	 * @param key key to remove
	 * @return true if {@code key} was present and removed
	 */
	@Override
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
	 * Rehashes the contiguous probe cluster after a deletion.
	 * @param initial cell just emptied
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
	 * Moves one live cell to its preferred slot if a preceding deletion opened that slot.
	 * @param cell occupied cell to reconsider
	 * @return true if the key moved to a different array cell
	 */
	private boolean rehashCell(final int cell){
		final int key=keys[cell];
		final int value=values[cell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(key);
		if(cell==dest){return false;}
		assert(keys[dest]==invalid);
		keys[cell]=invalid;
		values[cell]=0;
		keys[dest]=key;
		values[dest]=value;

		return true;
	}

	/** Chooses a new negative empty-cell sentinel and rewrites old sentinels in the key array. */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		Vector.changeAll(keys, old, x);
	}

	/**
	 * Locates a key by probing from its hash cell until the key or an empty cell is found.
	 * @return the cell index containing {@code key}, or -1 when absent or when {@code key} is the current sentinel
	 */
	@Override
	public int findCell(final int key){
		if(key==invalid){return -1;}
		final int hash=key&MASK;
		//Small IDs already name their home cell; only larger masked keys need division.
		final int initial=(hash<modulus ? hash : hash%modulus);
		return Vector.findKeyScalar(keys, key, initial, invalid);
	}

	/**
	 * Finds the cell containing a valid key or the first empty cell in that key's probe sequence.
	 * @param key non-sentinel key to find or insert
	 * @return a cell whose key is {@code key} or {@link #invalid}
	 */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		final int initial=(int)((key&MASK)%modulus);
		return Vector.findKeyOrInvalid(keys, key, initial, invalid);
	}

	/** Rebuilds the table at a larger prime capacity after the load threshold is exceeded. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}

	/**
	 * Rebuilds the backing arrays near the first usable prime table capacity at least {@code size2}.
	 * Requested sizes near the int-array ceiling are capped to reserve the {@link #extra}
	 * cells that participate in probing after the logical prime modulus.
	 * @param size2 requested minimum table capacity before ceiling caps
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2);
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra);
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;

		final int size3=(int)(newPrime+extra);
		// Near the int-array ceiling, pack more entries before a resize that can no longer double safely.
		final float lf=(size3<0x60000000 ? loadFactor : Math.max(loadFactor, 0.80f));
		sizeLimit=(int)(modulus*lf);
		final int[] oldK=keys;
		final int[] oldV=values;
		keys=KillSwitch.allocInt1D(size3);
		values=KillSwitch.allocInt1D(size3);
		Arrays.fill(keys, invalid);

//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i], v=oldV[i];
			if(k!=invalid){
				set(k, v);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies all live keys into a dense new array.
	 * @return owned snapshot array in raw table iteration order
	 */
	@Override
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

	/** Returns the borrowed raw key array; entries equal to {@link #invalid()} are empty cells. */
	@Override
	public int[] keys(){return keys;}

	/** Returns the borrowed raw value array, index-aligned with {@link #keys()}. */
	@Override
	public int[] values(){return values;}

	/** Returns the current empty-cell sentinel for the raw key array. */
	@Override
	public int invalid(){return invalid;}

	/** Returns the number of live mappings. */
	@Override
	public int size(){return size;}

	/** Returns true when there are no live mappings. */
	@Override
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Key table; cells equal to {@link #invalid} are empty. */
	private int[] keys;
	/** Value table aligned with {@link #keys}; empty cells are kept at 0. */
	private int[] values;
	/** Number of live mappings. */
	private int size=0;
	/** Negative empty-cell sentinel for the keys array. */
	private int invalid;
	/** Logical prime modulus; cells between this value and {@code keys.length} are usable overflow slots. */
	private int modulus;
	/** Maximum live entries before the next resize. */
	private int sizeLimit;
	/** Effective load factor used to derive {@link #sizeLimit}. */
	private final float loadFactor;

	/** Deterministic source for negative empty-cell sentinels. */
	private static final Random randy=shared.Shared.random(1);

}
