package map;

import java.util.Arrays;
import java.util.HashMap;

import shared.KillSwitch;
import shared.Random;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;

/**
 * High-capacity int-to-int map that mixes int keys into a power-of-two table.
 * Unlike {@link IntHashMap2}, this variant can grow to the JVM's safe int-array
 * ceiling. Empty cells hold a rerollable negative sentinel, values in empty cells
 * are zero, get() returns -1 on misses, and put()/set() return zero on fresh
 * insertions. Positive increment overflow is capped at {@link Integer#MAX_VALUE};
 * negative underflow wraps through signed int truncation. Instances are unsynchronized.
 *
 * @author Isla
 * @date November 2, 2025
 */
public final class IntHashMap3 implements IntHashMapInterface{

	/** Serialization version. */
	private static final long serialVersionUID=1L;

	/**
	 * Runs destructive timing and maximum-capacity comparisons.
	 * @param args Optional element count, repeat count, and whether to remove loaded keys
	 */
	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;
		boolean remove=args.length>2 ? Boolean.parseBoolean(args[2]) : true;

		// Generate random keys with collisions
		System.err.println("Generating "+size+" random keys...");
		Random randy=shared.Shared.random(12345);
		int[] keys=new int[size];
		for(int i=0; i<size; i++){
			keys[i]=(int)Math.sqrt(randy.nextLong()&Long.MAX_VALUE);
		}
		System.err.println("Keys generated.");

		bench(keys, repeats, remove);
	}

	/**
	 * Times repeated load/lookup/remove cycles, then runs an unconditional legacy
	 * capacity stress even when the requested timing workload is tiny.
	 * @param keys Randomized keys to insert and query
	 * @param repeats Positive number of small timing cycles per implementation
	 * @param remove Whether timing cycles remove all inserted keys after lookup
	 */
	private static void bench(int[] keys, int repeats, boolean remove){
		final int size=keys.length;
		System.gc();
		Timer t=new Timer();
		{
			System.err.println("\n*** Warmup ***");
			Shared.printMemory();
			t.start();
			HashMap<Integer, Integer> map0=new HashMap<Integer, Integer>();
			IntHashMap map1=new IntHashMap();
			IntHashMap2 map2=new IntHashMap2();
			IntHashMap3 map3=new IntHashMap3();
			long sum=0;
			int warmupCycles=32000;
			{
				for(int i=0; i<warmupCycles; i++){
					map0.put(i, 2*i);
					map1.put(i, 2*i);
					map2.put(i, 2*i);
					map3.put(i, 2*i);
				}
				for(int i=0; i<warmupCycles; i++){
					sum+=map0.get(i);
					sum+=map1.get(i);
					sum+=map2.get(i);
					sum+=map3.get(i);
				}
				if(remove){
					for(int i=0; i<warmupCycles; i++){
						map0.remove(i);
						map1.remove(i);
						map2.remove(i);
						map3.remove(i);
					}
				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+map0.size()+", sum="+sum);
			Shared.printMemory();
			map0=null;
			map1=null;
			map2=null;
			map3=null;
			System.gc();
		}

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
			System.err.println("\n*** IntHashMap3 ***");
			Shared.printMemory();
			t.start();
			IntHashMap3 map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new IntHashMap3();
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

		// Legacy dev capacity stress: this runs unconditionally after the small timing loops,
		// so main() is not a smoke test even with tiny args.
		keys=new int[Shared.MAX_ARRAY_LEN];
		for(int i=0; i<keys.length; i++){keys[i]=i;}
		Random randy=Shared.threadLocalRandom(1);
		for(int i=0; i<keys.length; i++){
			int loc=randy.nextInt(keys.length);
			int temp=keys[i];
			keys[i]=keys[loc];
			keys[loc]=temp;
		}
		{
			long max=0;
			t.start();
			try{
				IntHashMap map=new IntHashMap();
				//for(long i=Integer.MAX_VALUE; i>=0; i--){
				for(int i : keys){
					map.put((int)i, (int)(i*2));
					max=map.size();
				}
			}catch(Throwable e){
				e.printStackTrace();
			}
			System.err.println("IntHashMap: capacity="+max);
			t.stopAndPrint();
		}
		{
			long max=0;
			t.start();
			try{
				IntHashMap2 map=new IntHashMap2();
				//for(long i=Integer.MAX_VALUE; i>=0; i--){
				for(int i : keys){
					map.put((int)i, (int)(i*2));
					max=map.size();
				}
			}catch(Throwable e){
				e.printStackTrace();
			}
			System.err.println("IntHashMap2: capacity="+max);
			t.stopAndPrint();
		}

		{
			long max=0;
			t.start();
			try{
				IntHashMap3 map=new IntHashMap3();
				//for(long i=Integer.MAX_VALUE; i>=0; i--){
				for(int i : keys){
					map.put((int)i, (int)(i*2));
					max=map.size();
				}
			}catch(Throwable e){
				e.printStackTrace();
			}
			System.err.println("IntHashMap3: capacity="+max);
			t.stopAndPrint();
		}

		{
			long max=0;
			t.start();
			try{
				HashMap<Integer, Integer> map=new HashMap<Integer, Integer>();
				//for(long i=Integer.MAX_VALUE; i>=0; i--){
				for(int i : keys){
					map.put((int)i, (int)(i*2));
					max=map.size();
				}
			}catch(Throwable e){
				e.printStackTrace();
			}
			System.err.println("HashMap: capacity="+max);
			t.stopAndPrint();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a map with the default initial size and load factor. */
	public IntHashMap3(){
		this(256);
	}

	/**
	 * Creates a map with default load factor 0.7.
	 * @param initialSize Positive requested size, rounded down by resize(long) [map/IntHashMap3#003]
	 */
	public IntHashMap3(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Creates a map with caller-selected initial size and load factor.
	 * The load factor is clamped to [0.25,0.90].
	 * @param initialSize Positive requested size, rounded down by resize(long) [map/IntHashMap3#003]
	 * @param loadFactor_ Requested load in (0,1), clamped to [0.25,0.90] after validation
	 */
	public IntHashMap3(int initialSize, double loadFactor_){
		invalid=randy.nextInt()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, (float)loadFactor_, 0.90f);
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
	 * Stores a key/value pair.
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous value, or 0 if key was not present (empty cells hold 0;
	 *         use contains() to disambiguate a stored 0). NOT -1. [map/IntHashMap3#002]
	 */
	public int put(int key, int value){
		return set(key, value);
	}

	/**
	 * Copies all live mappings from another map, overwriting matching keys here.
	 * @param map Source map to read without modifying
	 */
	public void putAll(IntHashMap3 map){
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
	 *         use contains() to disambiguate a stored 0). NOT -1. [map/IntHashMap3#002]
	 */
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
	public void incrementAll(IntHashMap3 map){
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
	public void setToMax(IntHashMap3 map){
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
	private int findCell(final int key){
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
		return Vector.findKeyOrInvalidScalar(keys, key, initial, invalid);
	}

	/** Requests doubled storage; resize(long) rounds down, caps, and rehashes live mappings. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L);
	}

	/**
	 * Replaces raw storage with a rounded-down power-of-two size and rehashes live mappings.
	 * The physical array is capped at {@link Shared#SAFE_ARRAY_LEN}.
	 * @param size2 Requested size, rounded down to logical capacity and capped physically
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;

		final long old=(keys==null ? 0 : keys.length);
		long size3=Long.highestOneBit(size2);//rounds DOWN to a power of 2 (growth path strips the +extra -> clean doubling)
		//RESOLVED non-bug [map/IntHashMap3#001] - mask is taken from the UNCAPPED pow2 here (capped to SAFE_ARRAY_LEN
		//next line), so at the 2^31 tier mask=(int)(2^31-1)=2147483647 exceeds the array maxIndex. This is SAFE, not
		//a bug: findCell uses initial=Tools.hash32plus(key)&mask, and hash32plus is CAPPED below 0x7FFFF800 (< array
		//length; see hash32plus javadoc). Since hash<0x7FFFF800, hash&mask==hash is ALWAYS a valid in-range index
		//regardless of the over-wide mask. The HASH (not the mask) guarantees the range. Verdict reached with Brian;
		//escalation RETRACTED. (My earlier "output-neutral" note over-thought it: the out-of-range hash&mask it
		//described can never occur with a capped hash - hash itself is always < array length.)
		mask=(int)(size3-1);
		size3=Math.min(size3+extra, Shared.SAFE_ARRAY_LEN);
		if(size3<=old || size3>Shared.SAFE_ARRAY_LEN){//note: size3>SAFE_ARRAY_LEN is dead (Math.min just capped it)
			throw new RuntimeException("Map hit capacity at "+size+":"
				+"\nkeys.length="+keys.length+"\nsize2="+size2+"\nsize3="+size3);
		}
		final float loadFactor2=(size3<Shared.SAFE_ARRAY_LEN ? loadFactor : Tools.max(loadFactor, 0.85f));
		sizeLimit=(int)((size3-extra)*loadFactor2);

		final int[] oldK=keys;
		final int[] oldV=values;
		keys=KillSwitch.allocInt1D((int)size3);
		values=KillSwitch.allocInt1D((int)size3);
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
	/** Power-of-two mask applied to capped mixed key hashes. */
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
