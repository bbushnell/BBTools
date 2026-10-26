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
 * Int-to-object map that mixes int keys into a power-of-two open-addressed table.
 * Empty cells hold a rerollable negative sentinel in the key array; value cells
 * may store null, so get(), put(), and remove() returning null does not distinguish absence from a stored
 * null. Instances are unsynchronized. The borrowed key array encodes probe structure;
 * external key writes can break lookup and rehash invariants. The borrowed value
 * array is direct typed storage, so value writes change occupied cells.
 *
 * @author Isla
 * @date November 2, 2025
 *
 * @param <V> Value type
 */
public final class IntObjectMap<V> implements Serializable{

	/** Serialization version. */
	private static final long serialVersionUID=1L;

	/**
	 * Runs a destructive timing comparison against {@link HashMap}.
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

	/**
	 * Times repeated load/lookup cycles on the supplied keys and values.
	 * @param keys Randomized keys to insert and query
	 * @param values Values to assign by source index
	 * @param repeats Positive number of map construction cycles per implementation
	 */
	private static void bench(int[] keys, ArrayList<String> values, int repeats){
		final int size=keys.length;
		System.gc();
		Timer t=new Timer();

		{
			System.err.println("\n*** IntObjectMap<String> ***");
			Shared.printMemory();
			t.start();
			IntObjectMap<String> map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new IntObjectMap<String>(String.class);
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
			System.err.println("\n*** HashMap<Integer, String> ***");
			Shared.printMemory();
			t.start();
			HashMap<Integer, String> map=null;
			long sum=0;
			for(int r=0; r<repeats; r++){
				map=new HashMap<Integer, String>();
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
	 * Creates a map with default initial size and load factor.
	 * @param type Component class of the VALUES (e.g. String.class, byte[].class, long[].class) — used to allocate a
	 *        real V[] backing array via reflection so values() returns a genuine typed array (a plain Object[] backing
	 *        CCEs at the call site under a concrete V). Pass a class literal; it is a cheap cached constant. [map/IntObjectMap#002]
	 */
	public IntObjectMap(Class<V> type){
		this(256, type);
	}

	/**
	 * Creates a map with default load factor 0.7.
	 * @param initialSize Positive requested logical minimum size, rounded up to the next power of two
	 * @param type Component class of the values (see {@link #IntObjectMap(Class)}).
	 */
	public IntObjectMap(int initialSize, Class<V> type){
		this(initialSize, 0.7f, type);
	}

	/**
	 * Creates a map with caller-selected initial size and load factor.
	 * The load factor is clamped to [0.25,0.90].
	 * @param initialSize Positive requested logical minimum size, rounded up to the next power of two
	 * @param loadFactor_ Requested load in (0,1), clamped to [0.25,0.90] after validation
	 * @param type Component class of the values (see {@link #IntObjectMap(Class)}).
	 */
	public IntObjectMap(int initialSize, float loadFactor_, Class<V> type){
		invalid=randy.nextInt()|MINMASK;
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

	/** Removes all mappings, resets values to null, and preserves the allocated raw arrays. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, null);
		size=0;
	}

	/**
	 * Returns the value mapped from a key.
	 * @param key Key to look up
	 * @return Mapped value, or null if absent or mapped to null
	 */
	public V get(int key){
		int cell=findCell(key);
		return cell<0 ? null : values[cell];
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
	 * Stores a key/value pair.
	 * @param key Key to insert or overwrite
	 * @param value Value to store; may be null
	 * @return Previous value, or null if absent or previously mapped to null
	 */
	public V put(int key, V value){
		return set(key, value);
	}

	/**
	 * Copies all live mappings from another map, overwriting matching keys here.
	 * @param map Source map to read without modifying
	 */
	public void putAll(IntObjectMap<V> map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				put(map.keys[i], map.values[i]);
			}
		}
	}

	/**
	 * Stores a key/value pair, rerolling the sentinel if it collides with {@code key}.
	 * @param key Key to insert or overwrite
	 * @param value Value to store; may be null
	 * @return Previous value, or null if absent or previously mapped to null
	 */
	public V set(int key, V value){
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
	 * Removes a key and closes the following probe cluster.
	 * @param key Key to remove
	 * @return Previous value, or null if absent or previously mapped to null
	 */
	public V remove(int key){
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

	/** Rerolls the empty-cell sentinel and rewrites all cells holding the old sentinel. */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				keys[i]=invalid;
			}
		}
	}

	/**
	 * Locates a key by probing from {@code Tools.hash32plus(key) & mask}.
	 * @param key Key to search for
	 * @return Raw cell index if found, or -1 if absent
	 */
	private int findCell(final int key){
		if(key==invalid){return -1;}

		final int limit=keys.length;
		final int hash=Tools.hash32plus(key);
		final int initial=hash&mask;

		for(int cell=initial; cell<limit; cell++){
			final int x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}

	/**
	 * Locates a key or the empty cell that would receive it.
	 * @param key Non-sentinel key to search for
	 * @return Raw cell containing {@code key}, or a raw cell holding {@link #invalid}
	 */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";

		//CLEVER [verified]: no separate occupancy array - empty cells hold the negative `invalid` sentinel and keys are
		//compared directly (x==key || x==invalid), one int compare per cell. initial=hash&mask is in range for ANY
		//mask: hash=Tools.hash32plus(key) is capped in [0, 0x7FFFF800) < SAFE_ARRAY_LEN, so even at the top tier
		//(mask=Integer.MAX_VALUE after the 2x-growth fix) hash&mask == hash is a valid index. set() guards key==invalid
		//via resetInvalid() (re-rolls the sentinel + remaps old cells), so a real key can never alias an empty cell.
		final int limit=keys.length;
		final int hash=Tools.hash32plus(key);
		final int initial=hash&mask;

		for(int cell=initial; cell<limit; cell++){
			final int x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}

	/** Doubles the logical power-of-two table size and rehashes live mappings. */
	private final void resize(){
		assert(size>=sizeLimit);
		//grow 2x: double the LOGICAL pow2 capacity (keys.length-extra==pow2), not keys.length (=pow2+extra).
		//keys.length*2L overshot the power-of-2 boundary -> rounded up to 4*pow2 (same family bug as the pow2 maps,
		//anchor [map/IntLongHashMap2#002]); fixed to true 2x growth. Output unchanged; ~2x less memory.
		resize((keys.length-extra)*2L);
	}

	/**
	 * Replaces raw storage with a rounded-up logical size and rehashes live mappings.
	 * Near-ceiling requests cap physical storage at {@link Shared#SAFE_ARRAY_LEN} and
	 * use terminal load factor 0.85.
	 * @param size2 Requested logical minimum size before power-of-two rounding and physical cap
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

		final int[] oldK=keys;
		@SuppressWarnings("unchecked")
		final V[] oldV=values;
		keys=KillSwitch.allocInt1D((int)size3);
		//Option A [map/IntObjectMap#002]: allocate a REAL V[] via reflection from the stored component type, so the
		//backing array's runtime type is genuinely V[] (e.g. String[], byte[][]) and values() returns a usable typed
		//array instead of an Object[] that CCEs at the call site. `type` is a cheap cached class literal.
		@SuppressWarnings("unchecked")
		V[] temp=(V[])Array.newInstance(type, (int)size3);
		values=temp;
		Arrays.fill(keys, invalid);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i];
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
	 * Borrows the live internal value array.
	 * Stored null values are allowed, so inspect parallel keys() entries against invalid() for occupancy.
	 * This is a GENUINE V[] (e.g. String[], byte[][]) because the backing is allocated via Array.newInstance(type,...)
	 * from the component class passed to the constructor [map/IntObjectMap#002] — so it can be used directly (assigned
	 * to V[], iterated, .length) without the ClassCastException a plain Object[] backing would cause at the call site.
	 * WARNING: live array. Use with caution.
	 * @return Internal value array (a real V[])
	 */
	public V[] values(){return values;}

	/**
	 * Returns the current empty-cell sentinel.
	 * @return Current negative empty-cell sentinel
	 */
	public int invalid(){return invalid;}

	/**
	 * Returns the number of mappings.
	 * @return Number of live mappings
	 */
	public int size(){return size;}

	/**
	 * Tests whether the map has no mappings.
	 * @return true when the map has no live mappings
	 */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Raw int key table; empty cells contain {@link #invalid}. */
	private int[] keys;
	/** Raw typed values, index-aligned with {@link #keys}; live cells may contain null. */
	private V[] values;
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
	/** Component class of V, for allocating a real V[] backing via reflection (Array.newInstance). [map/IntObjectMap#002] */
	private final Class<V> type;

	/** Physical overflow cells above the logical power-of-two mask range. */
	private static final int extra=10;
	/** Mask that forces sentinel candidates negative. */
	private static final int MINMASK=Integer.MIN_VALUE;

	/** Reproducible source for negative empty-cell sentinel candidates. */
	private static final Random randy=shared.Shared.random(1);

}
