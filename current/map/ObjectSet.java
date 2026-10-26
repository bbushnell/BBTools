package map;

import java.io.Serializable;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import shared.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Hash set with Object keys.
 * Uses open addressing with linear probing for collision resolution.
 * Caches hash codes to avoid expensive equals() calls during probing.
 * Uses power-of-2 sizing for fast modulo via bitwise AND.
 *
 * Parallel arrays retain keys without allocating per-entry objects. Null keys
 * cannot be stored; null lookups/removals behave as absent keys. Equality and
 * hashes must stay stable while a key is stored. Not thread-safe; synchronize
 * accesses that overlap a mutation. Equal keys can have different identities:
 * add/addIfAbsent preserve the existing representative, while addAndReturnOld replaces it.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date January 26, 2026
 * @param <T> Key type - must properly implement hashCode() and equals()
 */
public final class ObjectSet<T> implements Serializable{

	/** Serialization version for this set layout. */
	private static final long serialVersionUID=1L;

	/** Runs a legacy presence comparison and allocation/timing benchmark on random strings.
	 * @param args Optional positional nonnegative key count and positive repetition count */
	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;
		ArrayList<String> list=randomStrings(size, 1, 50);
		test(list);
		bench(list, repeats);
	}

	/** Generates reproducible uppercase A-Z strings without ensuring uniqueness.
	 * @param size Nonnegative number of strings
	 * @param minLen Nonnegative minimum string length
	 * @param maxLen Inclusive maximum length, at least minLen; range must fit an int
	 * @return Newly allocated list and strings */
	private static ArrayList<String> randomStrings(int size, int minLen, int maxLen){
		Shared.printMemory();
		ArrayList<String> strings=new ArrayList<String>(size);
		Random randy=shared.Shared.random(12345);
		int range=maxLen-minLen+1;
		for(int i=0; i<size; i++){
			int len=randy.nextInt(range)+minLen;
			StringBuilder sb=new StringBuilder(len);
			for(int j=0; j<len; j++){
				sb.append((char)('A'+randy.nextInt(26)));
			}
			strings.add(sb.toString());
		}
		Shared.printMemory();
		System.err.println("Generated "+size+" random strings");
		return strings;
	}

	/**
	 * Compares insertion/presence with HashSet, reporting to stderr and exiting with
	 * status one on a presence mismatch. The size comparison requires assertions.
	 * @param list Non-null list of non-null keys with stable equality/hash codes; not modified
	 */
	private static void test(ArrayList<? extends Object> list){
		System.err.println("\n*** Testing Correctness ***");

		// Build HashSet
		HashSet<Object> hashSet=new HashSet<Object>();
		for(int i=0; i<list.size(); i++){
			hashSet.add(list.get(i));
		}

		// Build ObjectSet
		ObjectSet<Object> objectSet=new ObjectSet<Object>(Object.class);
		for(int i=0; i<list.size(); i++){
			objectSet.add(list.get(i));
		}

		System.err.println("HashSet size: "+hashSet.size());
		System.err.println("ObjectSet size: "+objectSet.size());
		assert(hashSet.size()==objectSet.size()) : "Size mismatch!";

		// Validate all values match
		int errors=0;
		for(int i=0; i<list.size(); i++){
			Object key=list.get(i);
			boolean h=hashSet.contains(key);
			boolean o=objectSet.contains(key);
			if(h!=o){
				System.err.println("ERROR at index "+i+": key='"+key+"', HashSet="+h+", ObjectSet="+o);
				errors++;
				if(errors>10){break;}
			}
		}

		if(errors==0){
			System.err.println("*** PASS: All presence checks match! ***");
		}else{
			System.err.println("*** FAIL: "+errors+" mismatches found! ***");
			System.exit(1);
		}
	}

	/**
	 * Benchmarks performance against {@code HashSet<Object>}.
	 * Reconstructs each set per repetition, inserts/looks up/removes all supplied keys,
	 * prints timings and requests GC. This is an ad hoc benchmark, not a throughput oracle.
	 * @param list Non-null list of non-null stable keys; not modified
	 * @param repeats Positive number of cycles per set implementation
	 */
	private static void bench(ArrayList<? extends Object> list, int repeats){
		System.gc();
		Timer t=new Timer();

		{
			System.err.println("\n*** ObjectSet<Object> ***");
			Shared.printMemory();
			t.start();
			ObjectSet<Object> set=null;
			for(int r=0; r<repeats; r++){
				set=new ObjectSet<Object>(Object.class);
				for(int i=0; i<list.size(); i++){
					set.add(list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					boolean x=set.contains(list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					set.remove(list.get(i));
				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+set.size());
			Shared.printMemory();
			set=null;
			System.gc();
		}

		{
			System.err.println("\n*** HashSet<Object> ***");
			Shared.printMemory();
			t.start();
			HashSet<Object> set=null;
			for(int r=0; r<repeats; r++){
				set=new HashSet<Object>();
				for(int i=0; i<list.size(); i++){
					set.add(list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					boolean x=set.contains(list.get(i));
				}
				for(int i=0; i<list.size(); i++){
					set.remove(list.get(i));
				}
			}
			t.stop("Time: \t");
			System.gc();
			System.err.println("Size: "+set.size());
			Shared.printMemory();
			set=null;
			System.gc();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an empty set with default logical capacity 256 and load factor 0.7.
	 * @param type Non-null reference component class assignable from every stored key;
	 * used to allocate correctly typed arrays returned by keys() [map/ObjectSet#002]
	 */
	public ObjectSet(Class<T> type){
		this(256, 0.7f, type);
	}

	/** Allocates a typed set with requested initial capacity and default load factor 0.7.
	 * @param initialSize Positive minimum logical size, rounded up to a power of two
	 * @param type Non-null reference component class assignable from every stored key */
	public ObjectSet(int initialSize, Class<T> type){
		this(initialSize, 0.7f, type);
	}

	/** Allocates typed storage for a new empty set.
	 * @param initialSize Positive minimum logical size, rounded up to a power of two
	 * @param loadFactor_ Requested fraction in (0,1), clamped to [0.25,0.90]
	 * @param type Non-null reference component class assignable from every stored key */
	public ObjectSet(int initialSize, float loadFactor_, Class<T> type){
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
	 */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, null);
		Arrays.fill(hashes, 0);
		size=0;
	}

	/**
	 * Checks if the set contains the given key.
	 * @param key Lookup key; null is treated as absent
	 * @return true if key is present
	 */
	public boolean contains(T key){
		return findCell(key)>=0;
	}

	/**
	 * Adds a key if absent, preserving the stored representative when an equal key exists.
	 * @param key Non-null key compatible with the constructor's component type
	 * @return true if the set did not already contain the key
	 */
	public boolean add(T key){
		//three insert variants share the findCellOrEmpty probe: add (bool added), addAndReturnOld (returns prior key,
		//overwrites), addIfAbsent (returns prior key, does NOT overwrite). All store the cached hash beside the key.
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		if(keys[cell]!=null){return false;} // Already present

		keys[cell]=key;
		hashes[cell]=hash;
		size++;
		if(size>sizeLimit){resize();}
		return true;
	}

	/**
	 * Finds the stored representative equal to the query key.
	 * @param key Lookup key; null is treated as absent
	 * @return Borrowed stored representative, or null when absent; do not mutate its hash/equality state
	 */
	public T get(T key){
		int cell=findCell(key);
		return cell<0 ? null : keys[cell];
	}

	/**
	 * Adds a key, replacing an equal stored representative with this object reference.
	 * Existing-key replacement preserves size; a new key increases size by one.
	 * @param key Non-null key compatible with the constructor's component type
	 * @return Previous stored representative, or null when newly inserted
	 */
	public T addAndReturnOld(T key){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		T old=keys[cell];
		if(old!=null){// Already present
			keys[cell]=key;
			return old;
		}

		keys[cell]=key;
		hashes[cell]=hash;
		size++;
		if(size>sizeLimit){resize();}
		return null;
	}

	/**
	 * Adds the specified key to the set, only if not present.
	 * @param key Non-null key compatible with the constructor's component type
	 * @return Existing stored representative, or null when newly inserted
	 */
	public T addIfAbsent(T key){
		assert(key!=null) : "Null keys not supported";
		final int hash=Tools.hash32plus(key.hashCode());
		final int cell=findCellOrEmpty(key, hash);
		T old=keys[cell];
		if(old!=null){return old;}// Already present

		keys[cell]=key;
		hashes[cell]=hash;
		size++;
		if(size>sizeLimit){resize();}
		return null;
	}

	/**
	 * Copies all entries from another set into this set.
	 * Existing equal representatives are retained; added keys are shared, not cloned.
	 * A distinct source is preserved; self-copy leaves the set unchanged.
	 * @param set Non-null source stable during this operation, with compatible key types
	 */
	public void addAll(ObjectSet<T> set){
		for(int i=0; i<set.keys.length; i++){
			if(set.keys[i]!=null){
				add(set.keys[i]);
			}
		}
	}

	/**
	 * Removes the specified key.
	 * @param key Key to remove; null is treated as absent
	 * @return true if the set contained the key
	 */
	public boolean remove(T key){
		if(key==null){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(keys[cell].equals(key));

		keys[cell]=null;
		hashes[cell]=0;
		size--;

		rehashFrom(cell);
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Repairs the contiguous probe cluster after deletion, wrapping at the array end.
	 * @param initial Cell just emptied; subsequent live entries may need rehoming */
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=keys.length;
		for(int cell=initial+1; cell<limit; cell++){
			final T key=keys[cell];
			if(key==null){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final T key=keys[cell];
			if(key==null){return;}
			rehashCell(cell);
		}
	}

	/** Moves one entry earlier in its probe chain, preserving its hash and the set size.
	 * @param cell Occupied cell index
	 * @return True if the entry moved and its old cell was cleared */
	private boolean rehashCell(final int cell){
		final T key=keys[cell];
		final int hash=hashes[cell];
		assert(key!=null);
		final int dest=findCellOrEmpty(key, hash);
		if(cell==dest){return false;}
		assert(keys[dest]==null);
		keys[cell]=null;
		hashes[cell]=0;
		keys[dest]=key;
		hashes[dest]=hash;

		return true;
	}

	/** Looks up a key without mutating storage, stopping at the first empty cell.
	 * @param key Lookup key; null is treated as absent
	 * @return Cell index containing an equal key, or -1 when absent */
	private int findCell(final T key){
		if(key==null){return -1;}

		final int limit=keys.length;
		final int hash=Tools.hash32plus(key.hashCode());
		final int initial=hash&mask;

		for(int cell=initial; cell<limit; cell++){
			final T x=keys[cell];
			if(x==null){return -1;}
			if(hashes[cell]==hash && x.equals(key)){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final T x=keys[cell];
			if(x==null){return -1;}
			if(hashes[cell]==hash && x.equals(key)){return cell;}
		}
		return -1;
	}

	/** Finds an equal key or the first empty slot in its probe sequence.
	 * @param key Non-null key to search for
	 * @param hash Precomputed Tools.hash32plus(key.hashCode()), not the raw hashCode()
	 * @return Matching or empty cell index
	 * @throws RuntimeException If no matching or empty cell exists */
	private int findCellOrEmpty(final T key, final int hash){
		assert(key!=null) : "Null keys not supported";

		//Cached hashes avoid equals() calls on mismatches. Tools.hash32plus caps the initial
		//index below SAFE_ARRAY_LEN even at the terminal, physically capped table size.
		final int limit=keys.length;
		final int initial=hash&mask;

		for(int cell=initial; cell<limit; cell++){
			final T x=keys[cell];
			if(x==null || (hashes[cell]==hash && x.equals(key))){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final T x=keys[cell];
			if(x==null || (hashes[cell]==hash && x.equals(key))){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}

	/** Doubles logical capacity before rehashing all keys. */
	private final void resize(){
		assert(size>=sizeLimit);
		//Double logical capacity, excluding slack; doubling the full array would round
		//ordinary growth up to four times the logical capacity. [map/IntLongHashMap2#002]
		resize((keys.length-extra)*2L);
	}

	/** Rebuilds both arrays, preserving keys, cached hashes and size.
	 * Rounds the logical request up to a power of two, adds slack and caps physical
	 * length at Shared.SAFE_ARRAY_LEN. The terminal allocation uses load factor 0.85.
	 * Borrowed array references are replaced; key objects themselves are retained.
	 * @param size2 Positive requested logical size, strictly greater than the entry count
	 * @throws RuntimeException If the capped allocation cannot grow the existing table */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;

		long size3=Long.highestOneBit(size2);
		if(size3<size2){size3<<=1;}
		mask=(int)(size3-1);
		size3=Math.min(size3+extra, Shared.SAFE_ARRAY_LEN);
		if((keys!=null && size3<=keys.length) || size3>Shared.SAFE_ARRAY_LEN){
			throw new RuntimeException("Set hit capacity at "+size);
		}
		final float loadFactor2=(size3<Shared.SAFE_ARRAY_LEN ? loadFactor : 0.85f);
		sizeLimit=(int)((size3-extra)*loadFactor2);

		@SuppressWarnings("unchecked")
		final T[] oldK=keys;
		final int[] oldH=hashes;
		//Option A [map/ObjectSet#002]: allocate a REAL T[] via reflection from the stored component type, so the
		//backing array's runtime type is genuinely T[] (e.g. String[], byte[][]) and keys() returns a usable typed
		//array instead of an Object[] that CCEs at the call site. `type` is a cheap cached class literal.
		@SuppressWarnings("unchecked")
		T[] tempK=(T[])Array.newInstance(type, (int)size3);
		keys=tempK;
		hashes=KillSwitch.allocInt1D((int)size3);

		if(size<1){return;}

		//Re-home every entry into the new array: mask changed, so each key gets a fresh probe start. Place key+cached
		//hash directly (size is unchanged - we're moving keys, not adding - so no size++ and no resize re-trigger).
		for(int i=0; i<oldK.length; i++){
			final T k=oldK[i];
			if(k!=null){
				final int h=oldH[i];
				int cell=findCellOrEmpty(k, h);
				keys[cell]=k;
				hashes[cell]=h;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows the live typed key array; null entries mark empty cells.
	 * Do not modify its entries or retain it across resizing. Stored keys must retain
	 * stable hashes/equality even when callers never write directly to the array.
	 * @return Set-owned array with the constructor-supplied component type
	 */
	public T[] keys(){return keys;}

	/** Returns occupied-key count.
	 * @return Number of distinct stored keys */
	public int size(){return size;}

	/** Returns whether no keys are stored.
	 * @return True when size() is zero */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Array of keys (null for empty cells) */
	private T[] keys;
	/** Array of cached hash codes (parallel to keys array, 0 for empty cells) */
	private int[] hashes;
	/** Number of entries in the set */
	private int size=0;
	/** Bit mask for fast modulo (always power of 2 minus 1) */
	private int mask;
	/** Size threshold for triggering resize (capacity * loadFactor) */
	private int sizeLimit;
	/** Requested load clamped to [0.25,0.90]; terminal allocation uses 0.85 instead. */
	private final float loadFactor;
	/** Component class of T, for allocating a real T[] backing via reflection (Array.newInstance). [map/ObjectSet#002] */
	private final Class<T> type;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Extra physical probe cells beyond the power-of-two logical capacity. */
	private static final int extra=10;

}
