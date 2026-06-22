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
 * memory-efficient replacement for HashSet<T> by:
 * - Using arrays instead of Entry objects
 * - Avoiding per-entry object overhead
 * - Better cache locality with linear probing
 * Thread-safety: Not thread-safe. External synchronization required for concurrent access.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date January 26, 2026
 * @param <T> Key type - must properly implement hashCode() and equals()
 */
public final class ObjectSet<T> implements Serializable {

	private static final long serialVersionUID = 1L;

	public static void main(String[] args){
		int size=args.length>0 ? Integer.parseInt(args[0]) : 1000000;
		int repeats=args.length>1 ? Integer.parseInt(args[1]) : 1;
		ArrayList<String> list=randomStrings(size, 1, 50);
		test(list);
		bench(list, repeats);
	}

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
	 * Tests correctness by comparing ObjectSet against HashSet.
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
	 * Benchmarks performance against HashSet<Object>.
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
	 * @param type Component class of the keys (e.g. String.class, byte[].class) — used to allocate a real T[] backing
	 *        array via reflection so keys() returns a genuine typed array (a plain Object[] backing CCEs at the call
	 *        site under a concrete T). Pass a class literal; it is a cheap cached constant. [map/ObjectSet#002]
	 */
	public ObjectSet(Class<T> type){
		this(256, 0.7f, type);
	}

	public ObjectSet(int initialSize, Class<T> type){
		this(initialSize, 0.7f, type);
	}

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
	 * Removes all entries from the set.
	 */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, null);
		Arrays.fill(hashes, 0);
		size=0;
	}

	/**
	 * Checks if the set contains the given key.
	 * @param key Key to check
	 * @return true if key is present
	 */
	public boolean contains(T key){
		return findCell(key)>=0;
	}

	/**
	 * Adds the specified key to the set.
	 * @param key Key to add
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
	 * Checks if the set contains the given key.
	 * @param key Key to check
	 * @return true if key is present
	 */
	public T get(T key){
		int cell=findCell(key);
		return cell<0 ? null : keys[cell];
	}

	/**
	 * Adds the specified key to the set, replacing if present.
	 * @param key Key to add
	 * @return Old key
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
	 * @param key Key to add
	 * @return Old key
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
	 * @param set Source set to copy from
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
	 * @param key Key to remove
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

	private int findCell(final T key){
		if(key==null){return -1;}

		final int limit=keys.length;
		final int hash=Tools.hash32plus(key.hashCode());
		final int initial=hash & mask;

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

	private int findCellOrEmpty(final T key, final int hash){
		assert(key!=null) : "Null keys not supported";

		//CLEVER [verified]: hashes[] caches each key's hash so the probe compares cheap ints first
		//(hashes[cell]==hash) and only calls equals() on a hash match. initial=hash & mask is in range for ANY mask:
		//hash=Tools.hash32plus(...) is CAPPED in [0, 0x7FFFF800) < SAFE_ARRAY_LEN, so even at the top tier
		//(mask=Integer.MAX_VALUE after the 2x-growth fix) hash & mask == hash is a valid index (see hash32plus javadoc).
		final int limit=keys.length;
		final int initial=hash & mask;

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

	private final void resize(){
		assert(size>=sizeLimit);
		//grow 2x: double the LOGICAL pow2 capacity (keys.length-extra==pow2), not keys.length (=pow2+extra).
		//keys.length*2L overshot the power-of-2 boundary -> rounded up to 4*pow2 (same family bug as the pow2 maps,
		//anchor [map/IntLongHashMap2#002]); fixed to true 2x growth. Output unchanged; ~2x less memory.
		resize((keys.length-extra)*2L);
	}

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
	 * Returns the internal key array (the live backing array; null entries mark empty cells).
	 * This is a GENUINE T[] (e.g. String[], byte[][]) because the backing is allocated via Array.newInstance(type,...)
	 * from the component class passed to the constructor [map/ObjectSet#002] — so it can be used directly (assigned to
	 * T[], iterated, .length) without the ClassCastException a plain Object[] backing would cause at the call site.
	 * WARNING: live array, contains nulls for empty cells. Use with caution.
	 * @return Internal key array (a real T[])
	 */
	public T[] keys(){return keys;}

	public int size(){return size;}

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
	/** Load factor (fraction of capacity before resize) */
	private final float loadFactor;
	/** Component class of T, for allocating a real T[] backing via reflection (Array.newInstance). [map/ObjectSet#002] */
	private final Class<T> type;

	/** Extra space beyond power-of-2 size to reduce wrap-around collisions */
	private static final int extra=10;

}