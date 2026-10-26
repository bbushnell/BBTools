package map;

import java.util.ArrayList;
import java.util.HashMap;

import shared.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntList;

/**
 * Shared base for mutable int-to-int maps backed by raw key/value tables.
 * Subclasses own storage, collision resolution, absent-value defaults, overflow policy,
 * and sentinel remapping; this class supplies membership wrappers, dense snapshots,
 * string views, whole-map aggregation, histogramming, and structural verification.
 * Instances are unsynchronized. The raw key array returned by {@link #keys()}
 * encodes probe structure; external key writes can break lookup and rehash invariants.
 * The raw value array returned by {@link #values()} is direct storage indexed by
 * the key table, so value writes change mappings at occupied key cells.
 *
 * @author Brian Bushnell
 * @date June 8, 2017
 */
public abstract class AbstractIntHashMap implements IntHashMapInterface{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs an empty abstract map base. */
	public AbstractIntHashMap(){}

	/*--------------------------------------------------------------*/
	/*----------------           Testing            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Exercises an initially empty int map against {@link HashMap} and then runs a
	 * destructive set/lookup/remove microbenchmark.
	 * @param set Non-null empty map to refill, drain, and finally leave populated
	 *        with the last benchmark cycle's keys and values
	 */
	public static final void test(AbstractIntHashMap set){
		Random randy2=Shared.threadLocalRandom();
		HashMap<Integer, Integer> set2=new HashMap<Integer, Integer>(20, 0.7f);
		ArrayList<Integer> klist=new ArrayList<Integer>();
		ArrayList<Integer> klist2=new ArrayList<Integer>();
		ArrayList<Integer> vlist=new ArrayList<Integer>();
		ArrayList<Integer> vlist2=new ArrayList<Integer>();

		for(int i=0; i<1000; i++){
			assert(!set.contains(i));
			assert(!set2.containsKey(i));
			klist.add(Integer.valueOf(i));
			vlist.add(Integer.valueOf(i*2+7));
		}
		for(int i=0; i<1000; i++){
			int r=randy2.nextInt();
			klist2.add(r);
			vlist2.add(randy2.nextInt()&Integer.MAX_VALUE);
		}


		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			set.put(k, v);
			set2.put(k, v);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			assert(set.size()==set2.size());
		}
		assert(!set.isEmpty());
		assert(!set2.isEmpty());
		assert(set.size()==set2.size());

		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
		}

		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			set.increment(k);
			assert(set.get(k)==v+1);
			set.increment(k, -1);
			assert(set.get(k)==v);
		}

		for(int i=0; i<klist.size(); i++){
			int k=klist.get(i), v=vlist.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			set.remove(k);
			set2.remove(k);
			assert(!set.containsKey(k));
			assert(!set2.containsKey(k));
			assert(set.size()==set2.size());
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		assert(set.size()==set2.size());


		for(int i=0; i<klist2.size(); i++){
			int k=klist2.get(i), v=vlist2.get(i);
			set.put(k, v);
			set2.put(k, v);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			assert(set.size()==set2.size());
		}
		assert(!set.isEmpty());
		assert(!set2.isEmpty());
		assert(set.size()==set2.size());

		for(int i=0; i<klist2.size(); i++){
			int k=klist2.get(i), v=vlist2.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
		}

		for(int i=0; i<klist2.size(); i++){
			int k=klist2.get(i), v=vlist2.get(i);
			assert(set.get(k)==v);
			assert(set2.get(k)==v);
			set.remove(k);
			set2.remove(k);
			assert(!set.containsKey(k));
			assert(!set2.containsKey(k));
			assert(set.size()==set2.size());
		}
		assert(set.isEmpty());
		assert(set2.isEmpty());
		assert(set.size()==set2.size());


		int count=4000000;
		int runs=32;
		IntList kil=new IntList(count);
		IntList vil=new IntList(count);
		for(int i=0; i<count; i++){
			kil.add(randy2.nextInt());
			vil.add(randy2.nextInt()&Integer.MAX_VALUE);
		}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("IntHashMap:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(int x : ll.array){
//					set.add(x);
//				}
				final int[] kila=kil.array;
				final int[] vila=vil.array;
				for(int z=0; z<count; z++){
					final int key=kila[z];
					final int value=vila[z];
					set.set(key, value);
					set.contains(key);
					set.remove(key);
					set.set(key, value);
				}
//				for(int x : ll.array){
//					set.remove(x);
//				}
//				set.clear();
//				assert(set.isEmpty());
//				System.err.println("Finished run "+i);
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();

//			System.err.println("HashMap:");
//			t.start();
//			for(int i=0; i<runs; i++){
//				for(int x : ll.array){
//					set2.add(x);
//				}
//				for(int x : ll.array){
//					set2.remove(x);
//				}
//				assert(set2.isEmpty());
////				System.err.println("Finished run "+i);
//			}
//			t.stop();
//			System.err.println(t);
//			Shared.printMemory();
		}
		t.stop();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Removes every mapping while preserving the map's allocated storage. */
	public abstract void clear();

	/**
	 * Tests whether the key is currently mapped.
	 * @param key Key to look up
	 * @return true if {@code key} is live in this map
	 */
	public final boolean contains(int key){
		return findCell(key)>=0;
	}

	/**
	 * Tests whether the key is currently mapped.
	 * @param key Key to look up
	 * @return true if {@code key} is live in this map
	 */
	public final boolean containsKey(int key){
		return findCell(key)>=0;
	}

	/**
	 * Returns the value mapped from a key.
	 * @param key Key to look up
	 * @return Mapped value, or the subclass-specific absent default
	 */
	public abstract int get(int key);

	/**
	 * Sets a key to a value, inserting the key if absent and replacing its old value if present.
	 *
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous mapped value, or the subclass-specific absent default
	 */
	public abstract int put(int key, int value);

	/**
	 * Sets a key to a value, inserting the key if absent and replacing its old value if present.
	 *
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous mapped value, or the subclass-specific absent default
	 */
	public abstract int set(int key, int value);

	/**
	 * Adds one to a key's value, inserting the key when absent.
	 * @param key Key to insert or update
	 * @return New value after the subclass-specific overflow policy
	 */
	public abstract int increment(int key);

	/**
	 * Adds a signed delta to a key's value, inserting the key when absent.
	 *
	 * @param key Key to insert or update
	 * @param incr Signed increment to apply
	 * @return New value after the subclass-specific overflow and underflow policy
	 */
	public abstract int increment(int key, int incr);

	/**
	 * Adds all live values from another int map into this map.
	 * Iterates the source's borrowed raw arrays and calls {@link #increment(int, int)}
	 * for every key not equal to the source's current sentinel. If {@code map} is
	 * this map, each value is doubled as its occupied cell is visited; distinct
	 * sources must remain structurally stable while their raw arrays are scanned.
	 * @param map Source map to scan for increments; may be this map
	 */
	public final void incrementAll(AbstractIntHashMap map){
		final int[] keys=map.keys();
		final int[] values=map.values();
		final int invalid=map.invalid();
		for(int i=0; i<keys.length; i++){
			int a=keys[i], b=values[i];
			if(a!=invalid){
				increment(a, b);
			}
		}
	}

	/**
	 * Removes a key if present and rehashes any dependent probe cluster cells.
	 * @param key Key to remove
	 * @return true if a live mapping was removed
	 */
	public abstract boolean remove(int key);

	/**
	 * Returns the number of live mappings.
	 * @return Number of live mappings
	 */
	public final int cardinality(){return size();}
	/**
	 * Returns the number of live mappings.
	 * @return Number of live mappings
	 */
	public abstract int size();
	/**
	 * Tests whether the map has no live mappings.
	 * @return true when the map has no live mappings
	 */
	public abstract boolean isEmpty();

	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns {@link #toStringListView()}.
	 * @return Key/value tuple list in raw table order
	 */
	@Override
	public final String toString(){
		return toStringListView();
	}

	/**
	 * Formats occupied cells in raw table order as triples containing the cell index, key, and value.
	 * @return List-like string containing {@code (cell, key, value)} triples
	 */
	public final String toStringSetView(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys(), values=values();
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys().length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+"("+i+", "+keys[i]+", "+values[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Formats mappings in raw table order as key/value tuples.
	 * @return List-like string containing {@code (key,value)} tuples
	 */
	public final String toStringListView(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys(), values=values();
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys().length; i++){
			if(keys[i]!=invalid){
				sb.append(comma);
				sb.append('(');
				sb.append(keys[i]);
				sb.append(',');
				sb.append(values[i]);
				sb.append(')');
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Copies all live keys into a dense new array.
	 * @return Owned key snapshot in raw table iteration order
	 */
	public final int[] toKeyArray(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys();
		int[] x=new int[size];
		int i=0;
		for(int v : keys){
			if(v!=invalid){
				x[i]=v;
				i++;
			}
		}
		return x;
	}

	/**
	 * Counts how often each stored value occurs.
	 * @return New binary map from stored value to the number of live cells with that value
	 */
	public final IntHashMapBinary toCountHistogram(){
		IntHashMapBinary counts=new IntHashMapBinary(Tools.mid(1, size(), 64));
		final int[] keys=keys();
		final int[] values=values();
		final int invalid=invalid();
		for(int i=0; i<keys.length; i++){
			int a=keys[i];
			int b=values[i];
			if(a!=invalid){
				counts.increment(b);
			}
		}
		return counts;
	}

	/**
	 * Checks raw array shape, empty-cell clearing, lookup placement, and live size.
	 * @return true if the current table satisfies the expected probe invariants
	 */
	public final boolean verify(){
		final int size=size(), invalid=invalid();
		final int[] keys=keys(), values=values();
		int numValues=0;
		int numFound=0;
		if(keys.length!=values.length){return false;}
		for(int i=0; i<keys.length; i++){
			final int key=keys[i];
			if(key==invalid){
				if(values[i]!=0){return false;}
			}else{
				numValues++;
				final int cell=findCell(key);//was findCell(i) - latent bug: findCell takes a KEY, so verify must look up keys[i], matching every Long* verify() sibling (e.g. LongHashMap:376). verify() is dead (no caller) -> output-neutral fix. [map/AbstractIntHashMap#001 FIXED]
				if(i==cell){
					numFound++;
				}else{
					return false;
				}
			}
		}
		return numValues==numFound && numValues==size;
	}

	/**
	 * Finds the raw cell holding a live key.
	 * @param key Key to locate
	 * @return Occupied cell index, or -1 if absent
	 */
	abstract int findCell(final int key);

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies all live keys into a dense new array.
	 * @return Owned key snapshot in raw table iteration order
	 */
	public abstract int[] toArray();

	/**
	 * Borrows the raw key array.
	 * @return Borrowed raw key array with {@link #invalid()} in empty cells
	 */
	public abstract int[] keys();

	/**
	 * Borrows the raw value array.
	 * @return Borrowed raw value array, index-aligned with {@link #keys()}
	 */
	public abstract int[] values();

	/**
	 * Returns the current empty-cell sentinel.
	 * @return Current empty-cell sentinel for {@link #keys()}
	 */
	public abstract int invalid();

	/** Mask to clear the sign bit from a raw int hash or key. */
	static final int MASK=Integer.MAX_VALUE;
	/** Mask to isolate the sign bit from a raw int key. */
	static final int MINMASK=Integer.MIN_VALUE;

	/** Physical overflow cells above the logical modulus for int-key prime tables. */
	static final int extra=10;

}
