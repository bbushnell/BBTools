package structures;

import java.util.ArrayList;
import java.util.HashMap;

import shared.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntList;

/**
 * Base class for primitive int-to-int hash maps.
 * Concrete implementations choose the table sizing and probing strategy; shared
 * helpers iterate the raw arrays from {@link #keys()} and {@link #values()} and
 * skip slots whose key is equal to {@link #invalid()}.
 *
 * @author Brian Bushnell
 * @date June 8, 2017
 */
public abstract class AbstractIntHashMap{
	
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
	
	/** Removes all mappings while retaining implementation-owned storage. */
	public abstract void clear();
	
	/**
	 * Tests whether the key is present.
	 * @param key Key to search for.
	 * @return true if the key is present and is not the current invalid sentinel.
	 */
	public final boolean contains(int key){
		return findCell(key)>=0;
	}
	
	/**
	 * Tests whether the key is present.
	 * @param key Key to search for.
	 * @return true if the key is present and is not the current invalid sentinel.
	 */
	public final boolean containsKey(int key){
		return findCell(key)>=0;
	}
	
	/**
	 * Returns the value mapped to a key.
	 * @param key Key to search for.
	 * @return Stored value, or the implementation's missing-value sentinel if absent.
	 */
	public abstract int get(int key);
	
	/**
	 * Associates a value with a key, overwriting any previous value.
	 *
	 * @param key Key to map.
	 * @param value Value to store.
	 * @return Previous value in the selected cell, or the implementation's
	 * default value for a newly inserted key.
	 */
	public abstract int put(int key, int value);
	
	/**
	 * Associates a value with a key, overwriting any previous value.
	 *
	 * @param key Key to map.
	 * @param value Value to store.
	 * @return Previous value in the selected cell, or the implementation's
	 * default value for a newly inserted key.
	 */
	public abstract int set(int key, int value);
	
	/**
	 * Increments a key's value by one.
	 * Absent keys are inserted as if their previous value were zero.
	 * @param key Key to increment.
	 * @return Updated value.
	 */
	public abstract int increment(int key);
	
	/**
	 * Increments a key's value by the specified amount.
	 * Absent keys are inserted as if their previous value were zero.
	 *
	 * @param key Key to increment.
	 * @param incr Increment amount.
	 * @return Updated value.
	 */
	public abstract int increment(int key, int incr);

	/**
	 * Adds every live value from another map into this map.
	 * @param map Source map to scan through its raw backing arrays.
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
	 * Removes a key from the map if present.
	 * @param key Key to remove.
	 * @return true if the key was present and removed.
	 */
	public abstract boolean remove(int key);
	
	/**
	 * Returns the number of live mappings.
	 * @return Number of live mappings.
	 */
	public final int cardinality(){return size();}
	/**
	 * Returns the number of live mappings.
	 * @return Number of live mappings.
	 */
	public abstract int size();
	/**
	 * Tests whether the map has no live mappings.
	 * @return true if empty.
	 */
	public abstract boolean isEmpty();
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Returns a string representation of this map in list view format.
	 * Shows key-value pairs as (key,value) tuples.
	 * @return String representation of the map
	 */
	@Override
	public final String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns a debug string containing (rawCell, key, value) triples.
	 * @return String containing live entries and their backing-array cells.
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
	 * Returns a string containing (key,value) pairs in backing-array order.
	 * @return String containing live entries.
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
	 * Copies live keys into a compact array.
	 * @return Newly allocated array containing exactly the live keys.
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
	 * Counts how many times each stored value occurs in this map.
	 * @return Histogram keyed by values from this map.
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
	 * Verifies raw-array lengths, empty-cell values, and live-key reachability.
	 * @return true if each live key can be found at its current cell and the
	 * live-cell count matches {@link #size()}.
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
				final int cell=findCell(key);
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
	 * Finds the backing-array cell containing a key.
	 * @param key Key to search for.
	 * @return Backing-array cell index, or -1 if absent.
	 */
	abstract int findCell(final int key);
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies live keys into a compact array.
	 * @return Newly allocated array containing exactly the live keys.
	 */
	public abstract int[] toArray();
	/**
	 * Returns the raw backing key array, including invalid slots.
	 * @return Implementation-owned key array.
	 */
	public abstract int[] keys();
	/**
	 * Returns the raw backing value array, including values for invalid slots.
	 * @return Implementation-owned value array.
	 */
	public abstract int[] values();
	/**
	 * Returns the sentinel key used for empty backing-array slots.
	 * @return Current invalid key sentinel.
	 */
	public abstract int invalid();
	
	static final int MASK=Integer.MAX_VALUE;
	static final int MINMASK=Integer.MIN_VALUE;
	
	static final int extra=10;
	
}
