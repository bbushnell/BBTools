package structures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date June 8, 2017
 *
 */
public abstract class AbstractIntHashMap{
	
	/**
	 * Comprehensive test method that validates hash map implementation correctness.
	 * Tests basic operations (put, get, remove, contains), increment operations,
	 * size tracking, and performance benchmarks against standard HashMap.
	 * Includes both sequential and random data testing patterns.
	 *
	 * @param set The hash map implementation to test
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
	
	/** Removes all key-value pairs from this map. */
	public abstract void clear();
	
	/**
	 * Tests whether the specified key is present in this map.
	 * @param key The key to search for
	 * @return true if the key exists, false otherwise
	 */
	public final boolean contains(int key){
		return findCell(key)>=0;
	}
	
	/**
	 * Tests whether the specified key is present in this map.
	 * Functionally identical to contains().
	 * @param key The key to search for
	 * @return true if the key exists, false otherwise
	 */
	public final boolean containsKey(int key){
		return findCell(key)>=0;
	}
	
	/**
	 * Returns the value associated with the specified key.
	 * @param key The key to look up
	 * @return The value associated with the key, or implementation-defined default if not found
	 */
	public abstract int get(int key);
	
	/**
	 * Set key to value.
	 * @param key
	 * @param value
	 * @return Old value.
	 */
	public abstract int put(int key, int value);
	
	/**
	 * Set key to value.
	 * @param key
	 * @param value
	 * @return Old value.
	 */
	public abstract int set(int key, int value);
	
	/**
	 * Increment key's value by 1.
	 * @param key
	 * @return New value.
	 */
	public abstract int increment(int key);
	
	/**
	 * Increment key's value.
	 * @param key
	 * @param incr
	 * @return New value.
	 */
	public abstract int increment(int key, int incr);

	/**
	 * Increments all keys in this map by their corresponding values from another map.
	 * Only processes valid (non-invalid) key-value pairs from the source map.
	 * @param map The source map containing increment values
	 */
	public final void incrementAll(AbstractIntHashMap map) {
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
	 * Remove this key.
	 * @param key
	 * @return true if the key was removed, false if it was not present.
	 */
	public abstract boolean remove(int key);
	
	/**
	 * Returns the number of key-value pairs in this map.
	 * Alias for size().
	 * @return The number of mappings
	 */
	public final int cardinality(){return size();}
	/** Returns the number of key-value pairs in this map.
	 * @return The number of mappings */
	public abstract int size();
	/** Returns true if this map contains no key-value pairs.
	 * @return true if empty, false otherwise */
	public abstract boolean isEmpty();
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns a detailed string representation showing internal structure.
	 * Includes array indices along with key-value pairs as (index, key, value) tuples.
	 * @return Detailed string representation with indices
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
	 * Returns a string representation of this map in list format.
	 * Shows only key-value pairs as (key,value) tuples without internal indices.
	 * @return String representation showing key-value pairs
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
	 * Returns an array containing all keys in this map.
	 * Invalid entries are excluded from the result.
	 * @return Array of all valid keys
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
	 * Creates a histogram of the values in this map.
	 * Each unique value becomes a key in the result, mapped to its frequency count.
	 * @return New IntHashMapBinary containing value frequency counts
	 */
	public final IntHashMapBinary toCountHistogram() {
		IntHashMapBinary counts=new IntHashMapBinary(Tools.mid(1, size(), 64));
		final int[] keys=keys();
		final int[] values=values();
		final int invalid=invalid();
		for(int i=0; i<keys.length; i++) {
			int a=keys[i];
			int b=values[i];
			if(a!=invalid){
				counts.increment(b);
			}
		}
		return counts;
	}
	
	/**
	 * Validates the internal consistency of this hash map.
	 * Checks that key and value arrays have matching lengths, invalid entries
	 * have zero values, and all keys can be found at their expected positions.
	 * @return true if the hash map structure is valid, false if corrupted
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
				final int cell=findCell(i);
				if(i==cell){
					numFound++;
				}else{
					return false;
				}
			}
		}
		return numValues==numFound && numValues==size;
	}
	
	abstract int findCell(final int key);
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the internal array representation of this map.
	 * Implementation-specific method for accessing underlying storage.
	 * @return Array representation of the map data
	 */
	public abstract int[] toArray();
	/** Returns the internal keys array.
	 * @return Array containing all keys (including invalid entries) */
	public abstract int[] keys();
	/** Returns the internal values array.
	 * @return Array containing all values */
	public abstract int[] values();
	/** Returns the sentinel value used to mark invalid/empty entries.
	 * @return The invalid key marker value */
	public abstract int invalid();
	
	/** Bit mask for positive integers (Integer.MAX_VALUE). */
	static final int MASK=Integer.MAX_VALUE;
	/** Bit mask for the minimum integer value (Integer.MIN_VALUE). */
	static final int MINMASK=Integer.MIN_VALUE;
	
	/** Additional capacity factor for internal array sizing. */
	static final int extra=10;
	
}
