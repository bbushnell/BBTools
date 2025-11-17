package structures;
import java.util.Arrays;

import shared.Tools;


/**
 * Array-backed integer-to-integer map for contiguous key ranges.
 * Optimized for dense mappings where keys fall within a known range.
 * Uses Integer.MIN_VALUE as sentinel value for unset entries.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class IntMap {
	
	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		
	}
	
	
	/**
	 * Constructs an IntMap for the specified key range.
	 * @param from Minimum key value (inclusive)
	 * @param to Maximum key value (inclusive)
	 */
	public IntMap(int from, int to){
		reset(from, to);
	}
	
	
	/**
	 * Retrieves the value associated with the specified key.
	 * @param key The key to look up
	 * @return The associated value, or INVALID if key not present
	 */
	public int get(int key){
		assert(key>=min && key<=max);
		return array[key-min];
	}
	
	
	/**
	 * Checks if the map contains a mapping for the specified key.
	 * @param key The key to check
	 * @return true if key is present, false otherwise
	 */
	public boolean containsKey(int key){
		assert(key>=min && key<=max);
		return array[key-min]!=INVALID;
	}
	
	
	/**
	 * Increments the value at the specified key by 1.
	 * @param key The key to increment
	 * @return The new value after incrementing
	 */
	public int increment(int key){
		return increment(key, 1);
	}
	
	
	/**
	 * Increments the value at the specified key by the given amount.
	 * Prevents integer overflow by capping at Integer.MAX_VALUE.
	 *
	 * @param key The key to increment
	 * @param value The amount to add
	 * @return The new value after incrementing
	 */
	public int increment(int key, int value){
		assert(key>=min && key<=max);
		int index=key-min;
		int old=array[index];
		int v2=(int)Tools.min(Integer.MAX_VALUE, old+(long)value);
		assert(array[index]!=INVALID);
		array[index]=v2;
		return v2;
	}
	
	
	/**
	 * Associates the specified value with the specified key.
	 * @param key The key to map
	 * @param value The value to associate (must not be INVALID)
	 * @return The previous value, or INVALID if key was not present
	 */
	public int put(int key, int value){
		assert(key>=min && key<=max);
		assert(value!=INVALID);
		int index=key-min;
		int old=array[index];
		array[index]=value;
		return old;
	}
	
	
	/**
	 * Removes the mapping for the specified key.
	 * @param key The key to remove
	 * @return The previous value, or INVALID if key was not present
	 */
	public int remove(int key){
		assert(key>=min && key<=max);
		int index=key-min;
		int old=array[index];
		array[index]=INVALID;
		return old;
	}
	
	
	/** Returns the number of key-value mappings in this map.
	 * @return The count of non-INVALID entries */
	public int size(){
		int sum=0;
		for(int i=0; i<array.length; i++){
			if(array[i]!=INVALID){sum++;}
		}
		return sum;
	}
	
	
	/** Returns an array containing all keys that have mappings.
	 * @return Array of keys in ascending order */
	public int[] keys(){
		int[] r=new int[size()];
		for(int i=0, j=0; j<r.length; i++){
			if(array[i]!=INVALID){
				r[j]=(min+i);
				j++;
			}
		}
		return r;
	}
	
	
	/** Returns an array containing all mapped values.
	 * @return Array of values corresponding to mapped keys */
	public int[] values(){
		int[] r=new int[size()];
		for(int i=0, j=0; j<r.length; i++){
			if(array[i]!=INVALID){
				r[j]=array[i];
				j++;
			}
		}
		return r;
	}
	
	
	/** Removes all mappings by setting all entries to INVALID */
	public void clear(){
		Arrays.fill(array, INVALID);
	}
	
	
	/**
	 * Resets the map to support a new key range.
	 * Reallocates the backing array if necessary and clears all entries.
	 * @param from New minimum key value (inclusive)
	 * @param to New maximum key value (inclusive)
	 */
	public void reset(int from, int to){
		min=from;
		max=to;
		assert(max>=min);
		assert(((long)max)-((long)min)<Integer.MAX_VALUE);
		
		int size=max-min+1;
		if(array==null || array.length<size){
			array=new int[size];
		}
		clear();
	}
	
	
	/** Minimum key value supported by this map */
	public int min;
	/** Maximum key value supported by this map */
	public int max;
	/** Backing array storing values indexed by (key - min) */
	public int[] array;
	
	/** Sentinel value indicating an unset or removed entry */
	private static final int INVALID=Integer.MIN_VALUE;
	
}
