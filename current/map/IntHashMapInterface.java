package map;

import java.io.Serializable;

public interface IntHashMapInterface extends Serializable{

		public void clear();
		
		public boolean contains(int key);
		
		public boolean containsKey(int key);
		
		public int get(int key);
		
		/**
		 * Associates the specified value with the specified key in this map.
		 * If the key was previously mapped, the old value is replaced.
		 *
		 * @param key The key to associate with the value
		 * @param value The value to associate with the key
		 * @return The previous value associated with the key, or implementation-defined default if none
		 */
		public int put(int key, int value);
		
		/**
		 * Associates the specified value with the specified key in this map.
		 * Functionally identical to put().
		 *
		 * @param key The key to associate with the value
		 * @param value The value to associate with the key
		 * @return The previous value associated with the key, or implementation-defined default if none
		 */
		public int set(int key, int value);
		
		/**
		 * Increments the value associated with the specified key by 1.
		 * If the key does not exist, it is created with value 1.
		 * @param key The key whose value should be incremented
		 * @return The new value after incrementing
		 */
		public int increment(int key);
		
		/**
		 * Increments the value associated with the specified key by the specified amount.
		 * If the key does not exist, it is created with the increment value.
		 *
		 * @param key The key whose value should be incremented
		 * @param incr The amount to increment by (may be negative)
		 * @return The new value after incrementing
		 */
		public int increment(int key, int incr);
		
		/**
		 * Removes the mapping for the specified key from this map if present.
		 * @param key The key whose mapping is to be removed
		 * @return true if the key was present and removed, false otherwise
		 */
		public boolean remove(int key);
		
		public int size();
		public boolean isEmpty();
		
		/*--------------------------------------------------------------*/
		/*----------------        String Methods        ----------------*/
		/*--------------------------------------------------------------*/
		
		/**
		 * Returns a string representation of this map in list view format.
		 * Shows key-value pairs as (key,value) tuples.
		 * @return String representation of the map
		 */
		public String toString();
		
//		int findCell(final int key);
		
		/*--------------------------------------------------------------*/
		/*----------------            Fields            ----------------*/
		/*--------------------------------------------------------------*/

		public abstract int[] toArray();
		public abstract int[] keys();
		public abstract int[] values();
		public abstract int invalid();
		
	}
