package map;

import java.io.Serializable;

/**
 * Contract for mutable int-to-int maps that expose raw storage to package utilities.
 * Entries equal to {@link #invalid()} in {@link #keys()} are empty cells. Instances are
 * unsynchronized. The borrowed key array encodes probe structure; external key writes can
 * break lookup and rehash invariants. The borrowed value array is direct int storage
 * index-aligned with {@link #keys()} and should only be edited with occupancy from keys.
 */
public interface IntHashMapInterface extends Serializable{

	/** Removes every live mapping while preserving reusable storage. */
	public void clear();

	/**
	 * Tests whether a key is currently mapped.
	 * @param key Key to look up
	 * @return true if {@code key} is live in this map
	 */
	public boolean contains(int key);

	/**
	 * Tests whether a key is currently mapped.
	 * @param key Key to look up
	 * @return true if {@code key} is live in this map
	 */
	public boolean containsKey(int key);

	/**
	 * Returns the int value mapped from a key.
	 * @param key Key to look up
	 * @return Mapped value, or the implementation-specific absent default
	 */
	public int get(int key);

	/**
	 * Sets a key to a value, inserting the key if absent and replacing its old value if present.
	 *
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous mapped value, or the implementation-specific absent default
	 */
	public int put(int key, int value);

	/**
	 * Sets a key to a value, inserting the key if absent and replacing its old value if present.
	 *
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous mapped value, or the implementation-specific absent default
	 */
	public int set(int key, int value);

	/**
	 * Adds one to a key's value, inserting the key if absent.
	 * @param key Key to insert or update
	 * @return New value after the implementation-specific overflow policy
	 */
	public int increment(int key);

	/**
	 * Adds a signed delta to a key's value, inserting the key if absent.
	 *
	 * @param key Key to insert or update
	 * @param incr Signed increment to apply
	 * @return New value after the implementation-specific overflow and underflow policy
	 */
	public int increment(int key, int incr);

	/**
	 * Removes a key if present.
	 * @param key Key to remove
	 * @return true if a live mapping was removed
	 */
	public boolean remove(int key);

	/**
	 * Returns the number of live mappings.
	 * @return Number of live mappings
	 */
	public int size();
	/**
	 * Tests whether the map has no live mappings.
	 * @return true when the map has no live mappings
	 */
	public boolean isEmpty();

	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Formats the live mappings.
	 * @return Implementation-defined string representation
	 */
	public String toString();

//		int findCell(final int key);

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies live keys into a dense array.
	 * @return Owned key snapshot
	 */
	public abstract int[] toArray();

	/**
	 * Borrows the raw key array.
	 * @return Borrowed raw key array; entries equal to {@link #invalid()} are empty cells
	 */
	public abstract int[] keys();

	/**
	 * Borrows the raw int value array.
	 * @return Borrowed raw int value array, index-aligned with {@link #keys()}
	 */
	public abstract int[] values();

	/**
	 * Returns the current empty-cell sentinel.
	 * @return Current empty-cell sentinel for {@link #keys()}
	 */
	public abstract int invalid();

}
