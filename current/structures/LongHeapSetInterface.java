package structures;

/**
 * Maintains a heap of unique values.
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public interface LongHeapSetInterface {
	
	/**
	 * Adds a value to the heap-set if not already present.
	 * If heap is at capacity, replaces the smallest value only if the new value is larger.
	 * Maintains uniqueness by checking set membership before insertion.
	 *
	 * @param key The value to add
	 * @return true if the value was added, false if it was a duplicate or too small
	 */
	public boolean add(long key);
	
	/**
	 * Increments the count associated with a key.
	 * For simple heap-set implementations, this typically behaves like add().
	 *
	 * @param key The key to increment
	 * @param incr The increment amount
	 * @return The new count value after incrementing
	 */
	public int increment(long key, int incr);
	
	/** Removes all elements from the heap-set.
	 * Resets both heap and set components to empty state. */
	public void clear();
	
	/** Returns the current number of elements in the heap-set.
	 * @return The number of unique values currently stored */
	public int size();
	
	/** Returns the maximum number of elements this heap-set can contain.
	 * @return The fixed capacity limit */
	public int capacity();
	
	/** Checks if there is space for additional elements without eviction.
	 * @return true if size is less than capacity, false otherwise */
	public boolean hasRoom();
	
	/**
	 * Provides access to the underlying heap structure.
	 * Allows direct manipulation of heap ordering and iteration over values.
	 * @return The internal LongHeap maintaining value ordering
	 */
	public LongHeap heap();
	
	/**
	 * Returns the smallest value in the heap-set without removing it.
	 * This is the value that would be evicted next if a larger value is added at capacity.
	 * @return The minimum value currently stored
	 */
	public long peek();

	/**
	 * Tests if a specific value is present in the heap-set.
	 * @param key The value to search for
	 * @return true if the value is found, false otherwise
	 */
	public boolean contains(long key);
	
}

