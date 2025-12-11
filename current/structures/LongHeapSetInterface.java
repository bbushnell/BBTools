package structures;

/**
 * Interface for heap-based collections that maintain unique values.
 * Combines heap and set functionality for efficient top-N selection with uniqueness constraints.
 * Implementations maintain a fixed-capacity heap where new values replace smaller values
 * when capacity is exceeded, while ensuring no duplicate values are stored.
 *
 * @author Brian Bushnell
 * @date July 6, 2016
 */
public interface LongHeapSetInterface {
	
	/**
	 * Adds a value if absent; at capacity, replaces the smallest only if the new value is larger.
	 * Maintains uniqueness by checking membership before insertion.
	 * @param key Value to add
	 * @return true if added; false if duplicate or too small
	 */
	public boolean add(long key);
	
	/**
	 * Increments the count associated with a key (or behaves like add in simple implementations).
	 * @param key Key to increment
	 * @param incr Amount to add
	 * @return Updated count
	 */
	public int increment(long key, int incr);
	
	/** Removes all elements, resetting heap and set state. */
	public void clear();
	
	/** Returns the current number of stored elements.
	 * @return Unique element count */
	public int size();
	
	/** Returns the fixed capacity of the heap-set.
	 * @return Maximum element count */
	public int capacity();
	
	/** Indicates whether there is space without eviction.
	 * @return true if size < capacity */
	public boolean hasRoom();
	
	/** Returns access to the underlying heap structure.
	 * @return Internal LongHeap maintaining ordering */
	public LongHeap heap();
	
	/** Returns the smallest value without removal (would evict next at capacity).
	 * @return Minimum stored value */
	public long peek();

	/**
	 * Tests whether a value is present in the heap-set.
	 * @param key Value to search for
	 * @return true if contained
	 */
	public boolean contains(long key);
	
}

