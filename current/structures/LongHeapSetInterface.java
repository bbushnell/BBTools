package structures;

/**
 * Interface for fixed-capacity heap-backed collections of unique long keys.
 * Implementations keep a min-heap as an eviction floor. When they are full, an
 * absent key must exceed the current minimum to be retained.
 *
 * @author Brian Bushnell
 * @date July 6, 2016
 */
public interface LongHeapSetInterface{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds or updates a key according to the implementation's semantics.
	 * Set implementations insert absent keys and reject retained duplicates.
	 * Count-tracking implementations may increment existing keys and report
	 * an implementation-specific count-update result. At capacity, new keys are
	 * retained only if they exceed the current eviction floor.
	 * @param key Key to add or update.
	 * @return Implementation-specific insertion or count-update result.
	 */
	public boolean add(long key);

	/**
	 * Increments the count for a key, or adds it with set semantics.
	 * Count-tracking implementations return the updated count. Set-only
	 * implementations may ignore incr and return 1 only for a successful add.
	 * @param key Key to increment.
	 * @param incr Amount to add.
	 * @return Updated count or set-style add result, depending on implementation.
	 */
	public int increment(long key, int incr);

	/** Removes all retained keys and backing heap state. */
	public void clear();

	/**
	 * Returns the current number of retained unique keys.
	 * @return Unique-key count.
	 */
	public int size();

	/**
	 * Returns the fixed capacity of the heap-backed collection.
	 * @return Maximum unique-key count.
	 */
	public int capacity();

	/** Indicates whether there is space without eviction.
	 * @return true if {@code size < capacity} */
	public boolean hasRoom();

	/**
	 * Returns access to the underlying heap structure.
	 * @return Internal {@link LongHeap} maintaining the eviction floor.
	 */
	public LongHeap heap();

	/** Returns the smallest value without removal.
	 * @return Minimum stored value, or {@link LongHeap#EMPTY} if empty */
	public long peek();

	/**
	 * Tests whether a value is present in the heap-set.
	 * @param key Value to search for
	 * @return true if contained
	 */
	public boolean contains(long key);

}

