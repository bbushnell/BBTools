package kmer;

/**
 * Stateful cursor over a short-key table. Call next before accessing a record.
 * Traversal order and filtering depend on the implementation; callers must keep
 * the underlying table stable while walking it. No removal operation is provided.
 * Accessors are valid only after a successful next(); no common exhausted value
 * or key is defined. The cursor provides no synchronization for shared access.
 */
public abstract class Walker{

	/** Advances the cursor; returns false when no further record is available. */
	public abstract boolean next();

	/**
	 * Returns the current k-mer key as a long-encoded value.
	 * Should be called after a successful next() operation.
	 * @return The current k-mer encoded as a long
	 */
	public abstract long kmer();

	/**
	 * Returns the stored integer value associated with the current k-mer.
	 * Counting-table walkers expose the count through this accessor.
	 * Should be called after a successful next() operation.
	 * @return The value associated with the current k-mer
	 */
	public abstract int value();

}
