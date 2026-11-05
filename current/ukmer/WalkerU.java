package ukmer;

/**
 * Stateful cursor over an occupied long-key table. Call next before accessing
 * a record. Traversal order and filtering depend on the implementation; callers
 * must keep the underlying table stable while walking it. No removal operation
 * is provided. Accessors are valid only after a successful next(); no common
 * exhausted value or key is defined. The cursor provides no synchronization for
 * shared access.
 */
public abstract class WalkerU{

	/**
	 * Advances to the next occupied table entry.
	 *
	 * @return true when {@link #kmer()} and {@link #value()} expose the next entry
	 */
	public abstract boolean next();

	/**
	 * Returns the current k-mer object (key) at the iterator position.
	 * Must be called after a successful next() call. The returned object may be
	 * reused by a later {@link #next()}, so callers retaining a key must copy it.
	 * @return Current k-mer object serving as the key
	 */
	public abstract Kmer kmer();

	/**
	 * Returns the current value associated with the k-mer at the iterator position.
	 * Must be called after a successful next() call.
	 * @return Current integer value associated with the k-mer
	 */
	public abstract int value();

}
