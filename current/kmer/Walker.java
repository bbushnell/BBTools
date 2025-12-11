package kmer;

public abstract class Walker {

	
	public abstract boolean next();
	
	/**
	 * Returns the current k-mer key as a long-encoded value.
	 * Should be called after a successful next() operation.
	 * @return The current k-mer encoded as a long
	 */
	public abstract long kmer();
	
	/**
	 * Returns the value associated with the current k-mer.
	 * Should be called after a successful next() operation.
	 * @return The value associated with the current k-mer
	 */
	public abstract int value();
	
}
