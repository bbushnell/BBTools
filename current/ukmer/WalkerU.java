package ukmer;

public abstract class WalkerU {
	
	public abstract boolean next();
	
	/**
	 * Returns the current k-mer object (key) at the iterator position.
	 * Must be called after a successful next() call.
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
