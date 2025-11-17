package kmer;

/**
 * Abstract base class for iterating through k-mer hash map data structures.
 * Provides iterator-like functionality for traversing k-mer collections
 * with access to both keys (k-mers) and their associated values.
 * @author Brian Bushnell
 */
public abstract class Walker {

	
	/** 
	 * Allows iteration through a hash map.
	 * Concurrent modification is not recommended.
	 */
	public abstract boolean next();
	
	/** Current object kmer (key) for kmer package */
	public abstract long kmer();
	
	/** Current value */
	public abstract int value();
	
}
