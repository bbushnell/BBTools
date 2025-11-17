package ukmer;

/**
 * Abstract iterator for traversing k-mer hash tables and data structures.
 * Provides unified iteration interface for k-mer table traversal and data extraction.
 * Supports sequential access to k-mer entries with key-value pair retrieval.
 * @author Brian Bushnell
 */
public abstract class WalkerU {
	
	/** 
	 * Allows iteration through a hash map.
	 * Concurrent modification is not recommended.
	 */
	public abstract boolean next();
	
	/** Current object kmer (key) for ukmer package */
	public abstract Kmer kmer();
	
	/** Current value */
	public abstract int value();
	
}
