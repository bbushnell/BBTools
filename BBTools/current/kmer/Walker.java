package kmer;

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
