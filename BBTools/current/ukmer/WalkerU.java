package ukmer;

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
