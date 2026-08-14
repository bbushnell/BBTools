package tax;

/**
 * Abstract batch reader for disk-backed hash tables.
 * Implementations vary the I/O strategy; the hash table
 * format and probe logic are shared.
 *
 * @author Brian Bushnell, Chloe
 */
public abstract class DiskAccessor {

	/**
	 * Look up multiple digitized keys and return their taxIDs.
	 * @param keys Array of digitized accession keys
	 * @return Array of taxIDs (same length as keys; -1 for not found)
	 */
	public abstract int[] getBatch(long[] keys);

	/** Name for benchmark display */
	public abstract String name();

	/** Close underlying file handles */
	public abstract void close() throws Exception;

	/*--------------------------------------------------------------*/
	/*----------------      Shared hash logic       ----------------*/
	/*--------------------------------------------------------------*/

	/** Compute initial partition and slot for a key */
	static long hash(long key){
		key^=(key>>>33);
		key*=0xff51afd7ed558ccdL;
		key^=(key>>>33);
		key*=0xc4ceb9fe1a85ec53L;
		key^=(key>>>33);
		return key&Long.MAX_VALUE;
	}

	static final long SLOT_SIZE=12;
	static final long EMPTY=0;
}
