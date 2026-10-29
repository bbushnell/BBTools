package tax;

/**
 * Batch lookup API for disk-backed accession tables.
 * Implementations vary the I/O strategy while preserving the table's digitized-key
 * contract: each queried key maps to one taxID or to -1 when absent.
 *
 * @author Brian Bushnell, Chloe
 */
public abstract class DiskAccessor{
	/** Constructs a batch accessor. */
	public DiskAccessor(){}

	/**
	 * Looks up multiple digitized keys and returns one taxID per input position.
	 * @param keys Array of digitized accession keys
	 * @return Array of taxIDs, with -1 for absent keys
	 */
	public abstract int[] getBatch(long[] keys);

	/**
	 * Returns the display name used by disk-accessor comparisons.
	 * @return Human-readable accessor name
	 */
	public abstract String name();

	/**
	 * Closes any underlying file handles.
	 * @throws Exception If handle closure fails
	 */
	public abstract void close() throws Exception;

	/*--------------------------------------------------------------*/
	/*----------------      Shared hash logic       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Computes the nonnegative MurmurHash3 finalizer value used to select a
	 * partition and initial linear-probe slot.
	 * @param key Digitized accession key
	 * @return Nonnegative 63-bit hash value
	 */
	static long hash(long key){
		key^=(key>>>33);
		key*=0xff51afd7ed558ccdL;
		key^=(key>>>33);
		key*=0xc4ceb9fe1a85ec53L;
		key^=(key>>>33);
		return key&Long.MAX_VALUE;
	}

	/** Bytes per disk slot: eight-byte key followed by four-byte taxID. */
	static final long SLOT_SIZE=12;
	/** On-disk key sentinel marking an empty slot and probe termination. */
	static final long EMPTY=0;
}
