package stream.bam;

import stream.HasID;

/**
 * Data shuttle for multithreaded BGZF compression/decompression.
 *
 * The same job object can be used for both reading and writing:
 * - For decompression: Producer fills compressed[], worker fills decompressed[]
 * - For compression: Producer fills decompressed[], worker fills compressed[]
 *
 * Jobs are identified by sequential IDs to maintain output order even when
 * workers complete out of order.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BgzfJob implements HasID, Comparable<BgzfJob> {

	/** Compressed BGZF block data (with header/footer for reading, without for writing) */
	public byte[] compressed;

	/** Decompressed block data (always 64KB or less) */
	public byte[] decompressed;

	/** Sequential job ID for maintaining output order */
	public final long id;

	/** Actual bytes used in compressed array */
	public int compressedSize;

	/** Actual bytes used in decompressed array */
	public int decompressedSize;

	/** Exception caught by worker thread during processing */
	public Exception error;

	/** Flag indicating this is the last job (no more jobs will be produced) */
	public final boolean lastJob;

	/** Poison pill marker for worker shutdown */
	public static final BgzfJob POISON_PILL = new BgzfJob(Long.MAX_VALUE, null, null, false);

	public BgzfJob(long id, byte[] raw, byte[] comp, boolean last) {
		this.id = id;
		decompressed=raw;
		compressed=comp;
		lastJob=last;
	}

	@Override
	public int compareTo(BgzfJob other) {
		return Long.compare(this.id(), other.id());
	}
	
	@Override
	public long id() {return id;}
	
	@Override
	public boolean poison() {return this==POISON_PILL;}
	
	@Override
	public boolean last() {return lastJob;}
	
	@Override
	public BgzfJob makePoison(long id_) {
		return new BgzfJob(id_, null, null, false);
	}
	
	@Override
	public BgzfJob makeLast(long id_){
		return new BgzfJob(id_, null, null, true);
	}

	/**
	 * Check if this job is the poison pill marker.
	 * Workers should exit their loop when receiving this.
	 */
	public boolean isPoisonPill() {
		return this == POISON_PILL;
	}

	/**
	 * Validate job state for debugging.
	 * Called with assert(!debugging || repOK()) pattern.
	 */
	public boolean repOK() {
		// At least one array should have data
		if (compressed == null && decompressed == null) return false;

		// Sizes should be non-negative and within array bounds
		if (compressed != null && (compressedSize < 0 || compressedSize > compressed.length)) {
			return false;
		}
		if (decompressed != null && (decompressedSize < 0 || decompressedSize > decompressed.length)) {
			return false;
		}

		// Decompressed data should never exceed BGZF max block size (64KB)
		if (decompressed != null && decompressed.length > 65536) {
			return false;
		}

		// ID should be non-negative (except for poison pill)
		if (this != POISON_PILL && id < 0) {
			return false;
		}

		return true;
	}
}
