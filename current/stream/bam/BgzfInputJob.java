package stream.bam;

import stream.HasID;

/**
 * Carries BGZF input metadata and mutable decompression results.
 * Input metadata fields are final, but compressed-array contents and result
 * fields remain mutable. Arrays are retained without copying. Callers coordinate
 * ownership and publication; this class supplies no synchronization for result writes.
 * Poison identity, the last-job flag, ID ordering, and payload checks are independent.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 14, 2025
 */
public class BgzfInputJob implements HasID, Comparable<BgzfInputJob>{

	/** Sequential job ID for maintaining output order */
	public final long id;

	/** Borrowed compressed bytes; may be null for control or predecoded jobs. */
	public final byte[] compressed;

	/** Expected CRC32 metadata, not validated by this container. */
	public final long expectedCrc;

	/** Expected output-size metadata; repOK does not compare it with the result size. */
	public final int expectedSize;

	/** Last-job flag, independent of poison identity. */
	public final boolean lastJob;

	/** Mutable decoded buffer supplied by processing code; only its used prefix is data. */
	public byte[] decompressed;

	/** Used byte count in decompressed, assigned by processing code. */
	public int decompressedSize;

	/** Mutable processing exception, not interpreted by repOK or the marker methods. */
	public Exception error;

	/** Shared identity marker with maximum ID and lastJob=false; do not use as a work buffer. */
	public static final BgzfInputJob POISON_PILL=new BgzfInputJob(Long.MAX_VALUE, null, 0, 0, false);

	/**
	 * Stores metadata without copying the input or invoking repOK.
	 * Result fields initially have their Java defaults: null buffer/error and zero size.
	 * @param id_ Ordering ID, retained as supplied
	 * @param compressed_ Borrowed compressed bytes, or null for a control/predecoded job
	 * @param expectedCrc_ Expected checksum metadata
	 * @param expectedSize_ Expected output-size metadata
	 * @param last_ Whether last() reports end of the sequence
	 */
	public BgzfInputJob(long id_, byte[] compressed_,
			long expectedCrc_, int expectedSize_, boolean last_){
		id=id_;
		compressed=compressed_;
		expectedCrc=expectedCrc_;
		expectedSize=expectedSize_;
		lastJob=last_;
	}

	/**
	 * Compares IDs without separate handling for poison or last flags.
	 * @param other Non-null job to compare
	 * @return Negative, zero, or positive according to the two IDs
	 */
	@Override
	public int compareTo(BgzfInputJob other){return Long.compare(this.id(), other.id());}

	/** Returns the stored ordering ID. */
	@Override
	public long id(){return id;}

	/** Returns true only for the shared POISON_PILL instance. */
	@Override
	public boolean poison(){return this==POISON_PILL;}

	/** Returns the stored last-job flag, independently of poison identity. */
	@Override
	public boolean last(){return lastJob;}

	/**
	 * Returns the shared identity marker rather than allocating a requested-ID marker.
	 * @param id_ Ignored; the singleton retains Long.MAX_VALUE
	 * @return The shared POISON_PILL, whose last-job flag is false
	 */
	//Historical queue-order rationale follows; this container review does not establish
	//its applicability to every caller or certify shutdown behavior. HasID documents this
	//fixed-ID factory variant; its consumer must support that marker policy.
	//Poison-by-IDENTITY (poison()/isPoisonPill() test `this==POISON_PILL`), unlike BgzfJob's poison-by-flag.
	//makePoison ignores id_ and returns the singleton (id=Long.MAX_VALUE), so a JobQueue post-last poison
	//(added as makePoison(job.id()+1)) sorts at MAX_VALUE instead of job.id()+1 - HARMLESS, because once a
	//`last` job is dequeued JobQueue sets lastSeen and the next take() returns null BEFORE that poison is ever
	//released (JobQueue.take L145). So the poison's id never matters on the consumer side; the singleton is
	//only meaningfully consumed on the WORKER inputQueue, where isPoisonPill() identity is exact.
	@Override
	public BgzfInputJob makePoison(long id_){return POISON_PILL;}

	/**
	 * Allocates an empty last marker with the requested ID and no payload or result.
	 * This fresh marker is not poison and fails repOK while both arrays remain null.
	 * @param id_ Ordering ID for the new marker
	 * @return New last-job instance with null arrays and zero checksum/size metadata
	 */
	@Override
	public BgzfInputJob makeLast(long id_){return new BgzfInputJob(id_, null, 0, 0, true);}

	/** Returns true only for the shared POISON_PILL instance, like poison(). */
	public boolean isPoisonPill(){return this==POISON_PILL;}

	/**
	 * Checks the stored payload representation, exempting the poison singleton.
	 * Other jobs need a nonnegative ID and at least one non-null array. A present
	 * decoded buffer must be at most 65536 bytes, with its used size in [0, length].
	 * When that buffer is null, decompressedSize is not checked. This does not
	 * validate compressed bytes, checksum, expectedSize, error, or last-job semantics;
	 * an empty factory marker or prototype can therefore fail this payload check.
	 * @return Whether these representation checks pass, or this is POISON_PILL
	 */
	public boolean repOK(){
		if(this==POISON_PILL){return true;}
		if(compressed==null && decompressed==null){return false;}
		if(decompressed!=null && (decompressedSize<0 || decompressedSize>decompressed.length)){return false;}
		if(decompressed!=null && decompressed.length>65536){return false;}
		if(id<0){return false;}
		return true;
	}
}
