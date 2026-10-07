package stream.bam;

import stream.HasID;

/**
 * Carries buffers, lengths and status for BGZF compression or decompression.
 * Processing code supplies compressed input and decoded output when reading,
 * or decoded input and compressed output when writing. Buffer layout depends on
 * the processing stage and can include metadata as well as payload.
 * Arrays are borrowed without copying, and buffers, lengths and error are mutable.
 * Callers coordinate ownership and publication; this class supplies no synchronization.
 * Callers also assign ordering IDs; comparison uses only those IDs.
 *
 * @author Chloe
 * @contributor Collei
 * @date October 18, 2025
 */
public class BgzfJob implements HasID, Comparable<BgzfJob>{

	/** Mutable compressed bytes; framing and used length are defined by the processing stage. */
	public byte[] compressed;

	/** Mutable decoded buffer or workspace; the input path can prepend eight metadata bytes. */
	public byte[] decompressed;

	/** Sequential job ID for maintaining output order */
	public final long id;

	/** Used length associated with compressed, assigned by processing code. */
	public int compressedSize;

	/** Used length associated with decompressed, interpreted by the processing stage. */
	public int decompressedSize;

	/** Mutable processing exception; marker methods and repOK do not interpret it. */
	public Exception error;

	/** Last-job flag, retained independently of the poison flag. */
	public final boolean lastJob;

	/** Poison status used by poison() and isPoisonPill(), independent of object identity. */
	public final boolean isPoison;

	/** Shared poison marker with maximum ID and lastJob=false; do not use as a work buffer. */
	public static final BgzfJob POISON_PILL=new BgzfJob(Long.MAX_VALUE, null, null, false, true);

	/**
	 * Retains buffers and metadata with the poison flag disabled.
	 * @param id Ordering ID
	 * @param raw Borrowed decoded bytes or workspace, possibly null
	 * @param comp Borrowed compressed bytes, possibly null
	 * @param last Whether this marks the last job
	 */
	public BgzfJob(long id, byte[] raw, byte[] comp, boolean last){this(id, raw, comp, last, false);}

	/**
	 * Retains arguments without copying buffers or invoking repOK.
	 * Used lengths start at zero and error starts null, regardless of buffer capacity.
	 * @param id Ordering ID, retained as supplied
	 * @param raw Borrowed decoded bytes or workspace, possibly null
	 * @param comp Borrowed compressed bytes, possibly null
	 * @param last Last-job flag
	 * @param poison Poison flag, stored independently of last
	 */
	public BgzfJob(long id, byte[] raw, byte[] comp, boolean last, boolean poison){
		this.id=id;
		decompressed=raw;
		compressed=comp;
		lastJob=last;
		isPoison=poison;
	}

	/**
	 * Compares IDs without separately ordering flags or payloads.
	 * @param other Non-null job to compare
	 * @return Negative, zero or positive according to the two ordering IDs
	 */
	@Override
	public int compareTo(BgzfJob other){return Long.compare(this.id(), other.id());}

	/** Returns the stored ordering ID. */
	@Override
	public long id(){return id;}

	/** Returns the stored poison flag; singleton identity is not required. */
	@Override
	public boolean poison(){return isPoison;}

	/** Returns the stored last-job flag, independently of poison status. */
	@Override
	public boolean last(){return lastJob;}

	/**
	 * Allocates an empty poison marker carrying the requested ID.
	 * @param id_ Ordering ID for the new marker
	 * @return New poison job with lastJob=false, null buffers and zero used lengths
	 */
	//Historical caller/order rationale follows; this container review does not certify
	//queue or shutdown behavior for every consumer. HasID now documents both factory forms.
	//Poison-by-FLAG (poison()/isPoisonPill() test the `isPoison` field), unlike BgzfInputJob's poison-by-
	//identity. makePoison HONORS id_ (a fresh poison stamped with the requested id), so a JobQueue post-last
	//poison sorts exactly at job.id()+1 as JobQueue expects - the write-side twin of BgzfInputJob's harmless
	//id-ignoring poison. The shared POISON_PILL singleton (L46) is the worker-shutdown marker on inputQueue.
	@Override
	public BgzfJob makePoison(long id_){return new BgzfJob(id_, null, null, false, true);}

	/**
	 * Allocates an empty last marker carrying the requested ID.
	 * The marker is not poison and fails repOK while both buffers remain null.
	 * @param id_ Ordering ID for the new marker
	 * @return New last job with poison=false, null buffers and zero used lengths
	 */
	@Override
	public BgzfJob makeLast(long id_){return new BgzfJob(id_, null, null, true, false);}

	/**
	 * Reports poison status by flag, like poison(), rather than singleton identity.
	 * @return The stored poison flag
	 */
	public boolean isPoisonPill(){return isPoison;}

	/**
	 * Checks stored sizes and IDs, accepting all poison jobs without further checks.
	 * Other jobs require at least one buffer and a nonnegative ID. For each non-null
	 * buffer, the associated used length must be in [0, capacity]. A decoded buffer
	 * must additionally have capacity at most 65536, including any metadata prefix.
	 * This cap conflicts with one input-side workspace layout, as noted below.
	 * No payload, checksum, error or last-job validation is performed. Null-buffer
	 * lengths are not checked, and an empty last marker fails the buffer-presence check.
	 * @return Whether these representation checks pass, or this is a poison job
	 */
	public boolean repOK(){
		// Poison pills are always OK
		if(isPoison){return true;}

		// At least one array should have data
		if(compressed==null && decompressed==null){return false;}

		// Sizes should be non-negative and within array bounds
		if(compressed!=null && (compressedSize<0 || compressedSize>compressed.length)){return false;}
		if(decompressed!=null && (decompressedSize<0 || decompressedSize>decompressed.length)){return false;}

		// Decompressed data should never exceed BGZF max block size (64KB)
		//TODO: Probable bug [stream/bam/BgzfJob#001] - BgzfInputStreamMT creates an
		//8+65536-byte decoded workspace with a metadata prefix, which this capacity check
		//rejects. Its DEBUG-gated check is currently disabled (STR-197). Reconcile the
		//representation separately; this is a source finding, not an observed runtime failure.
		if(decompressed!=null && decompressed.length>65536){return false;}

		// ID should be non-negative (except for poison pill)
		if(id<0){return false;}

		return true;
	}
}
