package stream.bam;

import stream.HasID;

/**
 * Data shuttle for BGZF decompression jobs.
 * Immutable metadata fields for efficient thread-safe access.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 14, 2025
 */
public class BgzfInputJob implements HasID, Comparable<BgzfInputJob> {

	/** Sequential job ID for maintaining output order */
	public final long id;

	/** Compressed BGZF block data */
	public final byte[] compressed;

	/** Expected CRC32 checksum of decompressed data */
	public final long expectedCrc;

	/** Expected size of decompressed data */
	public final int expectedSize;

	/** Flag indicating this is the last job */
	public final boolean lastJob;

	/** Decompressed block data (filled by worker) */
	public byte[] decompressed;

	/** Actual bytes decompressed (filled by worker) */
	public int decompressedSize;

	/** Exception caught by worker thread during processing */
	public Exception error;

	/** Poison pill marker for worker shutdown */
	public static final BgzfInputJob POISON_PILL=new BgzfInputJob(Long.MAX_VALUE, null, 0, 0, false);

	public BgzfInputJob(long id_, byte[] compressed_, 
			long expectedCrc_, int expectedSize_, boolean last_){
		id=id_;
		compressed=compressed_;
		expectedCrc=expectedCrc_;
		expectedSize=expectedSize_;
		lastJob=last_;
	}

	@Override
	public int compareTo(BgzfInputJob other){
		return Long.compare(this.id(), other.id());
	}

	@Override
	public long id(){return id;}

	@Override
	public boolean poison(){return this==POISON_PILL;}

	@Override
	public boolean last(){return lastJob;}

	@Override
	public BgzfInputJob makePoison(long id_){
		return POISON_PILL;
	}

	@Override
	public BgzfInputJob makeLast(long id_){
		return new BgzfInputJob(id_, null, 0, 0, true);
	}

	public boolean isPoisonPill(){
		return this==POISON_PILL;
	}

	public boolean repOK(){
		if(this==POISON_PILL){return true;}
		if(compressed==null && decompressed==null){return false;}
		if(decompressed!=null && (decompressedSize<0 || decompressedSize>decompressed.length)){
			return false;
		}
		if(decompressed!=null && decompressed.length>65536){
			return false;
		}
		if(id<0){return false;}
		return true;
	}
}