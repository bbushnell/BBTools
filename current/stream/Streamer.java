package stream;

import structures.ListNum;

/**
 * Unified interface for multithreaded sequence file readers.
 * Implementations use ordered job queues to parallelize decompression and parsing.
 * 
 * @author Brian Bushnell
 * @date October 30, 2025
 */
public interface Streamer {
	
	/** Source file */
	public String fname();
	
	/** Initialize and start background reading/parsing threads */
	public void start();
	
	/** Emergency shutdown - prefer poisoning via exhausting stream */
	public void close();
	
	/** True if the reads from this stream have their mate set */
	public boolean paired();
	
	/** 0 for R1 (or paired), 1 for R2 */
	public int pairnum();
	
	/** Number of reads processed */
	public long readsProcessed();
	
	/** Number of bases processed */
	public long basesProcessed();
	
	public void setSampleRate(float rate, long seed);

	/** 
	 * Returns next ordered batch of reads, or null when exhausted.
	 * Blocks if data not yet ready. Thread-safe for single consumer.
	 */
	public ListNum<Read> nextList();
	
	/** 
	 * Returns next ordered batch of SamLines (SAM/BAM only).
	 * May return null or throw UnsupportedOperationException for FASTA/FASTQ.
	 */
	public ListNum<SamLine> nextLines();

	/** 
	 * Returns true if more data may be available.
	 * May return false positives but must eventually return false.
	 * Used for pre-allocation optimizations, not correctness.
	 */
	public boolean hasMore();

	/** True if there was an error */
	public boolean errorState();

}