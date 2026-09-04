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
	 * Blocks if data not yet ready. Thread-safe for a single consumer, and
	 * also safe for multiple concurrent consumer threads calling this
	 * directly (each call atomically claims the next available batch) --
	 * verified 2026-09-04 across every implementation reachable via
	 * StreamerFactory: host-driven classes synchronize the whole read+advance
	 * step, and worker-thread classes hand off through a thread-safe queue
	 * whose terminal/poison marker is re-injected so every concurrent caller
	 * sees end-of-stream, not just the first. See template/A_SampleStreamerMT.java
	 * for the intended multi-consumer usage pattern.
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

	//TODO:  Remove these eventually
	public default void returnList(ListNum<Read> ln) {}
	public default void returnList(long id, boolean b) {}

}