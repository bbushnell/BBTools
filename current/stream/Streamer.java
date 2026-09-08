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

	/**
	 * Deterministic positional subsampler: the keep/drop decision for the record at position
	 * recordNum (0-based in the file) under the given seed and rate. A pure function of its
	 * arguments (SplitMix64-style mix), so it is thread-safe with no shared state, reproducible
	 * across runs and thread counts, and mate-safe: R1 and R2 streamers sampling the same
	 * positions with the same seed keep exactly the same subset. Replaces the shared-PRNG
	 * scheme, whose keep-decisions occurred in worker-scheduling order: under MT that desynced
	 * twin files, broke mate pairing, and killed the consumer on PairStreamer's numericID
	 * assert, hanging the JVM (replicated via stream.sh samplerate=0.5 threadsin=4, 2026-09-05).
	 */
	public static boolean sampleKeep(long recordNum, long seed, float rate){
		long x=recordNum*0x9E3779B97F4A7C15L+seed;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		x^=(x>>>31);
		return (x>>>40)<(long)(rate*0x1p24f);//top 24 bits vs rate scaled to 2^24
	}

	/** Resolve a user-supplied sampling seed: negative means "pick a random seed". */
	public static long resolveSampleSeed(long seed){
		return seed>=0 ? seed : new java.util.Random().nextLong();
	}

}