package stream;

import java.util.ArrayList;

import structures.ListNum;

/**
 * Unified interface for multithreaded sequence file writers.
 * Implementations use ordered job queues to parallelize formatting and compression.
 * 
 * @author Brian Bushnell
 * @date October 30, 2025
 */
public interface Writer {
	
	/** Initialize and start background formatting/writing threads */
	public void start();
	
//	/** Emergency shutdown - prefer poisonAndWait() for clean termination */
//	@Deprecated
//	public void close();
	
	/** Number of reads written */
	public long readsWritten();
	
	/** Number of bases written */
	public long basesWritten();
	
	/**
	 * Submit ordered batch of reads for writing.
	 * Blocks if queue is full. Thread-safe for a single producer, and also
	 * safe for multiple concurrent producer threads submitting in ANY
	 * arrival order (the ids used overall must still be an ascending, dense
	 * sequence -- see below) -- verified 2026-09-04 for every implementation
	 * reachable via WriterFactory: each is backed by OrderedQueueSystem2 or
	 * JobQueue, both of which reorder by id internally regardless of which
	 * thread submitted which id or in what order. The low-level `*ST`/`*ZT`
	 * classes constructible directly (not returned by WriterFactory) are NOT
	 * all safe for multiple producers -- see each class's own javadoc before
	 * using one outside a single-producer context.
	 * MUST USE ASCENDING NUMBERS FOR PROPER BEHAVIOUR.
	 */
	public void add(ArrayList<Read> reads, long id);
	
	/**
	 * Submit ordered batch of reads for writing.
	 * Blocks if queue is full. Thread-safe for a single producer, and also
	 * safe for multiple concurrent producer threads submitting in ANY
	 * arrival order (the ids used overall must still be an ascending, dense
	 * sequence -- see below) -- verified 2026-09-04 for every implementation
	 * reachable via WriterFactory: each is backed by OrderedQueueSystem2 or
	 * JobQueue, both of which reorder by id internally regardless of which
	 * thread submitted which id or in what order. The low-level `*ST`/`*ZT`
	 * classes constructible directly (not returned by WriterFactory) are NOT
	 * all safe for multiple producers -- see each class's own javadoc before
	 * using one outside a single-producer context.
	 * MUST USE ASCENDING NUMBERS FOR PROPER BEHAVIOUR.
	 */
	public void addReads(ListNum<Read> reads);
	
	/** 
	 * Submit ordered batch of SamLines for writing (SAM/BAM only).
	 * May throw UnsupportedOperationException for FASTA/FASTQ.
	 */
	public void addLines(ListNum<SamLine> lines);
	
	/** Signal no more data coming */
	public void poison();
	
	/** Wait for all queued data to be written
	 * @return errorState */
	public boolean waitForFinish();
	
	/** Convenience: poison and wait
	 * @return errorState */
	public boolean poisonAndWait();

	/**
	 * Force an immediate, non-blocking shutdown after an EXTERNAL failure (i.e. the caller
	 * detected a problem elsewhere in its own pipeline, not necessarily inside this Writer).
	 * Unlike poison()/poisonAndWait() -- the normal graceful-drain path -- this must return
	 * promptly and must NOT wait for any pending backlog to finish writing; whatever is still
	 * queued may be abandoned. Sets errorState(). Must be idempotent (safe to call more than
	 * once, including after finish() has already completed). 2026-09-03: added after verifying
	 * (by deliberately crashing SamWriter's worker and writer threads via stream.sh, not by
	 * inference) that every OQS2-backed Writer already recovers from an INTERNAL thread death
	 * without hanging -- this method exists for the DISTINCT case of an external caller wanting
	 * to abandon a still-healthy Writer early, which poison()'s normal drain does not address.
	 */
	public void finishError();

	public String fname();
	
//	//Can never be unset
//	public void setErrorState(boolean b);
	
	public boolean errorState();

	public boolean finishedSuccessfully();
	
}