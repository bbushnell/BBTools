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
	
	/** Emergency shutdown - prefer poisonAndWait() for clean termination */
	public void close();
	
	/** Number of reads written */
	public long readsWritten();
	
	/** Number of bases written */
	public long basesWritten();
	
	/** 
	 * Submit ordered batch of reads for writing.
	 * Blocks if queue is full. Thread-safe for single producer.
	 */
	public void add(ArrayList<Read> reads, long id);
	
	/** 
	 * Submit ordered batch of reads for writing.
	 * Blocks if queue is full. Thread-safe for single producer.
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
	
//	//Can never be unset
//	public void setErrorState(boolean b);
	
	public boolean errorState();

	public boolean finishedSuccessfully();
	
}