package processor;

import java.io.PrintStream;

import parse.Parser;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * Read-processing operations separated from a driver's streams, batching and threads.
 * Drivers finish argument parsing and initialization before copying worker instances,
 * then accumulate statistics after workers finish. Implementations define which
 * configuration can be shared and which processing state needs an independent copy;
 * this interface alone does not make a processor safe for concurrent calls.
 * The recursive type parameter keeps accumulation within one processor type.
 * 
 * @param <P> The implementing class type, which implements {@code Processor<P>}
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 16, 2025
 */
public interface Processor<P extends Processor<P>> extends Cloneable{
	
	/**
	 * Creates a worker instance under the implementation's copying contract.
	 * This need not be a deep copy or reset accumulated counters. Consult the
	 * implementation for borrowed configuration, shared state and initialization order.
	 * @return New processor instance with implementation-defined copied state
	 */
	public P clone();
	
	/*--------------------------------------------------------------*/
	/*----------------         Parsing              ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses a command-line argument.
	 * @param arg The full argument string
	 * @param a The argument key (before '=')
	 * @param b The argument value (after '=')
	 * @return true if the argument was recognized and processed
	 */
	public boolean parse(String arg, String a, String b);
	
	/**
	 * Transfers generic options after both the driver and processor finish parsing.
	 * Implementation-specific argument handling may also update derived parser settings.
	 * @param parser The parser containing parsed values
	 */
	public void setFromParser(Parser parser);
	
	/**
	 * Finalizes settings after all parsing is complete.
	 * Called once before processing begins.
	 */
	public void postParse();
	
	/*--------------------------------------------------------------*/
	/*----------------      Processing Method       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Mutates and filters a read or linked pair, including implementation-owned statistics.
	 * The driver retains responsibility for routing selected reads to output streams.
	 * @param r1 The first read (never null)
	 * @param r2 The second read (may be null for single-end data)
	 * @return Bit flags indicating which reads passed filters (0=none, 1=r1, 2=r2, 3=both)
	 */
	public int processReadPair(Read r1, Read r2);
	
	/**
	 * Processes one SAM record under the implementation's conversion/attachment contract.
	 * This operation need not provide the paired-read context of processReadPair.
	 * @param sl The SAM line to process
	 * @return true if the line passed all filters
	 */
	public boolean processSamLine(SamLine sl);
	
	/*--------------------------------------------------------------*/
	/*----------------      Utility Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Suggests a worker count based on configured operations and shared-state constraints.
	 * The driver must also enforce required serial-processing policies when applicable.
	 * @return Recommended worker count, not a measured throughput guarantee
	 */
	public int recommendedWorkers();
	
	/**
	 * Accumulates statistics from another processor instance.
	 * Call after the source worker finishes; this is not a worker-completion barrier.
	 * Callers must avoid merging copied initial totals or a worker's results twice.
	 * @param other The processor whose statistics should be added to this one
	 */
	public void add(P other);

	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints processing statistics to the specified stream.
	 * @param stream The output stream for statistics
	 */
	public void printStats(PrintStream stream);

	/**
	 * Generates a ByteBuilder containing formatted statistics.
	 * @return A ByteBuilder with statistics text
	 */
	public ByteBuilder toStats();
	
}
