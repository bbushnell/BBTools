package processor;

import java.io.PrintStream;
import shared.Parser;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * Interface for read processing operations.
 * Implements the curiously recurring template pattern to ensure type-safe
 * accumulation of statistics between processor instances.
 * 
 * @param <P> The implementing class type, which must extend Processor<P>
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 16, 2025
 */
public interface Processor<P extends Processor<P>> extends Cloneable {
	
	/**
	 * Creates a deep copy of this processor.
	 * @return A new instance with copied state
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
	 * Sets fields from a Parser after it has finished parsing.
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
	 * Processes a read pair.
	 * @param r1 The first read (never null)
	 * @param r2 The second read (may be null for single-end data)
	 * @return Bit flags indicating which reads passed filters (0=none, 1=r1, 2=r2, 3=both)
	 */
	public int processReadPair(Read r1, Read r2);
	
	/**
	 * Processes a SAM line directly.
	 * @param sl The SAM line to process
	 * @return true if the line passed all filters
	 */
	public boolean processSamLine(SamLine sl);
	
	/*--------------------------------------------------------------*/
	/*----------------      Utility Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculates the recommended number of worker threads based on configured operations.
	 * @return The optimal number of workers
	 */
	public int recommendedWorkers();
	
	/**
	 * Accumulates statistics from another processor instance.
	 * Used to combine results from multiple worker threads.
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