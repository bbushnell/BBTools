package stream;

import structures.ListNum;

public interface ConcurrentReadStreamInterface extends Runnable{
	
	/** Start this in a new thread. */
	public void start();
	
	/** Fetch the next list of reads.  Returns an empty list when end of input is reached. */
	public ListNum<Read> nextList();
	
	/** Re-calls returnList(returnList(ln.id, ln.list.isEmpty()) */
	public void returnList(ListNum<Read> ln);
	
	/** When the nextList() caller is done processing a list, it MUST be returned using this method.
	 * The 'poison' flag should be set to false normally.  When a consumer thread receives an empty list from nextList(),
	 * it should be returned with the poison flag set to true, then the consumer should terminate.
	 * This will return a list to the 'full' queue, allowing another thread to pull the empty list and terminate.  */
	public void returnList(long listNumber, boolean poison);
	
	/** This must be called (indirectly, via Thread.start()) before reads will be generated. */
	@Override
	public void run();
	
	/** Indicate to producer threads that no more reads are desired, and interrupts them. */
	public void shutdown();

	/** Reset state to allow production of reads from the beginning of the input files.
	 * Does not work with stdin (may cause strange behavior). */
	public void restart();
	
	/** Calls shutdown, then shuts down all threads and closes all associated files. */
	public void close();
	
	/** Returns true for paired-end stream, false for single-end stream. */
	public boolean paired();
	
	/** Name of the file or resource underlying this stream */
	public String fname();
	
	/** Returns the underlying read object producer(s), such as ReadInputStreams.  Optional method for things such as error messages. */
	public Object[] producers();
	
	/** Return true if this stream or its producers have detected an error. */
	public boolean errorState();
	
	/**
	 * Set the read sampling rate.  Optional method.
	 * @param rate Fraction of reads to use, 0-1.
	 * @param seed Random number generator seed when positive.  If negative, a random seed will be used.
	 */
	public void setSampleRate(float rate, long seed);
	
	/**
	 * @return Number of bases read by this stream.
	 */
	public long basesIn();
	
	/**
	 * @return Number of reads read by this stream.
	 */
	public long readsIn();
	
	/**
	 * @return Value of verbose field.
	 */
	public boolean verbose();
	
//	public String classname(); //e.g. getClass().getName()
	
}
