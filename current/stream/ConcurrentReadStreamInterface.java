package stream;

import structures.ListNum;

/**
 * Legacy producer/consumer API for numbered read batches, implemented by the cris family.
 * Consumers obtain batches through nextList and return spent batches through returnList.
 * Concrete implementations define ordering, terminal representation, supported consumer
 * concurrency, and lifecycle synchronization; this interface does not supply those mechanisms.
 * Predates the newer {@link Streamer} interface. {@link ConcurrentReadInputStream}
 * provides the common abstract implementation and an adapter for returning ListNum objects.
 *
 * @author Brian Bushnell
 */
public interface ConcurrentReadStreamInterface extends Runnable{
	
	/** Starts input processing according to the implementation's lifecycle and thread policy. */
	public void start();
	
	/**
	 * Obtains the next numbered batch, potentially waiting for input to become available.
	 * Follow the implementation's null/empty/terminal convention and return each supplied
	 * batch as required by its ownership protocol. No universal ordering policy is imposed here.
	 * @return Next batch or implementation-defined terminal result
	 */
	public ListNum<Read> nextList();
	
	/**
	 * Returns a batch the consumer has finished using so its storage may be reused.
	 * Consumers must not rely on retained mutable list contents after returning the batch.
	 * ConcurrentReadInputStream ignores null and forwards a nonnull batch's ID with
	 * ListNum.isEmpty() as the terminal flag; other implementations define their own adapter.
	 * @param ln Batch to return, including terminal batches according to the implementation
	 */
	public void returnList(ListNum<Read> ln);
	
	/**
	 * Returns a numbered batch or batch slot under the implementation's ownership protocol.
	 * @param listNumber Identifier of the batch being returned
	 * @param poison Whether this return represents a terminal signal under that protocol
	 */
	public void returnList(long listNumber, boolean poison);
	
	/** Executes stream processing through Runnable; ordinary callers initiate processing with start(). */
	@Override
	public void run();
	
	/**
	 * Requests shutdown. Handling of pending batches, resource release and completion
	 * waiting belongs to the concrete implementation.
	 */
	public void shutdown();

	/**
	 * Requests state reset for reading again from the beginning of a repeatable source.
	 * Reset/start sequencing is implementation-specific. Standard input is not a
	 * supported rewindable source for this legacy contract.
	 */
	public void restart();
	
	/**
	 * Closes the stream through its concrete lifecycle; completion waiting and pending
	 * batch handling are implementation-specific. Query errorState for reported status:
	 * ReadWrite's closure helpers call it after this void-returning method.
	 */
	public void close();
	
	/** @return Whether the implementation reports paired-read processing */
	public boolean paired();
	
	/** @return Input name or stream identifier supplied by the implementation */
	public String fname();
	
	/**
	 * Exposes producer objects associated with the stream.
	 * @return Implementation-defined array; element types, count and ownership are not specified here
	 */
	public Object[] producers();
	
	/**
	 * Reports the implementation's current error status, including any sources it consults.
	 * Closure helpers query this after close; this interface does not define status storage.
	 * @return True when the implementation reports an error
	 */
	public boolean errorState();
	
	/**
	 * Configures input subsampling. Sampling units, timing of configuration changes and
	 * seed interpretation depend on the concrete stream.
	 * @param rate Requested sampling fraction between 0 and 1
	 * @param seed Random seed; special values and reproducibility rules are implementation-specific
	 */
	public void setSampleRate(float rate, long seed);
	
	/**
	 * Reports a base count at the implementation's counting stage.
	 * Snapshot visibility and when totals become final are implementation-specific.
	 * @return Reported number of input bases
	 */
	public long basesIn();
	
	/**
	 * Reports a read count using the implementation's units and counting stage.
	 * Treatment of mates versus fragments and snapshot/final-total semantics are not
	 * standardized by this interface; consult the concrete stream's contract.
	 * @return Reported input-read count in the implementation's units
	 */
	public long readsIn();
	
	/** @return Whether the implementation reports verbose diagnostics enabled */
	public boolean verbose();
	
}
