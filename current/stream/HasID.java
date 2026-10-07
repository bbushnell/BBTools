package stream;

/**
 * Interface for jobs that can be queued and ordered by ID.
 * Queues such as {@link JobQueue} compare jobs through {@link #id()}.
 * Provides methods for job identification, poison pill detection, and completion signaling.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 23, 2025
 */
public interface HasID{

	/** Returns the ordering identifier; reusable markers can have fixed sentinel IDs. */
	public long id();
	
	/** Returns true if this is a poison pill message signaling thread shutdown */
	public boolean poison();
	
	/** Returns true if this is the last job in the sequence */
	public boolean last();

	/**
	 * Supplies a poison marker of the implementing job type.
	 * Allocation and ID handling belong to the concrete implementation: most create
	 * a new marker with the requested ID, while {@link stream.bam.BgzfInputJob} reuses
	 * a singleton with Long.MAX_VALUE and ignores the argument. Callers must use a
	 * marker policy compatible with their consumer rather than assume a fresh instance.
	 * Queue wrappers use poison markers for input-side and output-side termination;
	 * this interface does not perform either operation itself.
	 * @param id Requested marker ID, which a fixed-sentinel implementation may ignore
	 * @return A job-type marker whose poison() returns true, possibly a shared instance
	 */
	public HasID makePoison(long id);
	/**
	 * Creates a new last-job marker carrying the supplied ID.
	 * Consumers use this to mark the end of an ordered sequence; the marker itself
	 * does not close resources or perform queue operations.
	 * @param id Ordering ID for the new marker
	 * @return New job-type marker whose last() returns true and whose id() returns id
	 */
	public HasID makeLast(long id);
	
}
