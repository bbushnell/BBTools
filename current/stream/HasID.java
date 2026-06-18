package stream;

/**
 * Interface for jobs that can be queued and ordered by ID.
 * Extends Comparable to support priority queue ordering.
 * Provides methods for job identification, poison pill detection, and completion signaling.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 23, 2025
 */
public interface HasID {

	/** Returns unique identifier for this job, used for ordering */
	public long id();
	
	/** Returns true if this is a poison pill message signaling thread shutdown */
	public boolean poison();
	
	/** Returns true if this is the last job in the sequence */
	public boolean last();

	/** Factory: returns a new poison-pill instance of this job type carrying the given id. Serves double duty in OrderedQueueSystem: poison() puts one on the INPUT queue to shut workers down, and setFinished() puts one on the OUTPUT queue to terminate the consumer. */
	public HasID makePoison(long id);
	/** Factory: returns a new 'last' marker instance of this job type carrying the given id. Used on the OUTPUT side to mark end-of-ordered-sequence (see OrderedQueueSystem.poison). */
	public HasID makeLast(long id);
	
}