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

	public HasID makePoison(long id);
	public HasID makeLast(long id);
	
}