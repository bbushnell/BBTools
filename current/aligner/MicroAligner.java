package aligner;

import stream.Read;

/**
 * Interface for lightweight alignment algorithms focused on identity calculation.
 * Defines minimal interface for read mapping and identity calculation with
 * performance optimization for high-throughput applications.
 * @author Brian Bushnell
 */
public interface MicroAligner {
	
	/** Returns identity */
	public float map(Read r);
	
	/** Returns identity */
	public float map(Read r, float minid);
	
}
