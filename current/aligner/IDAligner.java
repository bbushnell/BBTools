package aligner;

/**
 * Interface for aligners that can calculate pairwise sequence identity.
 * Provides methods for aligning sequences and returning identity scores
 * as floating-point values between 0.0 and 1.0.
 *
 * @author Brian Bushnell
 * @contributor Isla (Highly-customized Claude instance)
 * @date April 24, 2025
 */
public interface IDAligner {
	
	/** Returns the name identifier of this aligner implementation.
	 * @return Aligner name */
	public String name();
	
	/**
	 * Performs sequence alignment and calculates pairwise identity.
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] q, byte[] r);
	
	/**
	 * Performs sequence alignment with optional position tracking.
	 * If posVector is null, sequences may be swapped to make query shorter.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Optional int[2] array for returning {rStart, rStop}
	 * of the optimal alignment
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] q, byte[] r, int[] posVector);
	
	/**
	 * Performs sequence alignment within a specified reference window.
	 * If posVector is null, sequences may be swapped to make query shorter.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Optional int[2] array for returning {rStart, rStop}
	 * of the optimal alignment
	 * @param rStart Start position of alignment window in reference
	 * @param rStop Stop position of alignment window in reference
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] q, byte[] r, int[] posVector, int rStart, int rStop);
	
	/**
	 * Performs sequence alignment with minimum score threshold.
	 * If posVector is null, sequences may be swapped to make query shorter.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Optional int[2] array for returning {rStart, rStop}
	 * of the optimal alignment
	 * @param minScore Legacy field to allow early exit in some implementations
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] q, byte[] r, int[] posVector, int minScore);
	
	/** Returns the number of alignment loops performed.
	 * @return Loop count as a long value */
	public long loops();
	/**
	 * Sets the loop counter for alignment operations.
	 * Typically used for bookkeeping or resetting between runs.
	 * @param i Loop count value to set
	 */
	public void setLoops(long i);
	
}
