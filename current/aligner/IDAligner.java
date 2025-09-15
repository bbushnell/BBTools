package aligner;

/**
 *Interface for aligners that can calculate pairwise identity.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 24, 2025
 */
public interface IDAligner {
	
	/**
	 * @return Aligner name.
	 */
	public String name();
	
	/**
	 * @param q Query sequence
	 * @param r Reference sequence
	 * @return Identity (0.0-1.0).
	 */
	public float align(byte[] q, byte[] r);
	
	/**
	 * @param q Query sequence
	 * @param r Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * If the posVector is null, sequences may be swapped so that the query is shorter.
	 * @return Identity (0.0-1.0).
	 */
	public float align(byte[] q, byte[] r, int[] posVector);
	
	/**
	 * @param q Query sequence
	 * @param r Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * If the posVector is null, sequences may be swapped so that the query is shorter.
	 * @param rStart Alignment window start.
	 * @param rStop Alignment window stop.
	 * @return Identity (0.0-1.0).
	 */
	public float align(byte[] q, byte[] r, int[] posVector, int rStart, int rStop);
	
	/**
	 * @param q Query sequence
	 * @param r Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * If the posVector is null, sequences may be swapped so that the query is shorter.
	 * @param minScore Legacy field to allow early exit in some classes
	 * @return Identity (0.0-1.0).
	 */
	public float align(byte[] q, byte[] r, int[] posVector, int minScore);
	
	public long loops();
	public void setLoops(long i);
	
}
