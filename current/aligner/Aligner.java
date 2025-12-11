package aligner;

/**
 * Interface for sequence alignment algorithms in BBTools.
 * Provides methods for dynamic programming-based alignment with various scoring options.
 * Implementations handle the core alignment computations and traceback operations.
 * @author Brian Bushnell
 */
public interface Aligner {

	/**
	 * Fills alignment matrix with minimum score constraint for efficiency.
	 * Skips areas that cannot achieve the minimum score to reduce computation.
	 *
	 * @param read Query sequence bytes
	 * @param ref Reference sequence bytes
	 * @param refStartLoc Starting position in reference
	 * @param refEndLoc Ending position in reference
	 * @param minScore Minimum alignment score threshold
	 * @return Array containing [rows, maxCol, maxState, maxScore, maxStart]
	 */
	int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore);
	
	/**
	 * Fills complete alignment matrix without score constraints.
	 * Equivalent to the original unlimited fill method.
	 *
	 * @param read Query sequence bytes
	 * @param ref Reference sequence bytes
	 * @param refStartLoc Starting position in reference
	 * @param refEndLoc Ending position in reference
	 * @return Array containing [rows, maxCol, maxState, maxScore, maxStart]
	 */
	int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc);
	
	/**
	 * Fills complete alignment matrix with optional minimum score parameter.
	 * Provides flexibility to specify minimum score or ignore it.
	 *
	 * @param read Query sequence bytes
	 * @param ref Reference sequence bytes
	 * @param refStartLoc Starting position in reference
	 * @param refEndLoc Ending position in reference
	 * @param minScore Optional minimum score threshold
	 * @return Array containing [rows, maxCol, maxState, maxScore, maxStart]
	 */
	int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore);

	/**
	 * Generates the alignment match string by tracing back through the matrix.
	 * Reconstructs the optimal alignment path from the endpoint coordinates.
	 *
	 * @param query Query sequence bytes
	 * @param ref Reference sequence bytes
	 * @param refStartLoc Starting position in reference
	 * @param refEndLoc Ending position in reference
	 * @param row Matrix row of alignment endpoint
	 * @param col Matrix column of alignment endpoint
	 * @param state Alignment state at endpoint
	 * @return Byte array representing the alignment match string
	 */
	byte[] traceback(byte[] query, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state);

	/**
	 * Computes alignment identity and detailed statistics via traceback.
	 * Populates extra array with alignment event counts if provided.
	 *
	 * @param query Query sequence bytes
	 * @param ref Reference sequence bytes
	 * @param refStartLoc Starting position in reference
	 * @param refEndLoc Ending position in reference
	 * @param row Matrix row of alignment endpoint
	 * @param col Matrix column of alignment endpoint
	 * @param state Alignment state at endpoint
	 * @param extra Output array for [match, sub, del, ins, N, clip] counts (if present)
	 * @return Alignment identity as a fraction
	 */
	float tracebackIdentity(byte[] query, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state, int[] extra);
	
	/**
	 * Computes alignment score from matrix coordinates.
	 * Returns score and reference boundaries of the optimal alignment.
	 *
	 * @param read Query sequence bytes
	 * @param ref Reference sequence bytes
	 * @param refStartLoc Starting position in reference
	 * @param refEndLoc Ending position in reference
	 * @param maxRow Matrix row with maximum score
	 * @param maxCol Matrix column with maximum score
	 * @param maxState Alignment state with maximum score
	 * @return Array containing [score, bestRefStart, bestRefStop]
	 */
	int[] score(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int maxRow, int maxCol,
			int maxState/*, final int maxScore, final int maxStart*/);

	/**
	 * Combined fill and scoring operation with minimum score constraint.
	 * Efficiently computes alignment while skipping low-scoring regions.
	 *
	 * @param read Query sequence bytes
	 * @param ref Reference sequence bytes
	 * @param refStartLoc Starting position in reference
	 * @param refEndLoc Ending position in reference
	 * @param minScore Minimum alignment score threshold
	 * @return Array containing [score, bestRefStart, bestRefStop]
	 */
	int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore);
	
	/**
	 * Calculates theoretical minimum alignment score for given length and identity.
	 * Used to determine score thresholds for alignment filtering.
	 *
	 * @param len Read length in bases
	 * @param identity Target alignment identity fraction
	 * @return Lowest possible alignment score achieving the specified identity
	 */
	int minScoreByIdentity(int len, float identity);

//	int maxRows();
//	int maxColumns();
	/**
	 * Returns the number of rows in the alignment matrix (typically read length + padding).
	 * Useful for diagnostics and buffer sizing.
	 * @return Row count in the alignment matrix
	 */
	int rows();
	/**
	 * Returns the number of columns in the alignment matrix (typically reference window length).
	 * Useful for diagnostics and buffer sizing.
	 * @return Column count in the alignment matrix
	 */
	int columns();
	
}