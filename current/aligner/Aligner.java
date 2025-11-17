package aligner;

/**
 * Interface for sequence alignment algorithms in BBTools.
 * Provides methods for dynamic programming-based alignment with various scoring options.
 * Implementations handle the core alignment computations and traceback operations.
 * @author Brian Bushnell
 */
public interface Aligner {

	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Will not fill areas that cannot match minScore */
	int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore);
	
	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Does not require a min score (ie, same as old method) */
	int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc);
	
	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Min score is optional */
	int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore);

	/** Generates the match string */
	byte[] traceback(byte[] query, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state);

	/** Generates identity;
	 * fills 'extra' with {match, sub, del, ins, N, clip} if present */
	float tracebackIdentity(byte[] query, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state, int[] extra);
	
	/** @return {score, bestRefStart, bestRefStop} */
	int[] score(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int maxRow, int maxCol,
			int maxState/*, final int maxScore, final int maxStart*/);

	/** Will not fill areas that cannot match minScore.
	 * @return {score, bestRefStart, bestRefStop}  */
	int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore);
	
	/** Lowest possible alignment score for a read with this length and this identity */
	int minScoreByIdentity(int len, float identity);

//	int maxRows();
//	int maxColumns();
	/** Gets the number of rows in the alignment matrix */
	int rows();
	/** Gets the number of columns in the alignment matrix */
	int columns();
	
}