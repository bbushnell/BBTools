package assemble;

import java.io.PrintStream;

/**
 * Holds shared constants for shaving/assembly operations including mode flags,
 * result codes, and branch classifications used throughout Tadpole traversal logic.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
public abstract class ShaveObject {
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print stream used for diagnostic messages during shaving operations. */
	static PrintStream outstream=System.err;
	
	/** Assembly mode for building contigs from input reads. */
	public static final int contigMode=0;
	/** Assembly mode for extending existing contigs. */
	public static final int extendMode=1;
	/** Assembly mode for error correction. */
	public static final int correctMode=2;
	/** Assembly mode for filling gaps between contigs. */
	public static final int insertMode=3;
	/** Assembly mode for discarding low-quality sequence. */
	public static final int discardMode=4;
	
	/** Exploration result code indicating the target was reached successfully. */
	/** Exploration result code indicating a circular path was detected. */
	/** Exploration result code when recursion depth is exceeded. */
	/**
	 * Exploration result code when the path exceeds the maximum length threshold.
	 */
	/**
	 * Exploration result code when the path is below the minimum length threshold.
	 */
	/** Exploration result code for a terminated path with no outgoing edges. */
	/** Exploration result code indicating traversal can continue. */
	public static final int KEEP_GOING=0, DEAD_END=1, TOO_SHORT=2, TOO_LONG=3, TOO_DEEP=4, LOOP=7, SUCCESS=8;
	/** Branch code for a dead-end branching point in the graph. */
	/** Branch code for a backward branching point in the graph. */
	/** Branch code for a forward branching point in the graph. */
	/** Bit mask used to identify branch-type result codes. */
	public static final int BRANCH_BIT=16, F_BRANCH=BRANCH_BIT|1, B_BRANCH=BRANCH_BIT|2, D_BRANCH=BRANCH_BIT|3;
	
	/**
	 * Tests whether a result code represents any branch type using the BRANCH_BIT mask.
	 * @param code Result code to test
	 * @return true if the code encodes a branch type
	 */
	public static final boolean isBranchCode(int code){return (code&BRANCH_BIT)==BRANCH_BIT;}
	
	/** Extension error code indicating an invalid seed sequence. */
	/** Extension error code indicating an ownership conflict. */
	public static final int BAD_OWNER=11, BAD_SEED=12/*, BRANCH=13*/;
	
	/** Traversal state constant for graph elements confirmed for retention */
	/** Traversal state constant for graph elements marked for deletion */
	/** Traversal state constant for successfully analyzed graph elements */
	public static final int STATUS_UNEXPLORED=0, STATUS_EXPLORED=1, STATUS_REMOVE=2, STATUS_KEEP=3;
	
	/** Human-readable names corresponding to each numeric result code. */
	public static final String[] codeStrings=new String[] {
			"KEEP_GOING", "DEAD_END", "TOO_SHORT", "TOO_LONG", "TOO_DEEP", "5",
			"6", "LOOP", "SUCCESS", "9", "10",
			"BAD_OWNER", "BAD_SEED", "BRANCH", "14", "15",
			"BRANCH", "F_BRANCH", "B_BRANCH", "D_BRANCH"
	};
	
	/** Maximum valid result code value (length of codeStrings). */
	public static final int MAX_CODE=codeStrings.length;
	
	/** Flag to enable performance monitoring output during assembly operations. */
	public static boolean printEventCounts=false;
	
	/** Enables verbose logging. */
	public static boolean verbose=false;
	/** Enables debug-level verbose logging. */
	public static boolean verbose2=false;
	
}
