package assemble;

import java.io.PrintStream;

/**
 * Holds constants for shaving.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 *
 */
public abstract class ShaveObject {
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	static PrintStream outstream=System.err;
	
	/** Assembly mode constant for standard contig assembly from input reads */
	public static final int contigMode=0;
	/** Assembly mode constant for extension of existing contig sequences */
	public static final int extendMode=1;
	/** Assembly mode constant for error correction processing operations */
	public static final int correctMode=2;
	/** Assembly mode constant for gap filling between existing contigs */
	public static final int insertMode=3;
	/** Assembly mode constant for low-quality sequence removal operations */
	public static final int discardMode=4;
	
	/** Explore codes */
	public static final int KEEP_GOING=0, DEAD_END=1, TOO_SHORT=2, TOO_LONG=3, TOO_DEEP=4, LOOP=7, SUCCESS=8;
	/** Branch codes */
	public static final int BRANCH_BIT=16, F_BRANCH=BRANCH_BIT|1, B_BRANCH=BRANCH_BIT|2, D_BRANCH=BRANCH_BIT|3;
	
	/**
	 * Tests whether the given code represents a branch type using bit mask detection.
	 * Uses BRANCH_BIT mask to identify codes that indicate graph branching conditions.
	 * @param code The result code to test for branch classification
	 * @return true if code represents any branch type (F_BRANCH, B_BRANCH, or D_BRANCH)
	 */
	public static final boolean isBranchCode(int code){return (code&BRANCH_BIT)==BRANCH_BIT;}
	
	/** Extend codes */
	public static final int BAD_OWNER=11, BAD_SEED=12/*, BRANCH=13*/;
	
	/** Traversal state constant for graph elements confirmed for retention */
	/** Traversal state constant for graph elements marked for deletion */
	/** Traversal state constant for successfully analyzed graph elements */
	/** Traversal state constant for initial unprocessed graph elements */
	public static final int STATUS_UNEXPLORED=0, STATUS_EXPLORED=1, STATUS_REMOVE=2, STATUS_KEEP=3;
	
	/**
	 * Human-readable names for all numeric result codes used in diagnostic output
	 */
	public static final String[] codeStrings=new String[] {
			"KEEP_GOING", "DEAD_END", "TOO_SHORT", "TOO_LONG", "TOO_DEEP", "5",
			"6", "LOOP", "SUCCESS", "9", "10",
			"BAD_OWNER", "BAD_SEED", "BRANCH", "14", "15",
			"BRANCH", "F_BRANCH", "B_BRANCH", "D_BRANCH"
	};
	
	/**
	 * Maximum valid code value for array bounds checking in diagnostic operations
	 */
	public static final int MAX_CODE=codeStrings.length;
	
	/**
	 * Configuration flag to enable performance monitoring output during assembly operations
	 */
	public static boolean printEventCounts=false;
	
	/** Verbose messages */
	public static boolean verbose=false;
	/** Debugging verbose messages */
	public static boolean verbose2=false;
	
}
