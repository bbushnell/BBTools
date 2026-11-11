package assemble;

import java.io.PrintStream;

/**
 * Holds shared constants for shaving/assembly operations including mode flags,
 * result codes, and branch classifications used throughout Tadpole traversal logic.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
public abstract class ShaveObject{
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print stream used for diagnostic messages during shaving operations. */
	static PrintStream outstream=System.err;
	
	/** Assembly mode for building contigs from input reads. */
	public static final int contigMode=0;
	/** Assembly mode for extending reads or sequences. */
	public static final int extendMode=1;
	/** Assembly mode for error correction. */
	public static final int correctMode=2;
	/** Assembly mode for estimating insert size from paired reads. */
	public static final int insertMode=3;
	/** Assembly mode for configured read filtering. */
	public static final int discardMode=4;
	
	/** Exploration result code indicating traversal can continue. */
	public static final int KEEP_GOING=0;
	/** Exploration result code for a terminated path with no outgoing edges. */
	public static final int DEAD_END=1;
	/** Exploration result code when the path is below the minimum length threshold. */
	public static final int TOO_SHORT=2;
	/** Exploration result code when the path exceeds the maximum length threshold. */
	public static final int TOO_LONG=3;
	/** Exploration result code when k-mer depth exceeds maxCount. */
	public static final int TOO_DEEP=4;
	/** Exploration result code indicating a circular path was detected. */
	public static final int LOOP=7;
	/** Exploration result code indicating the target was reached successfully. */
	public static final int SUCCESS=8;

	/** Bit mask used to identify branch-type result codes. */
	public static final int BRANCH_BIT=16;
	/** Branch code for a forward branching point in the graph. */
	public static final int F_BRANCH=BRANCH_BIT|1;
	/** Branch code for a backward branching point in the graph. */
	public static final int B_BRANCH=BRANCH_BIT|2;
	/** Branch code for a point branching in both directions in the graph. */
	public static final int D_BRANCH=BRANCH_BIT|3;
	
	/**
	 * Tests whether a result code represents any branch type using the BRANCH_BIT mask.
	 * @param code Result code to test
	 * @return true if the code encodes a branch type
	 */
	public static final boolean isBranchCode(int code){return (code&BRANCH_BIT)==BRANCH_BIT;}
	
	/** Extension error code indicating an ownership conflict. */
	public static final int BAD_OWNER=11;
	/** Extension error code indicating an invalid seed sequence. */
	public static final int BAD_SEED=12/*, BRANCH=13*/;
	
	/** Traversal state for graph elements not yet analyzed. */
	public static final int STATUS_UNEXPLORED=0;
	/** Traversal state for successfully analyzed graph elements. */
	public static final int STATUS_EXPLORED=1;
	/** Traversal state for graph elements marked for deletion. */
	public static final int STATUS_REMOVE=2;
	/** Traversal state for graph elements confirmed for retention. */
	public static final int STATUS_KEEP=3;
	
	/** Human-readable names corresponding to each numeric result code. */
	public static final String[] codeStrings=new String[] {
			"KEEP_GOING", "DEAD_END", "TOO_SHORT", "TOO_LONG", "TOO_DEEP", "5",
			"6", "LOOP", "SUCCESS", "9", "10",
			"BAD_OWNER", "BAD_SEED", "BRANCH", "14", "15",
			"BRANCH", "F_BRANCH", "B_BRANCH", "D_BRANCH"
	};
	
	/** Exclusive upper bound for valid indexes into codeStrings. */
	public static final int MAX_CODE=codeStrings.length;
	
	/** Flag to enable diagnostic event-count output during assembly operations. */
	public static boolean printEventCounts=false;
	
	/** Enables verbose logging. */
	public static boolean verbose=false;
	/** Enables debug-level verbose logging. */
	public static boolean verbose2=false;
	
}
