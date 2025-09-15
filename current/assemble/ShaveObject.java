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
	
	public static final int contigMode=0;
	public static final int extendMode=1;
	public static final int correctMode=2;
	public static final int insertMode=3;
	public static final int discardMode=4;
	
	/** Explore codes */
	public static final int KEEP_GOING=0, DEAD_END=1, TOO_SHORT=2, TOO_LONG=3, TOO_DEEP=4, LOOP=7, SUCCESS=8;
	/** Branch codes */
	public static final int BRANCH_BIT=16, F_BRANCH=BRANCH_BIT|1, B_BRANCH=BRANCH_BIT|2, D_BRANCH=BRANCH_BIT|3;
	
	public static final boolean isBranchCode(int code){return (code&BRANCH_BIT)==BRANCH_BIT;}
	
	/** Extend codes */
	public static final int BAD_OWNER=11, BAD_SEED=12/*, BRANCH=13*/;
	
	public static final int STATUS_UNEXPLORED=0, STATUS_EXPLORED=1, STATUS_REMOVE=2, STATUS_KEEP=3;
	
	public static final String[] codeStrings=new String[] {
			"KEEP_GOING", "DEAD_END", "TOO_SHORT", "TOO_LONG", "TOO_DEEP", "5",
			"6", "LOOP", "SUCCESS", "9", "10",
			"BAD_OWNER", "BAD_SEED", "BRANCH", "14", "15",
			"BRANCH", "F_BRANCH", "B_BRANCH", "D_BRANCH"
	};
	
	public static final int MAX_CODE=codeStrings.length;
	
	public static boolean printEventCounts=false;
	
	/** Verbose messages */
	public static boolean verbose=false;
	/** Debugging verbose messages */
	public static boolean verbose2=false;
	
}
