package aligner;

import java.util.Arrays;

import shared.Tools;
import structures.IntList;

/**
 * Presentation-only simplified version of the quantum alignment algorithm.
 * This implementation is decomposed into readable functions for publication
 * and educational purposes. For production use, see {@link QuantumAligner}.
 * Uses sparse dynamic programming with bandwidth optimization and quantum
 * teleportation across deletion gaps.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date April 24, 2025
 */
public class QuantumAlignerConcise2 implements IDAligner{

	/**
	 * Program entry point that delegates to Test framework.
	 * Uses reflection to avoid redundant testing code.
	 * @param args Command-line arguments for alignment testing
	 * @throws Exception If class cannot be instantiated or testing fails
	 */
	public static <C extends IDAligner> void main(String[] args) throws Exception {
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public QuantumAlignerConcise2() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the name identifier for this aligner implementation */
	@Override
	public final String name() {return "QuantumConcise2";}
	/**
	 * Aligns two sequences and returns identity score.
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and returns identity score with position information.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional array to store reference start/stop coordinates [rStart, rStop]
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold.
	 * Note: minScore parameter is ignored in this implementation.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional array to store reference coordinates
	 * @param minScore Minimum score threshold (ignored)
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns query to a window of the reference sequence.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional array to store reference coordinates
	 * @param rStart Reference window start coordinate
	 * @param rStop Reference window end coordinate
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculates optimal bandwidth for sparse alignment based on sequence similarity.
	 * Tests initial bases for mismatches to estimate required search width.
	 * Uses adaptive bandwidth with minimum of 2 and maximum based on sequence length.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @return Bandwidth value for sparse alignment
	 */
	private static int decideBandwidth(byte[] query, byte[] ref) {
		int subs=0, bandwidth=Tools.min(query.length/4+2, Math.max(query.length, ref.length)/32, 12);
		bandwidth=Math.max(2, bandwidth);
		for(int i=0, minlen=Math.min(query.length, ref.length); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);}
		return Math.min(subs+1, bandwidth);//At least 1
	}

	/**
	 * Core alignment algorithm using sparse dynamic programming with quantum teleportation.
	 * Implements bandwidth-limited alignment with score-based pruning and bridge building
	 * across deletion gaps. Uses bit-packed scores to track position, deletions, and score.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional array to store reference coordinates [rStart, rStop]
	 * @return Identity score between 0.0 and 1.0
	 */
	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector) {
		final int qLen=query.length, rLen=ref.length;
		final int addRight=2;// Cells to add to the right of score band
		final int bandWidth=decideBandwidth(query, ref);
		final int topBand=2*bandWidth;// Fully explored top rows
		final long scoreBand=(bandWidth+1L)<<SCORE_SHIFT;// Low score pruning cutoff
		long[] prev=new long[rLen+1], curr=new long[rLen+1];// Score arrays
		for(int j=0; j<=rLen; j++){prev[j]=j;}// Initialize prev scores to position
		// Create IntLists for tracking active positions
		IntList activeList=new IntList(rLen+addRight), nextList=new IntList(rLen+addRight);
		for(int j=0; j<=rLen; j++) {activeList.add(j);}	// Initialize active list
		int maxPos=0; // Best scoring position
		long maxScore=BAD, prevRowScore=BAD; //Initialize score outside of loop
		for(int i=1; i<=qLen; i++){// Fill alignment matrix using the sparse loop
			prevRowScore=maxScore; maxScore=BAD; maxPos=0;// Swap row best scores
			maxPos=processRow(i, rLen, addRight, query[i-1], ref, prev, curr, maxScore, 
					maxPos, prevRowScore, activeList, nextList, scoreBand, topBand);
			maxScore=curr[maxPos];// Grab maxScore from the processed row 
			long[] temp=prev; prev=curr; curr=temp;// Swap rows
			IntList tempL=activeList; activeList=nextList; nextList=tempL;// Swap lists
		}
		return postprocess(maxScore, maxPos, qLen, rLen, posVector);
	}
	
	/**
	 * Processes a single row of the alignment matrix using sparse computation.
	 * Only computes cells within the active bandwidth, with optional bridge building
	 * to catch long deletions. Updates active position lists for next iteration.
	 *
	 * @param i Current row index (1-based)
	 * @param rLen Reference sequence length
	 * @param addRight Additional cells to explore to the right
	 * @param q Query base at position i-1
	 * @param ref Reference sequence array
	 * @param prev Previous row scores
	 * @param curr Current row scores
	 * @param maxScore Maximum score in current row
	 * @param maxPos Position of maximum score
	 * @param prevRowScore Maximum score from previous row
	 * @param activeList Current active positions
	 * @param nextList Next row active positions
	 * @param scoreBand Score difference threshold for pruning
	 * @param topBand Number of fully explored top rows
	 * @return Position of highest score in processed row
	 */
	private static int processRow(int i, int rLen, int addRight, byte q, byte[] ref, 
			long[] prev, long[] curr, long maxScore, int maxPos, long prevRowScore, 
			IntList activeList, IntList nextList, long scoreBand, int topBand) {
		curr[0]=i*INS;// First column
		while(activeList.lastElement()>rLen) {activeList.pop();}// Remove any excess sites
		if(BUILD_BRIDGES){//Optionally race to catch up with long deletions
			int extra=((i&15)<1 ? 40 : 5), last=activeList.lastElement();
			int eLimit=Math.min(last+extra, rLen);//Extra horizontal cells to explore
			for(int e=last+1; e<eLimit; e++) {activeList.add(e);}
		}
		if(activeList.lastElement()<rLen) {prev[rLen]=BAD;}// Clear potential stale value
		nextList.clear().add(0);// Nonempty lists simplify logic
		for(int idx=1; idx<activeList.size; idx++){// Process only active positions
			final int j=activeList.get(idx);
			final long maxValue=processCell(i, j, rLen, addRight, q, ref[j-1], prev, curr, 
					maxScore, maxPos, prevRowScore, nextList, scoreBand, topBand);
			final boolean better=((maxValue&SCORE_MASK)>maxScore);
			maxScore=better ? maxValue : maxScore;// Track best score in row
			maxPos=better ? j : maxPos;
		}
		return maxPos;
	}
	
	/**
	 * Processes a single cell of the alignment matrix with quantum teleportation.
	 * Computes match/substitution, insertion, and deletion scores using bit-packed values.
	 * Implements quantum teleportation by allowing jumps across unexplored deletion gaps.
	 * Updates active position list based on score thresholds and match extension.
	 *
	 * @param i Current row index
	 * @param j Current column index
	 * @param rLen Reference sequence length
	 * @param addRight Additional positions to add to active list
	 * @param q Query base at position i-1
	 * @param r Reference base at position j-1
	 * @param prev Previous row scores
	 * @param curr Current row scores
	 * @param maxScore Maximum score in current row
	 * @param maxPos Position of maximum score
	 * @param prevRowScore Maximum score from previous row
	 * @param nextList Active position list for next row
	 * @param scoreBand Score threshold for position pruning
	 * @param topBand Number of fully explored rows
	 * @return Computed score for current cell
	 */
	private static long processCell(int i, int j, int rLen, int addRight, byte q, byte r, 
			long[] prev, long[] curr, long maxScore, int maxPos, long prevRowScore, 
			IntList nextList, long scoreBand, int topBand) {
		// Branchless score calculation. q and r are bases at query[i-1] and ref[j-1].
		final boolean isMatch=(q==r && q!='N');
		final boolean hasN=(q=='N' || r=='N');
		final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);// +1 match, 0 N, -1 sub
		final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];// Read adjacent scores
		final long diagScore=pj1+scoreAdd;// Match/Sub
		final long upScore=pj+INS;// Insertion
		final long leftScore1=cj1+DEL_INCREMENT; //Deletion; adjust both score and del counter
		// Allows quantum teleport across the unexplored gap from a long deletion
		final long leftScore=Math.max(leftScore1, maxScore+DEL_INCREMENT*(j-maxPos));
		// Find max using conditional expressions
		final long maxDiagUp=Math.max(diagScore, upScore);
		final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
		final long scoreDif=prevRowScore-maxValue;// Determines whether j is within score band
		final int lastPositionAdded=nextList.array[nextList.size-1];
		final boolean add=j<=rLen && (scoreDif<scoreBand || i<topBand);
		final boolean live=(EXTEND_MATCH && isMatch & lastPositionAdded<j+1);
		curr[j]=(add || live ? maxValue : BAD);// Update or prune current cell
		prev[j-1]=BAD;//Clear previous row
		if(add) {// Conditionally add multiple positions
			final int from=Math.max(lastPositionAdded+1, j);
			final int to=Math.min(j+addRight, rLen);
			for(int k=from; k<=to; k++) {nextList.addUnchecked(k);}
		}else if(live) {nextList.addUnchecked(j+1);}// Extend diagonal only
		return maxValue;
	}
	
	/**
	 * Extracts alignment statistics from bit-packed score and calculates identity.
	 * Solves system of equations to determine match, substitution, insertion, and
	 * deletion counts from the final alignment score and position information.
	 *
	 * @param maxScore Highest bit-packed score from final row
	 * @param maxPos Reference position of highest score
	 * @param qLen Query sequence length
	 * @param rLen Reference sequence length
	 * @param posVector Optional array to store reference coordinates
	 * @return Identity score as fraction of aligned bases that match
	 */
	private static float postprocess(long maxScore, int maxPos, int qLen, int rLen, int[] posVector) {
		// Extract alignment information
		final int originPos=(int)(maxScore&POSITION_MASK);
		final int deletions=(int)((maxScore&DEL_MASK) >> POSITION_BITS);
		final int refAlnLength=(maxPos-originPos);
		final int rawScore=(int)(maxScore >> SCORE_SHIFT);
		// Solve the system of equations to calculate operation counts:
		// 1. M + S + I = qLen
		// 2. M + S + D = refAlnLength
		// 3. M - S - I - D = Score
		final int insertions=Math.max(0, qLen+deletions-refAlnLength);
		final float matches=((rawScore+qLen+deletions)/2f);
		final float substitutions=Math.max(0, qLen-matches-insertions);
		//Generate results
		final float identity=matches/(matches+substitutions+insertions+deletions);
		if(posVector!=null){posVector[0]=originPos; posVector[1]=maxPos-1;} //Store rstart, rstop
		return identity;
	}

	/**
	 * Aligns query to a specified window of the reference sequence.
	 * Creates a subarray of the reference for the alignment window and adjusts
	 * coordinates in the position vector to reflect the original reference.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional array to store reference coordinates
	 * @param refStart Start coordinate of alignment window (inclusive)
	 * @param refEnd End coordinate of alignment window (inclusive)
	 * @return Identity score between 0.0 and 1.0
	 */
	public static final float alignStatic(final byte[] query, final byte[] ref, 
			final int[] posVector, int refStart, int refEnd) {
		refStart=Math.max(refStart, 0);
		refEnd=Math.min(refEnd, ref.length-1);
		final int rlen=refEnd-refStart+1;
		final byte[] region=(rlen==ref.length ? ref : Arrays.copyOfRange(ref, refStart, refEnd));
		final float id=alignStatic(query, region, posVector);
		if(posVector!=null) {
			posVector[0]+=refStart;
			posVector[1]+=refStart;
		}
		return id;
	}
	
	public long loops() {return -1;}
	public void setLoops(long x) {}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Bit field definitions
	private static final int POSITION_BITS=21;
	private static final int DEL_BITS=21;
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;
	// Masks
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	private static final long DEL_MASK=((1L << DEL_BITS)-1) << POSITION_BITS;
	private static final long SCORE_MASK=~(POSITION_MASK | DEL_MASK);
	// Scoring constants
	private static final long MATCH=1L << SCORE_SHIFT;
	private static final long SUB=(-1L) << SCORE_SHIFT;
	private static final long INS=(-1L) << SCORE_SHIFT;
	private static final long DEL=(-1L) << SCORE_SHIFT;
	private static final long N_SCORE=0L;
	private static final long BAD=Long.MIN_VALUE/2;
	private static final long DEL_INCREMENT=(1L<<POSITION_BITS)+DEL;
	// Run modes
	private static final boolean EXTEND_MATCH=true;
	private static final boolean BUILD_BRIDGES=true;

}
