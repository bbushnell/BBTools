package idaligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Glocal aligner optimized for computing approximate nucleotide identity (ANI).
 * Uses integer arithmetic and avoids traceback to maximize performance.
 * Implements a space-efficient dynamic programming algorithm using only two arrays.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date April 19, 2024
 */
public class GlocalAlignerInt implements IDAligner{

	/**
	 * Program entry point that delegates to Test class for standardized testing.
	 * Uses reflection to determine the calling class automatically.
	 * @param args Command-line arguments
	 * @throws Exception If reflection fails or testing encounters errors
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

	public GlocalAlignerInt() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the aligner name identifier */
	@Override
	public final String name() {return "GlocalInt";}
	/**
	 * Aligns two sequences and returns identity score.
	 * @param a First sequence
	 * @param b Second sequence
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and returns identity score with position information.
	 *
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional array to store alignment start/stop positions
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold.
	 * Note: minScore parameter is currently ignored in this implementation.
	 *
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional array to store alignment start/stop positions
	 * @param minScore Minimum score threshold (ignored)
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns sequences within a specified reference window.
	 *
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional array to store alignment start/stop positions
	 * @param rStart Start position of reference window
	 * @param rStop Stop position of reference window
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Core alignment method using dynamic programming to compute ANI.
	 * Implements space-efficient algorithm with integer scoring and position tracking.
	 * Automatically swaps sequences when posVector is null to ensure query <= ref length.
	 *
	 * @param query Query sequence (may be swapped with ref if shorter)
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of optimal alignment
	 * @return Identity score between 0.0 and 1.0
	 */
	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector) {
		// Swap to ensure query is not longer than ref
		if(posVector==null && query.length>ref.length) {
			byte[] temp=query;
			query=ref;
			ref=temp;
		}

		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;

		// Create arrays for current and previous rows
		int[] prev=new int[rLen+1], curr=new int[rLen+1];

		// Initialize first row with starting position in the lower bits
		for(int j=0; j<=rLen; j++) {prev[j]=j;}

		// Fill alignment matrix
		for(int i=1; i<=qLen; i++) {
			// First column-gap in reference
			curr[0]=i*INS;
			mloops+=rLen;

			for(int j=1; j<=rLen; j++) {
				byte q=query[i-1];
				byte r=ref[j-1];

				// Branchless score calculation
				boolean isMatch=(q==r);
				boolean hasN=(q=='N' || r=='N');
				int scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : MISMATCH);

				// Cache array accesses
				final int pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				int diagScore=pj1+scoreAdd;// Match/Sub
				int upScore=pj+INS;
				int leftScore=cj1+DEL;

				// Find max using conditional expressions
				int maxDiagUp=diagScore >= upScore ? diagScore : upScore;
				int maxValue=(maxDiagUp & SCORE_MASK) >= leftScore ? maxDiagUp : leftScore;

				curr[j]=maxValue;
			}

			// Swap rows
			int[] temp=prev;
			prev=curr;
			curr=temp;
		}
		
		// Find best score outside of main loop
		int maxScore=Integer.MIN_VALUE;
		int maxPos=0;
		for(int j=1; j<=rLen; j++){
		    int score=prev[j] & SCORE_MASK;
		    if(score>maxScore){
		        maxScore=score;
		        maxPos=j;
		    }
		}
		
		// Extract alignment information
		final int bestScore=prev[maxPos];
		final int originPos=bestScore & POSITION_MASK;
		final int endPos=maxPos;
		if(posVector!=null){
		    posVector[0]=originPos;
		    posVector[1]=endPos-1;
		}

		// Calculate alignment statistics
		final int refAlnLength=(endPos-originPos);
		final int rawScore=bestScore >> SCORE_SHIFT;

		// Calculate net gaps
		final int netGaps=Math.abs(qLen-refAlnLength);
		final float matches, insertions, deletions;

		// Apply the formulas we derived
		if(qLen>refAlnLength){
		    // More insertions than deletions case
		    matches=(rawScore+qLen)/2f;
		    insertions=netGaps;
		    deletions=0;
		}else{
		    // More deletions than insertions case
		    matches=(rawScore+refAlnLength)/2f;
		    insertions=0;
		    deletions=netGaps;
		}

		// Calculate mismatches directly
		float mismatches=Math.min(qLen, refAlnLength)-matches;
		mismatches=Math.max(mismatches, 0);

		// Calculate identity
		return matches/(matches+mismatches+insertions+deletions);
	}
	
	/**
	 * Aligns query to a bounded region of the reference by copying that slice
	 * and delegating to the core aligner. Returned positions (when provided) are
	 * offset by refStart to map back to the original reference coordinates.
	 * refStart/refEnd are clamped to valid indices; refEnd is treated as an
	 * exclusive upper bound for the copied slice.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of optimal alignment
	 * @param refStart Start position of alignment window
	 * @param refEnd End position of alignment window (exclusive for copying)
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
	
	private static AtomicLong loops=new AtomicLong(0);
	public long loops() {return loops.get();}
	public void setLoops(long x) {loops.set(x);}
	public static String output=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Position will use the lower 15 bits (sufficient for 32kbp)
	private static final int POSITION_BITS=15;
	private static final int POSITION_MASK=(1 << POSITION_BITS)-1;
	private static final int SCORE_MASK=~POSITION_MASK;
	private static final int SCORE_SHIFT=POSITION_BITS;

	// Equal weighting for operations
	private static final int MATCH=1 << SCORE_SHIFT;
	private static final int MISMATCH=(-1)*(1 << SCORE_SHIFT);
	private static final int INS=(-1)*(1 << SCORE_SHIFT);
	private static final int DEL=(-1)*(1 << SCORE_SHIFT);
	private static final int N_SCORE=0;
	private static final int BAD=Integer.MIN_VALUE/2;

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static final boolean GLOBAL=false;


}
