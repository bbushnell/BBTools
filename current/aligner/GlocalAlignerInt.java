package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Identity is approximate.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 19, 2024
 */
public class GlocalAlignerInt implements IDAligner{

	/** Main() passes the args and class to Test to avoid redundant code */
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

	@Override
	public final String name() {return "GlocalInt";}
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * If the posVector is null, sequences may be swapped so that the query is shorter.
	 * @return Identity (0.0-1.0).
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
	 * Lightweight wrapper for aligning to a window of the reference.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * If the posVector is null, sequences may be swapped so that the query is shorter.
	 * @param rStart Alignment window start.
	 * @param to Alignment window stop.
	 * @return Identity (0.0-1.0).
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
