package idaligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Performs sequence alignment to calculate Average Nucleotide Identity (ANI).
 * Uses a space-efficient algorithm with only 2 arrays and no traceback.
 * Calculates exact alignment scores and reference coordinates without full
 * traceback matrix reconstruction. Limited to sequences up to 2Mbp with
 * 21 position bits.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date April 23, 2025
 */
public class GlocalAligner implements IDAligner{

	/**
	 * Program entry point that delegates testing to the Test class.
	 * Uses reflection to determine the calling class for testing purposes.
	 *
	 * @param <C> IDAligner implementation type
	 * @param args Command-line arguments for testing
	 * @throws Exception If reflection or testing fails
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

	public GlocalAligner() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the aligner name "GlocalFull" */
	@Override
	public final String name() {return "GlocalFull";}
	/**
	 * Aligns two sequences and returns identity score.
	 * @param a First sequence
	 * @param b Second sequence
	 * @return Identity score (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and returns identity score with position information.
	 *
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional array for returning alignment coordinates
	 * @return Identity score (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold.
	 * Note: minScore parameter is ignored in current implementation.
	 *
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional array for returning alignment coordinates
	 * @param minScore Minimum score threshold (currently unused)
	 * @return Identity score (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns sequences within a specified reference region.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional array for returning alignment coordinates
	 * @param rStart Reference region start position
	 * @param rStop Reference region end position
	 * @return Identity score (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Core alignment algorithm that calculates ANI between two sequences.
	 * Uses dynamic programming with space-efficient row-by-row processing.
	 * Automatically swaps sequences to ensure query is not longer than reference
	 * when position vector is not needed. Implements bit-packed scoring to track
	 * position and deletion information without full traceback.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} coordinates.
	 * If null, sequences may be swapped for efficiency
	 * @return Identity score (0.0-1.0)
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
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));
		
//		// Banding parameters
//		final int bandWidth=rLen;
//		// Initialize band limits for use outside main loop
//		final int bandStart=1, bandEnd=rLen-1;

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}
		
		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
//			// Calculate band boundaries 
//			bandStart=Math.max(1, Math.min(i-bandWidth, rLen-bandWidth));
//			bandEnd=Math.min(rLen, i+bandWidth);
//			
//			//Clear stale data to the left of the band
//			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=i*INS;
			
			//Cache the query
			final byte q=query[i-1];
//			curr[0]=i*MATCH/1024;//This is just for pretty pictures

			for(int j=1; j<=rLen; j++){
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL_INCREMENT;

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				//Changing this conditional to max or removing the mask causes a slowdown.
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				
				curr[j]=maxValue;
			}
			if(viz!=null) {viz.print(curr, 1, rLen, rLen);}
			mloops+=rLen;

			// Swap rows
			long[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		loops.addAndGet(mloops);
		return postprocess(prev, qLen, rLen, posVector);
	}
	
	/**
	 * Extracts alignment statistics from the final scoring row and calculates identity.
	 * Uses bit field extraction to determine position, deletions, and raw score.
	 * Solves system of equations to calculate matches, substitutions, and insertions
	 * without full traceback reconstruction.
	 *
	 * @param prev Final row of alignment scores containing bit-packed information
	 * @param qLen Query sequence length
	 * @param rLen Reference sequence length
	 * @param posVector Optional array for returning reference start/stop coordinates
	 * @return Identity score calculated as matches/(matches+substitutions+insertions+deletions)
	 */
	private static final float postprocess(long[] prev, int qLen, int rLen, int[] posVector) {
		long maxScore=prev[rLen];// Find best score outside of main loop
		int maxPos=rLen;
		for(int j=1; !GLOBAL && j<=rLen; j++){
			long score=prev[j];
			if(score>maxScore){
				maxScore=score;
				maxPos=j;
			}
		}

		// Extract alignment information
		final long bestScore=prev[maxPos];
		final int originPos=(int)(bestScore&POSITION_MASK);
		final int endPos=maxPos;
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
		}

		// The bit field tracks deletion events 
		final int deletions=(int)((bestScore & DEL_MASK) >> POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(bestScore >> SCORE_SHIFT);

		// Solve the system of equations:
		// 1. M + S + I = qLen
		// 2. M + S + D = refAlnLength
		// 3. M - S - I - D = Score
		
		// Calculate operation counts
		final int insertions=Math.max(0, qLen+deletions-refAlnLength);
		final float matches=((rawScore+qLen+deletions)/2f);
		final float substitutions=Math.max(0, qLen-matches-insertions);
	    final float identity=matches/(matches+substitutions+insertions+deletions);

		if(PRINT_OPS) {
			System.err.println("originPos="+originPos);
			System.err.println("endPos="+endPos);
			System.err.println("qLen="+qLen);
			System.err.println("matches="+matches);
			System.err.println("refAlnLength="+refAlnLength);
			System.err.println("rawScore="+rawScore);
			System.err.println("deletions="+deletions);
			System.err.println("matches="+matches);
			System.err.println("substitutions="+substitutions);
			System.err.println("insertions="+insertions);
			System.err.println("identity="+identity);
		}
		
		return identity;
	}

	/**
	 * Aligns query to a specified window of the reference sequence.
	 * Extracts the reference region and adjusts returned coordinates to
	 * account for the window offset.
	 *
	 * @param query Query sequence
	 * @param ref Full reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} coordinates
	 * @param refStart Window start position (inclusive)
	 * @param refEnd Window end position (inclusive)
	 * @return Identity score (0.0-1.0) for alignment within the window
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
	private static final boolean PRINT_OPS=false;
	public static boolean GLOBAL=false;

}
