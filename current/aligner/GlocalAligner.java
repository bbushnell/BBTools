package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback. 
 *Limited to length 2Mbp with 21 position bits.
 *
 *@author Brian Bushnell
 *@contributor Isla
 *@date April 23, 2025
 */
public class GlocalAligner implements IDAligner{

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

	/** Default constructor for GlocalAligner instances */
	public GlocalAligner() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "GlocalFull";}
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
	 * Use alignment information to calculate identity and starting coordinate.
	 * @param prev Most recent score row
	 * @param qLen Query length
	 * @param rLen Reference length
	 * @param posVector Optional array for returning reference start/stop coordinates.
	 * @return Identity
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
	
	/** Thread-safe counter tracking total alignment loops performed */
	private static AtomicLong loops=new AtomicLong(0);
	/** Returns the total number of alignment loops performed across all threads */
	public long loops() {return loops.get();}
	/** Sets the loop counter to a specific value.
	 * @param x New loop count value */
	public void setLoops(long x) {loops.set(x);}
	/** Optional output file for alignment visualization */
	public static String output=null;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Bit field definitions
	/** Number of bits used to encode position information (21 bits) */
	private static final int POSITION_BITS=21;
	/** Number of bits used to encode deletion count information (21 bits) */
	private static final int DEL_BITS=21;
	/** Bit shift amount for score encoding (42 bits) */
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;

	// Masks
	/** Bit mask for extracting position information from packed scores */
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	/** Bit mask for extracting deletion count from packed scores */
	private static final long DEL_MASK=((1L << DEL_BITS)-1) << POSITION_BITS;
	/** Bit mask for extracting raw score from packed scores */
	private static final long SCORE_MASK=~(POSITION_MASK | DEL_MASK);

	// Scoring constants
	/** Score increment for matching bases (+1 in upper bits) */
	private static final long MATCH=1L << SCORE_SHIFT;
	/** Score penalty for substitutions (-1 in upper bits) */
	private static final long SUB=(-1L) << SCORE_SHIFT;
	/** Score penalty for insertions (-1 in upper bits) */
	private static final long INS=(-1L) << SCORE_SHIFT;
	/** Score penalty for deletions (-1 in upper bits) */
	private static final long DEL=(-1L) << SCORE_SHIFT;
	/** Score for ambiguous base matches (0) */
	private static final long N_SCORE=0L;
	/** Invalid score marker for out-of-bounds calculations */
	private static final long BAD=Long.MIN_VALUE/2;
	/** Combined deletion penalty and position increment for tracking deletions */
	private static final long DEL_INCREMENT=(1L<<POSITION_BITS)+DEL;

	// Run modes
	/** Debug flag for printing alignment operation details */
	private static final boolean PRINT_OPS=false;
	/** Flag determining global vs local alignment mode */
	public static boolean GLOBAL=false;

}
