package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Shared;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback. 
 *Limited to length 2Mbp with 21 position bits.
 *
 *@author Brian Bushnell
 *@contributor Isla
 *@date May 5, 2025
 */
public class GlocalPlusAligner implements IDAligner{

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

	/** Default constructor for GlocalPlusAligner instance */
	public GlocalPlusAligner() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final String name() {return "Glocal+";}
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
	public static final float alignStatic(byte[] query0, byte[] ref0, int[] posVector) {
		// Swap to ensure query is not longer than ref
		if(posVector==null && query0.length>ref0.length) {
			byte[] temp=query0;
			query0=ref0;
			ref0=temp;
		}
		final int qLen=query0.length;
		final int rLen=ref0.length;
		long mloops=0;

		final long[] query=Factory.encodeLong(query0, (byte)15);
		final long[] ref=Factory.encodeLong(ref0, (byte)31);

		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));
		
		// Banding parameters
		final int bandWidth=rLen;
		// Initialize band limits for use outside main loop
		final int bandStart=1, bandEnd=rLen;

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}
		
		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){

			// Clear first column score
			curr[0]=i*INS;
			
			//Cache the query
			final long q=query[i-1];
			
			if(Shared.SIMD) {
				shared.SIMDAlign.alignBandVector(q, ref, bandStart, bandEnd, prev, curr);
			}else {

				// Process only cells within the band
				for(int j=bandStart; j<=bandEnd; j++){
					final long r=ref[j-1];

					// Branchless score calculation
					final boolean isMatch=(q==r);
					final boolean hasN=((q|r)>=15);
					final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

					// Read adjacent scores
					final long pj1=prev[j-1], pj=prev[j];
					final long diagScore=pj1+scoreAdd;// Match/Sub
					final long upScore=pj+INS;
					//				final long leftScore=cj1+DEL_INCREMENT;

					// Find max using conditional expressions
					final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
					//Changing this conditional to max or removing the mask causes a slowdown.
					//				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;

					curr[j]=maxDiagUp;
				}
			}
			
			//Tail loop for deletions
			long leftCell=curr[bandStart-1];
			for(int j=bandStart; j<=bandEnd; j++){
				final long maxDiagUp=curr[j];
				final long leftScore=leftCell+DEL_INCREMENT;
				leftCell=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				curr[j]=leftCell;
			}
			
			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			mloops+=(bandEnd-bandStart+1);

			// Swap rows
			long[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		loops.addAndGet(mloops);
		return postprocess(prev, qLen, bandStart, bandEnd, posVector);
	}

	/**
	 * Extracts final alignment results from the DP matrix and calculates identity.
	 * Finds optimal score, extracts position information from bit fields,
	 * and solves system of equations to determine match/substitution/indel counts.
	 *
	 * @param prev Final row of the dynamic programming matrix
	 * @param qLen Length of the query sequence
	 * @param bandStart Starting column of the alignment band
	 * @param bandEnd Ending column of the alignment band
	 * @param posVector Optional array to receive alignment positions
	 * @return Calculated identity score between 0.0 and 1.0
	 */
	private static final float postprocess(long[] prev, int qLen, int bandStart, int bandEnd, int[] posVector) {
		// Find best score outside of main loop
		long maxScore=Long.MIN_VALUE;
		int maxPos=bandEnd;
		if(GLOBAL){
			maxScore=prev[bandEnd];
		}else{
			for(int j=bandStart; j<=bandEnd; j++){
				long score=prev[j];
				if(score>maxScore){
					maxScore=score;
					maxPos=j;
				}
			}
		}

		// Extract alignment information
		final int originPos=(int)(maxScore&POSITION_MASK);
		final int endPos=maxPos;
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
		}

		// The bit field tracks deletion events 
		final int deletions=(int)((maxScore & DEL_MASK) >> POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(maxScore >> SCORE_SHIFT);

		// Solve the system of equations:
		// 1. M + S + I = qLen
		// 2. M + S + D = refAlnLength
		// 3. Score = M - S - I - D
		
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
	
	/** Thread-safe counter for tracking total alignment operations performed */
	private static AtomicLong loops=new AtomicLong(0);
	/**
	 * Returns the total number of alignment loops performed across all instances
	 */
	public long loops() {return loops.get();}
	/** Sets the loop counter to a specific value.
	 * @param x New loop count value */
	public void setLoops(long x) {loops.set(x);}
	/** Optional output file path for visualization debugging */
	public static String output=null;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Bit field definitions
	/** Number of bits allocated for position information in packed long values */
	private static final int POSITION_BITS=21;
	/** Number of bits allocated for deletion count in packed long values */
	private static final int DEL_BITS=21;
	/**
	 * Bit shift amount to position score in the upper bits of packed long values
	 */
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;

	// Masks
	/** Bit mask for extracting position information from packed long values */
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	/** Bit mask for extracting deletion count from packed long values */
	private static final long DEL_MASK=((1L << DEL_BITS)-1) << POSITION_BITS;
	/** Bit mask for extracting score information from packed long values */
	private static final long SCORE_MASK=~(POSITION_MASK | DEL_MASK);

	// Scoring constants
	/** Scoring value for matches in the upper bits of packed long values */
	private static final long MATCH=1L << SCORE_SHIFT;
	/** Scoring penalty for substitutions in the upper bits of packed long values */
	private static final long SUB=(-1L) << SCORE_SHIFT;
	/** Scoring penalty for insertions in the upper bits of packed long values */
	private static final long INS=(-1L) << SCORE_SHIFT;
	/** Scoring penalty for deletions in the upper bits of packed long values */
	private static final long DEL=(-1L) << SCORE_SHIFT;
	/** Neutral scoring value for alignments involving ambiguous bases (N) */
	private static final long N_SCORE=0L;
	/** Sentinel value representing invalid or uninitialized alignment scores */
	private static final long BAD=Long.MIN_VALUE/2;
	/**
	 * Combined value for deletion penalty plus position increment used in DP calculations
	 */
	private static final long DEL_INCREMENT=(1L<<POSITION_BITS)+DEL;

	// Run modes
	/** Debug flag for printing alignment operation counts and statistics */
	private static final boolean PRINT_OPS=false;
	/** Flag controlling global vs local alignment mode */
	public static boolean GLOBAL=false;

}
