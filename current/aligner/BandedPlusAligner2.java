package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.PreParser;
import shared.Shared;

/**
 * Advanced traceback-free banded alignment with exact ANI calculation using mathematical constraint solving.
 * Restricts alignment to adaptive diagonal bands with 64-bit bit-packed encoding supporting sequences up to 2Mbp.
 * Uses SIMD vectorization when available and falls back to optimized scalar processing.
 *
 * @author Brian Bushnell
 * @contributor Isla (Highly-customized Claude instance)
 * @date April 24, 2025
 */
public class BandedPlusAligner2 implements IDAligner{

	/**
	 * Program entry point that delegates to Test framework for algorithm validation.
	 * Uses reflection to determine the calling class and passes it to Test.testAndPrint().
	 * @param args Command-line arguments for testing configuration
	 * @throws Exception If reflection or testing fails
	 */
	public static <C extends IDAligner> void main(String[] args) throws Exception {
		args=new PreParser(args, System.err, null, false, true, false).args;
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public BandedPlusAligner2() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the algorithm name for identification purposes */
	@Override
	public final String name() {return "Banded+2";}
	/**
	 * Aligns two sequences and returns identity without position information.
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and returns identity with optional position tracking.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional int[2] array for returning {rStart, rStop} coordinates
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold (currently ignored).
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional int[2] array for returning {rStart, rStop} coordinates
	 * @param minScore Minimum score threshold (not currently used)
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns query to a specific window of the reference sequence.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional int[2] array for returning {rStart, rStop} coordinates
	 * @param rStart Reference window start position
	 * @param rStop Reference window stop position
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculates optimal band width by analyzing early mismatch patterns.
	 * Uses square-root scaling with sequence length and early termination when mismatches exceed threshold.
	 *
	 * @param query Query sequence for bandwidth analysis
	 * @param ref Reference sequence for bandwidth analysis
	 * @return Optimal band width for alignment efficiency
	 */
	private static int decideBandwidth(byte[] query, byte[] ref) {
		int subs=0, qLen=query.length, rLen=ref.length;
		int bandwidth=Math.min(60+(int)Math.sqrt(rLen), 4+Math.max(qLen, rLen)/8);
		for(int i=0, minlen=Math.min(qLen, rLen); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);}
		return Math.min(subs+1, bandwidth);
	}

	/**
	 * Core alignment method implementing traceback-free banded alignment with mathematical constraint solving.
	 * Uses 64-bit bit-packed encoding and dual rolling arrays achieving O(n) space complexity.
	 * Automatically swaps sequences when beneficial and uses SIMD acceleration when available.
	 *
	 * @param query0 Query sequence
	 * @param ref0 Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} coordinates
	 * @return Identity score calculated using exact operation counts
	 */
	public static final float alignStatic(byte[] query0, byte[] ref0, int[] posVector) {
		// Swap to ensure query is not longer than ref
		if(posVector==null && query0.length>ref0.length) {
			byte[] temp=query0;
			query0=ref0;
			ref0=temp;
		}

		final long[] query=Factory.encodeLong(query0, (byte)15);
		final long[] ref=Factory.encodeLong(ref0, (byte)31);
		
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));
		
		// Banding parameters
		final int bandWidth=decideBandwidth(query0, ref0);
		// Initialize band limits for use outside main loop
		int bandStart=1, bandEnd=rLen-1;

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}
		
		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
			// Calculate band boundaries 
			bandStart=Math.max(1, Math.min(i-bandWidth, rLen-bandWidth));
			bandEnd=Math.min(rLen, i+bandWidth);
			
			//Clear stale data to the left of the band
			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=i*INS;
			
			//Cache the query
			final long q=query[i-1];
			
			if(Shared.SIMD) {
				shared.SIMDAlign.alignBandVectorDel(q, ref, bandStart, bandEnd, prev, curr);
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
//				Test.print(curr, "Row "+i);
			}
			
			if(!Shared.SIMD) {//Tail loop for deletions
				long leftCell=curr[bandStart-1];
				for(int j=bandStart; j<=bandEnd; j++){
					final long maxDiagUp=curr[j];
					final long leftScore=leftCell+DEL_INCREMENT;
					leftCell=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
					curr[j]=leftCell;
				}
			}
//			Test.print(curr, "Row "+i+"D");
			
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
	 * Extracts alignment results using mathematical constraint solving for exact operation counts.
	 * Uses the constraint system: M+S+I=qLen, M+S+D=refAlnLength, Score=M-S-I-D.
	 * Calculates identity without requiring traceback matrices.
	 *
	 * @param prev Final score row containing alignment results
	 * @param qLen Query sequence length
	 * @param bandStart Beginning of score band for result extraction
	 * @param bandEnd End of score band for result extraction
	 * @param posVector Optional array for returning reference coordinates
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
	 * Aligns query to a specific window of the reference sequence with coordinate adjustment.
	 * Extracts reference region and adjusts returned coordinates to global reference positions.
	 *
	 * @param query Query sequence
	 * @param ref Full reference sequence
	 * @param posVector Optional int[2] for returning global {rStart, rStop} coordinates
	 * @param refStart Window start position in reference
	 * @param refEnd Window end position in reference
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
