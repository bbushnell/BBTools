package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Tools;

/**
 * Sequence aligner for calculating average nucleotide identity (ANI) using dynamic programming.
 * Uses a drifting band approach to restrict computation to a diagonal region that moves
 * toward areas of highest scoring alignment. Optimized for high-identity sequences by
 * starting with a wide band that narrows as alignment progresses, widening again when
 * sequence identity drops. Avoids traceback by encoding position and match information
 * directly in the scoring matrix. Limited to sequences under 2Mbp due to 21-bit position encoding.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date April 24, 2025
 */
public class DriftingAlignerM implements IDAligner{

	/**
	 * Entry point that delegates to Test class for standardized aligner testing.
	 * Uses reflection to determine the actual aligner class from the call stack.
	 * @param args Command line arguments passed to Test.testAndPrint()
	 * @throws Exception from Test.testAndPrint() or class reflection
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

	public DriftingAlignerM() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the name identifier for this aligner implementation.
	 * @return "DrifitingM" (note: contains typo from original source) */
	@Override
	public final String name() {return "DrifitingM";}
	/**
	 * Aligns two sequences and returns identity score.
	 * @param a First sequence
	 * @param b Second sequence
	 * @return Identity fraction (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and returns identity score with position information.
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional int[2] array to receive alignment start/end positions
	 * @return Identity fraction (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold.
	 * Note: minScore parameter is currently ignored in this implementation.
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional int[2] array to receive alignment start/end positions
	 * @param minScore Minimum score threshold (ignored)
	 * @return Identity fraction (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns sequences within a specific reference window.
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional int[2] array to receive alignment start/end positions
	 * @param rStart Reference window start position
	 * @param rStop Reference window end position
	 * @return Identity fraction (0.0-1.0)
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Determines initial alignment bandwidth based on sequence lengths and substitution rate.
	 * Tests alignment quality by counting mismatches in the prefix and adjusts bandwidth
	 * accordingly. Uses adaptive sizing between 8-40 bp based on sequence length.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @return Initial bandwidth for alignment banding
	 */
	private static int decideBandwidth(byte[] query, byte[] ref) {
		int bandwidth=Tools.mid(8, 1+Math.max(query.length, ref.length)/16, 40);
		int subs=0;
		for(int i=0, minlen=Math.min(query.length, ref.length); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);
	}

	/**
	 * Core static alignment method implementing drifting band dynamic programming.
	 * Uses a diagonal band that starts wide, narrows as alignment progresses, and widens
	 * when sequence identity drops. Encodes position and match count in score values to
	 * avoid traceback. Band center drifts toward highest-scoring regions.
	 * @param query Query sequence (may be swapped with ref if longer and posVector is null)
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of optimal alignment
	 * @return Identity fraction (0.0-1.0) based on matches/(matches+subs+indels)
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
		
		//Create a visualizer if an output file is defined
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, MATCH_BITS));
		
		// Banding parameters
		final int bandWidth0=decideBandwidth(query, ref);
		final int maxDrift=2;

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}

		// Initialize band limits for use outside main loop
		int bandStart=0, bandEnd=rLen-1;
		int center=0;
		
		// Best scoring position
		int maxPos=0;
		long maxScore=BAD;
		int missingScore=0;
//		long prevRowScore=BAD;

		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
			// Calculate band boundaries
			// Bonus bandwidth due to low alignment quality
			final int scoreBonus=32-Integer.numberOfLeadingZeros(missingScore);
			// Bonus bandwidth near the top row
			final int bandWidth=bandWidth0+Math.max(20+bandWidth0*8-maxDrift*i, scoreBonus);
			// Center drift for this round
			final int drift=Tools.mid(-maxDrift, maxPos-center, maxDrift);
			// New band center
			center=center+1+drift;
			bandStart=Math.max(1, center-bandWidth);
			bandEnd=Math.min(rLen, center+bandWidth);
			
			//Clear stale data to the left of the band
			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=i*INS;
			
			//Cache the query
			final byte q=query[i-1];
			
			//Swap row best scores
//			prevRowScore=maxScore; //Not needed
			maxScore=BAD;
			maxPos=0;
			
			// Process only cells within the band
			for(int j=bandStart; j<=bandEnd; j++){
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final long scoreAdd=isMatch ? MATCH_INCREMENT : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL;

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				//This mask and conditional is no longer needed in the match version.
				//final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				final long maxValue=Math.max(maxDiagUp, leftScore);
				
				// Write score to current cell
				curr[j]=maxValue;
				
				// Track best score in row
				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			mloops+=(bandEnd-bandStart+1);
			final int score=(int)(maxScore>>SCORE_SHIFT);
			missingScore=i-score;//How much score is missing compared to a perfect match

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
	 * Extracts alignment statistics from final DP row and calculates identity.
	 * Decodes position and match information from bit-packed scores, solves system
	 * of equations to determine indel counts, and computes final identity percentage.
	 * @param prev Final row of alignment matrix with bit-packed scores
	 * @param qLen Query sequence length
	 * @param bandStart Start position of scoring band in final row
	 * @param bandEnd End position of scoring band in final row
	 * @param posVector Optional array for returning reference start/stop coordinates
	 * @return Identity fraction calculated as matches/(matches+substitutions+indels)
	 */
	private static final float postprocess(long[] prev, int qLen, int bandStart, int bandEnd, int[] posVector) {
		// Find best score outside of main loop
		long maxScore=Long.MIN_VALUE;
		int maxPos=0;
		for(int j=bandStart; j<=bandEnd; j++){
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

		// Extract tracked matches
		final int matches=(int)((bestScore & MATCH_MASK) >> POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(bestScore >> SCORE_SHIFT);

		// Solve the system of equations:
		// 1. M + S + I = qLen
		// 2. M + S + D = refAlnLength
		// 3. Score = M - S - I - D

		// From equations 1 and 2:
		// I - D = qLen - refAlnLength
		final int iMinusD=qLen-refAlnLength;

		// From equation 2:
		// S + D = refAlnLength - M
		final int sPlusD=refAlnLength-matches;
		
		final int deletions=Math.max(0, (2*matches-rawScore-qLen));

		// Now we can calculate the rest:
		final int substitutions=sPlusD-deletions;
		final int insertions=iMinusD+deletions;
		final float identity=matches/(float)(matches+substitutions+insertions+deletions);

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
	 * Aligns query to a specific window of the reference sequence.
	 * Extracts reference region and delegates to main alignment method, then adjusts
	 * returned coordinates to account for window offset.
	 * @param query Query sequence
	 * @param ref Full reference sequence
	 * @param posVector Optional int[2] for returning adjusted {rStart, rStop} coordinates
	 * @param refStart Window start position in reference (inclusive)
	 * @param refEnd Window end position in reference (inclusive)
	 * @return Identity fraction (0.0-1.0)
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
	private static final int MATCH_BITS=21;
	private static final int SCORE_SHIFT=POSITION_BITS+MATCH_BITS;

	// Masks
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	private static final long MATCH_MASK=((1L << MATCH_BITS)-1) << POSITION_BITS;
	private static final long SCORE_MASK=~(POSITION_MASK | MATCH_MASK);

	// Scoring constants
	private static final long MATCH=1L << SCORE_SHIFT;
	private static final long SUB=(-1L) << SCORE_SHIFT;
	private static final long INS=(-1L) << SCORE_SHIFT;
	private static final long DEL=(-1L) << SCORE_SHIFT;
	private static final long N_SCORE=0L;
	private static final long BAD=Long.MIN_VALUE/2;
	private static final long MATCH_INCREMENT=MATCH+(1L<<POSITION_BITS);

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static final boolean GLOBAL=false;

}
