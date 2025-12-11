package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Shared;

/**
 * Glocal alignment implementation that returns Average Nucleotide Identity (ANI).
 * Uses dynamic programming with only 2 arrays to avoid traceback overhead.
 * Calculates alignment start and end positions without full traceback reconstruction.
 * Limited to sequences of 2MB with 21-bit position encoding.
 *
 * Uses bit-packed score fields to track position, deletions, and score in single long values.
 * Supports both global and local alignment modes via GLOBAL flag.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 5, 2025
 */
public class GlocalPlusAligner4 implements IDAligner{

	/**
	 * Program entry point that delegates to Test class for standardized testing.
	 * Uses reflection to determine the actual class type and passes it to Test framework.
	 * @param args Command-line arguments passed to Test.testAndPrint
	 * @throws Exception If class reflection fails or Test execution throws
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

	public GlocalPlusAligner4() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the display name "Glocal+4" for this aligner implementation */
	@Override
	public final String name() {return "Glocal+4";}
	/**
	 * Aligns two sequences and returns identity score.
	 * @param a First sequence to align
	 * @param b Second sequence to align
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and returns identity score with alignment positions.
	 *
	 * @param a First sequence to align
	 * @param b Second sequence to align
	 * @param pos Output array for alignment start/end positions [start, end]
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold.
	 * Note: minScore parameter is currently ignored in this implementation.
	 *
	 * @param a First sequence to align
	 * @param b Second sequence to align
	 * @param pos Output array for alignment start/end positions [start, end]
	 * @param minScore Minimum score threshold (currently unused)
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns sequences within a specified reference window.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Output array for alignment start/end positions [start, end]
	 * @param rStart Reference window start position
	 * @param rStop Reference window end position
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main alignment algorithm using dynamic programming with bit-packed scoring.
	 * Encodes sequences into long arrays for efficient comparison.
	 * Uses 2-array approach to minimize memory usage.
	 * Optionally swaps query/reference to ensure query is shorter for efficiency.
	 *
	 * Bit packing format in score values:
	 * - Lower 21 bits: alignment start position
	 * - Next 21 bits: deletion count
	 * - Upper bits: alignment score
	 *
	 * @param query0 Query sequence bytes
	 * @param ref0 Reference sequence bytes
	 * @param posVector Optional output array for alignment positions [start, end].
	 * If null, sequences may be swapped for optimization.
	 * @return Identity score between 0.0 and 1.0
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
		final long[] ref=Factory.encodeLong(ref0, (byte)31, 4);

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
	 * Extracts final alignment statistics from completed DP matrix.
	 * Finds optimal score position and unpacks bit fields to calculate identity.
	 * Solves system of equations to determine match/substitution/indel counts.
	 *
	 * System of equations solved:
	 * 1. M + S + I = qLen
	 * 2. M + S + D = refAlnLength
	 * 3. Score = M - S - I - D
	 *
	 * @param prev Final DP row containing scores
	 * @param qLen Query sequence length
	 * @param bandStart Start of scoring band
	 * @param bandEnd End of scoring band
	 * @param posVector Output array for alignment positions [start, end]
	 * @return Identity score calculated as M/(M+S+I+D)
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
	 * Aligns query to a window of the reference sequence.
	 * Extracts reference subsequence and delegates to main alignment method.
	 * Adjusts returned positions to account for window offset.
	 *
	 * @param query Query sequence
	 * @param ref Full reference sequence
	 * @param posVector Output array for alignment positions [start, end]
	 * @param refStart Window start position (inclusive)
	 * @param refEnd Window end position (inclusive)
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
