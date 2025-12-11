package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Shared;
import shared.Tools;

/**
 * Banded sequence aligner optimized for high-identity alignments.
 * Uses dynamic programming with restricted bands to achieve linear time complexity
 * for sequences with low divergence. Computes approximate nucleotide identity (ANI)
 * without full traceback, using only two arrays for memory efficiency.
 *
 * The banding strategy adapts bandwidth based on observed substitution rate,
 * making it particularly efficient for comparing highly similar sequences
 * where full alignment matrices would be wasteful.
 *
 * @author Brian Bushnell
 * @contributor Isla (Highly-customized Claude instance)
 * @date April 19, 2025
 */
public class BandedPlusAlignerInt implements IDAligner{

	/**
	 * Program entry point that delegates to Test framework.
	 * Uses reflection to determine the calling class and pass it to the test harness.
	 * @param args Command-line arguments for testing
	 * @throws Exception if class reflection or testing fails
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

	public BandedPlusAlignerInt() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the aligner name for identification and logging */
	@Override
	public final String name() {return "BandedInt+";}
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
	 * Adaptively determines optimal bandwidth based on sequence divergence.
	 * Scans the initial portion of both sequences to count substitutions,
	 * then sets bandwidth to accommodate the observed divergence rate.
	 * This optimization reduces computation for highly similar sequences.
	 *
	 * @param query Encoded query sequence
	 * @param ref Encoded reference sequence
	 * @return Bandwidth for the banded alignment (minimum of observed subs + 1 or max bandwidth)
	 */
	private static int decideBandwidth(int[] query, int[] ref) {
		int bandwidth=Math.min(100, 4+Math.max(query.length, ref.length)/8);
		int subs=0;
		for(int i=0, minlen=Math.min(query.length, ref.length); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);
	}

	/**
	 * Performs banded sequence alignment to compute approximate nucleotide identity.
	 * Uses space-efficient dynamic programming with adaptive bandwidth based on
	 * sequence divergence. Optionally swaps sequences to ensure query is shorter.
	 *
	 * @param query0 Query sequence bases
	 * @param ref0 Reference sequence bases
	 * @param posVector Optional int[2] array to receive alignment coordinates [start, end]
	 * @return Nucleotide identity as float between 0.0 and 1.0
	 */
	public static final float alignStatic(byte[] query0, byte[] ref0, int[] posVector) {
		// Swap to ensure query is not inter than ref
		if(posVector==null && query0.length>ref0.length) {
			byte[] temp=query0;
			query0=ref0;
			ref0=temp;
		}

		final int[] query=Factory.encodeInt(query0, (byte)15);
		final int[] ref=Factory.encodeInt(ref0, (byte)31);

		assert(ref.length<=POSITION_MASK) : "Ref is too int: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, 0));
		
		// Banding parameters
		final int bandWidth=decideBandwidth(query, ref);

		// Create two arrays for current and previous rows
		int[] prev=new int[rLen+1], curr=new int[rLen+1];
		Arrays.fill(curr, BAD);

		// Initialize first row with starting position in the lower bits
		for(int j=0; j<=rLen; j++){prev[j]=j;}

		// Initialize band limits for use outside main loop
		int bandStart=0, bandEnd=rLen-1;
		
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
			final int q=query[i-1];
			
			if(Shared.SIMD) {
				shared.SIMDAlign.alignBandVectorInt(q, ref, bandStart, bandEnd, prev, curr, MATCH, N_SCORE, SUB, INS);
			}else {

				// Process only cells within the band
				for(int j=bandStart; j<=bandEnd; j++){
					final int r=ref[j-1];

					// Branchless score calculation
					final boolean isMatch=(q==r);
					final boolean hasN=((q|r)>=15);
					final int scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

					// Read adjacent scores
					final int pj1=prev[j-1], pj=prev[j];
					final int diagScore=pj1+scoreAdd;// Match/Sub
					final int upScore=pj+INS;
					//				final int leftScore=cj1+DEL_INCREMENT;

					// Find max using conditional expressions
					final int maxDiagUp=Math.max(diagScore, upScore);//This is fine

					curr[j]=maxDiagUp;
				}
			}
			
			//Tail loop for deletions
			int leftCell=curr[bandStart-1];
			for(int j=bandStart; j<=bandEnd; j++){
				final int maxDiagUp=curr[j];
				final int leftScore=leftCell+DEL;
				leftCell=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				curr[j]=leftCell;
			}
			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			mloops+=(bandEnd-bandStart+1);

			// Swap rows
			int[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		loops.addAndGet(mloops);
		return postprocess(prev, qLen, bandStart, bandEnd, posVector);
	}
	
	/**
	 * Extracts alignment results from the final DP array row.
	 * Finds the optimal alignment position and computes identity score
	 * using match/mismatch/gap counting. Includes safety checks for
	 * position value corruption.
	 *
	 * @param prev Final row of the DP matrix
	 * @param qLen Query sequence length
	 * @param bandStart Start position of the alignment band
	 * @param bandEnd End position of the alignment band
	 * @param posVector Optional array to receive alignment coordinates
	 * @return Nucleotide identity between 0.0 and 1.0
	 */
	private static final float postprocess(int[] prev, int qLen, int bandStart, int bandEnd, int[] posVector) {
		
		// Find best score outside of main loop
		int maxScore=Integer.MIN_VALUE;
		int maxPos=0;
		for(int j=bandStart; j<=bandEnd; j++){
		    int score=prev[j] & SCORE_MASK;
		    if(score>maxScore){
		        maxScore=score;
		        maxPos=j;
		    }
		}
		
		// Add safeguards to identity calculation
		int bestScore=prev[maxPos];
		int originPos=bestScore & POSITION_MASK;

		// Sanity check on position values
		if(originPos<0 || originPos>maxPos){
		    assert(false) : originPos+", "+maxPos;
		    return 0; // Return safe value if position is corrupted
		}

		int endPos=maxPos;
		int refAlnLength=(endPos-originPos);
		int rawScore=bestScore >> SCORE_SHIFT;

		if(posVector!=null){
		    posVector[0]=originPos;
		    posVector[1]=endPos-1;
		}

		// Calculate net gaps
		int netGaps=Math.abs(qLen-refAlnLength);
		assert(netGaps>=0);
		float matches, insertions, deletions;

		// Apply our improved formulas
		if(qLen>refAlnLength){
		    // More insertions than deletions
		    matches=(rawScore+qLen)/2;
		    insertions=netGaps;
		    deletions=0;
		}else{
		    // More deletions than insertions
		    matches=(rawScore+refAlnLength)/2;
		    insertions=0;
		    deletions=netGaps;
		}

		// Apply safeguards
		matches=Math.max(0, matches);
		float mismatches=Math.max(0, Math.min(qLen, refAlnLength)-matches);
		float operations=matches+mismatches+insertions+deletions;
		assert(operations>=0);

		// Final safety check
		if(operations<=0){return 0f;}

		float identity=(float)matches/operations;

		// Ensure return value is valid
		return Tools.mid(0f, identity, 1f);
	}
	
	/**
	 * Aligns query to a specific window of the reference sequence.
	 * Extracts the reference region and adjusts alignment coordinates
	 * to account for the window offset.
	 *
	 * @param query Query sequence to align
	 * @param ref Complete reference sequence
	 * @param posVector Optional array to receive adjusted alignment coordinates
	 * @param refStart Start position of alignment window in reference
	 * @param refEnd End position of alignment window in reference
	 * @return Nucleotide identity between 0.0 and 1.0
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
	private static final int SUB=(-1) << SCORE_SHIFT;
	private static final int INS=(-1) << SCORE_SHIFT;
	private static final int DEL=(-1) << SCORE_SHIFT;
	private static final int N_SCORE=0;
	private static final int BAD=Integer.MIN_VALUE/2;


}
