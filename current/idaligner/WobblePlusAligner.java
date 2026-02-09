package idaligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Tools;
import structures.RingBuffer;

/**
 * High-performance sequence aligner optimized for calculating average nucleotide identity (ANI).
 * Uses dynamic banding algorithm with position tracking to avoid traceback computation.
 * Limited to sequences up to 2Mbp with 21-bit position encoding.
 * Band width adapts to sequence identity and alignment quality for optimal performance.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 7, 2025
 */
public class WobblePlusAligner implements IDAligner{

	/**
	 * Program entry point that delegates to Test class for standardized testing.
	 * Uses reflection to determine the calling class for polymorphic testing.
	 * @param args Command-line arguments passed to Test.testAndPrint()
	 * @throws Exception If class reflection or testing fails
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

	public WobblePlusAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the aligner name for identification */
	@Override
	public final String name() {return "Wobble+";}
	/**
	 * Aligns two sequences and returns identity without position information.
	 * @param a First sequence
	 * @param b Second sequence
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and optionally returns alignment coordinates.
	 *
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional array to receive {rStart, rStop} coordinates
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold.
	 * Currently ignores minScore parameter and performs full alignment.
	 *
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Optional array to receive {rStart, rStop} coordinates
	 * @param minScore Minimum score threshold (unused)
	 * @return Identity score between 0.0 and 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns query to a specific window of the reference sequence.
	 *
	 * @param a Query sequence
	 * @param b Reference sequence
	 * @param pos Optional array to receive {rStart, rStop} coordinates
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
	 * Calculates optimal band width for alignment based on sequence similarity.
	 * Tests initial alignment quality to determine if narrow banding is appropriate.
	 * Returns smaller bandwidth for high-identity sequences with few indels.
	 *
	 * @param query Query sequence for testing
	 * @param ref Reference sequence for testing
	 * @return Optimal band width between 7 and 24 bases
	 */
	private static int decideBandwidth(byte[] query, byte[] ref) {
		int bandwidth=Tools.mid(7, 1+Math.max(query.length, ref.length)/24, 24);
		int subs=0;
		for(int i=0, minlen=Math.min(query.length, ref.length); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);
	}

	/**
	 * Core alignment method implementing dynamic banded Smith-Waterman algorithm.
	 * Encodes sequences, initializes scoring matrices, and performs alignment with adaptive
	 * band width that responds to sequence identity. Tracks position information without
	 * requiring traceback by encoding coordinates in score values.
	 *
	 * @param query0 Query sequence to align
	 * @param ref0 Reference sequence to align against
	 * @param posVector Optional int[2+] array for returning alignment coordinates.
	 * [0]=rStart, [1]=rStop, [2]=rawScore (optional), [3]=deletions (optional)
	 * @return Identity score between 0.0 and 1.0
	 */
	public static final float alignStatic(byte[] query0, byte[] ref0, int[] posVector) {
		// Swap to ensure query is not longer than ref
		if(posVector==null && query0.length>ref0.length) {
			byte[] temp=query0;
			query0=ref0;
			ref0=temp;
		}

		final byte[] query=Factory.encodeByte(query0, (byte)15);
		final byte[] ref=Factory.encodeByte(ref0, (byte)31);
		
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));
		
		// Banding parameters
		final int bandWidth0=decideBandwidth(query, ref);
		final int maxDrift=2, ringSize=(bandWidth0*5)/4;
		final RingBuffer ring=new RingBuffer(ringSize);

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}

		// Initialize band limits for use outside main loop
		int bandStart=1, bandEnd=rLen-1;
		int center=0;
		
		// Best scoring position
		int maxPos=0;
		long maxScore=2*SUB;
		
		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
			// Calculate bonus bandwidth due to low local alignment quality
			final int oldMaxScore=(int)(ring.getOldestUnchecked()>>SCORE_SHIFT);
			final int recentMissingScore=(oldMaxScore+ringSize)-(int)(maxScore>>SCORE_SHIFT);
			final int scoreBonus=Math.max(0, Math.min(ringSize*2, recentMissingScore*2));
			
			// Bonus bandwidth near the top row
			final int bandWidth=bandWidth0+Math.max(10+bandWidth0*8-maxDrift*i, scoreBonus);
			final int quarterBand=bandWidth/4;
			// Center drift for this round
			final int drift=Tools.mid(-1, maxPos-center, maxDrift);
			// New band center
			center=center+1+drift;
			bandStart=Math.max(bandStart, center-bandWidth+quarterBand);
			bandEnd=Math.min(rLen, center+bandWidth+quarterBand);
			
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
				final boolean isMatch=(q==r);
				final boolean hasN=((q|r)>=15);
				final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL_INCREMENT;

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				// Changing this conditional to max or removing the mask causes a slowdown.
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				
				// Write score to current cell
				curr[j]=maxValue;
				
				// Track best score in row
				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			mloops+=(bandEnd-bandStart+1);
			
			// Swap rows
			long[] temp=prev;
			prev=curr;
			curr=temp;
			ring.add(maxScore);
		}
		if(viz!=null) {viz.shutdown();}// Terminate visualizer
		if(GLOBAL) {maxPos=rLen;maxScore=prev[rLen-1]+DEL_INCREMENT;}//The last cell may be empty 
		loops.addAndGet(mloops);
		return postprocess(maxScore, maxPos, qLen, rLen, posVector);
	}
	
	/**
	 * Extracts alignment information from encoded score and calculates final identity.
	 * Decodes position and deletion count from bit fields in maxScore value.
	 * Solves system of equations to determine match, substitution, and indel counts.
	 *
	 * @param maxScore Encoded score containing position, deletions, and raw score
	 * @param maxPos Highest-scoring position in reference
	 * @param qLen Query sequence length
	 * @param rLen Reference sequence length
	 * @param posVector Optional array to receive alignment coordinates and statistics
	 * @return Final identity score calculated from operation counts
	 */
	private static float postprocess(long maxScore, int maxPos, int qLen, int rLen, int[] posVector) {
		// For conversion to global alignments
		if(GLOBAL && maxPos<rLen) {
			int dif=rLen-maxPos;
			maxPos+=dif;
			maxScore+=(dif*DEL_INCREMENT);
		}
		
		// Extract alignment information
		final int originPos=(int)(maxScore&POSITION_MASK);
		final int endPos=maxPos;

		// Calculate alignment statistics
		final int deletions=(int)((maxScore & DEL_MASK) >> POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(maxScore >> SCORE_SHIFT);
		
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
			if(posVector.length>2) {posVector[2]=rawScore;}
			if(posVector.length>3) {posVector[3]=deletions;}
		}
		
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
	 * Aligns query to a specific window of the reference sequence.
	 * Extracts reference region and adjusts returned coordinates to full reference space.
	 *
	 * @param query Query sequence to align
	 * @param ref Full reference sequence
	 * @param posVector Optional array to receive alignment coordinates in full reference space
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
		assert(posVector[1]>0) : id+", "+Arrays.toString(posVector)+", "+refStart;
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
	private static final long DEL_INCREMENT=DEL+(1L<<POSITION_BITS);

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static final boolean GLOBAL=false;

}
