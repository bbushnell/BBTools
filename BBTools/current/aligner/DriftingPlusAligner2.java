package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback.
 *Limited to length 2Mbp with 21 position bits.
 *Center of band drifts toward highest score.
 *Band starts wide and narrows to allow glocal alignments.
 *Band dynamically widens in response to low sequence identity.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 24, 2025
 */
public class DriftingPlusAligner2 implements IDAligner{

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

	public DriftingPlusAligner2() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "Drifting+2";}
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
	
	/** Tests for high-identity indel-free alignments needing low bandwidth */
	private static int decideBandwidth(byte[] query, byte[] ref) {
		int subs=0, qLen=query.length, rLen=ref.length;
		int bandwidth=Tools.mid(8, 1+Math.max(qLen, rLen)/16, 40+(int)Math.sqrt(rLen)/4);
		for(int i=0, minlen=Math.min(qLen, rLen); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);}
		return Math.min(subs+1, bandwidth);
	}

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

		final long[] query=Factory.encodeLong(query0, (byte)15);
		final long[] ref=Factory.encodeLong(ref0, (byte)31);
		
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));
		
		// Banding parameters
		final int bandWidth0=decideBandwidth(query0, ref0);
		final int maxDrift=2;
		// Initialize band limits for use outside main loop
		int bandStart=0, bandEnd=rLen-1;
		int center=0;

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}
		
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
			final int quarterBand=bandWidth/4;
			// Center drift for this round
			final int drift=Tools.mid(-maxDrift, maxPos-center, maxDrift);
			// New band center
			center=center+1+drift;
			bandStart=Math.max(1, center-bandWidth+quarterBand);
			bandEnd=Math.min(rLen, center+bandWidth+quarterBand);
			
			//Clear stale data to the left of the band
			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=i*INS;
			
			//Cache the query
			final long q=query[i-1];
			
			//Swap row best scores
//			prevRowScore=maxScore; //Not needed
			maxScore=BAD;
			maxPos=0;
			int posFromSimd=0;
			if(Shared.SIMD) {
				posFromSimd=shared.SIMDAlign.alignBandVectorAndReturnMaxPos(q, ref, bandStart, bandEnd, prev, curr);
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

					// Find max using conditional expressions
					final long maxDiagUp=Math.max(diagScore, upScore);//This is fine

					// Write score to current cell
					curr[j]=maxDiagUp;
				}
			}
			
			//Tail loop for deletions
			long leftCell=curr[bandStart-1];
			if(Shared.SIMD) {
				for(int j=bandStart; j<=bandEnd; j++){//Faster for LWIDTH=4
					final long maxDiagUp=curr[j];
					final long leftScore=leftCell+DEL_INCREMENT;
					leftCell=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
					curr[j]=leftCell;
				}
//				shared.SIMDAlign.processDeletionsTailVectorUnrolled(//Very slow
//				        curr, bandStart, bandEnd, DEL_INCREMENT);
				maxPos=posFromSimd;
				maxScore=curr[posFromSimd];
			}else {
				for(int j=bandStart; j<=bandEnd; j++){
					final long maxDiagUp=curr[j];
					final long leftScore=leftCell+DEL_INCREMENT;
					leftCell=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
					curr[j]=leftCell;

					// Track best score in row
					final boolean better=((leftCell&SCORE_MASK)>maxScore);
					maxScore=better ? leftCell : maxScore;
					maxPos=better ? j : maxPos;
				}
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
	 * Use alignment information to calculate identity and starting coordinate.
	 * @param prev Most recent score row
	 * @param qLen Query length
	 * @param bandStart Beginning of score band for the previous row
	 * @param bandEnd End of score band for the previous row
	 * @param posVector Optional array for returning reference start/stop coordinates.
	 * @return Identity
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

		// The bit field tracks deletion events 
		final int deletions=(int)((bestScore & DEL_MASK) >> POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(bestScore >> SCORE_SHIFT);

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
	public static final boolean GLOBAL=false;

}
