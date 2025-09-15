package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

/**
 * PRESENTATION-ONLY VERSION - DO NOT USE
 * This class contains simplified code for publication purposes.
 * For optimal implementation, see {@link BandedPlusAligner}
 * 
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback.
 *Limited to length 2Mbp with 21 position bits.
 *Restricts alignment to a fixed band around the diagonal.
 *
 *@author Brian Bushnell
 *@contributor Isla
 *@date April 24, 2025
 */
public class BandedAlignerConcise implements IDAligner{

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

	public BandedAlignerConcise() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "Banded";}
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
		int bandwidth=Math.min(60+(int)Math.sqrt(rLen), 4+Math.max(qLen, rLen)/8);
		for(int i=0, minlen=Math.min(qLen, rLen); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);}
		return Math.min(subs+1, bandwidth);
	}

	/** @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * @return Identity (0.0-1.0). */
	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector) {
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length, rLen=ref.length; //Cache array lengths
		long[] prev=new long[rLen+1], curr=new long[rLen+1]; // Score arrays for current/previous rows
		for(int j=0; j<=rLen; j++){prev[j]=j;}// Initialize top row with starting position in low bits
		Arrays.fill(curr, BAD); // Create current row to BAD
		final int bandWidth=decideBandwidth(query, ref);
		int bandStart=0, bandEnd=rLen-1; // Initialize band limits outside main loop
		for(int i=1; i<=qLen; i++){// Fill alignment matrix
			bandEnd=Math.min(rLen, i+bandWidth);// Calculate band boundaries 
			bandStart=Math.max(1, Math.min(i-bandWidth, rLen-bandWidth));
			curr[bandStart-1]=BAD; curr[0]=i*INS;// Clear stale data to the left of the band
			final byte q=query[i-1];//Cache the query
			// Process only cells within the band
			for(int j=bandStart; j<=bandEnd; j++){// Process only cells within the band
				//...Inner loop unchanged
				
				final byte r=ref[j-1];
				// Branchless score calculation
				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);// +1 match, 0 N, -1 sub
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];// Read adjacent scores
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;// Insertion
				final long leftScore=cj1+DEL_INCREMENT; //Deletion; adjust both score and del counter
				// Find max using conditional or arithmetic expressions
				final long maxDiagUp=Math.max(diagScore, upScore);
				//Changing this conditional to max or removing the mask causes a slowdown.
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				curr[j]=maxValue;// Store the maximum score
			}
			long[] temp=prev; prev=curr; curr=temp;// Swap rows
		}
		return postprocess(prev, qLen, bandStart, bandEnd, posVector);
	}

	/** Use alignment information to calculate identity and starting coordinate.
	 * @param prev Most recent score row
	 * @param qLen Query length
	 * @param bandStart Beginning of score band for the previous row
	 * @param bandEnd End of score band for the previous row
	 * @param posVector Optional array for returning reference start/stop coordinates.
	 * @return Identity  */
	private static final float postprocess(long[] prev, int qLen, int bandStart, int bandEnd, int[] posVector) {
		long maxScore=prev[bandEnd];
		int maxPos=bandEnd;
		for(int j=bandStart; j<=bandEnd; j++){// Find best score outside of main loop
			long score=prev[j];
			if(score>maxScore){maxScore=score; maxPos=j;}
		}
		// Extract alignment information
		final int originPos=(int)(maxScore&POSITION_MASK);
		final int deletions=(int)((maxScore&DEL_MASK) >> POSITION_BITS);
		final int refAlnLength=(maxPos-originPos);
		final int rawScore=(int)(maxScore >> SCORE_SHIFT);
		// Solve the system of equations to calculate operation counts:
		// 1. M + S + I = qLen
		// 2. M + S + D = refAlnLength
		// 3. M - S - I - D = Score
		final int insertions=Math.max(0, qLen+deletions-refAlnLength);
		final float matches=((rawScore+qLen+deletions)/2f);
		final float substitutions=Math.max(0, qLen-matches-insertions);
		//Generate results
	    final float identity=matches/(matches+substitutions+insertions+deletions);
		if(posVector!=null){posVector[0]=originPos; posVector[1]=maxPos-1;} //Store rstart, rstop
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

}
