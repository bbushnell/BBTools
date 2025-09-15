package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Tools;

/**
 * PRESENTATION-ONLY VERSION - DO NOT USE
 * This class contains simplified code for publication purposes.
 * For functional implementation, see {@link DriftingAligner}
 * 
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 24, 2025
 */
public class DriftingAlignerConcise implements IDAligner{

	/** Prints a warning */
	public static void main(String[] args) {
		throw new RuntimeException("Nonworking code demo; use WobbleAligner");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public DriftingAlignerConcise() {throw new RuntimeException("Nonworking code demo; use DriftingAligner");}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "Drifting";}
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

	/** This concise function should work, but use DriftingAligner instead.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * @return Identity (0.0-1.0). */
	public static final float alignStaticTrue(byte[] query, byte[] ref, int[] posVector) {
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length, rLen=ref.length; //Cache array lengths
		long[] prev=new long[rLen+1], curr=new long[rLen+1]; // Score arrays for current/previous rows
		for(int j=0; j<=rLen; j++){prev[j]=j;}// Initialize top row with starting position in low bits
		Arrays.fill(curr, BAD); // Initialize current row to BAD
		final int bandWidth0=decideBandwidth(query, ref), maxDrift=2;// Banding parameters
		int center=0, maxPos=0, missingScore=0;// Initialize center-tracking variables outside the loop
		long maxScore=BAD;// Best scoring position
		for(int i=1; i<=qLen; i++){// Fill alignment matrix
			// Bonus bandwidth from low total alignment quality
			final int scoreBonus=32-Integer.numberOfLeadingZeros(missingScore);
			// Bonus bandwidth near the top row
			final int bandWidth=bandWidth0+Math.max(20+bandWidth0*8-maxDrift*i, scoreBonus);
			final int rShift=bandWidth/4;// Shift band to the right to accommodate deletions
			final int drift=Tools.mid(-1, maxPos-center, maxDrift); // Band drift for this round
			center=center+1+drift;// New band center
			final int bandStart=Math.max(1, center-bandWidth+rShift);
			final int bandEnd=Math.min(rLen, center+bandWidth+rShift);
			curr[bandStart-1]=BAD; curr[0]=i*INS;// Clear stale data to the left of the band
			maxScore=BAD; maxPos=0;// Reset max score
			final byte q=query[i-1];// Cache the query
			for(int j=bandStart; j<=bandEnd; j++){// Process only cells within the band
				final byte r=ref[j-1];// Read the reference base
				// Branchless score calculation, with special case for 'N'
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
				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;// Track best score in row
				maxPos=better ? j : maxPos;
			}
			missingScore=i-(int)(maxScore>>SCORE_SHIFT);//How much score is missing vs a perfect match
			long[] temp=prev; prev=curr; curr=temp;// Swap rows
		}
		return postprocess(maxScore, maxPos, qLen, rLen, posVector);// Calculate identity and rStart
	}
	
	static long[] prev=null, curr=null;
	static int qLen=0, rLen=0;
	static long maxValue=0;

	 /** WARNING! This function is nonfunctional demo code.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * @return Identity (0.0-1.0). */
	@SuppressWarnings("unused")
	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector) {
		/* ... Prior initialization remains the same ... */
		final int bandWidth0=decideBandwidth(query, ref), maxDrift=2;// Banding parameters
		int center=0, maxPos=0, missingScore=0;//Initialize center-tracking variables
		long maxScore=BAD;// Best scoring position
		for(int i=1; i<=qLen; i++){// Fill alignment matrix
			// Bonus bandwidth from low total alignment quality
			final int scoreBonus=32-Integer.numberOfLeadingZeros(missingScore);
			// Bonus bandwidth near the top row
			final int bandWidth=bandWidth0+Math.max(20+bandWidth0*8-maxDrift*i, scoreBonus);
			final int rShift=bandWidth/4;// Shift band to the right to accommodate deletions
			final int drift=Tools.mid(-1, maxPos-center, maxDrift);// Band drift this round
			center=center+1+drift;// New band center; cannot drift left
			final int bandStart=Math.max(1, center-bandWidth+rShift);
			final int bandEnd=Math.min(rLen, center+bandWidth+rShift);
			curr[bandStart-1]=BAD; curr[0]=i*INS;// Clear stale data to the left of the band
			maxScore=BAD; maxPos=0;// Reset max score
			final byte q=query[i-1];// Cache the query
			for(int j=bandStart; j<=bandEnd; j++){// Process only cells within the band
				/* ... Prior inner loop remains the same ... */
				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;// Track best score in row
				maxPos=better ? j : maxPos;
			}
			missingScore=i-(int)(maxScore>>SCORE_SHIFT);// Score difference versus a perfect match
			long[] temp=prev; prev=curr; curr=temp;// Swap rows
		}
		return postprocess(maxScore, maxPos, qLen, rLen, posVector);// Calculate identity and rStart
	}

	/** Use alignment information to calculate identity and starting coordinate.
	 * @param maxScore Highest score in last row
	 * @param maxPos Highest-scoring position in last row
	 * @param qLen Query length
	 * @param rLen Reference length
	 * @param posVector Optional array for returning reference start/stop coordinates.
	 * @return Identity */
	private static float postprocess(long maxScore, int maxPos, int qLen, int rLen, int[] posVector) {
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
