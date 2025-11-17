package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import structures.IntList;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback.
 *Limited to length 2Mbp with 21 position bits.
 *
 *@author Brian Bushnell
 *@contributor Isla
 *@date April 24, 2024
 */
public class IntListAligner implements IDAligner{
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a new IntListAligner instance */
	public IntListAligner() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "IntList";}
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
	    final int refLen=ref.length;
	    long mloops=0;

	    // Banding parameters
	    final int bandWidth=Math.min(200, 1+Math.max(qLen, refLen)/4);
	    final int offset=qLen-refLen;

	    // Create arrays for current and previous rows
	    long[] prev=new long[refLen+1], curr=new long[refLen+1];
	    
	    // Create IntLists for tracking active positions
	    IntList activePositions = new IntList(bandWidth*2+3);
	    IntList nextActivePositions = new IntList(bandWidth*2+3);

	    // Initialize first row with starting position in the lower bits
	    for(int j=0; j<=refLen; j++){prev[j]=j;}
	    
	    // Initialize active positions for the first band
	    for(int j=0; j<=Math.min(refLen, bandWidth+1); j++){activePositions.add(j);}

	    // Best score tracking
	    long maxScore = Long.MIN_VALUE;
	    int maxPos = 0; // Position within the region
	    
	    // Fill alignment matrix
	    for(int i=1; i<=qLen; i++){
	        // First column
	        curr[0]=i*INS;
	        
	        // Clear next positions list
	        nextActivePositions.clear();
	        
	        // Always add first column to next row
	        nextActivePositions.add(0);

	        // Process only positions in the active list
	        for(int idx=0; idx<activePositions.size; idx++){
	            int j = activePositions.array[idx];
	            if(j == 0) continue; // Skip first column, already processed
	            
	            final byte q=query[i-1];
	            final byte r=ref[j-1];

	            // Branchless score calculation
	            final boolean isMatch=(q==r);
	            final boolean hasN=(q=='N' || r=='N');
	            final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : MISMATCH);

	            // Calculate scores
	            final long pj1=prev[j-1], pj=prev[j];
	            final long cj1=(j>1 ? curr[j-1] : 0);
	            final long diagScore=pj1+scoreAdd;// Match/Sub
	            final long upScore=pj+INS;
	            final long leftScore=cj1+DEL_INCREMENT2;

	            // Find max using conditional expressions
	            final long maxDiagUp=diagScore>=upScore ? diagScore : upScore;
	            final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;

	            curr[j]=maxValue;
	            
	            // Calculate band constraint once
	            final int diagCenter = i-offset;
	            final boolean jInBand = Math.abs(j-diagCenter) <= bandWidth;
	            final int lastPos = (nextActivePositions.size==0) ? -1 : nextActivePositions.array[nextActivePositions.size-1];

	            // ONLY add j if it's within the band (this was missing before)
	            if(jInBand && lastPos < j) {
	                nextActivePositions.add(j);
	            }

	            // j-1 (only if within band and > last position)
	            if(j-1 > 0 && Math.abs((j-1)-diagCenter) <= bandWidth && 
	               (nextActivePositions.size==0 || nextActivePositions.array[nextActivePositions.size-1] < j-1)) {
	                nextActivePositions.add(j-1);
	            }

	            // j+1 (only if within band and > last position)
	            if(j+1 <= refLen && Math.abs((j+1)-diagCenter) <= bandWidth && 
	               (nextActivePositions.size==0 || nextActivePositions.array[nextActivePositions.size-1] < j+1)) {
	                nextActivePositions.add(j+1);
	            }
	        	mloops++;
				
	            // Track best score in last row
	            final boolean better=(i == qLen && (maxValue&SCORE_MASK)>maxScore);
	            maxScore=better ? maxValue : maxScore;
	            maxPos=better ? j : maxPos;
	        }

	        // Swap rows
	        long[] temp=prev;
	        prev=curr;
	        curr=temp;
	        
	        // Swap position lists
	        IntList tempList = activePositions;
	        activePositions = nextActivePositions;
	        nextActivePositions = tempList;
	    }

	    // Extract alignment information
	    final int originPos=(int)(maxScore&POSITION_MASK);
	    final int endPos=maxPos;
	    if(posVector!=null){
	        posVector[0]=originPos;
	        posVector[1]=endPos-1;
	    }

	    // Calculate alignment statistics
	    final int deletions=(int)((maxScore & DEL_MASK) >> POSITION_BITS);
	    final int refAlnLength=(endPos-originPos);
	    final long rawScore=maxScore >> SCORE_SHIFT;

	    // Calculate operation counts
	    final int insertions=Math.max(0, qLen+deletions-refAlnLength);
	    final float matches=((rawScore+qLen+deletions)/2f);
	    final float mismatches=Math.max(0, qLen-matches-insertions);

	    // Calculate identity
	    return matches/(matches+mismatches+insertions+deletions);
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
	
	/** Thread-safe counter for tracking total alignment loops performed */
	private static AtomicLong loops=new AtomicLong(0);
	/** Gets the total number of alignment loops performed */
	public long loops() {return loops.get();}
	/** Sets the alignment loop counter.
	 * @param x New loop count value */
	public void setLoops(long x) {loops.set(x);}
	/** Optional output string for debugging or logging purposes */
	public static String output=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Bit field definitions
	/** Number of bits used to encode position information in packed long values */
	private static final int POSITION_BITS=21;
	/** Number of bits used to encode deletion count in packed long values */
	private static final int DEL_BITS=21;
	/** Bit shift value for extracting score from packed long values */
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;
	
	// Masks
	/** Bit mask for extracting position information from packed long values */
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	/** Bit mask for extracting deletion count from packed long values */
	private static final long DEL_MASK=((1L << DEL_BITS)-1) << POSITION_BITS;
	/** Bit mask for extracting alignment score from packed long values */
	private static final long SCORE_MASK=~(POSITION_MASK | DEL_MASK);

	// Scoring constants
	/** Score value for matching bases in alignment */
	private static final long MATCH=1L << SCORE_SHIFT;
	/** Score penalty for mismatched bases in alignment */
	private static final long MISMATCH=(-1L) << SCORE_SHIFT;
	/** Score penalty for insertion operations in alignment */
	private static final long INS=(-1L) << SCORE_SHIFT;
	/** Score penalty for deletion operations in alignment */
	private static final long DEL=(-1L) << SCORE_SHIFT;
	/** Score value for alignments involving ambiguous N bases */
	private static final long N_SCORE=0L;
	/** Score value representing invalid or bad alignment states */
	private static final long BAD=Long.MIN_VALUE/2;
	/** Increment value for tracking deletion counts in packed scoring */
	private static final long DEL_INCREMENT=1L<<POSITION_BITS;
	/** Combined deletion increment and penalty for efficient calculation */
	private static final long DEL_INCREMENT2=DEL_INCREMENT+DEL;

}
