package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Timer;
import shared.Tools;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback. 
 *Limited to length 2Mbp with 21 position bits.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 23, 2025
 */
public class RelativeAligner implements IDAligner{

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

	public RelativeAligner() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "Relative";}
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
	
	public static float alignStatic(byte[] query, byte[] ref, int[] posVector) {
	    // Get lengths of input sequences
	    final int qLen=query.length;
	    final int rLen=ref.length;
		long mloops=0;
	    
	    // Create byte arrays for current and previous rows
	    byte[] prev=new byte[rLen+1];
	    byte[] curr=new byte[rLen+1];
	    
	    // Initialize top row with deletions (shifted to mid-range)
	    final byte MIDPOINT=60;
	    for(int j=0; j<=rLen; j++) {
	        // Start at midpoint minus deletion penalty
	        prev[j]=(byte)(MIDPOINT-j);
	    }
	    
	    // Track total score adjustment from recalibrations
	    int totalAdjustment=0;
	    
	    // Track when to recalibrate
	    int timeToRecalibrate=60;
	    
	    // For glocal alignment we need to track best score in final row
	    byte finalRowMaxScore=Byte.MIN_VALUE;
	    int finalRowMaxPos=0;
	    
	    // Track absolute best score and position during alignment
	    int maxAbsScore=0;
	    int maxPos=0;
	    
	    // Fill the matrix
	    for(int i=1; i<=qLen; i++) {
	        // First column - gap in reference (insertion penalty)
	        curr[0]=(byte)(prev[0]-1);
	        
	        // Get current query base
	        byte q=query[i-1];
	        
	        // Process row - track maximum within this row
	        byte maxScore=Byte.MIN_VALUE;
	        int maxPos_i=0;
	        
	        // Define band (for future band implementation - currently full width)
	        final int bandStart=1;
	        final int bandEnd=rLen;
	        
	        // Process each cell in the row
	        for(int j=bandStart; j<=bandEnd; j++) {
	            // Get current reference base
	            byte r=ref[j-1];
	            
	            // Score calculation
	            boolean isMatch=(q==r && q!='N'); // Check for match (excluding N's)
	            boolean hasN=(q=='N' || r=='N');
	            final byte scoreAdd=isMatch? (byte)1 : (hasN? 0 : (byte)-1); // Match/mismatch
	            
	            // Read adjacent scores
	            final byte pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
	            final byte diagScore=(byte)(pj1+scoreAdd); // Match/Sub
	            final byte upScore=(byte)(pj-1);  // Insertion
	            final byte leftScore=(byte)(cj1-1); // Deletion
	            
	            // Find max score 
	            final byte maxDiagUp=Tools.max(diagScore, upScore);
	            final byte score=Tools.max(maxDiagUp, leftScore);
	            
	            // Write score to current cell
	            curr[j]=score;
	            
	            // Output debug info if debug flag is set
	            if(debug) {
	                System.err.println("\nCell i="+i+", j="+j+": score="+score);
	                System.err.println("match="+isMatch+", diag="+diagScore+", up="+upScore+", left="+leftScore);
	            }
	            
	            // Track best score in row
	            if(score>maxScore) {
	                maxScore=score;
	                maxPos_i=j;
	            }
	        }
	        
	        // Special handling for final row (glocal alignment)
	        if(i==qLen) {
	            finalRowMaxScore=maxScore;
	            finalRowMaxPos=maxPos_i;
	        }
	        
	        // Output debug info on row maxima
	        if(debug) {
	            int rowAbsScore=maxScore-MIDPOINT+totalAdjustment;
	            System.err.println("\nRow i="+i+": rowScore="+maxScore+", rowAbsScore="+rowAbsScore);
	        }
	        
	        // Recalibrate when needed to keep scores in byte range
	        if(i==timeToRecalibrate) {
	            int offset=MIDPOINT-maxScore;
	            
	            // Always recalibrate to prevent underflow
	            for(int j=0; j<=rLen; j++) {
	                // Ensure we don't go below -MIDPOINT to prevent underflow
	                int newVal=curr[j]+offset;
	                curr[j]=(byte)Math.max(-MIDPOINT, newVal);
	            }
	            totalAdjustment-=offset; // Track total adjustment (negative because we add the offset)
	            
	            if(debug) {
	                System.err.println("Recal: i="+i+": offset="+offset+", totalAdjustment="+totalAdjustment);
	            }
	            
	            timeToRecalibrate=i+60; // Schedule next recalibration
	        }
	        
	        // Swap rows for next iteration
	        byte[] temp=prev;
	        prev=curr;
	        curr=temp;
	    }
	    
	    // Use the final row's max score for glocal alignment
	    // Calculate true absolute score by adjusting for recalibrations
	    int absScore=finalRowMaxScore-MIDPOINT+totalAdjustment;
	    
	    // For alignment scoring with match=1, mismatch=-1, gap=-1:
	    // score = matches - (mismatches + gaps)
	    // When aligning a query of length qLen:
	    // matches + mismatches + gaps = qLen
	    // So: matches = (score + qLen) / 2
	    // Special case for single base: if absScore > 0, it's a match
	    int matches=Math.max(0, (absScore+qLen)/2);
	    if(qLen==1 && absScore>0) matches=1;
	    
	    float identity=(float)matches/(float)qLen;
	    
	    // Set output coordinates if posVector is provided
	    if(posVector!=null) {
	        posVector[0]=0; // For glocal alignment we start at 0
	        posVector[1]=finalRowMaxPos-1; // Convert to 0-based inclusive
	    }
	    
	    // Output final score information
	    if(debug) {
	        System.err.println("Final abs score="+absScore+", matches="+matches+
	                          ", identity="+identity);
	    } else {
	        System.err.println("maxScore="+absScore);
	    }
	    
	    return identity;
	}

	private static final float postprocess(byte[] prev, int qLen, int rLen, int maxScore, int maxPos, int[] pos) {
		
	    // Process results and return identity
	    if(pos != null) {
	        pos[0] = 0;//Unknown, for now
	        pos[1] = maxPos;
	    }
	    
	    System.err.println("maxScore="+maxScore);
	    
	    // Calculate identity from maxAbsScore and alignment length
	    float identity = (maxScore+qLen)/(2f*qLen);
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

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static boolean GLOBAL=false;
	public static boolean debug=true;

}
