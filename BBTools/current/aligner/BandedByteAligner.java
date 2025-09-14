package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.PreParser;
import shared.Shared;

/**
 *Aligns two sequences to return approximate ANI.
 *Uses only 2 arrays.
 *Requires traceback.
 *Limited to length 2Mbp with 21 position bits.
 *Restricts alignment to a fixed band around the diagonal.
 *
 *@author Brian Bushnell
 *@contributor Isla
 *@date May 4, 2025
 */
public class BandedByteAligner implements IDAligner{

	/** Main() passes the args and class to Test to avoid redundant code */
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

	public BandedByteAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "BandedByte";}
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

		final byte[] query=Factory.encodeByte(query0, (byte)15);
		final byte[] ref=Factory.encodeByte(ref0, (byte)31);
		
		final int qLen=query.length;
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, 0, 0));
		
		// Banding parameters
		final int bandWidth=decideBandwidth(query, ref);
		// Initialize band limits for use outside main loop
		int bandStart=1, bandEnd=rLen;

		// Create arrays for current and previous rows
		byte[] prev=new byte[rLen+1], curr=new byte[rLen+1];
		Arrays.fill(prev, (byte)60);
		Arrays.fill(curr, BAD);
		
		int timeSinceBalance=0;
		int absScore=0;
		
		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
			// Calculate band boundaries 
			bandStart=Math.max(1, Math.min(i-bandWidth, rLen-bandWidth));
			bandEnd=Math.min(rLen, i+bandWidth);
			
            if(debug) {
                System.err.println("\nRow "+i+": bw="+bandWidth+"; j is "+bandStart+" to "+bandEnd);
            }
			
			//Clear stale data to the left of the band
			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=0;//(byte)(i*INS);//TODO: This needs to be changed for relative scoring
			
			//Cache the query
			final byte q=query[i-1];
			
			if(Shared.SIMD) {
				shared.SIMDAlignByte.alignBandVector(q, ref, bandStart, bandEnd, prev, curr);
			}else {

				// Process only cells within the band
				for(int j=bandStart; j<=bandEnd; j++){
					final byte r=ref[j-1];

					// Branchless score calculation
					final boolean isMatch=(q==r);
					final boolean hasN=((q|r)>=15);
					final byte scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

					// Read adjacent scores
					final byte pj1=prev[j-1], pj=prev[j];
					final byte diagScore=(byte)(pj1+scoreAdd);// Match/Sub
					final byte upScore=(byte)(pj+INS);

					// Find max using conditional expressions
					final byte maxDiagUp=(byte)Math.max(diagScore, upScore);//This is fine

					curr[j]=maxDiagUp;
		            
		            // Output debug info if debug flag is set
		            if(debug) {
		                System.err.println("\nCell i="+i+", j="+j+": maxDiagUp="+maxDiagUp);
		                System.err.println("pj1="+pj1+", pj="+pj+", scoreAdd="+scoreAdd+", INS="+INS);
		                System.err.println("match="+isMatch+", diag="+diagScore+", up="+upScore);
		            }
				}
			}
			
			if(Shared.SIMD) {
				shared.SIMDAlignByte.processDeletionsTailVector(curr, bandStart, bandEnd);
			}else {
				//Tail loop for deletions
				byte leftCell=curr[bandStart-1];
				for(int j=bandStart; j<=bandEnd; j++){
					final byte maxDiagUp=curr[j];
					final byte leftScore=(byte)(leftCell+DEL);
					leftCell=(byte)Math.max(maxDiagUp, leftScore);
					curr[j]=leftCell;


					// Output debug info if debug flag is set
					if(debug) {
						System.err.println("\nCell i="+i+", j="+j+": score="+leftCell);
						System.err.println("maxDiagUp="+maxDiagUp+", leftScore="+leftScore);
					}
				}
			}
	        
	        // Output debug info on row maxima
	        if(debug) {//TODO
	        	System.err.println("\nRow "+i+" scores: "+Arrays.toString(curr));
//	            int rowAbsScore=maxScore-MIDPOINT+totalAdjustment;
//	            System.err.println("\nRow i="+i+": rowScore="+maxScore+", rowAbsScore="+rowAbsScore);
	        }
	        
	        timeSinceBalance++;
	        if(timeSinceBalance>=60) {
	        	absScore=rebalance(curr, bandStart, bandEnd, absScore);
	        	timeSinceBalance=0;
	        }
			
			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			mloops+=(bandEnd-bandStart+1);

			// Swap rows
			byte[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		loops.addAndGet(mloops);
		return postprocess(prev, qLen, bandStart, bandEnd, posVector, absScore);
	}
	
	private static int rebalance(byte[] curr, int bandStart, int bandEnd, int oldMax) {
		int max=-127;
		for(int j=bandStart; j<=bandEnd; j++) {
			max=Math.max(max, curr[j]);
		}
		byte delta=(byte)(60-max);//Negative if max is above 60
		for(int j=bandStart; j<=bandEnd; j++) {
			curr[j]=(byte)Math.max(-60, curr[j]+delta);
		}
		return oldMax-delta;
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
	private static final float postprocess(byte[] prev, int qLen, 
			int bandStart, int bandEnd, int[] posVector, int absScore) {
    	absScore=rebalance(prev, bandStart, bandEnd, absScore);
        
		// Find best score outside of main loop
		long maxScore=Long.MIN_VALUE;
		int maxPos=bandEnd;
		for(int j=bandStart; j<=bandEnd; j++){
			long score=prev[j];
			if(score>maxScore){
				maxScore=score;
				maxPos=j;
			}
		}
		maxScore=absScore;
        if(debug) {
            System.err.println("\npostprocess: prev="+Arrays.toString(prev)+
            		", qLen="+qLen+", bandStart="+bandStart+", bandEnd="+bandEnd);
            System.err.println("maxScore="+maxScore+", maxPos="+maxPos);
        }

		// Extract alignment information
		final int originPos=0;//Unknown
		final int endPos=maxPos;
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
		}

		float identity=(maxScore+qLen)/(2f*qLen);
		
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

	// Scoring constants
	private static final byte MATCH=1;
	private static final byte SUB=-1;
	private static final byte INS=-1;
	private static final byte DEL=-1;
	private static final byte N_SCORE=0;
	private static final byte BAD=Byte.MIN_VALUE/2;

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static final boolean debug=false;
//	public static boolean GLOBAL=false; //Cannot handle global alignments

}
