package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Timer;
import shared.Tools;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Banded for speed.
 *Identity is approximate.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 19, 2025
 */
public class BandedAlignerInt implements IDAligner{

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

	public BandedAlignerInt() {}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "BandedInt";}
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
		int bandwidth=Math.min(100, 4+Math.max(query.length, ref.length)/8);
		int subs=0;
		for(int i=0, minlen=Math.min(query.length, ref.length); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);
	}

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
		final int rLen=ref.length;
		long mloops=0;
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, 0));
		
		// Banding parameters
		final int bandWidth=decideBandwidth(query, ref);

		// Create two arrays for current and previous rows
		int[] prev=new int[rLen+1], curr=new int[rLen+1];
		Arrays.fill(curr, BAD);

		{// Initialize first row with starting position in the lower bits
			final int mult=(GLOBAL ? DEL : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}

		// Initialize band limits for use outside main loop
		int bandStart=0, bandEnd=rLen-1;
		
		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
			// Calculate band boundaries 
			bandStart=Tools.max(1, Tools.min(i-bandWidth, rLen-bandWidth));
			bandEnd=Math.min(rLen, i+bandWidth);
			
			//Clear stale data to the left of the band
			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=i*INS;
			
			//Cache the query
			final byte q=query[i-1];
			
			// Process only cells within the band
			for(int j=bandStart; j<=bandEnd; j++){
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final int scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final int pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final int diagScore=pj1+scoreAdd;// Match/Sub
				final int upScore=pj+INS;
				final int leftScore=cj1+DEL;

				// Find max using conditional expressions
				final int maxDiagUp=Math.max(diagScore, upScore);//This is fine
				//Changing this conditional to max or removing the mask causes a slowdown.
				final int maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				
				curr[j]=maxValue;
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

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static final boolean GLOBAL=false;


}
