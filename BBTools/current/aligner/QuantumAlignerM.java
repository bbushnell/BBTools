package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Timer;
import shared.Tools;
import structures.IntList;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback.
 *Limited to length 2Mbp with 21 position bits.
 *
 *Counts matches instead of dels.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 24, 2025
 */
public class QuantumAlignerM implements IDAligner{

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

	public QuantumAlignerM() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "QuantumM";}
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
		int bandwidth=Tools.min(query.length/2, Math.max(query.length, ref.length)/16, 24);
		bandwidth=Math.max(8, bandwidth);
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
		
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, MATCH_BITS));

		// Matrix exploration limits
		final int topWidth=decideBandwidth(query, ref);
		final int sideWidth0=1;//Set to >1 if you want a sideband.  Do NOT set >rLen.
		final int sideWidthMax=Tools.min(qLen, rLen);
		final int rightExtend=(LOOP_VERSION ? Math.max(5, topWidth/2-2) : 2);
		final long scoreWidth0=((long)(Math.max(8, topWidth/2)))<<SCORE_SHIFT;

		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];

		// Create IntLists for tracking active positions
		IntList activeList = new IntList(rLen+3);
		IntList nextList = new IntList(rLen+3);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}
		
		final int sparseStart=DENSE_TOP ? topWidth : 1;
		if(DENSE_TOP) {//Optionally use a dense strategy for aligning the top band
			long[][] arrays=alignDense(query, ref, prev, curr, viz, topWidth, rLen);
			curr=arrays[0];
			prev=arrays[1];
		}

		// Initialize active list to all but first column
		for(int j=0; j<=rLen; j++) {activeList.add(j);}
		
		//Prefill next list
		for(int j=0; (j<=sideWidth0 || j<=topWidth*2) && j<qLen; j++) {nextList.add(j);}

		int maxPos=0; // Best scoring position
		long maxScore=BAD;
		long prevRowScore=BAD;

		// Fill alignment matrix using the sparse loop
		for(int i=sparseStart; i<=qLen; i++){
			// First column
			curr[0]=i*INS;
			
			//Remove potential excess sites
			//This allows simplifying branch structure in the inner loop
			while(activeList.lastElement()>rLen) {activeList.pop();}
			
			if(BUILD_BRIDGES){//Race to catch up with long deletions
				int extra=((i&15)<1 ? 40 : 5);
				int last=activeList.lastElement();
				int eLimit=Math.min(last+extra, rLen);
				for(int e=last+1; e<eLimit; e++) {
					activeList.add(e);
				}
			}
			mloops+=activeList.size()-1;
//			if(i<50) {System.err.println("Row "+i+": "+activeList.size()+" D");}
			
			//Clear the potential stale value in the last cell of prev.
			//This action does not get seen by the visualizer
			if(activeList.lastElement()<rLen) {prev[rLen]=BAD;}
			
			//Swap row best scores
			prevRowScore=maxScore;
			maxScore=BAD;
			maxPos=0;

//			// Clear next positions list and add first column
//			nextList.clear();
//			nextList.add(1);
			
			//Moving the sideband test outside the inner loop is faster
			final int sideWidth=Tools.mid(sideWidth0, topWidth*2-i, sideWidthMax);
			assert(nextList.size()>=sideWidth);
			nextList.size=sideWidth;
			assert(nextList.lastElement()+1==sideWidth) : nextList+", "+sideWidth;
			
			//Allows skipping topband test
			final long scoreWidth=scoreWidth0+MATCH*(Math.max(0, topWidth-i));
			
			//Cache the query
			final byte q=query[i-1];

			// Process only positions in the active list
			for(int idx=1; idx<activeList.size; idx++){
				int j = activeList.array[idx];
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');
				final long scoreAdd=isMatch ? MATCH_INCREMENT : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore1=cj1+DEL;
				
				//Allows a long deletion - this is the quantum teleportation feature allowing
				//jumps between high-scoring regions, across an unexplored gap
				final long leftScore=Math.max(leftScore1, maxScore+INS*(j-maxPos));

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				//This mask and conditional is no longer needed in the match-counting version.
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
//				final long maxValue=Math.max(maxDiagUp, leftScore);

				final long scoreDif=prevRowScore-maxValue;
				final int last=nextList.array[nextList.size-1];
				//Eliminating to topWidth test increases speed
				final boolean add=j<=rLen && (/*i<topWidth ||*/ j<sideWidth || scoreDif<scoreWidth);
				final boolean live=(EXTEND_MATCH && isMatch & last<j+1);

				//Important: Injecting "BAD" into these cells clears stale values.
				//Update current cell
				curr[j]=(add || live ? maxValue : BAD);
				//Clear previous row
				//Required for correctness but has little impact in practice
				prev[j-1]=BAD;

				// Conditionally add positions
				if(add) {
					if(LOOP_VERSION) {
						final int from = Math.max(last+1, j);
						final int to = Math.min(j+rightExtend, rLen);
						for(int k=from; k<=to; k++) {nextList.addUnchecked(k);}
					}else {
						//Loop-free version is much faster
						final int jp2=j+2, jp3=j+3;
						if(last==jp2 && jp2<rLen) {//Common Case
							nextList.addUnchecked(jp3);
						}else {//Rare Case
							final int jp1=j+1;
							int tail=last;
							
							//Bounds unchecked version
							if(last<j) {nextList.addUnchecked(j); tail=j;}
							if(tail<jp1) {nextList.addUnchecked(jp1); tail=jp1;}
							if(tail<jp2) {nextList.addUnchecked(jp2); tail=jp2;}
							if(tail<jp3) {nextList.addUnchecked(jp3);}
						}
					}
				}
				else if(live) {//Extend from matching cells; for finding deletions
					nextList.addUnchecked(j+1);
				}

				// Track best score in row
				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			if(viz!=null) {viz.print(curr, activeList, rLen);}

			// Swap rows
			long[] temp=prev;
			prev=curr;
			curr=temp;

			// Swap position lists
			IntList tempList = activeList;
			activeList = nextList;
			nextList = tempList;
		}
		
		// Terminate visualizer
		if(viz!=null) {viz.shutdown();}

		loops.addAndGet(mloops);
		return postprocess(maxScore, maxPos,  qLen, posVector);
	}
	
	
	private static float postprocess(long maxScore, int maxPos, int qLen, int[] posVector) {
		// Extract alignment information
		final int originPos=(int)(maxScore&POSITION_MASK);
		final int endPos=maxPos;
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
		}

		// Calculate alignment statistics
		final int matches=(int)((maxScore & MATCH_MASK) >> POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(maxScore >> SCORE_SHIFT);

		// Solve the system of equations:
		// 1. M + S + I = qLen
		// 2. M + S + D = refAlnLength
		// 3. Score = M - S - I - D

		// From equations 1 and 2:
		// I - D = qLen - refAlnLength
		final int iMinusD=qLen-refAlnLength;

		// From equation 2:
		// S + D = refAlnLength - M
		final int sPlusD=refAlnLength-matches;
		
		final int deletions=Math.max(0, (2*matches-rawScore-qLen));

		// Calculate operation counts
		final int substitutions=sPlusD-deletions;
		final int insertions=iMinusD+deletions;
		final float identity=matches/(float)(matches+substitutions+insertions+deletions);

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
	
	
	// Process the first topWidth rows using a dense approach
	private static final long[][] alignDense(byte[] query, byte[] ref, long[] prev, long[] curr, Visualizer viz, int topWidth, int rLen) {
		long mloops=0;
		for(int i=1; i<topWidth; i++) {
			// First column should stay at zero in dense section
			// Oops! This is no longer true, now it is i*INS 
			curr[0]=i*INS;
			
			//Cache the query
			final byte q=query[i-1];

			// Process all columns in top rows
			for(int j=1; j<=rLen; j++) {
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean hasN=(q=='N' || r=='N');
				final boolean isMatch=(q==r && q!='N');
				final long scoreAdd=isMatch ? MATCH_INCREMENT : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL;

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				//This mask and conditional is no longer needed in the match version.
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
//				final long maxValue=Math.max(maxDiagUp, leftScore);

				// Update current cell
				curr[j]=maxValue;
			}

			if(viz!=null) {viz.print(curr, null, rLen);}

			// Swap rows
			long[] temp=prev;
			prev=curr;
			curr=temp;

			// Count loops for analysis
			mloops+=rLen;
		}
		loops.addAndGet(mloops);
		
		return new long[][] {curr, prev};
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
	private static final int MATCH_BITS=21;
	private static final int SCORE_SHIFT=POSITION_BITS+MATCH_BITS;

	// Masks
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	private static final long MATCH_MASK=((1L << MATCH_BITS)-1) << POSITION_BITS;
	private static final long SCORE_MASK=~(POSITION_MASK | MATCH_MASK);

	// Scoring constants
	private static final long MATCH=1L << SCORE_SHIFT;
	private static final long SUB=(-1L) << SCORE_SHIFT;
	private static final long INS=(-1L) << SCORE_SHIFT;
	private static final long DEL=(-1L) << SCORE_SHIFT;
	private static final long N_SCORE=0L;
	private static final long BAD=Long.MIN_VALUE/2;
	private static final long MATCH_INCREMENT=MATCH+(1L<<POSITION_BITS);

	// Run modes
	private static final boolean EXTEND_MATCH=true;
	private static final boolean LOOP_VERSION=false;
	private static final boolean BUILD_BRIDGES=true;
	private static final boolean DENSE_TOP=true;
	private static final boolean PRINT_OPS=false;
	public static final boolean GLOBAL=false;

}
