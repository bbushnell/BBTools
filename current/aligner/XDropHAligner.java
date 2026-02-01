package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import parse.PreParser;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * X-drop heuristic aligner for calculating Average Nucleotide Identity (ANI).
 * Uses a sparse dynamic programming approach with only two arrays to avoid traceback.
 * Calculates exact alignment identity and reference start/stop positions without full traceback.
 * Limited to sequences up to 2Mbp due to 21-bit position encoding.
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 30, 2025
 */
public class XDropHAligner implements IDAligner{

	/**
	 * Program entry point that delegates to Test class for standardized testing.
	 * Uses reflection to determine the actual class type from stack trace.
	 * @param <C> Class type extending IDAligner
	 * @param args Command-line arguments passed to testing framework
	 * @throws Exception If class reflection or testing fails
	 */
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

	public XDropHAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the aligner name identifier */
	@Override
	public final String name() {return "XDrop-H";}
	/**
	 * Aligns two sequences and returns identity.
	 * @param a First sequence
	 * @param b Second sequence
	 * @return Identity score from 0.0 to 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	/**
	 * Aligns two sequences and returns identity with position information.
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Array to store reference start/stop positions
	 * @return Identity score from 0.0 to 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	/**
	 * Aligns two sequences with minimum score threshold.
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Array to store reference start/stop positions
	 * @param minScore Minimum score threshold (currently unused)
	 * @return Identity score from 0.0 to 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	/**
	 * Aligns sequences within a specified reference window.
	 * @param a First sequence
	 * @param b Second sequence
	 * @param pos Array to store reference start/stop positions
	 * @param rStart Reference window start position
	 * @param rStop Reference window stop position
	 * @return Identity score from 0.0 to 1.0
	 */
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculates optimal bandwidth for alignment based on sequence characteristics.
	 * Tests for high-identity indel-free alignments requiring minimal bandwidth.
	 * Uses adaptive bandwidth based on sequence length and observed differences.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @return Bandwidth value (minimum 1)
	 */
	private static int decideBandwidth(byte[] query, byte[] ref) {
		int bandwidth=Tools.min(query.length/4+2, Math.max(query.length, ref.length)/32, 12);
		bandwidth=Math.max(2, bandwidth)+3;
		int subs=0;
		for(int i=0, minlen=Math.min(query.length, ref.length); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);//At least 1
	}

	/**
	 * Core static alignment method using X-drop heuristic with sparse DP matrix.
	 * Swaps sequences if query is longer than reference when position vector is null.
	 * Uses bit-packed scoring to track alignment score, position, and deletion count.
	 * Implements adaptive bandwidth and X-drop termination for efficiency.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional array for returning reference start/stop coordinates
	 * @return Identity score from 0.0 to 1.0
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
		
		//Create a visualizer if an output file is defined
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));

		// Matrix exploration limits
		final int bandWidth=(int)(1.4*decideBandwidth(query, ref));
//		final int topWidth=Math.min(query.length, bandWidth*2);
		//Lower insPad allows later top entrance
		final int insPad=-16-rLen/128-(int)Math.sqrt(rLen)-(5*Math.max(0, rLen-qLen))/4;
//		final int sideWidth0=1;//Set to >1 if you want a sideband.  Do NOT set >rLen.
//		final int sideWidthMax=Tools.min(qLen, rLen);
		final int rightExtend=(LOOP_VERSION ? 3 : 1);
		final long scoreWidth0=(bandWidth+1L)<<SCORE_SHIFT;//aka sandbarWalkThreshold
		
//		System.err.println("BW="+bandWidth+", topW="+topWidth+", scoreW="+(scoreWidth0>>SCORE_SHIFT));
		
		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];

		// Create IntLists for tracking active positions
		IntList activeList = new IntList(rLen+4);
		IntList nextList = new IntList(rLen+4);//aka upcomingRapids

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}

		// Initialize active list to all but first column
		for(int j=0; j<=rLen; j++) {activeList.addUnchecked(j);}
		
		//Prefill next list
		nextList.addUnchecked(0);

		int maxPos=0; // Best scoring position
		long currentMaxScore=BAD;
		long overallMaxScore=BAD;

		// Fill alignment matrix using the sparse loop
		for(int i=1; i<=qLen; i++){
			// First column
			curr[0]=i*INS;
			
			//Remove potential excess sites
			//This allows simplifying branch structure in the inner loop
			while(activeList.lastElement()>rLen) {activeList.pop();}
			mloops+=activeList.size()-1;
			if(activeList.size<2) {return 0;}//X-Drop
			
			//Clear the potential stale value in the last cell of prev.
			//This action does not get seen by the visualizer
			if(activeList.lastElement()<rLen) {prev[rLen]=BAD;}
			
			//Swap row best scores
			overallMaxScore=Math.max(overallMaxScore, currentMaxScore);
			currentMaxScore=BAD;
			maxPos=0;
			
			//Moving the sideband test outside the inner loop is faster
			final int sideWidth=1;// 0? //Tools.mid(sideWidth0, topWidth*2-i, sideWidthMax);
			assert(nextList.size()>=sideWidth || rLen<sideWidth) : "\nsize="+nextList.size+", sideW="+sideWidth;
//					+", sideW0="+sideWidth0+", qLen="+qLen+", rLen="+rLen+", "+(topWidth*2-i)+"\n"+nextList;
			nextList.size=sideWidth;//Lists are supposed to be monotonically increasing
			assert(nextList.size<=nextList.array.length) : nextList.size+", "+nextList.array.length;
			assert(nextList.lastElement()+1==sideWidth || rLen<sideWidth) : nextList+", "+sideWidth;
			
			//Allows skipping topband test
			final long scoreWidth=scoreWidth0;
			
			final byte q=query[i-1];//Cache the query
			final int forcedInsOffset=insPad-i;//j threshold for insertion penalty
			
			// Process only positions in the active list
			for(int idx=1; idx<activeList.size; idx++){
				int j = activeList.array[idx];
				
				assert(j>0) : idx+", "+j;//This can fail with a super-short ref...
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean hasN=(q=='N' || r=='N');
				final boolean isMatch=(q==r && q!='N');
				final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL_INCREMENT;

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				final long insScoreMod=Math.max(0L, (j+forcedInsOffset)/2)<<SCORE_SHIFT;
				//Changing this conditional to max or removing the mask causes a slowdown.
				final long maxValue0=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				//insScoreMod prevents false paths when ANI<60%, but eliminates pretty clouds
				final long maxValue=maxValue0-insScoreMod;

				final long scoreDif=overallMaxScore-maxValue;
				final int last=nextList.array[nextList.size-1];
				//Eliminating to topWidth test increases speed
				final boolean add=j<=rLen && (j<sideWidth || scoreDif<scoreWidth);
//				final boolean live=(EXTEND_MATCH && isMatch & last<j+1);
				
				//Important: Injecting "BAD" into these cells clears stale values.
				//Update current cell
				curr[j]=(add ? maxValue : BAD);
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
						final int jp2=j+2, jp1=j+1;
						if(last==jp1 && jp1<rLen) {//Common Case
							nextList.addUnchecked(jp2);
						}else {//Rare Case
							int tail=last;
							
							//Bounds unchecked version
							if(last<j) {nextList.addUnchecked(j); tail=j;}
							if(tail<jp1) {nextList.addUnchecked(jp1); tail=jp1;}
							if(tail<jp2) {nextList.addUnchecked(jp2); tail=jp2;}
//							if(tail<jp3) {nextList.addUnchecked(jp3);}
						}
					}
				}

				// Track best score in row
				// The mask is necessary for speed,
				// but either > or >= are OK.
				final boolean better=((maxValue&SCORE_MASK)>currentMaxScore);
				currentMaxScore=better ? maxValue : currentMaxScore;
				maxPos=better ? j : maxPos;
			}
			if(viz!=null) {viz.print(curr, activeList, rLen);}
			
			// Swap rows
			long[] temp=prev;
			prev=curr;
			curr=temp;

			// Swap position lists
			IntList tempList=activeList;
			activeList=nextList;
			nextList=tempList;
		}
		if(viz!=null) {viz.shutdown();}// Terminate visualizer
		if(GLOBAL) {maxPos=rLen;currentMaxScore=prev[rLen-1]+DEL_INCREMENT;}//The last cell may be empty 
		loops.addAndGet(mloops);
		return postprocess(currentMaxScore, maxPos, qLen, rLen, posVector);
	}
	
	/**
	 * Converts alignment score to identity and extracts position information.
	 * Solves system of equations to determine match/substitution/indel counts.
	 * Handles global alignment adjustments if GLOBAL mode is enabled.
	 * @param maxScore Highest score from final alignment row
	 * @param maxPos Position of highest score in reference
	 * @param qLen Query sequence length
	 * @param rLen Reference sequence length
	 * @param posVector Optional array for storing reference coordinates and statistics
	 * @return Final identity score from 0.0 to 1.0
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
		
		if(posVector!=null){//TODO: Enforce this as being an int[>=4], not int[2].
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
	 * Converts score array to readable string format for debugging.
	 * Extracts score bits from packed long values.
	 * @param array Array of packed score values
	 * @return Formatted string representation of scores
	 */
	private static ByteBuilder toScore(long[] array) {
		ByteBuilder bb=new ByteBuilder();
		bb.append('[');
		for(int i=0; i<array.length; i++) {
			bb.append(array[i]>>SCORE_SHIFT);
			bb.append(',');
		}
		bb.set(bb.length-1, ']');
		return bb;
	}

	/**
	 * Wrapper method for aligning to a specific window of the reference sequence.
	 * Extracts reference region and adjusts returned coordinates accordingly.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Array for storing alignment coordinates
	 * @param refStart Alignment window start position (inclusive)
	 * @param refEnd Alignment window end position (exclusive)
	 * @return Identity score from 0.0 to 1.0
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
	private static final boolean EXTEND_MATCH=true;
	private static final boolean LOOP_VERSION=false;
	private static final boolean PRINT_OPS=false;
//	private static final boolean debug=false;
	// This will force full-length alignment, but it will only be optimal
	// if the global alignment is within the glocal bandwidth.
	// Better to use Banded/Glocal for arbitrary global alignments.
	public static final boolean GLOBAL=false;

}
