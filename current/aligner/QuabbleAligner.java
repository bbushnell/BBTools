package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.PreParser;
import shared.Tools;
import structures.IntList;
import structures.RingBuffer;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback.
 *Limited to length 2Mbp with 21 position bits.
 *
 *@author Brian Bushnell
 *@contributor Isla
 *@date April 24, 2025
 */
public class QuabbleAligner implements IDAligner{

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

	/** Default constructor for QuabbleAligner instance */
	public QuabbleAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "Quabble";}
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
		int qLen=query.length, rLen=ref.length;
		int maxLen=Math.max(qLen, rLen), minLen=Math.min(qLen, rLen);
		int bandwidth=Tools.min(qLen/4+2, maxLen/32, (int)Tools.log2(maxLen+256)+2);
		bandwidth=Math.max(2, bandwidth)+2;
		int subs=0;
		for(int i=0; i<minLen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);//At least 1
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
		
		//Create a visualizer if an output file is defined
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));

		// Matrix exploration limits
		final int bandWidth=decideBandwidth(query, ref);
		final int topWidth=Math.min(query.length, bandWidth*2);
		//Lower insPad allows later top entrance
		final int insPad=-16-rLen/128-(int)Math.sqrt(rLen)-(5*Math.max(0, rLen-qLen))/4;
		final int denseWidth=DENSE_TOP ? Tools.min(topWidth, query.length-1) : 0;
		final int sideWidth0=1;//Set to >1 if you want a sideband.  Do NOT set >rLen.
		final int sideWidthMax=Tools.min(qLen, rLen);
		final int rightExtend=(LOOP_VERSION ? Math.max(5, bandWidth-2) : 2);
		final long scoreWidth0=(bandWidth+1L)<<SCORE_SHIFT;
		final int ringSize=Math.max(4, bandWidth/2);
		final RingBuffer ring=new RingBuffer(ringSize);
		
		// Create arrays for current and previous rows
		long[] prev=new long[rLen+1], curr=new long[rLen+1];

		// Create IntLists for tracking active positions
		IntList activeList = new IntList(rLen+4);
		IntList nextList = new IntList(rLen+4);

		{// Initialize first row with starting position in the lower bits
			final long mult=(GLOBAL ? DEL_INCREMENT : 1);
			for(int j=0; j<=rLen; j++){prev[j]=j*mult;}
		}
		
		final int sparseStart=1+denseWidth;
		if(denseWidth>0) {//Optionally use a dense strategy for aligning the top band
			long[][] arrays=alignDense(query, ref, prev, curr, viz, denseWidth+1, rLen);
			curr=arrays[0];
			prev=arrays[1];
		}

		// Initialize active list to all but first column
		for(int j=0; j<=rLen; j++) {activeList.addUnchecked(j);}
		
		//Prefill next list
		for(int j=0; (j<=sideWidth0 || j<=topWidth*2) && j<qLen; j++) {nextList.addUnchecked(j);}

		int bridgeTime=BRIDGE_PERIOD;
		int maxPos=0; // Best scoring position
		long maxScore=BAD;
		long prevRowScore=BAD;

		// Fill alignment matrix using the sparse loop
		for(int i=sparseStart; i<=qLen; i++){
			curr[0]=i*INS;// Fill first column
			
			//Remove potential excess sites
			//This allows simplifying branch structure in the inner loop
			while(activeList.lastElement()>rLen) {activeList.pop();}

			final byte q=query[i-1];//Cache the query
			final boolean nextMatch=(q==ref[Math.min(rLen-1, maxPos)]);
			if(BUILD_BRIDGES){//Race to catch up with long deletions
				int extra=((i&15)<1 ? 40 : 5);
				int last=activeList.lastElement();
				int eLimit=Math.min(last+extra, rLen);
				for(int e=last+1; e<eLimit; e++) {
					activeList.addUnchecked(e);
				}
			}
			mloops+=activeList.size()-1;
			
			//Clear the potential stale value in the last cell of prev.
			//This action does not get seen by the visualizer
			if(activeList.lastElement()<rLen) {prev[rLen]=BAD;}
			
			final int oldMaxScore=(int)(ring.getOldestUnchecked()>>SCORE_SHIFT);
			final int recentMissingScore=(oldMaxScore+ringSize)-(int)(maxScore>>SCORE_SHIFT)-1;
			final int scoreBonus=Math.max(0, Math.min(ringSize*2, recentMissingScore*2));
//			System.err.println("oldMaxScore="+oldMaxScore+", maxScore="
//					+(maxScore>>SCORE_SHIFT)+", scoreBonus="+scoreBonus);
			//Allows skipping topband test
//			final long scoreWidth=scoreWidth0+MATCH*((Math.max(0, topWidth-i))+scoreBonus);
			
			
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
			final long scoreWidth=scoreWidth0+MATCH*(Math.max(scoreBonus, topWidth-i));
			
			final int forcedInsOffset=insPad-i;//j threshold for insertion penalty
			
			// Process only positions in the active list
			for(int idx=1; idx<activeList.size; idx++){
				int j = activeList.array[idx];
				assert(j>0) : idx+", "+j;
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean hasN=(q=='N' || r=='N');
				final boolean isMatch=(q==r && q!='N');
				final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				// Read adjacent scores
				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore1=cj1+DEL_INCREMENT;
				
				//Allows a long deletion - this is the quantum teleportation feature allowing
				//jumps between high-scoring regions, across an unexplored gap
				final long leftScore=Math.max(leftScore1, maxScore+DEL_INCREMENT*(j-maxPos));

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				final long insScoreMod=Math.max(0L, (j+forcedInsOffset)/2)<<SCORE_SHIFT;
				//Changing this conditional to max or removing the mask causes a slowdown.
				final long maxValue0=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;
				//insScoreMod prevents false paths when ANI<60%, but eliminates pretty clouds
				final long maxValue=maxValue0-insScoreMod;

				final long scoreDif=prevRowScore-maxValue;
				final int last=nextList.array[nextList.size-1];
				//Eliminating to topWidth test increases speed
				final boolean add=j<=rLen && (/*i<topWidth ||*/ j<sideWidth || scoreDif<scoreWidth);
//						|| (GLOBAL && j>=i && qLen-i<100)); This is only marginally useful 
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
				// The mask is necessary for speed,
				// but either > or >= are OK.
				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			if(viz!=null) {viz.print(curr, activeList, rLen);}
			ring.add(maxScore);

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
		if(GLOBAL) {maxPos=rLen;maxScore=prev[rLen-1]+DEL_INCREMENT;}//The last cell may be empty 
		loops.addAndGet(mloops);
		return postprocess(maxScore, maxPos, qLen, rLen, posVector);
	}
	
	// Process the first topWidth rows using a dense approach
	/**
	 * Processes the first topWidth rows using dense dynamic programming approach.
	 * Fills all cells in the top portion of the alignment matrix.
	 * Used as initialization before switching to sparse algorithm.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param prev Previous row scores
	 * @param curr Current row scores
	 * @param viz Optional visualizer for debugging
	 * @param topWidth Number of top rows to process densely
	 * @param rLen Reference sequence length
	 * @return Array containing [current_row, previous_row] after processing
	 */
	private static final long[][] alignDense(byte[] query, byte[] ref, long[] prev, 
			long[] curr, Visualizer viz, int topWidth, int rLen) {
		long mloops=0;
		for(int i=1; i<topWidth; i++) {
			// First column penalizes clipping
			curr[0]=i*INS;
			
			//Cache the query
			final byte q=query[i-1];

			// Process all columns in top rows
			for(int j=1; j<=rLen; j++) {
				final byte r=ref[j-1];

				// Branchless score calculation
				final boolean hasN=(q=='N' || r=='N');
				final boolean isMatch=(q==r && q!='N');
				final long scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				// Calculate scores
				final long pj1=prev[j-1], pj=prev[j];
				final long cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;// Match/Sub
				final long upScore=pj+INS;
				final long leftScore=cj1+DEL_INCREMENT;

				// Find max using conditional expressions
				final long maxDiagUp=Math.max(diagScore, upScore);//This is fine
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;

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
	 * Use alignment information to calculate identity and starting coordinate.
	 * @param maxScore Highest score in last row
	 * @param maxPos Highest-scoring position in last row
	 * @param qLen Query length
	 * @param rLen Reference length
	 * @param posVector Optional array for returning reference start/stop coordinates.
	 * @return Identity
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
		assert(posVector[1]>0) : id+", "+Arrays.toString(posVector)+", "+refStart;
		if(posVector!=null) {
			posVector[0]+=refStart;
			posVector[1]+=refStart;
		}
		return id;
	}

	/** Thread-safe counter for total inner loop iterations across all alignments */
	private static AtomicLong loops=new AtomicLong(0);
	/**
	 * Returns the total number of inner loop iterations performed across all alignments
	 */
	public long loops() {return loops.get();}
	/** Sets the loop counter to a specific value.
	 * @param x New loop count value */
	public void setLoops(long x) {loops.set(x);}
	/** Optional output file path for visualization debugging */
	public static String output=null;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Bit field definitions
	/** Number of bits allocated for position information in score encoding */
	private static final int POSITION_BITS=21;
	/** Number of bits allocated for deletion count in score encoding */
	private static final int DEL_BITS=21;
	/** Bit shift amount for score component in packed long values */
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;

	// Masks
	/** Bit mask for extracting position information from packed scores */
	private static final long POSITION_MASK=(1L << POSITION_BITS)-1;
	/** Bit mask for extracting deletion count from packed scores */
	private static final long DEL_MASK=((1L << DEL_BITS)-1) << POSITION_BITS;
	/** Bit mask for extracting raw alignment score from packed values */
	private static final long SCORE_MASK=~(POSITION_MASK | DEL_MASK);

	// Scoring constants
	/** Score increment for sequence matches */
	private static final long MATCH=1L << SCORE_SHIFT;
	/** Score penalty for substitutions */
	private static final long SUB=(-1L) << SCORE_SHIFT;
	/** Score penalty for insertions */
	private static final long INS=(-1L) << SCORE_SHIFT;
	/** Score penalty for deletions */
	private static final long DEL=(-1L) << SCORE_SHIFT;
	/** Score for ambiguous bases (N characters) */
	private static final long N_SCORE=0L;
	/** Sentinel value indicating invalid or uninitialized alignment scores */
	private static final long BAD=Long.MIN_VALUE/2;
	/** Combined deletion penalty and position increment for packed scoring */
	private static final long DEL_INCREMENT=DEL+(1L<<POSITION_BITS);

	// Run modes
	/** Whether to extend alignment exploration from matching cells */
	private static final boolean EXTEND_MATCH=true;
	/** Whether to use loop-based or branchless version for position updates */
	private static final boolean LOOP_VERSION=false;
	/** Whether to build bridges across gaps to catch long deletions */
	private static final boolean BUILD_BRIDGES=true;
	/** Frequency of bridge building operations */
	private static final int BRIDGE_PERIOD=16;
	/** Whether to use dense alignment for top rows */
	private static final boolean DENSE_TOP=false;
	/** Whether to print detailed operation statistics for debugging */
	private static final boolean PRINT_OPS=false;
//	private static final boolean debug=false;
	// This will force full-length alignment, but it will only be optimal
	// if the global alignment is within the glocal bandwidth.
	// Better to use Banded/Glocal for arbitrary global alignments.
	/** Whether to force global alignment mode instead of glocal */
	public static final boolean GLOBAL=false;

}
