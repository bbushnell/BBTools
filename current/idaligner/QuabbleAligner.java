package idaligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import parse.PreParser;
import shared.Tools;
import structures.IntList;
import structures.RingBuffer;

/**
 * Fast sequence aligner using sparse dynamic programming for ANI calculation.
 * Uses only two arrays and avoids traceback, calculating reference start and
 * stop positions directly from packed scores. Limited to sequences up to
 * 2Mbp with 21 position bits and tuned for high-identity alignments via
 * adaptive bandwidth selection and bridge building.
 * @author Brian Bushnell
 * @contributor Isla
 * @date April 24, 2025
 */
public class QuabbleAligner implements IDAligner{

	/**
	 * Program entry point that delegates to the Test harness.
	 * Uses reflection to determine the calling aligner class and passes it to Test.
	 * @param args Command-line arguments forwarded to Test.testAndPrint()
	 * @throws Exception If class lookup or test execution fails
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

	@Override
	public final float align(byte[] a, byte[] b, AlignmentStats stats){
		float id=alignStatic(a, b, null);
		if(stats!=null){
			int[] pos=new int[4];
			id=alignStatic(a, b, pos);
			stats.setFromPos(pos, id);
			stats.qLen=a.length;
			stats.rLen=b.length;
		}
		return id;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Determines the bandwidth for alignment based on sequence lengths and an
	 * early mismatch scan, choosing a narrow band for high-identity, low-indel
	 * alignments.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @return Bandwidth for alignment, at least 2 plus a small safety margin
	 */
	private static int decideBandwidth(byte[] query, byte[] ref, int[] pos) {
		int qLen=query.length, rLen=ref.length;
		int maxLen=Math.max(qLen, rLen);
		int bandwidth=Tools.min(qLen/4+2, maxLen/32, (int)Tools.log2(maxLen+256)+2);
		bandwidth=Math.max(2, bandwidth)+2;
		return IDAlignerStatics.decideBandwidth(query, ref, pos, bandwidth);
	}

	/**
	 * Main static alignment method using sparse dynamic programming with packed
	 * scores, optional dense initialization, and a ring buffer–based heuristic to
	 * encourage recent improvements.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2+] for returning reference start/stop and stats
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
		if(posVector==null) {posVector=new int[2];}
		final int bandWidth=decideBandwidth(query, ref, posVector);
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
		return Tracer.postprocess(maxScore, maxPos, qLen, rLen, posVector, null);
	}

	// Process the first topWidth rows using a dense approach
	/**
	 * Processes the first topWidth rows using a dense dynamic programming
	 * approach that fills all cells in the top portion of the alignment matrix
	 * before switching to the sparse algorithm.
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
	 * Lightweight wrapper for aligning to a window of the reference sequence.
	 * Extracts the specified region, calls the core alignStatic, asserts a
	 * successful alignment, and adjusts returned coordinates back to the full
	 * reference.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional array for returning reference coordinates
	 * @param refStart Start position of alignment window (inclusive)
	 * @param refEnd End position of alignment window (inclusive)
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
	private static final boolean BUILD_BRIDGES=true;
	private static final int BRIDGE_PERIOD=16;
	private static final boolean DENSE_TOP=false;
	private static final boolean PRINT_OPS=false;
//	private static final boolean debug=false;
	// This will force full-length alignment, but it will only be optimal
	// if the global alignment is within the glocal bandwidth.
	// Better to use Banded/Glocal for arbitrary global alignments.
	public static final boolean GLOBAL=false;

}
