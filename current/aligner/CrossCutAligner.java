package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 *Aligns two sequences to return ANI.
 *Uses 3 scoring arrays and avoids traceback.
 *Gives an exact identity plus rstart and rstop.
 *Limited to length 2Mbp with 21 position bits.
 *Iterates over diagonals that span bottom left to top right.
 *This avoids all inter-loop data dependencies.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date May 2, 2025
 */
public class CrossCutAligner implements IDAligner{

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

	public CrossCutAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "CrossCut";}
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
	public static final float alignStatic(byte[] query0, byte[] ref0, int[] posVector) {
		if(posVector==null && query0.length>ref0.length) { // Swap sequences if query is longer than ref
			byte[] temp=query0; // Store query in temp
			query0=ref0; // Set query to ref
			ref0=temp; // Set ref to original query
		}

		final int qLen=query0.length; // Store query length
		final int rLen=ref0.length; // Store reference length
		long mloops=0;
		
		//Create a visualizer if an output file is defined
		Visualizer viz=(output==null ? null : new Visualizer(output, POSITION_BITS, DEL_BITS));

		final long[] ref=Factory.encodeLong(ref0, (byte)(15+32));
		final long[] revQuery=Factory.encodeLong(query0, (byte)(15+16));
		Tools.reverseInPlace(revQuery);

		final long[] bottom=new long[rLen+1]; // Array for bottom row

		if(qLen==0 || rLen==0) { // Handle empty sequences
			return postprocess(bottom, qLen, rLen, posVector); // Process and return
		}

		final int maxDiagLen=Math.max(qLen, rLen)+1; // Max diagonal length
		long[] diag_km2=new long[maxDiagLen]; // Array for diagonal k-2
		long[] diag_km1=new long[maxDiagLen]; // Array for diagonal k-1
		long[] diag_k=new long[maxDiagLen]; // Array for current diagonal k

		// Process all diagonals
		for(int k=2; k<=qLen+rLen; k++) {
			final int km1=k-1, km2=k-2;
			if(debug) {System.err.println("\nLooping because k="+k+"<="+(qLen+rLen));}

			// Min and max matrix coordinates for this diagonal
			int minRow=Math.max(1, k-rLen); // Min row value for this diagonal
			int maxRow=Math.min(qLen, km1); // Max row value for this diagonal

			// 1. Handle top row cell (if this diagonal intersects top row)
			if(minRow==1) {//Currently structured for maximal clarity, but can be simplified
			    final int row=1; // Top row
			    final int col=km1; // Column in top row
			    if(col>=1 && col<=rLen) { // Ensure valid column
			    	//Note revQuery[qLen-row] is a constant, queryBase0
			        handleTop(revQuery[qLen-row], ref[col-1], qLen, col, diag_k, diag_km1, bottom);
			    }
			}

			// 2. Handle left column cell
			if(k<=rLen+1) {//Currently structured for maximal clarity, but can be simplified
			    final int col=1; // Left column
			    final int row=km1; // Row in left column
			    if(row>=1 && row<=qLen) { // Ensure valid row
			    	//Note ref[col-1] is a constant, refBase0
			        handleLeft(revQuery[qLen-row], ref[col-1], qLen, row, diag_k, diag_km1, bottom);
			    }
			}

			// 3. Process inner cells of the diagonal
			if(debug) {System.err.println("Main loop start");}

			final int innerMinCol=Math.max(2, k-qLen);
			final int innerMaxCol=Math.min(rLen, k-2);

			// Check if bottom row is part of this diagonal
			final boolean processesBottomRow=(qLen>=minRow && qLen<=maxRow);
			final int bottomRowCol=(processesBottomRow) ? k-qLen : -1;

			if(Shared.SIMD) {
				shared.SIMDAlign.processCrossCutDiagonalSIMD(revQuery, ref, k, 
				        innerMinCol, innerMaxCol, qLen,
				        diag_km2, diag_km1, diag_k);
			}else {
				for(int col=innerMinCol; col<=innerMaxCol; col++) {
					final int row=k-col; // Calculate row from k and col

					assert(row>=2 && row<=qLen); // Verify row bounds

					// Calculate cell value
					final long q=revQuery[qLen-row]; // Get query base
					final long r=ref[col-1]; // Get ref base

					final boolean isMatch=(q==r); // Check if match
					final boolean hasN=(q|r)>15; // Check if has N
					long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB); // Score to add

					// Get adjacent cell values from previous diagonals
					long diagValue=diag_km2[col-2]; // Diagonal from (row-1, col-1)
					long upValue=diag_km1[col-1]; // Up from (row-1, col)
					long leftValue=diag_km1[col-2]; // Left from (row, col-1)

					// Calculate scores directly
					long diagScore=diagValue+scoreAdd; // Add score to diagonal
					long upScore=upValue+INS; // Add insertion penalty to up
					long leftScore=leftValue+DEL_INCREMENT; // Add deletion penalty to left

					// Find max using conditional expressions
					final long maxDiagUp=Math.max(diagScore, upScore); // Max of diagonal and up
					final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore?maxDiagUp:leftScore; // Compare with left

					// Store at diagonal array position (col-1)
					diag_k[col-1]=maxValue;

					if(debug) {
						System.err.println("Cell ("+row+","+col+") calculation:");
						System.err.println("  q="+q+", r="+r+", match="+isMatch);
						System.err.println("  diagValue="+diagValue+", upValue="+upValue+", leftValue="+leftValue);
						System.err.println("  diagScore="+diagScore+", upScore="+upScore+", leftScore="+leftScore);
						System.err.println("*Setting ("+row+","+col+") (diag_k["+(col-1)+"]) to "+maxValue+" ("+(maxValue>>SCORE_SHIFT)+")");
					}
				}
			}

			// Update bottom row outside the loop if it was processed in this diagonal
			if(processesBottomRow && bottomRowCol>=innerMinCol && bottomRowCol<=innerMaxCol) {
			    bottom[bottomRowCol]=diag_k[bottomRowCol-1];
			    if(debug) {System.err.println("Setting bottom["+bottomRowCol+"] to "+
			    		bottom[bottomRowCol]+" ("+(bottom[bottomRowCol]>>SCORE_SHIFT)+")");}
			}

			if(debug) {System.err.println("Main loop stop");}

			// Log diagonal arrays
			if(debug) {
				System.err.println("diag_km2: "+Arrays.toString(diag_km2));
				System.err.println("diag_km1: "+Arrays.toString(diag_km1));
				System.err.println("diag_k:   "+Arrays.toString(diag_k));
				System.err.println("Bottom:   "+Arrays.toString(bottom));
			}

			// Count loops for performance tracking
			mloops+=maxRow-minRow+1; // Count cells processed
			if(viz!=null) {viz.print(diag_k, innerMinCol-1, innerMaxCol-1, maxDiagLen);}

			// Rotate diagonals
			long[] temp=diag_km2; // Store oldest diagonal
			diag_km2=diag_km1; // Shift km1 to km2
			diag_km1=diag_k; // Shift k to km1
			diag_k=temp; // Reuse oldest array for new diagonal
		}
		if(viz!=null) {viz.shutdown();}
		loops.addAndGet(mloops);
		return postprocess(bottom, qLen, rLen, posVector); // Process and return final result
	}
	
	private static final long handleTop(long q, long r, int qLen, int col,
			long[] diag_k, long[] diag_km1, long[] bottom) {
		final int row=1;//This function is only for row 1

		// Get adjacent cell values
		final long diagValue=col-1; // Diagonal from top edge
		final long upValue=col; // Up from top edge
		final long leftValue=(col==1)?INS:diag_km1[col-2]; // Left from prev diagonal or edge

		// Calculate cell value and store
		final long maxValue=calculateCellValue(q, r, diagValue, upValue, leftValue);
		diag_k[col-1]=maxValue;

		if(debug) {
			System.err.println("Handle top:");
			System.err.println("q="+q+", r="+r+", match="+(q==r));
			System.err.println("*Setting ("+1+","+col+") (diag_k["+(col-1)+"]) to "+maxValue+" ("+(maxValue>>SCORE_SHIFT)+")");
		}

		// Update bottom row if query length is 1
		// This is such a rare case it's annoying to pollute the code with it...
		// Eventually it's best to special-case queries <2bp prior to the loop
		if(qLen==1) {
			bottom[col]=maxValue;
			if(debug) {System.err.println("Handled top: bottom["+col+"]="+maxValue);}
		}
		return maxValue;
	}
	
	private static final long handleLeft(long q, long r, int qLen, int row,
			long[] diag_k, long[] diag_km1, long[] bottom) {
		final int col=1;//This function is only for column 1
		
        // Get adjacent cell values
        final long leftValue=row*INS; // Left from left edge
        final long diagValue=leftValue-INS; // Diagonal from left edge
        final long upValue=(row==1)?1:diag_km1[0]; // Up from prev diagonal or edge

        // Calculate cell value and store
        final long maxValue=calculateCellValue(q, r, diagValue, upValue, leftValue);
        diag_k[0]=maxValue;

        if(debug) {
            System.err.println("Handle left:");
            System.err.println("q="+q+", r="+r+", match="+(q==r));
            System.err.println("*Setting ("+row+","+col+") (diag_k[0]) to "+
            		maxValue+" ("+(maxValue>>SCORE_SHIFT)+")");
        }

        // Update bottom row if this is the last row
        if(row==qLen) {
            bottom[col]=maxValue;
            if(debug) {System.err.println("Handled side: bottom["+col+"]="+maxValue);}
        }
		return maxValue;
	}

	/**
	 * Calculate score for a cell in the alignment matrix
	 * 
	 * @param q Query base
	 * @param r Reference base
	 * @param diagValue Diagonal cell value
	 * @param upValue Up cell value
	 * @param leftValue Left cell value
	 * @return Maximum score for this cell
	 */
	private static long calculateCellValue(long q, long r, long diagValue, long upValue, long leftValue) {
	    final boolean isMatch=(q==r); // Check if match
	    final boolean hasN=(q|r)>15; // Check if has N
	    final long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB); // Score to add

	    // Calculate scores directly
	    final long diagScore=diagValue+scoreAdd; // Add score to diagonal
	    final long upScore=upValue+INS; // Add insertion penalty to up
	    final long leftScore=leftValue+DEL_INCREMENT; // Add deletion penalty to left

	    // Find max using conditional expressions
	    final long maxDiagUp=Math.max(diagScore, upScore); // Max of diagonal and up
	    return (maxDiagUp&SCORE_MASK)>=leftScore?maxDiagUp:leftScore; // Compare with left
	}

	private static final float postprocess(long[] bottom, int qLen, int rLen, int[] posVector) {

		if(debug) { // If debug mode is on
			System.out.println("\n--- POSTPROCESS ---"); // Print section header
			System.out.println("qLen="+qLen+", rLen="+rLen); // Print sequence lengths
			System.out.println("bottom="+Arrays.toString(bottom)); // Print scores
			System.out.print("shifted: [");
			for(int i=0; i<bottom.length; i++) {
				if(i>0) {System.out.print(", ");}
				System.out.print(bottom[i]>>SCORE_SHIFT);
			}
			System.out.println("]");
		}

		long maxScore=Long.MIN_VALUE; // Initialize max score to minimum possible value
		int maxPos=0; // Initialize position of max score
		for(int j=1; j<=rLen; j++) { // Check each position in bottom row
			long score=bottom[j]; // Get score at this position
			if(debug) { // If debug mode is on
				System.out.println("j="+j+", score="+score+ 
						" (Score="+(score>>SCORE_SHIFT)+ 
						", Del="+((score&DEL_MASK)>>POSITION_BITS)+ 
						", Pos="+(score&POSITION_MASK)+")"); // Print score details
			}

			if(score>maxScore) { // If this score is better than current max
				maxScore=score; // Update max score
				maxPos=j; // Update position of max score
				if(debug) { // If debug mode is on
					System.out.println("New max at j="+j); // Print new max position
				}
			}
		}

		if(maxPos<qLen && rLen>=qLen && (maxScore>>SCORE_SHIFT)<0) { // Special case for mismatches
			long fullScore=bottom[qLen]; // Get score at query length position
			if(fullScore!=Long.MIN_VALUE) { // If it's a valid score
				maxScore=fullScore; // Use this score instead
				maxPos=qLen; // Update max position
				if(debug) { // If debug mode is on
					System.out.println("Using full query alignment at j="+maxPos); // Print info
				}
			}
		}

		final int originPos=(int)(maxScore&POSITION_MASK); // Get starting position
		final int endPos=maxPos; // Get ending position
		if(posVector!=null) { // If position vector was provided
			posVector[0]=originPos; // Store starting position
			posVector[1]=endPos-1; // Store ending position
		}

		final int deletions=(int)((maxScore&DEL_MASK)>>POSITION_BITS); // Get number of deletions
		final int refAlnLength=(endPos-originPos); // Calculate reference alignment length
		final long rawScore=maxScore>>SCORE_SHIFT; // Extract raw score

		final int insertions=Math.max(0, qLen+deletions-refAlnLength); // Calculate insertions
		final float matches=((rawScore+qLen+deletions)/2f); // Calculate matches
		final float substitutions=Math.max(0, qLen-matches-insertions); // Calculate substitutions

		final float identity=matches/(matches+substitutions+insertions+deletions); // Calculate identity

		if(PRINT_OPS || debug) {
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

		return identity; // Return identity score
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

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static final boolean debug=false;

}
