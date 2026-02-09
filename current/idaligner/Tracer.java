package idaligner;

import dna.AminoAcid;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Helper class for IDAligner.
 * Performs traceback and postprocessing using standard 21-bit packed fields:
 * score, deletions, position.  Assumes one contiguous score block per row.
 *
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 13, 2025
 */
public class Tracer{
	
	public static byte[] traceback(LongList trace, byte[] query, byte[] ref, int finalRow, int finalCol, ByteBuilder bb){
		if(bb==null) {bb=new ByteBuilder((query.length*5)/4);}
		else {bb.clear();}
		int r=finalRow, c=finalCol;
		
		// Find the last header(for finalRow)
		int currHeaderIdx=trace.size-1;
		// Scan back. Valid headers are negative. Negative scores start with 1111(due to 2s complement), 
		// headers start with 1000(Row index < 2M). So headers are > scores if treated as unsigned, 
		// but as signed longs, headers(10...) > scores(11...). 
		// Actually, simpler: verify the Row Index matches 'r'.
		while(true){
			long val=trace.get(currHeaderIdx);
			if(val<0){
				int row=(int)((val>>>42)&POSITION_MASK);
				if(row==r){break;}
			}
			currHeaderIdx--;
		}

		// Cache pointers to current and previous row headers
		int prevHeaderIdx=currHeaderIdx-((int)(trace.get(currHeaderIdx)&POSITION_MASK));
		
		while(r>0 && c>0){
			final long currVal=getTraceScore(trace, currHeaderIdx, c);
			
			final byte q=query[r-1];
			final byte refBase=ref[c-1];
			final boolean match=(q==refBase && AminoAcid.isFullyDefined(q));
			final boolean hasN=!(AminoAcid.isFullyDefined(q) && AminoAcid.isFullyDefined(refBase));
			final long scoreAdd=match ? MATCH : (hasN ? N_SCORE : SUB);
			
			// To find neighbors, we need:
			// Left(r, c-1): Same header(currHeaderIdx)
			// Up(r-1, c): Previous header(prevHeaderIdx)
			// Diag(r-1, c-1): Previous header(prevHeaderIdx)
			
			final long leftVal=getTraceScore(trace, currHeaderIdx, c-1);
			final long upVal=getTraceScore(trace, prevHeaderIdx, c);
			final long diagVal=getTraceScore(trace, prevHeaderIdx, c-1);
			
			final long fromLeft=leftVal+DEL_INCREMENT;
			final long fromUp=upVal+INS;
			final long fromDiag=diagVal+scoreAdd;
			final long maxDiagUp=Math.max(fromDiag, fromUp);
			
			if(currVal==maxDiagUp &&(maxDiagUp&SCORE_MASK)>=(fromLeft&SCORE_MASK)){
				if(fromDiag>=fromUp){
					bb.append(match ? 'm' : hasN ? 'N' : 'S');
					r--; c--;
				}else{
					bb.append('I');
					r--;
				}
			}else if(currVal==fromLeft){
				bb.append('D');
				c--;
			}else{
				if(r>c){bb.append('I'); r--;}
				else{bb.append('D'); c--;}
			}
			
			// If we moved up a row, shift the header pointers
			if(r < ((trace.get(currHeaderIdx)>>>42)&POSITION_MASK)){
				currHeaderIdx=prevHeaderIdx;
				// Decode the 'dist' from the new current header to find the new previous
				int dist=(int)(trace.get(currHeaderIdx)&POSITION_MASK);
				prevHeaderIdx=currHeaderIdx-dist;
			}
		}
		
		while(r>0){bb.append('I'); r--;}
		while(c>0){bb.append('D'); c--;}
		return bb.reverse().toBytes();
	}
	
	/**
	 * Traceback without looking at Query or Ref. 
	 * Deduces alignment state purely from score derivatives.
	 */
	public static byte[] traceback(LongList trace, int finalRow, int finalCol, ByteBuilder bb){
		int r=finalRow, c=finalCol;

		// 1. Find the header for the final row (You already have this loop)
		int currHeaderIdx=trace.size - 1;
		while(true){
			long val=trace.get(currHeaderIdx);
			if(val<0){
				int row=(int)((val>>>42)&POSITION_MASK);
				if(row==r){ break; }
			}
			currHeaderIdx--;
		}
		int prevHeaderIdx=currHeaderIdx-((int)(trace.get(currHeaderIdx)&POSITION_MASK));

		// 2. Get the Start Column of THIS row from the header
		long header=trace.get(currHeaderIdx);
		int rowStartCol=(int)((header>>>21)&POSITION_MASK);

		// 3. Calculate offset to the final column
		int offset=finalCol-rowStartCol;

		// 4. Retrieve the actual final score
		long finalScoreVal=trace.get(currHeaderIdx+1+offset);

		// 5. Extract the propagated rStart!
		final int rStart=(int)(finalScoreVal & POSITION_MASK);

		// 6. Size the builder
		// Ref Length + Query Length is a safe upper bound for the Cigar string
		int refLen=finalCol-rStart+1;
		int queryLen=finalRow+1; 

		if(bb==null) {
			bb=new ByteBuilder(refLen+queryLen+8);
		}else{
			bb.clear();
			bb.ensureExtra(refLen+queryLen+8);
		}

		// 2. The Loop
		while(r>0 && c>0){
			final long currVal=getTraceScore(trace, currHeaderIdx, c);
			final long diagVal=getTraceScore(trace, prevHeaderIdx, c-1);
			final long upVal=getTraceScore(trace, prevHeaderIdx, c);
			final long leftVal=getTraceScore(trace, currHeaderIdx, c-1);
			
			final long fromLeft=leftVal+DEL_INCREMENT; // Unmasked check
			final long fromUp=upVal+INS;
			
			// Calculate the derivative to check if Diag is valid
			final long delta=currVal-diagVal;
			final boolean isDiagValid=(delta==MATCH || delta==SUB || delta==N_SCORE);

			// Replicate Priority: (Diag/Up) >= Left
			if((currVal&SCORE_MASK)>=(fromLeft&SCORE_MASK)){

				// Tie-break: Diag >= Up
				// If Diag is valid and score is sufficient, take it.
				// Note: If we went Diag, currVal IS fromDiag, so currVal >= fromUp is the check.
				if(isDiagValid && currVal>=fromUp){
					if(delta==MATCH){bb.append('m');}
					else if(delta==N_SCORE){bb.append('N');}
					else{bb.append('S');} // delta == SUB
					r--; c--;
				}else{
					// Must have been Up
					bb.append('I');
					r--;
				}
			}else{
				// Must have been Left
				bb.append('D');
				c--;
			}
			
			// Move headers if we went up a row
			if(r < ((trace.get(currHeaderIdx)>>>42)&POSITION_MASK)){
				currHeaderIdx=prevHeaderIdx;
				if(currHeaderIdx < 0) break;
				int dist=(int)(trace.get(currHeaderIdx)&POSITION_MASK);
				prevHeaderIdx=currHeaderIdx-dist;
			}
		}
		
		while(r>0){bb.append('I'); r--;}
		while(c>0){bb.append('D'); c--;}
		
		return bb.reverse().toBytes();
	}
	
	private static long getTraceScore(LongList trace, int headerIdx, int c){
		long header=trace.get(headerIdx);
		int startCol=(int)((header>>>21)&POSITION_MASK);
		int relativeCol=c-startCol;
		if(relativeCol<0){return BAD;}
		
		// Bounds check: don't read past the block.
		// Next header is at headerIdx + scoreCount + 1.
		// Wait, we don't know the scoreCount easily without scanning or looking at the NEXT header.
		// But in traceback, we only query valid coordinates. 
		// If 'c' is outside the band, relativeCol will be large.
		// Let's rely on relativeCol being consistent with what was stored.
		// For safety, checking if we hit the next header(which is negative) is good.
		
		int idx=headerIdx+1+relativeCol;
		if(idx>=trace.size){return BAD;}
		long val=trace.get(idx);
		if(val<0 && (val&0x4000000000000000L)==0){
			return BAD;// Hit a header (Header starts 10..., Neg Score 11...)
		}
		return val;
	}
	
	/**
	 * Processes alignment results to calculate identity and extract alignment coordinates.
	 * Solves system of equations to determine matches, substitutions, insertions, and deletions.
	 *
	 * @param maxScore Highest score in last row of alignment matrix
	 * @param maxPos Highest-scoring position in last row
	 * @param qLen Query sequence length
	 * @param rLen Reference sequence length
	 * @param posVector Optional array for returning reference start/stop coordinates
	 * @return Identity percentage as float
	 */
	static float postprocess(long maxScore, int maxPos, int qLen, int rLen, 
			int[] posVector, AlignmentStats stats){
		// For conversion to global alignments
		if(GLOBAL && maxPos<rLen){
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
		
		// Solve the system of equations:
		// 1. M+S+I=qLen
		// 2. M+S+D=refAlnLength
		// 3. Score=M-S-I-D
		
		// Calculate operation counts
		final int insertions=Math.max(0, qLen+deletions-refAlnLength);
		final float matches=((rawScore+qLen+deletions)/2f);
		final float substitutions=Math.max(0, qLen-matches-insertions);
		final float identity=matches/(matches+substitutions+insertions+deletions);
		
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
			if(posVector.length>2){posVector[2]=rawScore;}
			if(posVector.length>3){posVector[3]=deletions;}
			if(posVector.length>4){posVector[4]=insertions;}
			if(posVector.length>5){posVector[5]=(int)matches;}
			if(posVector.length>6){posVector[6]=(int)substitutions;}
		}
		
		if(stats!=null){
			stats.clear();
			
			stats.qLen=qLen;
			stats.rLen=rLen;
			stats.identity=identity;
			stats.matches=(int)matches;
			stats.subs=(int)substitutions;
			stats.ins=insertions;
			stats.dels=deletions;
//			stats.ns=;
			stats.score=rawScore;
			stats.rStart=originPos;
			stats.rStop=endPos-1;
			
			if(PRINT_OPS) {stats.printOps();}
		}else if(PRINT_OPS){
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
	 * Inverts operations when query/ref were swapped.
	 * I -> D, D -> I, m -> m, S -> S, N -> N
	 */
	static void invertMatchString(byte[] ops){
		for(int i=0; i<ops.length; i++){
			if(ops[i]=='I'){ops[i]='D';}
			else if(ops[i]=='D'){ops[i]='I';}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean GLOBAL=false;
	public static boolean PRINT_OPS=false;
	
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
	
}
