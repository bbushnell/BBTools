package idaligner;

import structures.ByteBuilder;
import structures.LongList;

/**
 * Affine 3-state traceback for ScrabbleAffine.
 * Walks backward through interleaved M/I/D trace data stored by
 * ScrabbleAffine.alignWithTrace (Mode 2) to produce a match string.
 *
 * Storage: one LongList with interleaved rows. Each row:
 * [header] [M0 I0 D0] [M1 I1 D1] ... [Mn In Dn]
 * Header (negative long): bit63=1, bits61-42=row, bits41-21=startCol, bits20-0=distToPrev.
 *
 * @author Ady
 * @date July 1, 2026
 */
public class TracerAffine{

	static final int STATE_M=0, STATE_I=1, STATE_D=2;

	/*--------------------------------------------------------------*/
	/*----------------       Header packing        ----------------*/
	/*--------------------------------------------------------------*/

	private static final int HEADER_POS_BITS=21;
	private static final long HEADER_POS_MASK=(1L<<HEADER_POS_BITS)-1;

	static long packHeader(int row, int startCol, int distToPrev){
		return (1L<<63)
			|((long)row<<42)
			|((long)startCol<<21)
			|(distToPrev&HEADER_POS_MASK);
	}

	private static int headerRow(long h){return (int)((h>>>42)&HEADER_POS_MASK);}
	private static int headerStartCol(long h){return (int)((h>>>21)&HEADER_POS_MASK);}
	private static int headerDist(long h){return (int)(h&HEADER_POS_MASK);}

	/** Headers have bit63=1, bit62=0 (pattern 10...). Negative scores have 11... */
	private static boolean isHeader(long v){
		return v<0 && (v&0x4000000000000000L)==0;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Trace access          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retrieve a cell value from the interleaved trace.
	 * @param state 0=M, 1=I, 2=D
	 */
	private static long getTraceValue(LongList trace, int headerIdx, int col, int state){
		if(headerIdx<0) return ScrabbleAffine.BAD;
		long header=trace.get(headerIdx);
		int startCol=headerStartCol(header);
		int relativeCol=col-startCol;
		if(relativeCol<0) return ScrabbleAffine.BAD;
		int idx=headerIdx+1+3*relativeCol+state;
		if(idx>=trace.size) return ScrabbleAffine.BAD;
		long val=trace.get(idx);
		if(isHeader(val)) return ScrabbleAffine.BAD;
		return val;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Traceback           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Produce a match string by walking backward through stored M/I/D data.
	 *
	 * @param trace     Interleaved trace from ScrabbleAffine.alignWithTrace
	 * @param query     Query sequence
	 * @param ref       Reference sequence
	 * @param finalRow  Last query row (1-based, typically query.length)
	 * @param finalCol  Best ref column (1-based)
	 * @param bb        Reusable ByteBuilder (allocated if null)
	 * @return match string bytes (m=match, S=sub, N=nocall, I=insertion, D=deletion)
	 */
	public static byte[] traceback(LongList trace, byte[] query, byte[] ref,
			int finalRow, int finalCol, ByteBuilder bb){
		if(bb==null){bb=new ByteBuilder(query.length+ref.length);}
		else{bb.clear();}

		int r=finalRow, c=finalCol;

		// Find header for the final row
		int currHeaderIdx=findHeader(trace, r);
		int prevHeaderIdx=currHeaderIdx-headerDist(trace.get(currHeaderIdx));

		// Determine starting state: best of M and I at the endpoint.
		// Glocal: ending in D means trailing deletions, not counted.
		long mVal=getTraceValue(trace, currHeaderIdx, c, STATE_M);
		long iVal=getTraceValue(trace, currHeaderIdx, c, STATE_I);
		int state=(ScrabbleAffine.scoreOf(mVal)>=ScrabbleAffine.scoreOf(iVal))
			? STATE_M : STATE_I;

		while(r>0 && c>0){
			boolean rowChanged=false;

			if(state==STATE_M){
				// Diagonal move from (r-1, c-1). Character depends on bases.
				final byte q=query[r-1], rb=ref[c-1];
				final boolean isMatch=(q==rb && q!='N' && rb!='N');
				final boolean hasN=(q=='N' || rb=='N');
				bb.append(isMatch ? (byte)'m' : hasN ? (byte)'N' : (byte)'S');

				// Predecessor at (r-1, c-1): pick highest score. M>I>D on ties.
				long pM=getTraceValue(trace, prevHeaderIdx, c-1, STATE_M);
				long pI=getTraceValue(trace, prevHeaderIdx, c-1, STATE_I);
				long pD=getTraceValue(trace, prevHeaderIdx, c-1, STATE_D);
				long best=ScrabbleAffine.scoreOf(pM); state=STATE_M;
				if(ScrabbleAffine.scoreOf(pI)>best){best=ScrabbleAffine.scoreOf(pI); state=STATE_I;}
				if(ScrabbleAffine.scoreOf(pD)>best){state=STATE_D;}

				r--; c--;
				rowChanged=true;
			}else if(state==STATE_I){
				// Up move from (r-1, c). Insertion.
				bb.append((byte)'I');
				long cell=getTraceValue(trace, currHeaderIdx, c, STATE_I);
				// time==1 → opened from M; time>1 → extended from I
				state=(ScrabbleAffine.timeOf(cell)==1) ? STATE_M : STATE_I;
				r--;
				rowChanged=true;
			}else{
				// Left move from (r, c-1). Deletion.
				bb.append((byte)'D');
				long cell=getTraceValue(trace, currHeaderIdx, c, STATE_D);
				// time==1 → opened from M; time>1 → extended from D
				state=(ScrabbleAffine.timeOf(cell)==1) ? STATE_M : STATE_D;
				c--;
			}

			if(rowChanged){
				currHeaderIdx=prevHeaderIdx;
				if(currHeaderIdx>=0){
					prevHeaderIdx=currHeaderIdx-headerDist(trace.get(currHeaderIdx));
				}
			}
		}

		// Remaining query bases are insertions (leading unmatched query)
		while(r>0){bb.append((byte)'I'); r--;}
		// No trailing D's — glocal mode, unaligned ref prefix is not part of alignment

		return bb.reverse().toBytes();
	}

	/** Scan backward from the end to find the header for a given row. */
	private static int findHeader(LongList trace, int targetRow){
		for(int i=trace.size-1; i>=0; i--){
			long val=trace.get(i);
			if(isHeader(val) && headerRow(val)==targetRow){
				return i;
			}
		}
		throw new RuntimeException("TracerAffine: header not found for row "+targetRow);
	}

}
