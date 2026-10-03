package prot;

import java.util.Arrays;

import shared.Tools;
import structures.LongList;

/**
 * Fast banded protein aligner returning identity only, no traceback -- a protein analogue of
 * {@link idaligner.ScrabbleAligner}'s score-only path, for consensus-graph member screening
 * where a full O(m*n) {@link AAAligner} pass is too slow to run on every candidate.
 *
 * <p>Semantics, verified against the ported fill (not assumed from the class doc it was copied
 * from): this is GLOCAL, not local Smith-Waterman -- row 0 initializes every reference column to
 * score 0 (free reference start) and the result is read from the best cell of the LAST query row
 * only (free reference end), so the entire query is always consumed. This matches
 * {@link AAAligner#alignGlocal}, not {@link AAAligner#align}, and is the correct comparison
 * baseline: a full-query identity is what "is this candidate close enough to be a family member"
 * actually needs, not a local HSP that could report 100% on a small partial match.
 *
 * <p>Phase 1 (this class, per Brian/Elly 2026-09-03): plain match/mismatch scoring on
 * Blosum62-encoded residues (0-19 or {@link Blosum62#X_CODE}), linear indel cost, banded DP,
 * two reused packed-{@code long} rows -- structurally the same algorithm as
 * {@link idaligner.ScrabbleAligner#alignStatic}, with {@code X_CODE} playing the role
 * nucleotide 'N' plays there (never a match, contributes zero). No traceback, no
 * {@link AAAlignment}, no BLOSUM62 substitution scores -- those are deferred to a later
 * BLOSUM-aware variant (which would report a score/max-score ratio, not identity, per Brian)
 * and a separate traced glocal mode for graph insertion (per Elly, matching
 * {@link AAAligner#alignGlocal}'s semantics), once this phase is validated.
 *
 * @author Eru
 */
public final class ScrabbleAmino {

	/** Utility class; no instances. */
	private ScrabbleAmino(){}

	/**
	 * Banded glocal alignment (full query consumed, free reference start/end), score/identity
	 * only. Query and reference are Blosum62-encoded residue arrays (0-19 standard,
	 * {@link Blosum62#X_CODE} for ambiguous). Unlike {@link idaligner.ScrabbleAligner}, query and
	 * reference are NEVER swapped: glocal is asymmetric (query is the one fully consumed), and
	 * this must match {@link AAAligner#alignGlocal}'s (member, pivot) calling convention, where
	 * the member can legitimately be longer than the not-yet-fully-grown pivot.
	 *
	 * @param query Encoded query (member) residues -- fully consumed.
	 * @param ref Encoded reference (pivot) residues -- free start and end.
	 * @return Identity (0.0-1.0).
	 */
	public static final float scoreOnly(byte[] query, byte[] ref){
		return scoreOnly(query, ref, 0);
	}

	/**
	 * As {@link #scoreOnly(byte[],byte[])}, but with an explicit minimum band half-width floor
	 * ({@code minBand}) that {@link #decideBandwidth}'s adaptive estimate can widen past but never
	 * narrow below, and that {@code dynamicBW} never shrinks under either -- per Brian's direction
	 * (2026-09-03): real family divergence (mean ~55% identity, some members 20-30%) is far
	 * outside the >90%-identity/low-indel envelope the adaptive heuristic assumes, so a sweepable
	 * hard floor is needed to test "how wide is wide enough" against real data, not a heuristic
	 * that narrows itself back down based on short local match runs. minBand=0 reproduces the
	 * original adaptive-only behavior exactly.
	 */
	public static final float scoreOnly(byte[] query, byte[] ref, int minBand){
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length, rLen=ref.length;

		final int bandWidth0=Math.max(decideBandwidth(query, ref), minBand);
		final int maxDrift=2, maxDynamic=Math.max((bandWidth0*12)/4, minBand);

		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);
		for(int j=0; j<=rLen; j++){prev[j]=j;}//free reference start (local): row 0 = 0 everywhere

		int bandStart=1, bandEnd=rLen-1;
		int center=0;
		int prevBandEnd=rLen;
		int maxPos=0;
		long maxScore=2*SUB;
		int dynamicBW=0, deltaBW=0;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			final boolean nextMatch=(q==ref[Math.min(rLen-1, maxPos)]);
			if(nextMatch){deltaBW=(deltaBW<0 ? Math.max(-maxDynamic, deltaBW*2) : -2);}
			else{deltaBW=Tools.mid(1, (maxDynamic-dynamicBW)/2, 8);}
			dynamicBW=Tools.mid(0, dynamicBW+deltaBW, maxDynamic);

			final int bandWidth=bandWidth0+Math.max(16+bandWidth0*12-maxDrift*i, dynamicBW);
			final int quarterBand=bandWidth/4;
			final int drift=Tools.mid(-1, maxPos-center, maxDrift);
			center=center+1+drift;
			bandStart=Math.max(bandStart, center-bandWidth+quarterBand);
			bandEnd=Math.min(rLen, center+bandWidth+quarterBand);

			//Correctness: prev[] cells beyond the previous row's band still hold values from two
			//rows ago (double buffering). If the band re-widened, clear exactly those stale cells
			//before the inner loop reads them (see ScrabbleAligner.alignStatic, same reasoning).
			if(bandEnd>prevBandEnd){Arrays.fill(prev, prevBandEnd+1, bandEnd+1, BAD);}

			curr[bandStart-1]=BAD;
			curr[0]=i*INDEL;
			maxScore=BAD;
			maxPos=-1;

			for(int j=bandStart; j<=bandEnd; j++){
				final byte r=ref[j-1];

				final boolean isMatch=(q==r && q!=Blosum62.X_CODE);
				final boolean hasX=(q==Blosum62.X_CODE || r==Blosum62.X_CODE);
				final long scoreAdd=isMatch ? MATCH : (hasX ? X_SCORE : SUB);

				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long diagScore=pj1+scoreAdd;
				final long upScore=pj+INDEL;
				final long leftScore=cj1+DEL_INCREMENT;

				final long maxDiagUp=Math.max(diagScore, upScore);
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore ? maxDiagUp : leftScore;

				curr[j]=maxValue;

				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}

			long[] tmp=prev; prev=curr; curr=tmp;
			prevBandEnd=bandEnd;
		}
		return postprocessIdentity(maxScore, maxPos, qLen);
	}

	/**
	 * Banded glocal alignment WITH traceback, for consensus-graph member insertion --
	 * {@link AAGraph#add(byte[], AAAlignment)}'s exact input contract. Same fill as
	 * {@link #scoreOnly}, but records a sparse packed trace (row headers + per-row score cells,
	 * same layout as {@link idaligner.ScrabbleAligner#alignAndTraceStatic}) and reconstructs the
	 * op path using the ACTUAL query/ref bytes -- never a byte-less packed-delta reconstruction,
	 * which {@link idaligner.Tracer} documents (Aug 21 2026) as non-injective and a real
	 * production bug once shipped. qStart is always 0 and every query residue appears as an 'm'
	 * or 'I' op, matching {@link AAAligner#alignGlocal}'s own documented contract exactly (the
	 * entire query is placed; only 'm'/'D'/'I' appear, never AAAligner's local-only 'S'/'N'
	 * distinction -- identity is tallied separately, matching AAAlignment's field model).
	 *
	 * @param query Encoded query (member) residues -- fully consumed.
	 * @param ref Encoded reference (pivot) residues -- free start and end.
	 * @return An {@link AAAlignment} usable directly by {@link AAGraph#add}.
	 */
	public static final AAAlignment alignWithTrace(byte[] query, byte[] ref){
		return alignWithTrace(query, ref, 0);
	}

	/**
	 * As {@link #alignWithTrace(byte[],byte[])}, but with the same explicit minimum band
	 * half-width floor {@link #scoreOnly(byte[],byte[],int)} takes -- see that overload's doc for
	 * why (real family divergence far outside the adaptive heuristic's >90%-identity envelope).
	 * {@code minBand=0} reproduces the original adaptive-only behavior exactly.
	 */
	public static final AAAlignment alignWithTrace(byte[] query, byte[] ref, int minBand){
		if(query.length==0 || ref.length==0){return null;}//matches AAAligner.alignGlocal's contract
		assert(ref.length<=POSITION_MASK) : "Ref is too long: "+ref.length+">"+POSITION_MASK;
		final int qLen=query.length, rLen=ref.length;

		final LongList trace=new LongList(qLen*4+16);
		int lastHeaderIdx=0;

		final int bandWidth0=Math.max(decideBandwidth(query, ref), minBand);
		final int maxDrift=2, maxDynamic=Math.max((bandWidth0*12)/4, minBand);

		long[] prev=new long[rLen+1], curr=new long[rLen+1];
		Arrays.fill(curr, BAD);
		for(int j=0; j<=rLen; j++){prev[j]=j;}
		{//Row 0 header + cells, so the traceback can read column 0's synthesized boundary too.
			trace.add(0x8000000000000000L);
			lastHeaderIdx=0;
			for(int j=0; j<=rLen; j++){trace.add(prev[j]);}
		}

		int bandStart=1, bandEnd=rLen-1;
		int center=0;
		int prevBandEnd=rLen;
		int maxPos=0;
		long maxScore=2*SUB;
		int dynamicBW=0, deltaBW=0;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			final boolean nextMatch=(q==ref[Math.min(rLen-1, maxPos)]);
			if(nextMatch){deltaBW=(deltaBW<0 ? Math.max(-maxDynamic, deltaBW*2) : -2);}
			else{deltaBW=Tools.mid(1, (maxDynamic-dynamicBW)/2, 8);}
			dynamicBW=Tools.mid(0, dynamicBW+deltaBW, maxDynamic);

			final int bandWidth=bandWidth0+Math.max(16+bandWidth0*12-maxDrift*i, dynamicBW);
			final int quarterBand=bandWidth/4;
			final int drift=Tools.mid(-1, maxPos-center, maxDrift);
			center=center+1+drift;
			bandStart=Math.max(bandStart, center-bandWidth+quarterBand);
			bandEnd=Math.min(rLen, center+bandWidth+quarterBand);

			if(bandEnd>prevBandEnd){Arrays.fill(prev, prevBandEnd+1, bandEnd+1, BAD);}

			if(true){
				final int dist=trace.size-lastHeaderIdx;
				assert(dist<=(int)POSITION_MASK);
				final long header=0x8000000000000000L|((long)i<<42)|((long)bandStart<<21)|dist;
				lastHeaderIdx=trace.size;
				trace.add(header);
			}

			curr[bandStart-1]=BAD;
			curr[0]=i*INDEL;
			maxScore=BAD;
			maxPos=-1;

			for(int j=bandStart; j<=bandEnd; j++){
				final byte r=ref[j-1];

				final boolean isMatch=(q==r && q!=Blosum62.X_CODE);
				final boolean hasX=(q==Blosum62.X_CODE || r==Blosum62.X_CODE);
				final long scoreAdd=isMatch ? MATCH : (hasX ? X_SCORE : SUB);

				final long pj1=prev[j-1], pj=prev[j], cj1=curr[j-1];
				final long maxDiagUp=Math.max(pj1+scoreAdd, pj+INDEL);
				final long maxValue=(maxDiagUp&SCORE_MASK)>=(cj1+DEL_INCREMENT) ? maxDiagUp : (cj1+DEL_INCREMENT);

				curr[j]=maxValue;

				final boolean better=((maxValue&SCORE_MASK)>maxScore);
				maxScore=better ? maxValue : maxScore;
				maxPos=better ? j : maxPos;
			}
			trace.add(curr, bandStart, bandEnd+1);

			long[] tmp=prev; prev=curr; curr=tmp;
			prevBandEnd=bandEnd;
		}

		return traceback(trace, query, ref, qLen, maxPos, maxScore);
	}

	/**
	 * Sequence-aware traceback from the best last-row cell back to row 0, producing an
	 * AAGraph-compatible {@link AAAlignment}. Adapted from {@link idaligner.Tracer#traceback}'s
	 * CURRENT (non-disabled) method -- same header/cell lookup scheme -- but emits AAAligner's
	 * op vocabulary ('m' for every diagonal move, identity tallied separately) instead of
	 * Tracer's m/S/N split, and treats {@link Blosum62#X_CODE} as the ambiguous residue instead
	 * of nucleotide 'N'.
	 */
	private static AAAlignment traceback(LongList trace, byte[] query, byte[] ref,
			int finalRow, int finalCol, long maxScore){
		int r=finalRow, c=finalCol;
		final int qStop=finalRow-1, tStop=finalCol-1;
		int tStart=finalCol-1;

		int identities=0, mismatches=0, gapOpens=0, length=0;
		boolean prevQGap=false, prevTGap=false;

		final byte[] ops=new byte[query.length+ref.length];
		int op=ops.length;

		int currHeaderIdx=trace.size-1;
		while(true){
			final long val=trace.get(currHeaderIdx);
			if(val<0){
				final int row=(int)((val>>>42)&POSITION_MASK);
				if(row==r){break;}
			}
			currHeaderIdx--;
		}
		int prevHeaderIdx=currHeaderIdx-((int)(trace.get(currHeaderIdx)&POSITION_MASK));
		int currBlockEnd=trace.size;

		while(r>0 && c>0){
			final long currVal=getTraceScore(trace, currHeaderIdx, currBlockEnd, c);
			final byte q=query[r-1], t=ref[c-1];
			final boolean match=(q==t && q!=Blosum62.X_CODE);
			final boolean hasX=(q==Blosum62.X_CODE || t==Blosum62.X_CODE);
			final long scoreAdd=match ? MATCH : (hasX ? X_SCORE : SUB);

			final long leftVal=getTraceScore(trace, currHeaderIdx, currBlockEnd, c-1);
			final long upVal=getTraceScore(trace, prevHeaderIdx, currHeaderIdx, c);
			final long diagVal=getTraceScore(trace, prevHeaderIdx, currHeaderIdx, c-1);

			final long fromLeft=leftVal+DEL_INCREMENT;
			final long fromUp=upVal+INDEL;
			final long fromDiag=diagVal+scoreAdd;
			final long maxDiagUp=Math.max(fromDiag, fromUp);

			if((maxDiagUp&SCORE_MASK)>=fromLeft && currVal==maxDiagUp){
				if(fromDiag>=fromUp){
					length++;
					if(match){identities++;}else{mismatches++;}
					prevQGap=false; prevTGap=false;
					tStart=c-1;
					ops[--op]='m';
					r--; c--;
				}else{
					length++;
					if(!prevQGap){gapOpens++;}
					prevQGap=true; prevTGap=false;
					ops[--op]='I';
					r--;
				}
			}else if(currVal==fromLeft){
				length++;
				if(!prevTGap){gapOpens++;}
				prevTGap=true; prevQGap=false;
				ops[--op]='D';
				c--;
			}else{
				//No predecessor reproduces the stored value -- the recorded band did not connect
				//this cell (band re-narrowed past it, not a genuine desync since maxDynamic can
				//discard cells a wider pass would have kept). Fall back to forced insertions,
				//matching Tracer.traceback's own end-of-loop fallback for exhausted columns.
				break;
			}

			if(r<((trace.get(currHeaderIdx)>>>42)&POSITION_MASK)){
				currBlockEnd=currHeaderIdx;
				currHeaderIdx=prevHeaderIdx;
				final int dist=(int)(trace.get(currHeaderIdx)&POSITION_MASK);
				prevHeaderIdx=currHeaderIdx-dist;
			}
		}
		//Glocal contract: the ENTIRE query is placed, starting at residue 0 (AAAligner.
		//glocalTraceback's own documented convention). Any remaining rows -- whether the walk
		//reached c==0 (free reference start, genuinely done) or the recorded band ran out first
		//-- become insertions, exactly like Tracer.traceback's identical fallback.
		while(r>0){
			length++;
			if(!prevQGap){gapOpens++;}
			prevQGap=true; prevTGap=false;
			ops[--op]='I';
			r--;
		}

		final byte[] match=Arrays.copyOfRange(ops, op, ops.length);
		final int rawScore=(int)(maxScore>>SCORE_SHIFT);
		return new AAAlignment(rawScore, 0, qStop, tStart, tStop,
			identities, mismatches, gapOpens, length, match);
	}

	/** Reads one row's packed score cell for column c, or a boundary/BAD sentinel. Adapted from
	 *  {@link idaligner.Tracer#getTraceScore} (package-private there; duplicated for the same
	 *  reason as {@link #postprocessIdentity}). */
	private static long getTraceScore(LongList trace, int headerIdx, int blockEnd, int c){
		final long header=trace.get(headerIdx);
		final int startCol=(int)((header>>>21)&POSITION_MASK);
		final int relativeCol=c-startCol;
		if(relativeCol<0){
			if(c==0){
				final int row=(int)((header>>>42)&POSITION_MASK);
				return row*INDEL;
			}
			return BAD;
		}
		final int idx=headerIdx+1+relativeCol;
		if(idx>=blockEnd || idx>=trace.size){return BAD;}
		final long val=trace.get(idx);
		if(val<0 && (val&0x4000000000000000L)==0){return BAD;}//hit a header, not a score cell
		return val;
	}

	/**
	 * Adaptive bandwidth from sequence lengths and an early-mismatch scan, identical in shape
	 * to {@link idaligner.ScrabbleAligner}'s heuristic (high-identity, low-indel assumption).
	 * CORRECTED 2026-09-03: real mmseqs protein-family members are NOT already 90/80/80
	 * identity-grouped (that claim conflated a different, unrelated identity-grouping sub-project
	 * with these families) -- a real-family canary measured mean identity vs the medoid at 0.5466,
	 * spanning 20-100%, which this heuristic alone cannot band correctly. See
	 * {@code records/SCRABBLEAMINO_CANARY_v1.md} and {@link #scoreOnly(byte[],byte[],int)}'s
	 * {@code minBand} floor, added specifically to compensate.
	 */
	private static int decideBandwidth(byte[] query, byte[] ref){
		int subs=0, qLen=query.length, rLen=ref.length;
		int bandwidth=Tools.mid(7, 1+Math.max(qLen, rLen)/32, 20+(int)Math.sqrt(rLen)/8);
		for(int i=0, minlen=Math.min(qLen, rLen); i<minlen && subs<bandwidth; i++){
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);
	}

	/**
	 * Unpacks the best cell's score/deletions and solves the same M+S+I=qLen / M+S+D=refLen /
	 * score=M-S-I-D system {@link idaligner.Tracer#postprocess} uses -- duplicated here (not
	 * called cross-package; that method is package-private in idaligner) rather than widening
	 * Tracer's visibility for a first cut. Revisit consolidation once this design is validated.
	 */
	private static float postprocessIdentity(long maxScore, int maxPos, int qLen){
		final int originPos=(int)(maxScore&POSITION_MASK);
		final int endPos=maxPos;
		final int deletions=(int)((maxScore&DEL_MASK)>>POSITION_BITS);
		final int refAlnLength=(endPos-originPos);
		final int rawScore=(int)(maxScore>>SCORE_SHIFT);

		final int insertions=Math.max(0, qLen+deletions-refAlnLength);
		final float matches=((rawScore+qLen+deletions)/2f);
		final float substitutions=Math.max(0, qLen-matches-insertions);
		return matches/(matches+substitutions+insertions+deletions);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	//Bit layout identical to idaligner.ScrabbleAligner/Tracer -- same packed-long scheme,
	//same field widths; protein sequences here are far shorter than the 2M-residue cap.
	private static final int POSITION_BITS=21;
	private static final int DEL_BITS=21;
	private static final int SCORE_SHIFT=POSITION_BITS+DEL_BITS;

	private static final long POSITION_MASK=(1L<<POSITION_BITS)-1;
	private static final long DEL_MASK=((1L<<DEL_BITS)-1)<<POSITION_BITS;
	private static final long SCORE_MASK=~(POSITION_MASK|DEL_MASK);

	//Linear scoring, uniform +-1 scale -- Elly's spec: no BLOSUM62, no separate gap-open
	//premium for this first cut. INDEL folds the (identical) insertion and deletion cost.
	private static final long MATCH=1L<<SCORE_SHIFT;
	private static final long SUB=(-1L)<<SCORE_SHIFT;
	private static final long INDEL=(-1L)<<SCORE_SHIFT;
	private static final long X_SCORE=0L;
	private static final long BAD=Long.MIN_VALUE/2;
	private static final long DEL_INCREMENT=INDEL+(1L<<POSITION_BITS);
}
