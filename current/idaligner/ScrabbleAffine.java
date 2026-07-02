package idaligner;

import java.util.concurrent.atomic.AtomicLong;

import structures.ByteBuilder;
import structures.LongList;

/**
 * Affine, 3-state variant of ScrabbleAligner (WORKING NAME -- rename to taste).
 *
 * Same trick as ScrabbleAligner -- pack the alignment-start coordinate into the
 * score word so rStart can be read back with no traceback -- but with three
 * separate score states (match/sub "M", insertion "I", deletion "D") so it can
 * do affine gaps. That forces 6 rolling arrays (3 states x {prev,curr}) instead
 * of 2: three independent max-tracked scores cannot share one 64-bit cell.
 *
 * Packed long per cell: score (high 22 bits) | timeInState (mid 21) | start (low 21).
 * "timeInState" replaces ScrabbleAligner's deletion-count field; it is the length
 * of the current run in this state and is fed to the cost hooks, so the scoring
 * model can be retuned (open/extend now, full tiers later) by editing only the
 * cost functions -- the recurrence never hardcodes a penalty.
 *
 * v1 scope (deliberately simplified -- see inline NOTEs):
 *  - FULL DP, not yet banded (ScrabbleAligner's adaptive band ports on top later).
 *  - Affine 2-tier costs: open=2, extend=1, via subCost/insCost/delCost(timeInState).
 *  - Returns correct score + rStart/rStop. identity() is a STUB pending a decision
 *    (parallel match/length counters, vs packing one) -- exact affine identity
 *    cannot be inverted from the score the way ScrabbleAligner does.
 *
 * @author Ady
 * @contributor Brian Bushnell
 * @date June 21, 2026
 */
public class ScrabbleAffine implements IDAligner{

	public ScrabbleAffine(){}

	/** Self-test with BBMap's exact scoring. */
	public static void main(String[] args){
		// BBMap scoring: first match=70, consecutive=100, sub=-127, insOpen=-395, delOpen=-472, delExt=-33
		selfTest("ACGT", "ACGT");      // 4bp perfect: 70+100+100+100=370
		selfTest("ACGT", "AGGT");      // 1 sub at pos 1: 70 - 127 + 100 + 100 = 143
		selfTest("ACT",  "ACGT");      // 1 del: 70+100 - 472 + 100 = -202... hmm, del is very expensive
		selfTest("CGT",  "ACGT");      // start!=0: start=1, 70+100+100 = 270
		selfTest("ACGT", "ACT");       // 1 ins: depends on path
		selfTest("AAATTT", "AAAGGTTT");// 2bp del: matches around a deletion
		// Gapped reference tests
		selfTest("AAATTT", "AAA-TTT");    // gap symbol = cheap deletion
		selfTest("AAATTT", "AAA---TTT");  // 3 gap symbols

		// Mode 2 (trace + traceback) tests
		System.out.println("\n--- Mode 2: alignWithTrace + TracerAffine ---");
		selfTestTrace("ACGT", "ACGT");        // expect mmmm
		selfTestTrace("ACGT", "AGGT");        // expect mSmm
		selfTestTrace("ACT",  "ACGT");        // 1 del: mmDm
		selfTestTrace("CGT",  "ACGT");        // start=1: mmm
		selfTestTrace("ACGT", "ACT");         // 1 ins
		selfTestTrace("AAATTT", "AAAGGTTT");  // 2bp del: mmmDDmmm
		selfTestTrace("AAATTT", "AAA-TTT");   // gap symbol
		selfTestTrace("AAATTT", "AAA---TTT"); // 3 gap symbols
	}
	private static void selfTest(String q, String r){
		int[] pos=new int[3];
		float id=alignStatic(q.getBytes(), r.getBytes(), pos);
		System.out.println(q+" vs "+r+"  ->  start="+pos[0]+" stop="+pos[1]+" score="+pos[2]+"  (id stub="+id+")");
	}
	private static void selfTestTrace(String q, String r){
		byte[] query=q.getBytes(), ref=r.getBytes();
		int[] pos=new int[3];
		LongList trace=new LongList();
		alignWithTrace(query, ref, pos, trace);
		int finalCol=pos[1]+1; // rStop is 0-based; finalCol is 1-based
		byte[] match=TracerAffine.traceback(trace, query, ref, query.length, finalCol, null);
		// Cross-check: Mode 1 and Mode 2 should agree on score
		int[] pos1=new int[3];
		alignStatic(query, ref, pos1);
		String scoreMatch=(pos[2]==pos1[2]) ? "OK" : "MISMATCH(mode1="+pos1[2]+")";
		System.out.println(q+" vs "+r+"  ->  start="+pos[0]+" stop="+pos[1]+
			" score="+pos[2]+" ["+scoreMatch+"]  match="+new String(match));
	}

	@Override public final String name(){return "ScrabbleAffine";}
	@Override public final float align(byte[] a, byte[] b){return alignStatic(a, b, null);}
	@Override public final float align(byte[] a, byte[] b, int[] pos){return alignStatic(a, b, pos);}
	@Override public final float align(byte[] a, byte[] b, int[] pos, int minScore){return alignStatic(a, b, pos);}

	@Override public long loops(){return loops.get();}
	@Override public void setLoops(long x){loops.set(x);}
	private static final AtomicLong loops=new AtomicLong(0);

	/*--------------------------------------------------------------*/
	/*--------  Cost hooks: the ONLY place the scoring lives  -------*/
	/*--------------------------------------------------------------*/

	// timeInState = length of the current run; 0 means "opening" this state.
	// Returns a positive penalty (subtracted from score). BBMap's exact constants.
	private static final long POINTS_SUB=127, POINTS_SUBR=147;
	private static final long POINTS_SUB2=51, POINTS_SUB3=25;

	private static long subCost(int timeInState){
		if(timeInState==0) return POINTS_SUB;
		if(timeInState==1) return POINTS_SUB2;
		return POINTS_SUB3;
	}
	private static long subCostAfterShortMatch(int timeInState){
		return timeInState==0 ? POINTS_SUBR : subCost(timeInState);
	}

	private static final long POINTS_INS=395, POINTS_INS2=39, POINTS_INS3=23, POINTS_INS4=8;
	private static long insCost(int timeInState){
		if(timeInState==0) return POINTS_INS;
		if(timeInState<5) return POINTS_INS2;
		if(timeInState<20) return POINTS_INS3;
		return POINTS_INS4;
	}

	private static final long POINTS_DEL=472, POINTS_DEL2=33, POINTS_DEL3=9;
	private static final long POINTS_DEL4=1, POINTS_DEL5=1;
	private static final int DEL_COST3_LIMIT=5, DEL_COST4_LIMIT=20, DEL_COST5_LIMIT=80;
	private static final int DEL_TIMESLIP=4;
	private static long delCost(int timeInState){
		if(timeInState==0) return POINTS_DEL;
		if(timeInState<DEL_COST3_LIMIT) return POINTS_DEL2;
		if(timeInState<DEL_COST4_LIMIT) return POINTS_DEL3;
		if(timeInState<DEL_COST5_LIMIT) return POINTS_DEL4;
		return ((timeInState & (DEL_TIMESLIP-1))==0) ? POINTS_DEL5 : 0;
	}

	private static final long POINTS_DEL_REF_N=10;

	/** Cost for traversing a gap symbol (representing GAPLEN deleted bases). Nearly free. */
	private static long gapSymbolCost(){return GAP_SYMBOL_PENALTY;}
	private static final long GAP_SYMBOL_PENALTY=2;

	/** The gap symbol byte used by GappedReference. */
	public static final byte GAP_SYMBOL='-';

	/*--------------------------------------------------------------*/
	/*----------------       Banding               ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether to use adaptive banding. Off = full DP (for verification). */
	public static boolean USE_BANDING=true;

	/** Decide initial bandwidth from sequence lengths + a quick identity scan. */
	private static int decideBandwidth(byte[] query, byte[] ref){
		final int qLen=query.length, rLen=ref.length;
		int bandwidth=shared.Tools.mid(7, 1+Math.max(qLen, rLen)/32, 20+(int)Math.sqrt(rLen)/8);
		int subs=0;
		for(int i=0, minlen=Math.min(qLen, rLen); i<minlen && subs<bandwidth; i++){
			if(query[i]!=ref[i]){subs++;}
		}
		return Math.min(subs+1, bandwidth);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Core DP             ----------------*/
	/*--------------------------------------------------------------*/

	public static final float alignStatic(byte[] query, byte[] ref, int[] posVector){
		// Swap to keep query no longer than ref (only when caller doesn't need positions)
		if(posVector==null && query.length>ref.length){
			byte[] tmp=query; query=ref; ref=tmp;
		}
		final int qLen=query.length, rLen=ref.length;
		assert(rLen<=POSITION_MASK) : "Ref too long: "+rLen+">"+POSITION_MASK;
		long mloops=0;

		// Banding parameters
		final boolean banded=USE_BANDING;
		final int bandWidth0=banded ? decideBandwidth(query, ref) : 0;
		final int maxDrift=2, maxDynamic=banded ? (bandWidth0*12)/4 : 0;
		int center=0, dynamicBW=0, deltaBW=0;

		// 6 rolling arrays: 3 states x {prev,curr}
		long[] prevM=new long[rLen+1], currM=new long[rLen+1];
		long[] prevI=new long[rLen+1], currI=new long[rLen+1];
		long[] prevD=new long[rLen+1], currD=new long[rLen+1];

		// Row 0: the alignment may BEGIN at any ref column (glocal) -> score 0, start=j.
		// Time=1 so that the first match scores MATCH_FIRST (70), not MATCH (100).
		// (time=0 in M means "in a match streak"; time>0 means "not consecutive match".)
		// I and D cannot exist before any query base is consumed.
		for(int j=0; j<=rLen; j++){
			prevM[j]=pack(0, 1, j);
			prevI[j]=BAD;
			prevD[j]=BAD;
		}

		long bestScore=BAD; int bestPos=0; long bestWord=BAD;
		int bandStart=1, bandEnd=rLen;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			// Adaptive band: widen on mismatches, narrow on matches
			if(banded){
				final boolean nextMatch=(bestPos>=0 && bestPos<rLen && q==ref[Math.min(rLen-1, bestPos)]);
				if(nextMatch){
					deltaBW=(deltaBW<0 ? Math.max(-maxDynamic, deltaBW*2) : -2);
				}else{
					deltaBW=shared.Tools.mid(1, (maxDynamic-dynamicBW)/2, 8);
				}
				dynamicBW=shared.Tools.mid(0, dynamicBW+deltaBW, maxDynamic);
				final int bandWidth=bandWidth0+Math.max(16+bandWidth0*12-maxDrift*i, dynamicBW);
				final int quarterBand=bandWidth/4;
				final int drift=shared.Tools.mid(-1, bestPos-center, maxDrift);
				center=center+1+drift;
				bandStart=Math.max(1, center-bandWidth+quarterBand);
				bandEnd=Math.min(rLen, center+bandWidth+quarterBand);
			}

			// Column 0: leading query bases with no ref consumed are insertions.
			currM[0]=BAD;
			currD[0]=BAD;
			currI[0]=pack(0 - insCost(0) - (long)(i-1)*insCost(1), i, 0); // open + (i-1) extends

			// Clear stale data at band edges
			if(banded && bandStart>1){
				currM[bandStart-1]=BAD;
				currI[bandStart-1]=BAD;
				currD[bandStart-1]=BAD;
			}

			long rowBestScore=BAD; int rowBestPos=0; long rowBestWord=BAD;

			for(int j=bandStart; j<=bandEnd; j++){
				final byte r=ref[j-1];
				final boolean isGap=(r==GAP_SYMBOL);

				if(isGap){
					// Gap symbol: force a deletion at near-zero cost.
					// M and I are dead (can't match or insert against a gap symbol).
					currM[j]=BAD;
					currI[j]=BAD;

					// D: open from M or extend from D, but with gapSymbolCost instead of delCost.
					final long openD=currM[j-1], extD=currD[j-1];
					final long dOpen=scoreOf(openD)-gapSymbolCost();
					final long dExt=scoreOf(extD)-gapSymbolCost();
					currD[j]=(dOpen>=dExt)
						? pack(dOpen, 1, startOf(openD))
						: pack(dExt, timeOf(extD)+1, startOf(extD));
				}else{
				final boolean isMatch=(q==r && r!='N' && q!='N');
				final boolean hasN=(q=='N' || r=='N');

				/*--- M: diagonal move from the best of {M,I,D}[i-1][j-1] ---*/
				final long dM=prevM[j-1], dI=prevI[j-1], dD=prevD[j-1];
				long diag=dM; int diagState=0;
				if(scoreOf(dI)>scoreOf(diag)){diag=dI; diagState=1;}
				if(scoreOf(dD)>scoreOf(diag)){diag=dD; diagState=2;}
				final long mWord;
				if(isMatch){
					// Consecutive match (came from M and prev was also a match = time 0) gets +100.
					// First match (from I/D, or from M after a sub run = time>0) gets +70.
					final boolean consecutiveMatch=(diagState==0 && timeOf(diag)==0);
					final long matchPts=consecutiveMatch ? MATCH : MATCH_FIRST;
					mWord=pack(scoreOf(diag)+matchPts, 0, startOf(diag));
				}else if(hasN){
					mWord=pack(scoreOf(diag)+N_SCORE, 0, startOf(diag));
				}else{// substitution: the sub-run continues only if we came from M
					final int subRun=(diagState==0 ? timeOf(diag) : 0);
					mWord=pack(scoreOf(diag)-subCost(subRun), subRun+1, startOf(diag));
				}

				/*--- I: insertion, up move (i-1,j): open from M, extend from I ---*/
				final long openI=prevM[j], extI=prevI[j];
				final long iOpen=scoreOf(openI)-insCost(0);
				final long iExt=scoreOf(extI)-insCost(timeOf(extI));
				final long iWord=(iOpen>=iExt)
					? pack(iOpen, 1, startOf(openI))
					: pack(iExt, timeOf(extI)+1, startOf(extI));

				/*--- D: deletion, left move (i,j-1): open from M, extend from D ---*/
				final long openD=currM[j-1], extD=currD[j-1];
				final long dOpen=scoreOf(openD)-delCost(0);
				final long dExt=scoreOf(extD)-delCost(timeOf(extD));
				// Extra penalty for deleting across N (discourages cross-scaffold alignment)
				final long nPenalty=(r=='N' ? POINTS_DEL_REF_N : 0);
				final long dWord=(dOpen-nPenalty>=dExt-nPenalty)
					? pack(dOpen-nPenalty, 1, startOf(openD))
					: pack(dExt-nPenalty, timeOf(extD)+1, startOf(extD));

				currM[j]=mWord; currI[j]=iWord; currD[j]=dWord;
				} // end non-gap branch

				// Best cell ending in this row (query fully consumed when i==qLen).
				// Ending in D would charge trailing deletions; trailing ref is free in
				// glocal, so we only consider M and I for the endpoint.
				// In a gap column, M and I are BAD so D could be best, but we don't
				// want to end on a gap symbol — the real endpoint is after the gap.
				final long cm=currM[j], ci=currI[j];
				final long cellBest=(scoreOf(cm)>=scoreOf(ci) ? cm : ci);
				if(scoreOf(cellBest)>rowBestScore){
					rowBestScore=scoreOf(cellBest); rowBestPos=j; rowBestWord=cellBest;
				}
			}
			mloops+=(bandEnd-bandStart+1);

			// Update row-best into the band-tracking best
			if(scoreOf(rowBestWord)>scoreOf(bestWord) || i==qLen){
				bestScore=rowBestScore; bestWord=rowBestWord;
			}
			bestPos=rowBestPos; // always track for band centering

			if(i==qLen){bestScore=rowBestScore; bestWord=rowBestWord;}

			// Swap all three states' rows
			long[] t;
			t=prevM; prevM=currM; currM=t;
			t=prevI; prevI=currI; currI=t;
			t=prevD; prevD=currD; currD=t;
		}
		loops.addAndGet(mloops);

		final int rStart=startOf(bestWord);
		final int rStop=bestPos-1;
		if(posVector!=null){
			posVector[0]=rStart;
			posVector[1]=rStop;
			if(posVector.length>2){posVector[2]=(int)bestScore;}
		}
		return identity(bestWord, bestScore, bestPos, qLen, rLen);
	}

	/*--------------------------------------------------------------*/
	/*----  Mode 2: trace-recording DP (cold path, no banding)  ----*/
	/*--------------------------------------------------------------*/

	/**
	 * Same DP as alignStatic but stores every row's M/I/D cells into an
	 * interleaved LongList for TracerAffine to walk backward through.
	 * Full DP (no adaptive banding) for correctness; only called for the
	 * 1-2 winning candidates that need a match string.
	 */
	public static final float alignWithTrace(byte[] query, byte[] ref,
			int[] posVector, LongList trace){
		final int qLen=query.length, rLen=ref.length;
		assert(rLen<=POSITION_MASK) : "Ref too long: "+rLen+">"+POSITION_MASK;

		// 6 rolling arrays: 3 states x {prev,curr}
		long[] prevM=new long[rLen+1], currM=new long[rLen+1];
		long[] prevI=new long[rLen+1], currI=new long[rLen+1];
		long[] prevD=new long[rLen+1], currD=new long[rLen+1];

		// Row 0 init: glocal, alignment starts at any ref column with score 0.
		// Time=1 so first match scores MATCH_FIRST (70), not MATCH (100).
		for(int j=0; j<=rLen; j++){
			prevM[j]=pack(0, 1, j);
			prevI[j]=BAD;
			prevD[j]=BAD;
		}

		// Store row 0 in trace
		trace.clear();
		int lastHeaderIdx=0;
		trace.add(TracerAffine.packHeader(0, 0, 0));
		for(int j=0; j<=rLen; j++){
			trace.add(prevM[j]);
			trace.add(prevI[j]);
			trace.add(prevD[j]);
		}

		long bestScore=BAD; int bestPos=0; long bestWord=BAD;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			// Column 0: leading query bases with no ref consumed are insertions.
			currM[0]=BAD;
			currD[0]=BAD;
			currI[0]=pack(0-insCost(0)-(long)(i-1)*insCost(1), i, 0);

			long rowBestScore=BAD; int rowBestPos=0; long rowBestWord=BAD;

			for(int j=1; j<=rLen; j++){
				final byte r=ref[j-1];
				final boolean isGap=(r==GAP_SYMBOL);

				if(isGap){
					currM[j]=BAD;
					currI[j]=BAD;
					final long openD=currM[j-1], extD=currD[j-1];
					final long dOpen=scoreOf(openD)-gapSymbolCost();
					final long dExt=scoreOf(extD)-gapSymbolCost();
					currD[j]=(dOpen>=dExt)
						? pack(dOpen, 1, startOf(openD))
						: pack(dExt, timeOf(extD)+1, startOf(extD));
				}else{
				final boolean isMatch=(q==r && r!='N' && q!='N');
				final boolean hasN=(q=='N' || r=='N');

				/*--- M: diagonal from best of {M,I,D}[i-1][j-1] ---*/
				final long dM=prevM[j-1], dI=prevI[j-1], dD=prevD[j-1];
				long diag=dM; int diagState=0;
				if(scoreOf(dI)>scoreOf(diag)){diag=dI; diagState=1;}
				if(scoreOf(dD)>scoreOf(diag)){diag=dD; diagState=2;}
				final long mWord;
				if(isMatch){
					final boolean consecutiveMatch=(diagState==0 && timeOf(diag)==0);
					final long matchPts=consecutiveMatch ? MATCH : MATCH_FIRST;
					mWord=pack(scoreOf(diag)+matchPts, 0, startOf(diag));
				}else if(hasN){
					mWord=pack(scoreOf(diag)+N_SCORE, 0, startOf(diag));
				}else{
					final int subRun=(diagState==0 ? timeOf(diag) : 0);
					mWord=pack(scoreOf(diag)-subCost(subRun), subRun+1, startOf(diag));
				}

				/*--- I: up move from (i-1,j) ---*/
				final long openI=prevM[j], extI=prevI[j];
				final long iOpen=scoreOf(openI)-insCost(0);
				final long iExt=scoreOf(extI)-insCost(timeOf(extI));
				final long iWord=(iOpen>=iExt)
					? pack(iOpen, 1, startOf(openI))
					: pack(iExt, timeOf(extI)+1, startOf(extI));

				/*--- D: left move from (i,j-1) ---*/
				final long openD=currM[j-1], extD=currD[j-1];
				final long dOpen=scoreOf(openD)-delCost(0);
				final long dExt=scoreOf(extD)-delCost(timeOf(extD));
				final long nPenalty=(r=='N' ? POINTS_DEL_REF_N : 0);
				final long dWord=(dOpen-nPenalty>=dExt-nPenalty)
					? pack(dOpen-nPenalty, 1, startOf(openD))
					: pack(dExt-nPenalty, timeOf(extD)+1, startOf(extD));

				currM[j]=mWord; currI[j]=iWord; currD[j]=dWord;
				} // end non-gap branch

				final long cm=currM[j], ci=currI[j];
				final long cellBest=(scoreOf(cm)>=scoreOf(ci) ? cm : ci);
				if(scoreOf(cellBest)>rowBestScore){
					rowBestScore=scoreOf(cellBest); rowBestPos=j; rowBestWord=cellBest;
				}
			}

			// Store row i in trace (before swap)
			int newHeaderIdx=trace.size;
			int dist=newHeaderIdx-lastHeaderIdx;
			trace.add(TracerAffine.packHeader(i, 0, dist));
			for(int j=0; j<=rLen; j++){
				trace.add(currM[j]);
				trace.add(currI[j]);
				trace.add(currD[j]);
			}
			lastHeaderIdx=newHeaderIdx;

			if(i==qLen){bestScore=rowBestScore; bestWord=rowBestWord; bestPos=rowBestPos;}

			// Swap all three states' rows
			long[] t;
			t=prevM; prevM=currM; currM=t;
			t=prevI; prevI=currI; currI=t;
			t=prevD; prevD=currD; currD=t;
		}

		final int rStart=startOf(bestWord);
		final int rStop=bestPos-1;
		if(posVector!=null){
			posVector[0]=rStart;
			posVector[1]=rStop;
			if(posVector.length>2){posVector[2]=(int)bestScore;}
		}
		return identity(bestWord, bestScore, bestPos, qLen, rLen);
	}

	/**
	 * v1 STUB. Exact affine identity cannot be inverted from the score (variable
	 * open/extend costs) and the deletion-count field is gone (now timeInState).
	 * NEXT: carry a parallel match/length counter along the chosen path (or pack
	 * one) and return matches/length. For now, a rough score-based estimate so
	 * the pipeline runs end to end. Do NOT trust this number yet.
	 */
	private static float identity(long bestWord, long bestScore, int endPos, int qLen, int rLen){
		// crude placeholder assuming no indels (match=+1, sub~=-2): matches ~ (score+qLen)/2
		final float est=(bestScore+qLen)/(2f*Math.max(1, qLen));
		return Math.max(0f, Math.min(1f, est)); // TODO: replace with exact matches/length
	}

	/*--------------------------------------------------------------*/
	/*----------------          Packing             ----------------*/
	/*--------------------------------------------------------------*/

	private static final int POSITION_BITS=19;
	private static final int TIME_BITS=18;
	private static final int SCORE_SHIFT=POSITION_BITS+TIME_BITS; // 37
	private static final long POSITION_MASK=(1L<<POSITION_BITS)-1;
	private static final long TIME_MASK=(1L<<TIME_BITS)-1;

	// BBMap's exact scoring constants. Match the 4SA for debuggability.
	private static final long MATCH=100;
	private static final long MATCH_FIRST=70;
	private static final long N_SCORE=0;
	static final long BAD=Long.MIN_VALUE/2;

	/** Pack a raw (unshifted) score, timeInState, and start into one cell word. */
	static long pack(long score, int timeInState, int start){
		return (score<<SCORE_SHIFT)
			| (((long)timeInState & POSITION_MASK)<<POSITION_BITS)
			| (start & POSITION_MASK);
	}
	static long scoreOf(long w){return w>>SCORE_SHIFT;}                       // arithmetic; sign-extends
	static int timeOf(long w){return (int)((w>>>POSITION_BITS) & TIME_MASK);}
	static int startOf(long w){return (int)(w & POSITION_MASK);}

}
