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
 *  - Adaptive band (from ScrabbleAligner) in both alignStatic (Mode 1) and
 *    alignWithTrace (Mode 2); USE_BANDING / TRACE_BANDING toggle full DP.
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
		if(timeInState<5) return POINTS_SUB2;
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
	// Public so bbmap2's startup cross-assert (verifyConstants) can pin all four homes of
	// this value together; scoreFromMatch's gap-aware recomputation went stale once already.
	public static final long GAP_SYMBOL_PENALTY=2;

	/** The gap symbol byte used by GappedReference. */
	public static final byte GAP_SYMBOL='-';

	/*--------------------------------------------------------------*/
	/*----------------       Banding               ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether to use adaptive banding. Off = full DP (for verification). */
	public static boolean USE_BANDING=true;

	/** Band-growth tunables (see alignStatic): the initial band shrinks over rows as
	 *  bandWidth0 + max(BAND_CONST + bandWidth0*BAND_MULT - maxDrift*i, dynamicBW). The default
	 *  BAND_MULT=12 makes the early band ~157 wide (near-full DP for the first ~74 rows) — the
	 *  dominant CPU cost. Lowering these tightens the band (faster) at the risk of missing indels
	 *  whose path drifts outside it; sweep against the accuracy guard before committing a change. */
	public static int BAND_CONST=16, BAND_MULT=4;

	/** Diagnostic counters (DB_PROFILE): decideBandwidth call count + total iterations. */
	public static boolean DB_PROFILE=false;
	public static final AtomicLong DB_CALLS=new AtomicLong(), DB_ITERS=new AtomicLong();

	/** Band the trace-recording DP (Mode 2), mirroring alignStatic's adaptive band. Default OFF:
	 *  the full-DP winner alignment is load-bearing for INDEL RECOVERY — when the gap-array heuristic
	 *  misses a large deletion, full DP still finds it by searching the whole window, but a narrow band
	 *  soft-clips it (VERIFY_TRACE: lambda deletion reads 91.8%->57.9%, 95% of divergences band-WORSE).
	 *  Banding is only safe with a full-DP fallback for deficient reads (see plan). Keep off pending that. */
	public static boolean TRACE_BANDING=false;

	/** Diagnostic (VERIFY_TRACE): run the UNBANDED trace DP alongside the banded one and count how
	 *  often their (start,stop,score) disagree — i.e. how often the band misses the full-DP optimum.
	 *  Doubles Mode 2 cost; off by default. A divergence is an accepted band/DP tradeoff, NOT a bug:
	 *  the same band already scored this window in alignStatic (scoreOnly) during candidate selection. */
	public static boolean VERIFY_TRACE=false;
	public static final AtomicLong VT_TOTAL=new AtomicLong(), VT_DIVERGENCES=new AtomicLong();
	/** Divergence breakdown: banded score strictly WORSE than full-DP (band missed the optimum)
	 *  vs equal score but different start/stop (tie-break / traceback edge). */
	public static final AtomicLong VT_WORSE=new AtomicLong(), VT_TIE_DIFFPOS=new AtomicLong();

	/** Decide initial bandwidth from sequence lengths + a quick identity scan. */
	private static int decideBandwidth(byte[] query, byte[] ref){
		final int qLen=query.length, rLen=ref.length;
		int bandwidth=shared.Tools.mid(7, 1+Math.max(qLen, rLen)/32, 20+(int)Math.sqrt(rLen)/8);
		int subs=0;
		int i=0, minlen=Math.min(qLen, rLen);
		for(; i<minlen && subs<bandwidth; i++){
			if(query[i]!=ref[i]){subs++;}
		}
		if(DB_PROFILE){DB_CALLS.incrementAndGet(); DB_ITERS.addAndGet(i);}
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
		// Use posVector[0] as an alignment hint if provided (expected start in ref window)
		int center=(posVector!=null && posVector[0]>0) ? posVector[0] : 0;
		int dynamicBW=0, deltaBW=0;

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
		long col0InsCost=0;
		// Valid range of the prev arrays: {0} union [prevLo,prevHi] was written by the
		// previous row (row 0 wrote everything). Cells OUTSIDE it hold garbage: two-rows-ago
		// values, row-0 free-start words, or unpacked-zero phantom match-streak cells. The
		// band is not monotonic (dynamicBW can outgrow center drift on either edge), so each
		// row must BAD-clear the part of its read range [bandStart-1,bandEnd] the previous
		// row never wrote. ScrabbleAffineSP clears both edges for the same reason.
		int prevLo=0, prevHi=rLen;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			// Adaptive band: widen on mismatches, narrow on matches
			if(banded){
				assert(bandEnd<=rLen) : "Pre-band bandEnd="+bandEnd+", rLen="+rLen+
					", center="+center+", i="+i+", qLen="+qLen;
				final boolean nextMatch=(bestPos>=0 && bestPos<rLen && q==ref[Math.min(rLen-1, bestPos)]);
				if(nextMatch){
					deltaBW=(deltaBW<0 ? Math.max(-maxDynamic, deltaBW*2) : -2);
				}else{
					deltaBW=shared.Tools.mid(1, (maxDynamic-dynamicBW)/2, 8);
				}
				dynamicBW=shared.Tools.mid(0, dynamicBW+deltaBW, maxDynamic);
				final int bandWidth=bandWidth0+Math.max(BAND_CONST+bandWidth0*BAND_MULT-maxDrift*i, dynamicBW);
				final int quarterBand=bandWidth/4;
				final int drift=shared.Tools.mid(-1, bestPos-center, maxDrift);
				center=shared.Tools.mid(1, center+1+drift, rLen);
				bandStart=Math.max(1, center-bandWidth+quarterBand);
				bandEnd=Math.min(rLen, center+bandWidth+quarterBand);
				assert(bandEnd<=rLen) : "bandEnd="+bandEnd+", rLen="+rLen+", center="+center;
			}

			// Column 0: leading query bases with no ref consumed are insertions.
			col0InsCost-=insCost(i-1);
			currM[0]=BAD;
			currD[0]=BAD;
			currI[0]=pack(col0InsCost, i, 0);

			// Clear stale prev cells this row can read but the previous row never wrote:
			// left exposure (band moved left past the old bandStart-1 clear) and right
			// exposure (band grew right past the old bandEnd). Col 0 is always valid.
			for(int j=Math.max(1, bandStart-1); j<prevLo; j++){
				prevM[j]=BAD; prevI[j]=BAD; prevD[j]=BAD;
			}
			for(int j=prevHi+1; j<=bandEnd; j++){
				prevM[j]=BAD; prevI[j]=BAD; prevD[j]=BAD;
			}

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
					// Consecutive match: came from M, prev was also a match (time 0), and prev wasn't N.
					// N sets time=0 (to preserve sub-run mapping) but shouldn't start a match streak.
					final boolean prevWasN=(diagState==0 && timeOf(diag)==0 && i>=2 && j>=2
						&& (query[i-2]=='N' || ref[j-2]=='N'));
					final boolean consecutiveMatch=(diagState==0 && timeOf(diag)==0 && !prevWasN);
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
			// This row wrote col 0, the bandStart-1 BAD clear, and [bandStart,bandEnd].
			prevLo=bandStart-1; prevHi=bandEnd;
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
	/*----  Mode 2: trace-recording DP (cold path, banded)  --------*/
	/*--------------------------------------------------------------*/

	/**
	 * Trace-recording DP that produces a match string via TracerAffine. Banded by
	 * default (TRACE_BANDING), using the SAME adaptive band as alignStatic so the
	 * traced alignment matches the score alignStatic (scoreOnly) used to pick this
	 * winner. Only called for the 1-2 winning candidates that need a match string,
	 * so ~once per read. When VERIFY_TRACE is on, also runs the unbanded DP and
	 * counts (start,stop,score) divergences — an accepted band/DP tradeoff, not a bug.
	 * Reads posVector[0] as the band-center hint (expected read-start in the window).
	 */
	public static final float alignWithTrace(byte[] query, byte[] ref,
			int[] posVector, LongList trace){
		if(VERIFY_TRACE){
			// Unbanded reference run first (trace overwritten by the banded run below).
			int[] refPos=new int[3];
			if(posVector!=null && posVector.length>0){refPos[0]=posVector[0];}
			alignWithTraceImpl(query, ref, refPos, trace, false, null);
			float id=bandedTraceWithFallback(query, ref, posVector, trace);
			VT_TOTAL.incrementAndGet();
			if(posVector!=null && posVector.length>2 &&
					(refPos[0]!=posVector[0] || refPos[1]!=posVector[1] || refPos[2]!=posVector[2])){
				VT_DIVERGENCES.incrementAndGet();
				if(posVector[2]<refPos[2]){VT_WORSE.incrementAndGet();}
				else if(posVector[2]==refPos[2]){VT_TIE_DIFFPOS.incrementAndGet();}
			}
			return id;
		}
		return bandedTraceWithFallback(query, ref, posVector, trace);
	}

	/** Diagnostic: reads whose banded trace pressed the band edge and fell back to full DP. */
	public static final AtomicLong VT_FALLBACKS=new AtomicLong();

	/**
	 * Banded trace with a full-DP safety net (BBMap's recipe): run the banded DP; if the best
	 * path pressed a real band edge, the band may have clipped the true alignment (typically a
	 * large indel the gap-array heuristic missed), so re-run unbanded. Clean diagonal reads never
	 * touch the edge and pay only the cheap banded pass; only edge-touching reads pay full DP.
	 * When TRACE_BANDING is off, always full DP.
	 */
	private static final float bandedTraceWithFallback(byte[] query, byte[] ref,
			int[] posVector, LongList trace){
		if(!TRACE_BANDING){return alignWithTraceImpl(query, ref, posVector, trace, false, null);}
		boolean[] edge={false};
		float id=alignWithTraceImpl(query, ref, posVector, trace, true, edge);
		if(edge[0]){
			if(DB_PROFILE){VT_FALLBACKS.incrementAndGet();}
			id=alignWithTraceImpl(query, ref, posVector, trace, false, null);
		}
		return id;
	}

	/**
	 * Same DP as alignStatic (identical recurrence and adaptive band) but stores each
	 * row's M/I/D cells into an interleaved LongList for TracerAffine. When banded, only
	 * the [bandStart,bandEnd] block of each row is stored, and the header records
	 * startCol=bandStart so TracerAffine indexes cells relative to it (BAD outside the
	 * band). Row 0 is stored full so row 1's band can read its diagonal predecessors.
	 */
	private static final float alignWithTraceImpl(byte[] query, byte[] ref,
			int[] posVector, LongList trace, final boolean banded, final boolean[] edgeFlag){
		final int qLen=query.length, rLen=ref.length;
		assert(rLen<=POSITION_MASK) : "Ref too long: "+rLen+">"+POSITION_MASK;

		// Banding parameters (mirror alignStatic exactly).
		final int bandWidth0=banded ? decideBandwidth(query, ref) : 0;
		final int maxDrift=2, maxDynamic=banded ? (bandWidth0*12)/4 : 0;
		int center=(posVector!=null && posVector[0]>0) ? posVector[0] : 0;
		int dynamicBW=0, deltaBW=0;

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

		// Store row 0 in trace (full width — one row, cheap; row 1's band reads into it).
		trace.clear();
		int lastHeaderIdx=0;
		trace.add(TracerAffine.packHeader(0, 0, 0));
		for(int j=0; j<=rLen; j++){
			trace.add(prevM[j]);
			trace.add(prevI[j]);
			trace.add(prevD[j]);
		}

		long bestScore=BAD; int bestPos=0; long bestWord=BAD;
		int bandStart=1, bandEnd=rLen;
		long col0InsCost=0;
		// Valid range of the prev arrays (same stale-fill guard as alignStatic).
		int prevLo=0, prevHi=rLen;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			// Adaptive band: widen on mismatches, narrow on matches (mirrors alignStatic).
			if(banded){
				final boolean nextMatch=(bestPos>=0 && bestPos<rLen && q==ref[Math.min(rLen-1, bestPos)]);
				if(nextMatch){
					deltaBW=(deltaBW<0 ? Math.max(-maxDynamic, deltaBW*2) : -2);
				}else{
					deltaBW=shared.Tools.mid(1, (maxDynamic-dynamicBW)/2, 8);
				}
				dynamicBW=shared.Tools.mid(0, dynamicBW+deltaBW, maxDynamic);
				final int bandWidth=bandWidth0+Math.max(BAND_CONST+bandWidth0*BAND_MULT-maxDrift*i, dynamicBW);
				final int quarterBand=bandWidth/4;
				final int drift=shared.Tools.mid(-1, bestPos-center, maxDrift);
				center=shared.Tools.mid(1, center+1+drift, rLen);
				bandStart=Math.max(1, center-bandWidth+quarterBand);
				bandEnd=Math.min(rLen, center+bandWidth+quarterBand);
			}else{
				bandStart=1; bandEnd=rLen;
			}

			// Column 0: leading query bases with no ref consumed are insertions.
			col0InsCost-=insCost(i-1);
			currM[0]=BAD;
			currD[0]=BAD;
			currI[0]=pack(col0InsCost, i, 0);

			// Clear stale prev cells this row can read but the previous row never wrote
			// (same guard as alignStatic; no-op when unbanded).
			for(int j=Math.max(1, bandStart-1); j<prevLo; j++){
				prevM[j]=BAD; prevI[j]=BAD; prevD[j]=BAD;
			}
			for(int j=prevHi+1; j<=bandEnd; j++){
				prevM[j]=BAD; prevI[j]=BAD; prevD[j]=BAD;
			}

			// Clear stale rolling-array data at the left band edge.
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
					final boolean prevWasN=(diagState==0 && timeOf(diag)==0 && i>=2 && j>=2
						&& (query[i-2]=='N' || ref[j-2]=='N'));
					final boolean consecutiveMatch=(diagState==0 && timeOf(diag)==0 && !prevWasN);
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

			// Edge-touch detection (fallback trigger): the row's best cell pressed a REAL band
			// edge (not the ref boundary), so the true path may want to leave the band.
			if(banded && edgeFlag!=null && rowBestPos>0 &&
					((bandStart>1 && rowBestPos<=bandStart) || (bandEnd<rLen && rowBestPos>=bandEnd))){
				edgeFlag[0]=true;
			}

			// Store row i in trace (before swap). Banded: only [bandStart,bandEnd],
			// header startCol=bandStart. Unbanded: full [0,rLen], startCol=0.
			final int storeStart=banded ? bandStart : 0;
			int newHeaderIdx=trace.size;
			int dist=newHeaderIdx-lastHeaderIdx;
			trace.add(TracerAffine.packHeader(i, storeStart, dist));
			for(int j=storeStart; j<=bandEnd; j++){
				trace.add(currM[j]);
				trace.add(currI[j]);
				trace.add(currD[j]);
			}
			lastHeaderIdx=newHeaderIdx;

			// Track the best cell for band centering EVERY row (mirrors alignStatic). The
			// final answer is the best cell of the last row (glocal: query fully consumed).
			if(scoreOf(rowBestWord)>scoreOf(bestWord) || i==qLen){
				bestScore=rowBestScore; bestWord=rowBestWord;
			}
			bestPos=rowBestPos; // ALWAYS: the band re-centers on this next row (drift=bestPos-center)
			if(i==qLen){bestScore=rowBestScore; bestWord=rowBestWord;}

			// Swap all three states' rows
			long[] t;
			t=prevM; prevM=currM; currM=t;
			t=prevI; prevI=currI; currI=t;
			t=prevD; prevD=currD; currD=t;
			// This row wrote col 0, the bandStart-1 BAD clear, and [bandStart,bandEnd].
			prevLo=bandStart-1; prevHi=bandEnd;
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
	public static final long BAD=Long.MIN_VALUE/2;

	/** Pack a raw (unshifted) score, timeInState, and start into one cell word. */
	public static long pack(long score, int timeInState, int start){
		return (score<<SCORE_SHIFT)
			| (((long)timeInState & POSITION_MASK)<<POSITION_BITS)
			| (start & POSITION_MASK);
	}
	public static long scoreOf(long w){return w>>SCORE_SHIFT;}
	public static int timeOf(long w){return (int)((w>>>POSITION_BITS) & TIME_MASK);}
	public static int startOf(long w){return (int)(w & POSITION_MASK);}

}
