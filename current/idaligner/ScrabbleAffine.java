package idaligner;

import java.util.concurrent.atomic.AtomicLong;

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

	/** v1 self-test: prints start/stop/score for a few hand-checkable cases. */
	public static void main(String[] args){
		selfTest("ACGT", "ACGT");      // perfect:     expect start=0 stop=3 score=4
		selfTest("ACGT", "AGGT");      // 1 sub:       expect start=0 stop=3 score=1  (3 - subOpen 2)
		selfTest("ACT",  "ACGT");      // 1 del:       expect start=0 stop=3 score=1  (3 - delOpen 2)
		selfTest("CGT",  "ACGT");      // start!=0:    expect start=1 stop=3 score=3  (no-traceback start)
		selfTest("ACGT", "ACT");       // 1 ins:       expect start=0 stop=2 score=1  (3 - insOpen 2)
		selfTest("AAATTT", "AAAGGTTT");// del open+ext: expect start=0 stop=7 score=3  (6 - delOpen2 - delExt1)
	}
	private static void selfTest(String q, String r){
		int[] pos=new int[3];
		float id=alignStatic(q.getBytes(), r.getBytes(), pos);
		System.out.println(q+" vs "+r+"  ->  start="+pos[0]+" stop="+pos[1]+" score="+pos[2]+"  (id stub="+id+")");
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
	// Returns a positive penalty (subtracted from score). Edit freely to retune.
	private static long subCost(int timeInState){return timeInState==0 ? 2 : 1;}
	private static long insCost(int timeInState){return timeInState==0 ? 2 : 1;}
	private static long delCost(int timeInState){return timeInState==0 ? 2 : 1;}

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

		// 6 rolling arrays: 3 states x {prev,curr}
		long[] prevM=new long[rLen+1], currM=new long[rLen+1];
		long[] prevI=new long[rLen+1], currI=new long[rLen+1];
		long[] prevD=new long[rLen+1], currD=new long[rLen+1];

		// Row 0: the alignment may BEGIN at any ref column (glocal) -> score 0, start=j.
		// I and D cannot exist before any query base is consumed.
		for(int j=0; j<=rLen; j++){
			prevM[j]=pack(0, 0, j);
			prevI[j]=BAD;
			prevD[j]=BAD;
		}

		long bestScore=BAD; int bestPos=0; long bestWord=BAD;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			// Column 0: leading query bases with no ref consumed are insertions.
			currM[0]=BAD;
			currD[0]=BAD;
			currI[0]=pack(0 - insCost(0) - (long)(i-1)*insCost(1), i, 0); // open + (i-1) extends

			long rowBestScore=BAD; int rowBestPos=0; long rowBestWord=BAD;

			for(int j=1; j<=rLen; j++){
				final byte r=ref[j-1];
				final boolean isMatch=(q==r && q!='N');
				final boolean hasN=(q=='N' || r=='N');

				/*--- M: diagonal move from the best of {M,I,D}[i-1][j-1] ---*/
				final long dM=prevM[j-1], dI=prevI[j-1], dD=prevD[j-1];
				long diag=dM; int diagState=0;
				if(scoreOf(dI)>scoreOf(diag)){diag=dI; diagState=1;}
				if(scoreOf(dD)>scoreOf(diag)){diag=dD; diagState=2;}
				final long mWord;
				if(isMatch){
					mWord=pack(scoreOf(diag)+MATCH, 0, startOf(diag));
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
				final long dWord=(dOpen>=dExt)
					? pack(dOpen, 1, startOf(openD))
					: pack(dExt, timeOf(extD)+1, startOf(extD));

				currM[j]=mWord; currI[j]=iWord; currD[j]=dWord;

				// Best cell ending in this row (query fully consumed when i==qLen).
				// Ending in D would charge trailing deletions; trailing ref is free in
				// glocal, so we only consider M and I for the endpoint.
				final long cellBest=(scoreOf(mWord)>=scoreOf(iWord) ? mWord : iWord);
				if(scoreOf(cellBest)>rowBestScore){
					rowBestScore=scoreOf(cellBest); rowBestPos=j; rowBestWord=cellBest;
				}
			}
			mloops+=rLen;

			if(i==qLen){bestScore=rowBestScore; bestPos=rowBestPos; bestWord=rowBestWord;}

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

	private static final int POSITION_BITS=21;
	private static final int TIME_BITS=21;
	private static final int SCORE_SHIFT=POSITION_BITS+TIME_BITS; // 42
	private static final long POSITION_MASK=(1L<<POSITION_BITS)-1;

	private static final long MATCH=1;     // raw score delta (+1); pack() does the shifting
	private static final long N_SCORE=0;
	private static final long BAD=Long.MIN_VALUE/2;

	/** Pack a raw (unshifted) score, timeInState, and start into one cell word. */
	private static long pack(long score, int timeInState, int start){
		return (score<<SCORE_SHIFT)
			| (((long)timeInState & POSITION_MASK)<<POSITION_BITS)
			| (start & POSITION_MASK);
	}
	private static long scoreOf(long w){return w>>SCORE_SHIFT;}                       // arithmetic; sign-extends
	private static int timeOf(long w){return (int)((w>>>POSITION_BITS) & POSITION_MASK);}
	private static int startOf(long w){return (int)(w & POSITION_MASK);}

}
