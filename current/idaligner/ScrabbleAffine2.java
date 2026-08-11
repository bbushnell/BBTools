package idaligner;

import java.util.concurrent.atomic.AtomicLong;

import structures.LongList;

/**
 * ScrabbleAffine with a shadow band that starts wide and contracts.
 *
 * The shadow band guarantees that a wider window never scores worse than a
 * narrower one. It starts at ±qLen/2 from center and contracts from each side
 * when the edge cells are significantly worse than the best cell in the row.
 * Once the shadow contracts to the "desired" band (the narrow adaptive band
 * from ScrabbleAffine), it disappears and the aligner behaves identically.
 *
 * @author Ady
 * @date July 3, 2026
 */
public class ScrabbleAffine2{

	public ScrabbleAffine2(){}

	public static long loops(){return loops.get();}
	public static void resetLoops(){loops.set(0);}
	private static final java.util.concurrent.atomic.AtomicLong loops=new java.util.concurrent.atomic.AtomicLong(0);

	/*--------------------------------------------------------------*/
	/*--------  Cost hooks (identical to ScrabbleAffine)  ----------*/
	/*--------------------------------------------------------------*/

	private static final long POINTS_SUB=127, POINTS_SUBR=147;
	private static final long POINTS_SUB2=51, POINTS_SUB3=25;
	private static long subCost(int timeInState){
		if(timeInState==0) return POINTS_SUB;
		if(timeInState==1) return POINTS_SUB2;
		return POINTS_SUB3;
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
	private static long gapSymbolCost(){return GAP_SYMBOL_PENALTY;}
	// Public for the startup cross-assert; see GappedReference.GAP_SYMBOL_COST (canonical home).
	public static final long GAP_SYMBOL_PENALTY=2;
	public static final byte GAP_SYMBOL='-';

	/*--------------------------------------------------------------*/
	/*----------------       Shadow Band Config    ----------------*/
	/*--------------------------------------------------------------*/

	// How much worse an edge cell must be (vs row best) to allow contraction.
	// ~3 mismatches vs matches: 3*(127+100) ≈ 681. Start at 400, tune later.
	public static long CONTRACTION_THRESHOLD=1000;

	/*--------------------------------------------------------------*/
	/*----------------       Desired Band (from SA1) ---------------*/
	/*--------------------------------------------------------------*/

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
		if(posVector==null && query.length>ref.length){
			byte[] tmp=query; query=ref; ref=tmp;
		}
		final int qLen=query.length, rLen=ref.length;
		assert(rLen<=POSITION_MASK) : "Ref too long: "+rLen+">"+POSITION_MASK;

		// Desired band parameters (same as ScrabbleAffine)
		final int bandWidth0=decideBandwidth(query, ref);
		final int maxDrift=2, maxDynamic=bandWidth0*3;
		int center=(posVector!=null && posVector[0]>0) ? posVector[0] : 0;
		int dynamicBW=0, deltaBW=0;

		// Shadow band: starts covering the full plausible alignment region.
		// The read could start anywhere from column 0 to rLen-qLen, and end
		// at that start + qLen + some indel allowance. Cover the whole ref.
		int shadowLeft=1;
		int shadowRight=rLen;
		boolean shadowActive=true;

		// 6 rolling arrays
		long[] prevM=new long[rLen+1], currM=new long[rLen+1];
		long[] prevI=new long[rLen+1], currI=new long[rLen+1];
		long[] prevD=new long[rLen+1], currD=new long[rLen+1];

		for(int j=0; j<=rLen; j++){
			prevM[j]=pack(0, 1, j);
			prevI[j]=BAD;
			prevD[j]=BAD;
		}

		long bestScore=BAD; int bestPos=0; long bestWord=BAD;
		// Valid range of the prev arrays (same stale-fill guard as ScrabbleAffine):
		// cells outside {0} union [prevLo,prevHi] were not written by the previous row
		// and hold garbage (two-rows-ago values, row-0 free-start words, or all-zero
		// phantom match-streak words). The desired band is not monotonic, so each row
		// BAD-clears the part of its read range the previous row never wrote.
		int prevLo=0, prevHi=rLen;

		for(int i=1; i<=qLen; i++){
			final byte q=query[i-1];

			// Compute desired band (same logic as ScrabbleAffine)
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
			center=shared.Tools.mid(1, center+1+drift, rLen);
			int desiredStart=Math.max(1, center-bandWidth+quarterBand);
			int desiredEnd=Math.min(rLen, center+bandWidth+quarterBand);

			// Effective band: union of shadow and desired (shadow only extends, never shrinks desired)
			int bandStart, bandEnd;
			if(shadowActive){
				bandStart=Math.min(shadowLeft, desiredStart);
				bandEnd=Math.max(shadowRight, desiredEnd);
			}else{
				bandStart=desiredStart;
				bandEnd=desiredEnd;
			}

			// Column 0
			currM[0]=BAD;
			currD[0]=BAD;
			currI[0]=pack(0 - insCost(0) - (long)(i-1)*insCost(1), i, 0);

			// Clear stale prev cells this row can read but the previous row never wrote.
			for(int j=Math.max(1, bandStart-1); j<prevLo; j++){
				prevM[j]=BAD; prevI[j]=BAD; prevD[j]=BAD;
			}
			for(int j=prevHi+1; j<=bandEnd; j++){
				prevM[j]=BAD; prevI[j]=BAD; prevD[j]=BAD;
			}

			// Clear stale data at band edges
			if(bandStart>1){
				currM[bandStart-1]=BAD;
				currI[bandStart-1]=BAD;
				currD[bandStart-1]=BAD;
			}

			long rowBestScore=BAD; int rowBestPos=0; long rowBestWord=BAD;
			long mloops=0;

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

					final long openI=prevM[j], extI=prevI[j];
					final long iOpen=scoreOf(openI)-insCost(0);
					final long iExt=scoreOf(extI)-insCost(timeOf(extI));
					final long iWord=(iOpen>=iExt)
						? pack(iOpen, 1, startOf(openI))
						: pack(iExt, timeOf(extI)+1, startOf(extI));

					final long openD=currM[j-1], extD=currD[j-1];
					final long dOpen=scoreOf(openD)-delCost(0);
					final long dExt=scoreOf(extD)-delCost(timeOf(extD));
					final long nPenalty=(r=='N' ? POINTS_DEL_REF_N : 0);
					final long dWord=(dOpen-nPenalty>=dExt-nPenalty)
						? pack(dOpen-nPenalty, 1, startOf(openD))
						: pack(dExt-nPenalty, timeOf(extD)+1, startOf(extD));

					currM[j]=mWord; currI[j]=iWord; currD[j]=dWord;
				}

				final long cm=currM[j], ci=currI[j];
				mloops++;
				final long cellBest=(scoreOf(cm)>=scoreOf(ci) ? cm : ci);
				if(scoreOf(cellBest)>rowBestScore){
					rowBestScore=scoreOf(cellBest); rowBestPos=j; rowBestWord=cellBest;
				}
			}

			loops.addAndGet(mloops);
			// Shadow contraction: contract from each side when edge is bad
			if(shadowActive){
				// Contract left
				while(shadowLeft<desiredStart){
					long edgeScore=scoreOf(currM[shadowLeft]);
					if(edgeScore<rowBestScore-CONTRACTION_THRESHOLD){
						shadowLeft++;
					}else{break;}
				}
				// Contract right
				while(shadowRight>desiredEnd){
					long edgeScore=scoreOf(currM[shadowRight]);
					if(edgeScore<rowBestScore-CONTRACTION_THRESHOLD){
						shadowRight--;
					}else{break;}
				}
				// Check if shadow has contracted to desired band
				if(shadowLeft>=desiredStart && shadowRight<=desiredEnd){
					shadowActive=false;
				}
			}

			// Update best
			if(scoreOf(rowBestWord)>scoreOf(bestWord) || i==qLen){
				bestScore=rowBestScore; bestWord=rowBestWord;
			}
			bestPos=rowBestPos;
			if(i==qLen){bestScore=rowBestScore; bestWord=rowBestWord;}

			// Swap rows
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
		return 0f; // identity stub
	}

	/*--------------------------------------------------------------*/
	/*----------------          Packing             ----------------*/
	/*--------------------------------------------------------------*/

	private static final int POSITION_BITS=19;
	private static final int TIME_BITS=18;
	private static final int SCORE_SHIFT=POSITION_BITS+TIME_BITS;
	private static final long POSITION_MASK=(1L<<POSITION_BITS)-1;
	private static final long TIME_MASK=(1L<<TIME_BITS)-1;

	private static final long MATCH=100;
	private static final long MATCH_FIRST=70;
	private static final long N_SCORE=0;
	static final long BAD=Long.MIN_VALUE/2;

	static long pack(long score, int timeInState, int start){
		return (score<<SCORE_SHIFT)
			| (((long)timeInState & POSITION_MASK)<<POSITION_BITS)
			| (start & POSITION_MASK);
	}
	static long scoreOf(long w){return w>>SCORE_SHIFT;}
	static int timeOf(long w){return (int)((w>>>POSITION_BITS) & TIME_MASK);}
	static int startOf(long w){return (int)(w & POSITION_MASK);}

}
