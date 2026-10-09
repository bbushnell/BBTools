package align2;

import java.util.Arrays;
import dna.AminoAcid;
import static align2.MultiStateAligner9PacBio.*;

/** Unlimited PacBio scoring with six packed score rows and one ancestry byte per cell.
 * One worker owns this object. Consume score/trace before the next fill; inputs
 * remain borrowed until then. Retains native raw-predecessor traceback ties,
 * boundary insertion tiers, N/N matches, inclusive endpoints and padding.
 * No compressed-reference construction or limited-fill policy is provided.
 * @author Brian Bushnell, Raiden
 */
public final class PacBioRollingAligner {

	public PacBioRollingAligner(){this(PacBioScoreParameters.DEFAULT);}
	public PacBioRollingAligner(PacBioScoreParameters costs_){
		if(costs_==null){throw new IllegalArgumentException("Rolling PacBio requires immutable instance penalties");}
		costs=costs_;
	}

	/** Returns native {query length, ending column, state, score}. */
	public int[] fillUnlimited(byte[] query_, byte[] reference_, int start, int stop){
		ready=false;
		if(query_==null || reference_==null || query_.length<1 || start<0 || stop<start || stop>=reference_.length){
			throw new IllegalArgumentException("Unlimited PacBio needs a nonempty query and inclusive reference window");
		}
		query=query_;reference=reference_;refStart=start;refStop=stop;
		rows=query.length;columns=stop-start+1;stride=Math.addExact(columns, 1);
		final long gain=(rows-1L)*POINTSoff_MATCH2+POINTSoff_MATCH;
		if(-4L*gain<=BADoff){throw new IllegalArgumentException("Query exceeds native PacBio packed-score range: "+rows);}
		final int floor=(int)(-2L*gain);
		final int cells=Math.multiplyExact(Math.addExact(rows, 1), stride);
		final boolean allocate=m0==null || m0.length<stride || ancestry==null || ancestry.length<cells;
		final long allocationBegin=allocate ? System.nanoTime() : 0;
		if(m0==null || m0.length<stride){
			m0=new int[stride];m1=new int[stride];d0=new int[stride];d1=new int[stride];i0=new int[stride];i1=new int[stride];
		}
		if(ancestry==null || ancestry.length<cells){ancestry=new byte[cells];}
		allocationNanos=allocate ? System.nanoTime()-allocationBegin : 0;
		Arrays.fill(m0, 0, stride, 0);Arrays.fill(d0, 0, stride, 0);Arrays.fill(i0, 0, stride, 0);
		int[] pm=m0, cm=m1, pd=d0, cd=d1, pi=i0, ci=i1;
		int boundary=0;
		for(int row=1; row<=rows; row++){
			// Native constructor uses row<5 / row<20, one residue before internal tiers.
			boundary+=row==1 ? costs.offIns : row<LIMIT_FOR_COST_3 ? POINTSoff_INS2
				: row<LIMIT_FOR_COST_4 ? POINTSoff_INS3 : costs.offIns4;
			cm[0]=cd[0]=ci[0]=boundary;
			final byte q=query[row-1], q0=row<2 ? (byte)'?' : query[row-2];
			assert(q!=MSA.GAPC) : "Native PacBio query cannot contain reference gap-compression symbols";
			final int offset=row*stride;
			for(int col=1; col<=columns; col++){
				final byte r=reference[start+col-1], r0=col<2 ? (byte)'!' : reference[start+col-2];
				final boolean match=q==r && r!='N', previousMatch=q0==r0 && r0!='N', gap=r==MSA.GAPC;
				final int ms=pm[col-1]&SCOREMASK, del=pd[col-1]&SCOREMASK, ins=pi[col-1]&SCOREMASK;
				final int streak=pm[col-1]&TIMEMASK;
				int value=floor;
				if(!gap){
					final int gainMS=match ? (previousMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH)
						: r=='N' || q=='N' ? POINTSoff_NOCALL
						: previousMatch ? (streak<=1 ? costs.offSubR : costs.offSub)
						: streak==0 ? costs.offSub : streak<LIMIT_FOR_COST_3 ? costs.offSub2 : costs.offSub3;
					final int a=ms+gainMS, b=del+(match ? POINTSoff_MATCH : costs.offSub), c=ins+(match ? POINTSoff_MATCH : costs.offSub);
					final boolean fromMS=a>=b && a>=c;
					value=pack(fromMS ? a : Math.max(b, c), fromMS && match==previousMatch ? streak+1 : 1);
				}
				cm[col]=value;
				// Native traceback compares raw predecessor scores when time<=1, not weighted transitions.
				int flags=(value&TIMEMASK)>1 ? 0 : ms>=del && ms>=ins ? 0 : del>=ins ? 1 : 2;
				final int leftMS=cm[col-1]&SCOREMASK, leftD=cd[col-1]&SCOREMASK, ds=cd[col-1]&TIMEMASK;
				value=floor;
				if(row>=BARRIER_D && row<=rows-BARRIER_D){
					final int extra=r=='N' ? POINTSoff_DEL_REF_N : gap ? POINTSoff_GAP : 0;
					final int a=leftMS+costs.offDel+extra;
					final int b=leftD+extra+(ds==0 ? costs.offDel : ds<LIMIT_FOR_COST_3 ? POINTSoff_DEL2
						: ds<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 : ds<LIMIT_FOR_COST_5 ? costs.offDel4 : (ds&MASK5)==0 ? costs.offDel5 : 0);
					value=pack(Math.max(a, b), a>=b ? 1 : ds+1);
				}
				cd[col]=value;
				if((value&TIMEMASK)>1 || leftD>leftMS){flags|=4;}
				final int upMS=pm[col]&SCOREMASK, upI=pi[col]&SCOREMASK, is=pi[col]&TIMEMASK;
				value=floor;
				if(!gap && !(row<BARRIER_I && col>1) && !(row>rows-BARRIER_I && col<columns-1)){
					final int a=upMS+costs.offIns;
					final int b=upI+(is==0 ? costs.offIns : is<LIMIT_FOR_COST_3 ? POINTSoff_INS2
						: is<LIMIT_FOR_COST_4 ? POINTSoff_INS3 : costs.offIns4);
					value=pack(Math.max(a, b), a>=b ? 1 : is+1);
				}
				ci[col]=value;
				if((value&TIMEMASK)>1 || upI>upMS){flags|=8;}
				ancestry[offset+col]=(byte)flags;
			}
			int[] swap=pm;pm=cm;cm=swap;swap=pd;pd=cd;cd=swap;swap=pi;pi=ci;ci=swap;
		}
		bestScore=Integer.MIN_VALUE;bestCol=bestState=-1;
		for(int state=0; state<3; state++){
			final int[] last=state==0 ? pm : state==1 ? pd : pi;
			for(int col=1; col<=columns; col++){
				final int value=last[col]&SCOREMASK;
				if(value>bestScore){bestScore=value;bestCol=col;bestState=state;}
			}
		}
		lastInsTime=pi[bestCol]&TIMEMASK;bestScore>>=SCOREOFFSET;ready=true;
		return new int[]{rows, bestCol, bestState, bestScore};
	}

	private static int pack(int score, int time){
		assert(score>=MINoff_SCORE && score<=MAXoff_SCORE && (score&SCOREMASK)==score) :
			"Native PacBio packed cells require shifted scores in range; score="+score;
		if(time>MAX_TIME){time=MAX_TIME-MASK5;}
		return score|time;
	}
	/** Same six/eight-field score vector as native score2, including suggested padding. */
	public int[] score(){
		assert(ready) : "Score consumes the latest successful fill, never an interrupted fill";
		int row=rows, col=bestCol, state=bestState, stateTime=0;
		final int stop=refStart+col-1;
		while(row>0 && col>0){
			final int prev=previous(row, col, state);
			if(state!=MSA.MODE_DEL){row--;}if(state!=MSA.MODE_INS){col--;}
			stateTime=state==prev ? stateTime+1 : 0;state=prev;
		}
		if(row>col){col-=row;}
		final int start=refStart+col;
		final int left=start<refStart ? refStart-start : start==refStart && state==MSA.MODE_INS ? stateTime : 0;
		final int right=stop>refStop ? stop-refStop : stop==refStop && bestState==MSA.MODE_INS ? lastInsTime : 0;
		return left>0 || right>0 ? new int[]{bestScore, start, stop, rows, bestCol, bestState, left, right}
			: new int[]{bestScore, start, stop, rows, bestCol, bestState};
	}
	/** Native expanded m/S/N/D/I/X/Y operations, including equal ambiguous bases as m. */
	public byte[] traceback(){
		assert(ready) : "Traceback consumes retained ancestry from a successful fill";
		int row=rows, col=bestCol, state=bestState, length=0, gaps=0;
		final byte[] reverse=new byte[Math.addExact(row, col)];
		while(row>0 && col>0){
			final int prev=previous(row, col, state);
			final byte op;
			if(state==MSA.MODE_MS){
				final byte q=query[row-1], r=reference[refStart+col-1];
				op=q==r ? (byte)'m' : !AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r) ? (byte)'N' : (byte)'S';
				row--;col--;
			}else if(state==MSA.MODE_DEL){
				op=reference[refStart+col-1]==MSA.GAPC ? MSA.GAPC : (byte)'D';
				if(op==MSA.GAPC){gaps++;}col--;
			}else{op=col>=columns ? (byte)'Y' : (byte)'I';row--;}
			reverse[length++]=op;state=prev;
		}
		while(row-->0){reverse[length++]='X';}
		final byte[] result=new byte[Math.addExact(length, Math.multiplyExact(gaps, MSA.GAPLEN-1))];
		int next=0;
		for(int p=length-1; p>=0; p--){
			final byte op=reverse[p];
			if(op==MSA.GAPC){Arrays.fill(result, next, next+MSA.GAPLEN, (byte)'D');next+=MSA.GAPLEN;}
			else{result[next++]=op;}
		}
		assert(next==result.length) : "Every compressed-reference gap expands to exactly GAPLEN operations";
		return result;
	}
	private int previous(int row, int col, int state){
		final int flags=ancestry[row*stride+col];
		return state==0 ? flags&3 : state==1 ? (flags&4)==0 ? 0 : 1 : (flags&8)==0 ? 0 : 2;
	}
	/** Primitive DP payload only: six int rows and byte ancestry, excluding headers/trace output. */
	public long workingBytes(){assert(ready) : "Working dimensions require a completed fill";return (rows+1L)*stride+24L*stride;}
	public long retainedBytes(){return ancestry==null ? 0 : ancestry.length+24L*m0.length;}
	public long allocationNanos(){return allocationNanos;}
	private final PacBioScoreParameters costs;
	private int[] m0, m1, d0, d1, i0, i1;
	private byte[] ancestry, query, reference;
	private int rows, columns, stride, refStart, refStop, bestCol, bestState, bestScore, lastInsTime;
	private boolean ready;
	private long allocationNanos;
	/** Native MultiStateAligner9PacBio BARRIER_I1/BARRIER_D1; unchanged by cost experiments. */
	private static final int BARRIER_I=1, BARRIER_D=1;
}
