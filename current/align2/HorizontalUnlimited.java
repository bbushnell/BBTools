package align2;

import java.util.Arrays;
import dna.AminoAcid;
import shared.Shared;
import static align2.MultiStateAligner11ts.*;

/** Limited and unlimited legacy scoring with six packed rows and one ancestry byte per cell.
 * @author Collei */
class HorizontalUnlimited {
	HorizontalUnlimited(boolean vector_){vector=vector_ && Shared.SIMD && simd.Vector.simd256;}
	int[] fill(byte[] read,byte[] ref,int start,int stop){
		prepare(read,ref,start,stop);
		return fillPrepared();
	}
	private void prepare(byte[] read,byte[] ref,int start,int stop){
		query=read;reference=ref;refStart=start;refStop=stop;rows=read.length;columns=stop-start+1;stride=columns+1;
		assert(rows>0 && rows<POINTSoff_INS_ARRAY.length && start>=0 && stop<ref.length && columns>0) : "Legacy insertion table and reference bounds must cover the requested fill";
		if(m0==null || m0.length<stride){
			m0=new int[stride];m1=new int[stride];d0=new int[stride];d1=new int[stride];i0=new int[stride];i1=new int[stride];
			refCodes=new int[stride];traceFlags=new int[stride];
		}
		final int cells=Math.multiplyExact(rows+1,stride);
		if(trace==null || trace.length<cells){trace=new byte[cells];}
		Arrays.fill(m0,0,stride,0);Arrays.fill(d0,0,stride,0);Arrays.fill(i0,0,stride,0);
		refCodes[0]='!';for(int col=1;col<=columns;col++){refCodes[col]=ref[start+col-1];}
	}
	private int[] fillPrepared(){
		final byte[] read=query;
		int[] pm=m0,cm=m1,pd=d0,cd=d1,pi=i0,ci=i1;
		final int subfloor=-2*((rows-1)*POINTSoff_MATCH2+POINTSoff_MATCH);
		assert(subfloor>BADoff && subfloor*2L>BADoff) : "Same score-range bound as legacy fillUnlimited";
		int boundary=0;
		for(int row=1;row<=rows;row++){
			boundary+=POINTSoff_INS_ARRAY[row];cm[0]=cd[0]=ci[0]=boundary;
			final int q=read[row-1],q0=row<2?'?':read[row-2];assert(q!=MSA.GAPC) : "Query cannot contain reference-gap compression symbols";
			int firstScalar=1;
			if(vector){firstScalar=simd.SIMDMSA.fillRow(pm,pd,pi,cm,ci,refCodes,traceFlags,q,q0,row,rows,columns,subfloor,BARRIER_I1,MSA.GAPC);}
			for(int col=firstScalar;col<=columns;col++){fillCell(pm,pd,pi,cm,ci,refCodes,traceFlags,q,q0,row,rows,columns,subfloor,col);}
			final boolean blocked=row<BARRIER_D1 || row>rows-BARRIER_D1;
			int left=cd[0];final int offset=row*stride;
			for(int col=1;col<=columns;col++){
				final int ms=cm[col-1]&SCOREMASK,del=left&SCOREMASK,streak=left&TIMEMASK;
				int value=subfloor;
				if(!blocked){
					final int extra=refCodes[col]=='N'?POINTSoff_DEL_REF_N:refCodes[col]==MSA.GAPC?POINTSoff_GAP:0;
					final int a=ms+POINTSoff_DEL+extra,b=del+DEL_COST[streak]+extra;
					final boolean fromMS=a>=b;int time=fromMS?1:streak+1;if(time>MAX_TIME){time=MAX_TIME-MASK5;}
					value=(fromMS?a:b)|time;
				}
				cd[col]=left=value;
				// Match legacy traceback's raw-score tie resolution,not just the fill's winning transition.
				final int prevD=(value&TIMEMASK)>1 || del>ms ? 4:0;
				trace[offset+col]=(byte)(traceFlags[col]|prevD);
			}
			int[] swap=pm;pm=cm;cm=swap;swap=pd;pd=cd;cd=swap;swap=pi;pi=ci;ci=swap;
		}
		lastMS=pm;lastD=pd;lastI=pi;
		int best=Integer.MIN_VALUE,bestCol=-1,bestState=-1;
		for(int state=0;state<3;state++){
			final int[] values=state==0?lastMS:state==1?lastD:lastI;
			for(int col=1;col<=columns;col++){int value=values[col]&SCOREMASK;if(value>best){best=value;bestCol=col;bestState=state;}}
		}
		return new int[]{rows,bestCol,bestState,best>>SCOREOFFSET};
	}

	/** Lenient limited fill: optimistic remaining matches,live spans,and compact ancestry. */
	int[] fillLimited(byte[] read,byte[] ref,int start,int stop,int minimum,int halfband){
		prepare(read,ref,start,stop);lastFillCells=0;
		final int target=minimum<<SCOREOFFSET;
		final int floor=-2*((rows-1)*POINTSoff_MATCH2+POINTSoff_MATCH);
		assert(target>floor && target<MAXoff_SCORE) : "Limited target must fit the legacy score packing";
		int[] pm=m0,cm=m1,pd=d0,cd=d1,pi=i0,ci=i1;
		int previousLo=1,previousHi=columns,boundary=0;
		for(int row=1;row<=rows;row++){
			boundary+=POINTSoff_INS_ARRAY[row];cm[0]=cd[0]=ci[0]=boundary;
			final int from=Math.max(previousLo,halfband>0?Math.max(1,row-halfband):1);
			final int to=Math.min(columns,Math.min(previousHi+1,halfband>0?row+2*halfband:columns));
			if(from>to){return null;}
			if(from>1){cm[from-1]=cd[from-1]=ci[from-1]=floor;}
			final int q=read[row-1],q0=row<2?'?':read[row-2];
			int first=from;
			if(vector){first=simd.SIMDMSALimited.fillRow(pm,pd,pi,cm,ci,refCodes,traceFlags,q,q0,row,rows,columns,floor,BARRIER_I1,MSA.GAPC,from,to,target);}
			for(int col=first;col<=to;col++){
				fillCell(pm,pd,pi,cm,ci,refCodes,traceFlags,q,q0,row,rows,columns,floor,col);
				final int cutoff=target-Math.min(rows-row,columns-col)*POINTSoff_MATCH2;
				if(cm[col]<cutoff){cm[col]=floor;}if(ci[col]<cutoff){ci[col]=floor;}
			}
			int nextLo=-1,nextHi=-1,left=cd[from-1];
			final boolean blocked=row<BARRIER_D1 || row>rows-BARRIER_D1;
			final int offset=row*stride;
			for(int col=from;col<=columns;col++){
				lastFillCells++;
				if(col>to){cm[col]=ci[col]=floor;traceFlags[col]=0;}
				final int ms=cm[col-1]&SCOREMASK,del=left&SCOREMASK,s=left&TIMEMASK;
				final int cutoff=target-Math.min(rows-row,columns-col)*POINTSoff_MATCH2;
				int value=floor;
				if(!blocked){
					final int extra=refCodes[col]=='N'?POINTSoff_DEL_REF_N:refCodes[col]==MSA.GAPC?POINTSoff_GAP:0;
					int x=ms+POINTSoff_DEL+extra,y=del+DEL_COST[s]+extra,t=x>=y?1:s+1;
					if(t>MAX_TIME){t=MAX_TIME-MASK5;}value=Math.max(x,y)|t;if(value<cutoff){value=floor;}
				}
				cd[col]=left=value;
				trace[offset+col]=(byte)(traceFlags[col]|((value&TIMEMASK)>1 || del>ms?4:0));
				final boolean live=cm[col]>=cutoff || ci[col]>=cutoff || value>=cutoff;
				if(live){if(nextLo<0){nextLo=col;}nextHi=col;}
				if(col>=to && halfband>0){break;}
				if(col>to && !live){break;}
			}
			// Preserve a possible leading-insertion path on column0 even when no interior cell survives.
			if(boundary>=target-Math.min(rows-row,columns)*POINTSoff_MATCH2 && (halfband<1 || row+1-halfband<=1)){
				nextLo=1;nextHi=Math.max(0,nextHi);
			}
			if(nextLo<0){return null;}
			if(nextHi<columns){cm[nextHi+1]=cd[nextHi+1]=ci[nextHi+1]=floor;}
			previousLo=nextLo;previousHi=nextHi;
			int[] swap=pm;pm=cm;cm=swap;swap=pd;pd=cd;cd=swap;swap=pi;pi=ci;ci=swap;
		}
		lastMS=pm;lastD=pd;lastI=pi;
		Arrays.fill(lastMS,1,previousLo,floor);Arrays.fill(lastMS,previousHi+1,columns+1,floor);
		Arrays.fill(lastD,1,previousLo,floor);Arrays.fill(lastD,previousHi+1,columns+1,floor);
		Arrays.fill(lastI,1,previousLo,floor);Arrays.fill(lastI,previousHi+1,columns+1,floor);
		int best=Integer.MIN_VALUE,bestCol=-1,bestState=-1;
		for(int state=0;state<3;state++){
			int[] values=state==0?lastMS:state==1?lastD:lastI;
			for(int col=previousLo;col<=previousHi;col++){int value=values[col]&SCOREMASK;if(value>best){best=value;bestCol=col;bestState=state;}}
		}
		return best<target?null:new int[]{rows,bestCol,bestState,best>>SCOREOFFSET};
	}

	static void fillCell(int[] pm,int[] pd,int[] pi,int[] cm,int[] ci,int[] refs,int[] flags,int q,int q0,int row,int rows,int columns,int floor,int col){
		final int r=refs[col],r0=refs[col-1];final boolean gap=r==MSA.GAPC,match=q==r && r!='N',previousMatch=q0==r0 && r0!='N';
		final int a=pm[col-1]&SCOREMASK,b=pd[col-1]&SCOREMASK,c=pi[col-1]&SCOREMASK,streak=pm[col-1]&TIMEMASK;
		int ms=floor;
		if(!gap){
			final int gain=match?(previousMatch?POINTSoff_MATCH2:POINTSoff_MATCH):q=='N'||r=='N'?POINTSoff_NOCALL:
				previousMatch?(streak<=1?POINTSoff_SUBR:POINTSoff_SUB):POINTSoff_SUB_ARRAY[streak+1];
			final int x=a+gain,y=b+(match?POINTSoff_MATCH:POINTSoff_SUB),z=c+(match?POINTSoff_MATCH:POINTSoff_SUB);
			final boolean fromMS=x>=y && x>=z;int time=fromMS?(match==previousMatch?streak+1:1):1;
			if(time>MAX_TIME){time=MAX_TIME-MASK5;}ms=Math.max(x,Math.max(y,z))|time;
		}
		cm[col]=ms;
		final int prevMS=(ms&TIMEMASK)>1?0:a>=b && a>=c?0:b>=c?1:2;
		final int upMS=pm[col]&SCOREMASK,upI=pi[col]&SCOREMASK,is=pi[col]&TIMEMASK;
		int ins=floor;
		if(!gap && !(row<BARRIER_I1 && col>1) && !(row>rows-BARRIER_I1 && col<columns-1)){
			final int x=upMS+POINTSoff_INS,y=upI+POINTSoff_INS_ARRAY[is+1];
			int time=x>=y?1:is+1;if(time>MAX_TIME){time=MAX_TIME-MASK5;}ins=Math.max(x,y)|time;
		}
		ci[col]=ins;flags[col]=prevMS|((ins&TIMEMASK)>1 || upI>upMS?8:0);
	}

	// Compressed fills include a cushion column (greflimit), while score() supplies
	// the actual compressed endpoint (greflimit-1). The retained rows cover both.
	boolean matches(byte[] read,byte[] ref,int start,int stop,int row){return query==read && reference==ref && refStart==start && stop<=refStop && stop>=start && row==rows;}
	int[] score(int maxRow,int maxCol,int maxState,int requestedStop){
		assert(maxRow==rows && maxCol>=0 && maxCol<=columns && maxState>=0 && maxState<3) : "Rolling score queries require the retained last row and a valid state/column";
		int row=maxRow,col=maxCol,state=maxState,stateTime=0;
		final int score=(maxState==0?lastMS:maxState==1?lastD:lastI)[maxCol]&SCOREMASK;
		final int bestStop=refStart+col-1;
		while(row>0 && col>0){
			final int prev=previous(row,col,state);
			if(state!=MSA.MODE_DEL){row--;}if(state!=MSA.MODE_INS){col--;}
			stateTime=state==prev?stateTime+1:0;state=prev;
		}
		if(row>col){col-=row;}final int bestStart=refStart+col;
		final int padLeft=bestStart<refStart?refStart-bestStart:bestStart==refStart && state==MSA.MODE_INS?stateTime:0;
		final int padRight=bestStop>requestedStop?bestStop-requestedStop:bestStop==requestedStop && maxState==MSA.MODE_INS?lastI[maxCol]&TIMEMASK:0;
		return padLeft>0 || padRight>0?new int[]{score>>SCOREOFFSET,bestStart,bestStop,maxRow,maxCol,maxState,padLeft,padRight}:
			new int[]{score>>SCOREOFFSET,bestStart,bestStop,maxRow,maxCol,maxState};
	}
	byte[] traceback(int row,int col,int state){
		assert(row==rows && col>0 && col<=columns && state>=0 && state<3) : "Traceback must address this fill's last row";
		byte[] reverse=new byte[row+col];int length=0,gaps=0;
		while(row>0 && col>0){
			final int prev=previous(row,col,state);byte op;
			if(state==MSA.MODE_MS){
				final byte q=query[row-1],r=reference[refStart+col-1];
				op=q==r && r!='N'?(byte)'m':!AminoAcid.isFullyDefined(q)||!AminoAcid.isFullyDefined(r)?(byte)'N':(byte)'S';row--;col--;
			}else if(state==MSA.MODE_DEL){op=reference[refStart+col-1]==MSA.GAPC?MSA.GAPC:(byte)'D';if(op==MSA.GAPC){gaps++;}col--;}
			else{op=col>=columns?(byte)'Y':(byte)'I';row--;}
			reverse[length++]=op;state=prev;
		}
		while(row-->0){reverse[length++]='X';}
		final byte[] result=new byte[length+gaps*(MSA.GAPLEN-1)];int next=0;
		for(int p=length-1;p>=0;p--){byte c=reverse[p];if(c==MSA.GAPC){Arrays.fill(result,next,next+MSA.GAPLEN,(byte)'D');next+=MSA.GAPLEN;}else{result[next++]=c;}}
		assert(next==result.length) : "Every compressed deletion expands to GAPLEN bases";return result;
	}
	int traceIndex(int row,int col){return row*stride+col;}
	private int previous(int row,int col,int state){int f=trace[traceIndex(row,col)];return state==0?f&3:state==1?((f&4)==0?0:1):((f&8)==0?0:2);}
	private static final int[] DEL_COST=new int[MAX_TIME+1];
	static{for(int s=0;s<DEL_COST.length;s++){DEL_COST[s]=s==0?POINTSoff_DEL:s<LIMIT_FOR_COST_3?POINTSoff_DEL2:s<LIMIT_FOR_COST_4?POINTSoff_DEL3:s<LIMIT_FOR_COST_5?POINTSoff_DEL4:(s&MASK5)==0?POINTSoff_DEL5:0;}}
	final boolean vector;
	long lastFillCells;
	byte[] query,reference,trace;int refStart,refStop,rows,columns,stride;
	private int[] m0,m1,d0,d1,i0,i1,refCodes,traceFlags;
	int[] lastMS,lastD,lastI;
}
