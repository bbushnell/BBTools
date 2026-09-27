package align2;

import java.util.Arrays;
import static align2.MultiStateAligner11ts.*;

/** Experimental unlimited fill with three rotating antidiagonals per state.
 * Slower than horizontal SIMD on the tested real HG001 workload:154.931s versus98.084s.
 * Retained for experiments;never selected automatically. Limited calls use horizontal SIMD.
 * @author Collei */
final class DiagonalUnlimited extends HorizontalUnlimited {
	DiagonalUnlimited(){super(true);}
	@Override int[] fillLimited(byte[] read,byte[] ref,int start,int stop,int minimum,int halfband){
		diagonalTrace=false;
		return super.fillLimited(read,ref,start,stop,minimum,halfband);
	}
	@Override int[] fill(byte[] read,byte[] ref,int start,int stop){
		diagonalTrace=true;
		query=read;reference=ref;refStart=start;refStop=stop;rows=read.length;columns=stop-start+1;stride=columns+1;
		assert(rows>0 && rows<POINTSoff_INS_ARRAY.length && start>=0 && stop<ref.length && columns>0) : "Legacy insertion table and reference bounds must cover diagonal fill";
		if(scores==null || scores[0].length<stride){
			scores=new int[9][stride];refs=new int[stride];diagonalMS=new int[stride];diagonalD=new int[stride];diagonalI=new int[stride];
		}
		// A preceding horizontal limited fill installs its own last-row arrays.
		lastMS=diagonalMS;lastD=diagonalD;lastI=diagonalI;
		if(reverse==null || reverse.length<rows+2){reverse=new int[rows+2];boundary=new int[rows+1];}
		if(offsets==null || offsets.length<rows+columns+1){offsets=new int[rows+columns+1];}
		final int cells=Math.multiplyExact(rows,columns);
		if(trace==null || trace.length<cells){trace=new byte[cells];}
		for(int i=0;i<rows;i++){reverse[i]=read[rows-1-i];assert(reverse[i]!=MSA.GAPC) : "Query cannot contain reference-gap compression symbols";}
		reverse[rows]='?';refs[0]='!';for(int col=1;col<=columns;col++){refs[col]=ref[start+col-1];}
		boundary[0]=0;for(int row=1;row<=rows;row++){boundary[row]=boundary[row-1]+POINTSoff_INS_ARRAY[row];}
		for(int[] score:scores){Arrays.fill(score,0,stride,0);}
		int[] olderM=scores[0],olderD=scores[1],olderI=scores[2],prevM=scores[3],prevD=scores[4],prevI=scores[5],curM=scores[6],curD=scores[7],curI=scores[8];
		// Diagonal1 contains only the two boundary cells; diagonal0 is the origin.
		prevM[0]=prevD[0]=prevI[0]=boundary[1];prevM[1]=prevD[1]=prevI[1]=0;
		final int floor=-2*((rows-1)*POINTSoff_MATCH2+POINTSoff_MATCH);
		int next=0;
		for(int diagonal=2;diagonal<=rows+columns;diagonal++){
			final int lo=Math.max(1,diagonal-rows),hi=Math.min(columns,diagonal-1);
			offsets[diagonal]=next-lo;next+=hi-lo+1;
			if(diagonal<=rows){curM[0]=curD[0]=curI[0]=boundary[diagonal];}
			if(diagonal<=columns){curM[diagonal]=curD[diagonal]=curI[diagonal]=0;}
			simd.SIMDMSADiagonal.fill(olderM,olderD,olderI,prevM,prevD,prevI,curM,curD,curI,refs,reverse,trace,offsets[diagonal],diagonal,rows,columns,lo,hi,floor,BARRIER_I1,BARRIER_D1,MSA.GAPC);
			if(diagonal>rows){final int col=diagonal-rows;lastMS[col]=curM[col];lastD[col]=curD[col];lastI[col]=curI[col];}
			int[] swap=olderM;olderM=prevM;prevM=curM;curM=swap;swap=olderD;olderD=prevD;prevD=curD;curD=swap;swap=olderI;olderI=prevI;prevI=curI;curI=swap;
		}
		assert(next==cells) : "Antidiagonals must partition the complete row-column grid";
		int best=Integer.MIN_VALUE,bestCol=-1,bestState=-1;
		for(int state=0;state<3;state++){
			final int[] values=state==0?lastMS:state==1?lastD:lastI;
			for(int col=1;col<=columns;col++){int value=values[col]&SCOREMASK;if(value>best){best=value;bestCol=col;bestState=state;}}
		}
		return new int[]{rows,bestCol,bestState,best>>SCOREOFFSET};
	}
	@Override int traceIndex(int row,int col){return diagonalTrace?offsets[row+col]+col:super.traceIndex(row,col);}
	private boolean diagonalTrace;
	private int[] diagonalMS,diagonalD,diagonalI;
	private int[][] scores;
	private int[] refs,reverse,boundary,offsets;
}
