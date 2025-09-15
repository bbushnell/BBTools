package aligner;

import java.util.Arrays;

import shared.KillSwitch;
import shared.Tools;

/**
 * Based on MSA9PBA, but reduced to a single matrix. */
public final class SingleStateAlignerPacBioAdapter {
	
	
	public SingleStateAlignerPacBioAdapter(int maxRows_, int maxColumns_, int qlen){
//		assert(maxColumns_>=200);
//		assert(maxRows_>=200);
		maxRows=maxRows_;
		maxColumns=maxColumns_;
		packed=new int[maxRows+1][maxColumns+1];

		insScoreArray=KillSwitch.allocInt1D(5);
		delScoreArray=KillSwitch.allocInt1D(5);
		matchScoreArray=KillSwitch.allocInt1D(5);
		subScoreArray=KillSwitch.allocInt1D(5);

		Arrays.fill(insScoreArray, POINTSoff_INS|MODE_INS);
		Arrays.fill(delScoreArray, POINTSoff_DEL|MODE_DEL);
		Arrays.fill(matchScoreArray, POINTSoff_MATCH|MODE_MATCH);
		Arrays.fill(subScoreArray, POINTSoff_SUB|MODE_SUB);
		
		insScoreArray[MODE_INS]=POINTSoff_INS2|MODE_INS;
		delScoreArray[MODE_DEL]=POINTSoff_DEL2|MODE_DEL;
		matchScoreArray[MODE_MATCH]=POINTSoff_MATCH2|MODE_MATCH;
		subScoreArray[MODE_SUB]=POINTSoff_SUB2|MODE_SUB;
		
		vertLimit=new int[maxRows+1];
		horizLimit=new int[maxColumns+1];
		Arrays.fill(vertLimit, BADoff);
		Arrays.fill(horizLimit, BADoff);
		
//		for(int i=0; i<maxColumns+1; i++){
//			scores[0][i]=0-i;
//		}
		
		for(int i=1; i<=maxRows; i++){
			for(int j=1; j<packed[i].length; j++){
				packed[i][j]|=BADoff;
			}
		}
		for(int i=0; i<=maxRows; i++){

			int prevScore=(i<2 ? 0 : packed[i-1][0]);
			int score=(i<2 ? (i*POINTSoff_INS) : prevScore+POINTSoff_INS2);

			packed[i][0]=score;
		}
		for(int i=0; i<packed[0].length; i++){
			packed[0][i]=packed[0][i]|(i<<STARTOFFSET);
			int x=packed[0].length-i;
			int qbases=qlen-x;
			if(qbases>0){
				packed[0][i]|=calcDelScoreOffset(qbases); //Forces consumption of query
			}
		}
	}
	
	
	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Will not fill areas that cannot match minScore */
	public final int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		return fillLimitedX(read, ref, refStartLoc, refEndLoc, minScore);
	}
	
	
	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Will not fill areas that cannot match minScore */
	private final int[] fillLimitedX(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		return fillUnlimited(read, ref, refStartLoc, refEndLoc, minScore);
	}
	
	
	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Does not require a min score (ie, same as old method) */
	private final int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, final int minScore){
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;
		
		//temporary, for finding a bug
		if(rows>maxRows || columns>maxColumns){
			throw new RuntimeException("rows="+rows+", maxRows="+maxRows+", cols="+columns+", maxCols="+maxColumns+"\n"+new String(read)+"\n");
		}
		
		assert(rows<=maxRows) : "Check that values are in-bounds before calling this function: "+rows+", "+maxRows;
		assert(columns<=maxColumns) : "Check that values are in-bounds before calling this function: "+columns+", "+maxColumns;
		
		assert(refStartLoc>=0) : "Check that values are in-bounds before calling this function: "+refStartLoc;
		assert(refEndLoc<ref.length) : "Check that values are in-bounds before calling this function: "+refEndLoc+", "+ref.length;

		final int refOffset=refStartLoc-1;
		for(int row=1; row<=rows; row++){

			final byte qBase=read[row-1];
			for(int col=1; col<=columns; col++){
				iterationsUnlimited++;
				
				final byte rBase=ref[refOffset+col];
				
				final boolean match=(qBase==rBase);//Note that qBase will never be 'N'

				final int valueFromDiag=packed[row-1][col-1];
				final int valueFromDel=packed[row][col-1];
				final int valueFromIns=packed[row-1][col];
				final int stateFromDiag=valueFromDiag&MODEMASK;
				final int scoreFromDiag=valueFromDiag&HIGHMASK;
				final int stateFromDel=valueFromDel&MODEMASK;
				final int scoreFromDel=valueFromDel&HIGHMASK;
				final int stateFromIns=valueFromIns&MODEMASK;
				final int scoreFromIns=valueFromIns&HIGHMASK;
				
//				final boolean prevMatch=stateFromDiag==MODE_MATCH;
//				final boolean prevSub=stateFromDiag==MODE_SUB;
//				final boolean prevDel=stateFromDel==MODE_DEL;
//				final boolean prevIns=stateFromIns==MODE_INS;
				
				//Old conditional code, replaced by faster arrays
//				final int diagScoreM=MODE_MATCH|(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
//				final int diagScoreS=MODE_SUB|(prevSub ? POINTSoff_SUB2 : POINTSoff_SUB);
//				final int delScore=scoreFromDel+(prevDel ? POINTSoff_DEL2 : POINTSoff_DEL)|MODE_DEL;
//				final int insScore=scoreFromIns+(prevIns ? POINTSoff_INS2 : POINTSoff_INS)|MODE_INS;
				
				final int diagScoreM=matchScoreArray[stateFromDiag];
				final int diagScoreS=subScoreArray[stateFromDiag];
				final int delScore=scoreFromDel+delScoreArray[stateFromDel];
				final int insScore=scoreFromIns+insScoreArray[stateFromIns];
				
				int diagScore=(match ? diagScoreM : diagScoreS);
				diagScore=scoreFromDiag+(rBase!='N' ? diagScore : POINTSoff_NOREF_MODE_SUB);
				
				int score=diagScore>=delScore ? diagScore : delScore;
				score=score>=insScore ? score : insScore;
				
				packed[row][col]=score;
			}
		}
		

		int maxCol=-1;
		int maxState=-1;
		int maxStart=-1;
		int maxScore=Integer.MIN_VALUE;
		
		for(int col=1; col<=columns; col++){
			int x=packed[rows][col];
			if((x&SCOREMASK)>maxScore){
				maxScore=x;
				maxCol=col;
				maxState=x&MODEMASK;
				maxStart=x&STARTMASK;
			}
		}
		maxScore>>=SCOREOFFSET;

//		System.err.println("Returning "+rows+", "+maxCol+", "+maxState+", "+maxScore+"; minScore="+minScore);
		return maxScore<minScore ? null : new int[] {rows, maxCol, maxState, maxScore, maxStart};
	}
	
	
	
	/** Generates the match string */
	public final byte[] traceback(int refStartLoc, int refEndLoc, int row, int col, int state){
//		assert(false);
		assert(refStartLoc<=refEndLoc) : refStartLoc+", "+refEndLoc;
		assert(row==rows);
		
		byte[] out=new byte[row+col-1]; //TODO if an out of bound crash occurs, try removing the "-1".
		int outPos=0;

//		assert(state==(packed[row][col]&MODEMASK));
		
		while(row>0 && col>0){
			state=(packed[row][col]&MODEMASK);
			if(state==MODE_MATCH){
				col--;
				row--;
				out[outPos]='m';
			}else if(state==MODE_SUB){
				col--;
				row--;
				out[outPos]='S';
			}else if(state==MODE_DEL){
				col--;
				out[outPos]='D';
			}else if(state==MODE_INS){
				row--;
				out[outPos]='I';
			}else{
				assert(false) : state;
			}
			outPos++;
		}
		
		assert(row==0 || col==0);
		if(col!=row){
			while(row>0){
				out[outPos]='X';
				outPos++;
				row--;
				col--;
			}
			if(col>0){
				//do nothing
			}
		}
		
		//Shrink and reverse the string
		byte[] out2=new byte[outPos];
		for(int i=0; i<outPos; i++){
			out2[i]=out[outPos-i-1];
		}
		out=null;
		
		return out2;
	}
	
	/** @return {score, bestRefStart, bestRefStop} */
	public final int[] score(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
			final int maxRow, final int maxCol, final int maxState/*, final int maxScore, final int maxStart*/){
		
		int row=maxRow;
		int col=maxCol;
		int state=maxState;

		assert(maxState>=0 && maxState<packed.length) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		assert(maxRow>=0 && maxRow<packed.length) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		assert(maxCol>=0 && maxCol<packed[0].length) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		
		int score=packed[maxRow][maxCol]&SCOREMASK; //Or zero, if it is to be recalculated
		
		if(row<rows){
			int difR=rows-row;
			int difC=columns-col;
			
			while(difR>difC){
				score+=POINTSoff_NOREF;
				difR--;
			}
			
			row+=difR;
			col+=difR;
			
		}
		
		assert(refStartLoc<=refEndLoc);
		assert(row==rows);

		
		final int bestRefStop=refStartLoc+col-1;
		
		while(row>0 && col>0){
			state=(packed[row][col]&MODEMASK);
			if(state==MODE_MATCH){
				col--;
				row--;
			}else if(state==MODE_SUB){
				col--;
				row--;
			}else if(state==MODE_DEL){
				col--;
			}else if(state==MODE_INS){
				row--;
			}else{
				assert(false) : state;
			}
		}
//		assert(false) : row+", "+col;
		if(row>col){
			col-=row;
		}
		
		final int bestRefStart=refStartLoc+col;
		
		score>>=SCOREOFFSET;
//		System.err.println("t2\t"+score+", "+maxScore+", "+maxStart+", "+bestRefStart);
		int[] rvec;
		if(bestRefStart<refStartLoc || bestRefStop>refEndLoc){ //Suggest extra padding in cases of overflow
			int padLeft=Tools.max(0, refStartLoc-bestRefStart);
			int padRight=Tools.max(0, bestRefStop-refEndLoc);
			rvec=new int[] {score, bestRefStart, bestRefStop, padLeft, padRight};
		}else{
			rvec=new int[] {score, bestRefStart, bestRefStop};
		}
		return rvec;
	}
	
	
	/** Will not fill areas that cannot match minScore.
	 * @return {score, bestRefStart, bestRefStop}  */
	public final int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		
		if(b-a>=maxColumns){
			System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
			assert(false) : refStartLoc+", "+refEndLoc;
			b=Tools.min(ref.length-1, a+maxColumns-1);
		}
		int[] max=fillLimited(read, ref, a, b, minScore);
//		return max==null ? null : new int[] {max[3], 0, max[1]};
		
		int[] score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2]/*, max[3], max[4]*/));
		
		return score;
	}
	
	public static final String toString(byte[] ref, int startLoc, int stopLoc){
		StringBuilder sb=new StringBuilder(stopLoc-startLoc+1);
		for(int i=startLoc; i<=stopLoc; i++){sb.append((char)ref[i]);}
		return sb.toString();
	}
	
	public static int calcDelScore(int len){
		if(len<=0){return 0;}
		int score=POINTS_DEL;
		if(len>1){
			score+=(len-1)*POINTS_DEL2;
		}
		return score;
	}
	
	private static int calcDelScoreOffset(int len){
		if(len<=0){return 0;}
		int score=POINTSoff_DEL;
		if(len>1){
			score+=(len-1)*POINTSoff_DEL2;
		}
		return score;
	}
	
	public static int calcInsScore(int len){
		if(len<=0){return 0;}
		int score=POINTS_INS;
		
		if(len>1){
			score+=(len-1)*POINTS_INS2;
		}
		return score;
	}
	
	private static int calcInsScoreOffset(int len){
		if(len<=0){return 0;}
		int score=POINTSoff_INS;
		
		if(len>1){
			score+=(len-1)*POINTSoff_INS2;
		}
		return score;
	}
	
	
	private final int maxRows;
	public final int maxColumns;

	private final int[][] packed;

	private final int[] vertLimit;
	private final int[] horizLimit;

	private final int[] insScoreArray;
	private final int[] delScoreArray;
	private final int[] matchScoreArray;
	private final int[] subScoreArray;

	public static final int MODEBITS=3;
	public static final int STARTBITS=9;
	public static final int LOWBITS=MODEBITS+STARTBITS;
	public static final int SCOREBITS=32-STARTBITS;
	public static final int MAX_START=((1<<STARTBITS)-1);
	public static final int MAX_SCORE=((1<<(SCOREBITS-1))-1)-2000;
	public static final int MIN_SCORE=0-MAX_SCORE; //Keeps it 1 point above "BAD".

//	public static final int MODEOFFSET=0; //Always zero.
	public static final int STARTOFFSET=MODEBITS;
	public static final int SCOREOFFSET=LOWBITS;

	public static final int MODEMASK=~((-1)<<MODEBITS);
	public static final int STARTMASK=(~((-1)<<STARTBITS))<<STARTOFFSET;
	public static final int SCOREMASK=(~((-1)<<SCOREBITS))<<SCOREOFFSET;
	public static final int HIGHMASK=SCOREMASK|STARTMASK;

	//For some reason changing MODE_DEL from 1 to 0 breaks everything
	private static final byte MODE_DEL=1;
	private static final byte MODE_INS=2;
	private static final byte MODE_SUB=3;
	private static final byte MODE_MATCH=4;
	
	public static final int POINTS_NOREF=-10;
	public static final int POINTS_MATCH=90;
	public static final int POINTS_MATCH2=100; //Note:  Changing to 90 substantially reduces false positives
	public static final int POINTS_SUB=-143;
	public static final int POINTS_SUB2=-54;
	public static final int POINTS_INS=-207;
	public static final int POINTS_INS2=-51;
	public static final int POINTS_INS3=-37;
	public static final int POINTS_DEL=-273;
	public static final int POINTS_DEL2=-38;
	
	
//	public static final int POINTS_NOREF=-10;
//	public static final int POINTS_NOCALL=-10;
//	public static final int POINTS_MATCH=90;
//	public static final int POINTS_MATCH2=100; //Note:  Changing to 90 substantially reduces false positives
//	public static final int POINTS_SUB=-143;
//	public static final int POINTS_SUB2=-54;
//	public static final int POINTS_INS=-207;
//	public static final int POINTS_INS2=-51;
//	public static final int POINTS_INS3=-37;
//	public static final int POINTS_DEL=-273;
//	public static final int POINTS_DEL2=-38;
	
	public static final int BAD=MIN_SCORE-1;
	
	
	public static final int POINTSoff_NOREF=(POINTS_NOREF<<SCOREOFFSET);
	public static final int POINTSoff_NOREF_MODE_SUB=POINTSoff_NOREF|MODE_SUB;
	public static final int POINTSoff_MATCH=(POINTS_MATCH<<SCOREOFFSET);
	public static final int POINTSoff_MATCH2=(POINTS_MATCH2<<SCOREOFFSET);
	public static final int POINTSoff_SUB=(POINTS_SUB<<SCOREOFFSET);
	public static final int POINTSoff_SUB2=(POINTS_SUB2<<SCOREOFFSET);
	public static final int POINTSoff_INS=(POINTS_INS<<SCOREOFFSET);
	public static final int POINTSoff_INS2=(POINTS_INS2<<SCOREOFFSET);
	public static final int POINTSoff_DEL=(POINTS_DEL<<SCOREOFFSET);
	public static final int POINTSoff_DEL2=(POINTS_DEL2<<SCOREOFFSET);
	public static final int BADoff=(BAD<<SCOREOFFSET);
	public static final int MAXoff_SCORE=MAX_SCORE<<SCOREOFFSET;
	public static final int MINoff_SCORE=MIN_SCORE<<SCOREOFFSET;
	
	private int rows;
	private int columns;

	public long iterationsLimited=0;
	public long iterationsUnlimited=0;

	public boolean verbose=false;
	public boolean verbose2=false;
	
}
