package align2;
import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.SiteScore;

/**
 * Modification of MultiStateAligner9ts to replace fixed affine steps with an array */
public final class MultiStateAligner11tsJNI extends MSA{
	
	static {
		Shared.loadJNI();
	}

	private native void fillUnlimitedJNI(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int[] result, long[] iterationsUnlimited, int[] packed, int[] POINTSoff_SUB_ARRAY, int[] POINTSoff_INS_ARRAY, int maxRows, int maxColumns);

	private native void fillLimitedXJNI(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] result, long[] iterationsLimited, int[] packed, int[] POINTSoff_SUB_ARRAY, int[] POINTSoff_INS_ARRAY, int maxRows, int maxColumns, int bandwidth, float bandwidthRatio, int[] vertLimit, int[] horizLimit, byte[] baseToNumber, int[] POINTSoff_INS_ARRAY_C);

	public static void main(String[] args){
		byte[] read=args[0].getBytes();
		byte[] ref=args[1].getBytes();
		byte[] original=ref;
		
		MultiStateAligner11tsJNI msa=new MultiStateAligner11tsJNI(read.length, ref.length);
		System.out.println("Initial: ");
		//printMatrix(msa.packed, read.length, ref.length, TIMEMASK, SCOREOFFSET);
		
		int[] max=msa.fillLimited(read, ref, 0, ref.length-1, 0, null);
		
		System.out.println("Max: "+Arrays.toString(max));
		
		System.out.println("Final: ");
		//printMatrix(msa.packed, read.length, ref.length, TIMEMASK, SCOREOFFSET);
		
		byte[] out=msa.traceback(read, ref,  0, ref.length-1, max[0], max[1], max[2], false);
		
		int[] score=null;
		score=msa.score(read, ref,  0, ref.length-1, max[0], max[1], max[2], false);
		
		System.out.println(new String(ref));
		System.out.println(new String(read));
		System.out.println(new String(out));
		System.out.println("Score: "+Arrays.toString(score));
	}
	
	public MultiStateAligner11tsJNI(int maxRows_, int maxColumns_){
		super(maxRows_, maxColumns_);
		
		packed=KillSwitch.allocInt1D(3*(maxRows+1)*(maxColumns+1));
		grefbuffer=KillSwitch.allocByte1D(maxColumns+2);
		vertLimit=KillSwitch.allocInt1D(maxRows+1);
		horizLimit=KillSwitch.allocInt1D(maxColumns+1);
		
		Arrays.fill(vertLimit, BADoff);
		Arrays.fill(horizLimit, BADoff);
		
		for(int matrix=0; matrix<3; matrix++){
			for(int i=1; i<=maxRows; i++){
				for(int j=0; j<maxColumns+1; j++){
					packed[(matrix)*(maxRows+1)*(maxColumns+1)+(i)*(maxColumns+1)+(j)]|=BADoff;
				}
			}
			for(int i=0; i<=maxRows; i++){
				int prevScore=(i<2 ? 0 : packed[(matrix)*(maxRows+1)*(maxColumns+1)+(i-1)*(maxColumns+1)+(0)]);
				int score=prevScore+POINTSoff_INS_ARRAY[i];
				packed[(matrix)*(maxRows+1)*(maxColumns+1)+(i)*(maxColumns+1)+(0)]=score;
			}
		}
	}
	
	@Override
	public final int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] gaps){
		if(gaps==null){return fillLimitedX(read, ref, refStartLoc, refEndLoc, minScore);}
		else{
			byte[] gref=makeGref(ref, gaps, refStartLoc, refEndLoc);
			
			if(verbose && greflimit>0 && greflimit<500){
				System.err.println(new String(gref, 0, greflimit));
			}
			
			assert(gref!=null) : "Excessively long read:\n"+new String(read);
			return fillLimitedX(read, gref, 0, greflimit, minScore);
		}
	}
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Will not fill areas that cannot match minScore */
	private final int[] fillLimitedX(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		if(verbose){System.err.println("fillLimitedX");}
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;
		
		final int halfband=(bandwidth<1 && bandwidthRatio<=0) ? 0 :
			Tools.max(Tools.min(bandwidth<1 ? 9999999 : bandwidth, bandwidthRatio<=0 ? 9999999 : 8+(int)(rows*bandwidthRatio)), (columns-rows+8))/2;
		
		if(minScore<1 || (columns+rows<90) || ((halfband<1 || halfband*3>columns) && (columns>read.length+Tools.min(170, read.length+20)))){
			return fillUnlimited(read, ref, refStartLoc, refEndLoc);
		}

		minScore-=120; //Increases quality trivially

		//Create arrays for passing values to and from native library
		int[] result = KillSwitch.allocInt1D(5);
		long[] iterationsLimitedArray = new long[1];

		//Put values into array for passing to native library
		iterationsLimitedArray[0] = iterationsLimited;

		fillLimitedXJNI(read,ref,refStartLoc,refEndLoc,minScore,result,iterationsLimitedArray,packed,POINTSoff_SUB_ARRAY,POINTSoff_INS_ARRAY,maxRows,maxColumns,bandwidth,bandwidthRatio,vertLimit,horizLimit,AminoAcid.baseToNumber,POINTSoff_INS_ARRAY_C);

		//Retrieve variables from native library that were updated there
		iterationsLimited = iterationsLimitedArray[0];

		if(result[4]==1){
			return null;
		}

		//return new int[] {rows, maxCol, maxState, maxScore};
		return new int[] {result[0], result[1], result[2], result[3]};
	}
	
	@Override
	public final int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int[] gaps){
		if(gaps==null){return fillUnlimited(read, ref, refStartLoc, refEndLoc);}
		else{
			byte[] gref=makeGref(ref, gaps, refStartLoc, refEndLoc);
			assert(gref!=null) : "Excessively long read:\n"+new String(read);
			return fillUnlimited(read, gref, 0, greflimit);
		}
	}
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Does not require a min score (ie, same as old method) */
	private final int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc){

		//Create arrays for passing values to and from native library
		int[] result = KillSwitch.allocInt1D(4);
		long[] iterationsUnlimitedArray = new long[1];
		iterationsUnlimitedArray[0] = iterationsUnlimited;

		fillUnlimitedJNI(read,ref,refStartLoc,refEndLoc,result,iterationsUnlimitedArray,packed,POINTSoff_SUB_ARRAY,POINTSoff_INS_ARRAY,maxRows,maxColumns);

		//Retrieve variables from native library that were updated there
		long myiterationsUnlimited = iterationsUnlimitedArray[0];

		//return new int[] {rows, maxCol, maxState, maxScore};
		return new int[] {result[0], result[1], result[2], result[3]};
	}
	
	@Override
	@Deprecated
	/** return new int[] {rows, maxC, maxS, max}; */
	public final int[] fillQ(byte[] read, byte[] ref, byte[] baseScores, int refStartLoc, int refEndLoc){
		assert(false) : "Needs to be redone to work with score cutoffs.  Not difficult.";
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;

		assert(rows<=maxRows) : "Check that values are in-bounds before calling this function: "+rows+", "+maxRows;
		assert(columns<=maxColumns) : "Check that values are in-bounds before calling this function: "+columns+", "+maxColumns;
		
		assert(refStartLoc>=0) : "Check that values are in-bounds before calling this function: "+refStartLoc;
		assert(refEndLoc<ref.length) : "Check that values are in-bounds before calling this function: "+refEndLoc+", "+ref.length;
		
		for(int row=1; row<=rows; row++){
			for(int col=1; col<=columns; col++){
				final boolean match=(read[row-1]==ref[refStartLoc+col-1]);
				final boolean prevMatch=(row<2 || col<2 ? false : read[row-2]==ref[refStartLoc+col-2]);
				{//Calculate match and sub scores
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromDel=packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromIns=packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int streak=(packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&TIMEMASK);
					{//Calculate match/sub score
						if(match){
							int scoreMS=scoreFromDiag+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
							int scoreD=scoreFromDel+POINTSoff_MATCH;
							int scoreI=scoreFromIns+POINTSoff_MATCH;
							
							int score;
							int time;
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? streak+1 : 1);
							}else if(scoreD>=scoreI){
								score=scoreD;
								time=1;
							}else{
								score=scoreI;
								time=1;
							}
							score+=(((int)baseScores[row-1])<<SCOREOFFSET); //modifier
							
							if(time>MAX_TIME){time=MAX_TIME-MASK5;}
							assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
							assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
							packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col)]=(score|time);
							assert((score&SCOREMASK)==score);
							assert((time&TIMEMASK)==time);
							
						}else{
							int scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) :
								POINTSoff_SUB_ARRAY[streak+1]);
							int scoreD=scoreFromDel+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
							int scoreI=scoreFromIns+POINTSoff_SUB;
							
							int score;
							int time;
							byte prevState;
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? 1 : streak+1);
								prevState=MODE_MS;
							}else if(scoreD>=scoreI){
								score=scoreD;
								time=1;
								prevState=MODE_DEL;
							}else{
								score=scoreI;
								time=1;
								prevState=MODE_INS;
							}
							
							if(time>MAX_TIME){time=MAX_TIME-MASK5;}
							assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
							assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
							packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col)]=(score|time);
							assert((score&SCOREMASK)==score);
							assert((time&TIMEMASK)==time);
						}
					}
				}
				
				{//Calculate DEL score
					final int streak=packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col-1)]&TIMEMASK;
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromDel=packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col-1)]&SCOREMASK;
					
					int scoreMS=scoreFromDiag+POINTSoff_DEL;
					int scoreD=scoreFromDel+(streak==0 ? POINTSoff_DEL :
						streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
							streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 :
								streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 :
									((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));
					
					int score;
					int time;
					byte prevState;
					if(scoreMS>=scoreD){
						score=scoreMS;
						time=1;
						prevState=MODE_MS;
					}else{
						score=scoreD;
						time=streak+1;
						prevState=MODE_DEL;
					}
					
					if(time>MAX_TIME){time=MAX_TIME-MASK5;}
					assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
					assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
					packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col)]=(score|time);
					assert((score&SCOREMASK)==score);
					assert((time&TIMEMASK)==time);
				}
				
				{//Calculate INS score
					final int streak=packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col)]&TIMEMASK;
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col)]&SCOREMASK;
					final int scoreFromIns=packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col)]&SCOREMASK;
					
					int scoreMS=scoreFromDiag+POINTSoff_INS;
					int scoreI=scoreFromIns+POINTSoff_INS_ARRAY[streak+1];
					
					int score;
					int time;
					byte prevState;
					if(scoreMS>=scoreI){
						score=scoreMS;
						time=1;
						prevState=MODE_MS;
					}else{
						score=scoreI;
						time=streak+1;
						prevState=MODE_INS;
					}
					
					if(time>MAX_TIME){time=MAX_TIME-MASK5;}
					assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
					assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
					packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col)]=(score|time);
					assert((score&SCOREMASK)==score);
					assert((time&TIMEMASK)==time);
				}
			}
		}

		int maxCol=-1;
		int maxState=-1;
		int maxScore=Integer.MIN_VALUE;
		
		for(int state=0; state<3; state++){
			for(int col=1; col<=columns; col++){
				int x=packed[(state)*(maxRows+1)*(maxColumns+1)+(rows)*(maxColumns+1)+(col)]&SCOREMASK;
				if(x>maxScore){
					maxScore=x;
					maxCol=col;
					maxState=state;
				}
			}
		}
		maxScore>>=SCOREOFFSET;

		return new int[] {rows, maxCol, maxState, maxScore};
	}

	@Override
	/** @return {score, bestRefStart, bestRefStop} */
	/** Generates the match string */
	public final byte[] traceback(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state, boolean gapped){
		if(gapped){
			final byte[] gref=grefbuffer;
			int gstart=translateToGappedCoordinate(refStartLoc, gref);
			int gstop=translateToGappedCoordinate(refEndLoc, gref);
			byte[] out=traceback2(read, gref, gstart, gstop, row, col, state);
			return out;
		}else{
			return traceback2(read, ref, refStartLoc, refEndLoc, row, col, state);
		}
	}
	
	@Override
	/** Generates the match string */
	public final byte[] traceback2(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state){
		assert(refStartLoc<=refEndLoc) : refStartLoc+", "+refEndLoc;
		assert(row==rows);
		
		byte[] out=new byte[row+col-1]; //TODO if an out of bound crash occurs, try removing the "-1".
		int outPos=0;
		
		int gaps=0;
		
		if(state==MODE_INS){
		}
		
		while(row>0 && col>0){
			final int time=packed[(state)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col)]&TIMEMASK;
			final byte prev;
				
			if(state==MODE_MS){
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromDel=packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromIns=packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					if(scoreFromDiag>=scoreFromDel && scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
					else if(scoreFromDel>=scoreFromIns){prev=MODE_DEL;}
					else{prev=MODE_INS;}
				}
				
				byte c=read[row-1];
				byte r=ref[refStartLoc+col-1];
				if(c==r){
					out[outPos]='m';
				}else{
					if(!AminoAcid.isFullyDefined(c)){
						out[outPos]='N';
					}else if(!AminoAcid.isFullyDefined(r)){
						out[outPos]='N';
					}else{
						out[outPos]='S';
					}
				}
				
				row--;
				col--;
			}else if(state==MODE_DEL){
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromDel=packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col-1)]&SCOREMASK;
					if(scoreFromDiag>=scoreFromDel){prev=MODE_MS;}
					else{prev=MODE_DEL;}
				}
				
				byte r=ref[refStartLoc+col-1];
				if(r==GAPC){
					out[outPos]='-';
					gaps++;
				}else{
					out[outPos]='D';
				}
				col--;
			}else{
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col)]&SCOREMASK;
					final int scoreFromIns=packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col)]&SCOREMASK;
					if(scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
					else{prev=MODE_INS;}
				}
				
				assert(state==MODE_INS) : state;
				if(col==0){
					out[outPos]='X';
				}else if(col>=columns){
					out[outPos]='Y';
				}else{
					out[outPos]='I';
				}
				row--;
			}

			state=prev;
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
		
		byte[] out2=new byte[outPos];
		for(int i=0; i<outPos; i++){
			out2[i]=out[outPos-i-1];
		}
		out=null;
		
		if(gaps==0){return out2;}
		
		byte[] out3=new byte[out2.length+gaps*(GAPLEN-1)];
		for(int i=0, j=0; i<out2.length; i++){
			byte c=out2[i];
			if(c!=GAPC){
				out3[j]=c;
				j++;
			}else{
				int lim=j+GAPLEN;
				for(; j<lim; j++){
					out3[j]='D';
				}
			}
		}
		return out3;
	}
	
	@Override
	/** @return {score, bestRefStart, bestRefStop} */
	public final int[] score(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
			final int maxRow, final int maxCol, final int maxState, boolean gapped){
		if(gapped){
			if(verbose){
				System.err.println("score():");
				System.err.println("origin="+grefRefOrigin+", "+refStartLoc+", "+refEndLoc+", "+maxRow+", "+maxCol);
			}
			final byte[] gref=grefbuffer;
			int gstart=translateToGappedCoordinate(refStartLoc, gref);
			int gstop=translateToGappedCoordinate(refEndLoc, gref);
			
			assert(translateFromGappedCoordinate(gstart, gref)==refStartLoc); //TODO: Remove slow assertions
			assert(translateFromGappedCoordinate(gstop, gref)==refEndLoc);
			
			assert(gstart==0) : gstart; //TODO: skip translation if this is always zero
			
			if(verbose){System.err.println("gstart, gstop: "+gstart+", "+gstop);}
			int[] out=score2(read, gref, gstart, gstop, maxRow, maxCol, maxState);
			if(verbose){System.err.println("got score "+Arrays.toString(out));}
			
			assert(out[1]==translateToGappedCoordinate(translateFromGappedCoordinate(out[1], gref), gref)) :
				"Verifying: "+out[1]+" -> "+translateFromGappedCoordinate(out[1], gref)+" -> "+
				translateToGappedCoordinate(translateFromGappedCoordinate(out[1], gref), gref);
			assert(out[2]==translateToGappedCoordinate(translateFromGappedCoordinate(out[2], gref), gref));
			
			out[1]=translateFromGappedCoordinate(out[1], gref);
			out[2]=translateFromGappedCoordinate(out[2], gref);
			if(verbose){System.err.println("returning score "+Arrays.toString(out));}
			return out;
		}else{
			return score2(read, ref, refStartLoc, refEndLoc, maxRow, maxCol, maxState);
		}
	}
	
	@Override
	/** @return {score, bestRefStart, bestRefStop, maxRow, maxCol, maxState}, <br>
	 * or {score, bestRefStart, bestRefStop, maxRow, maxCol, maxState, padLeft, padRight} <br>
	 * if more padding is needed */
	public final int[] score2(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
			final int maxRow, final int maxCol, final int maxState){
		int row=maxRow;
		int col=maxCol;
		int state=maxState;

		assert(maxState>=0 && maxState<3) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		//assert(maxRow>=0 && maxRow<packed[0].length) :
		//	maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		//assert(maxCol>=0 && maxCol<packed[0][0].length) :
		//	maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		
		int score=packed[(maxState)*(maxRows+1)*(maxColumns+1)+(maxRow)*(maxColumns+1)+(maxCol)]&SCOREMASK; //Or zero, if it is to be recalculated
		
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
		
		if(verbose){System.err.println("Scoring.");}

		int stateTime=0;
		
		while(row>0 && col>0){
			if(verbose){System.err.println("state="+state+", row="+row+", col="+col);}

			final int time=packed[(state)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col)]&TIMEMASK;
			final byte prev;
			
			if(state==MODE_MS){
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromDel=packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromIns=packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col-1)]&SCOREMASK;
					if(scoreFromDiag>=scoreFromDel && scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
					else if(scoreFromDel>=scoreFromIns){prev=MODE_DEL;}
					else{prev=MODE_INS;}
				}
				row--;
				col--;
			}else if(state==MODE_DEL){
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col-1)]&SCOREMASK;
					final int scoreFromDel=packed[(MODE_DEL)*(maxRows+1)*(maxColumns+1)+(row)*(maxColumns+1)+(col-1)]&SCOREMASK;
					if(scoreFromDiag>=scoreFromDel){prev=MODE_MS;}
					else{prev=MODE_DEL;}
				}
				col--;
			}else{
				assert(state==MODE_INS);
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[(MODE_MS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col)]&SCOREMASK;
					final int scoreFromIns=packed[(MODE_INS)*(maxRows+1)*(maxColumns+1)+(row-1)*(maxColumns+1)+(col)]&SCOREMASK;
					if(scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
					else{prev=MODE_INS;}
				}
				row--;
			}
			
			if(col<0){
				if(verbose){
					System.err.println("Warning, column went below 0 at row="+row);
				}
				break; //prevents an out of bounds access
			}

			if(state==prev){stateTime++;}else{stateTime=0;}
			state=prev;
			
			if(verbose){System.err.println("state2="+state+", time="+time+", stateTime="+stateTime+", row2="+row+", col2="+col+"\n");}
		}
		if(row>col){
			col-=row;
		}
		
		final int bestRefStart=refStartLoc+col;
		
		score>>=SCOREOFFSET;
		
		if(verbose){
			System.err.println("bestRefStart="+bestRefStart+", refStartLoc="+refStartLoc);
			System.err.println("bestRefStop="+bestRefStop+", refEndLoc="+refEndLoc);
		}
		
		int padLeft=0;
		int padRight=0;
		if(bestRefStart<refStartLoc){
			padLeft=Tools.max(0, refStartLoc-bestRefStart);
		}else if(bestRefStart==refStartLoc && state==MODE_INS){
			padLeft=stateTime;
		}
		if(bestRefStop>refEndLoc){
			padRight=Tools.max(0, bestRefStop-refEndLoc);
		}else if(bestRefStop==refEndLoc && maxState==MODE_INS){
			padRight=packed[(maxState)*(maxRows+1)*(maxColumns+1)+(maxRow)*(maxColumns+1)+(maxCol)]&TIMEMASK;
		}
		
		int[] rvec;
		if(padLeft>0 || padRight>0){ //Suggest extra padding in cases of overflow
			rvec=new int[] {score, bestRefStart, bestRefStop, maxRow, maxCol, maxState, padLeft, padRight};
		}else{
			rvec=new int[] {score, bestRefStart, bestRefStop, maxRow, maxCol, maxState};
		}
		return rvec;
	}
	
	/**
	 * Fills grefbuffer
	 * @param ref
	 * @param a
	 * @param b
	 * @param gaps
	 * @return gref
	 */
	private final byte[] makeGref(byte[] ref, int[] gaps, int refStartLoc, int refEndLoc){
		assert(gaps!=null && gaps.length>0);
		
		assert(refStartLoc<=gaps[0]) : refStartLoc+", "+refEndLoc+", "+Arrays.toString(gaps);
		assert(refEndLoc>=gaps[gaps.length-1]);
		
		final int g0_old=gaps[0];
		final int gN_old=gaps[gaps.length-1];
		gaps[0]=Tools.min(gaps[0], refStartLoc);
		gaps[gaps.length-1]=Tools.max(gN_old, refEndLoc);
		grefRefOrigin=gaps[0];

		if(verbose){System.err.println("\ngaps2: "+Arrays.toString(gaps));}
		
		byte[] gref=grefbuffer;
		
		int gpos=0;
		for(int i=0; i<gaps.length; i+=2){
			int x=gaps[i];
			int y=gaps[i+1];
			
			for(int r=x; r<=y; r++, gpos++){
				//TODO: if out of bounds, use an 'N'
				assert(gpos<gref.length) :
					"\ngpos="+gpos+", gref.length="+gref.length+/*", read.length="+read.length+*/", gaps2="+Arrays.toString(gaps)+
					"\ni="+i+", r="+r+", x="+x+", y="+y+
					"\nGapTools.calcGrefLen("+gaps[0]+", "+gaps[gaps.length-1]+", gaps)="+GapTools.calcGrefLen(gaps[0], gaps[gaps.length-1], gaps)+
					"\nGapTools.calcGrefLen("+gaps[0]+", "+gaps[gaps.length-1]+", gaps)="+GapTools.calcGrefLen(gaps[0], gaps[gaps.length-1], gaps)+
					"\n"+refStartLoc+", "+refEndLoc+", "+greflimit+", "+GREFLIMIT2_CUSHION+"\n"+new String(gref)+"\n"/*+new String(read)+"\n"*/;
				gref[gpos]=ref[r];
			}
			
			if(i+2<gaps.length){
				int z=gaps[i+2];
				assert(z>y);
				int gap=z-y-1;
				assert(gap>=MINGAP) : gap+"\t"+MINGAP;
				if(gap<MINGAP){
					assert(false) : "TODO - just fill in normally";
				}else{
					int rem=gap%GAPLEN;
					int lim=y+GAPBUFFER+rem;
					
					int div=(gap-GAPBUFFER2)/GAPLEN;
					if(verbose){
						System.err.println("div = "+div);
					}
					assert(div>0);
					
					for(int r=y+1; r<=lim; r++, gpos++){
						gref[gpos]=ref[r];
					}
					for(int g=0; g<div; g++, gpos++){
						gref[gpos]=GAPC;
					}
					for(int r=z-GAPBUFFER; r<z; r++, gpos++){
						gref[gpos]=ref[r];
					}
				}
			}
		}
		
		greflimit=gpos;
		
		assert(gref[gpos-1]==ref[refEndLoc]);
		{
			final int lim=Tools.min(gref.length, greflimit+GREFLIMIT2_CUSHION);
			if(lim>gref.length){
				System.err.println("gref buffer overflow: "+lim+" > "+gref.length);
				return null;
			}
			for(int i=greflimit, r=refEndLoc+1; i<lim; i++, r++){
				gref[i]=(r<ref.length ? ref[r] : (byte)'N');
				greflimit2=i;
			}
		}
		
		if(verbose){
			System.err.println("gref:\n"+new String(gref));
		}
		
		gaps[0]=g0_old;
		gaps[gaps.length-1]=gN_old;

		if(verbose){
			System.err.println("\ngaps3: "+Arrays.toString(gaps));
		}
		
		return gref;
	}
	
	private final int translateFromGappedCoordinate(int point, byte[] gref){
		if(verbose){System.err.println("translateFromGappedCoordinate("+point+"), gro="+grefRefOrigin+", grl="+greflimit);}
		if(point<=0){return grefRefOrigin+point;}
		for(int i=0, j=grefRefOrigin; i<greflimit2; i++){
			byte c=gref[i];
			assert(point>=i) : "\n"+grefRefOrigin+"\n"+point+"\n"+new String(gref)+"\n";

			if(i==point){
				if(verbose){System.err.println(" -> "+j);}
				return j;
			}
			
			j+=(c==GAPC ? GAPLEN : 1);
		}

		System.err.println(grefRefOrigin);
		System.err.println(point);
		System.err.println(new String(gref));
		
		throw new RuntimeException("Out of bounds.");
	}
	
	private final int translateToGappedCoordinate(int point, byte[] gref){
		if(verbose){System.err.println("translateToGappedCoordinate("+point+"), gro="+grefRefOrigin+", grl="+greflimit);}
		if(point<=grefRefOrigin){return point-grefRefOrigin;}
		for(int i=0, j=grefRefOrigin; i<greflimit2; i++){
			assert(point>=j) : "\n"+grefRefOrigin+"\n"+point+"\n"+new String(gref)+"\n";
			byte c=gref[i];

			if(j==point){
				if(verbose){System.err.println(" -> "+i);}
				return i;
			}
			
			j+=(c==GAPC ? GAPLEN : 1);
		}

		System.err.println(grefRefOrigin);
		System.err.println(point);
		System.err.println(new String(gref));
		
		throw new RuntimeException("Out of bounds.");
	}

	/** Calculates score based on an array from Index */
	private final int calcAffineScore(int[] locArray){
		int score=0;
		int lastLoc=-2; //Last true location
		int lastValue=-1;
		int timeInMode=0;
		
		for(int i=0; i<locArray.length; i++){
			int loc=locArray[i];
			
			if(loc>0){//match
				if(loc==lastValue){//contiguous match
					score+=POINTS_MATCH2;
				}else if(loc==lastLoc || lastLoc<0){//match after a sub, or first match
					score+=POINTS_MATCH;
				}else if(loc<lastLoc){//deletion
					assert(lastLoc>=0);
					score+=POINTS_MATCH;
					score+=POINTS_DEL;
					int dif=lastLoc-loc+1;
					if(dif>MINGAP){
						int rem=dif%GAPLEN;
						int div=(dif-GAPBUFFER2)/GAPLEN;
						score+=(div*POINTS_GAP);
						assert(rem+GAPBUFFER2<dif);
						dif=rem+GAPBUFFER2;
						assert(dif>LIMIT_FOR_COST_4); //and probably LIMIT_FOR_COST_5
					}
					if(dif>LIMIT_FOR_COST_5){
						score+=((dif-LIMIT_FOR_COST_5+MASK5)/TIMESLIP)*POINTS_DEL5;
						dif=LIMIT_FOR_COST_5;
					}
					if(dif>LIMIT_FOR_COST_4){
						score+=(dif-LIMIT_FOR_COST_4)*POINTS_DEL4;
						dif=LIMIT_FOR_COST_4;
					}
					if(dif>LIMIT_FOR_COST_3){
						score+=(dif-LIMIT_FOR_COST_3)*POINTS_DEL3;
						dif=LIMIT_FOR_COST_3;
					}
					if(dif>1){
						score+=(dif-1)*POINTS_DEL2;
					}
					timeInMode=1;
				}else if(loc>lastLoc){//insertion
					assert(lastLoc>=0);
					score+=(POINTS_MATCH+POINTS_INS_ARRAY_C[lastLoc-loc]);
					timeInMode=1;
				}else{
					assert(false);
				}
				lastLoc=loc;
			}else{//substitution
				if(lastValue<0 && timeInMode>0){//contiguous
					timeInMode++;
					score+=(POINTS_SUB_ARRAY[timeInMode]);
				}else{
					score+=POINTS_SUB;
					timeInMode=1;
				}
			}
			lastValue=loc;
		}
		assert(score<=maxQuality(locArray.length));
		return score;
	}

	@Override
	public final int calcAffineScore(final int[] locArray, final byte[] baseScores, final byte bases[]){
		int score=0;
		int lastLoc=-3; //Last true location
		int lastValue=-1;
		int timeInMode=0;
		
		for(int i=0; i<locArray.length; i++){
			final int loc=locArray[i];
			
			if(loc>0){//match
				if(loc==lastValue){//contiguous match
					score+=(POINTS_MATCH2+baseScores[i]);
				}else if(loc==lastLoc || lastLoc<0){//match after a sub, or first match
					score+=(POINTS_MATCH+baseScores[i]);
				}else if(loc<lastLoc){//deletion
					assert(lastLoc>=0);
					score+=(POINTS_MATCH+baseScores[i]);
					score+=POINTS_DEL;
					int dif=lastLoc-loc+1;
					if(dif>MINGAP){
						int rem=dif%GAPLEN;
						int div=(dif-GAPBUFFER2)/GAPLEN;
						score+=(div*POINTS_GAP);
						assert(rem+GAPBUFFER2<dif);
						dif=rem+GAPBUFFER2;
						assert(dif>LIMIT_FOR_COST_4); //and probably LIMIT_FOR_COST_5
					}
					if(dif>LIMIT_FOR_COST_5){
						score+=((dif-LIMIT_FOR_COST_5+MASK5)/TIMESLIP)*POINTS_DEL5;
						dif=LIMIT_FOR_COST_5;
					}
					if(dif>LIMIT_FOR_COST_4){
						score+=(dif-LIMIT_FOR_COST_4)*POINTS_DEL4;
						dif=LIMIT_FOR_COST_4;
					}
					if(dif>LIMIT_FOR_COST_3){
						score+=(dif-LIMIT_FOR_COST_3)*POINTS_DEL3;
						dif=LIMIT_FOR_COST_3;
					}
					if(dif>1){
						score+=(dif-1)*POINTS_DEL2;
					}
					timeInMode=1;
				}else if(loc>lastLoc){//insertion
					assert(lastLoc>=0);
					score+=(POINTS_MATCH+baseScores[i]+POINTS_INS_ARRAY_C[Tools.min(loc-lastLoc, 5)]);
					timeInMode=1;
				}else{
					assert(false);
				}
				lastLoc=loc;
			}else if(loc==-1){//substitution
				if(lastValue<0 && timeInMode>0){//contiguous
					timeInMode++;
					score+=(POINTS_SUB_ARRAY[timeInMode]);
				}else{
					score+=POINTS_SUB;
					timeInMode=1;
				}
			}else{
				assert(loc==-2) : "\ni="+i+", loc="+loc+", score="+score+", lastLoc="+lastLoc+", lastValue="+lastValue
					+", time="+timeInMode+", length="+locArray.length+"\nbases=\n"+new String(bases)
					+"\nlocs[]=\n"+Arrays.toString(locArray)+"\n"+new String(bases)+"\n"
					+"If this happens please ensure that the reference has a startpad of Ns longer than readlength.";//N (no-call or no-ref)
				timeInMode=0;
				score+=POINTS_NOCALL;
			}
			lastValue=loc;
		}
		assert(score<=maxQuality(locArray.length));
		return score;
	}
	
	@Override
	public final int calcAffineScore(int[] locArray, byte[] baseScores, byte[] bases, int minContig){
		assert(minContig>1) : minContig;
		
		int contig=0;
		int maxContig=0;
		
		int score=0;
		int lastLoc=-3; //Last true location
		int lastValue=-1;
		int timeInMode=0;
		
		for(int i=0; i<locArray.length; i++){
			int loc=locArray[i];
			
			if(loc>0){//match
				if(loc==lastValue){//contiguous match
					contig++;
					score+=(POINTS_MATCH2+baseScores[i]);
				}else if(loc==lastLoc || lastLoc<0){//match after a sub, or first match
					maxContig=Tools.max(maxContig, contig);
					contig=1;
					score+=(POINTS_MATCH+baseScores[i]);
				}else if(loc<lastLoc){//deletion
					maxContig=Tools.max(maxContig, contig);
					contig=0;
					assert(lastLoc>=0);
					score+=(POINTS_MATCH+baseScores[i]);
					score+=POINTS_DEL;
					int dif=lastLoc-loc+1;
					if(dif>MINGAP){
						int rem=dif%GAPLEN;
						int div=(dif-GAPBUFFER2)/GAPLEN;
						score+=(div*POINTS_GAP);
						assert(rem+GAPBUFFER2<dif);
						dif=rem+GAPBUFFER2;
						assert(dif>LIMIT_FOR_COST_4); //and probably LIMIT_FOR_COST_5
					}
					if(dif>LIMIT_FOR_COST_5){
						score+=((dif-LIMIT_FOR_COST_5+MASK5)/TIMESLIP)*POINTS_DEL5;
						dif=LIMIT_FOR_COST_5;
					}
					if(dif>LIMIT_FOR_COST_4){
						score+=(dif-LIMIT_FOR_COST_4)*POINTS_DEL4;
						dif=LIMIT_FOR_COST_4;
					}
					if(dif>LIMIT_FOR_COST_3){
						score+=(dif-LIMIT_FOR_COST_3)*POINTS_DEL3;
						dif=LIMIT_FOR_COST_3;
					}
					if(dif>1){
						score+=(dif-1)*POINTS_DEL2;
					}
					timeInMode=1;
				}else if(loc>lastLoc){//insertion
					maxContig=Tools.max(maxContig, contig);
					contig=0;
					assert(lastLoc>=0);
					score+=(POINTS_MATCH+baseScores[i]+POINTS_INS_ARRAY_C[Tools.min(loc-lastLoc, 5)]);
					timeInMode=1;
				}else{
					assert(false);
				}
				lastLoc=loc;
			}else if(loc==-1){//substitution
				if(lastValue<0 && timeInMode>0){//contiguous
					timeInMode++;
					score+=(POINTS_SUB_ARRAY[timeInMode]);
				}else{
					score+=POINTS_SUB;
					timeInMode=1;
				}
			}else{
				assert(loc==-2) : loc+"\n"+Arrays.toString(locArray)+"\n"+Arrays.toString(baseScores)+"\n"+new String(bases)+"\n"
						+"If this happens please ensure that the reference has a startpad of Ns longer than readlength.";//N (no-call or no-ref)
				timeInMode=0;
				score+=POINTS_NOCALL;
			}
			lastValue=loc;
		}
		assert(score<=maxQuality(locArray.length));
		if(Tools.max(contig, maxContig)<minContig){score=Tools.min(score, -50*locArray.length);}
		return score;
	}
	
	@Override
	public final int scoreNoIndels(byte[] read, byte[] ref, final int refStart){
		return scoreNoIndels(read, ref, refStart, null);
	}
	@Override
	public final int scoreNoIndels(byte[] read, byte[] ref, final int refStart, final SiteScore ss){
		int score=0;
		int mode=-1;
		int timeInMode=0;
		
		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=read.length;
		final int refStop=refStart+read.length;
		boolean semiperfect=true;
		int norefs=0;
		
		if(refStart<0){
			readStart=0-refStart;
			score+=POINTS_NOREF*readStart;
			norefs+=readStart;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			readStop-=dif;
			score+=POINTS_NOREF*dif;
			norefs+=dif;
		}
		
		for(int i=readStart; i<readStop; i++){
			byte c=read[i];
			byte r=ref[refStart+i];
			
			if(c==r && c!='N'){
				if(mode==MODE_MS){
					timeInMode++;
					score+=POINTS_MATCH2;
				}else{
					timeInMode=0;
					score+=POINTS_MATCH;
				}
				mode=MODE_MS;
			}else if(c<0 || c=='N'){
				score+=POINTS_NOCALL;
				semiperfect=false;
			}else if(r<0 || r=='N'){
				score+=POINTS_NOREF;
				norefs++;
			}else{
				if(mode==MODE_SUB){timeInMode++;}
				else{timeInMode=0;}
				
				score+=(POINTS_SUB_ARRAY[timeInMode+1]);
				mode=MODE_SUB;
				semiperfect=false;
			}
		}
		
		return score;
	}
	
	@Override
	public final byte[] genMatchNoIndels(byte[] read, byte[] ref, final int refStart){
		if(read==null || ref==null){return null;}
		
		final byte[] match=new byte[read.length];
		
		for(int i=0, j=refStart; i<read.length; i++, j++){
			byte c=read[i];
			byte r=(j<0 || j>=ref.length) ? (byte)'N' : ref[j];

			if(c=='N' || r=='N'){match[i]='N';}
			else if(c==r){match[i]='m';}
			else{match[i]='S';}
			
		}
		
		return match;
	}
	
	@Override
	public final int scoreNoIndels(byte[] read, byte[] ref, byte[] baseScores, final int refStart){
		return scoreNoIndels(read, ref, baseScores, refStart, null);
	}
	@Override
	public final int scoreNoIndels(byte[] read, byte[] ref, byte[] baseScores, final int refStart, SiteScore ss){
		int score=0;
		int mode=-1;
		int timeInMode=0;
		int norefs=0;
		
		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=read.length;
		final int refStop=refStart+read.length;
		boolean semiperfect=true;
		
		if(refStart<0){
			readStart=0-refStart;
			score+=POINTS_NOREF*readStart;
			norefs+=readStart;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			readStop-=dif;
			score+=POINTS_NOREF*dif;
			norefs+=dif;
		}
		
		for(int i=readStart; i<readStop; i++){
			byte c=read[i];
			byte r=ref[refStart+i];
			
			if(c==r && c!='N'){
				if(mode==MODE_MS){
					timeInMode++;
					score+=POINTS_MATCH2;
				}else{
					timeInMode=0;
					score+=POINTS_MATCH;
				}
				score+=baseScores[i];
				mode=MODE_MS;
			}else if(c<0 || c=='N'){
				score+=POINTS_NOCALL;
				semiperfect=false;
			}else if(r<0 || r=='N'){
				score+=POINTS_NOREF;
				norefs++;
			}else{
				if(mode==MODE_SUB){timeInMode++;}
				else{timeInMode=0;}
				
				score+=(POINTS_SUB_ARRAY[timeInMode+1]);
				mode=MODE_SUB;
				semiperfect=false;
			}
		}
		
		return score;
	}
	
	@Override
	public final int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, byte[] baseScores, final int refStart, byte[][] matchReturn){
		int score=0;
		int mode=-1;
		int timeInMode=0;
		
		assert(refStart<=ref.length) : refStart+", "+ref.length;
		
		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=read.length;
		final int refStop=refStart+read.length;
		if(refStart<0 || refStop>ref.length){return -99999;}
		if(refStart<0){
			readStart=0-refStart;
			score+=POINTS_NOREF*readStart;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			System.err.println(new String(read)+"\ndif="+dif+", ref.length="+ref.length+", refStop="+refStop);
			readStop-=dif;
			score+=POINTS_NOREF*dif;
		}
		assert(refStart+readStop<=ref.length) : "readStart="+readStart+", readStop="+readStop+
		", refStart="+refStart+", refStop="+refStop+", ref.length="+ref.length+", read.length="+read.length;
		
		assert(matchReturn!=null);
		assert(matchReturn.length==1);
		if(matchReturn[0]==null || matchReturn[0].length!=read.length){
			assert(matchReturn[0]==null || matchReturn[0].length<read.length) : matchReturn[0].length+"!="+read.length;
			matchReturn[0]=new byte[read.length];
		}
		final byte[] match=matchReturn[0];
		
		for(int i=readStart; i<readStop; i++){
			byte c=read[i];
			byte r=ref[refStart+i];
			
			assert(r!='.' && c!='.');
			
			if(c==r && c!='N'){
				if(mode==MODE_MS){
					timeInMode++;
					score+=POINTS_MATCH2;
				}else{
					timeInMode=0;
					score+=POINTS_MATCH;
				}
				score+=baseScores[i];
				match[i]='m';
				mode=MODE_MS;
			}else if(c<0 || c=='N'){
				score+=POINTS_NOCALL;
				match[i]='N';
			}else if(r<0 || r=='N'){
				score+=POINTS_NOREF;
				match[i]='N';
			}else{
				match[i]='S';
				if(mode==MODE_SUB){timeInMode++;}
				else{timeInMode=0;}
				
				score+=(POINTS_SUB_ARRAY[timeInMode+1]);
				mode=MODE_SUB;
			}
		}
		
		return score;
	}
	
	@Override
	public final int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, final int refStart, byte[][] matchReturn){
		int score=0;
		int mode=-1;
		int timeInMode=0;
		
		assert(refStart<=ref.length) : refStart+", "+ref.length;
		
		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=read.length;
		final int refStop=refStart+read.length;
		if(refStart<0 || refStop>ref.length){return -99999;}
		if(refStart<0){
			readStart=0-refStart;
			score+=POINTS_NOREF*readStart;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			System.err.println(new String(read)+"\ndif="+dif+", ref.length="+ref.length+", refStop="+refStop);
			readStop-=dif;
			score+=POINTS_NOREF*dif;
		}
		assert(refStart+readStop<=ref.length) : "readStart="+readStart+", readStop="+readStop+
		", refStart="+refStart+", refStop="+refStop+", ref.length="+ref.length+", read.length="+read.length;
		
		assert(matchReturn!=null);
		assert(matchReturn.length==1);
		if(matchReturn[0]==null || matchReturn[0].length!=read.length){
			assert(matchReturn[0]==null || matchReturn[0].length<read.length) : matchReturn[0].length+"!="+read.length;
			matchReturn[0]=new byte[read.length];
		}
		final byte[] match=matchReturn[0];
		
		for(int i=readStart; i<readStop; i++){
			byte c=read[i];
			byte r=ref[refStart+i];
			
			assert(r!='.' && c!='.');
			
			if(c==r && c!='N'){
				if(mode==MODE_MS){
					timeInMode++;
					score+=POINTS_MATCH2;
				}else{
					timeInMode=0;
					score+=POINTS_MATCH;
				}
				match[i]='m';
				mode=MODE_MS;
			}else if(c<0 || c=='N'){
				score+=POINTS_NOCALL;
				match[i]='N';
			}else if(r<0 || r=='N'){
				score+=POINTS_NOREF;
				match[i]='N';
			}else{
				match[i]='S';
				if(mode==MODE_SUB){timeInMode++;}
				else{timeInMode=0;}
				
				if(AFFINE_ARRAYS){
					score+=(POINTS_SUB_ARRAY[timeInMode+1]);
				}else{
					if(timeInMode==0){score+=POINTS_SUB;}
					else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTS_SUB2;}
					else{score+=POINTS_SUB3;}
				}
				mode=MODE_SUB;
			}
		}
		
		return score;
	}
	
	@Override
	public final int maxQuality(int numBases){
		return POINTS_MATCH+(numBases-1)*(POINTS_MATCH2);
	}

	@Override
	public final int maxQuality(byte[] baseScores){
		return POINTS_MATCH+(baseScores.length-1)*(POINTS_MATCH2)+Tools.sumInt(baseScores);
	}

	@Override
	public final int maxImperfectScore(int numBases){
		int maxQ=maxQuality(numBases);
		int maxI=maxQ+Tools.min(POINTS_DEL, POINTS_INS-POINTS_MATCH2);
		assert(maxI<(maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB)));
		return maxI;
	}

	@Override
	public final int maxImperfectScore(byte[] baseScores){
		int maxQ=maxQuality(baseScores);
		int maxI=maxQ+Tools.min(POINTS_DEL, POINTS_INS-POINTS_MATCH2);
		assert(maxI<(maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB)));
		return maxI;
	}

	@Override
	public int calcDelScore(int len, boolean approximateGaps){
		if(len<=0){return 0;}
		int score=POINTS_DEL;
		
		if(approximateGaps && len>MINGAP){
			int rem=len%GAPLEN;
			int div=(len-GAPBUFFER2)/GAPLEN;
			score+=(div*POINTS_GAP);
			assert(rem+GAPBUFFER2<len);
			len=rem+GAPBUFFER2;
			assert(len>LIMIT_FOR_COST_4); //and probably LIMIT_FOR_COST_5
		}
		
		if(len>LIMIT_FOR_COST_5){
			score+=((len-LIMIT_FOR_COST_5+MASK5)/TIMESLIP)*POINTS_DEL5;
			len=LIMIT_FOR_COST_5;
		}
		if(len>LIMIT_FOR_COST_4){
			score+=(len-LIMIT_FOR_COST_4)*POINTS_DEL4;
			len=LIMIT_FOR_COST_4;
		}
		if(len>LIMIT_FOR_COST_3){
			score+=(len-LIMIT_FOR_COST_3)*POINTS_DEL3;
			len=LIMIT_FOR_COST_3;
		}
		if(len>1){
			score+=(len-1)*POINTS_DEL2;
		}
		return score;
	}
	
	private static int calcDelScoreOffset(int len){
		if(len<=0){return 0;}
		int score=POINTSoff_DEL;
		
		if(len>LIMIT_FOR_COST_5){
			score+=((len-LIMIT_FOR_COST_5+MASK5)/TIMESLIP)*POINTSoff_DEL5;
			len=LIMIT_FOR_COST_5;
		}
		if(len>LIMIT_FOR_COST_4){
			score+=(len-LIMIT_FOR_COST_4)*POINTSoff_DEL4;
			len=LIMIT_FOR_COST_4;
		}
		if(len>LIMIT_FOR_COST_3){
			score+=(len-LIMIT_FOR_COST_3)*POINTSoff_DEL3;
			len=LIMIT_FOR_COST_3;
		}
		if(len>1){
			score+=(len-1)*POINTSoff_DEL2;
		}
		return score;
	}

	@Override
	public int calcInsScore(int len){
		if(len<=0){return 0;}
		if(AFFINE_ARRAYS){
			return POINTS_INS_ARRAY_C[len];
		}else{
			int score=POINTS_INS;

			if(len>LIMIT_FOR_COST_4){
				score+=(len-LIMIT_FOR_COST_4)*POINTS_INS4;
				len=LIMIT_FOR_COST_4;
			}
			if(len>LIMIT_FOR_COST_3){
				score+=(len-LIMIT_FOR_COST_3)*POINTS_INS3;
				len=LIMIT_FOR_COST_3;
			}
			if(len>1){
				score+=(len-1)*POINTS_INS2;
			}
			return score;
		}
	}
	
	private static int calcInsScoreOffset(int len){
		if(len<=0){return 0;}
		if(AFFINE_ARRAYS){
			return POINTSoff_INS_ARRAY_C[len];
		}else{
			int score=POINTSoff_INS;
			if(len>LIMIT_FOR_COST_4){
				score+=(len-LIMIT_FOR_COST_4)*POINTSoff_INS4;
				len=LIMIT_FOR_COST_4;
			}
			if(len>LIMIT_FOR_COST_3){
				score+=(len-LIMIT_FOR_COST_3)*POINTSoff_INS3;
				len=LIMIT_FOR_COST_3;
			}
			if(len>1){
				score+=(len-1)*POINTSoff_INS2;
			}
			return score;
		}
	}
	
	private final int[]packed;
	private final byte[] grefbuffer;
	private int greflimit=-1;
	private int greflimit2=-1;
	private int grefRefOrigin=-1;
	
	@Override
	/**DO NOT MODIFY*/
	public final byte[] getGrefbuffer(){
		return grefbuffer;
	}

	public final int[] vertLimit;
	public final int[] horizLimit;

	@Override
	public CharSequence showVertLimit(){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<=rows; i++){sb.append(vertLimit[i]>>SCOREOFFSET).append(",");}
		return sb;
	}
	
	@Override
	public CharSequence showHorizLimit(){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<=columns; i++){sb.append(horizLimit[i]>>SCOREOFFSET).append(",");}
		return sb;
	}
	
	public static float minIdToMinRatio(double minid){
		if(minid>1){minid=minid/100;}
		assert(minid>0 && minid<=1) : "Min identity should be between 0 and 1.  Values above 1 will be assumed to be percent and divided by 100.";
		double matchdif=POINTS_MATCH-POINTS_MATCH2;
		double match=POINTS_MATCH2;
		double sub=-POINTS_MATCH2+0.5*(matchdif+POINTS_SUB)+0.5*POINTS_SUB2;
		double del=0.1*(matchdif+POINTS_DEL)+0.2*POINTS_DEL2+0.4*POINTS_DEL3+0.3*POINTS_DEL4;
		double ins=-POINTS_MATCH2+0.4*(matchdif+POINTS_INS)+0.3*(POINTS_INS2)+0.3*(POINTS_INS3);
		double badAvg=.7*sub+.2*del+.1*ins;
		double badFraction=1-minid;
		double minratio=(match+badFraction*badAvg)/match;
		assert(minratio<=1);
		minratio=Tools.max(0.1, minratio);
		return (float)minratio;
	}
	
	public static final int TIMEBITS=11;
	public static final int SCOREBITS=32-TIMEBITS;
	public static final int MAX_TIME=((1<<TIMEBITS)-1);
	public static final int MAX_SCORE=((1<<(SCOREBITS-1))-1)-2000;
	public static final int MIN_SCORE=0-MAX_SCORE; //Keeps it 1 point above "BAD".
	
	public static final int SCOREOFFSET=TIMEBITS;
	
	public static final int TIMEMASK=~((-1)<<TIMEBITS);
	public static final int SCOREMASK=(~((-1)<<SCOREBITS))<<SCOREOFFSET;
	
	public static final int POINTS_NOREF=0;
	public static final int POINTS_NOCALL=0;
	public static final int POINTS_MATCH=70;
	public static final int POINTS_MATCH2=100; //Note:  Changing to 90 substantially reduces false positives
	public static final int POINTS_COMPATIBLE=50;
	public static final int POINTS_SUB=-127;
	public static final int POINTS_SUBR=-147; //increased penalty if prior match streak was at most 1
	public static final int POINTS_SUB2=-51;
	public static final int POINTS_SUB3=-25;
	public static final int POINTS_MATCHSUB=-10;
	public static final int POINTS_INS=-395;
	public static final int POINTS_INS2=-39;
	public static final int POINTS_INS3=-23;
	public static final int POINTS_INS4=-8;
	public static final int POINTS_DEL=-472;
	public static final int POINTS_DEL2=-33;
	public static final int POINTS_DEL3=-9;
	public static final int POINTS_DEL4=-1;
	public static final int POINTS_DEL5=-1;
	public static final int POINTS_DEL_REF_N=-10;
	public static final int POINTS_GAP=0-GAPCOST;

	public static final int TIMESLIP=4;
	public static final int MASK5=TIMESLIP-1;
	static{assert(Integer.bitCount(TIMESLIP)==1);}
	
	private static final int BARRIER_I1=2;
	private static final int BARRIER_D1=3;

	public static final int LIMIT_FOR_COST_3=5;
	public static final int LIMIT_FOR_COST_4=20;
	public static final int LIMIT_FOR_COST_5=80;
	
	public static final int BAD=MIN_SCORE-1;
	
	public static final int POINTSoff_NOREF=(POINTS_NOREF<<SCOREOFFSET);
	public static final int POINTSoff_NOCALL=(POINTS_NOCALL<<SCOREOFFSET);
	public static final int POINTSoff_MATCH=(POINTS_MATCH<<SCOREOFFSET);
	public static final int POINTSoff_MATCH2=(POINTS_MATCH2<<SCOREOFFSET);
	public static final int POINTSoff_COMPATIBLE=(POINTS_COMPATIBLE<<SCOREOFFSET);
	public static final int POINTSoff_SUB=(POINTS_SUB<<SCOREOFFSET);
	public static final int POINTSoff_SUBR=(POINTS_SUBR<<SCOREOFFSET);
	public static final int POINTSoff_SUB2=(POINTS_SUB2<<SCOREOFFSET);
	public static final int POINTSoff_SUB3=(POINTS_SUB3<<SCOREOFFSET);
	public static final int POINTSoff_MATCHSUB=(POINTS_MATCHSUB<<SCOREOFFSET);
	public static final int POINTSoff_INS=(POINTS_INS<<SCOREOFFSET);
	public static final int POINTSoff_INS2=(POINTS_INS2<<SCOREOFFSET);
	public static final int POINTSoff_INS3=(POINTS_INS3<<SCOREOFFSET);
	public static final int POINTSoff_INS4=(POINTS_INS4<<SCOREOFFSET);
	public static final int POINTSoff_DEL=(POINTS_DEL<<SCOREOFFSET);
	public static final int POINTSoff_DEL2=(POINTS_DEL2<<SCOREOFFSET);
	public static final int POINTSoff_DEL3=(POINTS_DEL3<<SCOREOFFSET);
	public static final int POINTSoff_DEL4=(POINTS_DEL4<<SCOREOFFSET);
	public static final int POINTSoff_DEL5=(POINTS_DEL5<<SCOREOFFSET);
	public static final int POINTSoff_GAP=(POINTS_GAP<<SCOREOFFSET);
	public static final int POINTSoff_DEL_REF_N=(POINTS_DEL_REF_N<<SCOREOFFSET);
	public static final int BADoff=(BAD<<SCOREOFFSET);
	public static final int MAXoff_SCORE=MAX_SCORE<<SCOREOFFSET;
	public static final int MINoff_SCORE=MIN_SCORE<<SCOREOFFSET;
	
	public static final boolean AFFINE_ARRAYS=true;
	public static final int[] POINTS_INS_ARRAY;
	public static final int[] POINTSoff_INS_ARRAY;
	public static final int[] POINTS_INS_ARRAY_C;
	public static final int[] POINTSoff_INS_ARRAY_C;

	public static final int[] POINTS_SUB_ARRAY;
	public static final int[] POINTSoff_SUB_ARRAY;
	public static final int[] POINTS_SUB_ARRAY_C;
	public static final int[] POINTSoff_SUB_ARRAY_C;
	
	static{
		POINTS_INS_ARRAY=new int[604];
		POINTSoff_INS_ARRAY=new int[604];
		POINTS_INS_ARRAY_C=new int[604];
		POINTSoff_INS_ARRAY_C=new int[604];
		
		for(int i=1; i<POINTS_INS_ARRAY.length; i++){
			int pts, ptsoff;
			if(i>LIMIT_FOR_COST_4){
				pts=POINTS_INS4;
				ptsoff=POINTSoff_INS4;
			}else if(i>LIMIT_FOR_COST_3){
				pts=POINTS_INS3;
				ptsoff=POINTSoff_INS3;
			}else if(i>1){
				pts=POINTS_INS2;
				ptsoff=POINTSoff_INS2;
			}else{
				pts=POINTS_INS;
				ptsoff=POINTSoff_INS;
			}
			POINTS_INS_ARRAY[i]=pts;
			POINTSoff_INS_ARRAY[i]=ptsoff;
			POINTS_INS_ARRAY_C[i]=Tools.max(MIN_SCORE, pts+POINTS_INS_ARRAY_C[i-1]);
			POINTSoff_INS_ARRAY_C[i]=Tools.max(MINoff_SCORE, ptsoff+POINTSoff_INS_ARRAY_C[i-1]);
		}
		
		POINTS_SUB_ARRAY=new int[604];
		POINTSoff_SUB_ARRAY=new int[604];
		POINTS_SUB_ARRAY_C=new int[604];
		POINTSoff_SUB_ARRAY_C=new int[604];
		
		for(int i=1; i<POINTS_SUB_ARRAY.length; i++){
			int pts, ptsoff;
			if(i>LIMIT_FOR_COST_3){
				pts=POINTS_SUB3;
				ptsoff=POINTSoff_SUB3;
			}else if(i>1){
				pts=POINTS_SUB2;
				ptsoff=POINTSoff_SUB2;
			}else{
				pts=POINTS_SUB;
				ptsoff=POINTSoff_SUB;
			}
			POINTS_SUB_ARRAY[i]=pts;
			POINTSoff_SUB_ARRAY[i]=ptsoff;
			POINTS_SUB_ARRAY_C[i]=Tools.max(MIN_SCORE, pts+POINTS_SUB_ARRAY_C[i-1]);
			POINTSoff_SUB_ARRAY_C[i]=Tools.max(MINoff_SCORE, ptsoff+POINTSoff_SUB_ARRAY_C[i-1]);
		}
	}
	
	@Override
	public final int POINTS_NOREF(){return POINTS_NOREF;}
	@Override
	public final int POINTS_NOCALL(){return POINTS_NOCALL;}
	@Override
	public final int POINTS_MATCH(){return POINTS_MATCH;}
	@Override
	public final int POINTS_MATCH2(){return POINTS_MATCH2;}
	@Override
	public final int POINTS_COMPATIBLE(){return POINTS_COMPATIBLE;}
	@Override
	public final int POINTS_SUB(){return POINTS_SUB;}
	@Override
	public final int POINTS_SUBR(){return POINTS_SUBR;}
	@Override
	public final int POINTS_SUB2(){return POINTS_SUB2;}
	@Override
	public final int POINTS_SUB3(){return POINTS_SUB3;}
	@Override
	public final int POINTS_MATCHSUB(){return POINTS_MATCHSUB;}
	@Override
	public final int POINTS_INS(){return POINTS_INS;}
	@Override
	public final int POINTS_INS2(){return POINTS_INS2;}
	@Override
	public final int POINTS_INS3(){return POINTS_INS3;}
	@Override
	public final int POINTS_INS4(){return POINTS_INS4;}
	@Override
	public final int POINTS_DEL(){return POINTS_DEL;}
	@Override
	public final int POINTS_DEL2(){return POINTS_DEL2;}
	@Override
	public final int POINTS_DEL3(){return POINTS_DEL3;}
	@Override
	public final int POINTS_DEL4(){return POINTS_DEL4;}
	@Override
	public final int POINTS_DEL5(){return POINTS_DEL5;}
	@Override
	public final int POINTS_DEL_REF_N(){return POINTS_DEL_REF_N;}
	@Override
	public final int POINTS_GAP(){return POINTS_GAP;}

	@Override
	public final int TIMESLIP(){return TIMESLIP;}
	@Override
	public final int MASK5(){return MASK5;}
	@Override
	public final int SCOREOFFSET(){return SCOREOFFSET();}
	
	@Override
	final int BARRIER_I1(){return BARRIER_I1;}
	@Override
	final int BARRIER_D1(){return BARRIER_D1;}

	@Override
	public final int LIMIT_FOR_COST_3(){return LIMIT_FOR_COST_3;}
	@Override
	public final int LIMIT_FOR_COST_4(){return LIMIT_FOR_COST_4;}
	@Override
	public final int LIMIT_FOR_COST_5(){return LIMIT_FOR_COST_5;}
	
	@Override
	public final int BAD(){return BAD;}
	
	private int rows;
	private int columns;
	
}
