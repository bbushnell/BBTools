package aligner;

import java.util.Arrays;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SiteScore;

/**
 * Based on MSA9ts, with transform scores tweaked for PacBio. */
public final class MultiStateAligner9PacBioAdapter_WithBarriers {
	
	
	public static void main(String[] args){
		byte[] read=args[0].getBytes();
		byte[] ref=args[1].getBytes();
		
		byte[] original=ref;
		
		MultiStateAligner9PacBioAdapter_WithBarriers msa=new MultiStateAligner9PacBioAdapter_WithBarriers(read.length, ref.length);
		
		System.out.println("Initial: ");
		for(int mode=0; mode<msa.packed.length; mode++){
			for(int row=0; row<msa.packed[mode].length; row++){
				System.out.println(toScorePacked(msa.packed[mode][row]));
			}
			System.out.println();
			for(int row=0; row<msa.packed[mode].length; row++){
				System.out.println(toTimePacked(msa.packed[mode][row]));
			}
			System.out.println();
		}
		
		int[] max=msa.fillLimited(read, ref, 0, ref.length-1, 0, null);
		
		System.out.println("Max: "+Arrays.toString(max));
		
		System.out.println("Initial: ");
		for(int mode=0; mode<msa.packed.length; mode++){
			for(int row=0; row<msa.packed[mode].length; row++){
				System.out.println(toScorePacked(msa.packed[mode][row]));
			}
			System.out.println();
			for(int row=0; row<msa.packed[mode].length; row++){
				System.out.println(toTimePacked(msa.packed[mode][row]));
			}
			System.out.println();
		}
		
		byte[] out=msa.traceback(read, ref,  0, ref.length-1, max[0], max[1], max[2], false);
		
		int[] score=null;
		score=msa.score(read, ref,  0, ref.length-1, max[0], max[1], max[2], false);
		
		System.out.println(new String(ref));
		System.out.println(new String(read));
		System.out.println(new String(out));
		System.out.println("Score: "+Arrays.toString(score));
	}
	
	
	public MultiStateAligner9PacBioAdapter_WithBarriers(int maxRows_, int maxColumns_){
//		assert(maxColumns_>=200);
//		assert(maxRows_>=200);
		maxRows=maxRows_;
		maxColumns=maxColumns_;
		packed=new int[3][maxRows+1][maxColumns+1];
		grefbuffer=new byte[maxColumns+2];

		vertLimit=new int[maxRows+1];
		horizLimit=new int[maxColumns+1];
		Arrays.fill(vertLimit, BADoff);
		Arrays.fill(horizLimit, BADoff);
		
//		for(int i=0; i<maxColumns+1; i++){
//			scores[0][i]=0-i;
//		}
		
		for(int matrix=0; matrix<packed.length; matrix++){
			for(int i=1; i<=maxRows; i++){
				for(int j=0; j<packed[matrix][i].length; j++){
					packed[matrix][i][j]|=BADoff;
				}
//				packed[matrix][i][0]|=MODE_INS;
			}
//			for(int i=0; i<maxRows+1; i++){
//				scores[matrix][i][0]=(i*POINTSoff_NOREF);
//			}
			for(int i=0; i<=maxRows; i++){
				
				int prevScore=(i<2 ? 0 : packed[matrix][i-1][0]);
				int score=(i<2 ? (i*POINTSoff_INS) :
					(i<LIMIT_FOR_COST_3 ? prevScore+POINTSoff_INS2 :
						(i<LIMIT_FOR_COST_4 ? prevScore+POINTSoff_INS3 : prevScore+POINTSoff_INS4)));
				
				packed[matrix][i][0]=score;
			}
//			for(int i=1; i<maxColumns+1; i++){
//				prevState[matrix][0][i]=MODE_DEL;
//			}
//			for(int i=0; i<=maxColumns; i++){
//				packed[matrix][0][i]|=MODE_MS;
//			}
		}
	}
	
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Will not fill areas that cannot match minScore */
	public final int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] gaps){
		if(gaps==null){return fillLimitedX(read, ref, refStartLoc, refEndLoc, minScore);}
		else{
			byte[] gref=makeGref(ref, gaps, refStartLoc, refEndLoc);
			assert(gref!=null) : "Excessively long read:\n"+new String(read);
			return fillLimitedX(read, gref, 0, greflimit, minScore);
		}
	}
	
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Will not fill areas that cannot match minScore */
	private final int[] fillLimitedX(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
//		minScore=0;
//		assert(minScore>0);
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;
		
		if(/*read.length<40 || */false || minScore<=0 || columns>read.length+Tools.min(100, read.length)){
//			assert(false) : minScore;
			return fillUnlimited(read, ref, refStartLoc, refEndLoc);
		}

//		final int BARRIER_I2=columns-BARRIER_I1;
		final int BARRIER_I2=rows-BARRIER_I1, BARRIER_I2b=columns-1;
		final int BARRIER_D2=rows-BARRIER_D1;
		
		minScore-=100; //Increases quality trivially
		
		assert(rows<=maxRows) : "Check that values are in-bounds before calling this function: "+rows+", "+maxRows+"\n"+
			refStartLoc+", "+refEndLoc+", "+rows+", "+maxRows+", "+columns+", "+maxColumns+"\n"+new String(read)+"\n";
		assert(columns<=maxColumns) : "Check that values are in-bounds before calling this function: "+columns+", "+maxColumns+"\n"+
			refStartLoc+", "+refEndLoc+", "+rows+", "+maxRows+", "+columns+", "+maxColumns+"\n"+new String(read)+"\n";
		
		assert(refStartLoc>=0) : "Check that values are in-bounds before calling this function: "+refStartLoc+"\n"+
			refStartLoc+", "+refEndLoc+", "+rows+", "+maxRows+", "+columns+", "+maxColumns+"\n"+new String(read)+"\n";
		assert(refEndLoc<ref.length) : "Check that values are in-bounds before calling this function: "+refEndLoc+", "+ref.length+"\n"+
			refStartLoc+", "+refEndLoc+", "+rows+", "+maxRows+", "+columns+", "+maxColumns+"\n"+new String(read)+"\n";
		
//		for(int x=0; x<packed.length; x++){
//			for(int y=1; y<rows+1; y++){
//				Arrays.fill(packed[x][y], 1, columns+1, BADoff);
//			}
//		}
		for(int x=0; x<packed.length; x++){

//			Arrays.fill(packed[x][1], 1, columns+1, BADoff);
			Arrays.fill(packed[x][rows], 1, columns+1, BADoff);
//			for(int y=1; y<rows+1; y++){
//				Arrays.fill(packed[x][y], 1, columns+1, BADoff);
//			}
		}
		
		int minGoodCol=1;
		int maxGoodCol=columns;
		
		final int minScore_off=(minScore<<SCOREOFFSET);
		final int maxGain=(read.length-1)*POINTSoff_MATCH2+POINTSoff_MATCH;
		final int floor=minScore_off-maxGain;
//		final int subfloor=Tools.max(BADoff, floor-200*POINTSoff_MATCH2);
		final int subfloor=floor-5*POINTSoff_MATCH2;
		assert(subfloor>BADoff); //TODO: Actually, it needs to be substantially more.
		assert(subfloor<minScore_off) : minScore_off+", "+floor+", "+BADoff+", "+subfloor;
		
		if(verbose2){
			System.out.println();
			System.out.println("minScore="+minScore+"\t"+minScore_off);
			System.out.println("maxGain="+(maxGain>>SCOREOFFSET)+"\t"+(maxGain));
			System.out.println("floor="+(floor>>SCOREOFFSET)+"\t"+(floor));
			System.out.println("subfloor="+(subfloor>>SCOREOFFSET)+"\t"+(subfloor));
			System.out.println("BADoff="+(BADoff>>SCOREOFFSET)+"\t"+(BADoff));
			System.out.println("maxGain="+(maxGain>>SCOREOFFSET)+"\t"+(maxGain));
			System.out.println();
		}
		
		vertLimit[rows]=minScore_off;
		for(int i=rows-1; i>=0; i--){
			vertLimit[i]=Tools.max(vertLimit[i+1]-POINTSoff_MATCH2, floor);
		}
		
		horizLimit[columns]=minScore_off;
		for(int i=columns-1; i>=0; i--){
			horizLimit[i]=Tools.max(horizLimit[i+1]-POINTSoff_MATCH2, floor);
		}
		
		for(int row=1; row<=rows; row++){
			
			int colStart=minGoodCol;
			int colStop=maxGoodCol;
			
			minGoodCol=-1;
			maxGoodCol=-2;
			
			final int vlimit=vertLimit[row];
			
			if(verbose2){
				System.out.println();
				System.out.println("row="+row);
				System.out.println("colStart="+colStart);
				System.out.println("colStop="+colStop);
				System.out.println("vlimit="+(vlimit>>SCOREOFFSET)+"\t"+(vlimit));
			}
			
			if(colStart<0 || colStop<colStart){break;}
			
			
			if(colStart>1){
				assert(row>0);
				packed[MODE_MS][row][colStart-1]=subfloor;
				packed[MODE_INS][row][colStart-1]=subfloor;
				packed[MODE_DEL][row][colStart-1]=subfloor;
			}
			
			
			for(int col=colStart; col<=columns; col++){

				
				if(verbose2){
					System.out.println("\ncol "+col);
				}

				final byte call0=(row<2 ? (byte)'?' : read[row-2]);
				final byte call1=read[row-1];
				final byte ref0=(col<2 ? (byte)'!' : ref[refStartLoc+col-2]);
				final byte ref1=ref[refStartLoc+col-1];
				
				final boolean gap=(ref1==GAPC);
				assert(call1!=GAPC);

//				final boolean match=(read[row-1]==ref[refStartLoc+col-1]);
//				final boolean prevMatch=(row<2 || col<2 ? false : read[row-2]==ref[refStartLoc+col-2]);
				final boolean match=(call1==ref1 && ref1!='N');
				final boolean prevMatch=(call0==ref0 && ref0!='N');
				
//				System.err.println("")
				
				iterationsLimited++;
				final int limit=Tools.max(vlimit, horizLimit[col]);
				final int limit3=Tools.max(floor, (match ? limit-POINTSoff_MATCH2 : limit-POINTSoff_SUB3));

				final int delNeeded=Tools.max(0, row-col-1);
				final int insNeeded=Tools.max(0, (rows-row)-(columns-col)-1);

				final int delPenalty=calcDelScoreOffset(delNeeded);
				final int insPenalty=calcInsScoreOffset(insNeeded);
				
				
				final int scoreFromDiag_MS=packed[MODE_MS][row-1][col-1]&SCOREMASK;
				final int scoreFromDel_MS=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
				final int scoreFromIns_MS=packed[MODE_INS][row-1][col-1]&SCOREMASK;
				
				final int scoreFromDiag_DEL=packed[MODE_MS][row][col-1]&SCOREMASK;
				final int scoreFromDel_DEL=packed[MODE_DEL][row][col-1]&SCOREMASK;

				final int scoreFromDiag_INS=packed[MODE_MS][row-1][col]&SCOREMASK;
				final int scoreFromIns_INS=packed[MODE_INS][row-1][col]&SCOREMASK;
				
//				if(scoreFromDiag_MS<limit3 && scoreFromDel_MS<limit3 && scoreFromIns_MS<limit3
//						&& scoreFromDiag_DEL<limit && scoreFromDel_DEL<limit && scoreFromDiag_INS<limit && scoreFromIns_INS<limit){
//					iterationsLimited--; //A "fast" iteration
//				}
				
				if(gap || (scoreFromDiag_MS<=limit3 && scoreFromDel_MS<=limit3 && scoreFromIns_MS<=limit3)){
					packed[MODE_MS][row][col]=subfloor;
				}else{//Calculate match and sub scores
					final int streak=(packed[MODE_MS][row-1][col-1]&TIMEMASK);

					{//Calculate match/sub score
						
						int score;
						int time;
						byte prevState;
						
						if(match){

							int scoreMS=scoreFromDiag_MS+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
							int scoreD=scoreFromDel_MS+POINTSoff_MATCH;
							int scoreI=scoreFromIns_MS+POINTSoff_MATCH;
							
//							byte prevState;
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? streak+1 : 1);
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
							
//							if(time>MAX_TIME){time=MAX_TIME-MASK5;}
//							assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
//							assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
////							packed[MODE_MS][row][col]=(score|prevState|time);
//							packed[MODE_MS][row][col]=(score|time);
//							assert((score&SCOREMASK)==score);
////							assert((prevState&MODEMASK)==prevState);
//							assert((time&TIMEMASK)==time);
							
						}else{
							
							int scoreMS;
							if(ref1!='N' && call1!='N'){
								scoreMS=scoreFromDiag_MS+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) :
									(streak==0 ? POINTSoff_SUB : streak<LIMIT_FOR_COST_3 ? POINTSoff_SUB2 : POINTSoff_SUB3));
							}else{
								scoreMS=scoreFromDiag_MS+POINTSoff_NOCALL;
							}
							
							int scoreD=scoreFromDel_MS+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
							int scoreI=scoreFromIns_MS+POINTSoff_SUB;
							
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? 1 : streak+1);
//								time=(prevMatch ? (streak==1 ? 3 : 1) : streak+1);
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
							
//							if(time>MAX_TIME){time=MAX_TIME-MASK5;}
//							assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
//							assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
////							packed[MODE_MS][row][col]=(score|prevState|time);
//							packed[MODE_MS][row][col]=(score|time);
//							assert((score&SCOREMASK)==score);
////							assert((prevState&MODEMASK)==prevState);
//							assert((time&TIMEMASK)==time);
						}
						
						final int limit2;
						if(delNeeded>0){
							limit2=limit-delPenalty;
						}else if(insNeeded>0){
							limit2=limit-insPenalty;
						}else{
							limit2=limit;
						}
						assert(limit2>=limit);
						
						if(verbose2){System.err.println("MS: \tlimit2="+(limit2>>SCOREOFFSET)+"\t, score="+(score>>SCOREOFFSET));}
						
						if(score>=limit2){
							maxGoodCol=col;
							if(minGoodCol<0){minGoodCol=col;}
						}else{
							score=subfloor;
						}
						
						if(time>MAX_TIME){time=MAX_TIME-MASK5;}
						assert(score>=MINoff_SCORE || score==BADoff) : "Score overflow - use MSA2 instead";
						assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
//						packed[MODE_MS][row][col]=(score|prevState|time);
						packed[MODE_MS][row][col]=(score|time);
						assert((score&SCOREMASK)==score);
//						assert((prevState&MODEMASK)==prevState);
						assert((time&TIMEMASK)==time);
					}
				}
				
				if((scoreFromDiag_DEL<=limit && scoreFromDel_DEL<=limit) || row<BARRIER_D1 || row>BARRIER_D2){
//					assert((scoreFromDiag_DEL<=limit && scoreFromDel_DEL<=limit)) : scoreFromDiag_DEL+", "+row;
					packed[MODE_DEL][row][col]=subfloor;
				}else{//Calculate DEL score
							
					final int streak=packed[MODE_DEL][row][col-1]&TIMEMASK;
					
					int scoreMS=scoreFromDiag_DEL+POINTSoff_DEL;
					int scoreD=scoreFromDel_DEL+(streak==0 ? POINTSoff_DEL :
						streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
							streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 :
								streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 :
									((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));
//					int scoreI=scoreFromIns+POINTSoff_DEL;
					
					
					if(ref1=='N'){
						scoreMS+=POINTSoff_DEL_REF_N;
						scoreD+=POINTSoff_DEL_REF_N;
					}else if(gap){
						scoreMS+=POINTSoff_GAP;
						scoreD+=POINTSoff_GAP;
					}
					
					//if(match){scoreMS=subfloor;}
					
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
					
					final int limit2;
					if(insNeeded>0){
						limit2=limit-insPenalty;
					}else if(delNeeded>0){
						limit2=limit-calcDelScoreOffset(time+delNeeded)+calcDelScoreOffset(time);
					}else{
						limit2=limit;
					}
					assert(limit2>=limit);
					if(verbose2){System.err.println("DEL: \tlimit2="+(limit2>>SCOREOFFSET)+"\t, score="+(score>>SCOREOFFSET));}
					
					if(score>=limit2){
						maxGoodCol=col;
						if(minGoodCol<0){minGoodCol=col;}
					}else{
						score=subfloor;
					}
					
					if(time>MAX_TIME){time=MAX_TIME-MASK5;}
					assert(score>=MINoff_SCORE || score==BADoff) : "Score overflow - use MSA2 instead";
					assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
//					packed[MODE_DEL][row][col]=(score|prevState|time);
					packed[MODE_DEL][row][col]=(score|time);
					assert((score&SCOREMASK)==score);
//					assert((prevState&MODEMASK)==prevState);
					assert((time&TIMEMASK)==time);
				}

//				if(gap || (scoreFromDiag_INS<=limit && scoreFromIns_INS<=limit) || col<BARRIER_I1 || col>BARRIER_I2){
				if(gap || (scoreFromDiag_INS<=limit && scoreFromIns_INS<=limit) || (row<BARRIER_I1 && col>1) || (row>BARRIER_I2 && col<BARRIER_I2b)){
					packed[MODE_INS][row][col]=subfloor;
				}else{//Calculate INS score
					
					final int streak=packed[MODE_INS][row-1][col]&TIMEMASK;
					
					int scoreMS=scoreFromDiag_INS+POINTSoff_INS;
//					int scoreD=scoreFromDel+POINTSoff_INS;
					int scoreI=scoreFromIns_INS+(streak==0 ? POINTSoff_INS :
						streak<LIMIT_FOR_COST_3 ? POINTSoff_INS2 :
							streak<LIMIT_FOR_COST_4 ? POINTSoff_INS3 : POINTSoff_INS4);
					
					
//					System.err.println("("+row+","+col+")\t"+scoreFromDiag+"+"+POINTSoff_INS+"="+scoreM+", "+
//							scoreFromSub+"+"+POINTSoff_INS+"="+scoreS+", "
//							+scoreD+", "+scoreFromIns+"+"+
//							(streak==0 ? POINTSoff_INS : streak<LIMIT_FOR_COST_3 ? POINTSoff_INS2 : POINTSoff_INS3)+"="+scoreI);
					
					//if(match){scoreMS=subfloor;}
					
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
					
					final int limit2;
					if(delNeeded>0){
						limit2=limit-delPenalty;
					}else if(insNeeded>0){
						limit2=limit-calcInsScoreOffset(time+insNeeded)+calcInsScoreOffset(time);
					}else{
						limit2=limit;
					}
					assert(limit2>=limit);

					if(verbose2){System.err.println("INS: \tlimit2="+(limit2>>SCOREOFFSET)+"\t, score="+(score>>SCOREOFFSET));}
					if(score>=limit2){
						maxGoodCol=col;
						if(minGoodCol<0){minGoodCol=col;}
					}else{
						score=subfloor;
					}
					
					if(time>MAX_TIME){time=MAX_TIME-MASK5;}
					assert(score>=MINoff_SCORE || score==BADoff) : "Score overflow - use MSA2 instead";
					assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
//					packed[MODE_INS][row][col]=(score|prevState|time);
					packed[MODE_INS][row][col]=(score|time);
					assert((score&SCOREMASK)==score);
//					assert((prevState&MODEMASK)==prevState);
					assert((time&TIMEMASK)==time);
				}
				
				
				if(col>=colStop){
					if(col>colStop && maxGoodCol<col){break;}
					if(row>1){
						packed[MODE_MS][row-1][col+1]=subfloor;
						packed[MODE_INS][row-1][col+1]=subfloor;
						packed[MODE_DEL][row-1][col+1]=subfloor;
					}
				}
			}
		}
		

		int maxCol=-1;
		int maxState=-1;
		int maxScore=Integer.MIN_VALUE;
		
		for(int state=0; state<packed.length; state++){
			for(int col=1; col<=columns; col++){
				int x=packed[state][rows][col]&SCOREMASK;
				if(x>maxScore){
					maxScore=x;
					maxCol=col;
					maxState=state;
				}
			}
		}
		
		assert(maxScore>=BADoff);
//		if(maxScore==BADoff){
//			return null;
//		}
//		if(maxScore<floor){
//			return null;
//		}
		if(maxScore<minScore_off){
			return null;
		}
		
		maxScore>>=SCOREOFFSET;

//		System.err.println("Returning "+rows+", "+maxCol+", "+maxState+", "+maxScore);
		return new int[] {rows, maxCol, maxState, maxScore};
	}
	
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Will not fill areas that cannot match minScore */
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
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;
		
		final int maxGain=(read.length-1)*POINTSoff_MATCH2+POINTSoff_MATCH;
		final int subfloor=0-2*maxGain;
		assert(subfloor>BADoff && subfloor*2>BADoff); //TODO: Actually, it needs to be substantially more.
//		final int BARRIER_I2=columns-BARRIER_I1;
		final int BARRIER_I2=rows-BARRIER_I1, BARRIER_I2b=columns-1;
		final int BARRIER_D2=rows-BARRIER_D1;
		
		//temporary, for finding a bug
		if(rows>maxRows || columns>maxColumns){
			throw new RuntimeException("rows="+rows+", maxRows="+maxRows+", cols="+columns+", maxCols="+maxColumns+"\n"+new String(read)+"\n");
		}
		
		assert(rows<=maxRows) : "Check that values are in-bounds before calling this function: "+rows+", "+maxRows;
		assert(columns<=maxColumns) : "Check that values are in-bounds before calling this function: "+columns+", "+maxColumns;
		
		assert(refStartLoc>=0) : "Check that values are in-bounds before calling this function: "+refStartLoc;
		assert(refEndLoc<ref.length) : "Check that values are in-bounds before calling this function: "+refEndLoc+", "+ref.length;
		
		for(int row=1; row<=rows; row++){

//			int minc=max(1, row-20);
//			int maxc=min(columns, row+20);
			
			for(int col=1; col<=columns; col++){
				iterationsUnlimited++;
				
//				final boolean match=(read[row-1]==ref[refStartLoc+col-1]);
//				final boolean prevMatch=(row<2 || col<2 ? false : read[row-2]==ref[refStartLoc+col-2]);
				
				final byte call0=(row<2 ? (byte)'?' : read[row-2]);
				final byte call1=read[row-1];
				final byte ref0=(col<2 ? (byte)'!' : ref[refStartLoc+col-2]);
				final byte ref1=ref[refStartLoc+col-1];
				
				final boolean match=(call1==ref1 && ref1!='N');
				final boolean prevMatch=(call0==ref0 && ref0!='N');
				
				final boolean gap=(ref1==GAPC);
				assert(call1!=GAPC);

				if(gap){
					packed[MODE_MS][row][col]=subfloor;
				}else{//Calculate match and sub scores

					final int scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
					final int streak=(packed[MODE_MS][row-1][col-1]&TIMEMASK);

					{//Calculate match/sub score
						
						if(match){

							int scoreMS=scoreFromDiag+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
							int scoreD=scoreFromDel+POINTSoff_MATCH;
							int scoreI=scoreFromIns+POINTSoff_MATCH;
							
							int score;
							int time;
//							byte prevState;
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? streak+1 : 1);
//								prevState=MODE_MS;
							}else if(scoreD>=scoreI){
								score=scoreD;
								time=1;
//								prevState=MODE_DEL;
							}else{
								score=scoreI;
								time=1;
//								prevState=MODE_INS;
							}
							
							if(time>MAX_TIME){time=MAX_TIME-MASK5;}
							assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
							assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
//							packed[MODE_MS][row][col]=(score|prevState|time);
							packed[MODE_MS][row][col]=(score|time);
							assert((score&SCOREMASK)==score);
//							assert((prevState&MODEMASK)==prevState);
							assert((time&TIMEMASK)==time);
							
						}else{
							
							int scoreMS;
							if(ref1!='N' && call1!='N'){
								scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) :
									(streak==0 ? POINTSoff_SUB : streak<LIMIT_FOR_COST_3 ? POINTSoff_SUB2 : POINTSoff_SUB3));
							}else{
								scoreMS=scoreFromDiag+POINTSoff_NOCALL;
							}
							
							int scoreD=scoreFromDel+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
							int scoreI=scoreFromIns+POINTSoff_SUB;
							
							int score;
							int time;
							byte prevState;
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? 1 : streak+1);
//								time=(prevMatch ? (streak==1 ? 3 : 1) : streak+1);
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
//							packed[MODE_MS][row][col]=(score|prevState|time);
							packed[MODE_MS][row][col]=(score|time);
							assert((score&SCOREMASK)==score);
//							assert((prevState&MODEMASK)==prevState);
							assert((time&TIMEMASK)==time);
						}
					}
				}
				
				if(row<BARRIER_D1 || row>BARRIER_D2){
					packed[MODE_DEL][row][col]=subfloor;
				}else{//Calculate DEL score
							
					final int streak=packed[MODE_DEL][row][col-1]&TIMEMASK;
					
					final int scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;
					
					int scoreMS=scoreFromDiag+POINTSoff_DEL;
					int scoreD=scoreFromDel+(streak==0 ? POINTSoff_DEL :
						streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
							streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 :
								streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 :
									((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));
//					int scoreI=scoreFromIns+POINTSoff_DEL;
					
					if(ref1=='N'){
						scoreMS+=POINTSoff_DEL_REF_N;
						scoreD+=POINTSoff_DEL_REF_N;
					}else if(gap){
						scoreMS+=POINTSoff_GAP;
						scoreD+=POINTSoff_GAP;
					}
					
					//if(match){scoreMS=subfloor;}
					
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
//					packed[MODE_DEL][row][col]=(score|prevState|time);
					packed[MODE_DEL][row][col]=(score|time);
					assert((score&SCOREMASK)==score);
//					assert((prevState&MODEMASK)==prevState);
					assert((time&TIMEMASK)==time);
				}
				
				//Calculate INS score
//				if(gap || col<BARRIER_I1 || col>BARRIER_I2){
				if(gap || (row<BARRIER_I1 && col>1) || (row>BARRIER_I2 && col<BARRIER_I2b)){
					packed[MODE_INS][row][col]=subfloor;
				}else{//Calculate INS score
					
					final int streak=packed[MODE_INS][row-1][col]&TIMEMASK;

					final int scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;
					
					int scoreMS=scoreFromDiag+POINTSoff_INS;
//					int scoreD=scoreFromDel+POINTSoff_INS;
					int scoreI=scoreFromIns+(streak==0 ? POINTSoff_INS :
						streak<LIMIT_FOR_COST_3 ? POINTSoff_INS2 :
							streak<LIMIT_FOR_COST_4 ? POINTSoff_INS3 : POINTSoff_INS4);
					
//					System.err.println("("+row+","+col+")\t"+scoreFromDiag+"+"+POINTSoff_INS+"="+scoreM+", "+
//							scoreFromSub+"+"+POINTSoff_INS+"="+scoreS+", "
//							+scoreD+", "+scoreFromIns+"+"+
//							(streak==0 ? POINTSoff_INS : streak<LIMIT_FOR_COST_3 ? POINTSoff_INS2 : POINTSoff_INS3)+"="+scoreI);
					
					//if(match){scoreMS=subfloor;}
					
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
//					packed[MODE_INS][row][col]=(score|prevState|time);
					packed[MODE_INS][row][col]=(score|time);
					assert((score&SCOREMASK)==score);
//					assert((prevState&MODEMASK)==prevState);
					assert((time&TIMEMASK)==time);
				}
			}
		}
		

		int maxCol=-1;
		int maxState=-1;
		int maxScore=Integer.MIN_VALUE;
		
		for(int state=0; state<packed.length; state++){
			for(int col=1; col<=columns; col++){
				int x=packed[state][rows][col]&SCOREMASK;
				if(x>maxScore){
					maxScore=x;
					maxCol=col;
					maxState=state;
				}
			}
		}
		maxScore>>=SCOREOFFSET;

//		System.err.println("Returning "+rows+", "+maxCol+", "+maxState+", "+maxScore);
		return new int[] {rows, maxCol, maxState, maxScore};
	}
	
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

//			int minc=max(1, row-20);
//			int maxc=min(columns, row+20);
			
			for(int col=1; col<=columns; col++){
				
				final boolean match=(read[row-1]==ref[refStartLoc+col-1]);
				final boolean prevMatch=(row<2 || col<2 ? false : read[row-2]==ref[refStartLoc+col-2]);

				{//Calculate match and sub scores

					final int scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
					final int streak=(packed[MODE_MS][row-1][col-1]&TIMEMASK);

					{//Calculate match/sub score
						
						if(match){

							int scoreMS=scoreFromDiag+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
							int scoreD=scoreFromDel+POINTSoff_MATCH;
							int scoreI=scoreFromIns+POINTSoff_MATCH;
							
							int score;
							int time;
//							byte prevState;
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? streak+1 : 1);
//								prevState=MODE_MS;
							}else if(scoreD>=scoreI){
								score=scoreD;
								time=1;
//								prevState=MODE_DEL;
							}else{
								score=scoreI;
								time=1;
//								prevState=MODE_INS;
							}
							score+=(((int)baseScores[row-1])<<SCOREOFFSET); //modifier
							
							if(time>MAX_TIME){time=MAX_TIME-MASK5;}
							assert(score>=MINoff_SCORE) : "Score overflow - use MSA2 instead";
							assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
//							packed[MODE_MS][row][col]=(score|prevState|time);
							packed[MODE_MS][row][col]=(score|time);
							assert((score&SCOREMASK)==score);
//							assert((prevState&MODEMASK)==prevState);
							assert((time&TIMEMASK)==time);
							
						}else{

							int scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) :
								(streak==0 ? POINTSoff_SUB : streak<LIMIT_FOR_COST_3 ? POINTSoff_SUB2 : POINTSoff_SUB3));
							int scoreD=scoreFromDel+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
							int scoreI=scoreFromIns+POINTSoff_SUB;
							
							int score;
							int time;
							byte prevState;
							if(scoreMS>=scoreD && scoreMS>=scoreI){
								score=scoreMS;
								time=(prevMatch ? 1 : streak+1);
//								time=(prevMatch ? (streak==1 ? 3 : 1) : streak+1);
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
//							packed[MODE_MS][row][col]=(score|prevState|time);
							packed[MODE_MS][row][col]=(score|time);
							assert((score&SCOREMASK)==score);
//							assert((prevState&MODEMASK)==prevState);
							assert((time&TIMEMASK)==time);
						}
					}
				}
				
				{//Calculate DEL score
							
					final int streak=packed[MODE_DEL][row][col-1]&TIMEMASK;
					
					final int scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;
					
					int scoreMS=scoreFromDiag+POINTSoff_DEL;
					int scoreD=scoreFromDel+(streak==0 ? POINTSoff_DEL :
						streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
							streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 :
								streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 :
									((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));
//					int scoreI=scoreFromIns+POINTSoff_DEL;
					
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
//					packed[MODE_DEL][row][col]=(score|prevState|time);
					packed[MODE_DEL][row][col]=(score|time);
					assert((score&SCOREMASK)==score);
//					assert((prevState&MODEMASK)==prevState);
					assert((time&TIMEMASK)==time);
				}
				
				{//Calculate INS score
					
					final int streak=packed[MODE_INS][row-1][col]&TIMEMASK;

					final int scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;
					
					int scoreMS=scoreFromDiag+POINTSoff_INS;
//					int scoreD=scoreFromDel+POINTSoff_INS;
					int scoreI=scoreFromIns+(streak==0 ? POINTSoff_INS :
						streak<LIMIT_FOR_COST_3 ? POINTSoff_INS2 :
							streak<LIMIT_FOR_COST_4 ? POINTSoff_INS3 : POINTSoff_INS4);
					
//					System.err.println("("+row+","+col+")\t"+scoreFromDiag+"+"+POINTSoff_INS+"="+scoreM+", "+
//							scoreFromSub+"+"+POINTSoff_INS+"="+scoreS+", "
//							+scoreD+", "+scoreFromIns+"+"+
//							(streak==0 ? POINTSoff_INS : streak<LIMIT_FOR_COST_3 ? POINTSoff_INS2 : POINTSoff_INS3)+"="+scoreI);
					
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
//					packed[MODE_INS][row][col]=(score|prevState|time);
					packed[MODE_INS][row][col]=(score|time);
					assert((score&SCOREMASK)==score);
//					assert((prevState&MODEMASK)==prevState);
					assert((time&TIMEMASK)==time);
				}
			}
		}
		

		int maxCol=-1;
		int maxState=-1;
		int maxScore=Integer.MIN_VALUE;
		
		for(int state=0; state<packed.length; state++){
			for(int col=1; col<=columns; col++){
				int x=packed[state][rows][col]&SCOREMASK;
				if(x>maxScore){
					maxScore=x;
					maxCol=col;
					maxState=state;
				}
			}
		}
		maxScore>>=SCOREOFFSET;

//		System.err.println("Returning "+rows+", "+maxCol+", "+maxState+", "+maxScore);
		return new int[] {rows, maxCol, maxState, maxScore};
	}

	
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
	
	
	/** Generates the match string */
	public final byte[] traceback2(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state){
//		assert(false);
		assert(refStartLoc<=refEndLoc) : refStartLoc+", "+refEndLoc;
		assert(row==rows);
		
		byte[] out=new byte[row+col-1]; //TODO if an out of bound crash occurs, try removing the "-1".
		int outPos=0;
		
		int gaps=0;
		
		if(state==MODE_INS){
			//TODO ? Maybe not needed.
		}
		
		while(row>0 && col>0){
			
//			byte prev0=(byte)(packed[state][row][col]&MODEMASK);

			final int time=packed[state][row][col]&TIMEMASK;
			final byte prev;
				
//			System.err.println("state="+state+", prev="+prev+", row="+row+", col="+col+", score="+scores[state][row][col]);
			
			if(state==MODE_MS){
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
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
//						out[outPos]='X';
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
					final int scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;
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
					final int scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;
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

//			assert(prev==prev0);
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
		
		
		//Shrink and reverse the string
		byte[] out2=new byte[outPos];
		for(int i=0; i<outPos; i++){
			out2[i]=out[outPos-i-1];
		}
		out=null;
		
		if(gaps==0){return out2;}
		
		//TODO Consider outputting this compressed.
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
	
	/** @return {score, bestRefStart, bestRefStop} */
	public final int[] score2(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
			final int maxRow, final int maxCol, final int maxState){
		
		int row=maxRow;
		int col=maxCol;
		int state=maxState;

		assert(maxState>=0 && maxState<packed.length) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		assert(maxRow>=0 && maxRow<packed[0].length) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		assert(maxCol>=0 && maxCol<packed[0][0].length) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc);
		
		int score=packed[maxState][maxRow][maxCol]&SCOREMASK; //Or zero, if it is to be recalculated
		
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
//			System.err.println("state="+state+", row="+row+", col="+col);
			

			
//			byte prev0=(byte)(packed[state][row][col]&MODEMASK);

			final int time=packed[state][row][col]&TIMEMASK;
			final byte prev;
			
			if(state==MODE_MS){
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
					if(scoreFromDiag>=scoreFromDel && scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
					else if(scoreFromDel>=scoreFromIns){prev=MODE_DEL;}
					else{prev=MODE_INS;}
				}
				row--;
				col--;
			}else if(state==MODE_DEL){
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
					final int scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;
					if(scoreFromDiag>=scoreFromDel){prev=MODE_MS;}
					else{prev=MODE_DEL;}
				}
				col--;
			}else{
				assert(state==MODE_INS);
				if(time>1){prev=(byte)state;}
				else{
					final int scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
					final int scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;
					if(scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
					else{prev=MODE_INS;}
				}
				row--;
			}
			
			if(col<0){
				System.err.println(row);
				break; //prevents an out of bounds access
			
			}

//			assert(prev==prev0);
			state=prev;

//			System.err.println("state2="+state+", row2="+row+", col2="+col+"\n");
		}
//		assert(false) : row+", "+col;
		if(row>col){
			col-=row;
		}
		
		final int bestRefStart=refStartLoc+col;
		
		score>>=SCOREOFFSET;
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
	public final int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] gaps){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		
		int[] score;
		
		if(gaps==null){
			if(verbose){
				System.err.println("no gaps");
			}
			if(b-a>=maxColumns){
				System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
				assert(false) : refStartLoc+", "+refEndLoc;
				b=Tools.min(ref.length-1, a+maxColumns-1);
			}
			int[] max=fillLimited(read, ref, a, b, minScore, gaps);
			score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2], false));
		}else{
			if(verbose){System.err.println("\ngaps: "+Arrays.toString(gaps)+"\n"+new String(read)+"\ncoords: "+refStartLoc+", "+refEndLoc);}
			int[] max=fillLimited(read, ref, a, b, minScore, gaps);
			if(verbose){System.err.println("max: "+Arrays.toString(max));}
//			score=(max==null ? null : score(read, grefbuffer, 0, greflimit, max[0], max[1], max[2], true));
			score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2], true));
		}
		return score;
	}
	
	public final int[] fillAndScoreLimited(byte[] read, SiteScore ss, int thresh, int minScore){
		return fillAndScoreLimited(read, ss.chrom, ss.start, ss.stop, thresh, minScore, ss.gaps);
	}
	/*
	public final int[] fillAndScoreLimited_Gapped(byte[] read, SiteScore ss, int thresh, int minScore){
		if(ss.gaps==null){return fillAndScoreLimited(read, ss.chrom, ss.start, ss.stop, thresh, minScore);}
		int[] gaps=ss.gaps;
		final int bound1=gaps[0]=Tools.min(ss.start, gaps[0]);
		final int bound2=gaps[gaps.length-1]=Tools.max(ss.stop, gaps[gaps.length-1]);
		
		//This block is no longer needed since the array is preallocated.
		int len=0;
		final int gb2=GAPBUFFER*2;
		for(int i=0; i<gaps.length; i+=2){
			int x=gaps[i];
			int y=gaps[i+1];
			len+=(y-x+1);
			if(i+2<gaps.length){
				int z=gaps[i+2];
				assert(z>y);
				int gap=z-y-1;
				if(gap<MINGAP){
					len+=gap;
				}else{
					len+=gb2;
					gap-=gb2;
					int div=gap/GAPLEN;
					int rem=gap%GAPLEN;
					len+=(div+rem);
				}
			}
		}
		byte[] gref=grefbuffer;
		assert(gref.length>=len) : ss+"\t"+len+"\t"+gref.length;
		
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		
		for(int i=0, j=0; i<gaps.length; i+=2){
			int x=gaps[i];
			int y=gaps[i+1];
			
			for(int r=x; r<=y; r++, j++){
				gref[j]=cha.get(r);
			}
			
			if(i+2<gaps.length){
				int z=gaps[i+2];
				assert(z>y);
				int gap=z-y-1;
				assert(gap>=MINGAP);
				if(gap<MINGAP){
					assert(false) : "TODO - just fill in normally";
				}else{
					int div=gap/GAPLEN;
					assert(div>0);
					int rem=gap%GAPLEN;
					int lim=y+GAPBUFFER;
					
					for(int r=y+1; r<=lim; r++, j++){
						gref[j]=cha.get(r);
					}
					for(int g=0; g<div; g++, j++){
						gref[j]=GAPC;
					}
					for(int r=z-GAPBUFFER; r<z; r++, j++){
						gref[j]=cha.get(r);
					}
				}
			}
		}
//		fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore)
		int[] scoreArray=fillAndScoreLimited(read, ref, 0, ref.length-1, minScore);
		//Need to remap coordinates.
		
//		{score, bestRefStart, bestRefStop}
		if(scoreArray==null){return null;}
		
		int rstart=scoreArray[1];
		int rstop=scoreArray[2];
		
		int rstart2=-9999;
		int rstop2=-9999;
		
		for(int i=0, j=bound1; i<=len; i++){
			byte refc=ref[i];

			if(i==rstart){rstart2=j;}
			if(i==rstop){
				rstop2=j;
				assert(rstart2>-9999);
				break;
			}
			
			if(refc!=GAPC){
				j++;
			}else{
				j+=GAPLEN;
			}
		}
		assert(rstart2>-9999 && rstop2>-9999);
		scoreArray[1]=rstart2;
		scoreArray[2]=rstop2;
		
		return scoreArray;
	}*/
	
	/**
	 * Fills grefbuffer
	 * @param ref
	 * @param gaps
	 * @param refStartLoc
	 * @param refEndLoc
	 * @return gref
	 */
	public final byte[] makeGref(byte[] ref, int[] gaps, int refStartLoc, int refEndLoc){
		assert(gaps!=null && gaps.length>0);
		
		assert(refStartLoc<=gaps[0]) : refStartLoc+", "+refEndLoc+", "+Arrays.toString(gaps);
		assert(refEndLoc>=gaps[gaps.length-1]);
		
		final int g0_old=gaps[0];
		final int gN_old=gaps[gaps.length-1];
		gaps[0]=Tools.min(gaps[0], refStartLoc);
		gaps[gaps.length-1]=Tools.max(gN_old, refEndLoc);
		grefRefOrigin=gaps[0];

		if(verbose){System.err.println("\ngaps2: "+Arrays.toString(gaps));}
		
//		grefRefOrigin=Tools.min(gaps[0], refStartLoc);
		
//		//This block is no longer needed since the array is preallocated.
//		int len=0;
//		final int gb2=GAPBUFFER*2;
//		for(int i=0; i<gaps.length; i+=2){
//			int x=gaps[i];
//			int y=gaps[i+1];
//			len+=(y-x+1);
//			if(i+2<gaps.length){
//				int z=gaps[i+2];
//				assert(z>y);
//				int gap=z-y-1;
//				if(gap<MINGAP){
//					len+=gap;
//				}else{
//					len+=gb2;
//					gap-=gb2;
//					int div=gap/GAPLEN;
//					int rem=gap%GAPLEN;
//					len+=(div+rem);
//				}
//			}
//		}
		byte[] gref=grefbuffer;
		
		int gpos=0;
		for(int i=0; i<gaps.length; i+=2){
			int x=gaps[i];
			int y=gaps[i+1];
			
			for(int r=x; r<=y; r++, gpos++){
				//TODO: if out of bounds, use an 'N'
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
		
		//Add a cushion to the end to clear out the prior data (especially GAPC) that was there
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
	
	
//	public final int[] translateScoreFromGappedCoordinate(int[] score){
////		{score, bestRefStart, bestRefStop}
//		int a=score[1];
//		int b=score[2];
//		int a2=-9999;
//		int b2=-9999;
//		for(int i=0, j=grefRefOrigin; i<grefbuffer.length; i++){
//			byte c=grefbuffer[i];
//
//			if(i==a){a2=j;}
//			if(i==b){
//				b2=j;
//				assert(a2!=-9999);
//				score[1]=a2;
//				score[2]=b2;
//				return score;
//			}
//
//			j+=(c==GAPC ? GAPLEN : 1);
////			if(c!=GAPC){j++;}
////			else{j+=GAPLEN;}
//		}
//		throw new RuntimeException("Out of bounds.");
//	}
	
	private final int translateFromGappedCoordinate(int point, byte[] gref){
		if(verbose){System.err.println("translateFromGappedCoordinate("+point+"), gro="+grefRefOrigin+", grl="+greflimit);}
		if(point<=0){return grefRefOrigin+point;}
		for(int i=0, j=grefRefOrigin; i<greflimit2; i++){
//			if(verbose){System.err.println("i="+i+", j="+j+", sym="+(char)gref[i]);}
			byte c=gref[i];
			assert(point>=i) : "\n"+grefRefOrigin+"\n"+point+"\n"+new String(gref)+"\n";

			if(i==point){
				if(verbose){System.err.println(" -> "+j);}
				return j;
			}
			
			j+=(c==GAPC ? GAPLEN : 1);
//			if(c!=GAPC){j++;}
//			else{j+=GAPLEN;}
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
//			if(verbose){System.err.println("i="+i+", j="+j+", sym="+(char)gref[i]);}
			assert(point>=j) : "\n"+grefRefOrigin+"\n"+point+"\n"+new String(gref)+"\n";
			byte c=gref[i];

			if(j==point){
				if(verbose){System.err.println(" -> "+i);}
				return i;
			}
			
			j+=(c==GAPC ? GAPLEN : 1);
//			if(c!=GAPC){j++;}
//			else{j+=GAPLEN;}
		}

		System.err.println(grefRefOrigin);
		System.err.println(point);
		System.err.println(new String(gref));
		
		throw new RuntimeException("Out of bounds.");
	}
	
	public final int[] fillAndScoreLimited(byte[] read, int chrom, int start, int stop, int thresh, int minScore, int[] gaps){
		return fillAndScoreLimited(read, Data.getChromosome(chrom).array, start-thresh, stop+thresh, minScore, gaps);
	}
	
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, byte[] baseScores){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		if(b-a>=maxColumns){
			System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
			b=Tools.min(ref.length-1, a+maxColumns-1);
		}
		int[] max=fillQ(read, ref, baseScores, a, b);
//		int[] score=score(read, ref,  a, b, max[0], max[1], max[2]);
//		return score;
		return null;
	}
	
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, SiteScore ss, int thresh, byte[] baseScores){
		return fillAndScoreQ(read, ss.chrom, ss.start, ss.stop, thresh, baseScores);
	}
	
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, int chrom, int start, int stop, int thresh, byte[] baseScores){
		return fillAndScoreQ(read, Data.getChromosome(chrom).array, start-thresh, stop+thresh, baseScores);
	}
	
//	public final int scoreNoIndels(byte[] read, SiteScore ss){
//
//		ChromosomeArray cha=Data.getChromosome(ss.chrom);
//		final int refStart=ss.start;
//
//		int score=0;
//		int mode=MODE_START;
//		int timeInMode=0;
//		if(refStart<0 || refStart+read.length>cha.maxIndex+1){return -99999;} //TODO: Partial match
//
//		for(int i=0; i<read.length; i++){
//			byte c=read[i];
//			byte r=cha.get(refStart+i);
//
//			if(c==r){
//				if(mode==MODE_MS){
//					timeInMode++;
//					score+=POINTSoff_MATCH2;
//				}else{
//					timeInMode=0;
//					score+=POINTSoff_MATCH;
//				}
//				mode=MODE_MS;
//			}else if(c<0 || c=='N'){
//				score+=POINTSoff_NOCALL;
//			}else if(r<0 || r=='N'){
//				score+=POINTSoff_NOREF;
//			}else{
//				if(mode==MODE_SUB){timeInMode++;}
//				else{timeInMode=0;}
//
//				if(timeInMode==0){score+=POINTSoff_SUB;}
//				else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTSoff_SUB2;}
//				else{score+=POINTSoff_SUB3;}
//			}
//		}
//
//		return score;
//	}
	

	public final static int scoreNoIndels(byte[] read, SiteScore ss){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, ss.start, ss);
	}

	public final static int scoreNoIndels(byte[] read, final int chrom, final int refStart){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, refStart, null);
	}
	
	public final static int scoreNoIndels(byte[] read, SiteScore ss, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, baseScores, ss.start, ss);
	}

	public final static int scoreNoIndels(byte[] read, final int chrom, final int refStart, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, baseScores, refStart, null);
	}
	

//	public final int scoreNoIndels(byte[] read, final int chrom, final int refStart){
//
//		ChromosomeArray cha=Data.getChromosome(chrom);
//
//		int score=0;
//		int mode=MODE_START;
//		int timeInMode=0;
//
//		//This block handles cases where the read runs outside the reference
//		//Of course, padding the reference with 'N' would be better, but...
//		int readStart=0;
//		int readStop=read.length;
//		final int refStop=refStart+read.length;
//		if(refStart<0){
//			readStart=0-refStart;
//			score+=POINTSoff_NOREF*readStart;
//		}
//		if(refStop>cha.maxIndex+1){
//			int dif=(cha.maxIndex+1-refStop);
//			readStop-=dif;
//			score+=POINTSoff_NOREF*dif;
//		}
//
////		if(refStart<0 || refStart+read.length>cha.maxIndex+1){return -99999;} //No longer needed.
//
//		for(int i=readStart; i<readStop; i++){
//			byte c=read[i];
//			byte r=cha.get(refStart+i);
//
//			if(c==r){
//				if(mode==MODE_MS){
//					timeInMode++;
//					score+=POINTSoff_MATCH2;
//				}else{
//					timeInMode=0;
//					score+=POINTSoff_MATCH;
//				}
//				mode=MODE_MS;
//			}else if(c<0 || c=='N'){
//				score+=POINTSoff_NOCALL;
//			}else if(r<0 || r=='N'){
//				score+=POINTSoff_NOREF;
//			}else{
//				if(mode==MODE_SUB){timeInMode++;}
//				else{timeInMode=0;}
//
//				if(timeInMode==0){score+=POINTSoff_SUB;}
//				else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTSoff_SUB2;}
//				else{score+=POINTSoff_SUB3;}
//			}
//		}
//
//		return score;
//	}
	


	/** Calculates score based on an array from Index */
	public static final int calcAffineScore(int[] locArray){
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
//						assert(false) : div;
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
					score+=POINTS_MATCH;
					score+=POINTS_INS;
					int dif=Tools.min(loc-lastLoc+1, 5);
					assert(dif>0);
					if(dif>LIMIT_FOR_COST_4){
						score+=(dif-LIMIT_FOR_COST_4)*POINTS_INS4;
						dif=LIMIT_FOR_COST_4;
					}
					if(dif>LIMIT_FOR_COST_3){
						score+=(dif-LIMIT_FOR_COST_3)*POINTS_INS3;
						dif=LIMIT_FOR_COST_3;
					}
					if(dif>1){
						score+=(dif-1)*POINTS_INS2;
					}
					timeInMode=1;
				}else{
					assert(false);
				}
				lastLoc=loc;
			}else{//substitution
				if(lastValue<0 && timeInMode>0){//contiguous
					if(timeInMode<LIMIT_FOR_COST_3){score+=POINTS_SUB2;}
					else{score+=POINTS_SUB3;}
					timeInMode++;
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

	/** Calculates score based on an array from Index */
	public static final int calcAffineScore(int[] locArray, byte[] baseScores){
		int score=0;
		int lastLoc=-2; //Last true location
		int lastValue=-1;
		int timeInMode=0;
		
		for(int i=0; i<locArray.length; i++){
			int loc=locArray[i];
			
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
//						assert(false) : div;
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
					score+=(POINTS_MATCH+baseScores[i]);
					score+=POINTS_INS;
					int dif=Tools.min(loc-lastLoc+1, 5);
					assert(dif>0);
					if(dif>LIMIT_FOR_COST_4){
						score+=(dif-LIMIT_FOR_COST_4)*POINTS_INS4;
						dif=LIMIT_FOR_COST_4;
					}
					if(dif>LIMIT_FOR_COST_3){
						score+=(dif-LIMIT_FOR_COST_3)*POINTS_INS3;
						dif=LIMIT_FOR_COST_3;
					}
					if(dif>1){
						score+=(dif-1)*POINTS_INS2;
					}
					timeInMode=1;
				}else{
					assert(false);
				}
				lastLoc=loc;
			}else{//substitution
				if(lastValue<0 && timeInMode>0){//contiguous
					if(timeInMode<LIMIT_FOR_COST_3){score+=POINTS_SUB2;}
					else{score+=POINTS_SUB3;}
					timeInMode++;
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
	

	public final static int scoreNoIndels(byte[] read, byte[] ref, final int refStart, final SiteScore ss){
		
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
		
//		if(refStart<0 || refStart+read.length>ref.length){return -99999;} //No longer needed.
		
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
				
				if(timeInMode==0){score+=POINTS_SUB;}
				else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTS_SUB2;}
				else{score+=POINTS_SUB3;}
				mode=MODE_SUB;
				semiperfect=false;
			}
		}
		
		if(semiperfect && ss!=null){ss.semiperfect=((ss.stop==ss.start+read.length-1) && (norefs<=read.length/2));}
		
		return score;
	}
	

	public final static int scoreNoIndels(byte[] read, byte[] ref, byte[] baseScores, final int refStart, SiteScore ss){
		
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
		
//		if(refStart<0 || refStart+read.length>ref.length){return -99999;} //No longer needed.
		
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
				
				if(timeInMode==0){score+=POINTS_SUB;}
				else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTS_SUB2;}
				else{score+=POINTS_SUB3;}
				mode=MODE_SUB;
				semiperfect=false;
			}
		}
		
		if(semiperfect && ss!=null){ss.semiperfect=((ss.stop==ss.start+read.length-1) && (norefs<=read.length/2));}
		assert(Read.CHECKSITE(ss, read, -1));
		
		return score;
	}
	
	
	public final static int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, byte[] baseScores, final int refStart, byte[][] matchReturn){
		int score=0;
		int mode=-1;
		int timeInMode=0;
		
		assert(refStart<=ref.length) : refStart+", "+ref.length;
		
		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=read.length;
		final int refStop=refStart+read.length;
		if(refStart<0){
			readStart=0-refStart;
			score+=POINTS_NOREF*readStart;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			System.err.println("dif="+dif+", ref.length="+ref.length+", refStop="+refStop);
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
		
//		if(refStart<0 || refStart+read.length>ref.length){return -99999;} //No longer needed.
		
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
//				match[i]='m';
				match[i]='N';
			}else{
				match[i]='S';
				if(mode==MODE_SUB){timeInMode++;}
				else{timeInMode=0;}
				
				if(timeInMode==0){score+=POINTS_SUB;}
				else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTS_SUB2;}
				else{score+=POINTS_SUB3;}
				mode=MODE_SUB;
			}
		}
		
		return score;
	}
	
	
	public final static int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, final int refStart, byte[][] matchReturn){
		int score=0;
		int mode=-1;
		int timeInMode=0;
		
		assert(refStart<=ref.length) : refStart+", "+ref.length;
		
		//This block handles cases where the read runs outside the reference
		//Of course, padding the reference with 'N' would be better, but...
		int readStart=0;
		int readStop=read.length;
		final int refStop=refStart+read.length;
		if(refStart<0){
			readStart=0-refStart;
			score+=POINTS_NOREF*readStart;
		}
		if(refStop>ref.length){
			int dif=(refStop-ref.length);
			System.err.println("dif="+dif+", ref.length="+ref.length+", refStop="+refStop);
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
		
//		if(refStart<0 || refStart+read.length>ref.length){return -99999;} //No longer needed.
		
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
//				match[i]='m';
				match[i]='N';
			}else{
				match[i]='S';
				if(mode==MODE_SUB){timeInMode++;}
				else{timeInMode=0;}
				
				if(timeInMode==0){score+=POINTS_SUB;}
				else if(timeInMode<LIMIT_FOR_COST_3){score+=POINTS_SUB2;}
				else{score+=POINTS_SUB3;}
				mode=MODE_SUB;
			}
		}
		
		return score;
	}
	
	public static final int maxQuality(int numBases){
		return POINTS_MATCH+(numBases-1)*(POINTS_MATCH2);
	}
	
	public static final int maxQuality(byte[] baseScores){
		return POINTS_MATCH+(baseScores.length-1)*(POINTS_MATCH2)+Tools.sumInt(baseScores);
	}
	
	public static final int maxImperfectScore(int numBases){
//		int maxQ=maxQuality(numBases);
////		maxImperfectSwScore=maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB);
//		int maxI=maxQ+POINTS_DEL;
//		maxI=Tools.max(maxI, maxQ+POINTS_INS-POINTS_MATCH2);
//		maxI=Tools.min(maxI, maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB));
		
		int maxQ=maxQuality(numBases);
		int maxI=maxQ+Tools.min(POINTS_DEL, POINTS_INS-POINTS_MATCH2);
		assert(maxI<(maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB)));
		return maxI;
	}
	
	public static final int maxImperfectScore(byte[] baseScores){
//		int maxQ=maxQuality(numBases);
////		maxImperfectSwScore=maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB);
//		int maxI=maxQ+POINTS_DEL;
//		maxI=Tools.max(maxI, maxQ+POINTS_INS-POINTS_MATCH2);
//		maxI=Tools.min(maxI, maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB));
		
		int maxQ=maxQuality(baseScores);
		int maxI=maxQ+Tools.min(POINTS_DEL, POINTS_INS-POINTS_MATCH2);
		assert(maxI<(maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB)));
		return maxI;
	}
	
	public static final String toString(int[] a){
		
		int width=7;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num : a){
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public static final String toTimePacked(int[] a){
		int width=7;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num_ : a){
			int num=num_&TIMEMASK;
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public static final String toScorePacked(int[] a){
		int width=7;

		String minString=" -";
		String maxString="  ";
		while(minString.length()<width){minString+='9';}
		while(maxString.length()<width){maxString+='9';}
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num_ : a){
			int num=num_>>SCOREOFFSET;
			String s=" "+num;
			if(s.length()>width){s=num>0 ? maxString : minString;}
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public static final String toString(byte[] a){
		
		int width=6;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num : a){
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0);
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public static final String toString(byte[] ref, int startLoc, int stopLoc){
		StringBuilder sb=new StringBuilder(stopLoc-startLoc+1);
		for(int i=startLoc; i<=stopLoc; i++){sb.append((char)ref[i]);}
		return sb.toString();
	}
	
	public static int calcDelScore(int len){
		if(len<=0){return 0;}
		int score=POINTS_DEL;
		
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
	
	public static int calcInsScore(int len){
		if(len<=0){return 0;}
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
	
	private static int calcInsScoreOffset(int len){
		if(len<=0){return 0;}
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
	
	
	public final int maxRows;
	public final int maxColumns;

	private final int[][][] packed;
	private final byte[] grefbuffer;
	private int greflimit=-1;
	private int greflimit2=-1;
	private int grefRefOrigin=-1;

	public static final int GAPBUFFER=Shared.GAPBUFFER;
	public static final int GAPBUFFER2=Shared.GAPBUFFER2;
	public static final int GAPLEN=Shared.GAPLEN;
	public static final int MINGAP=Shared.MINGAP;
	public static final int GAPCOST=Shared.GAPCOST*2;
	public static final byte GAPC=Shared.GAPC;
	
	private static final int GREFLIMIT2_CUSHION=128; //Tools.max(GAPBUFFER2, GAPLEN);
	
	
	/**DO NOT MODIFY*/
	public final byte[] getGrefbuffer(){
		return grefbuffer;
	}

	public final int[] vertLimit;
	public final int[] horizLimit;

	CharSequence showVertLimit(){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<=rows; i++){sb.append(vertLimit[i]>>SCOREOFFSET).append(",");}
		return sb;
	}
	CharSequence showHorizLimit(){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<=columns; i++){sb.append(horizLimit[i]>>SCOREOFFSET).append(",");}
		return sb;
	}

//	public static final int MODEBITS=2;
	public static final int TIMEBITS=12;
	public static final int SCOREBITS=32-TIMEBITS;
	public static final int MAX_TIME=((1<<TIMEBITS)-1);
	public static final int MAX_SCORE=((1<<(SCOREBITS-1))-1)-2000;
	public static final int MIN_SCORE=0-MAX_SCORE; //Keeps it 1 point above "BAD".

//	public static final int MODEOFFSET=0; //Always zero.
//	public static final int TIMEOFFSET=0;
	public static final int SCOREOFFSET=TIMEBITS;

//	public static final int MODEMASK=~((-1)<<MODEBITS);
//	public static final int TIMEMASK=(~((-1)<<TIMEBITS))<<TIMEOFFSET;
	public static final int TIMEMASK=~((-1)<<TIMEBITS);
	public static final int SCOREMASK=(~((-1)<<SCOREBITS))<<SCOREOFFSET;
	
	private static final byte MODE_MS=0;
	private static final byte MODE_DEL=1;
	private static final byte MODE_INS=2;
	private static final byte MODE_SUB=3;
	
//	public static final int POINTS_NOREF=-1;
//	public static final int POINTS_MATCH=10;
//	public static final int POINTS_SUB=-15;
//	public static final int POINTS_SUB2=-10;
//	public static final int POINTS_SUB3=-2;
//	public static final int POINTS_INS=-30;
//	public static final int POINTS_INS2=-3;
//	public static final int POINTS_INS3=-1;
//	public static final int POINTS_DEL=-30;
//	public static final int POINTS_DEL2=-2;
//	public static final int POINTS_DEL3=-1;
	
//	public static final int POINTS_NOREF=-5;
//	public static final int POINTS_MATCH=10;
//	public static final int POINTS_SUB=-13;
//	public static final int POINTS_SUB2=-7;
//	public static final int POINTS_SUB3=-3;
//	public static final int POINTS_INS=-21;
//	public static final int POINTS_INS2=-2;
//	public static final int POINTS_INS3=-1;
//	public static final int POINTS_DEL=-20;
//	public static final int POINTS_DEL2=-2;
//	public static final int POINTS_DEL3=-1;
	
	public static final int POINTS_NOREF=-8;
	public static final int POINTS_NOCALL=-8;
	public static final int POINTS_MATCH=90;
	public static final int POINTS_MATCH2=100; //Note:  Changing to 90 substantially reduces false positives
	public static final int POINTS_COMPATIBLE=50;
	public static final int POINTS_SUB=-141;
	public static final int POINTS_SUBR=-159; //increased penalty if prior match streak was at most 1
	public static final int POINTS_SUB2=-49;
	public static final int POINTS_SUB3=-27;
	public static final int POINTS_MATCHSUB=-10;
	public static final int POINTS_INS=-204;
	public static final int POINTS_INS2=-42;
	public static final int POINTS_INS3=-25;
	public static final int POINTS_INS4=-8;
	public static final int POINTS_DEL=-287;
	public static final int POINTS_DEL2=-39;
	public static final int POINTS_DEL3=-21;
	public static final int POINTS_DEL4=-12;
	public static final int POINTS_DEL5=-8;
	public static final int POINTS_DEL_REF_N=-10;
	public static final int POINTS_GAP=0-GAPCOST;

	public static final int TIMESLIP=4;
	public static final int MASK5=TIMESLIP-1;
	static{assert(Integer.bitCount(TIMESLIP)==1);}
	
	//TODO:  Consider removing these barriers entirely for PacBio reads.  Would make code faster, too.
	private static final int BARRIER_I1=0;
	private static final int BARRIER_D1=0;
	

	public static final int LIMIT_FOR_COST_3=5;
	public static final int LIMIT_FOR_COST_4=25;
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
	
	private int rows;
	private int columns;

	public long iterationsLimited=0;
	public long iterationsUnlimited=0;

	public boolean verbose=false;
	public boolean verbose2=false;
	
}
