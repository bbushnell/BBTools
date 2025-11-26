package aligner;

import java.util.Arrays;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.Tools;
import stream.Read;
import stream.SiteScore;

/**
 * Based on MSA9ts, with transform scores tweaked for PacBio. */
public final class MultiStateAligner9PacBioAdapter3 {
	
	
	/**
	 * Constructs a new PacBio-optimized sequence aligner.
	 * Initializes dynamic programming matrices and precomputes scoring arrays
	 * for efficient alignment within the specified size constraints.
	 *
	 * @param maxRows_ Maximum number of query sequence positions (read length)
	 * @param maxColumns_ Maximum number of reference sequence positions
	 */
	public MultiStateAligner9PacBioAdapter3(int maxRows_, int maxColumns_){
//		assert(maxColumns_>=200);
//		assert(maxRows_>=200);
		maxRows=maxRows_;
		maxColumns=maxColumns_;
		packed=new int[3][maxRows+1][];
		
		for(int i=0; i<3; i++){
			packed[i][0]=new int[maxColumns+1];
			packed[i][1]=new int[maxColumns+1];
			packed[i][2]=new int[maxColumns+1];
			for(int j=3; j<maxRows+1; j++){
				packed[i][j]=packed[i][j-2];
			}
		}
		
		insScoreArray=new int[maxRows+1];
		delScoreArray=new int[maxColumns+1];
		matchScoreArray=new int[maxRows+1];
		subScoreArray=new int[maxRows+1];
		for(int i=0; i<insScoreArray.length; i++){
			insScoreArray[i]=(i==0 ? POINTSoff_INS :
				i<LIMIT_FOR_COST_3 ? POINTSoff_INS2 :
					i<LIMIT_FOR_COST_4 ? POINTSoff_INS3 : POINTSoff_INS4);
		}
		for(int i=0; i<delScoreArray.length; i++){
			delScoreArray[i]=(i==0 ? POINTSoff_DEL :
				i<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
					i<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 : POINTSoff_DEL4);
		}
		for(int i=0; i<matchScoreArray.length; i++){
			matchScoreArray[i]=(i==0 ? POINTSoff_MATCH2 : POINTSoff_MATCH);
		}
		for(int i=0; i<subScoreArray.length; i++){
			subScoreArray[i]=(i==0 ? POINTSoff_SUB : i<LIMIT_FOR_COST_3 ? POINTSoff_SUB2 : POINTSoff_SUB3);
		}
		
		vertLimit=new int[maxRows+1];
		horizLimit=new int[maxColumns+1];
		Arrays.fill(vertLimit, BADoff);
		Arrays.fill(horizLimit, BADoff);
		
		for(int matrix=0; matrix<packed.length; matrix++){
			for(int i=1; i<=3; i++){
				for(int j=0; j<packed[matrix][i].length; j++){
					packed[matrix][i][j]|=BADoff;
				}
			}
			
			//Initializes the rows; not needed.
			for(int i=1; i<=2; i++){
				
				int score=calcInsScoreOffset(i);
				packed[matrix][i][0]=score;
			}
		}
	}
	
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Will not fill areas that cannot match minScore */
	public final int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		return fillLimitedX(read, ref, refStartLoc, refEndLoc, minScore);
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
			
			{
				int score=calcInsScoreOffset(row);
				packed[0][row][0]=score;
				packed[1][row][0]=score;
				packed[2][row][0]=score;
			}
			
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

//				final boolean match=(read[row-1]==ref[refStartLoc+col-1]);
//				final boolean prevMatch=(row<2 || col<2 ? false : read[row-2]==ref[refStartLoc+col-2]);
				final boolean match=(call1==ref1);
				final boolean prevMatch=(call0==ref0);
				
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
				
				if((scoreFromDiag_MS<=limit3 && scoreFromDel_MS<=limit3 && scoreFromIns_MS<=limit3)){
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
							
						}else{
							
							int scoreMS;
							if(ref1!='N' && call1!='N'){
								scoreMS=scoreFromDiag_MS+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) : subScoreArray[streak]);
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
						
						assert(time<=MAX_TIME);//if(time>MAX_TIME){time=MAX_TIME-3;}
						assert(score>=MINoff_SCORE || score==BADoff) : "Score overflow - use MSA2 instead";
						assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
//						packed[MODE_MS][row][col]=(score|prevState|time);
						packed[MODE_MS][row][col]=(score|time);
						assert((score&SCOREMASK)==score);
//						assert((prevState&MODEMASK)==prevState);
						assert((time&TIMEMASK)==time);
					}
				}
				
				if((scoreFromDiag_DEL<=limit && scoreFromDel_DEL<=limit)){
//					assert((scoreFromDiag_DEL<=limit && scoreFromDel_DEL<=limit)) : scoreFromDiag_DEL+", "+row;
					packed[MODE_DEL][row][col]=subfloor;
				}else{//Calculate DEL score
							
					final int streak=packed[MODE_DEL][row][col-1]&TIMEMASK;
					
					int scoreMS=scoreFromDiag_DEL+POINTSoff_DEL;
					int scoreD=scoreFromDel_DEL+delScoreArray[streak];
//					int scoreI=scoreFromIns+POINTSoff_DEL;
					
					
					if(ref1=='N'){
						scoreMS+=POINTSoff_DEL_REF_N;
						scoreD+=POINTSoff_DEL_REF_N;
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
					
					assert(time<=MAX_TIME);//if(time>MAX_TIME){time=MAX_TIME-3;}
					assert(score>=MINoff_SCORE || score==BADoff) : "Score overflow - use MSA2 instead";
					assert(score<=MAXoff_SCORE) : "Score overflow - use MSA2 instead";
//					packed[MODE_DEL][row][col]=(score|prevState|time);
					packed[MODE_DEL][row][col]=(score|time);
					assert((score&SCOREMASK)==score);
//					assert((prevState&MODEMASK)==prevState);
					assert((time&TIMEMASK)==time);
				}
				
				if(scoreFromDiag_INS<=limit && scoreFromIns_INS<=limit){
					packed[MODE_INS][row][col]=subfloor;
				}else{//Calculate INS score
					
					final int streak=packed[MODE_INS][row-1][col]&TIMEMASK;
					
					int scoreMS=scoreFromDiag_INS+POINTSoff_INS;
//					int scoreD=scoreFromDel+POINTSoff_INS;
					int scoreI=scoreFromIns_INS+insScoreArray[streak];
					
					
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
					
					assert(time<=MAX_TIME);//if(time>MAX_TIME){time=MAX_TIME-3;}
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
	 * Does not require a min score (ie, same as old method) */
	private final int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc){
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;
		
		final int maxGain=(read.length-1)*POINTSoff_MATCH2+POINTSoff_MATCH;
		final int subfloor=0-2*maxGain;
		assert(subfloor>BADoff && subfloor*2>BADoff); //TODO: Actually, it needs to be substantially more.
		
		//temporary, for finding a bug
		if(rows>maxRows || columns>maxColumns){
			throw new RuntimeException("rows="+rows+", maxRows="+maxRows+", cols="+columns+", maxCols="+maxColumns+"\n"+new String(read)+"\n");
		}
		
		assert(rows<=maxRows) : "Check that values are in-bounds before calling this function: "+rows+", "+maxRows;
		assert(columns<=maxColumns) : "Check that values are in-bounds before calling this function: "+columns+", "+maxColumns;
		
		assert(refStartLoc>=0) : "Check that values are in-bounds before calling this function: "+refStartLoc;
		assert(refEndLoc<ref.length) : "Check that values are in-bounds before calling this function: "+refEndLoc+", "+ref.length;
		
		for(int row=1; row<=rows; row++){

			{
				int score=calcInsScoreOffset(row);
				packed[0][row][0]=score;
				packed[1][row][0]=score;
				packed[2][row][0]=score;
			}
			
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
				
				final boolean match=(call1==ref1);
				final boolean prevMatch=(call0==ref0);

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
							
							assert(time<=MAX_TIME);//if(time>MAX_TIME){time=MAX_TIME-3;}
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
								scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) : subScoreArray[streak]);
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
							
							assert(time<=MAX_TIME);//if(time>MAX_TIME){time=MAX_TIME-3;}
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
					int scoreD=scoreFromDel+delScoreArray[streak];
//					int scoreI=scoreFromIns+POINTSoff_DEL;
					
					if(ref1=='N'){
						scoreMS+=POINTSoff_DEL_REF_N;
						scoreD+=POINTSoff_DEL_REF_N;
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
					
					assert(time<=MAX_TIME);//if(time>MAX_TIME){time=MAX_TIME-3;}
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
					int scoreI=scoreFromIns+insScoreArray[streak];
					
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
					
					assert(time<=MAX_TIME);//if(time>MAX_TIME){time=MAX_TIME-3;}
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
	
	
	/** Generates the match string */
	public final byte[] traceback(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state){
//		assert(false);
		assert(refStartLoc<=refEndLoc) : refStartLoc+", "+refEndLoc;
		assert(row==rows);
		
		byte[] out=new byte[row+col-1]; //TODO if an out of bound crash occurs, try removing the "-1".
		int outPos=0;
		
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
				out[outPos]='D';
				
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
		
		return out2;
	}
	
	/** @return {score, bestRefStart, bestRefStop} */
	public final int[] score(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
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
	public final int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		
		int[] score;
		
		if(b-a>=maxColumns){
			System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
			assert(false) : refStartLoc+", "+refEndLoc;
			b=Tools.min(ref.length-1, a+maxColumns-1);
		}
		int[] max=fillLimited(read, ref, a, b, minScore);
		score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2]));
		
		return score;
	}
	
	

	/**
	 * Scores alignment without allowing insertions or deletions.
	 * Computes match/mismatch score for direct sequence comparison
	 * using chromosome data from the provided site score.
	 *
	 * @param read Query sequence as byte array
	 * @param ss Site score containing chromosome and position information
	 * @return Alignment score without indel penalties
	 */
	public final static int scoreNoIndels(byte[] read, SiteScore ss){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, ss.start, ss);
	}

	/**
	 * Scores alignment without allowing insertions or deletions.
	 * Computes match/mismatch score for direct sequence comparison
	 * at the specified chromosome and reference position.
	 *
	 * @param read Query sequence as byte array
	 * @param chrom Chromosome identifier
	 * @param refStart Starting position in reference sequence
	 * @return Alignment score without indel penalties
	 */
	public final static int scoreNoIndels(byte[] read, final int chrom, final int refStart){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, refStart, null);
	}
	
	/**
	 * Scores alignment without allowing insertions or deletions using quality scores.
	 * Incorporates base quality information to weight match/mismatch penalties
	 * and improve alignment scoring accuracy.
	 *
	 * @param read Query sequence as byte array
	 * @param ss Site score containing chromosome and position information
	 * @param baseScores Quality scores for each base in the read
	 * @return Quality-weighted alignment score without indel penalties
	 */
	public final static int scoreNoIndels(byte[] read, SiteScore ss, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, baseScores, ss.start, ss);
	}

	/**
	 * Scores alignment without allowing insertions or deletions using quality scores.
	 * Incorporates base quality information to weight match/mismatch penalties
	 * at the specified chromosome and reference position.
	 *
	 * @param read Query sequence as byte array
	 * @param chrom Chromosome identifier
	 * @param refStart Starting position in reference sequence
	 * @param baseScores Quality scores for each base in the read
	 * @return Quality-weighted alignment score without indel penalties
	 */
	public final static int scoreNoIndels(byte[] read, final int chrom, final int refStart, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, baseScores, refStart, null);
	}

	/**
	 * Scores alignment without allowing insertions or deletions.
	 * Computes match/mismatch score for direct sequence comparison
	 * and updates semiperfect status in the provided site score.
	 *
	 * @param read Query sequence as byte array
	 * @param ref Reference sequence as byte array
	 * @param refStart Starting position in reference sequence
	 * @param ss Site score object to update with semiperfect status (may be null)
	 * @return Alignment score without indel penalties
	 */
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
	

	/**
	 * Scores alignment without allowing insertions or deletions using quality scores.
	 * Incorporates base quality information and updates semiperfect alignment status.
	 * Most comprehensive no-indel scoring method with quality weighting.
	 *
	 * @param read Query sequence as byte array
	 * @param ref Reference sequence as byte array
	 * @param baseScores Quality scores for each base in the read
	 * @param refStart Starting position in reference sequence
	 * @param ss Site score object to update with semiperfect status (may be null)
	 * @return Quality-weighted alignment score without indel penalties
	 */
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
	
	
	/**
	 * Scores alignment without indels and generates match string representation.
	 * Creates detailed alignment visualization showing matches, mismatches, and no-calls
	 * while incorporating quality score weighting in the final score.
	 *
	 * @param read Query sequence as byte array
	 * @param ref Reference sequence as byte array
	 * @param baseScores Quality scores for each base in the read
	 * @param refStart Starting position in reference sequence
	 * @param matchReturn Single-element array to receive generated match string
	 * @return Quality-weighted alignment score without indel penalties
	 */
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
	
	
	/**
	 * Scores alignment without indels and generates match string representation.
	 * Creates detailed alignment visualization without quality score weighting
	 * for simpler alignment analysis and debugging.
	 *
	 * @param read Query sequence as byte array
	 * @param ref Reference sequence as byte array
	 * @param refStart Starting position in reference sequence
	 * @param matchReturn Single-element array to receive generated match string
	 * @return Alignment score without indel penalties
	 */
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
	
	/**
	 * Calculates theoretical maximum alignment score for perfect match.
	 * Assumes first base gets POINTS_MATCH and subsequent bases get POINTS_MATCH2
	 * to determine the highest possible score for a sequence of given length.
	 *
	 * @param numBases Length of sequence in bases
	 * @return Maximum possible alignment score for perfect match
	 */
	public static final int maxQuality(int numBases){
		return POINTS_MATCH+(numBases-1)*(POINTS_MATCH2);
	}
	
	/**
	 * Calculates theoretical maximum alignment score including quality bonuses.
	 * Incorporates base quality scores in addition to match scoring to determine
	 * the highest possible score achievable for the given sequence.
	 *
	 * @param baseScores Quality scores for each base
	 * @return Maximum possible quality-weighted alignment score
	 */
	public static final int maxQuality(byte[] baseScores){
		return POINTS_MATCH+(baseScores.length-1)*(POINTS_MATCH2)+Tools.sumInt(baseScores);
	}
	
	/**
	 * Calculates maximum alignment score allowing for single imperfection.
	 * Determines the highest score possible when one indel or substitution
	 * is allowed in an otherwise perfect alignment.
	 *
	 * @param numBases Length of sequence in bases
	 * @return Maximum alignment score with single imperfection
	 */
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
	
	/**
	 * Calculates maximum quality-weighted alignment score allowing single imperfection.
	 * Incorporates base quality scores while allowing one indel or substitution
	 * to determine realistic score thresholds for near-perfect alignments.
	 *
	 * @param baseScores Quality scores for each base
	 * @return Maximum quality-weighted alignment score with single imperfection
	 */
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
	
	/**
	 * Converts integer array to formatted string representation.
	 * Formats each element with consistent spacing for aligned output
	 * suitable for debugging and matrix visualization.
	 *
	 * @param a Integer array to format
	 * @return Formatted string with aligned columns
	 */
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
	
	/**
	 * Extracts and formats time values from packed integers.
	 * Isolates the time component from packed score/state/time values
	 * for debugging dynamic programming matrix state transitions.
	 *
	 * @param a Array of packed integers containing time information
	 * @return Formatted string showing time values only
	 */
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
	
	/**
	 * Extracts and formats score values from packed integers.
	 * Isolates the score component from packed score/state/time values
	 * with appropriate handling of overflow conditions.
	 *
	 * @param a Array of packed integers containing score information
	 * @return Formatted string showing score values only
	 */
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
	
	/**
	 * Converts byte array to formatted string representation.
	 * Formats each element with consistent spacing for debugging
	 * sequence and alignment data visualization.
	 *
	 * @param a Byte array to format
	 * @return Formatted string with aligned columns
	 */
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
	
	/**
	 * Extracts subsequence from byte array and converts to string.
	 * Creates human-readable representation of reference sequence region
	 * for debugging and alignment visualization.
	 *
	 * @param ref Reference sequence as byte array
	 * @param startLoc Starting position (inclusive)
	 * @param stopLoc Ending position (inclusive)
	 * @return String representation of specified subsequence
	 */
	public static final String toString(byte[] ref, int startLoc, int stopLoc){
		StringBuilder sb=new StringBuilder(stopLoc-startLoc+1);
		for(int i=startLoc; i<=stopLoc; i++){sb.append((char)ref[i]);}
		return sb.toString();
	}
	
	/**
	 * Calculates cumulative deletion penalty for specified length.
	 * Applies tiered penalty system with increasing costs for longer deletions
	 * to model realistic biological indel probabilities.
	 *
	 * @param len Length of deletion in bases
	 * @return Total deletion penalty score
	 */
	public static int calcDelScore(int len){
		if(len<=0){return 0;}
		int score=POINTS_DEL;
		
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
	
	/**
	 * Calculates deletion penalty with bit-shifted offset encoding.
	 * Internal version of deletion scoring used within dynamic programming
	 * matrix computations for efficient integer operations.
	 *
	 * @param len Length of deletion in bases
	 * @return Bit-shifted deletion penalty score
	 */
	private static int calcDelScoreOffset(int len){
		if(len<=0){return 0;}
		int score=POINTSoff_DEL;
		
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
	
	/**
	 * Calculates match bonus with bit-shifted offset encoding.
	 * Applies decreasing bonuses for longer match streaks using
	 * internal bit-shifted representation for matrix computations.
	 *
	 * @param len Length of match streak in bases
	 * @return Bit-shifted match bonus score
	 */
	private static int calcMatchScoreOffset(int len){
		if(len<=0){return 0;}
		int score=POINTSoff_MATCH;
		
		if(len>1){
			score+=(len-1)*POINTSoff_MATCH2;
		}
		return score;
	}
	
	/**
	 * Calculates substitution penalty with bit-shifted offset encoding.
	 * Applies tiered penalty system for consecutive substitutions
	 * using internal representation for dynamic programming efficiency.
	 *
	 * @param len Length of substitution streak in bases
	 * @return Bit-shifted substitution penalty score
	 */
	private static int calcSubScoreOffset(int len){
		if(len<=0){return 0;}
		int score=POINTSoff_SUB;
		
		if(len>LIMIT_FOR_COST_3){
			score+=(len-LIMIT_FOR_COST_3)*POINTSoff_SUB3;
			len=LIMIT_FOR_COST_3;
		}
		if(len>1){
			score+=(len-1)*POINTSoff_SUB2;
		}
		return score;
	}
	
	/**
	 * Calculates cumulative insertion penalty for specified length.
	 * Applies tiered penalty system with increasing costs for longer insertions
	 * to model realistic biological indel probabilities.
	 *
	 * @param len Length of insertion in bases
	 * @return Total insertion penalty score
	 */
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
	
	/**
	 * Calculates insertion penalty with bit-shifted offset encoding.
	 * Internal version of insertion scoring used within dynamic programming
	 * matrix computations for efficient integer operations.
	 *
	 * @param len Length of insertion in bases
	 * @return Bit-shifted insertion penalty score
	 */
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
	
	
	/** Maximum number of rows (read positions) in alignment matrix */
	private final int maxRows;
	/** Maximum number of columns (reference positions) in alignment matrix */
	private final int maxColumns;

	/**
	 * Three-dimensional matrix storing packed scores and states for match/insertion/deletion modes
	 */
	private final int[][][] packed;

	/** Vertical score limits for pruning low-scoring alignment paths */
	private final int[] vertLimit;
	/** Horizontal score limits for pruning low-scoring alignment paths */
	private final int[] horizLimit;

	/** Precomputed insertion penalty scores indexed by gap length */
	private final int[] insScoreArray;
	/** Precomputed deletion penalty scores indexed by gap length */
	private final int[] delScoreArray;
	/** Precomputed match bonus scores indexed by streak length */
	private final int[] matchScoreArray;
	/** Precomputed substitution penalty scores indexed by streak length */
	private final int[] subScoreArray;

	/** Returns formatted representation of vertical score limits for debugging */
	CharSequence showVertLimit(){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<=rows; i++){sb.append(vertLimit[i]>>SCOREOFFSET).append(",");}
		return sb;
	}
	/** Returns formatted representation of horizontal score limits for debugging */
	CharSequence showHorizLimit(){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<=columns; i++){sb.append(horizLimit[i]>>SCOREOFFSET).append(",");}
		return sb;
	}

//	public static final int MODEBITS=2;
	/** Number of bits allocated for time/streak information in packed integers */
	public static final int TIMEBITS=12;
	/** Number of bits allocated for score information in packed integers */
	public static final int SCOREBITS=32-TIMEBITS;
	/** Maximum time/streak value that can be stored in packed format */
	public static final int MAX_TIME=((1<<TIMEBITS)-1);
	/** Maximum score value that can be stored in packed format */
	public static final int MAX_SCORE=((1<<(SCOREBITS-1))-1)-2000;
	/** Minimum score value that can be stored in packed format */
	public static final int MIN_SCORE=0-MAX_SCORE; //Keeps it 1 point above "BAD".

//	public static final int MODEOFFSET=0; //Always zero.
//	public static final int TIMEOFFSET=0;
	/** Bit offset for score component in packed integers */
	public static final int SCOREOFFSET=TIMEBITS;

//	public static final int MODEMASK=~((-1)<<MODEBITS);
//	public static final int TIMEMASK=(~((-1)<<TIMEBITS))<<TIMEOFFSET;
	/** Bit mask for extracting time component from packed integers */
	public static final int TIMEMASK=~((-1)<<TIMEBITS);
	/** Bit mask for extracting score component from packed integers */
	public static final int SCOREMASK=(~((-1)<<SCOREBITS))<<SCOREOFFSET;
	
	/** Alignment state constant for match/substitution mode */
	private static final byte MODE_MS=0;
	/** Alignment state constant for deletion mode */
	private static final byte MODE_DEL=1;
	/** Alignment state constant for insertion mode */
	private static final byte MODE_INS=2;
	/** Alignment state constant for substitution mode */
	private static final byte MODE_SUB=3;
	
	/** Penalty score for alignment against missing reference sequence */
	public static final int POINTS_NOREF=-10;
	/** Penalty score for ambiguous base calls in query sequence */
	public static final int POINTS_NOCALL=-10;
	/** Bonus score for initial base match */
	public static final int POINTS_MATCH=90;
	/** Bonus score for subsequent consecutive base matches */
	public static final int POINTS_MATCH2=100; //Note:  Changing to 90 substantially reduces false positives
	/** Score for compatible base pairs (e.g., ambiguous codes) */
	public static final int POINTS_COMPATIBLE=50;
	/** Penalty score for initial substitution */
	public static final int POINTS_SUB=-143;
	/** Increased penalty for substitution after short match streak */
	public static final int POINTS_SUBR=-161; //increased penalty if prior match streak was at most 1
	/** Penalty score for second consecutive substitution */
	public static final int POINTS_SUB2=-54;
	/** Penalty score for third and subsequent consecutive substitutions */
	public static final int POINTS_SUB3=-35;
	/** Penalty for match immediately followed by substitution */
	public static final int POINTS_MATCHSUB=-10;
	/** Penalty score for initial insertion */
	public static final int POINTS_INS=-207;
	/** Penalty score for second consecutive insertion */
	public static final int POINTS_INS2=-51;
	/** Penalty score for third consecutive insertion */
	public static final int POINTS_INS3=-37;
	/** Penalty score for fourth and subsequent consecutive insertions */
	public static final int POINTS_INS4=-15;
	/** Penalty score for initial deletion */
	public static final int POINTS_DEL=-273;
	/** Penalty score for second consecutive deletion */
	public static final int POINTS_DEL2=-38;
	/** Penalty score for third consecutive deletion */
	public static final int POINTS_DEL3=-27;
	/** Penalty score for fourth and subsequent consecutive deletions */
	public static final int POINTS_DEL4=-15;
	/** Additional penalty for deletion against ambiguous reference base */
	public static final int POINTS_DEL_REF_N=-10;
	

	/** Indel length threshold for transitioning to third-tier penalty */
	public static final int LIMIT_FOR_COST_3=5;
	/** Indel length threshold for transitioning to fourth-tier penalty */
	public static final int LIMIT_FOR_COST_4=30;
	
	/** Sentinel value indicating invalid or impossible alignment score */
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
	public static final int POINTSoff_DEL_REF_N=(POINTS_DEL_REF_N<<SCOREOFFSET);
	public static final int BADoff=(BAD<<SCOREOFFSET);
	public static final int MAXoff_SCORE=MAX_SCORE<<SCOREOFFSET;
	public static final int MINoff_SCORE=MIN_SCORE<<SCOREOFFSET;
	
	/** Current number of rows in active alignment computation */
	private int rows;
	/** Current number of columns in active alignment computation */
	private int columns;

	/** Counter for alignment iterations with score-based pruning */
	public long iterationsLimited=0;
	/** Counter for alignment iterations without score-based pruning */
	public long iterationsUnlimited=0;

	/** Flag to enable verbose output for alignment debugging */
	public boolean verbose=false;
	/** Flag to enable detailed verbose output for alignment debugging */
	public boolean verbose2=false;
	
}
