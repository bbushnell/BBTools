#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "align2_MultiStateAligner11tsJNI.h"

// C doesn't have min() or max() so we define our own
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; }) 

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; }) 

#define GAPLEN 128
#define GAPCOST max(1,GAPLEN/64)

#define GAPC '-'

#define MODE_MS 0
#define MODE_DEL 1
#define MODE_INS 2
#define MODE_SUB 3

#define AFFINE_ARRAYS 1

#define TIMEBITS 11
#define SCOREBITS (32-TIMEBITS)
#define MAX_TIME ((1<<TIMEBITS)-1)
#define MAX_SCORE (((1<<(SCOREBITS-1))-1)-2000)
#define MIN_SCORE (0-MAX_SCORE) //Keeps it 1 point above "BAD".

#define SCOREOFFSET TIMEBITS

#define TIMEMASK (~((-1)<<TIMEBITS))
#define SCOREMASK ((~((-1)<<SCOREBITS))<<SCOREOFFSET)

#define POINTS_NOREF 0 //default -110
#define POINTS_NOCALL 0
#define POINTS_MATCH 70 //default 50
#define POINTS_MATCH2 100 //Note:  Changing to 90 substantially reduces false positives
#define POINTS_COMPATIBLE 50
#define POINTS_SUB -127 //default -133
#define POINTS_SUBR -147 //increased penalty if prior match streak was at most 1 (I have no idea why this improves things)
#define POINTS_SUB2 -51 //default -47
#define POINTS_SUB3 -25
#define POINTS_MATCHSUB -10
#define POINTS_INS -395 //default -251
#define POINTS_INS2 -39 //default -61
#define POINTS_INS3 -23 //default -20
#define POINTS_INS4 -8 //default -20
#define POINTS_DEL -472 //default -239
#define POINTS_DEL2 -33 //default -30
#define POINTS_DEL3 -9 //default -7
#define POINTS_DEL4 -1 //default -1
#define POINTS_DEL5 -1 //default -1
#define POINTS_DEL_REF_N -10 //default -10
#define POINTS_GAP (0-GAPCOST) //default -10

#define TIMESLIP 4
#define MASK5 (TIMESLIP-1)

#define BARRIER_I1 2
#define BARRIER_D1 3

#define LIMIT_FOR_COST_3 5
#define LIMIT_FOR_COST_4 20
#define LIMIT_FOR_COST_5 80

#define BAD (MIN_SCORE-1)

#define POINTSoff_NOREF (POINTS_NOREF<<SCOREOFFSET)
#define POINTSoff_NOCALL (POINTS_NOCALL<<SCOREOFFSET)
#define POINTSoff_MATCH (POINTS_MATCH<<SCOREOFFSET)
#define POINTSoff_MATCH2 (POINTS_MATCH2<<SCOREOFFSET)
#define POINTSoff_COMPATIBLE (POINTS_COMPATIBLE<<SCOREOFFSET)
#define POINTSoff_SUB (POINTS_SUB<<SCOREOFFSET)
#define POINTSoff_SUBR (POINTS_SUBR<<SCOREOFFSET)
#define POINTSoff_SUB2 (POINTS_SUB2<<SCOREOFFSET)
#define POINTSoff_SUB3 (POINTS_SUB3<<SCOREOFFSET)
#define POINTSoff_MATCHSUB (POINTS_MATCHSUB<<SCOREOFFSET)
#define POINTSoff_INS (POINTS_INS<<SCOREOFFSET)
#define POINTSoff_INS2 (POINTS_INS2<<SCOREOFFSET)
#define POINTSoff_INS3 (POINTS_INS3<<SCOREOFFSET)
#define POINTSoff_INS4 (POINTS_INS4<<SCOREOFFSET)
#define POINTSoff_DEL (POINTS_DEL<<SCOREOFFSET)
#define POINTSoff_DEL2 (POINTS_DEL2<<SCOREOFFSET)
#define POINTSoff_DEL3 (POINTS_DEL3<<SCOREOFFSET)
#define POINTSoff_DEL4 (POINTS_DEL4<<SCOREOFFSET)
#define POINTSoff_DEL5 (POINTS_DEL5<<SCOREOFFSET)
#define POINTSoff_GAP (POINTS_GAP<<SCOREOFFSET)
#define POINTSoff_DEL_REF_N (POINTS_DEL_REF_N<<SCOREOFFSET)
#define BADoff (BAD<<SCOREOFFSET)
#define MAXoff_SCORE (MAX_SCORE<<SCOREOFFSET)
#define MINoff_SCORE (MIN_SCORE<<SCOREOFFSET)
	
void fillUnlimited(
	jbyte * read,
	jbyte * ref,
	jsize read_length,
	jsize ref_length,
	jint refStartLoc,
	jint refEndLoc,
	jint * result,
	jlong * iterationsUnlimited,
	jint * packed,
	jint * POINTSoff_SUB_ARRAY,
	jint * POINTSoff_INS_ARRAY,
	jint maxRows,
	jint maxColumns
	) {

	const jint rows=read_length;
	const jint columns=refEndLoc-refStartLoc+1;

	const jint maxGain=(read_length-1)*POINTSoff_MATCH2+POINTSoff_MATCH;
	const jint subfloor=0-2*maxGain;
	const jint BARRIER_I2=rows-BARRIER_I1, BARRIER_I2b=columns-1;
	const jint BARRIER_D2=rows-BARRIER_D1;
	
	const int sizeXY=(maxRows+1)*(maxColumns+1);
	const int idxMS=MODE_MS*sizeXY;
	const int idxDEL=MODE_DEL*sizeXY;
	const int idxINS=MODE_INS*sizeXY;

	//temporary, for finding a bug
	if(rows>maxRows || columns>maxColumns){
		printf("error\n"); exit(0);
	}

	for(int row=1; row<=rows; row++){
		const int tmp1=(row-1)*(maxColumns+1);
		const int tmp2=(row)*(maxColumns+1);
		for(int col=1; col<=columns; col++){
			(*iterationsUnlimited)++;

			const jbyte call0=(row<2 ? (jbyte)'?' : read[row-2]);
			const jbyte call1=read[row-1];
			const jbyte ref0=(col<2 ? (jbyte)'!' : ref[refStartLoc+col-2]);
			const jbyte ref1=ref[refStartLoc+col-1];

			const jboolean match=(call1==ref1 && ref1!='N');
			const jboolean prevMatch=(call0==ref0 && ref0!='N');
			
			const jboolean gap=(ref1==GAPC);

			if(gap){
				packed[idxMS+tmp2+col]=subfloor;
			}else{//Calculate match and sub scores

				const jint scoreFromDiag=packed[idxMS+tmp1+col-1]&SCOREMASK;
				const jint scoreFromDel=packed[idxDEL+tmp1+col-1]&SCOREMASK;
				const jint scoreFromIns=packed[idxINS+tmp1+col-1]&SCOREMASK;
				const jint streak=(packed[idxMS+tmp1+col-1]&TIMEMASK);

				{//Calculate match/sub score
					
					if(match){

						const jint scoreMS=scoreFromDiag+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
						const jint scoreD=scoreFromDel+POINTSoff_MATCH;
						const jint scoreI=scoreFromIns+POINTSoff_MATCH;

						jint score;
						jint time;
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

						if(time>MAX_TIME){time=MAX_TIME-MASK5;}
						packed[idxMS+tmp2+col]=(score|time);
						
					}else{
						
						jint scoreMS;
						if(ref1!='N' && call1!='N'){
							scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) : 
								POINTSoff_SUB_ARRAY[streak+1]);
						}else{
							scoreMS=scoreFromDiag+POINTSoff_NOCALL;
						}

						const jint scoreD=scoreFromDel+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
						const jint scoreI=scoreFromIns+POINTSoff_SUB;

						jint score;
						jint time;
						jbyte prevState;
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
						packed[idxMS+tmp2+col]=(score|time);
					}
				}
			}

			if(row<BARRIER_D1 || row>BARRIER_D2){
				packed[idxDEL+tmp2+col]=subfloor;
			}else{//Calculate DEL score
						
				const jint streak=packed[idxDEL+tmp2+col-1]&TIMEMASK;
				
				const jint scoreFromDiag=packed[idxMS+tmp2+col-1]&SCOREMASK;
				const jint scoreFromDel=packed[idxDEL+tmp2+col-1]&SCOREMASK;

				jint scoreMS=scoreFromDiag+POINTSoff_DEL;
				jint scoreD=scoreFromDel+(streak==0 ? POINTSoff_DEL : 
					streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
						streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 : 
							streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 : 
								((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));

				if(ref1=='N'){
					scoreMS+=POINTSoff_DEL_REF_N;
					scoreD+=POINTSoff_DEL_REF_N;
				}else if(gap){
					scoreMS+=POINTSoff_GAP;
					scoreD+=POINTSoff_GAP;
				}
				
				jint score;
				jint time;
				jbyte prevState;
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
				packed[idxDEL+tmp2+col]=(score|time);
			}

			//Calculate INS score
			if(gap || (row<BARRIER_I1 && col>1) || (row>BARRIER_I2 && col<BARRIER_I2b)){
				packed[idxINS+tmp2+col]=subfloor;
			}else{//Calculate INS score
				
				const jint streak=packed[idxINS+tmp1+col]&TIMEMASK;

				const jint scoreFromDiag=packed[idxMS+tmp1+col]&SCOREMASK;
				const jint scoreFromIns=packed[idxINS+tmp1+col]&SCOREMASK;

				const jint scoreMS=scoreFromDiag+POINTSoff_INS;
				const jint scoreI=scoreFromIns+POINTSoff_INS_ARRAY[streak+1];

				jint score;
				jint time;
				jbyte prevState;
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
				packed[idxINS+tmp2+col]=(score|time);
			}
		}
	}

	jint maxCol=-1;
	jint maxState=-1;
	jint maxScore=INT_MIN;
	
	const int tmp=rows*(maxColumns+1);
	for(int state=0; state<3; state++){
		for(int col=1; col<=columns; col++){
			const int x=packed[(state)*sizeXY+tmp+col]&SCOREMASK;
			if(x>maxScore){
				maxScore=x;
				maxCol=col;
				maxState=state;
			}
		}
	}
	maxScore>>=SCOREOFFSET;

	result[0]=rows;
	result[1]=maxCol;
	result[2]=maxState;
	result[3]=maxScore;
	return;
}

jint calcDelScoreOffset(jint len){
        if(len<=0){return 0;}
        jint score=POINTSoff_DEL;

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

jint calcInsScoreOffset(jint len,
			jint * POINTSoff_INS_ARRAY_C
			){
        if(len<=0){return 0;}
        if(AFFINE_ARRAYS){
                return POINTSoff_INS_ARRAY_C[len];
        }else{
                jint score=POINTSoff_INS;
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

void fillLimitedX(
	jbyte * read,
	jbyte * ref,
	jsize read_length,
	jsize ref_length,
	jint refStartLoc,
	jint refEndLoc,
	jint minScore,
	jint * result,
	jlong * iterationsLimited,
	jint * packed,
	jint * POINTSoff_SUB_ARRAY,
	jint * POINTSoff_INS_ARRAY,
	jint maxRows,
	jint maxColumns,
	jint bandwidth,
	jfloat bandwidthRatio,
	jint * vertLimit,
	jint * horizLimit,
	jbyte * baseToNumber,
	jint * POINTSoff_INS_ARRAY_C
	) {

	const jint rows=read_length;
	const jint columns=refEndLoc-refStartLoc+1;

	const int sizeXY=(maxRows+1)*(maxColumns+1);
	const int idxMS=MODE_MS*sizeXY;
	const int idxDEL=MODE_DEL*sizeXY;
	const int idxINS=MODE_INS*sizeXY;

	const jint halfband=(bandwidth<1 && bandwidthRatio<=0) ? 0 : 
		max(min(bandwidth<1 ? 9999999 : bandwidth, bandwidthRatio<=0 ? 9999999 : 8+(jint)(rows*bandwidthRatio)), (columns-rows+8))/2;
	
	const jint BARRIER_I2=rows-BARRIER_I1, BARRIER_I2b=columns-1;
	const jint BARRIER_D2=rows-BARRIER_D1;

	const int tmp=rows*(maxColumns+1);
	for(int x=0; x<3; x++){
		for(int i=1; i<columns+1; i++) {
			packed[x*sizeXY+tmp+i]=BADoff;
		}
	}
	
	jint minGoodCol=1;
	jint maxGoodCol=columns;
	
	const jint minScore_off=(minScore<<SCOREOFFSET);
	const jint maxGain=(read_length-1)*POINTSoff_MATCH2+POINTSoff_MATCH;
	const jint floor=minScore_off-maxGain;
	const jint subfloor=floor-5*POINTSoff_MATCH2;

	vertLimit[rows]=minScore_off;
	jboolean prevDefined=0;
	for(int i=rows-1; i>=0; i--){
		jbyte c=read[i];
		//if(AminoAcid.isFullyDefined(c)){
		if(baseToNumber[c]>=0){
			vertLimit[i]=max(vertLimit[i+1]-(prevDefined ? POINTSoff_MATCH2 : POINTSoff_MATCH), floor);
			prevDefined=1;
		}else{
			vertLimit[i]=max(vertLimit[i+1]-POINTSoff_NOCALL, floor);
			prevDefined=0;
		}
	}

	horizLimit[columns]=minScore_off;
	prevDefined=0;
	for(int i=columns-1; i>=0; i--){
		jbyte c=ref[refStartLoc+i];
		if(baseToNumber[c]>=0){
			horizLimit[i]=max(horizLimit[i+1]-(prevDefined ? POINTSoff_MATCH2 : POINTSoff_MATCH), floor);
			prevDefined=1;
		}else{
			horizLimit[i]=max(horizLimit[i+1]-(prevDefined && c==GAPC ? POINTSoff_DEL : POINTSoff_NOREF), floor);
			prevDefined=0;
		}
	}

	for(int row=1; row<=rows; row++){
		const jint colStart=(halfband<1 ? minGoodCol : max(minGoodCol, row-halfband));
		const jint colStop=(halfband<1 ? maxGoodCol : min(maxGoodCol, row+halfband*2-1));
		
		minGoodCol=-1;
		maxGoodCol=-2;
		
		const jint vlimit=vertLimit[row];
		
		if(colStart<0 || colStop<colStart){break;}
		
		if(colStart>1){
			const int tmp3=(row)*(maxColumns+1)+(colStart-1);
			packed[idxMS+tmp3]=subfloor;
			packed[idxINS+tmp3]=subfloor;
			packed[idxDEL+tmp3]=subfloor;
		}

		const int tmp1=(row-1)*(maxColumns+1);
		const int tmp2=(row)*(maxColumns+1);
		for(int col=colStart; col<=columns; col++){
			const jbyte call0=(row<2 ? (jbyte)'?' : read[row-2]);
			const jbyte call1=read[row-1];
			const jbyte ref0=(col<2 ? (jbyte)'!' : ref[refStartLoc+col-2]);
			const jbyte ref1=ref[refStartLoc+col-1];
			
			const jboolean gap=(ref1==GAPC);

			const jboolean match=(call1==ref1 && ref1!='N');
			const jboolean prevMatch=(call0==ref0 && ref0!='N');
			
			(*iterationsLimited)++;
			const jint limit=max(vlimit, horizLimit[col]);
			const jint limit3=max(floor, (match ? limit-POINTSoff_MATCH2 : limit-POINTSoff_SUB3));

			const jint delNeeded=max(0, row-col-1);
			const jint insNeeded=max(0, (rows-row)-(columns-col)-1);

			const jint delPenalty=calcDelScoreOffset(delNeeded);
			const jint insPenalty=calcInsScoreOffset(insNeeded,POINTSoff_INS_ARRAY_C);
			
			const jint scoreFromDiag_MS=packed[idxMS+tmp1+col-1]&SCOREMASK;
			const jint scoreFromDel_MS=packed[idxDEL+tmp1+col-1]&SCOREMASK;
			const jint scoreFromIns_MS=packed[idxINS+tmp1+col-1]&SCOREMASK;
			
			const jint scoreFromDiag_DEL=packed[idxMS+tmp2+col-1]&SCOREMASK;
			const jint scoreFromDel_DEL=packed[idxDEL+tmp2+col-1]&SCOREMASK;

			const jint scoreFromDiag_INS=packed[idxMS+tmp1+col]&SCOREMASK;
			const jint scoreFromIns_INS=packed[idxINS+tmp1+col]&SCOREMASK;

			if(gap || (scoreFromDiag_MS<=limit3 && scoreFromDel_MS<=limit3 && scoreFromIns_MS<=limit3)){
				packed[idxMS+tmp2+col]=subfloor;
			}else{//Calculate match and sub scores
				const jint streak=(packed[idxMS+tmp1+col-1]&TIMEMASK);
				{//Calculate match/sub score
					jint score;
					jint time;
					jbyte prevState;
					
					if(match){
						const jint scoreMS=scoreFromDiag_MS+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
						const jint scoreD=scoreFromDel_MS+POINTSoff_MATCH;
						const jint scoreI=scoreFromIns_MS+POINTSoff_MATCH;
						
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
						jint scoreMS;
						if(ref1!='N' && call1!='N'){
							scoreMS=scoreFromDiag_MS+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) : 
								POINTSoff_SUB_ARRAY[streak+1]);
						}else{
							scoreMS=scoreFromDiag_MS+POINTSoff_NOCALL;
						}
						
						const jint scoreD=scoreFromDel_MS+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
						const jint scoreI=scoreFromIns_MS+POINTSoff_SUB;
						
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
					}
					
					jint limit2;
					if(delNeeded>0){
						limit2=limit-delPenalty;
					}else if(insNeeded>0){
						limit2=limit-insPenalty;
					}else{
						limit2=limit;
					}
					
					if(score>=limit2){
						maxGoodCol=col;
						if(minGoodCol<0){minGoodCol=col;}
					}else{
						score=subfloor;
					}
					
					if(time>MAX_TIME){time=MAX_TIME-MASK5;}
					packed[idxMS+tmp2+col]=(score|time);
				}
			}

			if((scoreFromDiag_DEL<=limit && scoreFromDel_DEL<=limit) || row<BARRIER_D1 || row>BARRIER_D2){
				packed[idxDEL+tmp2+col]=subfloor;
			}else{//Calculate DEL score
				const jint streak=packed[idxDEL+tmp2+col-1]&TIMEMASK;
				
				jint scoreMS=scoreFromDiag_DEL+POINTSoff_DEL;
				jint scoreD=scoreFromDel_DEL+(streak==0 ? POINTSoff_DEL : 
					streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
						streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 : 
							streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 : 
								((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));
				
				if(ref1=='N'){
					scoreMS+=POINTSoff_DEL_REF_N;
					scoreD+=POINTSoff_DEL_REF_N;
				}else if(gap){
					scoreMS+=POINTSoff_GAP;
					scoreD+=POINTSoff_GAP;
				}
				
				jint score;
				jint time;
				jbyte prevState;
				if(scoreMS>=scoreD){
					score=scoreMS;
					time=1;
					prevState=MODE_MS;
				}else{
					score=scoreD;
					time=streak+1;
					prevState=MODE_DEL;
				}
	
				jint limit2;
				if(insNeeded>0){
					limit2=limit-insPenalty;
				}else if(delNeeded>0){
					limit2=limit-calcDelScoreOffset(time+delNeeded)+calcDelScoreOffset(time);
				}else{
					limit2=limit;
				}
				
				if(score>=limit2){
					maxGoodCol=col;
					if(minGoodCol<0){minGoodCol=col;}
				}else{
					score=subfloor;
				}
				
				if(time>MAX_TIME){time=MAX_TIME-MASK5;}
				packed[idxDEL+tmp2+col]=(score|time);
			}

			if(gap || (scoreFromDiag_INS<=limit && scoreFromIns_INS<=limit) || (row<BARRIER_I1 && col>1) || (row>BARRIER_I2 && col<BARRIER_I2b)){
				packed[idxINS+tmp2+col]=subfloor;
			}else{//Calculate INS score
				const jint streak=packed[idxINS+tmp1+col]&TIMEMASK;
				
				const jint scoreMS=scoreFromDiag_INS+POINTSoff_INS;
				const jint scoreI=scoreFromIns_INS+POINTSoff_INS_ARRAY[streak+1];
				
				jint score;
				jint time;
				jbyte prevState;
				if(scoreMS>=scoreI){
					score=scoreMS;
					time=1;
					prevState=MODE_MS;
				}else{
					score=scoreI;
					time=streak+1;
					prevState=MODE_INS;
				}
				
				jint limit2;
				if(delNeeded>0){
					limit2=limit-delPenalty;
				}else if(insNeeded>0){
					limit2=limit-calcInsScoreOffset(time+insNeeded,POINTSoff_INS_ARRAY_C)+calcInsScoreOffset(time,POINTSoff_INS_ARRAY_C);
				}else{
					limit2=limit;
				}

				if(score>=limit2){
					maxGoodCol=col;
					if(minGoodCol<0){minGoodCol=col;}
				}else{
					score=subfloor;
				}
				
				if(time>MAX_TIME){time=MAX_TIME-MASK5;}
				packed[idxINS+tmp2+col]=(score|time);
			}
			
			if(col>=colStop){
				if(col>colStop && (maxGoodCol<col || halfband>0)){break;}
				if(row>1){
					const int tmp3=tmp1+(col+1);
					packed[idxMS+tmp3]=subfloor;
					packed[idxINS+tmp3]=subfloor;
					packed[idxDEL+tmp3]=subfloor;
				}
			}
		}
	}
	
	jint maxCol=-1;
	jint maxState=-1;
	jint maxScore=INT_MIN;
	
	//const int tmp=rows*(maxColumns+1);
	for(int state=0; state<3; state++){
		for(int col=1; col<=columns; col++){
			const int x=packed[(state)*sizeXY+tmp+col]&SCOREMASK;
			if(x>maxScore){
				maxScore=x;
				maxCol=col;
				maxState=state;
			}
		}
	}
	
	if(maxScore<minScore_off){
		result[0]=rows;
		result[1]=maxCol;
		result[2]=maxState;
		result[3]=maxScore;
		result[4]=1;
		return;
	}
	
	maxScore>>=SCOREOFFSET;
	result[0]=rows;
	result[1]=maxCol;
	result[2]=maxState;
	result[3]=maxScore;
	result[4]=0;
	return;
}
	

JNIEXPORT void JNICALL Java_align2_MultiStateAligner11tsJNI_fillUnlimitedJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray read,
							jbyteArray ref,
							jint refStartLoc,
							jint refEndLoc,
							jintArray result,
							jlongArray iterationsUnlimited,
							jintArray packed,
							jintArray POINTSoff_SUB_ARRAY,
							jintArray POINTSoff_INS_ARRAY,
							jint maxRows,
							jint maxColumns
                                                        ) {
  // Get the size of the read and the reference arrays
  jsize read_length=(*env)->GetArrayLength(env, read);
  jsize ref_length=(*env)->GetArrayLength(env, ref);

  // Copy arrays from Java
  jint * jPOINTSoff_INS_ARRAY=(jint*)(*env)->GetPrimitiveArrayCritical(env, POINTSoff_INS_ARRAY, NULL);
  jint * jPOINTSoff_SUB_ARRAY=(jint*)(*env)->GetPrimitiveArrayCritical(env, POINTSoff_SUB_ARRAY, NULL);
  jint * jpacked=(jint*)(*env)->GetPrimitiveArrayCritical(env, packed, NULL);
  jbyte * jread=(jbyte*)(*env)->GetPrimitiveArrayCritical(env, read, NULL);
  jbyte * jref=(jbyte*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
  jint * jresult=(jint*)(*env)->GetPrimitiveArrayCritical(env, result, NULL);
  jlong * jiterationsUnlimited=(jlong*)(*env)->GetPrimitiveArrayCritical(env, iterationsUnlimited, NULL);

  // Using pointers for variables that need to be passed back to Java so they can be updated by the called functions
  jlong * iterationsUnlimitedPointer=&jiterationsUnlimited[0];

  // Call the fillUnlimited function in C; the 4 return values will be in jresult[]
  fillUnlimited(jread,jref,read_length,ref_length,refStartLoc,refEndLoc,jresult,iterationsUnlimitedPointer,jpacked,jPOINTSoff_SUB_ARRAY,jPOINTSoff_INS_ARRAY,maxRows,maxColumns);

  // Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
  (*env)->ReleasePrimitiveArrayCritical(env, result, jresult, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, iterationsUnlimited, jiterationsUnlimited, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, read, jread, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, packed, jpacked, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, POINTSoff_SUB_ARRAY, jPOINTSoff_SUB_ARRAY, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, POINTSoff_INS_ARRAY, jPOINTSoff_INS_ARRAY, JNI_ABORT);

  return;
}

JNIEXPORT void JNICALL Java_align2_MultiStateAligner11tsJNI_fillLimitedXJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray read,
							jbyteArray ref,
							jint refStartLoc,
							jint refEndLoc,
							jint minScore,
							jintArray result,
							jlongArray iterationsLimited,
							jintArray packed,
							jintArray POINTSoff_SUB_ARRAY,
							jintArray POINTSoff_INS_ARRAY,
							jint maxRows,
							jint maxColumns,
							jint bandwidth,
							jfloat bandwidthRatio,
							jintArray vertLimit,
							jintArray horizLimit,
                                                        jbyteArray baseToNumber,
							jintArray POINTSoff_INS_ARRAY_C
                                                        ) {
  // Get the size of the read and the reference arrays
  jsize read_length=(*env)->GetArrayLength(env, read);
  jsize ref_length=(*env)->GetArrayLength(env, ref);

  // Copy arrays from Java
  jint * jPOINTSoff_INS_ARRAY=(jint*)(*env)->GetPrimitiveArrayCritical(env, POINTSoff_INS_ARRAY, NULL);
  jint * jPOINTSoff_SUB_ARRAY=(jint*)(*env)->GetPrimitiveArrayCritical(env, POINTSoff_SUB_ARRAY, NULL);
  jint * jpacked=(jint*)(*env)->GetPrimitiveArrayCritical(env, packed, NULL);
  jbyte * jread=(jbyte*)(*env)->GetPrimitiveArrayCritical(env, read, NULL);
  jbyte * jref=(jbyte*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
  jint * jresult=(jint*)(*env)->GetPrimitiveArrayCritical(env, result, NULL);
  jlong * jiterationsLimited=(jlong*)(*env)->GetPrimitiveArrayCritical(env, iterationsLimited, NULL);
  jint * jvertLimit=(jint*)(*env)->GetPrimitiveArrayCritical(env, vertLimit, NULL);
  jint * jhorizLimit=(jint*)(*env)->GetPrimitiveArrayCritical(env, horizLimit, NULL);
  jbyte * jbaseToNumber=(jbyte*)(*env)->GetPrimitiveArrayCritical(env, baseToNumber, NULL);
  jint * jPOINTSoff_INS_ARRAY_C=(jint*)(*env)->GetPrimitiveArrayCritical(env, POINTSoff_INS_ARRAY_C, NULL);

  // Using pointers for variables that need to be passed back to Java so they can be updated by the called functions
  jlong * iterationsLimitedPointer=&jiterationsLimited[0];

  // Call the fillLimitedX function in C; the 5 return values will be in jresult[]
  fillLimitedX(jread,jref,read_length,ref_length,refStartLoc,refEndLoc,minScore,jresult,iterationsLimitedPointer,jpacked,jPOINTSoff_SUB_ARRAY,jPOINTSoff_INS_ARRAY,maxRows,maxColumns,bandwidth,bandwidthRatio,jvertLimit,jhorizLimit,jbaseToNumber,jPOINTSoff_INS_ARRAY_C);

  // Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
  (*env)->ReleasePrimitiveArrayCritical(env, result, jresult, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, iterationsLimited, jiterationsLimited, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, read, jread, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, packed, jpacked, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, POINTSoff_SUB_ARRAY, jPOINTSoff_SUB_ARRAY, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, POINTSoff_INS_ARRAY, jPOINTSoff_INS_ARRAY, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, vertLimit, jvertLimit, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, horizLimit, jhorizLimit, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, baseToNumber, jbaseToNumber, JNI_ABORT);
  (*env)->ReleasePrimitiveArrayCritical(env, POINTSoff_INS_ARRAY_C, jPOINTSoff_INS_ARRAY_C, JNI_ABORT);

  return;
}

