package icecream;

import aligner.AlignmentResult;

public final class IceCreamAlignerJava extends IceCreamAligner {

	IceCreamAlignerJava(){}
	
	/**
	 * @param query
	 * @param ref
	 * @param rstart
	 * @param rstop
	 * @param minScore Quit early if score drops below this
	 * @param minRatio Don't return results if max score is less than this fraction of max possible score
	 * @return
	 */
	@Override
	public AlignmentResult alignForward(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio) {
		final int arrayLength=rstop-rstart+1;

		//Stack allocated; faster.
		final int[] array1=new int[arrayLength+1], array2=new int[arrayLength+1];

		final int minPassingScore=(int)(query.length*minRatio*pointsMatch);
		final int minPassingScore3=minPassingScore-query.length*pointsMatch;

		int[] prev=array1, next=array2;
		int maxScore=-999999;
		int maxQpos=-1;
		int maxRpos=-1;

		for(int i=0; i<=arrayLength-query.length; i++){prev[i]=0;}
		for(int i=arrayLength-query.length, score=0; i<=arrayLength; i++, score+=pointsDel) {
			prev[i]=score;
		}

		int currentRstart=rstart;
		for(int qpos=0; qpos<query.length; qpos++){
			prev[0]=pointsIns*qpos;

			final byte q=query[qpos];
			int maxScoreThisPass=-9999;
			final int remainingBases=(query.length-qpos);
			final int remainingPoints=remainingBases*pointsMatch;
			final int minViableScore=minPassingScore3-remainingPoints;
			//			int minViableRstart=rstop+1;
			for(int rpos=currentRstart, apos=1+currentRstart-rstart; rpos<=rstop; rpos++, apos++){
				final byte r=ref[rpos];
				final boolean match=(q==r);
				final int vScore=prev[apos]+pointsIns;
				final int hScore=next[apos-1]+pointsDel;
				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

				//Slow branchy code
				//				final int score=Tools.max(vScore, hScore, dScore);
				//				next[apos]=score;
				//				if(score>=maxScoreThisPass){
				//					maxScoreThisPass=score;
				//					if(score>=maxScore){
				//						maxScore=score;
				//						maxQpos=qpos;
				//						maxRpos=rpos;
				//					}
				//				}

				//Should be branchless conditional moves
				int score=(dScore>=vScore ? dScore : vScore);
				score=(hScore>score ? hScore : score);
				next[apos]=score;
				maxScoreThisPass=(score>maxScoreThisPass ? score : maxScoreThisPass);

				//				minViableRstart=((score<minViableScore || rpos>minViableRstart) ? minViableRstart : rpos);
			}
			iters+=arrayLength;
			if(maxScoreThisPass<minScore){//Aggressive early exit
				return null;
			}

			{//Safe early exit
				if(maxScoreThisPass<minViableScore){return null;}
			}

			int[] temp=prev;
			prev=next;
			next=temp;
		}

		maxScore=-999999;
		for(int rpos=rstart, apos=1; rpos<=rstop; rpos++, apos++){//Grab high score from last iteration
			int score=prev[apos];
			if(score>=maxScore){
				maxScore=score;
				maxQpos=query.length;
				maxRpos=rpos;
			}
		}
		
		maxScore=(maxScore-pointsSub*query.length)/(pointsMatch-pointsSub);//Rescale 0 to length
		final float ratio=maxScore/(float)query.length;
		
		if(ratio<minRatio){return null;}
		return new AlignmentResult(maxScore, maxQpos, maxRpos, query.length, ref.length, rstart, rstop, ratio);
	}

	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @param rstop
	 * @param minScore Quit early if score drops below this
	 * @param minRatio Don't return results if max score is less than this fraction of max possible score
	 * @return
	 */
	@Override
	public AlignmentResult alignForwardShort(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio) {
		
		final int arrayLength=query.length;

		//Stack allocated; faster.
		final int[] array1=new int[arrayLength+1], array2=new int[arrayLength+1];

		final int minPassingScore=(int)(query.length*minRatio);
		final int minPassingScore3=minPassingScore-query.length;

		int[] prev=array1, next=array2;
		int maxScore=-999999;
		int maxQpos=-1;
		int maxRpos=-1;

		for(int i=0; i<prev.length; i++) {
			prev[i]=pointsIns*i;
		}

		for(int rpos=rstart; rpos<=rstop; rpos++){
			if(rstop-rpos<arrayLength){prev[0]=next[0]+pointsDel;}

			final byte r=ref[rpos];
			int score=0;
			
			//In Java, the unsplit loop is faster
			//-XX:+UseSuperWord seems to have no impact on speed.
			//Might be worth packing the query and ref into int arrays and trying again.
			for(int qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
				final byte q=query[qpos];
				final boolean match=(q==r);
				final int vScore=prev[apos]+pointsIns;
				final int hScore=next[apos-1]+pointsDel;
				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
				
				score=(dScore>=vScore ? dScore : vScore);
				score=(hScore>score ? hScore : score);

				next[apos]=score;
			}

//			//Inner DV loop
//			for(int qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
//				final int q=query[qpos];
//				final boolean match=(q==r);
//				final int vScore=prev[apos]+pointsIns;
//				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
//
//				score=(dScore>=vScore ? dScore : vScore);
//
//				next[apos]=score;
//			}
//			
//			//Inner I loop
//			for(int qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
//				final int hScore=next[apos-1]+pointsDel;
//
//				score=next[apos];
//				score=(hScore>score ? hScore : score);
//
//				next[apos]=score;
//			}
			
			itersShort+=arrayLength;

			if(score>=maxScore){
				maxScore=score;
				maxQpos=arrayLength-1;
				maxRpos=rpos;
			}

			int[] temp=prev;
			prev=next;
			next=temp;
		}

		maxScore=(maxScore-pointsSub*query.length)/(pointsMatch-pointsSub);//Rescale 0 to length
		final float ratio=maxScore/(float)query.length;
		
		if(ratio<minRatio){return null;}
		
		return new AlignmentResult(maxScore, maxQpos, maxRpos, query.length, ref.length, rstart, rstop, ratio);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	long iters(){return iters;}

	@Override
	long itersShort(){return itersShort;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long iters = 0;
	long itersShort = 0;
	
	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int pointsMatch = 1;
	public static final int pointsSub = -1;
	public static final int pointsDel = -2;
	public static final int pointsIns = -2;
	
}
