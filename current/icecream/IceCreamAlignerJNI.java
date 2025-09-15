package icecream;

import aligner.AlignmentResult;
import shared.KillSwitch;
import shared.Shared;

public final class IceCreamAlignerJNI extends IceCreamAligner {

	static {
		Shared.loadJNI();
	}

	IceCreamAlignerJNI(){}
	
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
	public AlignmentResult alignForward(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio) {
		
		final int qlen=query.length;
		final int rlen=rstop-rstart+1;
		final int[] retVec=KillSwitch.allocInt1D(4);
		
//		alignForwardJNI(qInt, rInt, retVec, qlen, rlen, minScore, minRatio);
//		alignForwardPseudo(qInt, rInt, retVec, qlen, rlen, minScore, minRatio);
		
		if((qlen+rlen+32)*2>Short.MAX_VALUE) {
			final int[] qInt=new int[query.length];
			final int[] rInt=new int[rlen];
			for(int i=0; i<query.length; i++){qInt[i]=query[i];}
			for(int i=0; i<rlen; i++){rInt[i]=ref[i+rstart];}

//			alignForwardPseudo(qInt, rInt, retVec, qlen, rlen, minScore, minRatio);
			alignForwardJNI(qInt, rInt, retVec, qlen, rlen, minScore, minRatio);
//			alignForwardShortJNI(qInt, rInt, retVec, qlen, rlen);
		}else{
			final short[] qInt=new short[query.length];
			final short[] rInt=new short[rlen];
			for(int i=0; i<query.length; i++){qInt[i]=query[i];}
			for(int i=0; i<rlen; i++){rInt[i]=ref[i+rstart];}
			
			alignForward16JNI(qInt, rInt, retVec, (short)qlen, (short)rlen, (short)minScore, minRatio);
//			alignForwardShort16JNI(qInt, rInt, retVec, (short)qlen, (short)rlen);
		}
		
		int maxScore=retVec[0];
		final int maxQpos=retVec[1];
		final int maxRpos=retVec[2]+rstart;
		final int jniIters=retVec[3];
		
		iters+=jniIters;
		
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
		
		final int qlen=query.length;
		final int rlen=rstop-rstart+1;
		final int[] retVec=KillSwitch.allocInt1D(4);
		
		if((qlen+rlen+32)*2>Short.MAX_VALUE) {
			final int[] qInt=new int[query.length];
			final int[] rInt=new int[rlen];
			for(int i=0; i<query.length; i++){qInt[i]=query[i];}
			for(int i=0; i<rlen; i++){rInt[i]=ref[i+rstart];}

//			alignForwardShortPseudo(qInt, rInt, retVec, qlen, rlen);
			alignForwardShortJNI(qInt, rInt, retVec, qlen, rlen);
		}else{
			final short[] qInt=new short[query.length];
			final short[] rInt=new short[rlen];
			for(int i=0; i<query.length; i++){qInt[i]=query[i];}
			for(int i=0; i<rlen; i++){rInt[i]=ref[i+rstart];}
			
			alignForwardShort16JNI(qInt, rInt, retVec, (short)qlen, (short)rlen);
		}
		int maxScore=retVec[0];
		final int maxQpos=retVec[1];
		final int maxRpos=retVec[2]+rstart;
		final int jniIters=retVec[3];
		
		itersShort+=jniIters;
		
		maxScore=(maxScore-pointsSub*query.length)/(pointsMatch-pointsSub);//Rescale 0 to length
		final float ratio=maxScore/(float)query.length;
		
		if(ratio<minRatio){return null;}
		
		return new AlignmentResult(maxScore, maxQpos, maxRpos, query.length, ref.length, rstart, rstop, ratio);
	}

	/*--------------------------------------------------------------*/
	/*----------------             JNI              ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void alignForwardPseudo(final int[] query, final int[] ref, final int[] retArray, 
			final int qlen, final int rlen, final int minScore, final float minRatio) {
		
		final int arrayLength=rlen;
		final int arrayLength2=rlen+1;

		//Stack allocated; faster.
		int[] array1=new int[arrayLength2];
		int[] array2=new int[arrayLength2];

		int[] prev=array1;
		int[] next=array2;
		int maxScore=-32000;
		int maxQpos=-1;
		int maxRpos=-1;
		int iters=0;

		final int minPassingScore=(int)(qlen*minRatio*pointsMatch);
		final int minPassingScore3=minPassingScore-qlen*pointsMatch;

		for(int i=0; i<=arrayLength-qlen; i++){prev[i]=0;}
		for(int i=arrayLength-qlen, score=0; i<=arrayLength; i++, score+=pointsDel) {
			prev[i]=score;
			next[i]=0;
		}
		
		for(int qpos=0; qpos<qlen; qpos++){
			prev[0]=pointsIns*qpos;

			final int q=query[qpos];
			int maxScoreThisPass=-32000;
			final int remainingBases=(qlen-qpos);
			final int remainingPoints=remainingBases*pointsMatch;
			final int minViableScore=minPassingScore3-remainingPoints;
			
//			for(int rpos=0, apos=1; rpos<rlen; rpos++, apos++){
//				final int r=ref[rpos];
//				final boolean match=(q==r);
//				final int vScore=prev[apos]+pointsIns;
//				final int hScore=next[apos-1]+pointsDel;
//				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
//
//				//Should be branchless conditional moves
//				int score=(dScore>=vScore ? dScore : vScore);
//				score=(hScore>score ? hScore : score);
//				next[apos]=score;
//				maxScoreThisPass=(score>maxScoreThisPass ? score : maxScoreThisPass);
//			}
			
			for(int rpos=0, apos=1; rpos<rlen; rpos++, apos++){
				final int r=ref[rpos];
				final boolean match=(q==r);
				final int vScore=prev[apos]+pointsIns;
				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

				//Should be branchless conditional moves
				int score=(dScore>=vScore ? dScore : vScore);
				next[apos]=score;
			}
			
			for(int apos=1; apos<arrayLength2; apos++){
				final int hScore=next[apos-1]+pointsDel;

				//Should be branchless conditional moves
				int score=next[apos];
				score=(hScore>score ? hScore : score);
				next[apos]=score;
				maxScoreThisPass=(score>maxScoreThisPass ? score : maxScoreThisPass);
			}
			iters+=arrayLength;
			
			//Aggressive early exit
			if(maxScoreThisPass<minScore){return;}

			//Safe early exit
			if(maxScoreThisPass<minViableScore){return;}

			int[] temp=prev;
			prev=next;
			next=temp;
		}

		maxScore=-32000;
		for(int apos=1; apos<arrayLength2; apos++){//Grab high score from last iteration
			int score=prev[apos];
			if(score>=maxScore){
				maxScore=score;
				maxQpos=qlen;
				maxRpos=apos-1;
			}
		}
		retArray[0]=maxScore;
		retArray[1]=maxQpos;
		retArray[2]=maxRpos;
		retArray[3]=iters;
	}
	
	private static void alignForwardShortPseudo(int[] query, int[] ref, int[] retArray, int qlen, int rlen) {
		
		final int arrayLength=qlen;
		final int arrayLength2=qlen+1;

		//Stack allocated; faster.
		int[] array1=new int[arrayLength2];
		int[] array2=new int[arrayLength2];

		int[] prev=array1;
		int[] next=array2;
		int maxScore=-999999;
		int maxQpos=-1;
		int maxRpos=-1;
		int itersShort=0;

		for(int i=0; i<arrayLength2; i++) {
			prev[i]=pointsIns*i;
			next[i]=0;//For C version initialization
		}

		for(int rpos=0; rpos<rlen; rpos++){
			if(rlen-rpos<arrayLength){prev[0]=next[0]+pointsDel;}
			
			final int r=ref[rpos];
			int score=0;
			
//			//Inner loop
//			for(int qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
//				final int q=query[qpos];
//				final boolean match=(q==r);
//				final int vScore=prev[apos]+pointsIns;
//				final int hScore=next[apos-1]+pointsDel;
//				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
//
//				score=(dScore>=vScore ? dScore : vScore);
//				score=(hScore>score ? hScore : score);
//
//				next[apos]=score;
//			}
			
			//Inner DV loop
			for(int qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
				final int q=query[qpos];
				final boolean match=(q==r);
				final int vScore=prev[apos]+pointsIns;
				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

				score=(dScore>=vScore ? dScore : vScore);

				next[apos]=score;
			}
			
			//Inner I loop
			for(int qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
				final int hScore=next[apos-1]+pointsDel;

				score=next[apos];
				score=(hScore>score ? hScore : score);

				next[apos]=score;
			}

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
		retArray[0]=maxScore;
		retArray[1]=maxQpos;
		retArray[2]=maxRpos;
		retArray[3]=itersShort;
	}

	private static native void alignForwardJNI(int[] query, int[] ref, int[] retArray, int qlen, int rlen, int minScore, float minRatio);
	private static native void alignForward16JNI(short[] query, short[] ref, int[] retArray, short qlen, short rlen, short minScore, float minRatio);
	private static native void alignForwardShortJNI(int[] query, int[] ref, int[] retArray, int qlen, int rlen);
	private static native void alignForwardShort16JNI(short[] query, short[] ref, int[] retArray, short qlen, short rlen);
	
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
