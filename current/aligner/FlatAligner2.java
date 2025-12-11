package aligner;


/**
 * Sequence aligner using flatter scoring weights compared to standard aligners.
 * Implements dynamic programming alignment with optimized scoring constants
 * designed to reduce bias towards matches over gaps.
 * @author Brian Bushnell
 */
public final class FlatAligner2 {

	/** Default constructor for FlatAligner2. */
	public FlatAligner2(){}
	
	/**
	 * Performs forward alignment of query sequence against reference using dynamic programming.
	 * Uses stack-allocated arrays for improved performance and implements early termination
	 * when alignment score drops below minimum thresholds.
	 *
	 * @param query Query sequence bytes to align
	 * @param ref Reference sequence bytes
	 * @param rstart Starting position in reference sequence
	 * @param rstop Ending position in reference sequence
	 * @param minScore Minimum score threshold for early termination
	 * @param minRatio Minimum ratio of actual score to maximum possible score
	 * @return AlignmentResult containing score and positions, or null if thresholds not met
	 */
	public AlignmentResult alignForward(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio) {
//		if(ref.length<16300){return ica16.alignForward(query, ref, qstart, rstart, rstop, minScore, minRatio);}
		final int arrayLength=rstop-rstart+1;
		//		if(array1==null || arrayLength+1>array1.length){
		//			array1=new int[(int)((arrayLength+2)*1.2)];
		//			array2=new int[(int)((arrayLength+2)*1.2)];
		//		}

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
			//			currentRstart=minViableRstart;
			//			System.err.println("qPos="+qpos+", maxScoreThisPass="+maxScoreThisPass+", maxScore="+maxScore+", minScore="+minScore);
			//			System.err.println(Arrays.toString(prev));
			//			System.err.println(Arrays.toString(next));
			if(maxScoreThisPass<minScore){//Aggressive early exit
				//				System.err.print('.');
				//				System.err.println("qPos="+qpos+", maxScoreThisPass="+maxScoreThisPass+", maxScore="+maxScore+", minScore="+minScore);
				return null;
			}

			{//Safe early exit
				//				if((maxScoreThisPass+remaining+query.length)/2<minPassingScore){return null;}
				//				if((maxScoreThisPass+remaining+query.length)<minPassingScore2){return null;}
				if(maxScoreThisPass<minViableScore){return null;}
			}

			int[] temp=prev;
			prev=next;
			next=temp;
		}
		//		System.err.print('*');

		maxScore=-999999;
		for(int rpos=rstart, apos=1; rpos<=rstop; rpos++, apos++){//Grab high score from last iteration
			int score=prev[apos];
			if(score>=maxScore){
				maxScore=score;
				maxQpos=query.length;
				maxRpos=rpos;
			}
		}


		//		maxScore=(maxScore+query.length)/2;//Rescale 0 to length
		maxScore=(maxScore-pointsSub*query.length)/(pointsMatch-pointsSub);//Rescale 0 to length
		final float ratio=maxScore/(float)query.length;
		//		System.err.println("maxqPos="+maxQpos+", maxScore="+maxScore+", minScore="+(minRatio*query.length));
		//					if(minRatio*query.length>maxScore){return null;}
		if(ratio<minRatio){return null;}
		//		System.err.print('^');
		return new AlignmentResult(maxScore, maxQpos, maxRpos, query.length, ref.length, rstart, rstop, ratio);
	}

	/**
	 * Optimized alignment method for shorter sequences.
	 * Uses query length as array dimension rather than reference span for efficiency.
	 * Continues processing until minimum passing score is achieved.
	 *
	 * @param query Query sequence bytes to align
	 * @param ref Reference sequence bytes
	 * @param rstart Starting position in reference sequence
	 * @param rstop Ending position in reference sequence
	 * @param minRatio Minimum ratio of actual score to maximum possible score
	 * @return AlignmentResult containing score and positions, or null if ratio threshold not met
	 */
	public AlignmentResult alignForwardShort(final byte[] query, final byte[] ref, final int rstart, final int rstop,
			final float minRatio) {
//		if(ref.length<16300){return ica16.alignForwardShort(query, ref, qstart, rstart, rstop, minScore, minRatio);}

		final int arrayLength=query.length;
		//		if(array1==null || arrayLength+1>array1.length){
		//			array1=new int[(int)((arrayLength+2)*1.2)];
		//			array2=new int[(int)((arrayLength+2)*1.2)];
		//		}

		//Stack allocated; faster.
		final int[] array1=new int[arrayLength+1], array2=new int[arrayLength+1];

		final int minPassingScore=(int)(pointsMatch*query.length*minRatio);

		int[] prev=array1, next=array2;
		int maxScore=-1;
		int maxQpos=-1;
		int maxRpos=-1;

		for(int i=0; i<prev.length; i++) {
			prev[i]=pointsIns*i;
		}
		
		for(int rpos=rstart; rpos<=rstop && maxScore<minPassingScore; rpos++){
			if(rstop-rpos<arrayLength){prev[0]=next[0]+pointsDel;}

			final byte r=ref[rpos];
			int score=0;
			for(int qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
				final byte q=query[qpos];
				final boolean match=(q==r);
				final int vScore=prev[apos]+pointsIns;
				final int hScore=next[apos-1]+pointsDel;
				final int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

				//				score=Tools.max(vScore, hScore, dScore);
				score=(dScore>=vScore ? dScore : vScore);
				score=(hScore>score ? hScore : score);

				next[apos]=score;
			}
			itersShort+=arrayLength;

			if(score>=maxScore){
				maxScore=score;
				maxQpos=arrayLength-1;
				maxRpos=rpos;
			}

			//			System.err.println("qPos="+qpos+", maxScoreThisPass="+maxScoreThisPass+", maxScore="+maxScore+", minScore="+minScore);
			//			System.err.println(Arrays.toString(prev));
			//			System.err.println(Arrays.toString(next));
			//			if(maxScoreThisPass<minScore){//Aggressive early exit
			////				System.err.print('.');
			////				System.err.println("qPos="+qpos+", maxScoreThisPass="+maxScoreThisPass+", maxScore="+maxScore+", minScore="+minScore);
			//				return null;
			//			}

			//			{//Safe early exit
			//				int remaining=query.length-qpos;
			////				if((maxScoreThisPass+remaining+query.length)/2<minPassingScore){return null;}
			////				if((maxScoreThisPass+remaining+query.length)<minPassingScore2){return null;}
			//				if((maxScoreThisPass+remaining)<minPassingScore3){return null;}
			//			}

			int[] temp=prev;
			prev=next;
			next=temp;
		}
		//		System.err.print('*');

		//		maxScore=-999999;
		//		for(int rpos=rstart, apos=1; rpos<=rstop; rpos++, apos++){//Grab high score from last iteration
		//			int score=prev[apos];
		//			if(score>=maxScore){
		//				maxScore=score;
		//				maxQpos=query.length;
		//				maxRpos=rpos;
		//			}
		//		}
		final float ratio=maxScore/(float)(pointsMatch*query.length);
//		System.err.println(ratio+"\t"+maxScore+"\t"+query.length);
		//		System.err.println("maxqPos="+maxQpos+", maxScore="+maxScore+", minScore="+(minRatio*query.length));
		//					if(minRatio*query.length>maxScore){return null;}
		if(ratio<minRatio){return null;}
		//		System.err.print('^');
		return new AlignmentResult(maxScore, maxQpos, maxRpos, query.length, ref.length, rstart, rstop, ratio);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	long iters(){return iters;}
	long itersShort(){return itersShort;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long iters = 0;
	long itersShort = 0;
	
	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int pointsMatch = 10;
	public static final int pointsSub = -9;
	public static final int pointsDel = -11;
	public static final int pointsIns = -11;
	
}
