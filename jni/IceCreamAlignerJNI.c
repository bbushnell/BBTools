#include <jni.h>
#include <stdlib.h>
#include "icecream_IceCreamAlignerJNI.h"

#define pointsMatch 1
#define pointsSub -1
#define pointsDel -2
#define pointsIns -2

void alignForwardJNI(const jint * query, const jint * ref, jint * retArray, 
		const jint qlen, const jint rlen, const jint minScore, const float minRatio) {
	
	const int arrayLength=rlen;
	const int arrayLength2=rlen+1;

	//Stack allocated; faster.
	int array1[arrayLength2];
	int array2[arrayLength2];

	int * prev=array1;
	int * next=array2;
	int maxScore=-999999;
	int maxQpos=-1;
	int maxRpos=-1;
	int iters=0;
	retArray[0]=-999999;

	const int minPassingScore=(int)(qlen*minRatio*pointsMatch);
	const int minPassingScore3=minPassingScore-qlen*pointsMatch;

	for(int i=0; i<=arrayLength-qlen; i++){next[i]=prev[i]=0;}
	for(int i=arrayLength-qlen, score=0; i<=arrayLength; i++, score+=pointsDel) {
		prev[i]=score;
		next[i]=0;
	}
	
	for(int qpos=0; qpos<qlen; qpos++){
		prev[0]=pointsIns*qpos;

		const int q=query[qpos];
		int maxScoreThisPass=-999999;
		const int remainingBases=(qlen-qpos);
		const int remainingPoints=remainingBases*pointsMatch;
		const int minViableScore=minPassingScore3-remainingPoints;
		
		for(int rpos=0, apos=1; rpos<rlen; rpos++, apos++){
			const int r=ref[rpos];
			const int match=(q==r);
			const int vScore=prev[apos]+pointsIns;
			const int dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

			//Should be branchless conditional moves
			int score=(dScore>=vScore ? dScore : vScore);
			next[apos]=score;
		}
		
		for(int apos=1; apos<arrayLength2; apos++){
			const int hScore=next[apos-1]+pointsDel;

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

		int * temp=prev;
		prev=next;
		next=temp;
	}

	maxScore=-999999;
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


void alignForward16JNI(const jshort * query, const jshort * ref, jint * retArray, 
		const jshort qlen, const jshort rlen, const jshort minScore, const float minRatio) {
	
	const jshort arrayLength=rlen;
	const jshort arrayLength2=rlen+1;

	//Stack allocated; faster.
	jshort array1[arrayLength2];
	jshort array2[arrayLength2];
	
	jshort * prev=array1;
	jshort * next=array2;
	jshort maxScore=-32000;
	jshort maxQpos=-1;
	jshort maxRpos=-1;
	jint iters=0;
	retArray[0]=-32000;

	const jshort minPassingScore=(int)(qlen*minRatio*pointsMatch);
	const jshort minPassingScore3=minPassingScore-qlen*pointsMatch;

	for(jshort i=0; i<=arrayLength-qlen; i++){next[i]=prev[i]=0;}
	for(jshort i=arrayLength-qlen, score=0; i<=arrayLength; i++, score+=pointsDel) {
		prev[i]=score;
		next[i]=0;
	}
	
	for(jshort qpos=0; qpos<qlen; qpos++){
		prev[0]=pointsIns*qpos;

		const jshort q=query[qpos];
		jshort maxScoreThisPass=-32000;
		const jshort remainingBases=(qlen-qpos);
		const jshort remainingPoints=remainingBases*pointsMatch;
		const jshort minViableScore=minPassingScore3-remainingPoints;
		
//		for(jshort rpos=0, apos=1; rpos<rlen; rpos++, apos++){
//			const jshort r=ref[rpos];
//			const jshort match=(q==r);
//			const jshort vScore=prev[apos]+pointsIns;
//			const jshort hScore=next[apos-1]+pointsDel;
//			const jshort dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
//
//			//Should be branchless conditional moves
//			jshort score=(dScore>=vScore ? dScore : vScore);
//			score=(hScore>score ? hScore : score);
//			next[apos]=score;
//			maxScoreThisPass=(score>maxScoreThisPass ? score : maxScoreThisPass);
//		}
		
		for(jshort rpos=0, apos=1; rpos<rlen; rpos++, apos++){
			const jshort r=ref[rpos];
			const jshort match=(q==r);
			const jshort vScore=prev[apos]+pointsIns;
			const jshort dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

			//Should be branchless conditional moves
			jshort score=(dScore>=vScore ? dScore : vScore);
			next[apos]=score;
		}
		
		for(jshort apos=1; apos<arrayLength2; apos++){
			const jshort hScore=next[apos-1]+pointsDel;

			//Should be branchless conditional moves
			jshort score=next[apos];
			score=(hScore>score ? hScore : score);
			next[apos]=score;
			maxScoreThisPass=(score>maxScoreThisPass ? score : maxScoreThisPass);
		}
		iters+=arrayLength;
		
		//Aggressive early exit
		if(maxScoreThisPass<minScore){return;}

		//Safe early exit
		if(maxScoreThisPass<minViableScore){return;}

		jshort * temp=prev;
		prev=next;
		next=temp;
	}

	maxScore=-32000;
	for(jshort apos=1; apos<arrayLength2; apos++){//Grab high score from last iteration
		jshort score=prev[apos];
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


void alignForwardShortJNI(jint * query, jint * ref, jint * retArray, const jint qlen, const jint rlen) {

	const jint arrayLength=qlen;
	const jint arrayLength2=qlen+1;
	
	//Stack allocated; faster.
	jint array1[arrayLength2];
	jint array2[arrayLength2];
	
	//jint * array1=(jint *)malloc(4*arrayLength2);
	//jint * array2=(jint *)malloc(4*arrayLength2);
	
	jint * prev=array1;
	jint * next=array2;
	jint maxScore=-999999;
	jint maxQpos=-1;
	jint maxRpos=-1;
	jint itersShort=0;

	for(jint i=0; i<arrayLength2; i++) {
		prev[i]=pointsIns*i;
		next[i]=0;
	}

	for(jint rpos=0; rpos<rlen; rpos++){
		if(rlen-rpos<arrayLength){prev[0]=next[0]+pointsDel;}

		const jint r=ref[rpos];
		
//		//Inner loop
//		for(jint qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
//			const jint q=query[qpos];
//			const jboolean match=(q==r);
//			const jint vScore=prev[apos]+pointsIns;
//			const jint hScore=next[apos-1]+pointsDel;
//			const jint dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
//
//			score=(dScore>=vScore ? dScore : vScore);
//			score=(hScore>score ? hScore : score);
//
//			next[apos]=score;
//		}
		
		//Inner DV loop
		for(jint qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
			const jint q=query[qpos];
			const jboolean match=(q==r);
			const jint vScore=prev[apos]+pointsIns;
			const jint dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

			const jint score=(dScore>=vScore ? dScore : vScore);

			next[apos]=score;
		}
		
		//Inner I loop
		for(jint qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
			const jint hScore=next[apos-1]+pointsDel;

			jint score=next[apos];
			score=(hScore>score ? hScore : score);

			next[apos]=score;
		}
		
		itersShort+=arrayLength;
		
		jint score=next[arrayLength];
		if(score>=maxScore){
			maxScore=score;
			maxQpos=arrayLength-1;
			maxRpos=rpos;
		}

		jint * temp=prev;
		prev=next;
		next=temp;
	}
	retArray[0]=maxScore;
	retArray[1]=maxQpos;
	retArray[2]=maxRpos;
	retArray[3]=itersShort;
	
	//free(array1);
	//free(array2);
}


void alignForwardShort16JNI(jshort * query, jshort * ref, jint * retArray, const jshort qlen, const jshort rlen) {

	const jshort arrayLength=qlen;
	const jshort arrayLength2=qlen+1;
	
	//Stack allocated; faster.
	jshort array1[arrayLength2];
	jshort array2[arrayLength2];
	
	//jshort * array1=(jshort *)malloc(4*arrayLength2);
	//jshort * array2=(jshort *)malloc(4*arrayLength2);
	
	jshort * prev=array1;
	jshort * next=array2;
	jshort maxScore=-32000;
	jshort maxQpos=-1;
	jshort maxRpos=-1;
	jint itersShort=0;

	for(jshort i=0; i<arrayLength2; i++) {
		prev[i]=pointsIns*i;
		next[i]=0;
	}

	for(jshort rpos=0; rpos<rlen; rpos++){
		if(rlen-rpos<arrayLength){prev[0]=next[0]+pointsDel;}

		const jshort r=ref[rpos];
		jshort score=0;
		
//		//Inner loop
//		for(jshort qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
//			const jshort q=query[qpos];
//			const jboolean match=(q==r);
//			const jshort vScore=prev[apos]+pointsIns;
//			const jshort hScore=next[apos-1]+pointsDel;
//			const jshort dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
//
//			score=(dScore>=vScore ? dScore : vScore);
//			score=(hScore>score ? hScore : score);
//
//			next[apos]=score;
//		}
		
		//Inner DV loop
		for(jshort qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
			const jshort q=query[qpos];
			const jboolean match=(q==r);
			const jshort vScore=prev[apos]+pointsIns;
			const jshort dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

			score=(dScore>=vScore ? dScore : vScore);

			next[apos]=score;
		}
		
		//Inner I loop
		for(jshort qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
			const jshort hScore=next[apos-1]+pointsDel;

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

		jshort * temp=prev;
		prev=next;
		next=temp;
	}
	retArray[0]=maxScore;
	retArray[1]=maxQpos;
	retArray[2]=maxRpos;
	retArray[3]=itersShort;
	
	//free(array1);
	//free(array2);
}


void alignForwardShort16JNI_352(jshort * query, jshort * ref, jint * retArray, const jshort rlen) {

	const jshort arrayLength=352;
	const jshort arrayLength2=352+1;
	
	//Stack allocated; faster.
	jshort array1[arrayLength2];
	jshort array2[arrayLength2];
	
	//jshort * array1=(jshort *)malloc(4*arrayLength2);
	//jshort * array2=(jshort *)malloc(4*arrayLength2);
	
	jshort * prev=array1;
	jshort * next=array2;
	jshort maxScore=-16000;
	jshort maxQpos=-1;
	jshort maxRpos=-1;
	jint itersShort=0;

	for(jshort i=0; i<arrayLength2; i++) {
		prev[i]=pointsIns*i;
		next[i]=0;
	}

	for(jshort rpos=0; rpos<rlen; rpos++){
		if(rlen-rpos<arrayLength){prev[0]=next[0]+pointsDel;}

		const jshort r=ref[rpos];
		jshort score=0;
		
//		//Inner loop
//		for(jshort qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
//			const jshort q=query[qpos];
//			const jboolean match=(q==r);
//			const jshort vScore=prev[apos]+pointsIns;
//			const jshort hScore=next[apos-1]+pointsDel;
//			const jshort dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];
//
//			score=(dScore>=vScore ? dScore : vScore);
//			score=(hScore>score ? hScore : score);
//
//			next[apos]=score;
//		}
		
		//Inner DV loop
		for(jshort qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
			const jshort q=query[qpos];
			const jboolean match=(q==r);
			const jshort vScore=prev[apos]+pointsIns;
			const jshort dScore=(match ? pointsMatch : pointsSub)+prev[apos-1];

			score=(dScore>=vScore ? dScore : vScore);

			next[apos]=score;
		}
		
		//Inner I loop
		for(jshort qpos=0, apos=1; qpos<arrayLength; qpos++, apos++){
			const jshort hScore=next[apos-1]+pointsDel;

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

		jshort * temp=prev;
		prev=next;
		next=temp;
	}
	retArray[0]=maxScore;
	retArray[1]=maxQpos;
	retArray[2]=maxRpos;
	retArray[3]=itersShort;
	
	//free(array1);
	//free(array2);
}

JNIEXPORT void JNICALL Java_icecream_IceCreamAlignerJNI_alignForwardJNI
(JNIEnv *env, jobject obj, jintArray query, jintArray ref, jintArray rvector, jint qlen, jint rlen, jint minScore, jfloat minRatio){

	// Copy arrays from Java
	jint * jquery = (jint*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
	jint * jref = (jint*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
	jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);

	alignForwardJNI(jquery, jref, jrvector, qlen, rlen, minScore, minRatio);

	// Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
	(*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);
}

JNIEXPORT void JNICALL Java_icecream_IceCreamAlignerJNI_alignForward16JNI
(JNIEnv *env, jobject obj, jshortArray query, jshortArray ref, jintArray rvector, jshort qlen, jshort rlen, jshort minScore, jfloat minRatio){

	// Copy arrays from Java
	jshort * jquery = (jshort*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
	jshort * jref = (jshort*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
	jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);

	alignForward16JNI(jquery, jref, jrvector, qlen, rlen, minScore, minRatio);

	// Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
	(*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);
}

JNIEXPORT void JNICALL Java_icecream_IceCreamAlignerJNI_alignForwardShortJNI
(JNIEnv *env, jobject obj, jintArray query, jintArray ref, jintArray rvector, jint qlen, jint rlen){

	// Copy arrays from Java
	jint * jquery = (jint*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
	jint * jref = (jint*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
	jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);

	alignForwardShortJNI(jquery, jref, jrvector, qlen, rlen);

	// Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
	(*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);
}

JNIEXPORT void JNICALL Java_icecream_IceCreamAlignerJNI_alignForwardShort16JNI
(JNIEnv *env, jobject obj, jshortArray query, jshortArray ref, jintArray rvector, jshort qlen, jshort rlen){

	// Copy arrays from Java
	jshort * jquery = (jshort*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
	jshort * jref = (jshort*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
	jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);
	
	if(qlen==352){
		alignForwardShort16JNI_352(jquery, jref, jrvector, rlen);
	}else{
		alignForwardShort16JNI(jquery, jref, jrvector, qlen, rlen);
	}

	// Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
	(*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
	(*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);
}
