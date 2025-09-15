#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include "jgi_BBMergeOverlapper.h"

// C doesn't have min() or max() so we define our own
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; }) 

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; }) 

#define BAD_MULT 6
#define GOOD_MULT_1 8
#define GOOD_MULT_2 400

jint mid(int x, int y, int z){return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);}

//Looks fine on cursory inspection - BB Oct 27 2015
jint mateByOverlap(jbyte * a_bases, const jint a_bases_length, jbyte * b_bases, const jint b_bases_length, jbyte * a_quality, jbyte * b_quality, jfloat * aprob, jfloat * bprob, jint * rvector, jint minOverlap0, const jint minOverlap, const jint minInsert0, jint margin, const jint maxMismatches0, const jint maxMismatches, const jint minq) {

	minOverlap0=min(max(1, minOverlap0), minOverlap);
	margin=max(margin, 0);

	const jbyte *abases=a_bases, *bbases=b_bases;
	jbyte *aqual=NULL; 
	jbyte *bqual=NULL;
	if(a_quality!=NULL){
		aqual=a_quality;
	}
	if(b_quality!=NULL){
		bqual=b_quality;
	}
	const jint alen=a_bases_length, blen=b_bases_length;

	jint bestOverlap=-1;
	jint bestGood=-1;
	jint bestBad=maxMismatches0;

	jboolean ambig=0;
	const jint maxOverlap=alen+blen-max(minOverlap, minInsert0);

	const jfloat probCorrect[71] = 
	{0.000f, 0.251f, 0.369f, 0.499f, 0.602f, 0.684f, 0.749f, 0.800f, 0.842f, 0.874f, 0.900f, 0.921f, 0.937f, 0.950f, 0.960f, 0.968f,
	 0.975f, 0.980f, 0.984f, 0.987f, 0.990f, 0.992f, 0.994f, 0.995f, 0.996f, 0.997f, 0.997f, 0.998f, 0.998f, 0.999f, 0.999f, 0.999f,
	 0.999f, 0.999f, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

	if(aqual!=NULL && bqual!=NULL){
		for(jint i=0; i<alen; i++){aprob[i]=probCorrect[aqual[i]];}
		for(jint i=0; i<blen; i++){bprob[i]=probCorrect[bqual[i]];}
	}else{
		for(jint i=0; i<alen; i++){aprob[i]=0.98f;}
		for(jint i=0; i<blen; i++){bprob[i]=0.98f;}
	}
	
	const jfloat minprob=probCorrect[mid(1, minq, 41)];
	
	for(jint overlap=max(minOverlap0, 0); overlap<maxOverlap; overlap++){
		jint good=0, bad=0;
		jint istart=(overlap<=alen ? 0 : overlap-alen);
		jint jstart=(overlap<=alen ? alen-overlap : 0);
		{
			const jint iters=min(overlap-istart, min(blen-istart, alen-jstart));
			const jint imax=istart+iters;
			const jint badlim=bestBad+margin;

			for(jint i=istart, j=jstart; i<imax && bad<=badlim; i++, j++){
				const jbyte ca1=abases[j], cb1=bbases[i];
				const jfloat pc=aprob[j]*bprob[j];

				if(pc<=minprob){//do nothing
				}else if(ca1==cb1){good++;}
				else{bad++;}
			}
		}

		if(bad*2<good){
			if(good>minOverlap){//Candidate
				if(bad<=bestBad){
					if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
						if(bestBad-bad<margin){ambig=1;}
						bestOverlap=overlap;
						bestBad=bad;
						bestGood=good;
					}else if(bad==bestBad){
						ambig=1;
					}

					if(ambig && bestBad<margin){
						rvector[2]=bestBad;
						rvector[4]=(ambig ? 1 : 0);
						return -1;
					}
				}
			}else if(bad<margin){
				ambig=1;
				rvector[2]=bestBad;
				rvector[4]=(ambig ? 1 : 0);
				return -1;
			}else{
			}
		}
	}
	
	if(!ambig && bestBad>maxMismatches-margin){bestOverlap=-1;}
	
	rvector[2]=bestBad;
	rvector[4]=(ambig ? 1 : 0);
	
	return (bestOverlap<0 ? -1 : alen+blen-bestOverlap);
}
//Fixed - BB Oct 26 2015
jfloat findBestRatio(jbyte * a_bases, const jint a_bases_length, jbyte * b_bases, const jint b_bases_length,
		const jint minOverlap0, const jint minOverlap, const jint minInsert, const jfloat maxRatio, const jfloat offset, const jfloat gIncr, const jfloat bIncr) {
	const jbyte *abases=a_bases, *bbases=b_bases;
	const jint alen=a_bases_length, blen=b_bases_length;
	
	jfloat bestRatio=maxRatio+0.0001f;
	const jint maxOverlap=alen+blen-max(minOverlap, minInsert);
//	const jfloat altBadlimit=max(maxRatio, 0.07f)*2f*alen+1;
	const jfloat halfmax=maxRatio*0.5f;
	const jbyte N='N';
	
	const jint largestInsertToTest=(alen+blen-minOverlap);
	const jint smallestInsertToTest=minInsert;
	for(jint insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
		const jint istart=(insert<=blen ? 0 : insert-blen);
		const jint jstart=(insert>=blen ? 0 : blen-insert);	
		const jint overlapLength=min(alen-istart, min(blen-jstart, insert));

//		const jfloat badlimit=(min(altBadlimit, bestRatio*overlapLength));
		const jfloat badlimit=bestRatio*overlapLength;
		jfloat good=0, bad=0;

		const jint imax=istart+overlapLength;
		for(jint i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
			const jbyte ca=abases[i], cb=bbases[j];
			
			if(ca==cb){
				if(ca!=N){good+=gIncr;}
			}else{bad+=bIncr;}
		}

		if(bad<=badlimit){
			if(bad==0 && good>minOverlap0 && good<minOverlap){
				return 100.0f;
			}

			jfloat ratio=(bad+offset)/overlapLength;

			if(ratio<bestRatio){
				bestRatio=ratio;
				if(good>=minOverlap && ratio<halfmax){return bestRatio;}
			}
		}
	}
	
	return bestRatio;
}
//Fixed - BB Oct 27 2015
//TODO: remove a_quality and b_quality
jfloat findBestRatio_WithQualities(jbyte * a_bases, const jint a_bases_length, jbyte * b_bases, const jint b_bases_length, 
		jbyte * a_quality, jbyte * b_quality, //TODO:  Not needed
		jfloat * aprob, jfloat * bprob, 
		const jint minOverlap0, const jint minOverlap, const jint minInsert, const jfloat maxRatio, const jfloat offset) {
	const jbyte *abases=a_bases, *bbases=b_bases;
	const jint alen=a_bases_length, blen=b_bases_length;
	
	jfloat bestRatio=maxRatio+0.0001f;
//	const jfloat altBadlimit=max(maxRatio, 0.07f)*2f*alen+1;
	const jfloat halfmax=maxRatio*0.5f;
	
	
	const jint largestInsertToTest=(alen+blen-minOverlap); //TODO: test speed with minOverlap0
	const jint smallestInsertToTest=minInsert;
	for(jint insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
	
	
		const jint istart=(insert<=blen ? 0 : insert-blen);
		const jint jstart=(insert>=blen ? 0 : blen-insert);
		const jint overlapLength=min(alen-istart, min(blen-jstart, insert));

//		const jfloat badlimit=(min(altBadlimit, bestRatio*overlapLength));
		const jfloat badlimit=bestRatio*overlapLength;
		jfloat good=0, bad=0;

		const jint imax=istart+overlapLength;
		for(jint i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
			const jbyte ca=abases[i], cb=bbases[j];
			const jfloat x=aprob[i]*bprob[j];

			if(ca==cb){good+=x;}
			else{bad+=x;}
		}

		if(bad<=badlimit){
			if(bad==0 && good>minOverlap0 && good<minOverlap){
				return 100.0f;
			}

			jfloat ratio=(bad+offset)/overlapLength;

			if(ratio<bestRatio){
				bestRatio=ratio;
				if(good>=minOverlap && ratio<halfmax){return bestRatio;}
			}
		}
	}
	
	return bestRatio;
}
//Fixed - BB Oct 27 2015
jint mateByOverlapRatio_WithQualities(jbyte * a_bases, const jint a_bases_length, jbyte * b_bases, const jint b_bases_length, 
		jbyte * a_quality, jbyte * b_quality,
		jfloat * aprob, jfloat * bprob, 
		jint * rvector, jint minOverlap0, jint minOverlap, jint minInsert0, jint minInsert, jfloat maxRatio, const jfloat margin, const jfloat offset) {
	minOverlap=max(4, max(minOverlap0, minOverlap));
	minOverlap0=mid(4, minOverlap0, minOverlap);
	
        const jbyte *abases=a_bases, *bbases=b_bases, *aqual=a_quality, *bqual=b_quality;
        const jint alen=a_bases_length, blen=b_bases_length;
	const jint minLength=min(alen, blen);
	
	{
		const jfloat probCorrect[71] = 
		{0.000f, 0.251f, 0.369f, 0.499f, 0.602f, 0.684f, 0.749f, 0.800f, 0.842f, 0.874f, 0.900f, 0.921f, 0.937f, 0.950f, 0.960f, 0.968f,
	 	0.975f, 0.980f, 0.984f, 0.987f, 0.990f, 0.992f, 0.994f, 0.995f, 0.996f, 0.997f, 0.997f, 0.998f, 0.998f, 0.999f, 0.999f, 0.999f,
	 	0.999f, 0.999f, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

		for(jint i=0; i<alen; i++){aprob[i]=probCorrect[aqual[i]];}
		for(jint i=0; i<blen; i++){bprob[i]=probCorrect[bqual[i]];}
	}
	{
		jfloat x=findBestRatio_WithQualities(a_bases, a_bases_length, b_bases, b_bases_length, a_quality, b_quality, aprob, bprob, minOverlap0, minOverlap, minInsert, maxRatio, offset);
		if(x>maxRatio){
			rvector[2]=minLength;
			rvector[4]=0;
			return -1;
		}
		maxRatio=min(maxRatio, x);
	}

	const jfloat altBadlimit=max(maxRatio, 0.07f)*2.0f*alen+1;
	const jfloat margin2=(margin+offset)/minLength;
	
	jint bestInsert=-1;
	jfloat bestBad=minLength;
	jfloat bestRatio=1;
	jboolean ambig=0;
	
	const jint largestInsertToTest=(alen+blen-minOverlap0);
	const jint smallestInsertToTest=minInsert0;
	for(jint insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
		jfloat good=0, bad=0;

		const jint istart=(insert<=blen ? 0 : insert-blen);
		const jint jstart=(insert>=blen ? 0 : blen-insert);

		const jint overlapLength=min(alen-istart, min(blen-jstart, insert));
		const jfloat badlimit=min(altBadlimit, min(bestRatio, maxRatio)*margin*overlapLength);

		const jint imax=istart+overlapLength;
		for(jint i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
			const jbyte ca=abases[i], cb=bbases[j];
			const jfloat x=aprob[i]*bprob[j];

			if(ca==cb){good+=x;}
			else{bad+=x;}
		}
		
		if(bad<=badlimit){
			if(bad==0 && good>minOverlap0 && good<minOverlap){
				rvector[2]=(jint)bestBad;
				rvector[4]=1;
				return -1;
			}

			jfloat ratio=(bad+offset)/overlapLength;

			if(ratio<bestRatio*margin){

				ambig=(ratio*margin>=bestRatio || good<minOverlap);
				if(ratio<bestRatio){
					bestInsert=insert;
					bestBad=bad;
					bestRatio=ratio;
				}
				if(ambig && bestRatio<margin2){
					rvector[2]=(jint)bestBad;
					rvector[4]=1;
					return -1;
				}
			}
		}
	}
	
	if(!ambig && bestRatio>maxRatio){bestInsert=-1;}
	
	rvector[2]=(jint)bestBad;
	rvector[4]=(ambig ? 1 : 0);
		
	return (bestInsert<0 ? -1 : bestInsert);
}
//Fixed - BB Oct 27 2015
jint mateByOverlapRatio(jbyte * a_bases, const jint a_bases_length, jbyte * b_bases, const jint b_bases_length, jint * rvector, jint minOverlap0, jint minOverlap, const jint minInsert0, const jint minInsert, jfloat maxRatio, const jfloat margin, const jfloat offset, const jfloat gIncr, const jfloat bIncr) {

	minOverlap=max(4, max(minOverlap0, minOverlap));
	minOverlap0=mid(4, minOverlap0, minOverlap);

        const jbyte *abases=a_bases, *bbases=b_bases;
        const jint alen=a_bases_length, blen=b_bases_length;
	const jint minLength=min(alen, blen);
	{
		jfloat x=findBestRatio(a_bases, a_bases_length, b_bases, b_bases_length, minOverlap0, minOverlap, minInsert, maxRatio, offset, gIncr, bIncr);
		if(x>=maxRatio){
			rvector[2]=minLength;
			rvector[4]=0;
			return -1;
		}
		maxRatio=min(maxRatio, x);
	}
	
	const jfloat altBadlimit=max(maxRatio, 0.07f)*2.0f*alen+1;
	const jfloat margin2=(margin+offset)/minLength;
	const jbyte N='N';
	
	jint bestInsert=-1;
	jfloat bestBad=minLength;
	jfloat bestRatio=1;
	jboolean ambig=0;
	
	const jint largestInsertToTest=(alen+blen-minOverlap0);
	const jint smallestInsertToTest=minInsert0;
	for(jint insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
		const jint istart=(insert<=blen ? 0 : insert-blen);
		const jint jstart=(insert>=blen ? 0 : blen-insert);
		const jint overlapLength=min(alen-istart, min(blen-jstart, insert));
		
		const jfloat badlimit=(min(altBadlimit, min(bestRatio, maxRatio)*margin*overlapLength));
		jfloat good=0, bad=0;
		
		const int imax=istart+overlapLength;
		for(jint i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
			const jbyte ca=abases[i], cb=bbases[j];
			
			if(ca==cb){
				if(ca!=N){good+=gIncr;}
			}else{bad+=bIncr;}
		}
		
		if(bad<=badlimit){
			if(bad==0 && good>minOverlap0 && good<minOverlap){
				rvector[2]=(jint)bestBad;
				rvector[4]=1;
				return -1;
			}

			jfloat ratio=(bad+offset)/overlapLength;
			
			if(ratio<bestRatio*margin){

				ambig=(ratio*margin>=bestRatio || good<minOverlap);
				if(ratio<bestRatio){
					bestInsert=insert;
					bestBad=bad;
					bestRatio=ratio;
				}
				if(ambig && bestRatio<margin2){
					rvector[2]=(int)bestBad;
					rvector[4]=1;
					return -1;
				}
			}
		}
	}
	
	if(!ambig && bestRatio>maxRatio){bestInsert=-1;}
	
	rvector[2]=(jint)bestBad;
	rvector[4]=(ambig ? 1 : 0);
		
	return (bestInsert<0 ? -1 : bestInsert);
}

JNIEXPORT jint JNICALL Java_jgi_BBMergeOverlapper_mateByOverlapJNI_WithQualities(
							JNIEnv *env,
							jobject obj,
							jbyteArray a_bases,
							jbyteArray b_bases,
							jbyteArray a_quality,
							jbyteArray b_quality,
							jfloatArray aprob,
							jfloatArray bprob,
							jintArray rvector,
							jint minOverlap0,
							jint minOverlap,
							jint minInsert0,
							jint minInsert,
							jfloat maxRatio,
							jfloat margin,
							jfloat offset
                                                        ) {
   jbyte * ja_quality = NULL;
   jbyte * jb_quality = NULL;
   jfloat * japrob = NULL;
   jfloat * jbprob = NULL;

   // Get the size of the read and the reference arrays
   const jint a_bases_length = (*env)->GetArrayLength(env, a_bases);
   const jint b_bases_length = (*env)->GetArrayLength(env, b_bases);

   // Copy arrays from Java
   jbyte * ja_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, a_bases, NULL);
   jbyte * jb_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, b_bases, NULL);
   if(a_quality!=NULL) {ja_quality = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, a_quality, NULL);}
   if(b_quality!=NULL) {jb_quality = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, b_quality, NULL);}
   if(aprob!=NULL) {japrob = (jfloat*)(*env)->GetPrimitiveArrayCritical(env, aprob, NULL);}
   if(bprob!=NULL) {jbprob = (jfloat*)(*env)->GetPrimitiveArrayCritical(env, bprob, NULL);}
   jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);

   const jint returnVal = mateByOverlapRatio_WithQualities(ja_bases, a_bases_length, jb_bases, b_bases_length, ja_quality, jb_quality, japrob, jbprob, jrvector, minOverlap0, minOverlap, minInsert0, minInsert, maxRatio, margin, offset);

   // Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
   (*env)->ReleasePrimitiveArrayCritical(env, a_bases, ja_bases, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, b_bases, jb_bases, JNI_ABORT);
   if(ja_quality!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, a_quality, ja_quality, JNI_ABORT);}
   if(jb_quality!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, b_quality, jb_quality, JNI_ABORT);}
   if(japrob!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, aprob, japrob, JNI_ABORT);}
   if(jbprob!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, bprob, jbprob, JNI_ABORT);}
   (*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);

   return returnVal;
}

JNIEXPORT jint JNICALL Java_jgi_BBMergeOverlapper_mateByOverlapJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray a_bases,
							jbyteArray b_bases,
							jbyteArray a_quality,
							jbyteArray b_quality,
							jfloatArray aprob,
							jfloatArray bprob,
							jintArray rvector,
							jint minOverlap0,
							jint minOverlap,
							jint minInsert0,
							jint margin,
							jint maxMismatches0,
							jint maxMismatches,
							jint minq
                                                        ) {
   jbyte * ja_quality = NULL;
   jbyte * jb_quality = NULL;
   jfloat * japrob = NULL;
   jfloat * jbprob = NULL;

   // Get the size of the read and the reference arrays
   const jint a_bases_length = (*env)->GetArrayLength(env, a_bases);
   const jint b_bases_length = (*env)->GetArrayLength(env, b_bases);

   // Copy arrays from Java
   jbyte * ja_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, a_bases, NULL);
   jbyte * jb_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, b_bases, NULL);
   if(a_quality!=NULL) {ja_quality = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, a_quality, NULL);}
   if(b_quality!=NULL) {jb_quality = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, b_quality, NULL);}
   if(aprob!=NULL) {japrob = (jfloat*)(*env)->GetPrimitiveArrayCritical(env, aprob, NULL);}
   if(bprob!=NULL) {jbprob = (jfloat*)(*env)->GetPrimitiveArrayCritical(env, bprob, NULL);}
   jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);

   const jint returnVal = mateByOverlap(ja_bases, a_bases_length, jb_bases, b_bases_length, ja_quality, jb_quality, japrob, jbprob, jrvector, minOverlap0, minOverlap, minInsert0, margin, maxMismatches0, maxMismatches, minq);

   // Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
   (*env)->ReleasePrimitiveArrayCritical(env, a_bases, ja_bases, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, b_bases, jb_bases, JNI_ABORT);
   if(ja_quality!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, a_quality, ja_quality, JNI_ABORT);}
   if(jb_quality!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, b_quality, jb_quality, JNI_ABORT);}
   if(japrob!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, aprob, japrob, JNI_ABORT);}
   if(jbprob!=NULL) {(*env)->ReleasePrimitiveArrayCritical(env, bprob, jbprob, JNI_ABORT);}
   (*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);

   return returnVal;
}

JNIEXPORT jint JNICALL Java_jgi_BBMergeOverlapper_mateByOverlapRatioJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray a_bases,
							jbyteArray b_bases,
							jintArray rvector,
							jint minOverlap0,
							jint minOverlap,
							jint minInsert0,
							jint minInsert,
							jfloat maxRatio,
							jfloat margin,
							jfloat offset,
							jfloat gIncr,
							jfloat bIncr
                                                        ) {
   // Get the size of the read and the reference arrays
   const jint a_bases_length = (*env)->GetArrayLength(env, a_bases);
   const jint b_bases_length = (*env)->GetArrayLength(env, b_bases);

   // Copy arrays from Java
   jbyte * ja_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, a_bases, NULL);
   jbyte * jb_bases = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, b_bases, NULL);
   jint * jrvector = (jint*)(*env)->GetPrimitiveArrayCritical(env, rvector, NULL);

   const jint returnVal = mateByOverlapRatio(ja_bases, a_bases_length, jb_bases, b_bases_length, jrvector, minOverlap0, minOverlap, minInsert0, minInsert, maxRatio, margin, offset, gIncr, bIncr);

   // Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
   (*env)->ReleasePrimitiveArrayCritical(env, a_bases, ja_bases, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, b_bases, jb_bases, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, rvector, jrvector, 0);

   return returnVal;
}

