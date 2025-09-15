#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include "align2_BandedAlignerJNI.h"

// C doesn't have min() or max() so we define our own
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; }) 

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; }) 

// I am not aware of any place where BBMap changes this value
// taken from the value in BandedAligner.java, so I'm using a define
// to true so the compiler can optimize the if statements out
#define penalizeOffCenter 1

// Need these prototypes since these functions can call each other
jint alignForward(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber
	);

jint alignForwardRC(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber,
	jbyte * baseToComplementExtended
	);

jint alignReverse(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber
	);

jint alignReverseRC(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber,
	jbyte * baseToComplementExtended
	);

jint lastOffsetFunc(jint * array, jint halfWidth){
	const jint center=halfWidth+1;
	jint minLoc=center;
	for(jint i=1; i<=halfWidth; i++){
		if(array[center+i]<array[minLoc]){
			minLoc=center+i;
		}
		if(array[center-i]<array[minLoc]){
			minLoc=center-i;
		}
	}
	return center-minLoc;
}

jint penalizeOffCenterFunc(jint * array, jint halfWidth, jint big){
	const jint center=halfWidth+1;
	jint edits=array[center];
	for(jint i=1; i<=halfWidth; i++){
		array[center+i]=min(big, array[center+i]+i);
		edits=min(edits, array[center+i]);
		array[center-i]=min(big, array[center-i]+i);
		edits=min(edits, array[center-i]);
	}
	return edits;
}

jint alignForward(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber
	) {

	if(query_length-qstart>ref_length-rstart){
		jint x=alignForward(ref, query, ref_length, query_length, rstart,
			qstart, maxEdits, exact, lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, baseToNumber);
		const jint temp = *lastQueryLoc;
		*lastQueryLoc = *lastRefLoc;
		*lastRefLoc = temp;
		return x;
	}

	const jint big=999;
	jint edits=0;
	jint row=0;
	*lastRow=-1;
	*lastEdits=0;
	*lastOffset=0;

	const jint width=min(maxWidth, (maxEdits*2)+1);
	const jint halfWidth=width/2;
	const jboolean inexact=!exact;

	jint qloc=qstart;
	jint rsloc=rstart-halfWidth;
	const jint xlines=query_length-qstart;
	const jint ylines=ref_length-rstart;
	const jint len=min(xlines, ylines);

	if(len<1){
		return 0;
	}

	jint arrayCurrentArray[maxWidth+2];
	jint arrayPrevArray[maxWidth+2];
	jint * arrayCurrent=arrayCurrentArray;
	jint * arrayPrev=arrayPrevArray;
	jint * arrayTemp;

	for(jint i=0; i<maxWidth+2; i++) {
		arrayCurrent[i]=big;
		arrayPrev[i]=big;
	}

	{
		const jbyte q=query[qloc];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		edits=big;
		jint mloc=1+(colStart-rsloc);
		for(jint col=colStart; col<colLimit; mloc++, col++){
			const jbyte r=ref[col];
			const jint score=(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		row++; qloc++; rsloc++;
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}
	
	for(row=1; row<len; row++, qloc++, rsloc++){
		arrayTemp=arrayCurrent;
		arrayCurrent=arrayPrev;
		arrayPrev=arrayTemp;
		const jbyte q=query[qloc];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		for(jint i=0; i<maxWidth+2; i++) {
			arrayCurrent[i]=big;
		}
		edits=big;
		jint mloc=1+(colStart-rsloc);
		const jint forceDiag=(row==len-1);
		for(jint col=colStart; col<colLimit; mloc++, col++){
			const jbyte r=ref[col];
			const jint scoreUp=arrayPrev[mloc+1]+1;
			const jint scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			const jint scoreLeft=arrayCurrent[mloc-1]+1;
			const jint score=(forceDiag || col==ref_length-1) ? scoreDiag : min(scoreUp, min(scoreDiag, scoreLeft));
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		if(edits>maxEdits){
			row++;
			break;
		}
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}

	*lastRow=row-1;
	*lastEdits=edits;
	*lastQueryLoc=qloc-1;
	*lastOffset=lastOffsetFunc(arrayCurrent, halfWidth);
	*lastRefLoc=rsloc+halfWidth-(*lastOffset)-1;
	while((*lastRefLoc)>=ref_length || (*lastQueryLoc)>=query_length){(*lastRefLoc)--;(*lastQueryLoc)--;}

	return edits;
}

jint alignForwardRC(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber,
	jbyte * baseToComplementExtended
	) {

	if(qstart+1>ref_length-rstart){
		jint x=alignReverseRC(ref, query, ref_length, query_length, rstart,
			qstart, maxEdits, exact, lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, baseToNumber, baseToComplementExtended);
		const jint temp = *lastQueryLoc;
		*lastQueryLoc = *lastRefLoc;
		*lastRefLoc = temp;
		return x;
	}

	const jint big=999;
	jint edits=0;
	jint row=0;
	*lastRow=-1;
	*lastEdits=0;
	*lastOffset=0;

	const jint width=min(maxWidth, (maxEdits*2)+1);
	const jint halfWidth=width/2;
	const jboolean inexact=!exact;

	jint qloc=qstart;
	jint rsloc=rstart-halfWidth;
	const jint xlines=qstart+1;
	const jint ylines=ref_length-rstart;
	const jint len=min(xlines, ylines);

	if(len<1){
		return 0;
	}

	jint arrayCurrentArray[maxWidth+2];
	jint arrayPrevArray[maxWidth+2];
	jint * arrayCurrent=arrayCurrentArray;
	jint * arrayPrev=arrayPrevArray;
	jint * arrayTemp;

	for(jint i=0; i<maxWidth+2; i++) {
		arrayCurrent[i]=big;
		arrayPrev[i]=big;
	}

	{
		const jbyte q=baseToComplementExtended[query[qloc]];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		edits=big;
		jint mloc=1+(colStart-rsloc);
		for(jint col=colStart; col<colLimit; mloc++, col++){
			const jbyte r=ref[col];
			const jint score=(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		row++; qloc--; rsloc++;
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}
	
	for(row=1; row<len; row++, qloc--, rsloc++){
		arrayTemp=arrayCurrent;
		arrayCurrent=arrayPrev;
		arrayPrev=arrayTemp;
		const jbyte q=baseToComplementExtended[query[qloc]];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		for(jint i=0; i<maxWidth+2; i++) {
			arrayCurrent[i]=big;
		}
		edits=big;
		jint mloc=1+(colStart-rsloc);
		const jint forceDiag=(row==len-1);
		for(jint col=colStart; col<colLimit; mloc++, col++){
			const jbyte r=ref[col];
			const jint scoreUp=arrayPrev[mloc+1]+1;
			const jint scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			const jint scoreLeft=arrayCurrent[mloc-1]+1;
			const jint score=(forceDiag || col==ref_length-1) ? scoreDiag : min(scoreUp, min(scoreDiag, scoreLeft));
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		if(edits>maxEdits){row++; break;}
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}

	*lastRow=row-1;
	*lastEdits=edits;
	*lastOffset=lastOffsetFunc(arrayCurrent, halfWidth);
	*lastQueryLoc=qloc+1;
	*lastRefLoc=rsloc+halfWidth-(*lastOffset)-1;
	while((*lastRefLoc)>=ref_length || (*lastQueryLoc)<0){(*lastRefLoc)--; (*lastQueryLoc)++;}

	return edits;
}

jint alignReverse(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber
	) {

	if(qstart>rstart){
		jint x=alignReverse(ref, query, ref_length, query_length, rstart,
			qstart, maxEdits, exact, lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, baseToNumber);
		const jint temp = *lastQueryLoc;
		*lastQueryLoc = *lastRefLoc;
		*lastRefLoc = temp;
		return x;
	}

	const jint big=999;
	jint edits=0;
	jint row=0;
	*lastRow=-1;
	*lastEdits=0;
	*lastOffset=0;

	const jint width=min(maxWidth, (maxEdits*2)+1);
	const jint halfWidth=width/2;
	const jboolean inexact=!exact;

	jint qloc=qstart;
	jint rsloc=rstart-halfWidth;
	const jint xlines=qstart+1;
	const jint ylines=rstart+1;
	const jint len=min(xlines, ylines);

	if(len<1){
		return 0;
	}

	jint arrayCurrentArray[maxWidth+2];
	jint arrayPrevArray[maxWidth+2];
	jint * arrayCurrent=arrayCurrentArray;
	jint * arrayPrev=arrayPrevArray;
	jint * arrayTemp;

	for(jint i=0; i<maxWidth+2; i++) {
		arrayCurrent[i]=big;
		arrayPrev[i]=big;
	}

	{
		const jbyte q=query[qloc];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		edits=big;
		jint mloc=1+width-(colLimit-rsloc);
		for(jint col=colLimit-1; col>=colStart; mloc++, col--){
			const jbyte r=ref[col];
			const jint score=(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		row++; qloc--; rsloc--;
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}
	
	for(row=1; row<len; row++, qloc--, rsloc--){
		arrayTemp=arrayCurrent;
		arrayCurrent=arrayPrev;
		arrayPrev=arrayTemp;
		const jbyte q=query[qloc];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		for(jint i=0; i<maxWidth+2; i++) {
			arrayCurrent[i]=big;
		}
		edits=big;
		jint mloc=1+width-(colLimit-rsloc);
		const jint forceDiag=(row==len-1);
		for(jint col=colLimit-1; col>=colStart; mloc++, col--){
			const jbyte r=ref[col];
			const jint scoreUp=arrayPrev[mloc+1]+1;
			const jint scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			const jint scoreLeft=arrayCurrent[mloc-1]+1;
			const jint score=(forceDiag || col==0) ? scoreDiag : min(scoreUp, min(scoreDiag, scoreLeft));
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		if(edits>maxEdits){row++; break;}
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}

	*lastRow=row-1;
	*lastEdits=edits;
	*lastOffset=lastOffsetFunc(arrayCurrent, halfWidth);
	*lastQueryLoc=qloc+1;
	*lastRefLoc=rsloc+halfWidth+(*lastOffset)+1;
	while((*lastRefLoc)<0 || (*lastQueryLoc)<0){(*lastRefLoc)++; (*lastQueryLoc)++;}

	return edits;
}

jint alignReverseRC(
	jbyte * query,
	jbyte * ref,
	jint query_length,
	jint ref_length,
	jint qstart,
	jint rstart,
	jint maxEdits,
	jboolean exact,
	jint * lastQueryLoc,
	jint * lastRefLoc,
	jint * lastRow,
	jint * lastEdits,
	jint * lastOffset,
	jint maxWidth,
	jbyte * baseToNumber,
	jbyte * baseToComplementExtended
	) {

	if(query_length-qstart>rstart+1){
		jint x=alignForwardRC(ref, query, ref_length, query_length, rstart,
			qstart, maxEdits, exact, lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, baseToNumber, baseToComplementExtended);
		const jint temp = *lastQueryLoc;
		*lastQueryLoc = *lastRefLoc;
		*lastRefLoc = temp;
		return x;
	}

	const jint big=999;
	jint edits=0;
	jint row=0;
	*lastRow=-1;
	*lastEdits=0;
	*lastOffset=0;

	const jint width=min(maxWidth, (maxEdits*2)+1);
	const jint halfWidth=width/2;
	const jboolean inexact=!exact;

	jint qloc=qstart;
	jint rsloc=rstart-halfWidth;
	const jint xlines=query_length-qstart;
	const jint ylines=rstart+1;
	const jint len=min(xlines, ylines);

	if(len<1){
		return 0;
	}

	jint arrayCurrentArray[maxWidth+2];
	jint arrayPrevArray[maxWidth+2];
	jint * arrayCurrent=arrayCurrentArray;
	jint * arrayPrev=arrayPrevArray;
	jint * arrayTemp;

	for(jint i=0; i<maxWidth+2; i++) {
		arrayCurrent[i]=big;
		arrayPrev[i]=big;
	}

	{
		const jbyte q=baseToComplementExtended[query[qloc]];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		edits=big;
		jint mloc=1+width-(colLimit-rsloc);
		for(jint col=colLimit-1; col>=colStart; mloc++, col--){
			const jbyte r=ref[col];
			const jint score=(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		row++; qloc++; rsloc--;
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}
	
	for(row=1; row<len; row++, qloc++, rsloc--){
		arrayTemp=arrayCurrent;
		arrayCurrent=arrayPrev;
		arrayPrev=arrayTemp;
		const jbyte q=baseToComplementExtended[query[qloc]];
		const jint colStart=max(0, rsloc);
		const jint colLimit=min(rsloc+width, ref_length);
		for(jint i=0; i<maxWidth+2; i++) {
			arrayCurrent[i]=big;
		}
		edits=big;
		jint mloc=1+width-(colLimit-rsloc);
		const jint forceDiag=(row==len-1);
		for(jint col=colLimit-1; col>=colStart; mloc++, col--){
			const jbyte r=ref[col];
			const jint scoreUp=arrayPrev[mloc+1]+1;
			const jint scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!(baseToNumber[q]>=0) || !(baseToNumber[r]>=0))) ? 0 : 1);
			const jint scoreLeft=arrayCurrent[mloc-1]+1;
			const jint score=(forceDiag || col==0) ? scoreDiag : min(scoreUp, min(scoreDiag, scoreLeft));
			arrayCurrent[mloc]=score;
			edits=min(edits, score);
		}
		if(edits>maxEdits){row++; break;}
	}
	if(penalizeOffCenter){
		edits=penalizeOffCenterFunc(arrayCurrent, halfWidth, big);
	}

	*lastRow=row-1;
	*lastEdits=edits;
	*lastOffset=lastOffsetFunc(arrayCurrent, halfWidth);
	*lastQueryLoc=qloc-1;
	*lastRefLoc=rsloc+halfWidth+(*lastOffset)+1;
	while((*lastRefLoc)<0 || (*lastQueryLoc)>=query_length){(*lastRefLoc)++; (*lastQueryLoc)--;}
	return edits;
}

// The other three JNICALL functions are almost identical
JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignForwardJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray query,
							jbyteArray ref,
							jint qstart,
							jint rstart,
							jint maxEdits,
							jboolean exact,
                                                        jint maxWidth,
                                                        jbyteArray baseToNumber,
							jintArray returnVals
                                                        ) {
   jint edits = 0;

   // Get the size of the read and the reference arrays
   jint ref_length = (*env)->GetArrayLength(env, ref);
   jint query_length = (*env)->GetArrayLength(env, query);

   // Copy arrays from Java
   jbyte * jbaseToNumber = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, baseToNumber, NULL);
   jbyte * jref = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
   jbyte * jquery = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
   jint * jreturnVals = (jint*)(*env)->GetPrimitiveArrayCritical(env, returnVals, NULL);

   // Using pointers for variables that need to be passed back to Java so the called functions can update them
   jint * lastQueryLoc = &jreturnVals[0];
   jint * lastRefLoc = &jreturnVals[1];
   jint * lastRow = &jreturnVals[2];
   jint * lastEdits = &jreturnVals[3];
   jint * lastOffset = &jreturnVals[4];

   // Call the fillLimitedX function in C; the 5 return values will be in jresult[]
   edits = alignForward(jquery, jref, query_length, ref_length, qstart, rstart, maxEdits, exact,
      lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, jbaseToNumber);

   // Release Java arrays; 0 copies the array back to Java, JNI_ABORT does not copy the current array values to Java
   (*env)->ReleasePrimitiveArrayCritical(env, baseToNumber, jbaseToNumber, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, returnVals, jreturnVals, 0);

   return edits;
}

JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignForwardRCJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray query,
							jbyteArray ref,
							jint qstart,
							jint rstart,
							jint maxEdits,
							jboolean exact,
                                                        jint maxWidth,
                                                        jbyteArray baseToNumber,
                                                        jbyteArray baseToComplementExtended,
							jintArray returnVals
                                                        ) {
   jint edits = 0;

   jint ref_length = (*env)->GetArrayLength(env, ref);
   jint query_length = (*env)->GetArrayLength(env, query);

   jbyte * jbaseToComplementExtended = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, baseToComplementExtended, NULL);
   jbyte * jbaseToNumber = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, baseToNumber, NULL);
   jbyte * jref = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
   jbyte * jquery = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
   jint * jreturnVals = (jint*)(*env)->GetPrimitiveArrayCritical(env, returnVals, NULL);

   jint * lastQueryLoc = &jreturnVals[0];
   jint * lastRefLoc = &jreturnVals[1];
   jint * lastRow = &jreturnVals[2];
   jint * lastEdits = &jreturnVals[3];
   jint * lastOffset = &jreturnVals[4];

   edits = alignForwardRC(jquery, jref, query_length, ref_length, qstart, rstart, maxEdits, exact,
      lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, jbaseToNumber, jbaseToComplementExtended);

   (*env)->ReleasePrimitiveArrayCritical(env, baseToComplementExtended, jbaseToComplementExtended, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, baseToNumber, jbaseToNumber, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, returnVals, jreturnVals, 0);

   return edits;
}

JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignReverseJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray query,
							jbyteArray ref,
							jint qstart,
							jint rstart,
							jint maxEdits,
							jboolean exact,
                                                        jint maxWidth,
                                                        jbyteArray baseToNumber,
							jintArray returnVals
                                                        ) {
   jint edits = 0;

   jint ref_length = (*env)->GetArrayLength(env, ref);
   jint query_length = (*env)->GetArrayLength(env, query);

   jbyte * jbaseToNumber = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, baseToNumber, NULL);
   jbyte * jref = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
   jbyte * jquery = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
   jint * jreturnVals = (jint*)(*env)->GetPrimitiveArrayCritical(env, returnVals, NULL);

   jint * lastQueryLoc = &jreturnVals[0];
   jint * lastRefLoc = &jreturnVals[1];
   jint * lastRow = &jreturnVals[2];
   jint * lastEdits = &jreturnVals[3];
   jint * lastOffset = &jreturnVals[4];

   edits = alignReverse(jquery, jref, query_length, ref_length, qstart, rstart, maxEdits, exact,
      lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, jbaseToNumber);

   (*env)->ReleasePrimitiveArrayCritical(env, baseToNumber, jbaseToNumber, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, returnVals, jreturnVals, 0);

   return edits;
}

JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignReverseRCJNI(
							JNIEnv *env,
							jobject obj,
							jbyteArray query,
							jbyteArray ref,
							jint qstart,
							jint rstart,
							jint maxEdits,
							jboolean exact,
                                                        jint maxWidth,
                                                        jbyteArray baseToNumber,
                                                        jbyteArray baseToComplementExtended,
							jintArray returnVals
                                                        ) {
   jint edits = 0;

   jint ref_length = (*env)->GetArrayLength(env, ref);
   jint query_length = (*env)->GetArrayLength(env, query);

   jbyte * jbaseToComplementExtended = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, baseToComplementExtended, NULL);
   jbyte * jbaseToNumber = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, baseToNumber, NULL);
   jbyte * jref = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, ref, NULL);
   jbyte * jquery = (jbyte*)(*env)->GetPrimitiveArrayCritical(env, query, NULL);
   jint * jreturnVals = (jint*)(*env)->GetPrimitiveArrayCritical(env, returnVals, NULL);

   jint * lastQueryLoc = &jreturnVals[0];
   jint * lastRefLoc = &jreturnVals[1];
   jint * lastRow = &jreturnVals[2];
   jint * lastEdits = &jreturnVals[3];
   jint * lastOffset = &jreturnVals[4];

   edits = alignReverseRC(jquery, jref, query_length, ref_length, qstart, rstart, maxEdits, exact,
      lastQueryLoc, lastRefLoc, lastRow, lastEdits, lastOffset, maxWidth, jbaseToNumber, jbaseToComplementExtended);

   (*env)->ReleasePrimitiveArrayCritical(env, baseToComplementExtended, jbaseToComplementExtended, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, baseToNumber, jbaseToNumber, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, ref, jref, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, query, jquery, JNI_ABORT);
   (*env)->ReleasePrimitiveArrayCritical(env, returnVals, jreturnVals, 0);

   return edits;
}

