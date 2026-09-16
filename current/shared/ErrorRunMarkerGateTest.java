package shared;

import java.util.Arrays;
import java.util.Random;

/** Differential gate checks against the separately retained pre-gate source.
 * @author Fischl
 */
public final class ErrorRunMarkerGateTest {
	public static void main(final String[] args){
		final Random random=new Random(20260910);int cases=0,excluded=0;
		final float[] ratios={1,1.1f,2,4,Float.MAX_VALUE,Float.POSITIVE_INFINITY,Float.NaN};
		for(int trial=0;trial<6000;trial++){
			final int k=2+random.nextInt(31),size=1+random.nextInt(96);final int[] counts=new int[size];
			final int mode=trial%6;final int flat=random.nextInt(101);
			for(int i=0;i<size;i++){
				counts[i]=mode==0?flat:mode==1?random.nextInt(5):mode==2?random.nextInt(101):mode==3?(i%13==0?Integer.MAX_VALUE:1):mode==4?random.nextInt(5)-2:20;
			}
			if(mode==5 && size>k+4){final int start=2+random.nextInt(size-k-3);counts[start-2]=counts[start-1]=100;counts[start+k]=counts[start+k+1]=100;}
			final float ratio=ratios[trial%ratios.length];final int predicateMode=trial%3;
			final byte[] a=new byte[size+k-1],qa=trial%2==0?null:new byte[a.length];
			for(int i=0;i<a.length;i++){a[i]=(byte)(i%17==0?'N':"ACGT".charAt(i&3));if(qa!=null){qa[i]=(byte)(i%94);}}
			final byte[] b=a.clone(),qb=qa==null?null:qa.clone(),before=a.clone(),qbefore=qa==null?null:qa.clone();final int[] saved=counts.clone();
			final ErrorRunMarker.Result actual=ErrorRunMarker.mark(a,qa,counts,size,k,ratio,new ErrorRunMarker.ErrorPredicate(){public boolean isError(int h,int l){return predicate(h,l,predicateMode);}});
			final ErrorRunMarkerReference.Result expected=ErrorRunMarkerReference.mark(b,qb,counts,size,k,ratio,new ErrorRunMarkerReference.ErrorPredicate(){public boolean isError(int h,int l){return predicate(h,l,predicateMode);}});
			require(Arrays.equals(a,b) && Arrays.equals(qa,qb) && Arrays.equals(counts,saved),"Gate changed sequence/quality/profile at trial "+trial);
			require(actual.marked==expected.marked && actual.exactKRuns==expected.exactKRuns && actual.exactKMinus1Runs==expected.exactKMinus1Runs && actual.edgeRunsMarked==expected.edgeRunsMarked && actual.skippedEndRuns==expected.skippedEndRuns && actual.skippedMergedRuns==expected.skippedMergedRuns,"Gate changed counters at trial "+trial);
			if(!ErrorRunMarker.hasPotentialRun(counts,size,Tools.max(1f,ratio))){
				excluded++;require(expected.marked==0 && expected.exactKRuns==0 && expected.exactKMinus1Runs==0 && expected.edgeRunsMarked==0 && expected.skippedEndRuns==0 && expected.skippedMergedRuns==0 && Arrays.equals(b,before) && Arrays.equals(qb,qbefore),"Excluded profile has an observable original effect");
			}cases++;
		}
		require(excluded>0 && excluded<cases,"Differential panel must exercise both gate paths");
		System.out.println("MARKER_GATE_TEST_OK cases="+cases+" excluded="+excluded+"; all counters and complete base/quality arrays equal frozen original");
	}
	private static boolean predicate(int high,int low,int mode){return mode==0?high>4L*low:mode==1?true:((high^low)&1)==0;}
	private static void require(final boolean ok,final String why){if(!ok){throw new AssertionError(why);}}
}
