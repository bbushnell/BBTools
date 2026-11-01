package test.shared;

import shared.ErrorRunMarker;

/** Behavioral checks for bounded low-depth run marking and future gate work. */
public class ErrorRunMarkerTest{

	/** Runs finite profiles and sampled-gate counterexamples; arguments are unused. */
	public static void main(String[] args){
		testExactK();
		testExactKMinusOne();
		testEdgeRuns();
		testMergedRun();
		testFlatProfile();
		testNarrowFlankProfiles();
		testSamplePhases();
		System.out.println("ErrorRunMarkerTest PASS: "+checks+" checks");
	}

	/** An isolated k-long plateau changes exactly one base and its quality. */
	private static void testExactK(){
		final int k=5, countSize=20, start=5;
		final int[] counts=profile(countSize, 100);
		setRange(counts, start, start+k, 1);
		final byte[] bases=bases(countSize+k-1), before=bases.clone();
		final byte[] quals=qualities(bases.length), qbefore=quals.clone();
		final ErrorRunMarker.Result result=mark(bases, quals, counts, k);
		check(result.marked==1 && result.exactKRuns==1, "exact-K run was not marked once");
		check(bases[start+k-1]=='N' && quals[start+k-1]==0, "exact-K position/quality incorrect");
		checkUnchangedExcept(before, bases, start+k-1);
		checkQualityUnchangedExcept(qbefore, quals, start+k-1);
	}

	/** An isolated k-1 plateau uses the deletion-signature position with null qualities. */
	private static void testExactKMinusOne(){
		final int k=5, countSize=20, start=5;
		final int[] counts=profile(countSize, 100);
		setRange(counts, start, start+k-1, 1);
		final byte[] bases=bases(countSize+k-1), before=bases.clone();
		final ErrorRunMarker.Result result=mark(bases, null, counts, k);
		check(result.marked==1 && result.exactKMinus1Runs==1, "exact-K-minus-one run was not marked once");
		check(bases[start+k-2]=='N', "exact-K-minus-one position incorrect");
		checkUnchangedExcept(before, bases, start+k-2);
	}

	/** Left and right one-sided plateaus mark their literal edge positions. */
	private static void testEdgeRuns(){
		final int k=5, countSize=12;
		final int[] left=profile(countSize, 100);
		setRange(left, 0, 3, 1);
		final byte[] leftBases=bases(countSize+k-1), leftBefore=leftBases.clone();
		final ErrorRunMarker.Result leftResult=mark(leftBases, null, left, k);
		check(leftResult.marked==1 && leftResult.edgeRunsMarked==1, "left edge run was not marked");
		check(leftBases[2]=='N', "left edge marker position incorrect");
		checkUnchangedExcept(leftBefore, leftBases, 2);

		final int[] right=profile(countSize, 100);
		setRange(right, countSize-3, countSize, 1);
		final byte[] rightBases=bases(countSize+k-1), rightBefore=rightBases.clone();
		final ErrorRunMarker.Result rightResult=mark(rightBases, null, right, k);
		check(rightResult.marked==1 && rightResult.edgeRunsMarked==1, "right edge run was not marked");
		check(rightBases[countSize-3+k-1]=='N', "right edge marker position incorrect");
		checkUnchangedExcept(rightBefore, rightBases, countSize-3+k-1);
	}

	/** A k+1 plateau is reported as merged without changing bases. */
	private static void testMergedRun(){
		final int k=5, countSize=25, start=8;
		final int[] counts=profile(countSize, 100);
		setRange(counts, start, start+k+1, 1);//Longer than K: report, do not mark.
		final byte[] bases=bases(countSize+k-1), before=bases.clone();
		final ErrorRunMarker.Result result=mark(bases, null, counts, k);
		check(result.marked==0, "merged low run was marked");
		check(result.skippedMergedRuns==1, "merged low run was not reported once");
		checkUnchangedExcept(before, bases, -1);
	}

	/** Uniform depth produces no changes to bases or qualities. */
	private static void testFlatProfile(){
		final int k=31, countSize=100;
		final int[] counts=profile(countSize, 20);//Above a normal minCountCorrect threshold.
		final byte[] bases=bases(countSize+k-1), before=bases.clone();
		final byte[] quals=qualities(bases.length), qbefore=quals.clone();
		final ErrorRunMarker.Result result=mark(bases, quals, counts, k);
		check(result.marked==0 && result.exactKRuns==0 && result.exactKMinus1Runs==0,
				"flat profile was marked");
		checkUnchangedExcept(before, bases, -1);
		checkQualityUnchangedExcept(qbefore, quals, -1);
	}

	/**
	 * Keep the low plateau above minCountCorrect and make each high flank only
	 * the two k-mers required by the marker.  Some phase shifts hide both short
	 * flanks from a nine-position sampled gate; those profiles must still mark
	 * correctly, so a future gate cannot be enabled without proving equivalence.
	 */
	private static void testNarrowFlankProfiles(){
		final int k=31, countSize=160;
		int sampledGateMisses=0;
		for(int phase=0; phase<9; phase++){
			final int start=40+phase;
			final int[] counts=profile(countSize, 20);
			setRange(counts, start-2, start, 100);
			setRange(counts, start, start+k, 20);
			setRange(counts, start+k, start+k+2, 100);
			final byte[] bases=bases(countSize+k-1), before=bases.clone();
			final ErrorRunMarker.Result result=mark(bases, null, counts, k);
			check(result.marked==1 && result.exactKRuns==1, "narrow-flank run failed at phase "+phase);
			check(bases[start+k-1]=='N', "narrow-flank position failed at phase "+phase);
			checkUnchangedExcept(before, bases, start+k-1);
			if(!sampledTransition(counts, k)){sampledGateMisses++;}
		}
		check(sampledGateMisses>0, "profiles did not produce a prospective sampled-gate miss");
	}

	/**
	 * Generate exact-K runs at every phase relative to a prospective sampled
	 * profile gate.  The marker must remain phase-independent; this test does
	 * not assert that a future sampled gate is safe, it preserves counterexamples
	 * for that separate optimization decision.
	 */
	private static void testSamplePhases(){
		final int k=31, countSize=100;
		for(int phase=0; phase<9; phase++){
			final int start=10+phase;
			final int[] counts=profile(countSize, 100);
			setRange(counts, start, start+k, 1);
			final byte[] bases=bases(countSize+k-1), before=bases.clone();
			final ErrorRunMarker.Result result=mark(bases, null, counts, k);
			check(result.marked==1 && result.exactKRuns==1, "sample-phase exact-K run failed at phase "+phase);
			check(bases[start+k-1]=='N', "sample-phase marker position failed at phase "+phase);
			checkUnchangedExcept(before, bases, start+k-1);
		}
	}

	/** Runs the marker with the fixed fourfold-depth predicate used by these profiles. */
	private static ErrorRunMarker.Result mark(final byte[] bases, final byte[] quals, final int[] counts, final int k){
		return ErrorRunMarker.mark(bases, quals, counts, counts.length, k, 4f,
			new ErrorRunMarker.ErrorPredicate(){
				@Override
				public boolean isError(final int high, final int low){return high>=low*4;}
			});
	}

	/** Creates an independent constant-depth profile. */
	private static int[] profile(final int size, final int value){
		final int[] out=new int[size];
		for(int i=0; i<size; i++){out[i]=value;}
		return out;
	}

	/** Overwrites the half-open interval [from,to) in the fixture profile. */
	private static void setRange(final int[] values, final int from, final int to, final int value){
		for(int i=from; i<to; i++){values[i]=value;}
	}

	/** Models an intentionally insufficient sampled gate to preserve its counterexamples. */
	private static boolean sampledTransition(final int[] counts, final int k){
		int prev=-1;
		final int incr=Math.min(9, Math.max(1, k/2));
		for(int i=0; i<counts.length; i+=incr){
			final int count=counts[i];
			if(count<3 || (prev>=0 && count>=prev*4)){return true;}
			if(prev>=0 && prev>=count*4){return true;}
			prev=count;
		}
		final int count=counts[counts.length-1];
		return count<3 || (prev>=0 && (count>=prev*4 || prev>=count*4));
	}

	/** Creates all-A bases so any N is attributable to marking. */
	private static byte[] bases(final int length){
		final byte[] out=new byte[length];
		for(int i=0; i<length; i++){out[i]='A';}
		return out;
	}

	/** Creates Phred-30 qualities so zeroed positions are observable. */
	private static byte[] qualities(final int length){
		final byte[] out=new byte[length];
		for(int i=0; i<length; i++){out[i]=30;}
		return out;
	}

	/** Checks every base; changed=-1 requires the entire array to remain unchanged. */
	private static void checkUnchangedExcept(final byte[] before, final byte[] after, final int changed){
		for(int i=0; i<before.length; i++){
			check(i==changed ? after[i]=='N' : before[i]==after[i], "unexpected base change at "+i);
		}
	}

	/** Checks every quality; changed=-1 requires the entire array to remain unchanged. */
	private static void checkQualityUnchangedExcept(final byte[] before, final byte[] after, final int changed){
		for(int i=0; i<before.length; i++){
			check(i==changed ? after[i]==0 : before[i]==after[i], "unexpected quality change at "+i);
		}
	}

	/** Counts and enforces an expectation even when JVM assertions are disabled. */
	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new RuntimeException(message);}
	}

	/** Number of expectations attempted by this process. */
	private static long checks;
}
