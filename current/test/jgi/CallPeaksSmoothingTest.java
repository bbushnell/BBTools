package test.jgi;

import java.util.Arrays;

import jgi.CallPeaks;

/**
 * Finite constant-histogram checks for CallPeaks smoothing arithmetic.
 * Every triangular tap clamps to a constant value, so the correct result is
 * that same value. Baseline expectations independently describe the known int
 * intermediate overflows; fixed expectations require constant preservation.
 * Inputs contain at most eight elements. No peak calling or large data is used.
 * @author Shinobu
 */
public final class CallPeaksSmoothingTest{

	/** Accepts only phase=baseline or phase=fixed; prints every checked result. */
	public static void main(final String[] args){
		if(args.length!=1 || (!args[0].equals("phase=baseline") && !args[0].equals("phase=fixed"))){
			throw new IllegalArgumentException("Set phase=baseline or phase=fixed.");
		}
		final boolean baseline=args[0].equals("phase=baseline");
		int mismatches=0;
		for(final boolean progressive : new boolean[] {false, true}){
			final String mode=progressive ? "progressive" : "fixed";
			mismatches+=checkInt("int_product_"+mode, 1, 1073741824, 2, progressive, baseline, 0, 0);
			mismatches+=checkInt("int_radius_"+mode, 1, 1, 46341, progressive, baseline, -1, 0);
			mismatches+=checkLong("long_radius_"+mode, 1, 100, 46341, progressive, baseline, -100, 0);
			mismatches+=checkInt("int_control_"+mode, 8, 17, 2, progressive, baseline, 17, 0);
			mismatches+=checkLong("long_control_"+mode, 8, 100, 2, progressive, baseline, 100, 0);
		}
		//The first radius increment occurs after index 6; only index 7 uses 46341.
		mismatches+=checkInt("int_radius_increment", 8, 1, 46340, true, baseline, -1, 7);
		mismatches+=checkLong("long_radius_increment", 8, 100, 46340, true, baseline, -100, 7);
		final int expected=baseline ? 8 : 0;
		if(mismatches!=expected){throw new AssertionError("Constant-preservation failures: expected "+expected+", got "+mismatches);}
		System.out.println("PASS "+args[0]+" cases=12 mismatches="+mismatches);
	}

	/** Checks an int overload, the complete result, and input preservation. */
	private static int checkInt(final String name, final int length, final int value, final int radius,
			final boolean progressive, final boolean baseline, final int buggyValue, final int firstBuggyIndex){
		assert(length>0 && radius>0 && firstBuggyIndex>=0 && firstBuggyIndex<length) :
			"Constant-histogram fixture needs positive dimensions and a valid expected-change index: "+name;
		final int[] input=new int[length];
		Arrays.fill(input, value);
		final int[] before=input.clone();
		final int[] result=progressive ? CallPeaks.smoothProgressive(input, radius) : CallPeaks.smooth(input, radius);
		if(result.length!=length || !Arrays.equals(input, before)){
			throw new AssertionError(name+": smoothing changed the input or result length.");
		}
		boolean mismatch=false;
		for(int i=0; i<length; i++){
			final int expected=baseline && i>=firstBuggyIndex ? buggyValue : value;
			if(result[i]!=expected){throw new AssertionError(name+"["+i+"]: expected "+expected+", got "+result[i]);}
			mismatch|=result[i]!=value;
		}
		System.out.println(name+"\t"+Arrays.toString(result));
		return mismatch ? 1 : 0;
	}

	/** Checks a long overload using the same independent constant oracle. */
	private static int checkLong(final String name, final int length, final long value, final int radius,
			final boolean progressive, final boolean baseline, final long buggyValue, final int firstBuggyIndex){
		assert(length>0 && radius>0 && firstBuggyIndex>=0 && firstBuggyIndex<length) :
			"Constant-histogram fixture needs positive dimensions and a valid expected-change index: "+name;
		final long[] input=new long[length];
		Arrays.fill(input, value);
		final long[] before=input.clone();
		final long[] result=progressive ? CallPeaks.smoothProgressive(input, radius) : CallPeaks.smooth(input, radius);
		if(result.length!=length || !Arrays.equals(input, before)){
			throw new AssertionError(name+": smoothing changed the input or result length.");
		}
		boolean mismatch=false;
		for(int i=0; i<length; i++){
			final long expected=baseline && i>=firstBuggyIndex ? buggyValue : value;
			if(result[i]!=expected){throw new AssertionError(name+"["+i+"]: expected "+expected+", got "+result[i]);}
			mismatch|=result[i]!=value;
		}
		System.out.println(name+"\t"+Arrays.toString(result));
		return mismatch ? 1 : 0;
	}
}
