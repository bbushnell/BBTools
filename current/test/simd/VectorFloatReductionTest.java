package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks floating reductions against literal scalar-contract expectations.
 * Large finite float sums must remain finite when their double result fits.
 * Maximum selection retains the first value on ties and ignores later NaNs,
 * matching the facade's scalar comparison loop. Run in an isolated JVM because
 * this test temporarily changes the process-wide SIMD flag.
 * @author Shinobu */
public final class VectorFloatReductionTest{

	/** Runs exact finite, NaN and signed-zero controls across dispatch boundaries.
	 * @param args Optional simd=t or simd=f selector */
	public static void main(final String[] args){
		boolean simd=false;
		for(final String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD checks require 256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			failures+=checkSum("null", null, 0); cases++;
			failures+=checkSum("empty", new float[0], 0); cases++;
			failures+=checkMax("null", null, 0); cases++;
			for(final int length : LENGTHS){
				final float[] values=new float[length];
				Arrays.fill(values, 1);
				failures+=checkSum("ones-"+length, values, length); cases++;
				failures+=checkMax("ones-"+length, values, 1); cases++;
				Arrays.fill(values, Float.MAX_VALUE);
				failures+=checkSum("large-finite-"+length, values,
						length*(double)Float.MAX_VALUE); cases++;
				Arrays.fill(values, -Float.MAX_VALUE);
				failures+=checkSum("large-negative-"+length, values,
						-length*(double)Float.MAX_VALUE); cases++;
				for(int i=0; i<length; i++){values[i]=i-(length>>>1);}
				final double rampSum=length*(length-1)*0.5-length*(length>>>1);
				failures+=checkSum("signed-ramp-"+length, values, rampSum); cases++;
				for(int i=0; i<length; i++){
					values[i]=(i&1)==0 ? Float.MAX_VALUE : -Float.MAX_VALUE;
				}
				failures+=checkSum("large-cancellation-"+length, values,
						(length&1)==0 ? 0 : (double)Float.MAX_VALUE); cases++;
				Arrays.fill(values, -1);
				values[0]=Float.NaN;
				failures+=checkMax("first-nan-"+length, values, Float.NaN); cases++;
				Arrays.fill(values, -1);
				values[0]=-0.0f;
				if(length>1){values[length-1]=+0.0f;}
				failures+=checkMax("first-negative-zero-"+length, values, -0.0f); cases++;
				values[0]=+0.0f;
				if(length>1){values[length-1]=-0.0f;}
				failures+=checkMax("first-positive-zero-"+length, values, +0.0f); cases++;
				if(length>1){
					Arrays.fill(values, -1);
					values[0]=1;
					values[length-1]=Float.NaN;
					failures+=checkMax("later-nan-"+length, values, 1); cases++;
				}
			}
			final float[] special=new float[33];
			special[16]=Float.NaN;
			failures+=checkSum("nan-second-block", special, Double.NaN); cases++;
			special[16]=Float.POSITIVE_INFINITY;
			failures+=checkSum("positive-infinity", special, Double.POSITIVE_INFINITY); cases++;
			special[16]=Float.NEGATIVE_INFINITY;
			failures+=checkSum("negative-infinity", special, Double.NEGATIVE_INFINITY); cases++;
			special[8]=Float.POSITIVE_INFINITY;
			failures+=checkSum("opposing-infinities", special, Double.NaN); cases++;
			Arrays.fill(special, -0.0f);
			failures+=checkSum("negative-zeros", special, 0.0); cases++;
			special[8]=+0.0f;
			failures+=checkSum("mixed-zeros", special, 0.0); cases++;
			System.out.println("FLOAT_REDUCTION_RESULT simd="+simd+" cases="+cases+
					" failures="+failures);
			if(failures>0){throw new AssertionError("Floating reduction failures="+failures);}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Checks a sum against an exact finite value or a specified nonfinite result.
	 * @param label Diagnostic case name
	 * @param values Input array, possibly null
	 * @param expected Literal or closed-form expected result
	 * @return One on mismatch, otherwise zero */
	private static int checkSum(final String label, final float[] values, final double expected){
		final double actual=Vector.sum(values);
		if(Double.doubleToLongBits(actual)==Double.doubleToLongBits(expected)){return 0;}
		System.out.println("FAIL sum "+label+" expected="+expected+" actual="+actual);
		return 1;
	}

	/** Checks scalar maximum-selection semantics, preserving signed-zero ties.
	 * NaN payload bits are not part of this test's contract.
	 * @param label Diagnostic case name
	 * @param values Input array, possibly null
	 * @param expected Expected selected value
	 * @return One on mismatch, otherwise zero */
	private static int checkMax(final String label, final float[] values, final float expected){
		final float actual=Vector.max(values);
		if(Float.floatToIntBits(actual)==Float.floatToIntBits(expected)){return 0;}
		System.out.println("FAIL max "+label+" expected="+expected+" actual="+actual);
		return 1;
	}

	/** Prevents construction of this standalone test. */
	private VectorFloatReductionTest(){}

	/** Lengths around the facade threshold, vector groups and unrolled-loop tails. */
	private static final int[] LENGTHS={1, 8, 15, 16, 17, 23, 24, 25, 31, 32, 33, 64, 65, 127, 128, 129};
}
