package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks sparse vector facades against literal bSet-indexed expectations.
 * Full-width sparse rows must still use bSet when source and destination
 * lengths match, including when sparse SIMD is disabled per call.
 * @author Sucrose */
public final class VectorSparseFmaTest{

	/** Runs sparse fma/addProduct checks with identity and swapped indices.
	 * @param args Optional simd=t or simd=f selector */
	public static void main(final String[] args){
		boolean simd=false;
		for(final String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(simd && !Vector.simd128){
			throw new IllegalStateException("Sparse fma checks require 128/256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD, oldSimdFma=Shared.SIMD_FMA;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			Shared.SIMD_FMA=true;
			final int[] identity=makeIdentityIndices();
			final int[] swapped=makeSwappedIndices();
			failures+=checkFma("identity-allow", identity, true); cases++;
			failures+=checkFma("identity-disallow", identity, false); cases++;
			failures+=checkFma("swapped-allow", swapped, true); cases++;
			failures+=checkFma("swapped-disallow", swapped, false); cases++;
			failures+=checkAddProduct("identity", identity); cases++;
			failures+=checkAddProduct("swapped", swapped); cases++;
			System.out.println("SPARSE_FMA_RESULT simd="+simd+" cases="+cases+
					" failures="+failures);
			if(failures>0){throw new AssertionError("Sparse fma/addProduct failures="+failures);}
		}finally{
			Shared.SIMD_FMA=oldSimdFma;
			Shared.SIMD=oldSimd;
		}
	}

	/** Compares one sparse fma against an independent bSet-indexed sum.
	 * @param label Diagnostic case name
	 * @param indices Indices into the weights array
	 * @param allowSimd Per-call sparse SIMD permission
	 * @return One on mismatch, otherwise zero */
	private static int checkFma(final String label, final int[] indices, final boolean allowSimd){
		final float[] values=makeValues();
		final float[] weights=makeWeights();
		final float expected=expectedFma(values, weights, indices);
		final float actual=Vector.fma(values, weights, indices, 1, allowSimd);
		if(Float.floatToIntBits(actual)==Float.floatToIntBits(expected)){return 0;}
		System.out.println("FAIL fma "+label+" expected="+expected+" actual="+actual);
		return 1;
	}

	/** Compares one sparse addProduct against independent per-lane updates.
	 * @param label Diagnostic case name
	 * @param indices Indices into the increment array
	 * @return One on mismatch, otherwise zero */
	private static int checkAddProduct(final String label, final int[] indices){
		final float[] actual=makeValues();
		final float[] expected=makeValues();
		final float[] increments=makeWeights();
		Vector.addProduct(actual, increments, indices, MULT, 1);
		for(int i=0; i<expected.length; i++){expected[i]+=increments[indices[i]]*MULT;}
		if(Arrays.equals(actual, expected)){return 0;}
		System.out.println("FAIL addProduct "+label+" expected="+Arrays.toString(expected)+
				" actual="+Arrays.toString(actual));
		return 1;
	}

	/** Computes the sparse dot-product contract without using the Vector facade.
	 * @param values Multipliers
	 * @param weights Indexed values
	 * @param indices Indices into {@code weights}
	 * @return Sum of {@code values[i]*weights[indices[i]]} */
	private static float expectedFma(final float[] values, final float[] weights, final int[] indices){
		float sum=0;
		for(int i=0; i<values.length; i++){sum+=values[i]*weights[indices[i]];}
		return sum;
	}

	/** Creates a 16-element identity sparse-index vector.
	 * @return New identity index vector */
	private static int[] makeIdentityIndices(){
		final int[] indices=new int[LENGTH];
		for(int i=0; i<indices.length; i++){indices[i]=i;}
		return indices;
	}

	/** Creates a 16-element nonidentity sparse-index vector.
	 * @return New vector with its first two indices swapped */
	private static int[] makeSwappedIndices(){
		final int[] indices=makeIdentityIndices();
		indices[0]=1;
		indices[1]=0;
		return indices;
	}

	/** Creates multiplicands whose first two lanes expose swapped-index use.
	 * @return New value array */
	private static float[] makeValues(){
		final float[] values=new float[LENGTH];
		for(int i=0; i<values.length; i++){values[i]=i+1;}
		return values;
	}

	/** Creates indexed weights whose first two lanes are distinguishable.
	 * @return New weight array */
	private static float[] makeWeights(){
		final float[] weights=new float[LENGTH];
		for(int i=0; i<weights.length; i++){weights[i]=100+3*i;}
		return weights;
	}

	/** Prevents construction of this standalone test. */
	private VectorSparseFmaTest(){}

	/** Full vector length that triggers Vector's equal-length sparse shortcuts. */
	private static final int LENGTH=16;

	/** Exact dyadic multiplier for addProduct checks. */
	private static final float MULT=1.5f;
}
