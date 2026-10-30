package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks that integer sums widen before addition, including vector block reductions.
 * Expected sums use independent long accumulation and literal boundary controls.
 * This finite fixture is not a throughput benchmark.
 * @author Shinobu */
public final class VectorIntegerSumTest{

	/** Runs scalar or hardware-supported SIMD checks without hiding failed cases.
	 * @param args simd=t or simd=f
	 */
	public static void main(final String[] args){
		boolean simd=false;
		for(String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(simd && !Vector.simd128){
			throw new IllegalStateException("SIMD test requires the optional Vector module and 128/256-bit support");
		}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			for(int length : LENGTHS){
				for(int pattern=0; pattern<5; pattern++){
					final int[] values=makeValues(length, pattern);
					failures+=check(values, 0, length-1, true, pattern); cases++;
					failures+=check(values, 0, length-1, false, pattern); cases++;
					if(length>2){
						failures+=check(values, 1, length-2, false, pattern); cases++;
					}
					failures+=check(values, length, length-1, false, pattern); cases++;
				}
			}
			final int[] literal=new int[16];
			Arrays.fill(literal, Integer.MAX_VALUE);
			final long actual=Vector.sum(literal);
			if(actual!=34359738352L){
				System.out.println("FAIL literal16max expected=34359738352 actual="+actual);
				failures++;
			}
			cases++;
			System.out.println("INTEGER_SUM_RESULT simd="+simd+" cases="+cases+" failures="+failures);
			if(failures!=0){throw new AssertionError("Integer sums must widen before addition; failed cases="+failures);}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Produces bounded signed and unsigned-magnitude controls.
	 * @param length Array length
	 * @param pattern Zero/small/large-positive/large-negative/mixed-sign selector
	 * @return New array for one case
	 */
	private static int[] makeValues(final int length, final int pattern){
		final int[] values=new int[length];
		for(int i=0; i<length; i++){
			values[i]=(pattern==0 ? 0 : pattern==1 ? i+1 : pattern==2 ? Integer.MAX_VALUE :
					pattern==3 ? Integer.MIN_VALUE : (i%3==0 ? Integer.MAX_VALUE : -1000000000));
		}
		return values;
	}

	/** Compares one facade operation with a scalar long oracle.
	 * @param values Input values
	 * @param from Inclusive start
	 * @param to Inclusive end, from-1 for an empty range
	 * @param whole Whether to use the whole-array overload
	 * @param pattern Diagnostic input-pattern selector
	 * @return One on mismatch, otherwise zero
	 */
	private static int check(final int[] values, final int from, final int to,
			final boolean whole, final int pattern){
		long expected=0;
		for(int i=from; i<=to; i++){expected+=(long)values[i];}
		final long actual=(whole ? Vector.sum(values) : Vector.sum(values, from, to));
		if(actual==expected){return 0;}
		System.out.println("FAIL length="+values.length+" pattern="+pattern+" from="+from+
				" to="+to+" whole="+whole+" expected="+expected+" actual="+actual);
		return 1;
	}

	/** Prevents construction of this standalone test. */
	private VectorIntegerSumTest(){}

	/** Lengths spanning scalar thresholds, vector widths and residual tails. */
	private static final int[] LENGTHS={0, 1, 7, 8, 9, 15, 16, 17, 23, 24, 25, 31,
			32, 33, 65, 79, 80, 81};
}
