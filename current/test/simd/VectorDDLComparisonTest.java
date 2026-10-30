package test.simd;

import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks exact scalar/SIMD parity for DynamicDemiLog bucket comparisons.
 * Includes empty and lane-boundary lengths plus unsigned short values that cross
 * Java's signed-short boundary.
 * @author Sucrose */
public final class VectorDDLComparisonTest{

	/** Runs finite DDL comparison fixtures.
	 * @param args simd=t or simd=f
	 */
	public static void main(final String[] args){
		final boolean simd=parseSimd(args);
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD test requires 256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0, oracleChecksum=0;
		try{
			Shared.SIMD=simd;
			for(final int length : LENGTHS){
				for(int pattern=0; pattern<PATTERNS; pattern++){
					final char[] a=makeValues(length, pattern, false);
					final char[] b=makeValues(length, pattern, true);
					final short[] sa=Vector.toShortArray(a);
					final short[] sb=Vector.toShortArray(b);
					final int[] expected=oracle(a, b);
					failures+=check("char-char", length, pattern, expected, Vector.compareDDL(a, b));
					cases++;
					failures+=check("short-char", length, pattern, expected, Vector.compareDDL(sa, b));
					cases++;
					failures+=check("short-short", length, pattern, expected, Vector.compareDDL(sa, sb));
					cases++;
					oracleChecksum=31*oracleChecksum+expected[0]*3+expected[1]*5+expected[2]*7+expected[3]*11;
				}
			}
			System.out.println("DDL_COMPARISON_RESULT simd="+simd+" cases="+cases+
					" failures="+failures+" oracle_checksum="+oracleChecksum);
			if(failures!=0){throw new AssertionError("DDL comparison failures="+failures);}
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Parses the SIMD selector.
	 * @param args Command-line arguments
	 * @return True for simd=t
	 */
	private static boolean parseSimd(final String[] args){
		boolean simd=false;
		for(final String arg : args){
			if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		return simd;
	}

	/** Builds one bounded DDL array.
	 * @param length Array length
	 * @param pattern Pattern selector
	 * @param second Whether to build the right operand
	 * @return New char array
	 */
	private static char[] makeValues(final int length, final int pattern, final boolean second){
		final char[] array=new char[length];
		for(int i=0; i<length; i++){
			final int selector=(i+pattern*3+(second ? 5 : 0))&7;
			final int value;
			if(pattern==0){value=0;}
			else if(pattern==1){value=selector==0 ? 0 : 1;}
			else if(pattern==2){value=VALUES[(i+selector)%VALUES.length];}
			else{value=VALUES[(i*3+selector)%VALUES.length];}
			array[i]=(char)value;
		}
		return array;
	}

	/** Computes the expected lower/equal/higher/both-empty tuple independently.
	 * @param a Left values
	 * @param b Right values
	 * @return int[]{lower, equal, higher, bothEmpty}
	 */
	private static int[] oracle(final char[] a, final char[] b){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<a.length; i++){
			final int ai=a[i], bi=b[i];
			if(ai==0 && bi==0){bothEmpty++;}
			else if(ai<bi){lower++;}
			else if(ai>bi){higher++;}
			else{equal++;}
		}
		return new int[] {lower, equal, higher, bothEmpty};
	}

	/** Compares one result vector.
	 * @param label Facade overload label
	 * @param length Input length
	 * @param pattern Pattern selector
	 * @param expected Expected tuple
	 * @param actual Actual tuple
	 * @return One on mismatch, otherwise zero
	 */
	private static int check(final String label, final int length, final int pattern,
			final int[] expected, final int[] actual){
		if(Arrays.equals(expected, actual)){return 0;}
		System.out.println("FAIL "+label+" length="+length+" pattern="+pattern+
				" expected="+Arrays.toString(expected)+" actual="+Arrays.toString(actual));
		return 1;
	}

	/** Prevents construction of this standalone test. */
	private VectorDDLComparisonTest(){}

	/** Lengths spanning scalar tails and 16-lane short vectors. */
	private static final int[] LENGTHS={0, 1, 15, 16, 17, 31, 32, 33, 2048};
	/** Values around empty buckets and signed-short boundaries. */
	private static final int[] VALUES={0, 1, 32767, 32768, 40000, 65535};
	/** Number of deterministic data patterns. */
	private static final int PATTERNS=4;
}
