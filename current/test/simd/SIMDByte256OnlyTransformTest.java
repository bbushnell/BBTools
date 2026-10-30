package test.simd;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

/** Checks edited SIMDByte256_only byte transforms directly.
 * @author Sucrose */
public final class SIMDByte256OnlyTransformTest{

	/** Runs direct historical-backend fixtures.
	 * @param args Unused
	 * @throws Exception If reflective backend calls fail
	 */
	public static void main(final String[] args) throws Exception{
		final Method addAndCapMin=method("addAndCapMin", byte[].class, byte.class, int.class);
		final Method dotDashXToN=method("dotDashXToN", byte[].class);
		final Method toUpperCase=method("toUpperCase", byte[].class);
		int cases=0, failures=0;
		failures+=checkDotDashXToN(dotDashXToN); cases++;
		failures+=checkToUpperCase(toUpperCase); cases++;
		failures+=checkAddAndCapMin(addAndCapMin); cases++;
		System.out.println("BYTE256_ONLY_RESULT cases="+cases+" failures="+failures);
		if(failures!=0){throw new AssertionError("SIMDByte256_only failures="+failures);}
	}

	/** Gets an accessible SIMDByte256_only method.
	 * @param name Method name
	 * @param types Parameter types
	 * @return Accessible method
	 * @throws Exception If the method lookup fails
	 */
	private static Method method(final String name, final Class<?>... types) throws Exception{
		final Method m=Class.forName("simd.SIMDByte256_only").getDeclaredMethod(name, types);
		m.setAccessible(true);
		return m;
	}

	/** Checks dot/dash/X/x/n conversion to N.
	 * @param method Backend method
	 * @return One on mismatch, otherwise zero
	 * @throws Exception If the reflective call fails
	 */
	private static int checkDotDashXToN(final Method method) throws Exception{
		final byte[] actual=".-XxnACGTumrwsyk?!0123456789.-XxnACGTACG".
				getBytes(StandardCharsets.US_ASCII);
		final byte[] expected=actual.clone();
		for(int i=0; i<expected.length; i++){
			if(expected[i]=='.' || expected[i]=='-' || expected[i]=='X' ||
					expected[i]=='x' || expected[i]=='n'){
				expected[i]='N';
			}
		}
		final boolean actualReturn=(Boolean)method.invoke(null, (Object)actual);
		if(actualReturn && Arrays.equals(expected, actual)){return 0;}
		System.out.println("FAIL dotDashXToN expected="+Arrays.toString(expected)+
				" actual="+Arrays.toString(actual)+" actualReturn="+actualReturn);
		return 1;
	}

	/** Checks ASCII lower-case conversion and punctuation preservation.
	 * @param method Backend method
	 * @return One on mismatch, otherwise zero
	 * @throws Exception If the reflective call fails
	 */
	private static int checkToUpperCase(final Method method) throws Exception{
		final byte[] actual="a.c-gtnxACGT?!0123456789a.c-gtnxACGT!?az".
				getBytes(StandardCharsets.US_ASCII);
		final byte[] expected=actual.clone();
		for(int i=0; i<expected.length; i++){
			if(expected[i]>='a' && expected[i]<='z'){expected[i]=(byte)(expected[i]-32);}
		}
		final boolean actualReturn=(Boolean)method.invoke(null, (Object)actual);
		if(actualReturn && Arrays.equals(expected, actual)){return 0;}
		System.out.println("FAIL toUpperCase expected="+Arrays.toString(expected)+
				" actual="+Arrays.toString(actual)+" actualReturn="+actualReturn);
		return 1;
	}

	/** Checks min capping, extreme thresholds and the empty-array sentinel.
	 * @param method Backend method
	 * @return Number of mismatches
	 * @throws Exception If a reflective call fails
	 */
	private static int checkAddAndCapMin(final Method method) throws Exception{
		int failures=0;
		final int[] lengths={0, 40, 40, 40, 40};
		final int[] caps={4, Integer.MAX_VALUE, Integer.MIN_VALUE, 255, -129};
		final byte[] deltas={(byte)-17, (byte)-128, (byte)127, (byte)-1, (byte)1};
		for(int i=0; i<caps.length; i++){
			final byte[] actual=makeSigned(lengths[i]);
			final byte[] expected=actual.clone();
			final byte expectedMin=scalarAddAndCapMin(expected, deltas[i], caps[i]);
			final byte actualMin=(Byte)method.invoke(null, actual, deltas[i], caps[i]);
			if(expectedMin!=actualMin || !Arrays.equals(expected, actual)){
				System.out.println("FAIL addAndCapMin index="+i+
						" expectedMin="+expectedMin+" actualMin="+actualMin+
						" expected="+Arrays.toString(expected)+
						" actual="+Arrays.toString(actual));
				failures++;
			}
		}
		return failures;
	}

	/** Applies the scalar addAndCapMin contract.
	 * @param array Bytes to modify
	 * @param delta Delta to add
	 * @param cap Minimum cap
	 * @return Minimum after addition, before capping
	 */
	private static byte scalarAddAndCapMin(final byte[] array, final byte delta, final int cap){
		int min=127;
		for(int i=0; i<array.length; i++){
			final int b=array[i]+delta;
			min=Math.min(min, b);
			array[i]=(byte)Math.max(cap, b);
		}
		return (byte)min;
	}

	/** Makes a signed byte array.
	 * @param length Array length
	 * @return New array
	 */
	private static byte[] makeSigned(final int length){
		final byte[] out=new byte[length];
		for(int i=0; i<length; i++){out[i]=(byte)(i*19-117);}
		return out;
	}

	/** Prevents construction of this standalone test. */
	private SIMDByte256OnlyTransformTest(){}
}
