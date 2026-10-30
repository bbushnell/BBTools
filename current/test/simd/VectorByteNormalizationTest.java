package test.simd;

import java.lang.reflect.Method;
import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks ASCII base normalization at scalar, vector and tail boundaries.
 * @author Sucrose */
public final class VectorByteNormalizationTest{

	/** Runs finite live and historical iupacToN fixtures.
	 * @param args simd=t or simd=f
	 * @throws Exception If reflective historical-backend calls fail
	 */
	public static void main(final String[] args) throws Exception{
		final boolean simd=parseSimd(args);
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD test requires 256-bit Vector support");
		}
		final Method oldIupacToN=method("iupacToN");
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			failures+=checkNull(oldIupacToN); cases++;
			failures+=checkVector(new byte[0], "empty"); cases++;
			failures+=checkHistorical(oldIupacToN, new byte[0], "empty"); cases++;
			for(final int length : LENGTHS){
				for(final int ascii : ASCII_VALUES){
					for(final int pos : positions(length)){
						final byte[] input=asciiCase(length, (byte)ascii, pos);
						final String label="ascii-"+ascii+"-"+pos;
						failures+=checkVector(input, label); cases++;
						failures+=checkHistorical(oldIupacToN, input, label); cases++;
					}
				}
				final byte[] input=canonicalCase(length);
				failures+=checkVector(input, "canonical"); cases++;
				failures+=checkHistorical(oldIupacToN, input, "canonical"); cases++;
			}
			System.out.println("BYTE_NORMALIZATION_RESULT simd="+simd+" cases="+
					cases+" failures="+failures+" ascii_values="+ASCII_VALUES.length);
			if(failures!=0){throw new AssertionError("Byte normalization failures="+failures);}
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

	/** Gets an accessible SIMDByte256_only method.
	 * @param name Method name
	 * @return Accessible method
	 * @throws Exception If the method lookup fails
	 */
	private static Method method(final String name) throws Exception{
		final Method m=Class.forName("simd.SIMDByte256_only").getDeclaredMethod(name,
				byte[].class);
		m.setAccessible(true);
		return m;
	}

	/** Checks null input policies.
	 * @param iupacToN Reflective historical iupacToN method
	 * @return One on mismatch, otherwise zero
	 * @throws Exception If a reflective call fails
	 */
	private static int checkNull(final Method iupacToN) throws Exception{
		final boolean vector=Vector.iupacToN(null);
		final boolean old=(Boolean)iupacToN.invoke(null, (Object)null);
		if(vector && old){return 0;}
		System.out.println("FAIL null vector="+vector+" old="+old);
		return 1;
	}

	/** Checks the public Vector facade.
	 * @param input Input bases
	 * @param label Case label
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkVector(final byte[] input, final String label){
		final byte[] expected=expectedLive(input), actual=input.clone();
		final boolean actualReturn=Vector.iupacToN(actual);
		if(actualReturn && Arrays.equals(expected, actual)){return 0;}
		System.out.println("FAIL vector "+label+" length="+input.length+
				" expected="+Arrays.toString(expected)+" actual="+
				Arrays.toString(actual)+" actualReturn="+actualReturn+
				" input="+Arrays.toString(input));
		return 1;
	}

	/** Checks the historical SIMDByte256_only backend directly.
	 * @param iupacToN Reflective historical iupacToN method
	 * @param input Input bases
	 * @param label Case label
	 * @return One on mismatch, otherwise zero
	 * @throws Exception If a reflective call fails
	 */
	private static int checkHistorical(final Method iupacToN, final byte[] input,
			final String label) throws Exception{
		final byte[] expected=expectedHistorical(input), actual=input.clone();
		final boolean actualReturn=(Boolean)iupacToN.invoke(null, (Object)actual);
		if(actualReturn && Arrays.equals(expected, actual)){return 0;}
		System.out.println("FAIL historical "+label+" length="+input.length+
				" expected="+Arrays.toString(expected)+" actual="+
				Arrays.toString(actual)+" actualReturn="+actualReturn+
				" input="+Arrays.toString(input));
		return 1;
	}

	/** Builds expected live ambiguity-letter-only normalization.
	 * @param input Input bases
	 * @return New expected output
	 */
	private static byte[] expectedLive(final byte[] input){
		final byte[] out=input.clone();
		for(int i=0; i<out.length; i++){out[i]=normalizeLive(out[i]);}
		return out;
	}

	/** Builds expected historical full ACGTU normalization.
	 * @param input Input bases
	 * @return New expected output
	 */
	private static byte[] expectedHistorical(final byte[] input){
		final byte[] out=input.clone();
		for(int i=0; i<out.length; i++){out[i]=normalizeHistorical(out[i]);}
		return out;
	}

	/** Applies the live iupacToN contract independently for nonnegative ASCII.
	 * @param b Byte to normalize
	 * @return N for IUPAC ambiguity letters, otherwise the original byte
	 */
	private static byte normalizeLive(final byte b){
		switch(b){
			case 'B':
			case 'b':
			case 'D':
			case 'd':
			case 'H':
			case 'h':
			case 'K':
			case 'k':
			case 'M':
			case 'm':
			case 'N':
			case 'n':
			case 'R':
			case 'r':
			case 'S':
			case 's':
			case 'V':
			case 'v':
			case 'W':
			case 'w':
			case 'X':
			case 'x':
			case 'Y':
			case 'y':
				return 'N';
			default:
				return b;
		}
	}

	/** Applies the historical scalar baseToACGTN contract independently.
	 * @param b Byte to normalize
	 * @return A/C/G/T for valid DNA or U/u, otherwise N
	 */
	private static byte normalizeHistorical(final byte b){
		switch(b){
			case 'A':
			case 'a':
				return 'A';
			case 'C':
			case 'c':
				return 'C';
			case 'G':
			case 'g':
				return 'G';
			case 'T':
			case 't':
			case 'U':
			case 'u':
				return 'T';
			default:
				return 'N';
		}
	}

	/** Builds one canonical sequence with case, U/u, IUPAC and nonletters.
	 * @param length Sequence length
	 * @return New sequence
	 */
	private static byte[] canonicalCase(final int length){
		final byte[] out=new byte[length];
		for(int i=0; i<length; i++){out[i]=CANONICAL_BASES[i%CANONICAL_BASES.length];}
		return out;
	}

	/** Builds an A-only sequence with one ASCII byte inserted.
	 * @param length Sequence length
	 * @param ascii ASCII byte
	 * @param pos ASCII byte position
	 * @return New sequence
	 */
	private static byte[] asciiCase(final int length, final byte ascii, final int pos){
		final byte[] out=new byte[length];
		Arrays.fill(out, (byte)'A');
		out[pos]=ascii;
		return out;
	}

	/** Returns selected ASCII insertion positions.
	 * @param length Sequence length
	 * @return Positions that cover initial, 64-bit, 256-bit and scalar-tail lanes
	 */
	private static int[] positions(final int length){
		return new int[] {0, Math.min(7, length-1), Math.min(8, length-1),
				Math.min(31, length-1), Math.min(32, length-1), length-1};
	}

	/** Lengths around scalar, 64-bit, 256-bit and 256+64-bit tails. */
	private static final int[] LENGTHS={1, 7, 8, 9, 31, 32, 33, 40, 41, 63, 64, 65};
	/** Baseline bytes for positions not under test. */
	private static final byte[] CANONICAL_BASES={'A', 'c', 'G', 't', 'U', 'u', 'M',
			'r', 'N', 'x', 'Z', '.', '-', '?', ' '};
	/** All nonnegative ASCII bytes; negative Java bytes remain outside the 128-entry table. */
	private static final int[] ASCII_VALUES=new int[128];
	static{
		for(int i=0; i<ASCII_VALUES.length; i++){ASCII_VALUES[i]=i;}
	}

	/** Prevents construction of this standalone test. */
	private VectorByteNormalizationTest(){}
}
