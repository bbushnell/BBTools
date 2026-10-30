package test.simd;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import shared.Shared;
import simd.Vector;
import structures.ByteBuilder;

/** Checks byte transform, quality-helper and fake-quality append parity.
 * The lowercase/IUPAC quality case records the current SIMDByte256#002 policy
 * divergence separately from ordinary uppercase read parity.
 * @author Sucrose */
public final class VectorByteTransformTest{

	/** Runs finite byte transform fixtures.
	 * @param args simd=t or simd=f
	 */
	public static void main(final String[] args){
		final boolean simd=parseSimd(args);
		if(simd && !Vector.simd128){
			throw new IllegalStateException("SIMD test requires 128/256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			for(final int length : LENGTHS){
				final byte[] bases=makeBases(length);
				failures+=checkUnary("uToT", simd, bases, a -> {Vector.uToT(a); return true;}); cases++;
				failures+=checkUnary("lowerCaseToN", simd, bases, Vector::lowerCaseToN); cases++;
				failures+=checkUnary("dotDashXToN", simd, bases, Vector::dotDashXToN); cases++;
				failures+=checkUnary("toUpperCase", simd, bases, Vector::toUpperCase); cases++;
				failures+=checkAdd(length, simd); cases++;
				failures+=checkAddAndCap(length, simd); cases++;
				failures+=checkQuality(length, simd, false); cases++;
				failures+=checkQuality(length, simd, true); cases++;
				failures+=checkAppendFake(length, simd); cases++;
			}
			if(simd){
				failures+=checkQualityDivergence(false); cases++;
				failures+=checkQualityDivergence(true); cases++;
			}
			failures+=checkAddAndCapExtreme(simd); cases++;
			failures+=checkAddAndCapGrid(simd); cases++;
			System.out.println("BYTE_TRANSFORM_RESULT simd="+simd+" cases="+cases+
					" failures="+failures);
			if(failures!=0){throw new AssertionError("Byte transform failures="+failures);}
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

	/** Compares one scalar and selected-mode unary transform result.
	 * @param label Transform label
	 * @param simd Whether the tested run should force SIMD
	 * @param input Input array
	 * @param transform Transform function
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkUnary(final String label, final boolean simd, final byte[] input,
			final ByteTransform transform){
		final byte[] expected=input.clone();
		final byte[] actual=input.clone();
		Shared.SIMD=false;
		final boolean expectedReturn=transform.apply(expected);
		Shared.SIMD=simd;
		final boolean actualReturn=transform.apply(actual);
		if(expectedReturn==actualReturn && Arrays.equals(expected, actual)){return 0;}
		System.out.println("FAIL "+label+" length="+input.length+
				" expectedReturn="+expectedReturn+" actualReturn="+actualReturn+
				" expected="+Arrays.toString(expected)+" actual="+Arrays.toString(actual));
		return 1;
	}

	/** Checks signed byte addition.
	 * @param length Array length
	 * @param simd Whether to force SIMD
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkAdd(final int length, final boolean simd){
		final byte[] expected=makeSigned(length), actual=expected.clone();
		Shared.SIMD=false;
		Vector.add(expected, (byte)17);
		Shared.SIMD=simd;
		Vector.add(actual, (byte)17);
		return checkArray("add", expected, actual);
	}

	/** Checks addAndCapMin output bytes and returned pre-cap minimum.
	 * @param length Array length
	 * @param simd Whether to force SIMD
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkAddAndCap(final int length, final boolean simd){
		final byte[] expected=makeSigned(length), actual=expected.clone();
		Shared.SIMD=false;
		final byte expectedMin=Vector.addAndCapMin(expected, (byte)29, 4);
		Shared.SIMD=simd;
		final byte actualMin=Vector.addAndCapMin(actual, (byte)29, 4);
		if(expectedMin==actualMin){return checkArray("addAndCapMin", expected, actual);}
		System.out.println("FAIL addAndCapMin-min length="+length+" expected="+
				expectedMin+" actual="+actualMin);
		return 1;
	}

	/** Checks addAndCapMin with caps outside the signed-byte range.
	 * @param simd Whether to force SIMD
	 * @return Number of mismatches
	 */
	private static int checkAddAndCapExtreme(final boolean simd){
		int failures=0;
		final int[] caps={Integer.MAX_VALUE, Integer.MIN_VALUE, 255, -129};
		final byte[] deltas={(byte)-128, (byte)127, (byte)-1, (byte)1};
		for(int i=0; i<caps.length; i++){
			final byte[] expected=makeSigned(40), actual=expected.clone();
			Shared.SIMD=false;
			final byte expectedMin=Vector.addAndCapMin(expected, deltas[i], caps[i]);
			Shared.SIMD=simd;
			final byte actualMin=Vector.addAndCapMin(actual, deltas[i], caps[i]);
			if(expectedMin==actualMin){
				failures+=checkArray("addAndCapMin-extreme-"+i, expected, actual);
			}else{
				System.out.println("FAIL addAndCapMin-extreme-min index="+i+
						" expected="+expectedMin+" actual="+actualMin);
				failures++;
			}
		}
		return failures;
	}

	/** Checks every signed-byte/delta pair against literal integer arithmetic.
	 * Selected caps include both int extremes and all clamp-region boundaries.
	 * A second all-127 pattern checks the pre-cap minimum's initial127 ceiling.
	 * Both patterns include 32-byte, eight-byte and scalar-tail work.
	 * @param simd Whether to force SIMD
	 * @return One on the first mismatch, otherwise zero */
	private static int checkAddAndCapGrid(final boolean simd){
		Shared.SIMD=simd;
		int arrays=0, values=0;
		for(int pattern=0; pattern<2; pattern++){
			final int length=(pattern==0 ? 297 : 41);
			final byte[] actual=new byte[length], expected=new byte[length];
			for(int delta=Byte.MIN_VALUE; delta<=Byte.MAX_VALUE; delta++){
				for(final int cap : GRID_CAPS){
					int expectedMin=127;
					for(int i=0; i<length; i++){
						final int input=(pattern==0 && i<256 ? i-128 : 127);
						actual[i]=(byte)input;
						final int sum=input+delta;
						expected[i]=(byte)Math.max(cap, sum);
						expectedMin=Math.min(expectedMin, sum);
					}
					final byte observedMin=Vector.addAndCapMin(actual, (byte)delta, cap);
					if(observedMin!=(byte)expectedMin || !Arrays.equals(actual, expected)){
						System.out.println("FAIL cap-grid pattern="+pattern+" delta="+delta+
								" cap="+cap+" expectedMin="+(byte)expectedMin+" actualMin="+observedMin);
						return 1;
					}
					arrays++;
					values+=length;
				}
			}
		}
		System.out.println("CAP_GRID_RESULT simd="+simd+" arrays="+arrays+" values="+values);
		return 0;
	}

	/** Checks ordinary uppercase read quality helpers.
	 * @param length Array length
	 * @param simd Whether to force SIMD
	 * @param cap Whether to test capQuality instead of applyQualOffset
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkQuality(final int length, final boolean simd, final boolean cap){
		final byte[] bases=new byte[length];
		final byte[] expected=new byte[length], actual=new byte[length];
		for(int i=0; i<length; i++){
			bases[i]=UPPERCASE_READ_BASES[i%UPPERCASE_READ_BASES.length];
			expected[i]=actual[i]=(byte)(i%60);
		}
		Shared.SIMD=false;
		if(cap){Vector.capQuality(expected, bases);}
		else{Vector.applyQualOffset(expected, bases, 5);}
		Shared.SIMD=simd;
		if(cap){Vector.capQuality(actual, bases);}
		else{Vector.applyQualOffset(actual, bases, 5);}
		return checkArray(cap ? "capQuality" : "applyQualOffset", expected, actual);
	}

	/** Checks the current lowercase/IUPAC/punctuation quality divergence.
	 * @param cap Whether to test capQuality instead of applyQualOffset
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkQualityDivergence(final boolean cap){
		final byte[] bases="nMR.-xacgtunMR.-x".getBytes(StandardCharsets.US_ASCII);
		final byte[] scalar=new byte[bases.length], vector=new byte[bases.length];
		Arrays.fill(scalar, (byte)20);
		Arrays.fill(vector, (byte)20);
		Shared.SIMD=false;
		if(cap){Vector.capQuality(scalar, bases);}
		else{Vector.applyQualOffset(scalar, bases, 5);}
		Shared.SIMD=true;
		if(cap){Vector.capQuality(vector, bases);}
		else{Vector.applyQualOffset(vector, bases, 5);}
		final byte[] expectedScalar=(cap ? QUALITY_CAP_SCALAR : QUALITY_OFFSET_SCALAR);
		final byte[] expectedVector=(cap ? QUALITY_CAP_VECTOR : QUALITY_OFFSET_VECTOR);
		if(Arrays.equals(expectedScalar, scalar) && Arrays.equals(expectedVector, vector)){return 0;}
		System.out.println("FAIL quality-divergence cap="+cap+
				" expectedScalar="+Arrays.toString(expectedScalar)+
				" actualScalar="+Arrays.toString(scalar)+
				" expectedVector="+Arrays.toString(expectedVector)+
				" actualVector="+Arrays.toString(vector));
		return 1;
	}

	/** Checks Vector.appendFake and existing ByteBuilder content preservation.
	 * @param length Array length
	 * @param simd Whether to force SIMD
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkAppendFake(final int length, final boolean simd){
		final byte[] bases=makeBases(length);
		final ByteBuilder expected=new ByteBuilder(1), actual=new ByteBuilder(1);
		expected.append("xy");
		actual.append("xy");
		Shared.SIMD=false;
		Vector.appendFake(bases, expected, 40, 2);
		Shared.SIMD=simd;
		Vector.appendFake(bases, actual, 40, 2);
		final byte[] e=Arrays.copyOf(expected.array, expected.length);
		final byte[] a=Arrays.copyOf(actual.array, actual.length);
		if(Arrays.equals(e, a)){return 0;}
		System.out.println("FAIL appendFake length="+length+" expected="+
				Arrays.toString(e)+" actual="+Arrays.toString(a));
		return 1;
	}

	/** Compares byte arrays and reports mismatches.
	 * @param label Test label
	 * @param expected Expected bytes
	 * @param actual Actual bytes
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkArray(final String label, final byte[] expected, final byte[] actual){
		if(Arrays.equals(expected, actual)){return 0;}
		System.out.println("FAIL "+label+" length="+expected.length+
				" expected="+Arrays.toString(expected)+" actual="+Arrays.toString(actual));
		return 1;
	}

	/** Makes a positive-ASCII base and punctuation array.
	 * @param length Array length
	 * @return New array
	 */
	private static byte[] makeBases(final int length){
		final byte[] out=new byte[length];
		for(int i=0; i<length; i++){out[i]=BASES[i%BASES.length];}
		return out;
	}

	/** Makes a signed byte array for arithmetic tests.
	 * @param length Array length
	 * @return New array
	 */
	private static byte[] makeSigned(final int length){
		final byte[] out=new byte[length];
		for(int i=0; i<length; i++){out[i]=(byte)(i*19-117);}
		return out;
	}

	/** Byte-array transform function. */
	private interface ByteTransform{
		/** Applies this transform in-place.
		 * @param array Array to transform
		 * @return Transform result
		 */
		boolean apply(byte[] array);
	}

	/** Prevents construction of this standalone test. */
	private VectorByteTransformTest(){}

	/** Lengths spanning 8-byte, 16-byte, 32-byte and tail thresholds. */
	private static final int[] LENGTHS={0, 1, 7, 8, 9, 15, 16, 17, 23, 24, 25, 31,
			32, 33, 39, 40, 41, 47, 48, 49};
	/** Caps spanning ordinary values, byte boundaries and int extremes. */
	private static final int[] GRID_CAPS={Integer.MIN_VALUE, -257, -256, -255, -129,
			-128, -127, -1, 0, 1, 126, 127, 128, 254, 255, 256, 257, Integer.MAX_VALUE};
	/** Uppercase read bases used by ordinary quality cases. */
	private static final byte[] UPPERCASE_READ_BASES={'A', 'C', 'G', 'T', 'N'};
	/** Positive ASCII bases and punctuation covered by transform fixtures. */
	private static final byte[] BASES={'A', 'C', 'G', 'T', 'U', 'N', 'a', 'c', 'g', 't',
			'u', 'n', 'M', 'R', 'W', 'S', 'Y', 'K', '.', '-', 'X', 'x', '!', '?'};
	/** Scalar applyQualOffset result for "nMR.-xacgtunMR.-x". */
	private static final byte[] QUALITY_OFFSET_SCALAR={0, 0, 0, 0, 0, 0, 25, 25, 25,
			25, 25, 0, 0, 0, 0, 0, 0};
	/** SIMD applyQualOffset result for 16 vector lanes plus one scalar tail. */
	private static final byte[] QUALITY_OFFSET_VECTOR={25, 25, 25, 25, 25, 25, 25,
			25, 25, 25, 25, 25, 25, 25, 25, 25, 0};
	/** Scalar capQuality result for "nMR.-xacgtunMR.-x". */
	private static final byte[] QUALITY_CAP_SCALAR={0, 0, 0, 0, 0, 0, 20, 20, 20,
			20, 20, 0, 0, 0, 0, 0, 0};
	/** SIMD capQuality result for 16 vector lanes plus one scalar tail. */
	private static final byte[] QUALITY_CAP_VECTOR={20, 20, 20, 20, 20, 20, 20, 20,
			20, 20, 20, 20, 20, 20, 20, 20, 0};
}
