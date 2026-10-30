package test.simd;

import java.lang.reflect.Method;
import java.util.Arrays;

import shared.Shared;
import simd.Vector;

/** Checks byte predicate parity across ASCII scalar tails and vector chunks.
 * @author Sucrose */
public final class VectorBytePredicateTest{

	/** Runs finite predicate fixtures.
	 * @param args simd=t or simd=f
	 * @throws Exception If reflective historical-backend calls fail
	 */
	public static void main(final String[] args) throws Exception{
		final boolean simd=parseSimd(args);
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD test requires 256-bit Vector support");
		}
		final Method oldIsProtein=method("isProtein");
		final Method oldIsNucleotide=method("isNucleotide");
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			Shared.SIMD=simd;
			failures+=checkNull(oldIsProtein, oldIsNucleotide); cases++;
			for(final int length : LENGTHS){
				for(final int ascii : ASCII_VALUES){
					for(final int pos : positions(length)){
						final byte[] input=asciiCase(length, (byte)ascii, pos);
						final String label="ascii-"+ascii+"-"+pos;
						failures+=checkVector(input, label);
						cases++;
						failures+=checkHistorical(oldIsProtein, oldIsNucleotide, input, label);
						cases++;
					}
				}
				failures+=checkVector(nucleotideCase(length), "nucleotide"); cases++;
				failures+=checkHistorical(oldIsProtein, oldIsNucleotide,
						asciiCase(length, (byte)'f', Math.max(0, length-1)),
						"historical-f-"+length);
				cases++;
				failures+=checkHistorical(oldIsProtein, oldIsNucleotide,
						asciiCase(length, (byte)'Z', 0), "historical-Z-"+length);
				cases++;
			}
			System.out.println("BYTE_PREDICATE_RESULT simd="+simd+" cases="+cases+
					" failures="+failures+" ascii_values="+ASCII_VALUES.length);
			if(failures!=0){throw new AssertionError("Byte predicate failures="+failures);}
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

	/** Gets an accessible SIMDByte256_only predicate.
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
	 * @param isProtein Reflective historical isProtein method
	 * @param isNucleotide Reflective historical isNucleotide method
	 * @return One on mismatch, otherwise zero
	 * @throws Exception If a reflective call fails
	 */
	private static int checkNull(final Method isProtein, final Method isNucleotide)
			throws Exception{
		final boolean vectorProtein=Vector.isProtein(null);
		final boolean vectorNucleotide=Vector.isNucleotide(null);
		final boolean oldProtein=(Boolean)isProtein.invoke(null, (Object)null);
		final boolean oldNucleotide=(Boolean)isNucleotide.invoke(null, (Object)null);
		if(vectorProtein && vectorNucleotide && !oldProtein && oldNucleotide){return 0;}
		System.out.println("FAIL null vectorProtein="+vectorProtein+
				" vectorNucleotide="+vectorNucleotide+" oldProtein="+oldProtein+
				" oldNucleotide="+oldNucleotide);
		return 1;
	}

	/** Checks the public Vector predicates.
	 * @param input Input sequence
	 * @param label Case label
	 * @return One on mismatch, otherwise zero
	 */
	private static int checkVector(final byte[] input, final String label){
		final boolean expectedProtein=expectedProtein(input);
		final boolean expectedNucleotide=expectedNucleotide(input);
		final boolean actualProtein=Vector.isProtein(input.clone());
		final boolean actualNucleotide=Vector.isNucleotide(input.clone());
		if(expectedProtein==actualProtein && expectedNucleotide==actualNucleotide){return 0;}
		System.out.println("FAIL vector "+label+" length="+input.length+
				" expectedProtein="+expectedProtein+" actualProtein="+actualProtein+
				" expectedNucleotide="+expectedNucleotide+
				" actualNucleotide="+actualNucleotide+
				" input="+Arrays.toString(input));
		return 1;
	}

	/** Checks the historical SIMDByte256_only predicates directly.
	 * @param isProtein Reflective isProtein method
	 * @param isNucleotide Reflective isNucleotide method
	 * @param input Input sequence
	 * @param label Case label
	 * @return One on mismatch, otherwise zero
	 * @throws Exception If a reflective call fails
	 */
	private static int checkHistorical(final Method isProtein, final Method isNucleotide,
			final byte[] input, final String label) throws Exception{
		final boolean expectedProtein=expectedProtein(input);
		final boolean expectedNucleotide=expectedNucleotide(input);
		final boolean actualProtein=(Boolean)isProtein.invoke(null, (Object)input.clone());
		final boolean actualNucleotide=(Boolean)isNucleotide.invoke(null, (Object)input.clone());
		if(expectedProtein==actualProtein && expectedNucleotide==actualNucleotide){return 0;}
		System.out.println("FAIL historical "+label+" length="+input.length+
				" expectedProtein="+expectedProtein+" actualProtein="+actualProtein+
				" expectedNucleotide="+expectedNucleotide+
				" actualNucleotide="+actualNucleotide+
				" input="+Arrays.toString(input));
		return 1;
	}

	/** Builds one valid nucleotide sequence.
	 * @param length Sequence length
	 * @return New sequence
	 */
	private static byte[] nucleotideCase(final int length){
		final byte[] out=new byte[length];
		for(int i=0; i<length; i++){out[i]=IUPAC_BASES[i%IUPAC_BASES.length];}
		return out;
	}

	/** Builds a nucleotide sequence with one ASCII byte inserted.
	 * @param length Sequence length
	 * @param ascii ASCII byte
	 * @param pos ASCII byte position
	 * @return New sequence
	 */
	private static byte[] asciiCase(final int length, final byte ascii, final int pos){
		final byte[] out=nucleotideCase(length);
		if(length>0){out[pos]=ascii;}
		return out;
	}

	/** Returns selected ASCII insertion positions.
	 * @param length Sequence length
	 * @return Positions that cover initial, 64-bit, 256-bit and scalar-tail lanes
	 */
	private static int[] positions(final int length){
		if(length<1){return new int[] {0};}
		return new int[] {0, Math.min(7, length-1), Math.min(8, length-1),
				Math.min(31, length-1), Math.min(32, length-1), length-1};
	}

	/** Applies the Vector.isProtein scalar table contract independently.
	 * @param array Input sequence
	 * @return True when any byte is an amino-acid-only symbol
	 */
	private static boolean expectedProtein(final byte[] array){
		for(final byte b : array){
			if(isAminoOnly(b) || b=='.' || b=='-'){return true;}
		}
		return false;
	}

	/** Applies the Vector.isNucleotide scalar table contract independently.
	 * @param array Input sequence
	 * @return True when all bytes are IUPAC nucleotide bases
	 */
	private static boolean expectedNucleotide(final byte[] array){
		for(final byte b : array){
			if(!isIupacBase(b)){return false;}
		}
		return true;
	}

	/** Tests one amino-acid symbol that is not an extended IUPAC base.
	 * @param b Byte to test
	 * @return True for E/F/I/J/L/O/P/Q/Z ignoring case, plus stop '*'
	 */
	private static boolean isAminoOnly(final byte b){
		final int x=b&~32;
		return x=='E' || x=='F' || x=='I' || x=='J' || x=='L' || x=='O' ||
				x=='P' || x=='Q' || x=='Z' || b=='*';
	}

	/** Tests one extended IUPAC base.
	 * @param b Byte to test
	 * @return True for an IUPAC base, ignoring case
	 */
	private static boolean isIupacBase(final byte b){
		switch(b&~32){
			case 'A':
			case 'B':
			case 'C':
			case 'D':
			case 'G':
			case 'H':
			case 'K':
			case 'M':
			case 'N':
			case 'R':
			case 'S':
			case 'T':
			case 'U':
			case 'V':
			case 'W':
			case 'X':
			case 'Y':
				return true;
			default:
				return false;
		}
	}

	/** Lengths around empty, scalar, 64-bit, 256-bit and 256+64-bit tails. */
	private static final int[] LENGTHS={0, 1, 7, 8, 9, 31, 32, 33, 40, 41, 63, 64, 65};
	/** Baseline bytes for positions not under test; all are extended IUPAC bases. */
	private static final byte[] IUPAC_BASES={'A', 'c', 'G', 't', 'U', 'n', 'R', 'y',
			'K', 'm', 'S', 'w', 'B', 'd', 'H', 'v'};
	/** All nonnegative ASCII bytes; negative Java bytes remain outside the 128-entry tables. */
	private static final int[] ASCII_VALUES=new int[128];
	static{
		for(int i=0; i<ASCII_VALUES.length; i++){ASCII_VALUES[i]=i;}
	}

	/** Prevents construction of this standalone test. */
	private VectorBytePredicateTest(){}
}
