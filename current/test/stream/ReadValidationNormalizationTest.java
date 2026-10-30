package test.stream;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import shared.Shared;
import simd.Vector;
import stream.Read;

/** Checks Read nucleotide normalization flag combinations.
 * @author Jean */
public final class ReadValidationNormalizationTest{

	/** Runs finite scalar and SIMD validation fixtures.
	 * @param args Command-line arguments, ignored */
	public static void main(final String[] args){
		final State oldState=new State();
		final boolean oldAminoIn=Shared.AMINO_IN;
		final boolean oldSimd=Shared.SIMD;
		int cases=0, failures=0;
		try{
			for(int mode=0; mode<2; mode++){
				final boolean simd=(mode==1);
				if(simd && !Vector.simd128){continue;}
				failures+=check("identity", simd, SEED, false, false, false, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("utot", simd, SEED, false, true, false, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("upper", simd, SEED, false, false, true, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("utot-upper", simd, SEED, false, true, true, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("lower", simd, SEED, false, false, false, true,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("lower-utot", simd, SEED, false, true, false, true,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("iupac", simd, SEED, true, false, false, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("iupac-utot", simd, SEED, true, true, false, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("iupac-upper", simd, SEED, true, false, true, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("iupac-utot-upper", simd, SEED, true, true, true, false,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("iupac-lower", simd, SEED, true, false, false, true,
						false, Read.IGNORE_JUNK, false); cases++;
				failures+=check("dotdash", simd, SEED, false, false, false, false,
						true, Read.IGNORE_JUNK, false); cases++;
				failures+=check("dotdash-iupac", simd, SEED, true, false, false, false,
						true, Read.IGNORE_JUNK, false); cases++;
				failures+=check("dotdash-iupac-lower", simd, SEED, true, false, false,
						true, true, Read.IGNORE_JUNK, false); cases++;
				failures+=check("quality-iupac", simd, AMBIGUOUS_SEED, true, false,
						false, false, false, Read.IGNORE_JUNK, true); cases++;
				failures+=check("quality-iupac-upper", simd, AMBIGUOUS_SEED, true,
						false, true, false, false, Read.IGNORE_JUNK, true); cases++;
				failures+=check("mode4-iupac-only", simd, AMBIGUOUS_SEED, false,
						false, false, false, false, Read.FIX_JUNK_AND_IUPAC, false); cases++;
				failures+=check("mode4-junk-iupac", simd, SEED, false, false, false,
						false, false, Read.FIX_JUNK_AND_IUPAC, false); cases++;
				failures+=check("quality-mode4-junk-iupac", simd, SEED, false,
						false, false, false, false, Read.FIX_JUNK_AND_IUPAC, true); cases++;
				failures+=check("fix-junk", simd, JUNK_SEED, false, false, false,
						false, false, Read.FIX_JUNK, false); cases++;
				failures+=check("fix-junk-iupac", simd, JUNK_SEED, true, false,
						false, false, false, Read.FIX_JUNK, false); cases++;
				failures+=check("quality-fix-junk-iupac", simd, JUNK_SEED, true,
						false, false, false, false, Read.FIX_JUNK, true); cases++;
				failures+=checkJunkFlag("flag-junk", simd, JUNK_SEED, false, false,
						false, false, false); cases++;
				failures+=checkJunkFlag("flag-junk-iupac", simd, JUNK_SEED, true,
						false, false, false, false); cases++;
			}
			System.out.println("READ_VALIDATION_NORMALIZATION_RESULT cases="+cases+
					" failures="+failures+" simd_bits="+Vector.SIMD_BITS);
			if(failures!=0){
				throw new AssertionError("Read normalization failures="+failures);
			}
		}finally{
			oldState.restore();
			Shared.AMINO_IN=oldAminoIn;
			Shared.SIMD=oldSimd;
		}
	}

	/** Checks one flag combination.
	 * @param label Case label
	 * @param simd Whether to enable Vector dispatch
	 * @param seed Source pattern
	 * @param iupac Whether Read.IUPAC_TO_N is set
	 * @param uToT Whether Read.U_TO_T is set
	 * @param upper Whether Read.TO_UPPER_CASE is set
	 * @param lower Whether Read.LOWER_CASE_TO_N is set
	 * @param dotDash Whether Read.DOT_DASH_X_TO_N is set
	 * @param junkMode Read.JUNK_MODE value
	 * @param quals Whether to assert quality normalization
	 * @return One on mismatch, otherwise zero
	 */
	private static int check(final String label, final boolean simd, final String seed,
			final boolean iupac, final boolean uToT, final boolean upper,
			final boolean lower, final boolean dotDash, final int junkMode,
			final boolean quals){
		configure(simd, iupac, uToT, upper, lower, dotDash, junkMode, quals);
		final byte[] input=repeat(seed, 96);
		final byte[] expectedBases=expectedBases(input, iupac, uToT, upper,
				lower, dotDash, junkMode);
		final byte[] expectedQuals=(quals ? expectedQuals(expectedBases) : null);
		final byte[] actualBases=input.clone();
		final byte[] actualQuals=(quals ? filledQuals(actualBases.length) : null);
		final Read read=new Read(actualBases, actualQuals, label, 1, false);
		read.validate(true);
		if(Arrays.equals(expectedBases, read.bases) && Arrays.equals(expectedQuals, read.quality)){
			return 0;
		}
		System.out.println("FAIL "+label+" simd="+simd+
				" expectedBases="+new String(expectedBases, StandardCharsets.US_ASCII)+
				" actualBases="+new String(read.bases, StandardCharsets.US_ASCII)+
				" expectedQuals="+Arrays.toString(expectedQuals)+
				" actualQuals="+Arrays.toString(read.quality));
		return 1;
	}

	/** Checks FLAG_JUNK base handling and junk-bit setting.
	 * @param label Case label
	 * @param simd Whether to enable Vector dispatch
	 * @param seed Source pattern
	 * @param iupac Whether Read.IUPAC_TO_N is set
	 * @param uToT Whether Read.U_TO_T is set
	 * @param upper Whether Read.TO_UPPER_CASE is set
	 * @param lower Whether Read.LOWER_CASE_TO_N is set
	 * @param dotDash Whether Read.DOT_DASH_X_TO_N is set
	 * @return One on mismatch, otherwise zero */
	private static int checkJunkFlag(final String label, final boolean simd,
			final String seed, final boolean iupac, final boolean uToT,
			final boolean upper, final boolean lower, final boolean dotDash){
		configure(simd, iupac, uToT, upper, lower, dotDash, Read.FLAG_JUNK, false);
		final byte[] input=repeat(seed, 96);
		final byte[] expectedBases=expectedBases(input, iupac, uToT, upper, lower,
				dotDash, Read.FLAG_JUNK);
		final Read read=new Read(input.clone(), null, label, 1, false);
		read.validate(false);
		if(read.junk() && Arrays.equals(expectedBases, read.bases)){return 0;}
		System.out.println("FAIL "+label+" simd="+simd+" junk="+read.junk()+
				" expectedBases="+new String(expectedBases, StandardCharsets.US_ASCII)+
				" actualBases="+new String(read.bases, StandardCharsets.US_ASCII));
		return 1;
	}

	/** Sets the full static validation policy touched by this test.
	 * @param simd Whether to enable Vector dispatch
	 * @param iupac Whether Read.IUPAC_TO_N is set
	 * @param uToT Whether Read.U_TO_T is set
	 * @param upper Whether Read.TO_UPPER_CASE is set
	 * @param lower Whether Read.LOWER_CASE_TO_N is set
	 * @param dotDash Whether Read.DOT_DASH_X_TO_N is set
	 * @param junkMode Read.JUNK_MODE value
	 * @param quals Whether to enable Read.CHANGE_QUALITY */
	private static void configure(final boolean simd, final boolean iupac,
			final boolean uToT, final boolean upper, final boolean lower,
			final boolean dotDash, final int junkMode, final boolean quals){
		Shared.AMINO_IN=false;
		Shared.SIMD=simd;
		Read.SKIP_SLOW_VALIDATION=false;
		Read.VALIDATE_IN_CONSTRUCTOR=false;
		Read.VALIDATE_VECTOR=true;
		Read.CHANGE_QUALITY=quals;
		Read.TO_UPPER_CASE=upper;
		Read.LOWER_CASE_TO_N=lower;
		Read.DOT_DASH_X_TO_N=dotDash;
		Read.IUPAC_TO_N=iupac;
		Read.U_TO_T=uToT;
		Read.JUNK_MODE=junkMode;
	}

	/** Builds expected bases for Brian's independent normalization flags.
	 * @param input Input bases
	 * @param iupac Whether to mask IUPAC ambiguity letters
	 * @param uToT Whether to convert U/u to T/t
	 * @param upper Whether to uppercase ASCII letters
	 * @param lower Whether to convert lowercase ASCII letters to N
	 * @param dotDash Whether to map dot, dash, X/x and n to N
	 * @param junkMode Read.JUNK_MODE value
	 * @return Expected final bases */
	private static byte[] expectedBases(final byte[] input, final boolean iupac,
			final boolean uToT, final boolean upper, final boolean lower,
			final boolean dotDash, final int junkMode){
		final byte[] out=input.clone();
		for(int i=0; i<out.length; i++){
			byte b=out[i];
			if(uToT){b=uToT(b);}
			if(upper){b=toUpperCase(b);}
			else if(lower){b=lowerCaseToN(b);}
			if(dotDash){b=dotDashToN(b);}
			if(iupac || junkMode==Read.FIX_JUNK_AND_IUPAC){
				b=iupacToN(b);
			}
			if((junkMode==Read.FIX_JUNK || junkMode==Read.FIX_JUNK_AND_IUPAC)
					&& !isIupacNucleotide(b)){
				b='N';
			}
			out[i]=b;
		}
		return out;
	}

	/** Builds expected capped qualities from the final bases.
	 * @param bases Expected final bases
	 * @return Expected quality values */
	private static byte[] expectedQuals(final byte[] bases){
		final byte[] quals=filledQuals(bases.length);
		for(int i=0; i<quals.length; i++){
			if(!isCalledBase(bases[i])){quals[i]=0;}
		}
		return quals;
	}

	/** Repeats a seed to a selected byte length.
	 * @param seed Repeated ASCII pattern
	 * @param length Output length
	 * @return Repeated byte array */
	private static byte[] repeat(final String seed, final int length){
		final byte[] pattern=seed.getBytes(StandardCharsets.US_ASCII);
		final byte[] out=new byte[length];
		for(int i=0; i<out.length; i++){out[i]=pattern[i%pattern.length];}
		return out;
	}

	/** Returns quality values inside Read.qMap's unchanged range.
	 * @param length Output length
	 * @return Filled quality array */
	private static byte[] filledQuals(final int length){
		final byte[] quals=new byte[length];
		Arrays.fill(quals, (byte)40);
		return quals;
	}

	/** Applies Read's U/u-to-T/t policy.
	 * @param b Input base
	 * @return T/t for U/u, otherwise the original byte */
	private static byte uToT(final byte b){
		return b=='U' ? (byte)'T' : b=='u' ? (byte)'t' : b;
	}

	/** Applies ASCII uppercase.
	 * @param b Input byte
	 * @return Uppercase ASCII byte, or the original byte */
	private static byte toUpperCase(final byte b){
		return (b>='a' && b<='z') ? (byte)(b-32) : b;
	}

	/** Applies DOT_DASH_X_TO_N nucleotide mapping.
	 * @param b Input base
	 * @return N for dot, dash, X/x and n; otherwise the original byte */
	private static byte dotDashToN(final byte b){
		return (b=='.' || b=='-' || b=='X' || b=='x' || b=='n') ? (byte)'N' : b;
	}

	/** Applies LOWER_CASE_TO_N nucleotide mapping.
	 * @param b Input byte
	 * @return N for lowercase ASCII, otherwise the original byte */
	private static byte lowerCaseToN(final byte b){
		return (b>='a' && b<='z') ? (byte)'N' : b;
	}

	/** Applies ambiguity-letter-only IUPAC masking.
	 * @param b Input base
	 * @return N for ambiguity letters, otherwise the original byte */
	private static byte iupacToN(final byte b){
		return isAmbiguousIupac(b) ? (byte)'N' : b;
	}

	/** Returns true for non-ACGTU extended IUPAC symbols masked to N.
	 * @param b Input base
	 * @return True for an IUPAC ambiguity letter */
	private static boolean isAmbiguousIupac(final byte b){
		switch(b){
			case 'B':
			case 'D':
			case 'H':
			case 'K':
			case 'M':
			case 'N':
			case 'R':
			case 'S':
			case 'V':
			case 'W':
			case 'X':
			case 'Y':
			case 'b':
			case 'd':
			case 'h':
			case 'k':
			case 'm':
			case 'n':
			case 'r':
			case 's':
			case 'v':
			case 'w':
			case 'x':
			case 'y':
				return true;
			default:
				return false;
		}
	}

	/** Returns true for the accepted nucleotide alphabet before junk repair.
	 * @param b Input base
	 * @return True for A/C/G/T/U or an IUPAC ambiguity letter */
	private static boolean isIupacNucleotide(final byte b){
		return isCalledBase(b) || isAmbiguousIupac(b);
	}

	/** Returns true for bases with retained quality after validation.
	 * @param b Input base
	 * @return True for A/C/G/T/U in either case */
	private static boolean isCalledBase(final byte b){
		switch(b){
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'U':
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'u':
				return true;
			default:
				return false;
		}
	}

	/** Snapshot of Read's static validation policy. */
	private static final class State{
		/** Captures the current static validation policy. */
		State(){}

		/** Captured TO_UPPER_CASE. */
		final boolean toUpperCase=Read.TO_UPPER_CASE;
		/** Captured LOWER_CASE_TO_N. */
		final boolean lowerCaseToN=Read.LOWER_CASE_TO_N;
		/** Captured DOT_DASH_X_TO_N. */
		final boolean dotDashXToN=Read.DOT_DASH_X_TO_N;
		/** Captured VALIDATE_VECTOR. */
		final boolean validateVector=Read.VALIDATE_VECTOR;
		/** Captured SKIP_SLOW_VALIDATION. */
		final boolean skipSlowValidation=Read.SKIP_SLOW_VALIDATION;
		/** Captured IUPAC_TO_N. */
		final boolean iupacToN=Read.IUPAC_TO_N;
		/** Captured U_TO_T. */
		final boolean uToT=Read.U_TO_T;
		/** Captured VALIDATE_IN_CONSTRUCTOR. */
		final boolean validateInConstructor=Read.VALIDATE_IN_CONSTRUCTOR;
		/** Captured CHANGE_QUALITY. */
		final boolean changeQuality=Read.CHANGE_QUALITY;
		/** Captured JUNK_MODE. */
		final int junkMode=Read.JUNK_MODE;

		/** Restores Read's static validation policy. */
		void restore(){
			Read.TO_UPPER_CASE=toUpperCase;
			Read.LOWER_CASE_TO_N=lowerCaseToN;
			Read.DOT_DASH_X_TO_N=dotDashXToN;
			Read.VALIDATE_VECTOR=validateVector;
			Read.SKIP_SLOW_VALIDATION=skipSlowValidation;
			Read.IUPAC_TO_N=iupacToN;
			Read.U_TO_T=uToT;
			Read.VALIDATE_IN_CONSTRUCTOR=validateInConstructor;
			Read.CHANGE_QUALITY=changeQuality;
			Read.JUNK_MODE=junkMode;
		}
	}

	/** Mixed called, ambiguous and junk symbols used by most flag combinations. */
	private static final String SEED="acgtuACGTUryswkmbdhvnRYKMSWBDHVNXx.-!z";
	/** Ambiguity-only input for FIX_JUNK_AND_IUPAC clean-read coverage. */
	private static final String AMBIGUOUS_SEED="acgtuACGTUryswkmbdhvnRYKMSWBDHVNXx";
	/** Ambiguity plus true junk for FIX_JUNK and FLAG_JUNK coverage. */
	private static final String JUNK_SEED="acgtuACGTUryswkmbdhvnRYKMSWBDHVNXx!z";

	/** Prevents construction of this standalone test. */
	private ReadValidationNormalizationTest(){}
}
