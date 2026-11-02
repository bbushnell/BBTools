package test.repeat;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import repeat.Palindrome;
import repeat.PalindromeFinder;

/** Exhaustive finite oracle for short uppercase nucleotide palindrome searches.
 * Enumerates intervals and arm lengths independently of the detector's center
 * expansion. Tests every A/C/G/T sequence of length 2 through 6, both as a whole
 * array and inside a padded subregion. Perfect and one-mismatch configurations
 * are reported separately. This does not certify long inputs, arbitrary tails,
 * undefined/lowercase symbols or concurrent finder use.
 * @author Sucrose
 */
public final class PalindromeOracleTest{

	/** Accepts mode=perfect|imperfect|both; default both. No files or worker threads are used. */
	public static void main(String[] args){
		String mode="both";
		if(args.length>1 || (args.length==1 && !args[0].startsWith("mode="))){
			throw new IllegalArgumentException("Expected mode=perfect|imperfect|both");
		}
		if(args.length==1){mode=args[0].substring(5);}
		require(mode.equals("perfect") || mode.equals("imperfect") || mode.equals("both"), "Unknown mode: "+mode);
		long cases=0, failures=0;
		final int first=(mode.equals("imperfect") ? 1 : 0), last=(mode.equals("perfect") ? 0 : 1);
		for(int mismatches=first; mismatches<=last; mismatches++){
			for(int arm=1; arm<=2; arm++){
				for(int[] bounds : LOOPS){
					final PalindromeOracleTest test=new PalindromeOracleTest(arm, bounds[0], bounds[1], mismatches);
					test.run();
					cases+=test.cases;
					failures+=test.failures;
					System.out.println("arm="+arm+" minloop="+bounds[0]+" maxloop="+bounds[1]
							+" maxmismatches="+mismatches+" cases="+test.cases+" failures="+test.failures);
				}
			}
		}
		final long expected=10912L*2*LOOPS.length*(last-first+1);
		require(cases==expected, "Finite enumeration count changed: "+cases+" expected="+expected);
		System.out.println("PalindromeOracleTest cases="+cases+" failures="+failures);
		if(failures>0){throw new AssertionError("Finite palindrome oracle mismatches="+failures);}
	}

	/** Uses unrestricted tails, at least one match, and a fixed minimum arm length. */
	private PalindromeOracleTest(int arm, int minLoop_, int maxLoop_, int mismatches){
		minArm=arm;
		minLoop=minLoop_;
		maxLoop=maxLoop_;
		maxMismatches=mismatches;
		finder=new PalindromeFinder(arm, minLoop, maxLoop, 1, mismatches, 0,
				Integer.MAX_VALUE, Integer.MAX_VALUE);
	}

	/** Enumerates 5456 sequences, each with two equivalent coordinate presentations. */
	private void run(){
		for(int length=2; length<=6; length++){
			final byte[] bases=new byte[length], padded=new byte[length+4];
			Arrays.fill(padded, (byte)'G');
			for(int code=0, limit=1<<(2*length); code<limit; code++){
				int value=code;
				for(int i=0; i<length; i++, value>>>=2){bases[i]=ALPHABET[value&3];}
				check(bases, 0, length-1);
				System.arraycopy(bases, 0, padded, 2, length);
				check(padded, 2, length+1);
			}
		}
	}

	/** Compares all four returned scalars; null and thrown assertions remain distinct failures. */
	private void check(byte[] bases, int from, int to){
		cases++;
		final Palindrome expected=oracle(bases, from, to);
		Palindrome actual=null;
		AssertionError thrown=null;
		try{actual=finder.longestPalindrome(bases, from, to);}
		catch(AssertionError failure){thrown=failure;}
		final boolean equal=(thrown==null && (expected==null ? actual==null : actual!=null
				&& expected.a==actual.a && expected.b==actual.b && expected.matches==actual.matches
				&& expected.mismatches==actual.mismatches));
		if(!equal){
			failures++;
			if(failures<=4){
				System.out.println("MISMATCH sequence="+new String(bases, StandardCharsets.US_ASCII)
						+" region="+from+".."+to+" expected="+expected+" actual="+actual+" assertion="+thrown);
			}
		}
	}

	/** Enumerates every inclusive interval and possible arm length, then applies the contract.
	 * Both innermost and outermost paired positions must match, as required by the
	 * detector's imperfect-stem contract. Interior pairs may consume the mismatch budget. */
	private Palindrome oracle(byte[] bases, int from, int to){
		Palindrome best=null;
		for(int a=from; a<to; a++){
			for(int b=a+1; b<=to; b++){
				for(int arm=minArm; 2*arm<=b-a+1; arm++){
					final int loop=b-a+1-2*arm;
					if(loop<minLoop || loop>maxLoop){continue;}
					if(!complement(bases[a], bases[b]) || !complement(bases[a+arm-1], bases[b-arm+1])){continue;}
					int bad=0;
					for(int offset=0; offset<arm; offset++){
						if(!complement(bases[a+offset], bases[b-offset])){bad++;}
					}
					final int good=arm-bad;
					if(bad>maxMismatches || good<1){continue;}
					if(better(a, b, good, bad, best)){best=new Palindrome(a, b, good, bad);}
				}
			}
		}
		return best;
	}

	/** Independent explicit nucleotide complement relation, with no production lookup table. */
	private static boolean complement(byte a, byte b){
		return (a=='A' && b=='T') || (a=='T' && b=='A') || (a=='C' && b=='G') || (a=='G' && b=='C');
	}

	/** Implements the documented comparison priority without calling production compareTo. */
	private static boolean better(int a, int b, int matches, int mismatches, Palindrome best){
		if(best==null){return true;}
		if(matches!=best.matches){return matches>best.matches;}
		if(mismatches!=best.mismatches){return mismatches<best.mismatches;}
		if(b-a!=best.b-best.a){return b-a>best.b-best.a;}
		return a<best.a;
	}

	/** Checks the corpus and argument contracts under both JVM assertion modes. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** One reused, single-caller production finder and its independently applied constraints. */
	private final PalindromeFinder finder;
	private final int minArm, minLoop, maxLoop, maxMismatches;
	/** Per-configuration comparison and failure totals. */
	private long cases, failures;
	/** Literal DNA alphabet and exact/ranged loop constraints; no random sampling. */
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
	private static final int[][] LOOPS={{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 3}, {2, 2}, {2, 4}, {3, 3}};
}
