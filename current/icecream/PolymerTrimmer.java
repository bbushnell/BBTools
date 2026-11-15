package icecream;

/**
 * Finds candidate homopolymer end spans using a running match/mismatch score.
 * Methods return lengths without mutating input. Comparisons are exact byte
 * comparisons, including case. The fraction setting controls score penalties and
 * scan tolerance; no separate measured-fraction check is applied to the result.
 * Configuration is static and unsynchronized: finish configuration before workers
 * trim concurrently, since a setter changes later calls throughout the JVM.
 * @author Brian Bushnell, Sucrose
 */
public class PolymerTrimmer{
	
	/**
	 * Updates shared configuration for minPolymer, minFraction, or polyerror.
	 * Parameter names are case-insensitive; polyerror is converted to 1-value
	 * before calling setMinFraction. minPolymer receives the parsed integer directly.
	 * @param arg Original argument, unused by this parser
	 * @param a Non-null parameter name
	 * @param b Value to parse for a recognized parameter
	 * @return Whether the name was recognized; malformed recognized values can fail parsing
	 */
	public static boolean parse(String arg, String a, String b){
		if(a.equalsIgnoreCase("minPolymer")){
			minPolymer=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minFraction")){
			float f=Float.parseFloat(b);
			setMinFraction(f);
		}else if(a.equalsIgnoreCase("polyerror")){
			float f=Float.parseFloat(b);
			setMinFraction(1-f);
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Delegates to the byte-symbol prefix scan after casting the character to byte.
	 * @param bases Non-null input, read without mutation
	 * @param symbol Symbol cast to the compared byte value
	 * @return Candidate prefix length, or zero
	 */
	public static int testLeft(byte[] bases, char symbol){return testLeft(bases, (byte)symbol);}
	
	/**
	 * Scans from the left while the running score is at least minScore. A match
	 * adds one and a mismatch subtracts penalty. The earliest position attaining the
	 * maximum positive score determines the prefix, including intervening mismatches.
	 * @param bases Non-null input, read without mutation; empty arrays return zero
	 * @param symbol Byte value compared exactly, without case conversion
	 * @return Prefix span if at least minPolymer positions long, otherwise zero
	 */
	public static int testLeft(byte[] bases, byte symbol){
		float score=0;
		float max=0;
		int maxPos=-1;
		for(int i=0; i<bases.length && score>=minScore; i++){
			byte b=bases[i];
			if(b==symbol){
				score++;
				if(score>max){
					max=score;
					maxPos=i;
				}
			}else{
				score-=penalty;
			}
		}
		int trim=maxPos+1;
		return (trim<minPolymer ? 0 : trim);
	}
	
	/**
	 * Delegates to the byte-symbol suffix scan after casting the character to byte.
	 * @param bases Non-null input, read without mutation
	 * @param symbol Symbol cast to the compared byte value
	 * @return Candidate suffix length, or zero
	 */
	public static int testRight(byte[] bases, char symbol){return testRight(bases, (byte)symbol);}
	
	/**
	 * Mirrors the prefix scan from the right. The first occurrence of the maximum
	 * positive score in that direction determines the suffix, including mismatches.
	 * @param bases Non-null input, read without mutation; empty arrays return zero
	 * @param symbol Byte value compared exactly, without case conversion
	 * @return Suffix span if at least minPolymer positions long, otherwise zero
	 */
	public static int testRight(byte[] bases, byte symbol){
		float score=0;
		float max=0;
		int maxPos=bases.length;
		for(int i=bases.length-1; i>=0 && score>=minScore; i--){
			byte b=bases[i];
			if(b==symbol){
				score++;
				if(score>max){
					max=score;
					maxPos=i;
				}
			}else{
				score-=penalty;
			}
		}
		int trim=bases.length-maxPos;
		return (trim<minPolymer ? 0 : trim);
	}
	
	/**
	 * Recomputes global penalty and scan cutoff from the requested fraction.
	 * Below one, penalty is 1/(1-f)-1 and minScore is -4*penalty. At one, the
	 * historical finite penalty99 and cutoff0 are retained. This configures a score
	 * heuristic, not a post-scan fraction calculation. Call before concurrent trimming.
	 * @param f Fraction setting, checked by assertion to lie in [0,1]
	 * @throws AssertionError If assertions are enabled and f is outside [0,1] or NaN
	 */
	public static void setMinFraction(float f){
		assert(f>=0 && f<=1) : f;
		minFraction=f;
		penalty=(f>=1 ? 99 : ((1f/(1-minFraction))-1));
		minScore=(f>=1 ? 0 : -4*penalty);
	}
	
	/** Minimum candidate span, counting all positions rather than only matching symbols. */
	static int minPolymer=5;
	/** Requested fraction setting used to derive the scoring parameters. */
	private static float minFraction=0.8f;
	/** Score subtracted for each nonmatching byte. */
	private static float penalty=(1f/(1-minFraction))-1;
	/** Running-score cutoff checked before examining the next byte. */
	private static float minScore=-4*penalty;
	
}
