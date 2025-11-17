package icecream;

/**
 * Trims byte arrays by detecting homogeneous symbol sequences with dynamic scoring.
 * Implements left and right side trimming of byte arrays by calculating a dynamic
 * score for symbol matching, allowing flexible polymer sequence removal.
 * Used in bioinformatics for sequence preprocessing and adapter removal.
 *
 * @author Brian Bushnell
 */
public class PolymerTrimmer {
	
	/**
	 * Parses command-line arguments for polymer trimming configuration.
	 * Supports minPolymer, minFraction, and polyerror parameters.
	 *
	 * @param arg The full argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the argument was recognized and parsed, false otherwise
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
	 * Tests for polymer trimming from the left side using char symbol.
	 * Convenience method that converts char to byte and calls main testLeft method.
	 *
	 * @param bases The sequence bases to analyze
	 * @param symbol The polymer symbol to detect as char
	 * @return Number of bases to trim from left end, or 0 if below minimum
	 */
	public static int testLeft(byte[] bases, char symbol){return testLeft(bases, (byte)symbol);}
	
	/**
	 * Tests for polymer trimming from the left side using dynamic scoring.
	 * Scans left to right, incrementing score for matching symbols and applying
	 * penalty for mismatches. Tracks maximum score position for optimal trim point.
	 *
	 * @param bases The sequence bases to analyze
	 * @param symbol The polymer symbol to detect
	 * @return Number of bases to trim from left end, or 0 if below minPolymer
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
	 * Tests for polymer trimming from the right side using char symbol.
	 * Convenience method that converts char to byte and calls main testRight method.
	 *
	 * @param bases The sequence bases to analyze
	 * @param symbol The polymer symbol to detect as char
	 * @return Number of bases to trim from right end, or 0 if below minimum
	 */
	public static int testRight(byte[] bases, char symbol){return testRight(bases, (byte)symbol);}
	
	/**
	 * Tests for polymer trimming from the right side using dynamic scoring.
	 * Scans right to left, incrementing score for matching symbols and applying
	 * penalty for mismatches. Tracks maximum score position for optimal trim point.
	 *
	 * @param bases The sequence bases to analyze
	 * @param symbol The polymer symbol to detect
	 * @return Number of bases to trim from right end, or 0 if below minPolymer
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
	 * Sets the minimum fraction of matching symbols required for polymer detection.
	 * Recalculates penalty and minimum score thresholds based on the fraction.
	 * Higher fractions require more homogeneous sequences for trimming.
	 * @param f Minimum fraction (0.0 to 1.0) of symbols that must match
	 */
	public static void setMinFraction(float f){
		assert(f>=0 && f<=1) : f;
		minFraction=f;
		penalty=(f>=1 ? 99 : ((1f/(1-minFraction))-1));
		minScore=(f>=1 ? 0 : -4*penalty);
	}
	
	/** Minimum length of polymer sequence required for trimming */
	static int minPolymer=5;
	/** Minimum fraction of matching symbols required for polymer detection */
	private static float minFraction=0.8f;
	/** Penalty score applied for non-matching symbols during dynamic scoring */
	private static float penalty=(1f/(1-minFraction))-1;
	/** Minimum score threshold for continuing polymer detection scan */
	private static float minScore=-4*penalty;
	
}
