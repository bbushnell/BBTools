package structures;

/**
 * Simple container for storing a pair of strings.
 * Used throughout BBTools for paired string data such as key-value associations,
 * paired identifiers, or dual string parameters.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class StringPair {
	
	/**
	 * Creates a new StringPair with the specified string values.
	 * @param a_ First string in the pair
	 * @param b_ Second string in the pair
	 */
	public StringPair(String a_, String b_){
		a=a_;
		b=b_;
	}
	
	@Override
	public String toString(){return "("+a+", "+b+")";}
	
	/** First string in the pair */
	public String a;
	/** Second string in the pair */
	public String b;
	
}
