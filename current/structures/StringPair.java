package structures;

/**
 * Mutable pair of public strings.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class StringPair{

	/** Stores two string references without copying them.
	 * @param a_ First value
	 * @param b_ Second value */
	public StringPair(String a_, String b_){
		a=a_;
		b=b_;
	}

	@Override
	public String toString(){return "("+a+", "+b+")";}

	/** First value. */
	public String a;
	/** Second value. */
	public String b;

}
