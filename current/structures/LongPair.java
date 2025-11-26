package structures;

/**
 * Simple pair structure for holding two long values with comparison support.
 * Provides basic arithmetic operations and lexicographic ordering.
 * Commonly used for coordinate pairs, ranges, and key-value associations in genomic data processing.
 *
 * @author Brian Bushnell
 * @date March 2014
 */
public class LongPair implements Comparable<LongPair>{

	/**
	 * Creates a LongPair with specified values.
	 * @param a_ First long value
	 * @param b_ Second long value
	 */
	public LongPair(long a_, long b_){
		a=a_;
		b=b_;
	}

	/** Creates a LongPair with default values (0, 0). */
	public LongPair(){}

	/** Returns the smaller of the two values.
	 * @return Minimum of a and b */
	public long min() {return Math.min(a, b);}
	/** Returns the larger of the two values.
	 * @return Maximum of a and b */
	public long max() {return Math.max(a, b);}
	/** Returns the sum of both values.
	 * @return a + b */
	public long sum() {return a+b;}
	
	@Override
	public int compareTo(LongPair other) {
		if(a!=other.a){return a>other.a ? 1 : -1;}
		return b>other.b ? 1 : b<other.b ? -1 : 0;
	}
	
	/** Second long value stored in this pair. */
	/** First long value stored in this pair. */
	public long a, b;
	
}
