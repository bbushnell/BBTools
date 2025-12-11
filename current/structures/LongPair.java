package structures;

/**
 * Pair of long values with comparison support.
 * Provides min/max/sum helpers, lexicographic ordering, and is commonly used for coordinate pairs, ranges, or key/value associations in genomic pipelines.
 * @author Brian Bushnell
 * @date March 2014
 */
public class LongPair implements Comparable<LongPair>{

	/**
	 * Creates a LongPair with the specified values.
	 * @param a_ First long value
	 * @param b_ Second long value
	 */
	public LongPair(long a_, long b_){
		a=a_;
		b=b_;
	}

	/** Creates a LongPair with default values of 0 for both elements. */
	public LongPair(){}

	public long min() {return Math.min(a, b);}
	public long max() {return Math.max(a, b);}
	public long sum() {return a+b;}
	
	/**
	 * Compares this pair to another using lexicographic ordering.
	 * Compares 'a' values first, then 'b' values if 'a' values are equal.
	 * @param other The LongPair to compare against
	 * @return Negative if this < other, positive if this > other, 0 if equal
	 */
	@Override
	public int compareTo(LongPair other) {
		if(a!=other.a){return a>other.a ? 1 : -1;}
		return b>other.b ? 1 : b<other.b ? -1 : 0;
	}
	
	/** Second long value stored in this pair. */
	public long a, b;
	
}
