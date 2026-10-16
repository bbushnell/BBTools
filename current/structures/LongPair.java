package structures;

/**
 * Mutable pair of public long values.
 * <p>
 * Current callers use {@code a}/{@code b} as paired counters, histogram entries
 * or count/key pairs.  Sorting is ascending by {@code a}, then {@code b},
 * which is inconsistent with inherited identity equality.
 *
 * @author Brian Bushnell
 * @date March 2014
 */
public class LongPair implements Comparable<LongPair>{

	/** Creates a pair with explicit values. */
	public LongPair(long a_, long b_){
		a=a_;
		b=b_;
	}

	/** Creates a zero-initialized pair. */
	public LongPair(){}

	/** Returns the smaller element. */
	public long min(){return Math.min(a, b);}

	/** Returns the larger element. */
	public long max(){return Math.max(a, b);}

	/** Returns {@code a+b}. */
	public long sum(){return a+b;}

	/**
	 * Orders pairs by {@code a}, then {@code b}, both ascending.
	 */
	@Override
	public int compareTo(LongPair other){
		if(a!=other.a){return a>other.a ? 1 : -1;}
		return b>other.b ? 1 : b<other.b ? -1 : 0;
	}

	/** Public mutable values. */
	public long a, b;

}
