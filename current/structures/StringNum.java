package structures;

import java.io.Serializable;

/**
 * Serializable pair of a String label and mutable long count.
 * <p>
 * Ordering is by count, then label; equality and hashing use both fields.
 * Mutation is unsynchronized; callers sharing an instance coordinate externally.
 * Because the count participates in equality, hashing and ordering, callers must
 * not mutate it while an instance is in a hash-based or sorted collection.
 *
 * @author Brian Bushnell
 */
public class StringNum implements Comparable<StringNum>, Serializable{

	private static final long serialVersionUID=1L;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a label/count pair.
	 * @param s_ Label; callers that sort or hash instances must supply non-null strings.
	 * @param n_ Initial count.
	 */
	public StringNum(String s_, long n_){
		s=s_;
		n=n_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/** Increments the count by 1 and returns the updated value. */
	public long increment(){
		return (n=n+1);
	}

	/** Adds x to the count and returns the updated value. */
	public long increment(long x){
		return (n=n+x);
	}

	/** Adds sn's count to this count. */
	public void add(StringNum sn){
		n+=sn.n;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Ordering           ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares counts in ascending order, then labels lexicographically. */
	@Override
	public int compareTo(StringNum o){
		if(n<o.n){return -1;}
		if(n>o.n){return 1;}
		return s.compareTo(o.s);
	}

	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns {@code s} and {@code n} separated by one tab. */
	@Override
	public String toString(){
		return s+"\t"+n;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Equality           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int hashCode(){
		return ((int)(n&Integer.MAX_VALUE))^(s.hashCode());
	}

	@Override
	public boolean equals(Object other){
		return other instanceof StringNum && equals((StringNum)other);
	}

	/** Tests for the same count and label, including null-safe label comparison. */
	public boolean equals(StringNum other){
		if(other==null){return false;}
		if(n!=other.n){return false;}
		if(s==other.s){return true;}
		if(s==null || other.s==null){return false;}
		return s.equals(other.s);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Label used for identity and as the secondary sort key. */
	public final String s;

	/** Mutable count. */
	public long n;

}
