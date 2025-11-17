package structures;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Thread-safe association of a string identifier with an atomic long counter.
 * Provides atomic increment operations and comparison capabilities for concurrent
 * counting applications such as k-mer frequency tracking and statistical analysis.
 * Comparison is first by numeric value, then by string lexicographically.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class AtomicStringNum implements Comparable<AtomicStringNum> {

	/**
	 * Creates a new AtomicStringNum with the specified string and initial count.
	 * @param s_ The string identifier
	 * @param n_ The initial numeric value
	 */
	public AtomicStringNum(String s_, long n_){
		s=s_;
		n=new AtomicLong(n_);
	}

	/** Atomically increments the numeric value by 1 and returns the new value.
	 * @return The incremented value */
	public long increment(){
		return n.incrementAndGet();
	}
	
	/**
	 * Atomically adds the specified value to the numeric counter.
	 * @param x The value to add
	 * @return The new value after addition
	 */
	public long increment(long x){
		return n.addAndGet(x);
	}
	
	/** Atomically adds the numeric value from another AtomicStringNum to this one.
	 * @param sn The AtomicStringNum whose value to add */
	public void add(AtomicStringNum sn) {
		n.addAndGet(sn.n.get());
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(AtomicStringNum o) {
		final long a=n.get(), b=o.n.get();
		if(a<b){return -1;}
		if(a>b){return 1;}
		return s.compareTo(o.s);
	}

	@Override
	public String toString(){
		return s+"\t"+n;
	}

	@Override
	public int hashCode(){
		return ((int)(n.get()&Integer.MAX_VALUE))^(s.hashCode());
	}
	
	@Override
	public boolean equals(Object other){
		return equals((AtomicStringNum)other);
	}
	
	/**
	 * Tests equality with another AtomicStringNum.
	 * Equality requires both identical numeric values and string content.
	 * @param other The AtomicStringNum to compare
	 * @return true if both string and numeric components are equal
	 */
	public boolean equals(AtomicStringNum other){
		if(other==null){return false;}
		if(n!=other.n){return false;}
		if(s==other.s){return true;}
		if(s==null || other.s==null){return false;}
		return s.equals(other.s);
	}
	
	/*--------------------------------------------------------------*/

	/** The string identifier associated with this counter */
	public final String s;
	/** The atomic long counter providing thread-safe numeric operations */
	public AtomicLong n;

}
