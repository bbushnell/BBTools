package structures;

import java.util.concurrent.atomic.AtomicLong;

/**
 * String label paired with an atomically updated {@code long} counter.
 * Increment operations are atomic for the counter value itself; comparisons,
 * equality and hashing read the current counter values separately.
 * Ordering is ascending by live counter value, then by string label.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class AtomicStringNum implements Comparable<AtomicStringNum>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an association with the specified label and initial count.
	 * @param s_ String label.
	 * @param n_ Initial counter value.
	 */
	public AtomicStringNum(String s_, long n_){
		s=s_;
		n=new AtomicLong(n_);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Atomically increments the counter by 1.
	 * @return New counter value.
	 */
	public long increment(){
		return n.incrementAndGet();
	}

	/**
	 * Atomically adds a delta to the counter.
	 * @param x Amount to add; negative values are permitted.
	 * @return New counter value.
	 */
	public long increment(long x){
		return n.addAndGet(x);
	}

	/**
	 * Atomically adds one snapshot of another counter's current value.
	 * @param sn Counter whose current value to add.
	 */
	public void add(AtomicStringNum sn){
		n.addAndGet(sn.n.get());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Ordering           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Compares live counter values, breaking ties by non-null label.
	 * @param o Other value to compare.
	 * @return A negative value, zero or a positive value as this value is less
	 * than, equal to or greater than {@code o}.
	 */
	@Override
	public int compareTo(AtomicStringNum o){
		final long a=n.get(), b=o.n.get();
		if(a<b){return -1;}
		if(a>b){return 1;}
		return s.compareTo(o.s);
	}

	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Formats the label and live counter value separated by a tab.
	 * @return Tab-delimited label and counter.
	 */
	@Override
	public String toString(){
		return s+"\t"+n;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Equality           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Hashes the non-null label and live counter value.
	 * @return Hash code derived from the current state.
	 */
	@Override
	public int hashCode(){
		return ((int)(n.get()&Integer.MAX_VALUE))^(s.hashCode());
	}

	/**
	 * Tests equality with another object.
	 * @param other Object to compare.
	 * @return True if {@code other} is an AtomicStringNum with the same label
	 * and live counter value.
	 */
	@Override
	public boolean equals(Object other){
		return other instanceof AtomicStringNum && equals((AtomicStringNum)other);
	}

	/**
	 * Tests equality with another AtomicStringNum.
	 * @param other AtomicStringNum to compare; may be null.
	 * @return True if both objects have the same label and live counter value.
	 */
	public boolean equals(AtomicStringNum other){
		if(other==null){return false;}
		if(n.get()!=other.n.get()){return false;}
		if(s==other.s){return true;}
		if(s==null || other.s==null){return false;}
		return s.equals(other.s);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** String label associated with this counter. */
	public final String s;

	/** Publicly replaceable atomic counter. */
	public AtomicLong n;

}
