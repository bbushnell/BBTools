package structures;

import java.io.Serializable;

public class StringNum implements Comparable<StringNum>, Serializable {

	private static final long serialVersionUID=1L;
	
	/**
	 * Constructs a StringNum with the specified string and numeric value.
	 * @param s_ The string component
	 * @param n_ The numeric value component
	 */
	public StringNum(String s_, long n_){
		s=s_;
		n=n_;
	}

	/** Increments the numeric value by 1 and returns the new value.
	 * @return The incremented numeric value */
	public long increment(){
		return (n=n+1);
	}
	
	/**
	 * Increments the numeric value by the specified amount and returns the new value.
	 * @param x The amount to add to the numeric value
	 * @return The incremented numeric value
	 */
	public long increment(long x){
		return (n=n+x);
	}
	
	/**
	 * Adds another StringNum's numeric value to this object's numeric value.
	 * The string component is unchanged.
	 * @param sn The StringNum whose numeric value to add
	 */
	public void add(StringNum sn) {
		n+=sn.n;
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(StringNum o) {
		if(n<o.n){return -1;}
		if(n>o.n){return 1;}
		return s.compareTo(o.s);
	}

	@Override
	public String toString(){
		return s+"\t"+n;
	}

	@Override
	public int hashCode(){
		return ((int)(n&Integer.MAX_VALUE))^(s.hashCode());
	}
	
	@Override
	public boolean equals(Object other){
		return equals((StringNum)other); //Possible bug: Unchecked cast may throw ClassCastException
	}
	
	/**
	 * Tests equality with another StringNum based on both string and numeric values.
	 * Returns true only if both components match exactly.
	 * @param other The StringNum to compare with
	 * @return true if both string and numeric values are equal, false otherwise
	 */
	public boolean equals(StringNum other){
		if(other==null){return false;}
		if(n!=other.n){return false;}
		if(s==other.s){return true;}
		if(s==null || other.s==null){return false;}
		return s.equals(other.s);
	}
	
	/*--------------------------------------------------------------*/

	/** The string component of this StringNum */
	public final String s;
	/** The numeric value component that can be modified via increment operations */
	public long n;

}
