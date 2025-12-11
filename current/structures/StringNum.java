package structures;

import java.io.Serializable;

public class StringNum implements Comparable<StringNum>, Serializable {

	private static final long serialVersionUID=1L;
	
	public StringNum(String s_, long n_){
		s=s_;
		n=n_;
	}

	public long increment(){
		return (n=n+1);
	}
	
	public long increment(long x){
		return (n=n+x);
	}
	
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

	/**
	 * Returns a tab-separated string representation of this StringNum.
	 * Format: "string\tnumber"
	 * @return Tab-delimited string and number
	 */
	@Override
	public String toString(){
		return s+"\t"+n;
	}

	/**
	 * Returns a hash code combining the string hash and lower 31 bits of the number.
	 * Uses XOR combination to distribute hash values effectively.
	 * @return Combined hash code of string and numeric components
	 */
	@Override
	public int hashCode(){
		return ((int)(n&Integer.MAX_VALUE))^(s.hashCode());
	}
	
	/**
	 * Tests equality with another object by casting to StringNum.
	 * @param other The object to compare with
	 * @return true if objects are equal, false otherwise
	 */
	@Override
	public boolean equals(Object other){
		return equals((StringNum)other); //Possible bug: Unchecked cast may throw ClassCastException
	}
	
	public boolean equals(StringNum other){
		if(other==null){return false;}
		if(n!=other.n){return false;}
		if(s==other.s){return true;}
		if(s==null || other.s==null){return false;}
		return s.equals(other.s);
	}
	
	/*--------------------------------------------------------------*/

	public final String s;
	public long n;

}
