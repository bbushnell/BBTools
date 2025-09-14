package structures;

public class StringNum implements Comparable<StringNum> {

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
		return equals((StringNum)other);
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
