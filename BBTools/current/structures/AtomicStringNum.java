package structures;

import java.util.concurrent.atomic.AtomicLong;

public class AtomicStringNum implements Comparable<AtomicStringNum> {

	public AtomicStringNum(String s_, long n_){
		s=s_;
		n=new AtomicLong(n_);
	}

	public long increment(){
		return n.incrementAndGet();
	}
	
	public long increment(long x){
		return n.addAndGet(x);
	}
	
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
	
	public boolean equals(AtomicStringNum other){
		if(other==null){return false;}
		if(n!=other.n){return false;}
		if(s==other.s){return true;}
		if(s==null || other.s==null){return false;}
		return s.equals(other.s);
	}
	
	/*--------------------------------------------------------------*/

	public final String s;
	public AtomicLong n;

}
