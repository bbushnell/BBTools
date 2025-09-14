package structures;

public class LongPair implements Comparable<LongPair>{

	public LongPair(long a_, long b_){
		a=a_;
		b=b_;
	}

	public LongPair(){}

	public long min() {return Math.min(a, b);}
	public long max() {return Math.max(a, b);}
	public long sum() {return a+b;}
	
	@Override
	public int compareTo(LongPair other) {
		if(a!=other.a){return a>other.a ? 1 : -1;}
		return b>other.b ? 1 : b<other.b ? -1 : 0;
	}
	
	public long a, b;
	
}
