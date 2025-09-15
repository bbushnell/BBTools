package structures;

import shared.Tools;

public class SeqCountM extends SeqCount {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public SeqCountM(SeqCount sq) {
		super(sq.bases);
		count=sq.count();
	}
	
	public SeqCountM(byte[] s, int start, int stop) {
		super(s, start, stop);
	}
	
	public SeqCountM(byte[] s) {
		super(s);
	}
	
	@Override
	public SeqCountM clone() {
		synchronized(this) {
			SeqCountM clone=(SeqCountM) super.clone();
			return clone;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
//	@Override
	public void add(SeqCount s) {
//		assert(equals(s));
		count+=s.count();
	}

//	@Override
	public void increment(int x) {
		count+=x;
	}

	@Override
	public int count() {return count;}
	
	@Override
	public int compareTo(SeqCount s) {
		if(count()!=s.count()) {return count()-s.count();}
		if(bases.length!=s.bases.length) {return bases.length-s.bases.length;}
		if(s.getClass()==SeqCountM.class) {
			SeqCountM scm=(SeqCountM)s;
			if(score!=scm.score) {return score>scm.score ? 1 : -1;}
		}
		return Tools.compare(bases, s.bases);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int count=1;
	public float score=-1;
	
}