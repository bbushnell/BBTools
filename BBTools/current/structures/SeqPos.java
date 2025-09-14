package structures;

import shared.Tools;

public class SeqPos implements Cloneable, Comparable<SeqPos>{

//	public SeqPos(byte[] seq_, int pos_) {this(seq_, pos_, 1);}
	public SeqPos(SeqPos sp) {this(sp.seq, sp.pos, sp.count, sp.hashcode, sp.gc, sp.score);}
//	public SeqPos(byte[] seq_, int pos_, int count_) {
////		synchronized(this) {
//			seq=seq_;
//			pos=pos_;
//			count=count_;
//			hashcode=Tools.hash(seq, 22);
//			gc=Tools.calcGC(seq);
//			assert(count>=0);
////		}
//	}
	SeqPos(byte[] seq_, int pos_, int count_, int code_, float gc_, float score_) {
//		synchronized(this) {
			seq=seq_;
			pos=pos_;
			count=count_;
			hashcode=code_;
			gc=gc_;
			score=score_;
			assert(count>=0);
//		}
	}
	
//	public void setFrom(SeqPos sp) {
//		synchronized(sp) {
//			synchronized(sp.seq) {synchronized(this) {
//			seq=sp.seq;
//			pos=sp.pos;
////			count=sp.count;
//			hashcode=sp.hashcode;
//			assert(count>0);
//			}}
//		}
//	}
	
	@Override
	public boolean equals(Object o) {
		return equals((SeqPos)o);
	}
	
	public boolean equals(SeqPos o) {
		if(pos!=o.pos || hashcode!=o.hashcode) {return false;}
		return Tools.equals(seq, o.seq);
	}
	
	@Override
	public int hashCode() {
		return hashcode;
	}
	
	@Override
	public SeqPos clone() {
		try {
			return (SeqPos) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public int compareTo(SeqPos o) {
		if(count!=o.count) {return o.count-count;}
		if(seq.length!=o.seq.length) {return o.seq.length-seq.length;}
		if(score!=o.score) {return o.score>score ? 1 : -1;}
//		if(gc!=o.gc) {
//			return Tools.absdif((int)(gc*100),  50)<Tools.absdif((int)(o.gc*100), 50) ? -1 : 1;
//		}//Neutral gc first
//		return Tools.compare(seq, o.seq);
		return 0;
	}
	
	public final byte[] seq() {return seq;}
	public final int pos() {return pos;}
	
	public final byte[] seq;
	public final int pos;
	public final int hashcode;
	public final int count;
	public final float gc;
	public final float score;
	
}
