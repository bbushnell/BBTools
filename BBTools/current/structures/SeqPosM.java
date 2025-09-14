package structures;

import shared.Tools;

public class SeqPosM implements Cloneable, Comparable<SeqPosM>{

	public SeqPosM(byte[] seq_, int pos_) {this(seq_, pos_, 1);}
	public SeqPosM(SeqPosM sp) {this(sp.seq, sp.pos, sp.count, sp.hashcode, sp.gc);}
	public SeqPosM(byte[] seq_, int pos_, int count_) {
//		synchronized(this) {
			Object o=(seq_==null ? this : seq_);
//			synchronized(o) {
				seq=seq_;
				pos=pos_;
				count=count_;
				hashcode=Tools.hash(seq, 22);
				gc=Tools.calcGC(seq);
				assert(count>=0);
//			}
//		}
	}
	private SeqPosM(byte[] seq_, int pos_, int count_, int code_, float gc_) {
//		synchronized(this) {
//			synchronized(seq_) {
			seq=seq_;
			pos=pos_;
			count=count_;
			hashcode=code_;
			assert(count>=0);
			gc=gc_;
//			}
//		}
	}
	
	public void setFrom(SeqPos sp) {
//		synchronized(sp) {
//			synchronized(sp.seq()) {synchronized(this) {
			seq=sp.seq();
			pos=sp.pos();
			count=sp.count;
			hashcode=sp.hashcode;
			assert(count>=0);
//			}}
//		}
	}
	
	@Override
	public boolean equals(Object o) {
		return equals((SeqPosM)o);
	}
	
	public boolean equals(SeqPosM o) {
		if(pos!=o.pos || hashcode!=o.hashcode) {return false;}
		return Tools.equals(seq, o.seq);
	}
	
	@Override
	public int hashCode() {
		return hashcode;
	}
	
	@Override
	public SeqPosM clone() {
		try {
			return (SeqPosM) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public int compareTo(SeqPosM o) {
		if(count!=o.count) {return o.count-count;}
		if(seq.length!=o.seq.length) {return o.seq.length-seq.length;}
//		if(pos!=o.pos) {return o.pos-pos;}
//		return Tools.compare(seq, o.seq);//Slow; not needed for deterministic behavior if using a stable sort
		return 0;
	}
	
	public final byte[] seq() {return seq;}
	public final int pos() {return pos;}
	public void setPos(int x) {pos=x;}
	
	public byte[] seq;
	public int pos;
	public int hashcode;
	public int count;
	public float gc;
	
}
