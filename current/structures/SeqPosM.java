package structures;

import shared.Tools;

/**
 * Represents a sequence position with metadata including count, GC content, and hash code.
 * This mutable version of SeqPos provides comparison capabilities based on count and length,
 * with optimized hashing and GC content calculation for sequence analysis.
 * @author Brian Bushnell
 */
public class SeqPosM implements Cloneable, Comparable<SeqPosM>{

	/**
	 * Creates a SeqPosM with count defaulted to 1.
	 * @param seq_ The sequence bytes
	 * @param pos_ Position within the sequence
	 */
	public SeqPosM(byte[] seq_, int pos_) {this(seq_, pos_, 1);}
	/**
	 * Copy constructor that creates a new SeqPosM from an existing one.
	 * Copies all fields including pre-computed hash code and GC content.
	 * @param sp The SeqPosM to copy
	 */
	public SeqPosM(SeqPosM sp) {this(sp.seq, sp.pos, sp.count, sp.hashcode, sp.gc);}
	/**
	 * Creates a SeqPosM with specified count and computes hash code and GC content.
	 * Hash code is computed using Tools.hash with seed 22 for consistent hashing.
	 *
	 * @param seq_ The sequence bytes
	 * @param pos_ Position within the sequence
	 * @param count_ Count value for this sequence position
	 */
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
	/**
	 * Private constructor that sets all fields directly without computation.
	 * Used internally when hash code and GC content are already known.
	 *
	 * @param seq_ The sequence bytes
	 * @param pos_ Position within the sequence
	 * @param count_ Count value
	 * @param code_ Pre-computed hash code
	 * @param gc_ Pre-computed GC content
	 */
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
	
	/**
	 * Updates this SeqPosM with values from a SeqPos object.
	 * Copies sequence, position, count, and hash code from the source.
	 * @param sp The SeqPos to copy values from
	 */
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
	
	/**
	 * Compares this SeqPosM with another for equality.
	 * First checks position and hash code for quick rejection,
	 * then performs full sequence comparison if needed.
	 *
	 * @param o The SeqPosM to compare with
	 * @return true if sequences and positions are equal
	 */
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
	
	/** Returns the sequence bytes */
	public final byte[] seq() {return seq;}
	/** Returns the position within the sequence */
	public final int pos() {return pos;}
	/** Sets the position within the sequence.
	 * @param x The new position value */
	public void setPos(int x) {pos=x;}
	
	/** The sequence bytes */
	public byte[] seq;
	/** Position within the sequence */
	public int pos;
	/** Pre-computed hash code for the sequence using Tools.hash with seed 22 */
	public int hashcode;
	/** Count or frequency associated with this sequence position */
	public int count;
	/** GC content of the sequence calculated by Tools.calcGC */
	public float gc;
	
}
