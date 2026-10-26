package map;

import shared.Tools;

/** Mutable sequence placement used as per-thread SeqMap lookup scratch.
 * Equality uses sequence content and position; ranking uses count and length.
 * Sequence storage is borrowed. Do not mutate a live hash key; SeqMap clones
 * its scratch placement before retaining it in a result set.
 * @author Brian Bushnell
 */
public class SeqPosM implements Cloneable, Comparable<SeqPosM>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Computes metadata for a sequence with occurrence count one. */
	public SeqPosM(byte[] seq_, int pos_){this(seq_, pos_, 1);}

	/** Copies metadata and borrows the same sequence array. */
	public SeqPosM(SeqPosM sp){this(sp.seq, sp.pos, sp.count, sp.hashcode, sp.gc);}

	/** Computes sequence hash and GC once; null is permitted for empty scratch.
	 * @param seq_ Borrowed sequence bytes, or null for uninitialized scratch
	 * @param pos_ Placement coordinate
	 * @param count_ Nonnegative occurrence count
	 */
	public SeqPosM(byte[] seq_, int pos_, int count_){
		seq=seq_;
		pos=pos_;
		count=count_;
		hashcode=Tools.hash(seq, 22);
		gc=Tools.calcGC(seq);
		assert(count>=0) : "Sequence occurrence counts must be nonnegative for descending count ranking: "+count;
	}

	/** Copies already-computed metadata without rescanning the borrowed sequence. */
	private SeqPosM(byte[] seq_, int pos_, int count_, int code_, float gc_){
		seq=seq_;
		pos=pos_;
		count=count_;
		hashcode=code_;
		assert(count>=0) : "Copied sequence occurrence count must remain nonnegative: "+count;
		gc=gc_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces all corresponding metadata and borrows the source sequence. */
	public void setFrom(SeqPos sp){
		seq=sp.seq();
		pos=sp.pos();
		count=sp.count;
		hashcode=sp.hashcode;
		gc=sp.gc;
		assert(count>=0) : "SeqMap candidate metadata requires a nonnegative occurrence count: "+count;
	}

	/** Compares a homogeneous, non-null SeqPosM key. Foreign/null keys are unsupported. */
	@Override
	public boolean equals(Object o){return equals((SeqPosM)o);}

	/** Tests position and sequence identity, independently of count and GC. */
	public boolean equals(SeqPosM o){
		if(pos!=o.pos || hashcode!=o.hashcode){return false;}
		return Tools.equals(seq, o.seq);
	}

	/** Returns the cached sequence hash; different positions may share a hash. */
	@Override
	public int hashCode(){return hashcode;}

	/** Returns a shallow metadata copy sharing the sequence array. */
	@Override
	public SeqPosM clone(){
		try{
			return (SeqPosM)super.clone();
		}catch(CloneNotSupportedException e){
			throw new RuntimeException(e);
		}
	}

	/** Ranks nonnegative counts and sequence lengths descending; stable sorting
	 * preserves input order for rank ties. Ranking is distinct from key identity.
	 */
	@Override
	public int compareTo(SeqPosM o){
		if(count!=o.count){return o.count-count;}
		if(seq.length!=o.seq.length){return o.seq.length-seq.length;}
		return 0;
	}

	/** Returns the borrowed sequence. */
	public final byte[] seq(){return seq;}
	/** Returns the placement coordinate. */
	public final int pos(){return pos;}
	/** Changes position; the object must not be a live hash key. */
	public final void setPos(int x){pos=x;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed sequence bytes; null is allowed only for uninitialized scratch. */
	public byte[] seq;
	/** Sequence offset relative to the first query region in SeqMap results. */
	public int pos;
	/** Cached sequence-content hash. */
	public int hashcode;
	/** Nonnegative occurrence count. */
	public int count;
	/** Cached GC fraction, potentially NaN when no bases are defined. */
	public float gc;
}
