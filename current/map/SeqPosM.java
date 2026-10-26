package map;

import shared.Tools;

/** Mutable placement used by SeqMap as per-thread lookup scratch and result entries.
 * Equality uses sequence content and position; ranking uses count and length.
 * Sequence storage is borrowed. Do not mutate a live hash key; SeqMap clones
 * its scratch placement before retaining it in a result set.
 * @author Brian Bushnell
 */
public class SeqPosM implements Cloneable, Comparable<SeqPosM>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Computes metadata for a sequence with occurrence count one.
	 * @param seq_ Borrowed sequence, or null for uninitialized scratch
	 * @param pos_ Placement coordinate */
	public SeqPosM(byte[] seq_, int pos_){this(seq_, pos_, 1);}

	/** Copies metadata and borrows the same sequence array.
	 * @param sp Non-null scratch placement to copy */
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

	/** Copies already-computed metadata without rescanning the borrowed sequence.
	 * @param seq_ Borrowed sequence, possibly null for uninitialized scratch
	 * @param pos_ Placement coordinate
	 * @param count_ Nonnegative occurrence count
	 * @param code_ Cached sequence hash consistent with seq_
	 * @param gc_ Cached GC fraction, possibly NaN for no defined bases */
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

	/** Replaces corresponding metadata and borrows the source sequence; its score is not stored.
	 * This mutable object must not be a live hash key while its identity is changed.
	 * @param sp Non-null source placement */
	public void setFrom(SeqPos sp){
		seq=sp.seq();
		pos=sp.pos();
		count=sp.count;
		hashcode=sp.hashcode;
		gc=sp.gc;
		assert(count>=0) : "SeqMap candidate metadata requires a nonnegative occurrence count: "+count;
	}

	/** Compares a homogeneous placement key; this implementation does not guard the cast.
	 * @param o SeqPosM to compare, or null
	 * @return True when position, cached hash and sequence contents match; false for null
	 * @throws ClassCastException If o is non-null and not a SeqPosM */
	@Override
	public boolean equals(Object o){return equals((SeqPosM)o);}

	/** Compares position and sequence contents, independently of count and GC.
	 * @param o Placement with a cached hash consistent with its sequence, or null
	 * @return True when position, cached hash and sequence contents match; false for null */
	public boolean equals(SeqPosM o){
		if(o==null){return false;}
		if(pos!=o.pos || hashcode!=o.hashcode){return false;}
		return Tools.equals(seq, o.seq);
	}

	/** Returns the cached sequence hash; different positions may share a hash.
	 * @return Cached hash, which callers must keep consistent with seq */
	@Override
	public int hashCode(){return hashcode;}

	/** Creates a shallow metadata copy sharing the sequence array.
	 * @return Distinct mutable placement with the same fields and sequence reference */
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
	 * Both placements require non-null sequences; score and GC do not affect ordering.
	 * @param o Non-null placement with a nonnegative count
	 * @return Negative if this ranks first, positive if o ranks first, or zero for a rank tie
	 */
	@Override
	public int compareTo(SeqPosM o){
		if(count!=o.count){return o.count-count;}
		if(seq.length!=o.seq.length){return o.seq.length-seq.length;}
		return 0;
	}

	/** Returns the borrowed sequence.
	 * @return Stored array reference, possibly null for uninitialized scratch */
	public final byte[] seq(){return seq;}
	/** Returns the placement coordinate.
	 * @return Projected sequence start in query coordinates for SeqMap results */
	public final int pos(){return pos;}
	/** Changes position; the object must not be a live hash key.
	 * @param x New coordinate, which may be negative for a sequence extending left of the query */
	public final void setPos(int x){pos=x;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed sequence bytes; null is allowed only for uninitialized scratch. */
	public byte[] seq;
	/** Placement coordinate; SeqMap results store starts projected onto the first query region. */
	public int pos;
	/** Cached sequence-content hash. */
	public int hashcode;
	/** Nonnegative occurrence count. */
	public int count;
	/** Cached GC fraction, potentially NaN when no bases are defined. */
	public float gc;
}
