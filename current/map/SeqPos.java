package map;

import shared.Tools;

/** Sequence placement with immutable metadata and a borrowed sequence array.
 * Equality uses sequence content and position; natural ordering ranks count,
 * length and score, so equal ranks need not identify the same placement.
 * Callers must not modify the sequence while its cached hash/GC are in use.
 * @author Brian Bushnell
 * @date 2013
 */
public class SeqPos implements Cloneable, Comparable<SeqPos>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies metadata and borrows the same sequence array.
	 * @param sp Non-null placement to copy */
	public SeqPos(SeqPos sp){this(sp.seq, sp.pos, sp.count, sp.hashcode, sp.gc, sp.score);}

	/** Constructs an indexed placement with precomputed sequence metadata.
	 * @param seq_ Borrowed sequence bytes; must be non-null for ordering and indexing
	 * @param pos_ Indexed kmer end position, or a caller-defined placement coordinate
	 * @param count_ Nonnegative occurrence count
	 * @param code_ Precomputed sequence hash
	 * @param gc_ Precomputed GC fraction
	 * @param score_ Optional score; SeqMap uses -1 for an unscored sequence
	 */
	SeqPos(byte[] seq_, int pos_, int count_, int code_, float gc_, float score_){
		seq=seq_;
		pos=pos_;
		count=count_;
		hashcode=code_;
		gc=gc_;
		score=score_;
		assert(count>=0) : "Sequence occurrence counts must be nonnegative for descending count ranking: "+count;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares a homogeneous placement key; this implementation does not guard the cast.
	 * @param o SeqPos to compare, or null
	 * @return True when position, cached hash and sequence contents match; false for null
	 * @throws ClassCastException If o is non-null and not a SeqPos */
	@Override
	public boolean equals(Object o){return equals((SeqPos)o);}

	/** Compares position and sequence contents, independently of ranking metadata.
	 * Cached hashes must agree before sequence bytes are compared.
	 * @param o Placement with a hash consistent with its sequence, or null
	 * @return True when position, cached hash and sequence contents match; false for null */
	public boolean equals(SeqPos o){
		if(o==null){return false;}
		if(pos!=o.pos || hashcode!=o.hashcode){return false;}
		return Tools.equals(seq, o.seq);
	}

	/** Returns the cached sequence hash; different positions may share a hash.
	 * @return Constructor-supplied sequence hash */
	@Override
	public int hashCode(){return hashcode;}

	/** Creates a shallow metadata copy sharing the borrowed sequence.
	 * @return Distinct placement with the same metadata and sequence reference */
	@Override
	public SeqPos clone(){
		try{
			return (SeqPos)super.clone();
		}catch(CloneNotSupportedException e){
			throw new RuntimeException(e);
		}
	}

	/** Ranks count, length and score descending. Signed zeros tie; NaN sorts before
	 * numeric scores, following reversed Float ordering. This ranking intentionally
	 * differs from the sequence/position identity used by equals and hashCode.
	 * Both placements require nonnegative counts and non-null sequences.
	 * @param o Non-null placement to rank against
	 * @return Negative if this ranks first, positive if o ranks first, or zero for a rank tie
	 */
	@Override
	public int compareTo(SeqPos o){
		if(count!=o.count){return o.count-count;}
		if(seq.length!=o.seq.length){return o.seq.length-seq.length;}
		if(score!=o.score){return Float.compare(o.score, score);}
		return 0;
	}

	/** Returns the borrowed sequence array; its contents must remain stable.
	 * @return Stored sequence reference, not a copy */
	public final byte[] seq(){return seq;}
	/** Returns the stored position.
	 * @return Indexed kmer end coordinate for SeqMap entries */
	public final int pos(){return pos;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed sequence; its contents must remain stable while indexed. */
	public final byte[] seq;
	/** Kmer end coordinate for SeqMap entries. */
	public final int pos;
	/** Cached sequence-content hash. */
	public final int hashcode;
	/** Nonnegative occurrence count. */
	public final int count;
	/** GC fraction; may be NaN for sequences with no defined bases. */
	public final float gc;
	/** Optional rank score. */
	public final float score;
}
