package structures;

import shared.Tools;

/**
 * Mutable SeqCount that accumulates occurrences and optional scores.
 * Allows incrementing counts for repeated sequences while retaining sequence bytes.
 * Used where SeqCount immutability is insufficient for aggregation.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
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
	
	/**
	 * Creates a synchronized clone of this SeqCountM.
	 * Preserves the sequence data, count, and score values.
	 * @return Deep copy of this SeqCountM
	 */
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
	/** Adds the count from another SeqCount to this instance (no sequence copy).
	 * @param s SeqCount whose count should be accumulated */
	public void add(SeqCount s) {
//		assert(equals(s));
		count+=s.count();
	}

//	@Override
	/** Increments the occurrence count by the specified amount.
	 * @param x Amount to add to the current count */
	public void increment(int x) {
		count+=x;
	}

	/** Returns the current occurrence count for this sequence */
	@Override
	public int count() {return count;}
	
	/**
	 * Compares SeqCounts for ordering, with enhanced score comparison.
	 * Primary comparison by count, then sequence length, then score (for SeqCountM),
	 * finally by lexicographic sequence comparison.
	 *
	 * @param s SeqCount to compare against
	 * @return Negative, zero, or positive for less than, equal, or greater than
	 */
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
	
	/** Number of occurrences observed for this sequence (mutable). */
	public int count=1;
	/** Optional score associated with this sequence; defaults to -1 when unset. */
	public float score=-1;
	
}