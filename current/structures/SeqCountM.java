package structures;

import shared.Tools;

/**
 * Mutable version of SeqCount that tracks variable occurrence counts and scores.
 * While the parent SeqCount is immutable with fixed count of 1, SeqCountM allows
 * incremental counting of sequence occurrences and maintains an optional score.
 * Used for applications requiring accumulation of sequence statistics.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class SeqCountM extends SeqCount {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a mutable SeqCountM from an existing SeqCount.
	 * Copies the sequence data and initializes count to the source count.
	 * @param sq Source SeqCount to copy from
	 */
	public SeqCountM(SeqCount sq) {
		super(sq.bases);
		count=sq.count();
	}
	
	/**
	 * Creates a SeqCountM from a subsequence of bases with initial count of 1.
	 * @param s Source sequence array
	 * @param start Starting index (inclusive)
	 * @param stop Ending index (exclusive)
	 */
	public SeqCountM(byte[] s, int start, int stop) {
		super(s, start, stop);
	}
	
	/** Creates a SeqCountM from a complete sequence array with initial count of 1.
	 * @param s Sequence bases to track */
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
	/**
	 * Adds the count from another SeqCount to this SeqCountM.
	 * Accumulates occurrence counts for the same sequence.
	 * @param s SeqCount whose count should be added to this one
	 */
	public void add(SeqCount s) {
//		assert(equals(s));
		count+=s.count();
	}

//	@Override
	/** Increments the count by the specified amount.
	 * @param x Amount to add to the current count */
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
	
	/** Number of occurrences observed for this sequence */
	public int count=1;
	/** Optional score associated with this sequence, defaults to -1 */
	public float score=-1;
	
}