package structures;

/**
 * SeqCount with mutable occurrence and score fields.
 * <p>
 * Equality and hashing are inherited from SeqCount and remain sequence-only.
 * Ordering delegates to SeqCount's count, length, score and bases comparison.
 * Mutation is unsynchronized; callers sharing an instance coordinate externally.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class SeqCountM extends SeqCount{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a mutable counter retaining and canonizing sq's base array, then copying its current count. */
	public SeqCountM(SeqCount sq){
		super(sq.bases);
		count=sq.count();
	}

	/** Creates a mutable counter from a copied half-open sequence slice. */
	public SeqCountM(byte[] s, int start, int stop){
		super(s, start, stop);
	}

	/** Creates a mutable counter retaining and canonizing the supplied bases. */
	public SeqCountM(byte[] s){
		super(s);
	}

	/** Returns a shallow clone that shares the final base array. */
	@Override
	public SeqCountM clone(){
		synchronized(this){
			SeqCountM clone=(SeqCountM) super.clone();
			return clone;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

//	@Override
	/**
	 * Adds another SeqCount's count to this one.
	 * Caller must supply the same sequence; only the count is read.
	 */
	public void add(SeqCount s){
//		assert(equals(s));
		count+=s.count();
	}

//	@Override
	/** Adds x to the occurrence count. */
	public void increment(int x){
		count+=x;
	}

	/** Returns the current occurrence count. */
	@Override
	public int count(){return count;}

	/** Returns the current comparison score. */
	@Override
	public float score(){return score;}

	/**
	 * Orders by count, length, score, then canonicalized bases.
	 */
	@Override
	public int compareTo(SeqCount s){
		return super.compareTo(s);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Mutable occurrence count. */
	public int count=1;

	/** Mutable optional score; -1 when unset. */
	public float score=-1;

}
