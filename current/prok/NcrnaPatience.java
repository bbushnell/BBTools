package prok;

/** Reusable, worker-owned stagnation budget for one ncRNA primary model walk.
 * Only span-eligible actual alignments participate. Strict raw-identity gains
 * reset the counter, including before the first pass; latching a pass alone
 * does not reset it. Winning-model selection is deliberately independent.
 * Reset between primary walks; never share between workers or nested active walks.
 * Failed-primary rescue does not use this budget and retains its first-pass exit.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPatience{

	/** -1 leaves stopping to the legacy caller; MAX walks the entire shortlist. */
	public void reset(final int limit_){
		validate(limit_);
		limit=limit_;
		bestIdentity=0;
		sinceImproved=0;
		passed=false;
	}

	/** Returns a stop only after an eligible primary pass has been latched in this walk.
	 * Ineligible attempts consume nothing, even when their identity is higher.
	 */
	public boolean observe(final float identity, final boolean currentPass, final boolean eligible){
		if(limit<0 || !eligible){return false;}
		assert(Float.isFinite(identity) && identity>=0 && identity<=1)
			: "Patience compares actual finite alignment identities in [0,1], not scores or missing results: "+identity;
		if(identity>bestIdentity){
			bestIdentity=identity;
			sinceImproved=0;
		}else if(sinceImproved<Integer.MAX_VALUE){sinceImproved++;}
		passed|=currentPass;
		return limit!=Integer.MAX_VALUE && passed && sinceImproved>=limit;
	}

	public static void validate(final int value){
		if(value<-1){throw new IllegalArgumentException("ncRNA patience must be -1 (legacy) or a nonnegative alignment count: "+value);}
	}

	private int limit=-1, sinceImproved;
	private float bestIdentity;
	private boolean passed;
}
