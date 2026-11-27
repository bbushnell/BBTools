package prok;

/** Query admission through the two already-loaded reference indexes. No bitset,
 * duplicate postings or shipped-profile changes. Empty key zero is not a word.
 * @author Brian, Ganyu
 */
public final class NcrnaDdlReferenceKeys{

	public NcrnaDdlReferenceKeys(final NcrnaDdlProfileResource.ConsensusPair pair_){
		if(pair_==null || !pair_.memberOnly.orderOnly() || !pair_.replacement.orderOnly()
			|| !pair_.memberOnly.sameOrderingQuery(pair_.replacement)
			|| pair_.memberOnly.models()!=pair_.replacement.models()){
			throw new IllegalArgumentException("Reference-key admission requires one validated member/consensus-replacement pair");
		}
		pair=pair_;
		bins=pair.memberOnly.buckets();
		seed=pair.memberOnly.seed();
	}

	/** Existing posting presence is checked before query-bin maximum reduction. */
	public boolean contains(final int bin, final int key){return pair.containsReferenceKey(bin, key);}
	/** Exact object binding prevents reusing a posting lookup from another reference axis. */
	public boolean matches(final NcrnaDdlProfileResource.ConsensusPair other){return pair==other;}
	/** No additional array payload; the small wrapper still has ordinary object overhead. */
	public long payloadBytes(){return 0;}

	public final int bins;
	public final long seed;
	private final NcrnaDdlProfileResource.ConsensusPair pair;
}
