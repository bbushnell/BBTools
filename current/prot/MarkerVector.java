package prot;

/**
 * The gene presence/count feature vector a genome bin produces against one
 * {@link MarkerSet}: for each selected single-copy marker family, how many of the
 * bin's proteins matched it (the copy count; 0 = absent).
 *
 * <p><b>This is the primary output of {@link MarkerVectorizer}.</b> It is the
 * fixed-length, comparably-indexed vector the completeness/contamination neural
 * net reads. The dimensions and order are defined by the marker set's selected
 * markers (family-id ascending), so every bin scored against the same marker set
 * yields a vector whose entry {@code i} always refers to the same family.</p>
 *
 * <p>The obvious CheckM1-style derived scalars — families present, present exactly
 * once, and multi-copy (&gt;1) — are computed on demand from {@link #counts}; the
 * raw count vector remains the load-bearing feature.</p>
 *
 * @author Eru
 */
public final class MarkerVector {

	/** counts[i] = number of bin proteins assigned to selected marker family i. */
	public final int[] counts;
	/** familyIds[i] = the marker family id at vector position i (parallel to counts). */
	public final int[] familyIds;
	/** representativeIds[i] = that family's representative id (parallel to counts). */
	public final String[] representativeIds;
	/** Domain of the marker set this vector was scored against. */
	public final String domain;
	/** Bin proteins that matched some marker family (assigned to a count). */
	public final int proteinsMatched;
	/** Bin proteins that matched no marker family above threshold. */
	public final int proteinsUnmatched;

	/**
	 * Constructs a marker vector.
	 * @param counts Per-family copy counts (defines the dimension).
	 * @param familyIds Family ids parallel to counts.
	 * @param representativeIds Family representative ids parallel to counts.
	 * @param domain Marker-set domain.
	 * @param proteinsMatched Count of bin proteins assigned to a family.
	 * @param proteinsUnmatched Count of bin proteins assigned to no family.
	 */
	public MarkerVector(final int[] counts, final int[] familyIds,
			final String[] representativeIds, final String domain,
			final int proteinsMatched, final int proteinsUnmatched){
		if(counts==null || familyIds==null || representativeIds==null){
			throw new RuntimeException("Null array in MarkerVector.");
		}
		if(counts.length!=familyIds.length || counts.length!=representativeIds.length){
			throw new RuntimeException("MarkerVector array length mismatch: counts="+
				counts.length+" familyIds="+familyIds.length+" reps="+representativeIds.length);
		}
		this.counts=counts;
		this.familyIds=familyIds;
		this.representativeIds=representativeIds;
		this.domain=domain;
		this.proteinsMatched=proteinsMatched;
		this.proteinsUnmatched=proteinsUnmatched;
	}

	/** Vector length = number of selected marker families. @return Dimension. */
	public final int dimension(){return counts.length;}

	/** Marker families with at least one copy (counts[i]&gt;0). @return Present count. */
	public final int familiesPresent(){
		int n=0;
		for(final int c : counts){if(c>0){n++;}}
		return n;
	}

	/** Marker families present exactly once (counts[i]==1). @return Single-copy count. */
	public final int familiesExactlyOnce(){
		int n=0;
		for(final int c : counts){if(c==1){n++;}}
		return n;
	}

	/** Marker families with more than one copy (counts[i]&gt;1). @return Multi-copy count. */
	public final int familiesMultiCopy(){
		int n=0;
		for(final int c : counts){if(c>1){n++;}}
		return n;
	}
}
