package prot;

/**
 * One protein family as recorded in a {@link MarkerSet}: its stable-within-run
 * id, representative sequence, the {@link CopyNumberDistribution} across the
 * marker set's domain, and whether it was selected as a single-copy marker.
 *
 * @author Eru
 */
public final class MarkerFamily {

	/** Family id — the clustering cluster id; stable within a run only. */
	public final int familyId;
	/** Representative protein sequence of the family. */
	public final ProteinSequence representative;
	/** Copy-number distribution across this domain's genomes. */
	public final CopyNumberDistribution dist;
	/** True if selected as a single-copy marker (exactly-once fraction &ge; threshold). */
	public final boolean selectedSingleCopy;

	/**
	 * Constructs a marker family record.
	 * @param familyId Cluster id (stable within run).
	 * @param representative Family representative.
	 * @param dist Copy-number distribution over the domain's genomes.
	 * @param selectedSingleCopy Whether it qualifies as a single-copy marker.
	 */
	public MarkerFamily(final int familyId, final ProteinSequence representative,
			final CopyNumberDistribution dist, final boolean selectedSingleCopy){
		if(representative==null){throw new RuntimeException("Null representative.");}
		if(dist==null){throw new RuntimeException("Null distribution.");}
		this.familyId=familyId;
		this.representative=representative;
		this.dist=dist;
		this.selectedSingleCopy=selectedSingleCopy;
	}

	/** Prevalence in this domain [0,1]. @return Fraction of genomes carrying it. */
	public final double prevalence(){return dist.prevalence();}

	/** Exactly-once fraction in this domain [0,1]. @return Single-copy fraction. */
	public final double fractionExactlyOnce(){return dist.fractionExactlyOnce();}
}
