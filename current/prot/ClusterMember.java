package prot;

/**
 * One member of a {@link ProteinCluster}: a protein sequence plus its measured
 * similarity to the cluster representative.
 *
 * <p>The representative is itself a member of its own cluster, recorded with
 * 100% identity and full coverage. All other members carry the identity and
 * coverage of their best BLOSUM62 local alignment against the representative,
 * as computed by {@link AAAligner}.</p>
 *
 * @author Eru
 */
public final class ClusterMember {

	/** The member protein sequence. */
	public final ProteinSequence seq;
	/** Percent identity to the representative in [0,100] (100 for the rep itself). */
	public final double identity;
	/**
	 * Fraction of this member's residues covered by the aligned span against the
	 * representative, in [0,1] (1.0 for the rep itself).
	 */
	public final double coverage;

	/**
	 * Constructs a cluster member.
	 * @param seq Member protein sequence.
	 * @param identity Percent identity to the representative [0,100].
	 * @param coverage Aligned fraction of this member [0,1].
	 */
	public ClusterMember(final ProteinSequence seq, final double identity, final double coverage){
		if(seq==null){throw new RuntimeException("Null member sequence.");}
		this.seq=seq;
		this.identity=identity;
		this.coverage=coverage;
	}

	/** Member id (convenience). @return The member sequence id. */
	public final String id(){return seq.id;}
}
