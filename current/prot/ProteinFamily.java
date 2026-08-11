package prot;

import java.util.ArrayList;

/**
 * A protein family produced by consensus-reassignment clustering
 * ({@link ProteinFamilyClusterer}): a consensus sequence that serves as the
 * family's search representative, plus the input members assigned to it.
 *
 * <p>Unlike a greedy {@link ProteinCluster} (whose representative is the longest
 * member), a family's representative is the {@link AAGraph} consensus of all its
 * members, so it fits divergent members rather than only the seed.</p>
 *
 * @author Eru
 */
public final class ProteinFamily {

	/** Family id, assigned in output order (0-based). */
	public final int id;
	/** Consensus sequence (the family's search representative). */
	public final ProteinSequence consensus;
	/** Input members assigned to this family. */
	public final ArrayList<ProteinSequence> members;

	/**
	 * Constructs a family.
	 * @param id Family id.
	 * @param consensus The AAGraph consensus representative.
	 * @param members The assigned member sequences.
	 */
	public ProteinFamily(final int id, final ProteinSequence consensus,
			final ArrayList<ProteinSequence> members){
		this.id=id;
		this.consensus=consensus;
		this.members=members;
	}

	/** Number of members. @return Member count. */
	public final int size(){return members.size();}
}
