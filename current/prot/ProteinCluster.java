package prot;

import java.util.ArrayList;

/**
 * A protein family/cluster: one representative sequence plus every member
 * assigned to it, each with its identity and coverage to the representative.
 *
 * <p>This is the object returned by the in-memory clustering API
 * ({@link ProteinClusterer#cluster}). The representative is always the first
 * member ({@code members.get(0)}), recorded at 100% identity and full coverage.
 * A singleton cluster has exactly one member (the representative itself).</p>
 *
 * <p>The cluster id is assigned in representative-creation order and is stable
 * within a single run. A cross-run stable-id / incremental-update policy is a
 * later decision and is deliberately not implemented here.</p>
 *
 * @author Eru
 */
public final class ProteinCluster {

	/** Cluster id, unique and stable within a single run (0-based). */
	public final int id;
	/** The representative protein sequence. */
	public final ProteinSequence representative;
	/** All members, representative first, then assignment order. */
	public final ArrayList<ClusterMember> members=new ArrayList<ClusterMember>();

	/**
	 * Creates a cluster seeded with its representative as the first member
	 * (identity 100%, coverage 1.0).
	 * @param id Cluster id.
	 * @param representative Representative sequence.
	 */
	public ProteinCluster(final int id, final ProteinSequence representative){
		if(representative==null){throw new RuntimeException("Null representative.");}
		this.id=id;
		this.representative=representative;
		members.add(new ClusterMember(representative, 100.0, 1.0));
	}

	/**
	 * Adds a non-representative member with its similarity to the representative.
	 * @param seq Member sequence.
	 * @param identity Percent identity to the representative [0,100].
	 * @param coverage Aligned fraction of the member [0,1].
	 */
	public final void add(final ProteinSequence seq, final double identity, final double coverage){
		members.add(new ClusterMember(seq, identity, coverage));
	}

	/** Number of members, including the representative. @return Member count. */
	public final int size(){return members.size();}

	/** True if the cluster holds only its representative. @return Singleton flag. */
	public final boolean isSingleton(){return members.size()==1;}
}
