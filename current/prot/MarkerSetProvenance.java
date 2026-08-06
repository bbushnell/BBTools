package prot;

import java.util.ArrayList;
import java.util.List;

/**
 * Provenance stamp for one {@link MarkerSet}: the parameters and source genome
 * ids that produced it, plus a caller-supplied build timestamp.
 *
 * <p>Per the freeze (§4e), the factory never reads the clock itself — the build
 * timestamp is passed in so a rebuild is reproducible. This is the basic
 * provenance slice; the full immutable build manifest (input/output checksums,
 * commit + dirty-tree state, taxonomy snapshot id, accession split lists) is a
 * later step and its remaining fields are left as documented slots.</p>
 *
 * @author Eru
 */
public final class MarkerSetProvenance {

	/** Caller-supplied build timestamp (never read from the clock); may be "NA". */
	public final String buildTimestamp;
	/** Single-copy selection threshold on the exactly-once fraction, in [0,1]. */
	public final double selectionThreshold;
	/** Clustering minimum percent identity used. */
	public final double minIdentity;
	/** Clustering minimum coverage used. */
	public final double minCoverage;
	/** Clustering seed k-mer length used. */
	public final int k;
	/** Source genome ids contributing to this domain's set (input order). */
	public final List<String> sourceGenomeIds;
	/** Taxonomy snapshot id (slot; "NA" until the reference-corpus step). */
	public final String taxonomyVersion;

	/**
	 * Constructs a provenance stamp.
	 * @param buildTimestamp Caller-supplied timestamp (nullable -&gt; "NA").
	 * @param selectionThreshold Exactly-once selection threshold [0,1].
	 * @param minIdentity Clustering min percent identity.
	 * @param minCoverage Clustering min coverage.
	 * @param k Clustering seed k-mer length.
	 * @param sourceGenomeIds Contributing genome ids (copied defensively).
	 * @param taxonomyVersion Taxonomy snapshot id (nullable -&gt; "NA").
	 */
	public MarkerSetProvenance(final String buildTimestamp, final double selectionThreshold,
			final double minIdentity, final double minCoverage, final int k,
			final List<String> sourceGenomeIds, final String taxonomyVersion){
		this.buildTimestamp=(buildTimestamp==null ? "NA" : buildTimestamp);
		this.selectionThreshold=selectionThreshold;
		this.minIdentity=minIdentity;
		this.minCoverage=minCoverage;
		this.k=k;
		this.sourceGenomeIds=new ArrayList<String>(
			sourceGenomeIds==null ? new ArrayList<String>() : sourceGenomeIds);
		this.taxonomyVersion=(taxonomyVersion==null ? "NA" : taxonomyVersion);
	}
}
