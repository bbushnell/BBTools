package prot;

import java.util.ArrayList;
import java.util.List;

/**
 * A per-domain, versioned, provenance-stamped set of protein families with their
 * per-domain prevalence and copy-number distributions, and a selected subset of
 * single-copy marker families.
 *
 * <p>The freeze (§4c) requires marker sets to be stored per domain (Bacteria and
 * Archaea separately) because the two domains have distinct single-copy marker
 * genes; therefore {@link MarkerFactory} returns one {@code MarkerSet} per domain
 * rather than a single combined object.</p>
 *
 * <p>Basic slice: the lineage hierarchy (domain -&gt; ... -&gt; species) with
 * nearest-ancestor fallback, the marker-set 4-tuple key, and cross-run stable
 * family ids are later steps and are not implemented here.</p>
 *
 * @author Eru
 */
public final class MarkerSet {

	/** Domain this set describes (e.g. "Bacteria"). */
	public final String domain;
	/** Version id for this marker set (caller-supplied). */
	public final String version;
	/** Provenance stamp (params, source genome ids, build timestamp). */
	public final MarkerSetProvenance provenance;
	/** Number of genomes of this domain (the prevalence denominator). */
	public final int genomeCount;
	/** All families observed in this domain (prevalence &gt; 0), family-id ascending. */
	public final List<MarkerFamily> families;

	/**
	 * Constructs a marker set.
	 * @param domain Domain label.
	 * @param version Version id.
	 * @param provenance Provenance stamp.
	 * @param genomeCount Genomes of this domain.
	 * @param families Observed families (family-id ascending).
	 */
	public MarkerSet(final String domain, final String version,
			final MarkerSetProvenance provenance, final int genomeCount,
			final List<MarkerFamily> families){
		this.domain=domain;
		this.version=version;
		this.provenance=provenance;
		this.genomeCount=genomeCount;
		this.families=families;
	}

	/**
	 * The families selected as single-copy markers.
	 * @return Newly-built list of selected families (family-id order preserved).
	 */
	public final List<MarkerFamily> selectedMarkers(){
		final ArrayList<MarkerFamily> out=new ArrayList<MarkerFamily>();
		for(final MarkerFamily f : families){if(f.selectedSingleCopy){out.add(f);}}
		return out;
	}

	/** Number of selected single-copy markers. @return Selected count. */
	public final int selectedCount(){
		int n=0;
		for(final MarkerFamily f : families){if(f.selectedSingleCopy){n++;}}
		return n;
	}
}
