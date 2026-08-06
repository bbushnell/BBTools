package prot;

import java.util.List;

/**
 * One genome's contribution to the marker factory: a labeled bag of that
 * genome's protein sequences. The domain (Bacteria / Archaea) and lineage label
 * drive per-domain, per-lineage marker statistics.
 *
 * <p>Protein sequence ids must be unique across the whole pool handed to
 * {@link MarkerFactory} (the clustering engine fails loud on duplicate ids). The
 * factory maps each clustered protein back to its genome by object identity, so
 * an id only needs to be pool-unique; it need not encode the genome.</p>
 *
 * @author Eru
 */
public final class GenomeProteins {

	/** Genome identifier (e.g. an accession); non-empty. */
	public final String id;
	/** Domain label ("Bacteria" or "Archaea"); non-empty (per-domain marker sets). */
	public final String domain;
	/** Lineage label (e.g. a TaxID or name); may be null when unused. */
	public final String lineage;
	/** This genome's protein sequences. */
	public final List<ProteinSequence> proteins;

	/**
	 * Constructs a labeled genome.
	 * @param id Genome identifier (non-empty).
	 * @param domain Domain label (non-empty).
	 * @param lineage Lineage label (nullable).
	 * @param proteins This genome's protein sequences (non-null).
	 */
	public GenomeProteins(final String id, final String domain, final String lineage,
			final List<ProteinSequence> proteins){
		if(id==null || id.length()==0){throw new RuntimeException("Empty genome id.");}
		if(domain==null || domain.length()==0){
			throw new RuntimeException("Empty domain for genome '"+id+"'.");
		}
		if(proteins==null){throw new RuntimeException("Null protein list for genome '"+id+"'.");}
		this.id=id;
		this.domain=domain;
		this.lineage=lineage;
		this.proteins=proteins;
	}

	/** Number of proteins in this genome. @return Protein count. */
	public final int size(){return proteins.size();}
}
