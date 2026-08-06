package prot;

import java.util.List;

/**
 * Turns a genome bin's proteins into a gene presence/count feature vector against
 * a {@link MarkerSet} — the feature-extraction step that feeds the MAG-QC
 * completeness/contamination neural net (Capability 4).
 *
 * <p><b>This is the load-bearing API.</b> {@link #vectorize(List, MarkerSet)}
 * takes a bin's proteins as in-memory {@link ProteinSequence} objects plus a
 * marker set and returns a {@link MarkerVector} — no disk round-trip; callable
 * directly from other BBTools code (e.g. a binner scoring an in-memory bin).</p>
 *
 * <p>Algorithm (basic version, correctness-first):</p>
 * <ul>
 * <li>The vector's fixed dimension and order are the marker set's selected
 * single-copy markers ({@link MarkerSet#selectedMarkers()}, family-id ascending),
 * so every bin scored against the same marker set is comparably indexed.
 * <li>Each bin protein is aligned against every selected family representative
 * with BLOSUM62 affine-gap local alignment ({@link AAAligner}) — the same
 * identity/coverage measure {@link ProteinClusterer} uses to build the families.
 * <li>A protein is assigned to its single best-matching family (highest percent
 * identity) that meets both the identity and coverage thresholds; ties go to the
 * earliest family in vector order, so the result is deterministic. That family's
 * count is incremented. A protein that meets no family's thresholds matches
 * nothing.
 * </ul>
 *
 * <p>Thresholds default to the ones recorded in the marker set's provenance (the
 * thresholds that defined the families), so a bin is scored consistently with how
 * the markers were built; set {@link #minIdentity}/{@link #minCoverage} to a
 * non-negative value to override, falling back to 90% / 0.8 when neither is
 * available.</p>
 *
 * <p>Deferred (reported, not built): real bin/{@code CallGenes} hookup (genome
 * FASTA -&gt; proteins), the completeness/contamination net itself (this only
 * produces its input), threshold tuning, k-mer-seeded / multithreaded matching
 * for speed, and multi-family (paralog) assignment.</p>
 *
 * @author Eru
 */
public final class MarkerVectorizer {

	/** Minimum percent identity to assign a protein to a family; -1 = inherit. */
	public double minIdentity=-1;
	/** Minimum aligned fraction of both protein and representative; -1 = inherit. */
	public double minCoverage=-1;
	/** Fallback identity threshold when neither override nor provenance supplies one. */
	public static final double DEFAULT_MIN_IDENTITY=90.0;
	/** Fallback coverage threshold when neither override nor provenance supplies one. */
	public static final double DEFAULT_MIN_COVERAGE=0.8;

	/**
	 * Builds the gene presence/count vector for one bin against one marker set.
	 *
	 * @param binProteins The bin's protein sequences (in-memory; non-null).
	 * @param markerSet The marker set defining the vector's families and order.
	 * @return The fixed-length count vector plus derived match statistics.
	 */
	public MarkerVector vectorize(final List<ProteinSequence> binProteins,
			final MarkerSet markerSet){
		if(binProteins==null){throw new RuntimeException("Null bin protein list.");}
		if(markerSet==null){throw new RuntimeException("Null marker set.");}

		final List<MarkerFamily> markers=markerSet.selectedMarkers();
		final int dim=markers.size();
		final int[] counts=new int[dim];
		final int[] familyIds=new int[dim];
		final String[] repIds=new String[dim];
		for(int i=0; i<dim; i++){
			final MarkerFamily f=markers.get(i);
			familyIds[i]=f.familyId;
			repIds[i]=f.representative.id;
		}

		final double idThresh=resolveIdentity(markerSet);
		final double covThresh=resolveCoverage(markerSet);

		int matched=0, unmatched=0;
		for(final ProteinSequence p : binProteins){
			if(p==null){throw new RuntimeException("Null protein in bin.");}
			int best=-1;
			double bestIdentity=-1;
			for(int i=0; i<dim; i++){
				final ProteinSequence rep=markers.get(i).representative;
				final AAAlignment aln=AAAligner.align(p.enc, rep.enc);
				if(aln==null){continue;}
				final double identity=aln.pident();
				if(identity<idThresh){continue;}
				final double memberCov=(aln.qStop-aln.qStart+1)/(double)p.length();
				final double repCov=(aln.tStop-aln.tStart+1)/(double)rep.length();
				if(memberCov<covThresh || repCov<covThresh){continue;}
				//Strict '>' makes the earliest family in vector order win ties (deterministic).
				if(identity>bestIdentity){bestIdentity=identity; best=i;}
			}
			if(best>=0){counts[best]++; matched++;}
			else{unmatched++;}
		}
		return new MarkerVector(counts, familyIds, repIds, markerSet.domain, matched, unmatched);
	}

	/** Resolves the identity threshold: override, else provenance, else default. */
	private double resolveIdentity(final MarkerSet ms){
		if(minIdentity>=0){return minIdentity;}
		if(ms.provenance!=null && ms.provenance.minIdentity>0){return ms.provenance.minIdentity;}
		return DEFAULT_MIN_IDENTITY;
	}

	/** Resolves the coverage threshold: override, else provenance, else default. */
	private double resolveCoverage(final MarkerSet ms){
		if(minCoverage>=0){return minCoverage;}
		if(ms.provenance!=null && ms.provenance.minCoverage>0){return ms.provenance.minCoverage;}
		return DEFAULT_MIN_COVERAGE;
	}
}
