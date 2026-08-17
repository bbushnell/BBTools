package prot;

/**
 * A single protein-search hit: one row of BLAST-tab {@code outfmt 6} output,
 * plus the raw score and a flag recording whether the E-value is rigorous.
 *
 * <p>Coordinates are stored 1-based inclusive (BLAST convention). Field order on
 * output is the frozen contract:
 * {@code query target pident length mismatch gapopen qstart qend tstart tend evalue bitscore}.</p>
 *
 * <p>This is the object returned by the in-memory search API — callers consume it
 * directly without any disk round-trip.</p>
 *
 * @author Eru
 */
public final class ProteinHit {

	/** Query sequence id. */
	public final String query;
	/** Target sequence id. */
	public final String target;
	/** Percent identity in [0,100]. */
	public final double pident;
	/** Alignment length (columns, gaps included). */
	public final int length;
	/** Mismatch count (aligned non-gap, non-identical columns). */
	public final int mismatch;
	/** Gap-open count (contiguous gap runs). */
	public final int gapopen;
	/** Query start (1-based, inclusive). */
	public final int qstart;
	/** Query end (1-based, inclusive). */
	public final int qend;
	/** Target start (1-based, inclusive). */
	public final int tstart;
	/** Target end (1-based, inclusive). */
	public final int tend;
	/** E-value; may be approximate (see evalueApproximate). */
	public final double evalue;
	/** Bitscore (rigorous for gapped BLOSUM62 11/1). */
	public final double bitscore;
	/** Raw BLOSUM62 alignment score. */
	public final int rawScore;
	/**
	 * True if the E-value omits the Karlin-Altschul edge-effect length correction
	 * (an approximation). The bitscore is rigorous regardless.
	 */
	public final boolean evalueApproximate;

	/**
	 * Builds a hit from an alignment result and its query/target metadata.
	 * @param query Query id.
	 * @param target Target id.
	 * @param aln Alignment result (0-based coords internally).
	 * @param evalue Computed E-value.
	 * @param evalueApproximate Whether the E-value is approximate.
	 */
	public ProteinHit(String query, String target, AAAlignment aln,
			double evalue, boolean evalueApproximate){
		this.query=query;
		this.target=target;
		this.pident=aln.pident();
		this.length=aln.length;
		this.mismatch=aln.mismatches;
		this.gapopen=aln.gapOpens;
		this.qstart=aln.qStart+1;
		this.qend=aln.qStop+1;
		this.tstart=aln.tStart+1;
		this.tend=aln.tStop+1;
		this.rawScore=aln.rawScore;
		this.bitscore=aln.bitscore();
		this.evalue=evalue;
		this.evalueApproximate=evalueApproximate;
	}

	/**
	 * Formats this hit as a BLAST-tab {@code outfmt 6} row (no trailing newline).
	 * Precision follows the frozen contract: pident 3 decimals, E-value 2
	 * significant figures in scientific notation, bitscore 1 decimal.
	 *
	 * <p>Kept on StringBuilder/String.format deliberately (assessed for the
	 * String.format-&gt;ByteBuilder conformance pass, not converted): pident is
	 * 100.000 on every exact/self hit -- routine, not an edge case -- and
	 * ByteBuilder.append(double,int)'s whole-number special-case would print "100"
	 * there (same defect confirmed empirically on ClusterProteins.row() and
	 * MarkerFactoryCLI.row()). {@link #formatEvalue} also has no ByteBuilder
	 * equivalent at all (scientific notation isn't one of its numeric formatters).
	 * Output-preservation wins over the stylistic conversion here too.</p>
	 * @return Tab-separated row.
	 */
	public final String toTsv(){
		final StringBuilder sb=new StringBuilder();
		sb.append(query).append('\t');
		sb.append(target).append('\t');
		sb.append(String.format("%.3f", pident)).append('\t');
		sb.append(length).append('\t');
		sb.append(mismatch).append('\t');
		sb.append(gapopen).append('\t');
		sb.append(qstart).append('\t');
		sb.append(qend).append('\t');
		sb.append(tstart).append('\t');
		sb.append(tend).append('\t');
		sb.append(formatEvalue(evalue)).append('\t');
		sb.append(formatBits(bitscore));
		return sb.toString();
	}

	/**
	 * Formats an E-value in scientific notation with 2 significant figures.
	 * @param e E-value.
	 * @return Formatted string.
	 */
	static final String formatEvalue(double e){
		if(e==0){return "0.0e+00";}
		return String.format("%.1e", e);
	}

	/**
	 * Formats a bitscore with 1 decimal place.
	 * @param b Bitscore.
	 * @return Formatted string.
	 */
	static final String formatBits(double b){
		return String.format("%.1f", b);
	}

	@Override
	public String toString(){return toTsv();}
}
