package prot;

/**
 * The CheckM1-style marker-counting completeness/contamination estimator: the
 * ORACLE/baseline of the MAG-QC product. It reads a {@link MarkerVector} (the gene
 * presence/count vector a bin produced against a marker set) and computes, by
 * direct counting, the two frozen quantities and every field that must travel with
 * them (freeze §4d).
 *
 * <p>This is deliberately the simple counting method, not the shipped estimator.
 * The freeze (§3d) locks the shipped completeness/contamination estimator as a
 * NEURAL NET (Brian's design, built on {@code ml.Trainer}); this class is the
 * fixed oracle that net is validated against (§6 invariant 5: correct basic
 * version + locked oracle before optimizing). The output <b>fields</b> are the
 * same for oracle and net; only the estimator differs.</p>
 *
 * <h3>Frozen formulas (freeze §4d)</h3>
 * <ul>
 * <li><b>completeness</b> = detected marker families / effective denominator.</li>
 * <li><b>contamination</b> (headline) = excess copies / effective denominator,
 * where excess copies = &Sigma; max(0, copies-1) over families — so 3 copies
 * contributes 2 and differs from 2 copies (contributes 1).</li>
 * <li><b>effective denominator</b> = expected single-copy marker families for the
 * lineage set; in this basic oracle that is the vector dimension (every position
 * is an expected family).</li>
 * </ul>
 *
 * <p>Also reported (secondary) is the simpler {@code multiCopyMarkers /
 * denominator}; the excess-copy value is the headline number.</p>
 *
 * <p><b>Scope of this oracle (basic per-family only):</b> no collocated-marker-set
 * weighting, no lineage-specific expected-marker adjustment, no taxonomy
 * assignment (so {@code assignmentConfidence} is NaN and {@code lineageTaxID}/
 * {@code rank} are "NA"). Those are the net's / later phases' job.</p>
 *
 * @author Eru
 */
public final class MagQC {

	/** Model+version stamp recorded on every result this oracle produces. */
	public static final String ORACLE_MODEL="oracle-checkm1-count-v1";

	/**
	 * Estimates completeness/contamination from a marker vector alone.
	 * @param mv The bin's gene presence/count vector.
	 * @return The populated QC result.
	 */
	public MagQCResult estimate(final MarkerVector mv){
		return estimate(mv, null);
	}

	/**
	 * Estimates completeness/contamination from a marker vector, using the marker
	 * set (when supplied) only for context fields (id/version); the counting is on
	 * the vector.
	 * @param mv The bin's gene presence/count vector.
	 * @param ms The marker set the vector was scored against, or null.
	 * @return The populated QC result.
	 */
	public MagQCResult estimate(final MarkerVector mv, final MarkerSet ms){
		if(mv==null){throw new RuntimeException("Null MarkerVector.");}

		final int denom=mv.dimension();//effective denominator = expected families
		final int detected=mv.familiesPresent();
		final int multiCopy=mv.familiesMultiCopy();

		//Excess copies = sum over families of max(0, copies-1).
		int excess=0;
		for(final int c : mv.counts){
			if(c<0){throw new RuntimeException("Negative copy count in vector: "+c);}
			if(c>1){excess+=(c-1);}
		}

		//Crash-loud sanity: detected can never exceed the denominator.
		assert(detected<=denom) : "detected="+detected+" > denom="+denom;
		assert(multiCopy<=detected) : "multiCopy="+multiCopy+" > detected="+detected;

		final String domain=mv.domain;
		final String markerSetId=(ms==null ? "NA" :
			(ms.version==null ? "NA" : ms.version));

		//Insufficient-evidence: no expected families -> report unknown, not a
		//fabricated percentage (freeze §4d).
		if(denom==0){
			return new MagQCResult(Double.NaN, Double.NaN, Double.NaN,
				0, detected, multiCopy, excess, 0, domain, markerSetId,
				"NA", "NA", Double.NaN, ORACLE_MODEL, "unknown", false);
		}

		final double completeness=100.0*detected/denom;
		final double contamination=100.0*excess/denom;//headline (excess-copy)
		final double contamMultiCopy=100.0*multiCopy/denom;//secondary

		return new MagQCResult(completeness, contamination, contamMultiCopy,
			denom, detected, multiCopy, excess, denom, domain, markerSetId,
			"NA", "NA", Double.NaN, ORACLE_MODEL, "unknown", true);
	}
}
