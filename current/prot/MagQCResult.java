package prot;

/**
 * The completeness/contamination QC estimate for a single genome bin, plus every
 * field the freeze (§4d) requires to travel with it. This is a plain immutable
 * result object: {@link MagQC} produces it from a {@link MarkerVector}, and both
 * the in-process caller and {@link MagQCCLI} read it.
 *
 * <p>Barbara's invariant (freeze §4d, §6.3): a completeness/contamination
 * percentage is <b>never</b> reported alone. The raw counts, the effective
 * denominator, the domain, the marker-set id, the lineage, an
 * {@link #assignmentConfidence}, and an {@link #oodStatus} all travel with the
 * two percentages so an alternate estimator remains auditable.</p>
 *
 * <p><b>Scale:</b> {@link #completeness}, {@link #contamination}, and
 * {@link #contaminationMultiCopy} are percentages on the frozen {@code [0,100]}
 * scale (§4d). Contamination is <b>not</b> clamped: a heavily duplicated bin can
 * legitimately exceed 100% (excess copies can outnumber expected families), and
 * hiding that would be fabrication — the true value is reported as-is.</p>
 *
 * <p><b>Insufficient evidence:</b> when {@link #effectiveDenominator} is 0 the
 * percentages are {@code NaN} and {@link #sufficientEvidence} is false — the
 * freeze forbids reporting a fabricated percentage when marker support is absent
 * (§4d "unknown / insufficient-evidence behavior").</p>
 *
 * @author Eru
 */
public final class MagQCResult {

	/** Completeness percentage in [0,100]: detected families / effective denom. */
	public final double completeness;
	/** HEADLINE contamination percentage: excess copies / effective denom (not clamped). */
	public final double contamination;
	/** Secondary contamination percentage: multi-copy families / effective denom. */
	public final double contaminationMultiCopy;

	/** Expected marker families for the lineage set (= vector dimension in the oracle). */
	public final int expectedMarkers;
	/** Detected marker families (families with at least one copy). */
	public final int detectedMarkers;
	/** Marker families observed in more than one copy. */
	public final int multiCopyMarkers;
	/** Excess copies: sum over families of max(0, copies-1) (3 copies contributes 2). */
	public final int excessCopies;
	/** Effective denominator: expected single-copy marker families (oracle: dimension). */
	public final int effectiveDenominator;

	/** Domain the marker set belongs to (Bacteria / Archaea); the domain assignment. */
	public final String domainAssignment;
	/** Marker-set id/version the estimate was produced against ("NA" if unknown). */
	public final String markerSetId;
	/** Lineage TaxID of the marker set ("NA" in the basic oracle — no lineage input). */
	public final String lineageTaxID;
	/** Lineage rank of the marker set ("NA" in the basic oracle). */
	public final String rank;

	/**
	 * P(taxonomy assignment correct). {@code NaN} in this oracle: the CheckM1-style
	 * counting method assigns no taxonomy, so it emits no confidence — that is the
	 * job of the §3d neural net, not this baseline. Reported (not omitted) and
	 * honestly {@code NaN} rather than a fabricated number.
	 */
	public final double assignmentConfidence;
	/** Model+version that produced {@link #assignmentConfidence}. */
	public final String assignmentConfidenceModel;

	/**
	 * Out-of-distribution status. Fixed to {@code "unknown"} because OD-9 (the
	 * OOD-validation + abstention protocol) is not closed; the freeze forbids
	 * populating this with anything else until then (§4d, §7 OD-9).
	 */
	public final String oodStatus;

	/** False when there was no marker support (denominator 0) and % are NaN. */
	public final boolean sufficientEvidence;

	/**
	 * Constructs a fully-populated QC result.
	 * @param completeness Completeness % [0,100].
	 * @param contamination Headline (excess-copy) contamination %.
	 * @param contaminationMultiCopy Secondary (multi-copy-family) contamination %.
	 * @param expectedMarkers Expected marker families.
	 * @param detectedMarkers Detected marker families.
	 * @param multiCopyMarkers Families in &gt;1 copy.
	 * @param excessCopies Sum of max(0, copies-1).
	 * @param effectiveDenominator Effective denominator (expected families).
	 * @param domainAssignment Domain of the marker set.
	 * @param markerSetId Marker-set id/version.
	 * @param lineageTaxID Lineage TaxID ("NA" in the oracle).
	 * @param rank Lineage rank ("NA" in the oracle).
	 * @param assignmentConfidence P(assignment correct); NaN in the oracle.
	 * @param assignmentConfidenceModel Model+version producing the confidence.
	 * @param oodStatus OOD status (fixed "unknown" until OD-9).
	 * @param sufficientEvidence Whether marker support existed.
	 */
	public MagQCResult(final double completeness, final double contamination,
			final double contaminationMultiCopy, final int expectedMarkers,
			final int detectedMarkers, final int multiCopyMarkers, final int excessCopies,
			final int effectiveDenominator, final String domainAssignment,
			final String markerSetId, final String lineageTaxID, final String rank,
			final double assignmentConfidence, final String assignmentConfidenceModel,
			final String oodStatus, final boolean sufficientEvidence){
		this.completeness=completeness;
		this.contamination=contamination;
		this.contaminationMultiCopy=contaminationMultiCopy;
		this.expectedMarkers=expectedMarkers;
		this.detectedMarkers=detectedMarkers;
		this.multiCopyMarkers=multiCopyMarkers;
		this.excessCopies=excessCopies;
		this.effectiveDenominator=effectiveDenominator;
		this.domainAssignment=(domainAssignment==null ? "NA" : domainAssignment);
		this.markerSetId=(markerSetId==null ? "NA" : markerSetId);
		this.lineageTaxID=(lineageTaxID==null ? "NA" : lineageTaxID);
		this.rank=(rank==null ? "NA" : rank);
		this.assignmentConfidence=assignmentConfidence;
		this.assignmentConfidenceModel=
			(assignmentConfidenceModel==null ? "NA" : assignmentConfidenceModel);
		this.oodStatus=(oodStatus==null ? "unknown" : oodStatus);
		this.sufficientEvidence=sufficientEvidence;
	}
}
