package prok;

/** Default-off observations of family seed hits and each alignment window's terminal outcome.
 * Coordinates are zero-based, inclusive, in the oriented strand; scaflen allows conversion
 * to genomic coordinates. ACCEPT means verifier acceptance before snapshot resolution and
 * final gene-path selection. It does not prove a returned or emitted call.
 * Implementations shared by workers must serialize writes. ncrnadiag= opts into this hook.
 * The TSV implementation emits exactly21fields under its type_v2 header.
 * Unavailable numbers are-1; unavailable text/list/trim flags are empty.
 * Explicit operation flags preserve real measured0 and attempted-trim false;
 * nonfinite measurements fail rather than becoming an absent observation.
 * @author G11, Raiden
 */
public interface NcrnaStageDiagSink {

	/** One event after family input/minLen guards, including zero hits. The supplied array
	 * is a private copy, so retaining or modifying it cannot change the caller's seed stream. */
	void seed(String family, String contig, int strand, int scaflen, int[] hitPositions);

	/** One event per normal alignWindow return (exceptions remain failed runs).
	 * Outcomes: REJECT_WLEN, REJECT_KHITS, REJECT_INDEX, REJECT_ID, REJECT_LEN,
	 * REJECT_RESCUE, REJECT_VERIFY, REJECT_MAXLEN, ACCEPT_IDPASS, ACCEPT_RESCUE, ACCEPT_HBM.
	 * REJECT_VERIFY is the single-alignment path's combined identity/HBM rejection.
	 * Unreached numeric fields are -1. Measured zero is retained; measuredFlags
	 * comes from actual primary/rescue/HBM/trim operations, never from score values.
	 * An unattempted trim has no text value; an attempted failed trim is false.
	 * Rejects may retain an aligned span, otherwise coordinates are -1.
	 * Accept coordinates normally include configured trim/voting/NN refinement.
	 * Glycine instead records raw discovery/selection bounds; deferred reporting
	 * endpoints appear in the final GFF only after selection. trimSucceeded
	 * reports only an additional alignment trim, not voting or NN activity. */
	void window(String family, String contig, int strand, int scaflen, int pass,
		int wStart, int wStop, String outcome, int khits, int shortlistSize,
		int modelIndex, String modelName, float bestId, float bestReId, float bestHbm,
		int orfStart, int orfStop, int measuredFlags, boolean trimSucceeded);

	/** Explicit current-window measurement presence; the TSV remains21columns. */
	int PRIMARY_ID=1, RESCUE_ID=2, HBM=4, TRIM=8;

	/** Actual generic-RNA DP input or retained output, after scavenger resolution.
	 * Coordinates retain the oriented, zero-based inclusive convention of window().
	 * Glycine retained rows still use raw selection bounds, before deferred reporting.
	 * Only primitive/immutable metadata crosses the observer boundary; no caller
	 * candidate or mutable list is exposed. Existing observers may ignore this event. */
	default void path(String family, String contig, int strand, int scaflen,
		String modelName, int start, int stop, boolean retained){}

	/** Selected verifier identity carried by the exact16S/23S/prok5S candidate through DP.
	 * Other families do not supply this measurement; the TSV omits it as numeric-1.
	 * Existing observers retain their old path callback. */
	default void path(String family,String contig,int strand,int scaflen,
		String modelName,int start,int stop,boolean retained,float acceptedIdentity){
		path(family,contig,strand,scaflen,modelName,start,stop,retained);
	}

	/** One initial alignment actually attempted for a shortlisted model. Identity
	 * is observed before length/identity acceptance; it is not a counterfactual
	 * for models outside the shortlist. Rescue realignments are not reported here. */
	default void model(String family,String contig,int strand,int scaflen,int pass,
		int wStart,int wStop,int modelIndex,String modelName,float identity){}

	/** Supplemental experimental local-alignment geometry, independent of the
	 * stable 21-column stage stream. Null means no positive alignment. */
	default void modelClip(String family,String contig,int strand,int pass,int wStart,int wStop,
		int model,idaligner.EndClippedAligner.Result result,int minLen){}

	/** Supplemental PacBio evidence and separately accounted Quantum comparison. */
	default void modelPacBio(String family,String contig,int strand,int pass,int wStart,int wStop,
		int model,Euk18sPacBioAligner.Result result,float cutoff,int maxLen){}
}
