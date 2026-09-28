package prok;

/** Default-off observations of family seed hits and each alignment window's terminal outcome.
 * Coordinates are zero-based, inclusive, in the oriented strand; scaflen allows conversion
 * to genomic coordinates. ACCEPT means verifier acceptance before snapshot resolution and
 * final gene-path selection. It does not prove a returned or emitted call.
 * Implementations shared by workers must serialize writes. ncrnadiag= opts into this hook.
 * @author G11, Raiden
 */
interface NcrnaStageDiagSink {

	/** One event after family input/minLen guards, including zero hits. The supplied array
	 * is a private copy, so retaining or modifying it cannot change the caller's seed stream. */
	void seed(String family, String contig, int strand, int scaflen, int[] hitPositions);

	/** One event per normal alignWindow return (exceptions remain failed runs).
	 * Outcomes: REJECT_WLEN, REJECT_KHITS, REJECT_INDEX, REJECT_ID, REJECT_LEN,
	 * REJECT_RESCUE, REJECT_VERIFY, REJECT_MAXLEN, ACCEPT_IDPASS, ACCEPT_RESCUE, ACCEPT_HBM.
	 * REJECT_VERIFY is the single-alignment path's combined identity/HBM rejection.
	 * Counts/model indices not reached are -1; uncomputed scores are0 (an attempted
	 * HBM maximum with no eligible score retains -999). bestReId is0 on the single-
	 * alignment path. Rejects may retain an aligned span, otherwise coordinates are -1.
	 * Accept coordinates include configured trim/voting/NN refinement. trimSucceeded
	 * reports only an additional alignment trim, not voting or NN activity. */
	void window(String family, String contig, int strand, int scaflen, int pass,
		int wStart, int wStop, String outcome, int khits, int shortlistSize,
		int modelIndex, String modelName, float bestId, float bestReId, float bestHbm,
		int orfStart, int orfStop, boolean trimSucceeded);
}
