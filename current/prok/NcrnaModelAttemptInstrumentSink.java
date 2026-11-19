package prok;

/** Opt-in development trace for ranked model-alignment attempts. */
public interface NcrnaModelAttemptInstrumentSink {

	/**
	 * Fires after each initial shortlist alignment. Ranks are one-based and sharedKmers is the
	 * exact index-kmer count used to order that model for this window, or -1 when a
	 * sketch-only caller deliberately has no native index. alignedLength is the
	 * inclusive reference span reported by QuantumAligner, or -1 when the short-window
	 * Scrabble path did not request positional output.
	 * In the explicitly enabled exhaustive mapping-oracle arm, rank is instead the
	 * one-based library traversal position; every model is attempted without ranking
	 * or early exit. Index scratch is refreshed, so sharedKmers still belongs to this
	 * window (the one-model direct-alignment control reports0 without counting).
	 */
	void modelAttempt(String contigName, int strand, int pass, int windowStart, int windowStop,
		int rank, int model, int sharedKmers, float identity, int alignedLength, boolean identityPass);

	/** Raw actual alignment endpoints before identity/span gates, trimming or path
	 * selection. Coordinates are oriented zero-based inclusive. Invalid/missing
	 * coordinates are explicit -1/-1. Stage distinguishes primary, edge and rescue;
	 * genuineClip/edgeFloor apply only to physical-edge attempts. No alignment is
	 * added to produce this observation. Existing observers retain their contract. */
	default void alignmentSpan(String contigName, int strand, int pass, int windowStart, int windowStop,
		int model, float identity, int alignedStart, int alignedStop, int stage,
		boolean genuineClip, boolean edgeFloor){}

	int PRIMARY=0, PHYSICAL_EDGE=1, RESCUE=2;

	/** Selected-format rescue uses one Quantum trace for identity, span and HBM.
	 * The original window and the narrower rescue operand are both oriented,
	 * zero-based inclusive intervals. hbmMeasured is authoritative; hbm is read
	 * only when that flag is true and is then finite, including a real zero.
	 * An unmeasured HBM value is omitted by any serialized observer. This callback
	 * is separate from primary attempts; existing primary-only observers receive
	 * no additional rows. Accepted means this attempt passed rescue verification,
	 * before final winning-model, boundary and path selection.
	 */
	default void rescueAttempt(String contigName, int strand, int pass, int windowStart, int windowStop,
		int candidateStart, int candidateStop, int rank, int model, float identity,
		int alignedLength, float hbm, boolean hbmMeasured, boolean accepted){}
}
