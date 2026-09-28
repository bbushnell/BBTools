package prok;

/** Opt-in development trace for ranked model-alignment attempts. */
public interface NcrnaModelAttemptInstrumentSink {

	/**
	 * Fires after each initial shortlist alignment. Ranks are one-based and sharedKmers is the
	 * exact index-kmer count used to order that model for this window. alignedLength is the
	 * inclusive reference span reported by QuantumAligner, or -1 when the short-window
	 * Scrabble path did not request positional output.
	 * In the explicitly enabled exhaustive mapping-oracle arm, rank is instead the
	 * one-based library traversal position; every model is attempted without ranking
	 * or early exit. Index scratch is refreshed, so sharedKmers still belongs to this
	 * window (the one-model direct-alignment control reports0 without counting).
	 */
	void modelAttempt(String contigName, int strand, int pass, int windowStart, int windowStop,
		int rank, int model, int sharedKmers, float identity, int alignedLength, boolean identityPass);
}
