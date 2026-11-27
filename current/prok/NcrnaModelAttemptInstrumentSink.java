package prok;

/** Opt-in development trace for ranked model-alignment attempts. */
public interface NcrnaModelAttemptInstrumentSink {

	/** Optional scheduled-window lifecycle. Begin precedes even length/seed rejects.
	 * Shortlist is a borrowed live prefix, in actual one-based verification order.
	 * alignedModel observes actual primary/rescue model verification, not endpoint
	 * refinement kernels. acceptedModel binds the FINAL returned proposal's model.
	 * Completed counts include every nested kernel request, including refinement;
	 * accepted is before candidate resolution and final path selection. No extra
	 * alignment or score computation is performed for these observations.
	 */
	default void windowStarted(String contig, int strand, int pass, int start, int stop){}
	/** Begins one actual ordinary ranking invocation before its index work and
	 * before rankingInputs/shortlistedModels. Stage is supplied by the call site,
	 * including when all models are skipped or a later gate rejects the query.
	 * No event is synthesized for a window rejected before ranking, an explicit
	 * one-model oracle, or a verifier route that performs no new ranking query.
	 * This identifies the query's entry stage, not every later verifier that may
	 * reuse its shortlist. Do not infer a stage from invocation number.
	 */
	default void rankingStarted(int stage){}
	default void shortlistedModels(int[] models, int size){}
	/** Both actual paired-ranking vectors before maximum-score merging or the
	 * ranked cap. Eligible is a borrowed live prefix on the complete model axis;
	 * nativeIndex.lastSharedCount(model) reads this same query's integer counts.
	 * All arguments are read-only and valid only during this callback. Never
	 * rerank, invoke the index or retain these mutable arrays without copying.
	 * Fires per actual ranking invocation, including a repeated edge-stage query;
	 * the consumer must retain invocation identity rather than infer a stage.
	 */
	default void rankingInputs(int[] eligible, int size, float[] memberHits,
		float[] consensusHits, TrnaKmerIndex nativeIndex){}
	default void alignedModel(int model, int stage){}
	default void acceptedModel(int model){}
	/** Actual primary comparison result, after span eligibility and before updating
	 * the winner. Tier -1 means the legacy identity-only comparison has no tiers. */
	default void primaryCompared(int model, float identity, int previousModel, float previousIdentity,
		int tier, int previousTier, boolean replaced){}
	/** Final accepted model/identity and refined oriented inclusive coordinates,
	 * before candidate resolution. This observes the returned proposal only. */
	default void proposalAccepted(int model, float identity, int start, int stop){}
	/** Borrow the actual terminal proposal object before resolution/path selection.
	 * Null means this observed window returned no proposal. No sequence capture. */
	default void proposalCompleted(Orf accepted){}
	/** Existing alignment statistics for explicit shared-only winner experiments.
	 * Reported packed MSDI is separate because integer fields can truncate halves. */
	default void winnerEvidence(int stage, int model, float identity, float matches, float columns, int modelLength, long reportedColumns,
		int previousModel, float previousIdentity, float previousMatches, float previousColumns, int previousModelLength,
		long previousReportedColumns, int tier, int previousTier, int rule, double offset, boolean replaced){}
	/** Physical ordinary primary alignment on the explicit shared winner route.
	 * Captures the input offset before the kernel and the exact returned position/
	 * score fields before span gates. Buffers are borrowed read-only, never a
	 * permission to reset the caller's hint or to run an extra alignment. With
	 * modelQuery=true the reference coordinates are window-relative; otherwise
	 * they are model-relative. An unknown input offset has incomingKnown=false;
	 * its numeric argument is unused and must be serialized as absent, not zero.
	 * Kernel state may depend on history: this event is not a cache-invariance
	 * promise. Cached logical requests and endpoint refinements do not fire it.
	 */
	default void primaryAlignmentState(int model, int kernel, boolean modelQuery,
		boolean incomingKnown, int incomingOffset, float identity, int[] positions){}
	/** Actual one-based primary/physical-edge stopping position, never inferred
	 * from an absent model. Other observers may ignore the optional event. */
	default void patienceStopped(int afterRank, int stage){}
	/** Opt-in score-patience terminal evidence, before candidate finalization.
	 * Ranked visits include pre-alignment length skips; eligible counts do not. */
	default void scoreWalkCompleted(int stage, int visitedRanks, int eligibleAlignments, int reason, boolean ceilingReached){}
	default void windowCompleted(boolean accepted, long requests, long physical, long modelBases){}

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
	int KERNEL_QUANTUM=0, KERNEL_SCRABBLE=1;

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
