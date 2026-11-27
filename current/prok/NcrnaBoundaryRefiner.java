package prok;

import idaligner.AlignmentStats;

/** Optional worker-owned refinement of a selected ncRNA model's endpoints.
 * Called once at the accepted-window boundary seam, outside the alignment loops
 * and before gene-path selection. No implementation is loaded by default.
 * Implementations may update start/stop, but must preserve the detection score,
 * model identity, source bases and any physically truncated endpoint.
 * @author Brian Bushnell, Raiden
 */
public interface NcrnaBoundaryRefiner{

	/** Whether this refiner owns the selected model. False retains the ordinary
	 * endpoint path; true replaces that path even when refinement makes no move. */
	boolean applies(int model);

	/** False preserves the original pre-trim replacement hook. True requests the
	 * ordinary scalar-NN anchor after unchanged alignment trimming, and runs only
	 * when trim succeeded, raw-end refinement is enabled or one end is fixed.
	 * Selected before calling; must not change during a worker's lifetime.
	 */
	default boolean afterAlignmentTrim(){return false;}

	/** Opt-in trace borrowers run synchronously before these worker-owned stats
	 * can be reused. Existing endpoint methods do not retain alignment objects. */
	default boolean needsAlignmentTrace(){return false;}

	/** Both origins are zero-based on the complete sense-oriented contig.
	 * accepting=true starts a new selected-model call and clears old trim state. */
	default void alignment(AlignmentStats stats, int origin, boolean accepting){}

	/** Covariance endpoints use no genomic-GC feature; existing networks do. */
	default boolean usesContigGc(){return true;}

	/** Only an implementation with an explicit trace contract may publish score windows. */
	default void setBoundaryTrace(NcrnaBoundaryTrace.Family trace){
		if(trace!=null){throw new IllegalArgumentException("This refiner has no boundary trace contract");}
	}

	/** Optional invocation diagnostics, separate from NN score-window replay.
	 * Implementations that do not publish proposals need no diagnostic state. */
	default void setProposalTrace(NcrnaBoundaryTrace.Family trace){}

	/** Refines zero-based inclusive coordinates on the complete sense-oriented
	 * contig. fixedMask uses bit 1 for the physical start and bit 2 for the stop.
	 * Identity is the selected alignment's score; quantumIdentity distinguishes
	 * full-model Quantum anchors from other alignment policies. Length limits are
	 * inclusive bases. Returns whether coordinates changed, not detection success.
	 */
	boolean refine(Orf orf, byte[] bases, int model, float gc, float identity,
		boolean quantumIdentity, int fixedMask, int minLen, int maxLen);

	/** Shared read-only factory installed before caller workers start. Each call
	 * validates its family and returns a new refiner with independent mutable state.
	 * Never mutate the family or its shared resources from this method.
	 */
	interface Factory{
		NcrnaBoundaryRefiner create(NcrnaFamily family);
	}
}
