package prok;

/** Opt-in RF00174 candidate observation immediately before model shortlisting.
 * A driver binds one sink to one CallGenesIteration input worker, retaining that
 * input's identity and current Read ordinal outside the family-local scavenger.
 * Both arrays are borrowed, read-only and valid only during the callback; a
 * consumer retaining sequence must copy it. The contig is strand-oriented, while
 * window is the caller's exact uppercase copy of the inclusive [start,stop] span.
 * Original-contig context is available for observation without copying genomes
 * into a run directory. This interface does not change candidates or acceptance.
 * @author Brian Bushnell, Raiden
 */
public interface NcrnaWindowCaptureSink{

	/** Observes a post-seed-gate event before index filtering or early acceptance.
	 * Repeated events and the two scavenger passes remain distinct.
	 */
	void capture(String contig, int strand, int pass, int start, int stop,
		int source, byte[] contigBases, byte[] window);

	/** Completes one actual scheduled window after all primary/rescue/edge/trim
	 * alignments. A rejected seed/length gate may have no preceding capture and
	 * zero work. The accepted object is borrowed before candidate resolution and
	 * final path selection; null means no accepted proposal. No work is replayed.
	 * The count is logical caller requests; it equals physical kernels only with
	 * alignment caching disabled, as the selected-results collector requires.
	 */
	default void completed(String contig, int strand, int pass, int start, int stop,
		Orf accepted, long alignments){}
}
