package prok;

/**
 * Opt-in capture callback for boundary-NN instrumentation (Citan/Brian, 2026-08-28): invoked
 * once per ACCEPTED ncRNA locus, immediately after trimToAlignmentExtent's alignment-extent
 * snap and before any boundary-NN refinement, with the real production state at that point --
 * winning model, original seed-derived candidate window, and post-trim pre-NN boundaries.
 * Deliberately thin/generic (raw state only, no sweep/label logic) so a future live-policy or
 * DAgger-style capture pass can reuse this same hook with different downstream consumption
 * (Citan's (3): "design the sink so a later NN-on capture can record actual scorer-visited
 * rows" -- this interface does not assume bootstrap-style sweeping is the only consumer).
 *
 * <p>NOT a hot-path interface: invoked only on already-accepted, rare ncRNA calls (bounded by
 * final call count, not scan-window count) -- a single null-check per accepted locus is the
 * only cost when instrumentation is off (the field holding this is null by default; see
 * NcrnaScavenger.instrumentSink).
 *
 * <p>THREAD SAFETY: NcrnaScavenger copies the needed window bytes via Arrays.copyOfRange
 * BEFORE invoking capture() -- windowCopy is a fresh, private array, never a live reference
 * into the mutable per-strand bases[] that GeneCaller reverse-complements in place between
 * strand passes (GeneCaller.java:330). Implementations may retain windowCopy freely.
 * Instrumentation mode itself requires t=1 (enforced by whatever wires this in, not by this
 * interface) -- concurrency safety here is only "the array you were given is exclusively
 * yours for as long as you want it", not a guarantee about concurrent capture() calls.
 *
 * @author G11
 */
public interface NcrnaBoundaryInstrumentSink {

	/**
	 * @param contigName Contig/sequence name (Orf.scafName).
	 * @param strand 0=forward pass, 1=reverse-complement pass (GeneCaller's convention; matches
	 *   Orf.strand at this pre-flip hook point -- GeneCaller.flip() has NOT run yet). All
	 *   coordinates passed to this method are in this per-strand-oriented frame, NOT final
	 *   contig-absolute GFF coordinates. A caller matching against a reference GFF (which IS
	 *   contig-absolute) must transform the reference INTO this frame for strand==1 using
	 *   PFeature.flip()'s formula (start,stop -> scaflen-1-stop, scaflen-1-start) -- that
	 *   formula is self-inverse, so the same one-line transform that converts a scavenged Orf
	 *   OUT of this frame (GeneCaller.java:290-292) converts a truth locus INTO it.
	 * @param modelIndex Index into the family's library[]/models[] of the model this locus was
	 *   actually scored/trimmed against in production -- unmodified by this hook, i.e. the
	 *   real production selection, not a re-derived one.
	 * @param origWStart @param origWStop The ORIGINAL seed-derived candidate window bounds
	 *   (alignWindow's own wStart/wStop, pre-trim), preserved separately from the post-trim
	 *   boundaries so P1 pad/window failures can be distinguished from trim/sweep failures
	 *   downstream (Citan's (5)).
	 * @param postTrimStart @param postTrimStop The Orf's boundaries immediately after
	 *   trimToAlignmentExtent's alignment-extent snap, before any boundary-NN adjustment.
	 * @param windowCopy A fresh, private byte[] copy (see class javadoc) spanning at least
	 *   [postTrimStart-CAPTURE_PAD, postTrimStop+CAPTURE_PAD] clamped to the contig -- generous
	 *   enough for the family-configured boundary candidate ranges and the
	 *   enrichment profile's own local +-2 radius plus k-mer window widths (k up to 11 in the
	 *   currently-staged tables).
	 * @param windowCopyOffset The frame-relative position that windowCopy[0] corresponds to --
	 *   origWStart/origWStop/postTrimStart/postTrimStop are in the SAME frame as windowCopy but
	 *   are NOT pre-translated into windowCopy-relative indices; subtract this offset to index
	 *   into windowCopy (e.g. windowCopy[postTrimStart-windowCopyOffset] is the base at
	 *   postTrimStart).
	 * @param trimSucceeded The STRUCTURAL eligibility signal for downstream bootstrap work. True
	 *   if the alignment-extent snap (ScrabbleAligner traceback against the extended window)
	 *   actually produced a valid rStart/rStop and postTrimStart/postTrimStop were updated from
	 *   it -- a real post-trim position, usable as a bootstrap candidate pending the driver's
	 *   own sweep-reachability check against truth. False means the snap failed (no valid
	 *   traceback, or the extended window was shorter than the model, or a degenerate
	 *   rStart/rStop) and postTrimStart/postTrimStop are simply the RAW pre-trim alignWindow
	 *   span, unchanged -- NOT a real post-trim position, must not be used as if it were one.
	 *   Citan, 2026-08-28: this field exists because the original code silently skipped capture
	 *   entirely for these loci (returning early before the hook), undercounting the accepted
	 *   denominator; capture now always fires, and this flag is what lets a downstream consumer
	 *   tell the two cases apart.
	 * @param nnInvoked A PRODUCTION-RUN FACT, not a bootstrap-suitability signal: true only if
	 *   trimSucceeded AND boundary5Net/boundary3Net were both non-null for this call -- i.e.
	 *   refineBoundaryNN was ACTUALLY invoked immediately after this capture, exactly
	 *   reproducing the original pre-instrumentation condition for running it. Never true when
	 *   trimSucceeded is false (refineBoundaryNN never ran on an untrimmed span, structurally).
	 *   Bootstrap capture runs are intentionally NN-off (no net staged yet), so
	 *   trimSucceeded=true with nnInvoked=false is the NORMAL, expected state during bootstrap
	 *   capture, NOT a degraded one -- do not gate candidate emission or the bootstrap-achievable
	 *   denominator on this field; use trimSucceeded plus the driver's own sweep-reachability
	 *   check instead. This field only tells you whether THIS SPECIFIC RUN happened to have a
	 *   net staged and invoke it -- useful for distinguishing "no net was staged this run" from
	 *   "a net was staged and did/didn't fire", not for judging trainability.
	 */
	void capture(String contigName, int strand, int modelIndex,
			int origWStart, int origWStop, int postTrimStart, int postTrimStop,
			byte[] windowCopy, int windowCopyOffset, boolean trimSucceeded, boolean nnInvoked);
}
