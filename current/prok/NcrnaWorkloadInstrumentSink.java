package prok;

/**
 * Opt-in capture callback for B4 seed-trigger/workload instrumentation (Citan's schema B4,
 * 2026-08-29: "the GATE/WORKLOAD metric = DEDUPLICATED SCHEDULED candidate windows, per shred
 * AND per Mbp -- NOT raw seed-hit occurrences, NOT bbduk contaminated-read counts"). Mirrors
 * NcrnaBoundaryInstrumentSink's opt-in pattern exactly (same null-field-by-default cost model,
 * same reason: a hook that never fires in any production run must cost nothing when off).
 *
 * <p>Deliberately raw/generic (coordinates only, no dedup/counting logic here) -- the driver
 * owns aggregation. In particular the "distinct physical shred" dedup Citan specified
 * (a shred is triggered if EITHER strand hits) cannot be done inside NcrnaScavenger: one
 * instance only ever sees ONE strand's bases[] per scavenge() call (GeneCaller reverse-
 * complements between passes), so it has no way to know whether the OTHER strand of the same
 * physical shred also hit. The driver, which calls scavenge() once per strand per shred, is the
 * only place that can correlate both calls back to one physical shred ID.
 *
 * <p>NOT a hot-path interface in the traditional sense (fires once per scavenge() call for
 * seedHits, once per scheduled window for scheduledWindow -- bounded by shred count and
 * candidate-window count, not by genome length), but IS on the main scan path unlike the
 * accepted-only boundary sink -- a single null-check per call site is still the only cost when
 * off (field null by default; see NcrnaScavenger.workloadSink).
 *
 * <p>THREAD SAFETY: hitPositions is defensively COPIED (Arrays.copyOf) by NcrnaScavenger before
 * this callback fires -- Citan's correction, 2026-08-29: the original javadoc here claimed no
 * copy was needed because findKmerHitPositions returns "a fresh array per call", which is true
 * of the array findKmerHitPositions returns but was never actually verified against IntList's
 * internal buffer-reuse behavior before being asserted -- an unverified claim about shared
 * production code, exactly the kind of thing this project's whole discipline exists to catch.
 * The copy removes the dependency on that internal detail entirely rather than relying on it.
 * Instrumentation mode requires t=1 (enforced by whatever wires this in, not by this interface),
 * same as the boundary sink.
 *
 * @author G11
 */
public interface NcrnaWorkloadInstrumentSink {

	/**
	 * Fires once per scavenge() invocation that reaches findKmerHitPositions -- including when
	 * hitPositions.length==0 -- immediately after findKmerHitPositions returns, before the
	 * length==0 early-return that would otherwise skip the rest of scavenge(). This is NOT
	 * unconditional across every scavenge() call (Citan's correction, 2026-08-29): scavenge()'s
	 * OWN earlier guard (`bases==null || bases.length&lt;minLen || library==null ||
	 * kmerSet==null`) returns before findKmerHitPositions is even called, so a shred shorter
	 * than the family's minLen (a real, common case -- e.g. every L100 shred against rnasep's
	 * minLen=200) produces NO seedHits call at all, correctly reflecting that the real production
	 * caller never scanned it. A driver MUST count totalShreds/totalBp from the physical input
	 * records directly, never by counting seedHits callbacks -- a guarded-out record contributes
	 * zero to every callback-derived metric while still being a real physical shred.
	 * @param contigName Contig/sequence name (the shred ID in this project's usage).
	 * @param strand 0=forward pass, 1=reverse-complement pass (GeneCaller's convention; same
	 *   frame as NcrnaBoundaryInstrumentSink.capture's strand parameter).
	 * @param hitPositions ALL conserved-kmer hit coordinates found on this strand, in the
	 *   per-strand-oriented frame findKmerHitPositions computed them in (see
	 *   NcrnaBoundaryInstrumentSink's strand javadoc for the frame-transform note if matching
	 *   against a contig-absolute reference). Length 0 is valid and expected for most shreds.
	 */
	void seedHits(String contigName, int strand, int[] hitPositions);

	/**
	 * Fires once per candidate window immediately BEFORE NcrnaScavenger calls alignWindow on it
	 * -- i.e. once per actual alignment REQUEST, after buildCandidateWindows -> collapse
	 * ByIntersection -> subtractClaimed have already run, matching schema B4's "actual alignment
	 * requests" gate metric exactly (NOT raw seed-hit count, NOT a pre-collapse window count).
	 * Fires for every window in both scavenge passes -- pass=1 for the primary windows loop,
	 * pass=2 for the findNearbyUnclaimed/pass2Windows loop (only reached when scavengePass2 is
	 * enabled, the production default) -- so a driver can report combined and per-pass totals.
	 * @param contigName Contig/sequence name (the shred ID).
	 * @param strand 0=forward, 1=reverse-complement (same frame as seedHits above).
	 * @param pass 1 or 2, identifying which scavenge() pass scheduled this window.
	 * @param wStart @param wStop The scheduled window's bounds (alignWindow's own wStart/wStop
	 *   argument for this call), in the same per-strand frame as seedHits' hitPositions.
	 */
	void scheduledWindow(String contigName, int strand, int pass, int wStart, int wStop);
}
