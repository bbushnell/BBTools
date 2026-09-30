package align2;

/** Per-worker counters for selective max-indel retries.
 * Each BBMapThread owns one instance; BBMapS merges these after its workers finish.
 * Methods are not synchronized and must not race with aggregation or reporting.
 * A retry counts one single-read attempt or one paired attempt, not two mate ends.
 * Reason flags can overlap and do not partition retries: a no-good-pair retry can
 * occur without any reason bit. Probe and observation counts also include work
 * that never becomes a full wide retry. These are routing diagnostics, not truth
 * accuracy measurements.
 * @author Collei */
public final class HybridMaxIndelStats{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Records entry into a full wide retry, excluding exploratory pseudo probes. */
	public void retryAttempted(){retries++;}
	/** Records selection of the wide attempt after its acceptance checks. */
	public void wideSelected(){wideSelected++;}
	/** Records restoration of the saved low attempt after a rejected wide attempt. */
	public void lowRestored(){lowRestored++;}
	/** Records MAPQ rejection after the other wide-acceptance conditions passed. */
	public void mapqRejected(){mapqRejected++;}
	/** Records a nonrejected single read that remained unmapped at low bounds. */
	public void singleNoAccepted(){singleNoAccepted++;}
	/** Records a wide pseudo probe, not a full alignment or final wide selection.
	 * @param gap Whether the probe found a candidate gap longer than 50 bases
	 * @param pairable Whether the probe found pairable candidate sites */
	public void wideProbe(final boolean gap, final boolean pairable){
		wideProbeAttempted++;
		if(gap){wideProbeGap++;}
		if(pairable){wideProbePairable++;}
		if(gap || pairable){wideProbeAccepted++;}
	}
	/** Records low-pair observations, whether or not they trigger a retry.
	 * @param pairable Whether the current top low sites form an acceptable pair
	 * @param halfSupported Caller-supplied half-read support, or accepted probe
	 * evidence when the experimental wide-probe route is enabled */
	public void pairPseudoObservation(final boolean pairable, final boolean halfSupported){
		if(!pairable){pairNotPairable++;}
		if(halfSupported){pairHalfSupported++;}
		if(!pairable && halfSupported){pairNotPairableHalfSupported++;}
	}
	/** Records overlapping HybridPairPolicy reason bits at paired retry entry.
	 * Zero is valid for the no-good-pair route; both means NO_SITE and HALF_ERRORS. */
	public void pairReasons(final int flags){
		if((flags&HybridPairPolicy.TERMINAL_ERRORS)!=0){pairTerminalErrors++;}
		final boolean noSite=(flags&HybridPairPolicy.NO_SITE)!=0;
		final boolean half=(flags&HybridPairPolicy.HALF_ERRORS)!=0;
		final boolean pseudo=(flags&HybridPairPolicy.PSEUDO_HALF)!=0;
		if(noSite){pairNoSite++;}
		if(half){pairHalfErrors++;}
		if(pseudo){pairPseudoHalf++;}
		if(noSite && half){pairBoth++;}
	}
	/** Records the outcome using the same reason flags as pairReasons.
	 * @param selected True if the wide attempt was retained; false if low was restored */
	public void pairOutcome(final int flags, final boolean selected){
		if((flags&HybridPairPolicy.TERMINAL_ERRORS)!=0){if(selected){pairWideTerminal++;}else{pairRestoredTerminal++;}}
		final boolean noSite=(flags&HybridPairPolicy.NO_SITE)!=0;
		final boolean half=(flags&HybridPairPolicy.HALF_ERRORS)!=0;
		final boolean pseudo=(flags&HybridPairPolicy.PSEUDO_HALF)!=0;
		if(selected){
			if(noSite){pairWideNoSite++;}if(half){pairWideHalfErrors++;}
			if(pseudo){pairWidePseudoHalf++;}
			if(noSite && half){pairWideBoth++;}
		}else{
			if(noSite){pairRestoredNoSite++;}if(half){pairRestoredHalfErrors++;}
			if(pseudo){pairRestoredPseudoHalf++;}
			if(noSite && half){pairRestoredBoth++;}
		}
	}
	/** Records one reused native match, rather than one read pair or retry. */
	public void matchCacheHit(){matchCacheHits++;}

	/** Adds all counters from a completed worker into this aggregate.
	 * The caller merges each worker once after completion; the source is not reset.
	 * @throws IllegalArgumentException If the source statistics are null */
	public void add(final HybridMaxIndelStats b){
		if(b==null){throw new IllegalArgumentException("Cannot add null hybrid max-indel statistics");}
		retries+=b.retries; wideSelected+=b.wideSelected; lowRestored+=b.lowRestored; mapqRejected+=b.mapqRejected;
		singleNoAccepted+=b.singleNoAccepted; pairNotPairable+=b.pairNotPairable;
		wideProbeAttempted+=b.wideProbeAttempted; wideProbeAccepted+=b.wideProbeAccepted;
		wideProbeGap+=b.wideProbeGap; wideProbePairable+=b.wideProbePairable;
		pairHalfSupported+=b.pairHalfSupported; pairNotPairableHalfSupported+=b.pairNotPairableHalfSupported;
		pairNoSite+=b.pairNoSite;
		pairTerminalErrors+=b.pairTerminalErrors; pairWideTerminal+=b.pairWideTerminal; pairRestoredTerminal+=b.pairRestoredTerminal;
		pairHalfErrors+=b.pairHalfErrors; pairPseudoHalf+=b.pairPseudoHalf; pairBoth+=b.pairBoth;
		pairWideNoSite+=b.pairWideNoSite; pairWideHalfErrors+=b.pairWideHalfErrors; pairWideBoth+=b.pairWideBoth;
		pairWidePseudoHalf+=b.pairWidePseudoHalf;
		pairRestoredNoSite+=b.pairRestoredNoSite; pairRestoredHalfErrors+=b.pairRestoredHalfErrors; pairRestoredBoth+=b.pairRestoredBoth;
		pairRestoredPseudoHalf+=b.pairRestoredPseudoHalf;
		matchCacheHits+=b.matchCacheHits;
	}

	/** Returns one newline-terminated row with a prefix and alternating key/value fields.
	 * Retains the established field names and order for existing diagnostic parsers. */
	public String toTsv(){return "hybrid_maxindel\tretries\t"+retries+
		"\twide_selected\t"+wideSelected+"\tlow_restored\t"+lowRestored+
		"\tmapq_rejected\t"+mapqRejected+
		"\tsingle_no_accepted\t"+singleNoAccepted+
		"\twide_probe_attempted\t"+wideProbeAttempted+"\twide_probe_accepted\t"+wideProbeAccepted+
		"\twide_probe_gap\t"+wideProbeGap+"\twide_probe_pairable\t"+wideProbePairable+
		"\tpair_not_pairable\t"+pairNotPairable+"\tpair_half_supported\t"+pairHalfSupported+
		"\tpair_not_pairable_half_supported\t"+pairNotPairableHalfSupported+
		"\tpair_no_site\t"+pairNoSite+
		"\tpair_half_errors\t"+pairHalfErrors+"\tpair_pseudo_half\t"+pairPseudoHalf+"\tpair_both\t"+pairBoth+
		"\tpair_terminal_errors\t"+pairTerminalErrors+"\tpair_wide_terminal\t"+pairWideTerminal+"\tpair_restored_terminal\t"+pairRestoredTerminal+
		"\tpair_wide_no_site\t"+pairWideNoSite+"\tpair_wide_half_errors\t"+pairWideHalfErrors+
		"\tpair_wide_pseudo_half\t"+pairWidePseudoHalf+
		"\tpair_wide_both\t"+pairWideBoth+"\tpair_restored_no_site\t"+pairRestoredNoSite+
		"\tpair_restored_half_errors\t"+pairRestoredHalfErrors+
		"\tpair_restored_pseudo_half\t"+pairRestoredPseudoHalf+"\tpair_restored_both\t"+pairRestoredBoth+
		"\tmatch_cache_hits\t"+matchCacheHits+"\n";}

	/** Number of full wide retries, counting each pair once. */
	public long retryCount(){return retries;}
	/** Number of full wide attempts selected. */
	public long wideSelectedCount(){return wideSelected;}
	/** Number of saved low attempts restored after a full wide retry. */
	public long lowRestoredCount(){return lowRestored;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long retries;
	private long wideSelected;
	private long lowRestored;
	private long mapqRejected;
	private long singleNoAccepted;
	private long wideProbeAttempted;
	private long wideProbeAccepted;
	private long wideProbeGap;
	private long wideProbePairable;
	private long pairNotPairable;
	private long pairHalfSupported;
	private long pairNotPairableHalfSupported;
	private long pairNoSite;
	private long pairHalfErrors;
	private long pairTerminalErrors, pairWideTerminal, pairRestoredTerminal;
	private long pairPseudoHalf;
	private long pairBoth;
	private long pairWideNoSite;
	private long pairWideHalfErrors;
	private long pairWidePseudoHalf;
	private long pairWideBoth;
	private long pairRestoredNoSite;
	private long pairRestoredHalfErrors;
	private long pairRestoredPseudoHalf;
	private long pairRestoredBoth;
	private long matchCacheHits;
}
