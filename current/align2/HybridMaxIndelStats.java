package align2;

/** Per-worker counters for selective max-indel retries.
 * @author Collei
 */
public final class HybridMaxIndelStats {

	public void retryAttempted(){retries++;}
	public void wideSelected(){wideSelected++;}
	public void lowRestored(){lowRestored++;}
	public void mapqRejected(){mapqRejected++;}
	public void singleNoAccepted(){singleNoAccepted++;}
	public void wideProbe(boolean gap,boolean pairable){
		wideProbeAttempted++;
		if(gap){wideProbeGap++;}
		if(pairable){wideProbePairable++;}
		if(gap || pairable){wideProbeAccepted++;}
	}
	public void pairPseudoObservation(boolean pairable,boolean halfSupported){
		if(!pairable){pairNotPairable++;}
		if(halfSupported){pairHalfSupported++;}
		if(!pairable && halfSupported){pairNotPairableHalfSupported++;}
	}
	public void pairReasons(int flags){
		final boolean noSite=(flags&HybridPairPolicy.NO_SITE)!=0;
		final boolean half=(flags&HybridPairPolicy.HALF_ERRORS)!=0;
		final boolean pseudo=(flags&HybridPairPolicy.PSEUDO_HALF)!=0;
		if(noSite){pairNoSite++;}
		if(half){pairHalfErrors++;}
		if(pseudo){pairPseudoHalf++;}
		if(noSite && half){pairBoth++;}
	}
	public void pairOutcome(int flags,boolean selected){
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
	public void matchCacheHit(){matchCacheHits++;}

	public void add(HybridMaxIndelStats b){
		if(b==null){throw new IllegalArgumentException("Cannot add null hybrid max-indel statistics");}
		retries+=b.retries;wideSelected+=b.wideSelected;lowRestored+=b.lowRestored;mapqRejected+=b.mapqRejected;
		singleNoAccepted+=b.singleNoAccepted;pairNotPairable+=b.pairNotPairable;
		wideProbeAttempted+=b.wideProbeAttempted;wideProbeAccepted+=b.wideProbeAccepted;
		wideProbeGap+=b.wideProbeGap;wideProbePairable+=b.wideProbePairable;
		pairHalfSupported+=b.pairHalfSupported;pairNotPairableHalfSupported+=b.pairNotPairableHalfSupported;
		pairNoSite+=b.pairNoSite;
		pairHalfErrors+=b.pairHalfErrors;pairPseudoHalf+=b.pairPseudoHalf;pairBoth+=b.pairBoth;
		pairWideNoSite+=b.pairWideNoSite;pairWideHalfErrors+=b.pairWideHalfErrors;pairWideBoth+=b.pairWideBoth;
		pairWidePseudoHalf+=b.pairWidePseudoHalf;
		pairRestoredNoSite+=b.pairRestoredNoSite;pairRestoredHalfErrors+=b.pairRestoredHalfErrors;pairRestoredBoth+=b.pairRestoredBoth;
		pairRestoredPseudoHalf+=b.pairRestoredPseudoHalf;
		matchCacheHits+=b.matchCacheHits;
	}

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
		"\tpair_wide_no_site\t"+pairWideNoSite+"\tpair_wide_half_errors\t"+pairWideHalfErrors+
		"\tpair_wide_pseudo_half\t"+pairWidePseudoHalf+
		"\tpair_wide_both\t"+pairWideBoth+"\tpair_restored_no_site\t"+pairRestoredNoSite+
		"\tpair_restored_half_errors\t"+pairRestoredHalfErrors+
		"\tpair_restored_pseudo_half\t"+pairRestoredPseudoHalf+"\tpair_restored_both\t"+pairRestoredBoth+
		"\tmatch_cache_hits\t"+matchCacheHits+"\n";}

	public long retryCount(){return retries;}
	public long wideSelectedCount(){return wideSelected;}
	public long lowRestoredCount(){return lowRestored;}

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
