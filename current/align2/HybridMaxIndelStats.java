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
	public void pairReasons(int flags){
		final boolean noSite=(flags&HybridPairPolicy.NO_SITE)!=0;
		final boolean half=(flags&HybridPairPolicy.HALF_ERRORS)!=0;
		if(noSite){pairNoSite++;}
		if(half){pairHalfErrors++;}
		if(noSite && half){pairBoth++;}
	}
	public void pairOutcome(int flags,boolean selected){
		final boolean noSite=(flags&HybridPairPolicy.NO_SITE)!=0;
		final boolean half=(flags&HybridPairPolicy.HALF_ERRORS)!=0;
		if(selected){
			if(noSite){pairWideNoSite++;}if(half){pairWideHalfErrors++;}
			if(noSite && half){pairWideBoth++;}
		}else{
			if(noSite){pairRestoredNoSite++;}if(half){pairRestoredHalfErrors++;}
			if(noSite && half){pairRestoredBoth++;}
		}
	}
	public void matchCacheHit(){matchCacheHits++;}

	public void add(HybridMaxIndelStats b){
		if(b==null){throw new IllegalArgumentException("Cannot add null hybrid max-indel statistics");}
		retries+=b.retries;wideSelected+=b.wideSelected;lowRestored+=b.lowRestored;mapqRejected+=b.mapqRejected;
		singleNoAccepted+=b.singleNoAccepted;pairNoSite+=b.pairNoSite;
		pairHalfErrors+=b.pairHalfErrors;pairBoth+=b.pairBoth;
		pairWideNoSite+=b.pairWideNoSite;pairWideHalfErrors+=b.pairWideHalfErrors;pairWideBoth+=b.pairWideBoth;
		pairRestoredNoSite+=b.pairRestoredNoSite;pairRestoredHalfErrors+=b.pairRestoredHalfErrors;pairRestoredBoth+=b.pairRestoredBoth;
		matchCacheHits+=b.matchCacheHits;
	}

	public String toTsv(){return "hybrid_maxindel\tretries\t"+retries+
		"\twide_selected\t"+wideSelected+"\tlow_restored\t"+lowRestored+
		"\tmapq_rejected\t"+mapqRejected+
		"\tsingle_no_accepted\t"+singleNoAccepted+"\tpair_no_site\t"+pairNoSite+
		"\tpair_half_errors\t"+pairHalfErrors+"\tpair_both\t"+pairBoth+
		"\tpair_wide_no_site\t"+pairWideNoSite+"\tpair_wide_half_errors\t"+pairWideHalfErrors+
		"\tpair_wide_both\t"+pairWideBoth+"\tpair_restored_no_site\t"+pairRestoredNoSite+
		"\tpair_restored_half_errors\t"+pairRestoredHalfErrors+"\tpair_restored_both\t"+pairRestoredBoth+
		"\tmatch_cache_hits\t"+matchCacheHits+"\n";}

	public long retryCount(){return retries;}
	public long wideSelectedCount(){return wideSelected;}
	public long lowRestoredCount(){return lowRestored;}

	private long retries;
	private long wideSelected;
	private long lowRestored;
	private long mapqRejected;
	private long singleNoAccepted;
	private long pairNoSite;
	private long pairHalfErrors;
	private long pairBoth;
	private long pairWideNoSite;
	private long pairWideHalfErrors;
	private long pairWideBoth;
	private long pairRestoredNoSite;
	private long pairRestoredHalfErrors;
	private long pairRestoredBoth;
	private long matchCacheHits;
}
