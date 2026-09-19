package align2;

/** Worker-local counters for tiered Quantum qualification and opt-in routing. */
public class TieredQuantumStats {

	void observe(final QuantumRanker.Result qr, final int bbScore,
			final int legacyScore, final int legacyStart, final int legacyStop,
			final int originalMinScore){
		attempts++;
		if(!qr.supported){unsupported++;return;}
		supported++;
		if(qr.uncertain){uncertain++;return;}
		assert(qr.match!=null) : "Tiered shadow requires Quantum traceback";
		final int indelBases=qr.insertions+qr.deletions;
		if(legacyScore!=Integer.MIN_VALUE){legacyResult++;}
		if(indelBases<=1){
			smallIndel++;
			final boolean scoreGood=(legacyScore!=Integer.MIN_VALUE && bbScore>=legacyScore);
			final boolean endpointsGood=(qr.rStart==legacyStart && qr.rStop==legacyStop);
			if(scoreGood){smallScoreAtLeastLegacy++;}
			if(endpointsGood){smallEndpointExact++;}
			if(scoreGood && endpointsGood){smallReusable++;}
		}else{
			longerIndel++;
			if(bbScore>originalMinScore){longerRaisesMinScore++;}
			longerEndpointSpanSum+=qr.rStop-qr.rStart+1L;
		}
	}

	void add(final TieredQuantumStats other){
		attempts+=other.attempts;unsupported+=other.unsupported;supported+=other.supported;
		legacyResult+=other.legacyResult;
		uncertain+=other.uncertain;smallIndel+=other.smallIndel;
		smallScoreAtLeastLegacy+=other.smallScoreAtLeastLegacy;
		smallEndpointExact+=other.smallEndpointExact;smallReusable+=other.smallReusable;
		longerIndel+=other.longerIndel;
		longerRaisesMinScore+=other.longerRaisesMinScore;
		longerEndpointSpanSum+=other.longerEndpointSpanSum;
	}

	String toTsv(){
		return "quantum_tiered_shadow\tattempts\t"+attempts+
				"\tsupported\t"+supported+"\tunsupported\t"+unsupported+
				"\tuncertain\t"+uncertain+"\tlegacy_result\t"+legacyResult+
				"\tsmall_indel\t"+smallIndel+
				"\tsmall_score_ge_legacy\t"+smallScoreAtLeastLegacy+
				"\tsmall_endpoint_exact\t"+smallEndpointExact+
				"\tsmall_reusable\t"+smallReusable+
				"\tlonger_indel\t"+longerIndel+
				"\tlonger_raises_minscore\t"+longerRaisesMinScore+
				"\tlonger_endpoint_span_sum\t"+longerEndpointSpanSum+'\n';
	}

	long attempts, unsupported, supported, uncertain, legacyResult;
	long smallIndel, smallScoreAtLeastLegacy, smallEndpointExact, smallReusable;
	long longerIndel, longerRaisesMinScore, longerEndpointSpanSum;
}
