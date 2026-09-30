package align2;

/**
 * Worker-local site counters for the explicit no-MSA Quantum mode.
 * Each gapArray/observe call records one scoring attempt, not one read or final
 * mapping. observe assigns the first applicable outcome, so its outcome buckets
 * are disjoint; unsupported subreasons additionally partition unsupported.
 * Merge only completed workers into a separate total; updates are unsynchronized.
 * @author Collei
 */
public class QuantumOnlyStats{

	/** Records rejection of a compressed-gap site without a Quantum call. */
	void gapArray(){attempts++; compressedGapArray++;}

	/** Consumes a nonnull result before the ranker reuses it. Scores must already
	 * share BBMap's units. Priority: unsupported, uncertain, missing trace, score.
	 * The supplied result is read immediately and never retained or modified. */
	void observe(final QuantumRanker.Result result, final int scaledScore,
			final int minimumScore){
		attempts++;
		if(!result.supported){
			unsupported++;
			if(result.failureCode==QuantumRanker.FAIL_EMPTY){emptyWindow++;}
			else if(result.failureCode==QuantumRanker.FAIL_QUERY_LONGER){queryLonger++;}
			else if(result.failureCode==QuantumRanker.FAIL_WINDOW_TOO_LONG){windowTooLong++;}
			else if(result.failureCode==QuantumRanker.FAIL_EDIT_BUDGET){editBudget++;}
			else{noPath++;}// Includes any failure code without a dedicated bucket.
			return;
		}
		if(result.uncertain){uncertain++;return;}
		if(result.match==null){missingTrace++;return;}
		if(scaledScore<minimumScore){belowThreshold++;return;}
		accepted++;
	}

	/** Adds a completed worker without clearing it or making concurrent updates safe. */
	void add(final QuantumOnlyStats other){
		attempts+=other.attempts; accepted+=other.accepted;
		compressedGapArray+=other.compressedGapArray; unsupported+=other.unsupported;
		emptyWindow+=other.emptyWindow; queryLonger+=other.queryLonger;
		windowTooLong+=other.windowTooLong; editBudget+=other.editBudget;
		noPath+=other.noPath; uncertain+=other.uncertain;
		missingTrace+=other.missingTrace; belowThreshold+=other.belowThreshold;
	}

	/** Emits one newline-terminated record with alternating key/value columns. */
	String toTsv(){
		return "quantum_only\tattempts\t"+attempts+"\taccepted\t"+accepted+
				"\tcompressed_gaps\t"+compressedGapArray+
				"\tunsupported\t"+unsupported+"\tempty_window\t"+emptyWindow+
				"\tquery_longer\t"+queryLonger+
				"\twindow_too_long\t"+windowTooLong+
				"\tedit_budget\t"+editBudget+"\tno_path\t"+noPath+
				"\tuncertain\t"+uncertain+"\tmissing_trace\t"+missingTrace+
				"\tbelow_threshold\t"+belowThreshold+'\n';
	}

	long attempts, accepted, compressedGapArray, unsupported;
	long emptyWindow, queryLonger, windowTooLong, editBudget, noPath;
	long uncertain, missingTrace, belowThreshold;
}
