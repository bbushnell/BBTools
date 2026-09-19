package align2;

/** Worker-local counters for the explicit no-MSA Quantum mode. */
public class QuantumOnlyStats {

	void gapArray(){attempts++;compressedGapArray++;}

	void observe(final QuantumRanker.Result result, final int scaledScore,
			final int minimumScore){
		attempts++;
		if(!result.supported){
			unsupported++;
			if(result.failureCode==QuantumRanker.FAIL_EMPTY){emptyWindow++;}
			else if(result.failureCode==QuantumRanker.FAIL_QUERY_LONGER){queryLonger++;}
			else if(result.failureCode==QuantumRanker.FAIL_WINDOW_TOO_LONG){windowTooLong++;}
			else if(result.failureCode==QuantumRanker.FAIL_EDIT_BUDGET){editBudget++;}
			else{noPath++;}
			return;
		}
		if(result.uncertain){uncertain++;return;}
		if(result.match==null){missingTrace++;return;}
		if(scaledScore<minimumScore){belowThreshold++;return;}
		accepted++;
	}

	void add(final QuantumOnlyStats other){
		attempts+=other.attempts;accepted+=other.accepted;
		compressedGapArray+=other.compressedGapArray;unsupported+=other.unsupported;
		emptyWindow+=other.emptyWindow;queryLonger+=other.queryLonger;
		windowTooLong+=other.windowTooLong;editBudget+=other.editBudget;
		noPath+=other.noPath;uncertain+=other.uncertain;
		missingTrace+=other.missingTrace;belowThreshold+=other.belowThreshold;
	}

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
