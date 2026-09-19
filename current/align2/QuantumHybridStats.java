package align2;

/** Worker-local routing counters for the broad Quantum/MSA hybrid. */
public class QuantumHybridStats {

	void ordinaryAttempt(){ordinaryAttempts++;}
	void accepted(){accepted++;}
	void fallbackUnsupported(){fallback++;unsupported++;}
	void fallbackUncertain(){fallback++;selectedTraceUncertain++;}
	void fallbackMissing(){fallback++;missingTrace++;}
	void fallbackThreshold(){fallback++;belowThreshold++;}
	void longTraceFallback(){fallback++;longTraceFallback++;}
	void geometryFallback(){fallback++;geometryFallback++;}
	void compressedEvaluated(){compressedEvaluated++;}
	void compressedSkipped(){compressedSkipped++;strictCompressedSkipped++;}
	void scoreAwareSkipped(){compressedSkipped++;scoreAwareSkipped++;}
	boolean adaptiveOpportunity(){
		final long total=strictCompressedSkipped+compressedEvaluated;
		return total>=64 && strictCompressedSkipped*2<=total;
	}

	void add(final QuantumHybridStats other){
		ordinaryAttempts+=other.ordinaryAttempts;accepted+=other.accepted;
		fallback+=other.fallback;unsupported+=other.unsupported;
		selectedTraceUncertain+=other.selectedTraceUncertain;
		missingTrace+=other.missingTrace;belowThreshold+=other.belowThreshold;
		longTraceFallback+=other.longTraceFallback;
		geometryFallback+=other.geometryFallback;
		compressedEvaluated+=other.compressedEvaluated;
		compressedSkipped+=other.compressedSkipped;
		strictCompressedSkipped+=other.strictCompressedSkipped;
		scoreAwareSkipped+=other.scoreAwareSkipped;
	}

	String toTsv(){
		return "quantum_hybrid\tordinary_attempts\t"+ordinaryAttempts+
			"\taccepted\t"+accepted+"\tfallback\t"+fallback+
			"\tunsupported\t"+unsupported+
			"\tselected_trace_uncertain\t"+selectedTraceUncertain+
			"\tmissing_trace\t"+missingTrace+
			"\tbelow_threshold\t"+belowThreshold+
			"\tlong_trace_fallback\t"+longTraceFallback+
			"\tgeometry_fallback\t"+geometryFallback+
			"\tcompressed_evaluated\t"+compressedEvaluated+
			"\tcompressed_skipped\t"+compressedSkipped+
			"\tstrict_compressed_skipped\t"+strictCompressedSkipped+
			"\tscore_aware_skipped\t"+scoreAwareSkipped+'\n';
	}

	long ordinaryAttempts, accepted, fallback, unsupported;
	long selectedTraceUncertain, missingTrace, belowThreshold, longTraceFallback, geometryFallback;
	long compressedEvaluated, compressedSkipped, strictCompressedSkipped;
	long scoreAwareSkipped;
}
