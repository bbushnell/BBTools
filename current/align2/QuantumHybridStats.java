package align2;

/**
 * Worker-local site-routing counters for the broad Quantum/MSA hybrid.
 * Counts scoreSlow decisions, not reads or final accepted mappings. Geometry
 * fallbacks do not attempt Quantum, so accepted+fallback is not ordinaryAttempts.
 * Compressed-site decisions have separate counters. These mutable counters also
 * drive the worker's adaptive policy; they are not merely end-of-run diagnostics.
 * Merge completed workers into a separate total; updates are unsynchronized.
 * @author Collei
 */
public class QuantumHybridStats{

	/** Records a non-compressed site eligible for a Quantum attempt. */
	void ordinaryAttempt(){ordinaryAttempts++;}
	/** Records a Quantum trace accepted for reuse in place of this MSA fill. */
	void accepted(){accepted++;}
	void fallbackUnsupported(){fallback++; unsupported++;}
	/** The caller found N in the selected trace; this is not Result.uncertain. */
	void fallbackUncertain(){fallback++; selectedTraceUncertain++;}
	void fallbackMissing(){fallback++; missingTrace++;}
	void fallbackThreshold(){fallback++; belowThreshold++;}
	/** More than one total inserted/deleted base, not necessarily one long gap. */
	void longTraceFallback(){fallback++; longTraceFallback++;}
	/** Records a site excluded by projected geometry before attempting Quantum. */
	void geometryFallback(){fallback++; geometryFallback++;}
	void compressedEvaluated(){compressedEvaluated++;}
	void compressedSkipped(){compressedSkipped++; strictCompressedSkipped++;}
	void scoreAwareSkipped(){compressedSkipped++; scoreAwareSkipped++;}
	/** Enables the relaxed opportunity test after 64 strict-skipped/evaluated sites
	 * when at most half were strict-skipped. Score-aware skips are excluded from
	 * this denominator. BBMapThread still applies its per-site score/hit conditions. */
	boolean adaptiveOpportunity(){
		final long total=strictCompressedSkipped+compressedEvaluated;
		return total>=64 && strictCompressedSkipped*2<=total;
	}

	/** Adds a completed worker's counters without clearing it; calling twice double-counts. */
	void add(final QuantumHybridStats other){
		ordinaryAttempts+=other.ordinaryAttempts; accepted+=other.accepted;
		fallback+=other.fallback; unsupported+=other.unsupported;
		selectedTraceUncertain+=other.selectedTraceUncertain;
		missingTrace+=other.missingTrace; belowThreshold+=other.belowThreshold;
		longTraceFallback+=other.longTraceFallback;
		geometryFallback+=other.geometryFallback;
		compressedEvaluated+=other.compressedEvaluated;
		compressedSkipped+=other.compressedSkipped;
		strictCompressedSkipped+=other.strictCompressedSkipped;
		scoreAwareSkipped+=other.scoreAwareSkipped;
	}

	/** Emits one newline-terminated record with alternating key/value columns. */
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
