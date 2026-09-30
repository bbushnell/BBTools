package align2;

/**
 * Routes total reference bases to separately calibrated paired-MAPQ tables.
 * BBMapS passes Data.numBases. Reference sizes between the two supported regimes
 * intentionally retain heuristic MAPQ, rather than extrapolating a calibration.
 * @author Collei
 */
public final class NeuralMapqPairedReferenceScale{
	private NeuralMapqPairedReferenceScale(){}

	/** Inclusive boundaries; missing/nonpositive sizes and the middle interval are unsupported. */
	public static int regime(final long referenceBases){
		if(referenceBases<1){return UNSUPPORTED;}
		if(referenceBases<=SMALL_MAX_BASES){return SMALL;}
		if(referenceBases>=LARGE_MIN_BASES){return LARGE;}
		return UNSUPPORTED;
	}

	public static final int UNSUPPORTED=0, SMALL=1, LARGE=2;
	public static final long SMALL_MAX_BASES=20_000_000L, LARGE_MIN_BASES=1_000_000_000L;
	/** Historical paired-V1 caps; current inference uses NeuralMapqLengthCaps instead. */
	public static final int LARGE_MAPQ_CAP=42, SMALL_MAPQ_CAP=43;
}
