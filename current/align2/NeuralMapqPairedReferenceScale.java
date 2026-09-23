package align2;

/** Conservative paired-V1 reference-size routing for separately calibrated LUTs. */
public final class NeuralMapqPairedReferenceScale {
	private NeuralMapqPairedReferenceScale(){}
	public static int regime(final long referenceBases){
		if(referenceBases<1){return UNSUPPORTED;}
		if(referenceBases<=SMALL_MAX_BASES){return SMALL;}
		if(referenceBases>=LARGE_MIN_BASES){return LARGE;}
		return UNSUPPORTED;
	}
	public static final int UNSUPPORTED=0,SMALL=1,LARGE=2;
	public static final long SMALL_MAX_BASES=20_000_000L,LARGE_MIN_BASES=1_000_000_000L;
	public static final int LARGE_MAPQ_CAP=42,SMALL_MAPQ_CAP=43;
}
