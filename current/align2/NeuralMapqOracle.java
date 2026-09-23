package align2;

import java.util.Arrays;

/** Loose synthetic-truth oracle used to train and grade neural MAPQ. */
public final class NeuralMapqOracle {

	private NeuralMapqOracle(){}

	public static boolean isCorrectLoose(final boolean mapped,
			final String mappedReference, final byte mappedStrand,
			final int mappedStart, final int mappedStop,
			final String truthReference, final byte truthStrand,
			final int truthStart, final int truthStop){
		return isCorrectLoose(mapped, mappedReference, mappedStrand, mappedStart,
				mappedStop, truthReference, truthStrand, truthStart, truthStop,
				NeuralMapqFeatureSchema.ORACLE_TOLERANCE);
	}

	/**
	 * Mirrors MakeRocCurve.isCorrectHitLoose: mapped, reference, and strand must
	 * agree, then either endpoint may fall within the requested tolerance.
	 */
	public static boolean isCorrectLoose(final boolean mapped,
			final String mappedReference, final byte mappedStrand,
			final int mappedStart, final int mappedStop,
			final String truthReference, final byte truthStrand,
			final int truthStart, final int truthStop, final int tolerance){
		if(tolerance<0){
			throw new IllegalArgumentException("Loose MAPQ truth tolerance must be nonnegative: "+tolerance);
		}
		if(!mapped || mappedReference==null || truthReference==null){return false;}
		if(mappedStrand!=truthStrand || !truthReference.equals(mappedReference)){return false;}
		return absdif(mappedStart, truthStart)<=tolerance ||
				absdif(mappedStop, truthStop)<=tolerance;
	}

	/** Allocation-free reference-name variant used by the feature exporter. */
	public static boolean isCorrectLoose(final boolean mapped,
			final byte[] mappedReference, final byte mappedStrand,
			final int mappedStart, final int mappedStop,
			final byte[] truthReference, final byte truthStrand,
			final int truthStart, final int truthStop, final int tolerance){
		if(tolerance<0){
			throw new IllegalArgumentException("Loose MAPQ truth tolerance must be nonnegative: "+tolerance);
		}
		if(!mapped || mappedReference==null || truthReference==null){return false;}
		if(mappedStrand!=truthStrand || !Arrays.equals(truthReference,mappedReference)){return false;}
		return absdif(mappedStart,truthStart)<=tolerance || absdif(mappedStop,truthStop)<=tolerance;
	}

	private static long absdif(final int a, final int b){
		final long difference=(long)a-b;
		return difference<0 ? -difference : difference;
	}
}
