package align2;

import java.util.Arrays;

/**
 * Loose and exact synthetic-truth predicates used to train and grade neural MAPQ.
 * Callers supply scaffold-relative coordinates using the same convention for
 * truth and mapping. This class does no SAM coordinate conversion or validation
 * of endpoint order. Reference names are compared exactly and case-sensitively.
 *
 * @author Collei
 */
public final class NeuralMapqOracle{

	private NeuralMapqOracle(){}

	/** Uses the raw feature schema's default endpoint tolerance. */
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
	 * The boundary is inclusive. Zero tolerance still accepts just one exact
	 * endpoint; it is not strict correctness (which requires both endpoints).
	 * Negative tolerance is rejected even when the read is unmapped.
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

	/** Allocation-free reference-name variant with the same loose endpoint rule. */
	public static boolean isCorrectLoose(final boolean mapped,
			final byte[] mappedReference, final byte mappedStrand,
			final int mappedStart, final int mappedStop,
			final byte[] truthReference, final byte truthStrand,
			final int truthStart, final int truthStop, final int tolerance){
		if(tolerance<0){
			throw new IllegalArgumentException("Loose MAPQ truth tolerance must be nonnegative: "+tolerance);
		}
		if(!mapped || mappedReference==null || truthReference==null){return false;}
		if(mappedStrand!=truthStrand || !Arrays.equals(truthReference, mappedReference)){return false;}
		return absdif(mappedStart, truthStart)<=tolerance || absdif(mappedStop, truthStop)<=tolerance;
	}

	/**
	 * Exact-placement target used by the strict MAPQ training recipe: reference,
	 * strand and both endpoints must match. Unlike loose tolerance zero, matching
	 * just one endpoint is insufficient. Callers validate coordinate conventions.
	 */
	public static boolean isCorrectStrict(final boolean mapped,
			final String mappedReference, final byte mappedStrand,
			final int mappedStart, final int mappedStop,
			final String truthReference, final byte truthStrand,
			final int truthStart, final int truthStop){
		return mapped && mappedReference!=null && truthReference!=null &&
				mappedStrand==truthStrand && truthReference.equals(mappedReference) &&
				mappedStart==truthStart && mappedStop==truthStop;
	}

	/** Allocation-free exact-placement variant with the same endpoint rule. */
	public static boolean isCorrectStrict(final boolean mapped,
			final byte[] mappedReference, final byte mappedStrand,
			final int mappedStart, final int mappedStop,
			final byte[] truthReference, final byte truthStrand,
			final int truthStart, final int truthStop){
		return mapped && mappedReference!=null && truthReference!=null &&
				mappedStrand==truthStrand && Arrays.equals(truthReference, mappedReference) &&
				mappedStart==truthStart && mappedStop==truthStop;
	}

	/** Widen before subtraction so opposite int extremes do not wrap into a near hit. */
	private static long absdif(final int a, final int b){
		final long difference=(long)a-b;
		return difference<0 ? -difference : difference;
	}
}
