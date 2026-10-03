package prot;

/**
 * Allocation-free acceptance predicate shared by production assignment and
 * consensus construction.
 *
 * <p>The caller supplies the profile-specific thresholds and the already
 * evaluated length decision.  Production uses schema-7 effective thresholds;
 * construction uses its fixed global floors and neutral thresholds for gates
 * that are intentionally disabled there.  Every comparison is inclusive.</p>
 *
 * @author Yoimiya
 */
public final class FamilyAcceptanceGate {

	private FamilyAcceptanceGate(){}
	/** Explicitly disables an integer gate; no other integer has this meaning. */
	public static final int PASS_ALL_INT=Integer.MIN_VALUE;
	/** Explicitly disables a double gate; NaN and infinities are always invalid. */
	public static final double PASS_ALL_DOUBLE=-Double.MAX_VALUE;
	/** Explicitly disables a float gate; NaN and infinities are always invalid. */
	public static final float PASS_ALL_FLOAT=-Float.MAX_VALUE;

	/** Applies the seven gates that precede the conditional HBM traceback score. */
	public static boolean passesBeforeHbm(final boolean enabled,
			final boolean lengthPass, final int rawScore, final double identity,
			final double normalizedR, final double mutualOverlap,
			final int kmerCount, final double kmerQueryDensity,
			final int minRawScore, final double minIdentity,
			final double minNormalizedR, final double minMutualOverlap,
			final int minKmerCount, final double minKmerQueryDensity){
		validate(minRawScore,minIdentity,minNormalizedR,minMutualOverlap,
			minKmerCount,minKmerQueryDensity,PASS_ALL_FLOAT);
		return enabled && lengthPass &&
			(minRawScore==PASS_ALL_INT || rawScore>=minRawScore) &&
			(minIdentity==PASS_ALL_DOUBLE || identity>=minIdentity) &&
			(minNormalizedR==PASS_ALL_DOUBLE || normalizedR>=minNormalizedR) &&
			(minMutualOverlap==PASS_ALL_DOUBLE || mutualOverlap>=minMutualOverlap) &&
			(minKmerCount==PASS_ALL_INT || kmerCount>=minKmerCount) &&
			(minKmerQueryDensity==PASS_ALL_DOUBLE ||
				kmerQueryDensity>=minKmerQueryDensity);
	}

	/** Applies all eight schema-7 gates, including HBM path-relative score. */
	public static boolean passes(final boolean enabled,
			final boolean lengthPass, final int rawScore, final double identity,
			final double normalizedR, final double mutualOverlap,
			final int kmerCount, final double kmerQueryDensity,
			final float hbmPathRelative, final int minRawScore,
			final double minIdentity, final double minNormalizedR,
			final double minMutualOverlap, final int minKmerCount,
			final double minKmerQueryDensity, final float minHbmPathRelative){
		validate(minRawScore,minIdentity,minNormalizedR,minMutualOverlap,
			minKmerCount,minKmerQueryDensity,minHbmPathRelative);
		return passesBeforeHbm(enabled,lengthPass,rawScore,identity,normalizedR,
			mutualOverlap,kmerCount,kmerQueryDensity,minRawScore,minIdentity,
			minNormalizedR,minMutualOverlap,minKmerCount,minKmerQueryDensity) &&
			(minHbmPathRelative==PASS_ALL_FLOAT ||
				hbmPathRelative>=minHbmPathRelative);
	}

	/** Rejects unset, nonfinite, or out-of-domain thresholds before comparison. */
	private static void validate(final int minRawScore, final double minIdentity,
			final double minNormalizedR, final double minMutualOverlap,
			final int minKmerCount, final double minKmerQueryDensity,
			final float minHbmPathRelative){
		if(minRawScore==Integer.MAX_VALUE ||
				(minIdentity!=PASS_ALL_DOUBLE && (!Double.isFinite(minIdentity) ||
					minIdentity<0 || minIdentity>100)) ||
				(minNormalizedR!=PASS_ALL_DOUBLE && !Double.isFinite(minNormalizedR)) ||
				(minMutualOverlap!=PASS_ALL_DOUBLE && (!Double.isFinite(minMutualOverlap) ||
					minMutualOverlap<0 || minMutualOverlap>1)) ||
				(minKmerCount!=PASS_ALL_INT && minKmerCount<0) ||
				(minKmerQueryDensity!=PASS_ALL_DOUBLE &&
					(!Double.isFinite(minKmerQueryDensity) || minKmerQueryDensity<0 ||
					minKmerQueryDensity>1)) ||
				(minHbmPathRelative!=PASS_ALL_FLOAT &&
					(!Float.isFinite(minHbmPathRelative) || minHbmPathRelative<0 ||
					minHbmPathRelative>1))){
			throw new IllegalArgumentException("Invalid or unset family-acceptance threshold");
		}
	}
}
