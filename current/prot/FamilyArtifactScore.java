package prot;

/**
 * Frozen authoritative score contract for the MAG-QC family-artifact assay.
 *
 * <p>Flat consensuses use the symmetric, length-adjusted BLOSUM62 ratio from
 * {@link ProteinGroupingScoreSweep}. Profile-HMM scores remain HMMER
 * full-sequence bit scores and are passed to the generic margin helpers without
 * pretending that their numeric scale is interchangeable with the flat score.</p>
 *
 * @author Elly
 */
public final class FamilyArtifactScore {

	/** Linear cost per gap residue for the flat-consensus alignment. */
	public static final int FLAT_GAP_COST=4;
	/** Gentle terminal-length penalty exponent. */
	public static final double FLAT_LENGTH_GAMMA=0.25;
	/** Pre-registered flat-consensus margin in BSR percentage points. */
	public static final double FLAT_MARGIN_MIN=5.0;
	/** Pre-registered HMMER full-sequence margin in bits. */
	public static final double HMM_MARGIN_MIN=5.0;

	private FamilyArtifactScore(){}

	/**
	 * Scores a query against one flat consensus in BSR percentage points [0,100].
	 * A pair with no positive local alignment scores zero.
	 */
	public static double flatScore(final byte[] query, final byte[] consensus){
		if(query==null || consensus==null){throw new IllegalArgumentException("Null protein sequence.");}
		final ProteinGroupingScoreSweep.Metrics m=ProteinGroupingScoreSweep.scorePair(
			query, consensus, FLAT_GAP_COST, FLAT_LENGTH_GAMMA);
		return m==null ? 0 : 100.0*m.adjustedScoreRatio;
	}

	/** Own-family score minus the strongest non-self score, in the scorer's native units. */
	public static double margin(final double ownScore, final double strongestOtherScore){
		if(!Double.isFinite(ownScore) || !Double.isFinite(strongestOtherScore)){
			throw new IllegalArgumentException("Artifact scores must be finite: own="+ownScore+
				", other="+strongestOtherScore);
		}
		return ownScore-strongestOtherScore;
	}

	/**
	 * True only when the own family is the unique top scorer. Deterministic output
	 * may hash-break a tie, but a tied biological assignment is still a failure.
	 */
	public static boolean uniqueOwnTop(final double ownScore, final double strongestOtherScore){
		return margin(ownScore, strongestOtherScore)>0;
	}

	/** True when the native-unit margin clears its pre-registered method-specific bar. */
	public static boolean clearsMargin(final double ownScore, final double strongestOtherScore,
			final double minimumMargin){
		if(!Double.isFinite(minimumMargin) || minimumMargin<0){
			throw new IllegalArgumentException("Margin threshold must be finite and nonnegative: "+minimumMargin);
		}
		return margin(ownScore, strongestOtherScore)>=minimumMargin;
	}
}
