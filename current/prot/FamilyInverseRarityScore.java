package prot;

/** Immutable inverse-posting-rarity rank score bound to an existing caller pivot.
 * The index supplies u=sum(meanPostingDepth/postingDepth) in the raw-hit pass.
 * Linear scoring changes only the hit contribution. Protected sort-grid modes
 * add a post-length bonus, retaining the already validated member medians,
 * combined25 cosine, KD similarity and physical-rank tie break.
 * It neither fits coefficients nor changes admission or acceptance gates.
 * @author Keqing
 */
public final class FamilyInverseRarityScore{

	/** Retain the bound metadata rather than reading or reconstructing its files. */
	public FamilyInverseRarityScore(final FamilyTailCascade.LoadedClassifier base, final double x_, final double y_){
		require(base!=null, "Inverse rarity needs an existing loaded classifier");
		validateCoefficients(x_, y_);
		config=base.config(); heuristic=config.heuristicRanking; mode=Mode.LINEAR;
		x=x_; y=y_; alpha=kappa=saturation=0; protectHead=false; arm=-1;
		validate();
	}

	/** Bind one of the two protected sort-grid finalists to the real caller order. */
	public FamilyInverseRarityScore(final FamilyTailCascade.LoadedClassifier base, final int sortGridArm){
		require(base!=null, "Inverse rarity needs an existing loaded classifier");
		config=base.config(); heuristic=config.heuristicRanking; arm=sortGridArm;
		x=y=0; protectHead=true;
		if(sortGridArm==319){
			mode=Mode.MEDIAN; alpha=4; kappa=1; saturation=0;
		}else if(sortGridArm==448){
			mode=Mode.POSITIVE; alpha=1.5; kappa=0; saturation=0.25;
		}else{throw new IllegalArgumentException("Unsupported sort-grid arm for caller ranking: "+sortGridArm);}
		validate();
	}

	private void validate(){
		require(heuristic!=null && "hitfourth_h020_c150".equals(heuristic.armName()),
			"Inverse rarity replaces hits only in the bound fourth-root, cosine1.5, KD0.2 pivot");
		require(boundTo(config), "Inverse rarity requires the same complete raw uint16 index, metadata and fixed floor5/top50 view");
		final double mean=config.index.meanPostingDepth();
		require(Double.isFinite(mean) && mean>0, "Inverse rarity requires positive finite mean posting depth, not an empty index: "+mean);
	}

	/** Score one raw-admitted family from the existing worker's measured features. */
	public double score(final int queryLength, final int ordinal, final int hits, final double rarity,
			final double cosine, final double queryMean){
		require(mode==Mode.LINEAR, "Sort-grid modes need a complete admitted pool, not one-at-a-time scoring");
		return scorePrepared(x, y, hits, rarity, lengthFactor(queryLength, ordinal), cosine,
			hydrophobicitySimilarity(queryMean, ordinal));
	}

	/** Prepared once per query/family and reusable across the finite coefficient grid. */
	public double lengthFactor(final int queryLength, final int ordinal){
		heuristic.validateQueryLength(queryLength);
		final double median=heuristic.familyMedianAt(ordinal);
		final double ratio=Math.min(queryLength, median)/Math.max(queryLength, median);
		assert(ratio>0 && ratio<=1) : "Validated positive query and member lengths define a ratio in (0,1]";
		return Math.sqrt(Math.sqrt(ratio));
	}

	/** Uses the identical member KD mean and native residue scale as the old pivot. */
	public double hydrophobicitySimilarity(final double queryMean, final int ordinal){
		return FamilyHeuristicRanking.hydrophobicitySimilarity(queryMean, heuristic.familyHydrophobicityAt(ordinal));
	}

	/** Validate even queries rejected before candidate scoring. */
	public void validateQueryLength(final int queryLength){heuristic.validateQueryLength(queryLength);}
	/** Recheck the original metadata binding using its existing verification policy. */
	public void verify(){heuristic.verify();}
	/** The query-median arm needs all admitted rarity ratios before any candidate is sorted. */
	public boolean needsQueryMedian(){return mode==Mode.MEDIAN;}
	/** Protected finalists pin the deployed baseline winner at rank one. */
	public boolean protectsBaselineRankOne(){return protectHead;}
	/** Sort-grid arm number, or -1 for the older linear X/Y grid. */
	public int sortGridArm(){return arm;}
	public double x(){return x;}
	public double y(){return y;}

	/** Score one candidate after this query's complete admitted pool has been observed. */
	public double scorePrepared(final int hits, final double rarity, final double lengthFactor,
			final double cosine, final double similarity, final double queryMedian){
		if(mode==Mode.LINEAR){return scorePrepared(x, y, hits, rarity, lengthFactor, cosine, similarity);}
		final double base=scorePrepared(hits, lengthFactor, cosine, similarity);
		final double ratio=rarityRatio(hits, rarity);
		final double bonus;
		if(mode==Mode.MEDIAN){
			final double centered=(ratio-queryMedian)/(ratio+queryMedian);
			final double evidence=hits/(hits+kappa);
			bonus=alpha*evidence*(centered<-1 ? -1 : centered>1 ? 1 : centered);
		}else{
			bonus=alpha*Math.min(1, ratio/saturation);
		}
		final double score=base+bonus;
		require(Double.isFinite(score), "Sort-grid inverse-rarity score overflowed");
		return score;
	}

	public static double rarityRatio(final int hits, final double rarity){
		if(hits<=0 || hits>1024 || !Double.isFinite(rarity) || rarity<=0){
			throw new IllegalArgumentException("Raw-admitted uint16 candidates need positive raw hits and inverse-rarity sum: hits="+hits+", rarity="+rarity);
		}
		return rarity/hits;
	}

	public static double median(final double[] values, final int size){
		require(values!=null && size>0 && size<=values.length, "Query-median rarity needs a nonempty scratch prefix");
		java.util.Arrays.sort(values, 0, size);
		return values[size/2];
	}

	/** Gate-only copies may retain this binding; feature, index and ranking changes may not. */
	boolean boundTo(final FamilyTailCascade.SearchConfig candidate){
		return candidate!=null && candidate.index==config.index && candidate.index.supportsInversePostingRarity() &&
			candidate.heuristicRanking==heuristic && heuristic!=null && heuristic.boundTo(candidate) &&
			candidate.scalarBounds==config.scalarBounds && candidate.combinedSpectra==config.combinedSpectra &&
			candidate.combinedNorms==config.combinedNorms && candidate.geometry==config.geometry &&
			candidate.width==config.width && candidate.k==config.k && candidate.seed==config.seed &&
			candidate.positionMode==config.positionMode && candidate.partitions==config.partitions && candidate.overlap==config.overlap &&
			candidate.blacklist==config.blacklist && candidate.top==50 && candidate.ddlReject==5 &&
			candidate.ddlAccept>candidate.width && !candidate.union && candidate.truth==null &&
			candidate.lengthRanking==null && candidate.lengthPolicy==null && candidate.memberScalarMeans==null &&
			!candidate.lengthScaledShortAdmission && candidate.familyAdmission==null && candidate.hitFloors==null;
	}

	/** Pure grid arithmetic. Inputs are measured once and remain resident in worker scratch.
	 * X=1,Y=0 uses the same operation order as the existing common pivot. */
	public static double scorePrepared(final double x, final double y, final int hits, final double rarity,
			final double lengthFactor, final double cosine, final double similarity){
		validateCoefficients(x, y);
		if(hits<=0 || hits>1024 || !Double.isFinite(rarity) || rarity<=0){
			throw new IllegalArgumentException("Raw-admitted uint16 candidates need positive raw hits and inverse-rarity sum: hits="+hits+", rarity="+rarity);
		}
		final double weighted=x*hits+y*rarity;
		return scorePrepared(weighted, lengthFactor, cosine, similarity);
	}

	/** Common pivot from an already-combined hit term. */
	public static double scorePrepared(final double weighted, final double lengthFactor,
			final double cosine, final double similarity){
		require(Double.isFinite(weighted) && weighted>0, "Prepared hit contribution must be positive and finite");
		require(Double.isFinite(lengthFactor) && lengthFactor>0 && lengthFactor<=1 &&
			Double.isFinite(cosine) && cosine>=0 && cosine<=1 &&
			Double.isFinite(similarity) && similarity>=0 && similarity<=1,
			"Prepared ranking coordinates must retain positive fourth-root length factor and measured cosine/KD in [0,1]");
		final double score=weighted*lengthFactor+1.5*cosine+0.2*similarity;
		require(Double.isFinite(score), "Inverse-rarity rank score overflowed; no nonfinite ordering is permitted");
		return score;
	}

	private static void validateCoefficients(final double x, final double y){
		if(!Double.isFinite(x) || !Double.isFinite(y) || x<0 || y<0 || (x==0 && y==0)){
			throw new IllegalArgumentException("Inverse-rarity coefficients must be finite and nonnegative, with at least one positive: X="+x+", Y="+y);
		}
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private final FamilyTailCascade.SearchConfig config;
	private final FamilyHeuristicScore heuristic;
	private final Mode mode;
	private final double x, y, alpha, kappa, saturation;
	private final boolean protectHead;
	private final int arm;
	private enum Mode{LINEAR, MEDIAN, POSITIVE}
}
