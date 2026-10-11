package var2;

import ml.CellNet;
import parse.Parse;
import shared.Tools;

/**
 * Mutable thresholds for variant emission and configuration for nearby counting.
 * passesFast checks a subset of the emission thresholds and accepts forced variants;
 * passesFilter applies the full implemented thresholds without a forced bypass.
 * Callers such as VarMap supply that bypass explicitly and calculate coverage first.
 * Nearby counting/flagging is performed by callers, not by this filter.
 *
 * Threshold fields are unsynchronized; configure before concurrent evaluation.
 * Inputs must contain consistent statistics. This class does not validate all
 * ranges or reject every NaN. Scoring delegates to Var and may update cached
 * features or use a mutable CellNet supplied by the caller.
 * @author Brian Bushnell
 * @contributor Isla
 */
public class VarFilter{

	/**
	 * Applies a recognized option in argument order, returning false for unknown keys.
	 * Most numeric options are stored without range validation. minIdentity and
	 * minMaxIdentity values above one are divided by 100 once; maximum mean
	 * identity is divided only when above one and at most 100. Allele fractions and
	 * pairing/strand ratios are taken literally. rarity can lower minAlleleFraction;
	 * lowcoveragepenalty changes global Var state. Assertion failures can leave an
	 * assigned value in place. Parsing does not count neighbors or score a variant.
	 * @param a Non-null, already lowercased option key
	 * @param b Value interpreted by the selected parser
	 * @param arg Original argument, unused
	 * @return Whether the option name was recognized
	 */
	public boolean parse(String a, String b, String arg){
		if(a.equals("minreads") || a.equals("minad") || a.equals("minalleledepth") || a.equals("mincount")){
			minAlleleDepth=Integer.parseInt(b);
		}else if(a.equals("maxreads") || a.equals("maxad") || a.equals("maxalleledepth")){
			maxAlleleDepth=Integer.parseInt(b);
		}else if(a.equals("mincov") || a.equals("mincoverage") || a.equals("mindepth")){
			minCov=Integer.parseInt(b);
		}else if(a.equals("maxcov") || a.equals("maxcoverage") || a.equals("maxdepth")){
			maxCov=Integer.parseInt(b);
		}else if(a.equals("minqualitymax") || a.equals("minmaxquality")){
			minMaxQuality=Integer.parseInt(b);
		}else if(a.equals("minedistmax") || a.equals("minmaxedist")){
			minMaxEdist=Integer.parseInt(b);
		}else if(a.equals("minmapqmax") || a.equals("minmaxmapq")){
			minMaxMapq=Integer.parseInt(b);
		}else if(a.equals("minidmax") || a.equals("minmaxid")){
			minMaxIdentity=Double.parseDouble(b);
			if(minMaxIdentity>1){minMaxIdentity/=100;}
		}

		else if(a.equals("minpairingrate") || a.equals("minpairrate")){
			minPairingRate=Double.parseDouble(b);
		}else if(a.equals("minstrandratio")){
			minStrandRatio=Double.parseDouble(b);
		}else if(a.equals("minscore")){
			minScore=Double.parseDouble(b);
		}else if(a.equals("maxscore")){
			maxScore=Double.parseDouble(b);
		}else if(a.equals("minquality") || a.equals("minavgquality") || a.equals("maq")){
			minAvgQuality=Double.parseDouble(b);
		}else if(a.equals("maxquality") || a.equals("maxavgquality")){
			maxAvgQuality=Double.parseDouble(b);
		}else if(a.equals("minedist") || a.equals("minavgedist") || a.equals("mae")){
			minAvgEdist=Double.parseDouble(b);
		}else if(a.equals("minavgmapq")){
			minAvgMapq=Double.parseDouble(b);
		}else if(a.equals("maxavgmapq")){
			maxAvgMapq=Double.parseDouble(b);
		}else if(a.equals("minallelefraction") || a.equals("minallelefrequency") || a.equals("maf")){
			minAlleleFraction=Double.parseDouble(b);
		}else if(a.equals("maxallelefraction") || a.equals("maxallelefrequency")){
			maxAlleleFraction=Double.parseDouble(b);
		}else if(a.equals("minidentity") || a.equals("mid") || a.equals("minid")){
			minIdentity=Double.parseDouble(b);
			if(minIdentity>1){minIdentity/=100;}
		}else if(a.equals("maxidentity") || a.equals("maxid")){
			maxIdentity=Double.parseDouble(b);
			if(maxIdentity>1 && maxIdentity<=100){maxIdentity/=100;}
		}else if(a.equals("lowcoveragepenalty") || a.equals("lowcovpenalty") || a.equals("covpenalty")){
			Var.lowCoveragePenalty=Double.parseDouble(b);
			assert(Var.lowCoveragePenalty>=0) : "Low coverage penalty must be at least 0.";
		}

		//For adjacent nearby 'rainbow' SNPs
		else if(a.equals("adjacentcount") || a.equals("nearbycount") ||
				a.equals("maxnearbycount") || a.equals("maxnearby")){
			maxNearbyCount=Parse.parseIntKMG(b);
			assert(maxNearbyCount>=1); //1 means at least 1 nearby var
			maxNearbyCount=Tools.max(maxNearbyCount, 1);
		}else if(a.equals("failadjacent") || a.equals("failnearby")){
			failNearby=Parse.parseBoolean(b);
		}else if(a.equals("flagadjacent") || a.equals("flagnearby")){
			flagNearby=Parse.parseBoolean(b);
		}else if(a.equals("penalizeadjacent") || a.equals("penalizenearby")){
			penalizeNearby=Parse.parseBoolean(b);
		}else if(a.equals("adjacentdist") || a.equals("nearbydist") ||
				a.equals("nearbyrange") || a.equals("maxnearbydist")){
			nearbyDist=Parse.parseIntKMG(b);
			assert(nearbyDist>0); //0 means itself
		}else if(a.equals("nearbygap") || a.equals("maxnearbygap") || a.equals("nearbyslop")){
			nearbyGap=Parse.parseIntKMG(b);
			assert(nearbyGap>=0); //0 means on top of each other
		}

		//Bootstrapping: gate which variants COUNT toward a neighbor's NVC, independent of emission filters
		else if(a.equals("nvcmincount") || a.equals("nvcminreads") || a.equals("nvcminad")){
			nvcMinCount=Parse.parseIntKMG(b);
		}else if(a.equals("nvcmaf") || a.equals("nvcminallelefraction")){
			nvcMaf=Double.parseDouble(b);
		}

		else if(a.equals("rarity")){
			rarity=Double.parseDouble(b);
			assert(rarity>=0 && rarity<=1);
			minAlleleFraction=Tools.min(minAlleleFraction, rarity);
		}

		else if(a.equals("clearfilters")){
			if(Parse.parseBoolean(b)){clear();}
		}else{
			return false;
		}
		return true;
	}

	/**
	 * Relaxes emission bounds and disables failNearby, resetting maxNearbyCount to -1
	 * and the two NVC contributor thresholds to permissive values. Retains rarity,
	 * nearbyDist, nearbyGap, flagNearby, penalizeNearby and countNearbyVars. Does not
	 * reset Var's global lowCoveragePenalty or any other global scoring configuration.
	 */
	public void clear(){
		minAlleleDepth=-1;
		maxAlleleDepth=Integer.MAX_VALUE;
		minCov=-1;
		maxCov=Integer.MAX_VALUE;

		minMaxQuality=0;
		minMaxEdist=0;
		minMaxMapq=0;
		minMaxIdentity=0;

		minPairingRate=0;
		minStrandRatio=0;
		minScore=0;
		minAvgQuality=0;
		minAvgEdist=0;
		minAvgMapq=0;
		minAlleleFraction=0;
		minIdentity=0;

		maxScore=Integer.MAX_VALUE;
		maxAvgQuality=Integer.MAX_VALUE;
		maxAvgMapq=Integer.MAX_VALUE;
		maxAlleleFraction=Integer.MAX_VALUE;
		maxIdentity=Integer.MAX_VALUE;

		maxNearbyCount=-1;
		failNearby=false;
		nvcMinCount=-1;
		nvcMaf=0;
	}

	/**
	 * Copies all per-instance thresholds, scoring policy and nearby/counting settings.
	 * Does not copy or change global Var scoring configuration. Self-copy is allowed;
	 * callers must coordinate concurrent mutation of either filter.
	 * @param filter Non-null source; source fields are not modified
	 */
	public void setFrom(VarFilter filter){
		//VAR2-022: CalcTrueQuality transfers parsed options through this method;
		//rarity and nearby/NVC policy must travel with the emission thresholds.
		minAlleleDepth=filter.minAlleleDepth;
		maxAlleleDepth=filter.maxAlleleDepth;
		minCov=filter.minCov;
		maxCov=filter.maxCov;

		minMaxQuality=filter.minMaxQuality;
		minMaxEdist=filter.minMaxEdist;
		minMaxMapq=filter.minMaxMapq;
		minMaxIdentity=filter.minMaxIdentity;

		minPairingRate=filter.minPairingRate;
		minStrandRatio=filter.minStrandRatio;
		minScore=filter.minScore;
		minAvgQuality=filter.minAvgQuality;
		minAvgEdist=filter.minAvgEdist;
		minAvgMapq=filter.minAvgMapq;
		minAlleleFraction=filter.minAlleleFraction;
		minIdentity=filter.minIdentity;

		maxScore=filter.maxScore;
		maxAvgQuality=filter.maxAvgQuality;
		maxAvgMapq=filter.maxAvgMapq;
		maxAlleleFraction=filter.maxAlleleFraction;
		maxIdentity=filter.maxIdentity;
		rarity=filter.rarity;

		maxNearbyCount=filter.maxNearbyCount;
		nearbyDist=filter.nearbyDist;
		nearbyGap=filter.nearbyGap;
		nvcMinCount=filter.nvcMinCount;
		nvcMaf=filter.nvcMaf;
		flagNearby=filter.flagNearby;
		failNearby=filter.failNearby;
		penalizeNearby=filter.penalizeNearby;
		countNearbyVars=filter.countNearbyVars;
	}

	/**
	 * Accepts forced variants immediately; otherwise checks supporting-read count,
	 * maximum base quality, maximum end distance and maximum mapping quality.
	 * Does not check total coverage, maximum identity, averages, AF, nearby count or score.
	 * Passing this method therefore does not imply passing passesFilter.
	 * @param v Non-null variant with accumulated statistics
	 * @return Whether these preliminary checks pass
	 */
	public boolean passesFast(Var v){
		if(v.forced()){return true;}
		final int count=v.alleleCount();
		if(count<minAlleleDepth || count>maxAlleleDepth){return false;}
		if(v.baseQMax<minMaxQuality){return false;}
		if(v.endDistMax<minMaxEdist){return false;}
		if(v.mapQMax<minMaxMapq){return false;}
		return true;
	}

	/**
	 * Tests accumulated statistics against the implemented inclusive bounds.
	 * Unlike passesFast, this method has no forced-variant bypass. It reads stored
	 * coverage without computing it; callers must populate it before filtering.
	 * Nearby rejection runs only when both considerNearby and failNearby are true,
	 * and requires a counted, nonnegative nearbyVarCount under assertions.
	 *
	 * The dataset pairing rate gates the proper-pair-fraction check when positive;
	 * that check compares to minPairingRate, not a dataset-normalized ratio. Mean
	 * identity uses stored per-mille sums. AF limits run only for positive coverage
	 * and use an existing revised-AF cache without computing a new revision here.
	 * The score stage runs only for an enabled score bound and calls Var.phredScore;
	 * a supplied network is used only through that stage. No universal finite-value
	 * check is made, so invalid statistics can escape individual comparisons.
	 * @param v Non-null variant with coverage and supporting-read statistics
	 * @param pairingRate Dataset proper-pair fraction; also used by scoring
	 * @param totalQualityAvg Dataset mean base quality for scoring
	 * @param totalMapqAvg Dataset mean mapping quality for scoring
	 * @param readLengthAvg Dataset mean read length in bases for scoring
	 * @param ploidy Sample ploidy for scoring
	 * @param map Optional scaffold context, interpreted by Var's scoring path
	 * @param net Optional network used when score bounds require scoring
	 * @param considerNearby Whether this call may enforce failNearby
	 * @return Whether all checks reached by this configuration pass
	 */
	public boolean passesFilter(Var v, double pairingRate, double totalQualityAvg, double totalMapqAvg,
			double readLengthAvg, int ploidy, ScafMap map, CellNet net, boolean considerNearby){

		//TIER 1: BASIC DEPTH AND COVERAGE FILTERING
		//These are fast checks that eliminate obviously bad variants early
		final int count=v.alleleCount();
		if(count<minAlleleDepth || count>maxAlleleDepth){return false;}
		final int cov=v.coverage();
		if(cov<minCov || cov>maxCov){return false;}

		//TIER 2: MAXIMUM QUALITY THRESHOLDS
		//Check that at least one supporting read had acceptable quality metrics
		if(v.baseQMax<minMaxQuality){return false;}
		if(v.endDistMax<minMaxEdist){return false;}
		if(v.mapQMax<minMaxMapq){return false;}
		if(v.idMax*0.001f<minMaxIdentity){return false;}

		//TIER 3: PROXIMITY-BASED FILTERING
		//Filter variants in regions with high variant density (potential sequencing errors)
		if(considerNearby && failNearby){
			assert(v.nearbyVarCount>=0) : "Nearby vars were not counted.";
			if(v.nearbyVarCount>maxNearbyCount){return false;}
		}

		//TIER 4: STATISTICAL QUALITY FILTERS
		//These use averaged metrics across all supporting reads
		//Optimization: Use multiplication instead of division for better performance
		//For positive counts and finite values, compare count*threshold with the sum

		if(pairingRate>0 && minPairingRate>0 && count*minPairingRate>v.properPairCount){return false;}
		if(minAvgQuality>0 && count*minAvgQuality>v.baseQSum){return false;}
		if(minAvgEdist>0 && count*minAvgEdist>v.endDistSum){return false;}
		if(minAvgMapq>0 && count*minAvgMapq>v.mapQSum){return false;}
		if(minIdentity>0 && count*minIdentity*1000>v.idSum){return false;}

		//Upper bounds on quality metrics (detect potential systematic biases)
		if(maxAvgQuality<Integer.MAX_VALUE && count*maxAvgQuality<v.baseQSum){return false;}
		if(maxAvgMapq<Integer.MAX_VALUE && count*maxAvgMapq<v.mapQSum){return false;}
		if(maxIdentity<Integer.MAX_VALUE && count*maxIdentity*1000<v.idSum){return false;}

		//TIER 5: STRAND BIAS DETECTION
		//This requires division but is crucial for detecting sequencing artifacts
		if(minStrandRatio>0 && v.strandRatio()<minStrandRatio){return false;}

		//TIER 6: ALLELE FRACTION FILTERING
		//Evaluates variant frequency within the population/sample
		//Uses the existing revised allele fraction if available, without revising it here
		if(minAlleleFraction>0 && v.coverage()>0){
			final double af=v.revisedAlleleFraction==-1 ? v.alleleFraction() : v.revisedAlleleFraction;
			if(af<minAlleleFraction){return false;}
		}
		if(maxAlleleFraction<Integer.MAX_VALUE && v.coverage()>0){
			final double af=v.revisedAlleleFraction==-1 ? v.alleleFraction() : v.revisedAlleleFraction;
			if(af>maxAlleleFraction){return false;}
		}

		//TIER 7: INTEGRATED SCORING
		//phredScore selects composite scoring or the cutoff-scaled network output.
		//The scoring path may update cached variant features or network state.
		if(minScore>0 || maxScore<Integer.MAX_VALUE){
			double phredScore=v.phredScore(pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map, net);
			if(phredScore<minScore || phredScore>maxScore){return false;}
		}

		return true;
	}

	/**
	 * Returns a legacy partial diagnostic summary with the two supplied context values.
	 * Omits several upper bounds, rarity and all nearby/counting options. The final
	 * minAlleleFraction and minIdentity entries currently lack a separator; this is
	 * not a complete or round-trippable configuration serialization.
	 * @param pairingRate Dataset pairing rate printed as context
	 * @param ploidy Sample ploidy printed as context
	 * @return Newly allocated partial summary
	 */
	public String toString(double pairingRate, int ploidy){
		StringBuilder sb=new StringBuilder();

		sb.append("pairingRate=").append(pairingRate).append("\n");
		sb.append("ploidy=").append(ploidy).append("\n");

		sb.append("minReads=").append(minAlleleDepth).append("\n");
		sb.append("maxReads=").append(maxAlleleDepth).append("\n");
		sb.append("minCov=").append(minCov).append("\n");
		sb.append("maxCov=").append(maxCov).append("\n");
		sb.append("minMaxQuality=").append(minMaxQuality).append("\n");
		sb.append("minMaxEdist=").append(minMaxEdist).append("\n");
		sb.append("minMaxMapq=").append(minMaxMapq).append("\n");
		sb.append("minMaxIdentity=").append(minMaxIdentity).append("\n");

		sb.append("minPairingRate=").append(minPairingRate).append("\n");
		sb.append("minStrandRatio=").append(minStrandRatio).append("\n");
		sb.append("minScore=").append(minScore).append("\n");
		sb.append("minAvgQuality=").append(minAvgQuality).append("\n");
		sb.append("minAvgEdist=").append(minAvgEdist).append("\n");
		sb.append("minAvgMapq=").append(minAvgMapq).append("\n");
		//TODO: Probable bug [VAR2-023] - no separator before minIdentity below.
		//Only commented diagnostic calls found in the scoped var2 search; retain output for now.
		sb.append("minAlleleFraction=").append(minAlleleFraction);
		sb.append("minIdentity=").append(minIdentity);

		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	//DEPTH FILTERING PARAMETERS
	/** Inclusive minimum supporting-read count; used by both filter methods. */
	public int minAlleleDepth=2;
	/** Inclusive maximum supporting-read count; used by both filter methods. */
	public int maxAlleleDepth=Integer.MAX_VALUE;
	/** Inclusive minimum total coverage; passesFilter only. */
	public int minCov=-1;
	/** Inclusive maximum total coverage; passesFilter only. */
	public int maxCov=Integer.MAX_VALUE;

	//MAXIMUM QUALITY THRESHOLDS (each maximum must meet its bound)
	/** Minimum maximum base quality among supporting reads. */
	public int minMaxQuality=15;
	/** Minimum maximum end distance in bases among supporting reads. */
	public int minMaxEdist=20;
	/** Minimum maximum mapping quality among supporting reads. */
	public int minMaxMapq=0;
	/** Minimum maximum alignment identity as a fraction; passesFilter only. */
	public double minMaxIdentity=0;

	//STATISTICAL FILTERING PARAMETERS
	/** Minimum supporting-read proper-pair fraction, enabled for positive dataset pairing. */
	public double minPairingRate=0.1;
	/** Minimum Var.strandRatio; enabled when positive. */
	public double minStrandRatio=0.1;
	/** Inclusive minimum Var.phredScore; a positive value enables scoring. */
	public double minScore=20;
	/** Inclusive maximum Var.phredScore; below MAX_VALUE enables scoring. */
	public double maxScore=Integer.MAX_VALUE;
	/** Minimum mean base quality; enabled when positive. */
	public double minAvgQuality=12;
	/** Maximum mean base quality; below MAX_VALUE enables this bound. */
	public double maxAvgQuality=Integer.MAX_VALUE;
	/** Minimum mean end distance in bases; enabled when positive. */
	public double minAvgEdist=10;
	/** Minimum mean mapping quality; enabled when positive. */
	public double minAvgMapq=0;
	/** Maximum mean mapping quality; below MAX_VALUE enables this bound. */
	public double maxAvgMapq=Integer.MAX_VALUE;
	/** Minimum cached-revised or raw AF; applied only with positive coverage and bound. */
	public double minAlleleFraction=0.1;
	/** Maximum cached-revised or raw AF; applied with positive coverage and bound below MAX_VALUE. */
	public double maxAlleleFraction=Integer.MAX_VALUE;
	/** Minimum mean alignment identity as a fraction; compared against per-mille sums. */
	public double minIdentity=0;
	/** Maximum mean alignment identity as a fraction; below MAX_VALUE enables it. */
	public double maxIdentity=Integer.MAX_VALUE;
	/** Var scoring parameter; parsing it may also lower the minimum AF. Not reset by clear. */
	public double rarity=1;

	//PROXIMITY-BASED FILTERING PARAMETERS
	/** Largest accepted nearby count when failNearby is enabled for the call. */
	public int maxNearbyCount=1;
	/** Distance in bases supplied by callers to nearby-variant counting. */
	public int nearbyDist=20;
	/** Gap in bases supplied by callers to nearby-variant counting. */
	public int nearbyGap=2;

	//NVC contributor thresholds are interpreted by VarMap separately from emission thresholds.
	/** Minimum supporting-read count for contributing to NVC; tested by VarMap. */
	public int nvcMinCount=-1;
	/** Minimum AF for contributing to NVC; tested by VarMap when positive. */
	public double nvcMaf=0;

	/** Request flagging during caller-owned nearby counting; not read by passesFilter. */
	public boolean flagNearby=false;
	/** Enable nearby-count rejection when considerNearby is true. */
	public boolean failNearby=false;
	/** Legacy parsed setting; no current Java consumer found beyond configuration. */
	public boolean penalizeNearby=false;
	/** Legacy field; CallVariants/CallVariants2 use their own counting switches. */
	public boolean countNearbyVars=true;
}
