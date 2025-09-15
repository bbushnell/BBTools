package var2;

import ml.CellNet;
import shared.Parse;
import shared.Tools;

/**
* Comprehensive filtering system for genetic variants based on multiple quality metrics.
* This is the core decision-making class that determines which variants are genuine
* versus sequencing artifacts or low-confidence calls.
* 
* Implements a multi-tier filtering approach:
* 1. Fast preliminary filters (read depth, max quality scores)
* 2. Comprehensive statistical analysis (pairing rates, strand bias, coverage)
* 3. Complex scoring algorithms (phred scores, allele fractions)
* 4. Proximity-based filtering (nearby variant detection)
* 
* The passesFilter() method is the primary integration point where all statistical
* evidence is evaluated. This is also the natural insertion point for machine learning
* models that could supplement or replace traditional statistical filters.
* 
* @author Brian Bushnell
* @contributor Isla
*/
public class VarFilter {
	
	/**
	* Parses command-line arguments to configure filtering parameters.
	* Handles a comprehensive set of filtering options with intelligent defaults
	* and automatic unit conversion (e.g., percentages to fractions).
	* 
	* @param a Argument key (lowercase)
	* @param b Argument value
	* @param arg Original argument string
	* @return true if argument was recognized and parsed
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
	 * Resets all filtering parameters to permissive defaults.
	 * Useful for starting with a clean slate or disabling all filters.
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
	}
	
	/**
	 * Copies all filtering parameters from another VarFilter.
	 * Useful for creating filtered copies or applying consistent parameters.
	 * @param filter Source VarFilter to copy parameters from
	 */
	public void setFrom(VarFilter filter){
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
	}
	
	/**
	 * Fast preliminary filtering using only simple numeric thresholds.
	 * This method provides rapid screening of variants before more expensive
	 * statistical calculations. Used for initial filtering in high-throughput scenarios.
	 * 
	 * Checks basic quality metrics that don't require division or complex calculations:
	 * - Read depth (allele count)
	 * - Maximum base quality observed
	 * - Maximum end distance (position within reads)
	 * - Maximum mapping quality
	 * 
	 * @param v Variant to evaluate
	 * @return true if variant passes preliminary filters
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
	* Comprehensive variant filtering using all available statistical evidence.
	* 
	* This is the core filtering method that integrates multiple lines of evidence
	* to determine variant quality. The method implements a layered filtering approach:
	* 
	* 1. DEPTH FILTERS: Checks read count and coverage depth
	* 2. QUALITY FILTERS: Evaluates base quality, mapping quality, alignment identity
	* 3. PROXIMITY FILTERS: Considers nearby variant density (potential sequencing errors)
	* 4. STATISTICAL FILTERS: Analyzes pairing rates, strand bias, average qualities
	* 5. ALLELE FRACTION FILTERS: Evaluates variant frequency in the population
	* 6. INTEGRATED SCORING: Uses phred-scaled composite scores
	* 
	* The method uses several optimization strategies:
	* - Early returns for forced variants (user-specified high-confidence calls)
	* - Fast rejection based on simple thresholds before expensive calculations
	* - Multiplication-based comparisons instead of division (count*threshold > sum)
	* - Conditional evaluation of expensive metrics only when thresholds are set
	* 
	* NEURAL NETWORK INTEGRATION POINT:
	* This method represents the natural insertion point for machine learning models.
	* A neural network could either:
	* 1. Replace this entire method with learned decision boundaries
	* 2. Supplement the statistical filters with additional evidence
	* 3. Provide a final confidence score alongside traditional filtering
	* 
	* The method already calculates all the statistical features that would be
	* useful for ML training: coverage, quality scores, strand bias, proximity metrics, etc.
	* 
	* @param v Variant to evaluate
	* @param pairingRate Overall proper pairing rate from sequencing run
	* @param totalQualityAvg Average base quality from the entire dataset
	* @param totalMapqAvg Average mapping quality from the entire dataset  
	* @param readLengthAvg Average read length from the sequencing run
	* @param ploidy Sample ploidy (typically 1 or 2)
	* @param map Scaffold mapping for coordinate-based calculations
	* @param considerNearby Whether to apply proximity-based filtering
	* @return true if variant passes all filtering criteria
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
		//(count*threshold>sum) is equivalent to (sum/count<threshold) but faster
		
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
		//Uses revised allele fraction if available (corrects for insertion length bias)
		if(minAlleleFraction>0 && v.coverage()>0){
			final double af=v.revisedAlleleFraction==-1 ? v.alleleFraction() : v.revisedAlleleFraction;
			if(af<minAlleleFraction){return false;}
		}
		if(maxAlleleFraction<Integer.MAX_VALUE && v.coverage()>0){
			final double af=v.revisedAlleleFraction==-1 ? v.alleleFraction() : v.revisedAlleleFraction;
			if(af>maxAlleleFraction){return false;}
		}
		
		//TIER 7: INTEGRATED SCORING
		//This is the most expensive calculation - a composite phred score that integrates
		//multiple lines of evidence using sophisticated statistical models
		if(minScore>0 || maxScore<Integer.MAX_VALUE){
			double phredScore=v.phredScore(pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map, null);
			if(net!=null){
				if(phredScore<minScore*0.5f){return false;}
				float[] vec=FeatureVectorMaker.toVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
				float output=net.applyInput(vec).feedForward();
				output=Tools.mid(0, output, 1);
				//Convert NN output to phred score multiplier: cutoff->1.0, 0.0->0.1, 1.0->1.9
				float multiplier=(output<=net.cutoff) ?
					0.1f+(0.9f*output/net.cutoff) :
					1.0f+(0.9f*(output-net.cutoff)/(1.0f-net.cutoff));
				
				phredScore*=multiplier;
			}
			if(phredScore<minScore || phredScore>maxScore){return false;}
		}
		
		return true;
	}
	
	/**
	 * Returns a formatted string representation of all filter parameters.
	 * Useful for logging, debugging, and reproducibility documentation.
	 * @param pairingRate Overall pairing rate for context
	 * @param ploidy Sample ploidy for context  
	 * @return Formatted parameter summary
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
		sb.append("minAlleleFraction=").append(minAlleleFraction);
		sb.append("minIdentity=").append(minIdentity);
		
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	//DEPTH FILTERING PARAMETERS
	/** Minimum number of reads supporting the variant allele */
	public int minAlleleDepth=2;
	/** Maximum number of reads supporting the variant allele */
	public int maxAlleleDepth=Integer.MAX_VALUE;
	/** Minimum total coverage at the variant position */
	public int minCov=-1;
	/** Maximum total coverage at the variant position */
	public int maxCov=Integer.MAX_VALUE;
	
	//MAXIMUM QUALITY THRESHOLDS (at least one read must exceed these)
	/** Minimum base quality score observed among supporting reads */
	public int minMaxQuality=15;
	/** Minimum end distance (position within read) observed among supporting reads */
	public int minMaxEdist=20;
	/** Minimum mapping quality observed among supporting reads */
	public int minMaxMapq=0;
	/** Minimum alignment identity observed among supporting reads */
	public double minMaxIdentity=0;
	
	//STATISTICAL FILTERING PARAMETERS
	/** Minimum proper pairing rate among supporting reads */
	public double minPairingRate=0.1;
	/** Minimum strand ratio (balance between + and - strands) */
	public double minStrandRatio=0.1;
	/** Minimum composite phred score for the variant */
	public double minScore=20;
	/** Maximum composite phred score for the variant */
	public double maxScore=Integer.MAX_VALUE;
	/** Minimum average base quality among supporting reads */
	public double minAvgQuality=12;
	/** Maximum average base quality among supporting reads */
	public double maxAvgQuality=Integer.MAX_VALUE;
	/** Minimum average end distance among supporting reads */
	public double minAvgEdist=10;
	/** Minimum average mapping quality among supporting reads */
	public double minAvgMapq=0;
	/** Maximum average mapping quality among supporting reads */
	public double maxAvgMapq=Integer.MAX_VALUE;
	/** Minimum allele fraction (variant frequency in the sample) */
	public double minAlleleFraction=0.1;
	/** Maximum allele fraction (variant frequency in the sample) */
	public double maxAlleleFraction=Integer.MAX_VALUE;
	/** Minimum average alignment identity among supporting reads */
	public double minIdentity=0;
	/** Maximum average alignment identity among supporting reads */
	public double maxIdentity=Integer.MAX_VALUE;
	/** Expected rarity of variants in the population (affects scoring) */
	public double rarity=1;
	
	//PROXIMITY-BASED FILTERING PARAMETERS
	/** Maximum number of nearby variants allowed before flagging/failing */
	public int maxNearbyCount=1;
	/** Distance threshold for considering variants "nearby" */
	public int nearbyDist=20;
	/** Minimum gap between variants to avoid proximity penalties */
	public int nearbyGap=2;
	
	/** Whether to flag variants with too many nearby variants */
	public boolean flagNearby=false;
	/** Whether to fail/reject variants with too many nearby variants */
	public boolean failNearby=false;
	/** Whether to penalize scores of variants with nearby variants */
	public boolean penalizeNearby=false;
	/** Whether to count nearby variants (enables proximity-based filtering) */
	public boolean countNearbyVars=true;
}