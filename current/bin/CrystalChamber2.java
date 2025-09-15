package bin;

import java.util.ArrayList;
import java.util.Random;

/**
 * Enhanced recrystallization-based bin refinement using iterative centroid clustering.
 * Implements binary splitting (k=2) with farthest-first centroid initialization and assignment stability detection.
 * The algorithm performs K-means-style clustering using Oracle similarity calculations instead of Euclidean distance,
 * making it suitable for biological sequence similarity rather than numerical feature clustering.
 * Uses fixed binary splits to avoid algorithmic complexity since recursive splitting can achieve multi-way partitions.
 * 
 * The recrystallization process iteratively assigns contigs to the nearest centroid (highest similarity),
 * then updates centroids to represent each cluster. Convergence is detected when assignments stabilize
 * between iterations, providing faster termination than traditional centroid movement thresholds.
 * 
 * "Perfect binary splits with better initialization - no need for complexity." - UMP45
 * 
 * @author UMP45
 * @author Brian Bushnell
 */
class CrystalChamber2 extends AbstractRefiner {
	
	/**
	 * Creates enhanced CrystalChamber2 refiner with specified Oracle for similarity calculations.
	 * Initializes recrystallization parameters for iterative centroid-based clustering with binary splitting.
	 * Sets reasonable defaults for convergence (50 iterations max), split quality thresholds (0.1 minimum improvement),
	 * and reproducible random seed for consistent centroid initialization across runs.
	 * Enables debugging output to track split success rates and similarity thresholds.
	 * @param oracle_ Oracle instance for contig similarity evaluation during clustering
	 */
	public CrystalChamber2(Oracle oracle_){
		oracle=oracle_;
		maxIterations=50; //Prevent infinite clustering loops
		convergenceThreshold=0.01f; //Traditional centroid movement threshold (unused in stability detection)
		minSplitImprovement=0.1f; //Minimum similarity difference to justify cluster separation
		random=new Random(12345); //Reproducible results for testing
		debug=true; //Enable debugging output for split analysis
		splitAttempts=0; //Initialize debugging counters
		successfulSplits=0;
	}
	
	/**
	 * Performs binary splitting refinement on input cluster using recrystallization algorithm.
	 * Attempts to separate the input cluster into two sub-clusters using iterative centroid-based clustering.
	 * First validates input size (minimum 4 contigs required for meaningful split), then performs recrystallization
	 * with k=2 binary splitting. Validates split quality using Oracle similarity calculations and split improvement thresholds.
	 * Returns null if input is too small, recrystallization fails, split quality is insufficient, or resulting clusters
	 * would be merged back together based on similarity scores.
	 * @param input Input cluster to refine through binary splitting
	 * @return ArrayList containing two refined clusters, or null if splitting failed or was not beneficial
	 */
	@Override
	ArrayList<Bin> refine(Bin input){
		if(input==null || input.numContigs()<4){return null;} //Require minimum size for binary splitting
		if(!input.isCluster()){return null;} //Only process cluster types, not bins
		
		Cluster cluster=(Cluster)input;
		ArrayList<Contig> contigs=new ArrayList<>(cluster.contigs); //Work with defensive copy
		
		if(contigs.size()<4){return null;} //Double-check minimum size after cast
		
		splitAttempts++; //Track total split attempts for debugging statistics
		
		int k=2; //Fixed binary splitting - recursive splits handle multi-way partitions
		
		ArrayList<Cluster> crystals=recrystallize(contigs, k); //Perform iterative clustering
		
		if(crystals==null || crystals.size()!=2){
			if(debug && splitAttempts%100==0) System.err.println("CrystalChamber2: Recrystallization failed, attempts="+splitAttempts);
			return null;
		}
		
		ArrayList<Bin> result=new ArrayList<Bin>(crystals);
		if(!isSplitBeneficial(input, result)){ //Validate split improves overall quality
			if(debug && splitAttempts%100==0) System.err.println("CrystalChamber2: Split not beneficial, attempts="+splitAttempts);
			return null;
		}
		
		float similarity=oracle.similarity(crystals.get(0), crystals.get(1), 1.0f); //Calculate inter-cluster similarity
		boolean wouldMerge=similarity>minSplitImprovement; //Check if clusters are too similar to justify split
		
		if(debug && splitAttempts%100==0){
			System.err.println("CrystalChamber2 DEBUG: similarity="+similarity+" threshold="+minSplitImprovement+" wouldMerge="+wouldMerge+" attempts="+splitAttempts+" successes="+successfulSplits);
		}
		
		if(wouldMerge){
			return null; //Reject split if clusters would be merged back together
		}
		
		successfulSplits++; //Track successful splits for debugging statistics
		
		cluster.clear(); //Clear original cluster to prevent stale references and memory leaks
		
		return result;
	}
	
	/**
	 * Determines optimal number of clusters based on cluster size heuristics.
	 * Uses simple rules to balance clustering effectiveness with computational complexity.
	 * Small clusters (size < 8) use binary splits to avoid over-fragmentation.
	 * Medium clusters (8-19) can support 3-way splits for better separation.
	 * Large clusters (20-49) use 4-way splits to handle complexity without excessive fragmentation.
	 * Very large clusters use up to 5-way splits or size/10, whichever is smaller.
	 * Note: Current implementation uses fixed k=2, making this method unused but preserved for future flexibility.
	 * @param contigs List of contigs to analyze for optimal k value determination
	 * @return Recommended number of clusters (2-5) based on input size heuristics
	 */
	private int determineOptimalK(ArrayList<Contig> contigs){
		int size=contigs.size();
		
		if(size<8) return 2; //Small clusters: binary split to avoid over-fragmentation
		if(size<20) return 3; //Medium clusters: 3-way split provides better separation
		if(size<50) return 4; //Large clusters: 4-way split balances complexity and effectiveness
		return Math.min(5, size/10); //Very large: cap at 5-way or 10% of size to prevent excessive fragmentation
	}
	
	/**
	 * Performs iterative centroid-based clustering to separate contigs into k clusters using Oracle similarity.
	 * Implements K-means-style algorithm adapted for biological sequence similarity rather than Euclidean distance.
	 * Uses farthest-first centroid initialization to maximize initial separation, then iteratively assigns
	 * contigs to nearest centroids and updates centroid representatives. Enhanced convergence detection uses
	 * assignment stability (no contig reassignments) rather than traditional centroid movement thresholds.
	 * Fails if any cluster becomes empty during iterations or if initialization cannot find sufficient diversity.
	 * Creates proper Cluster objects with correct contig.cluster pointer updates to maintain data integrity.
	 * @param contigs List of contigs to partition into clusters using similarity-based assignments
	 * @param k Number of clusters to create (typically k=2 for binary splitting)
	 * @return ArrayList of k clusters with proper contig assignments, or null if clustering fails
	 */
	private ArrayList<Cluster> recrystallize(ArrayList<Contig> contigs, int k){
		if(contigs.size()<k){return null;} //Cannot create more clusters than contigs available
		
		ArrayList<Centroid> centroids=initializeCentroids(contigs, k); //Farthest-first initialization for maximum separation
		if(centroids==null){return null;}
		
		ArrayList<ArrayList<Contig>> assignments=new ArrayList<>(k);
		for(int i=0; i<k; i++){
			assignments.add(new ArrayList<Contig>());
		}
		
		ArrayList<ArrayList<Contig>> previousAssignments=null; //Store previous iteration for stability detection
		for(int iter=0; iter<maxIterations; iter++){
			for(ArrayList<Contig> list : assignments){
				list.clear(); //Clear previous iteration assignments
			}
			
			for(Contig contig : contigs){
				int bestCentroid=findNearestCentroid(contig, centroids); //Find centroid with highest similarity
				assignments.get(bestCentroid).add(contig);
			}
			
			boolean hasEmpty=false;
			for(ArrayList<Contig> list : assignments){
				if(list.isEmpty()){
					hasEmpty=true;
					break;
				}
			}
			if(hasEmpty){return null;} //Initialization failure - some centroids attract no contigs
			
			if(previousAssignments!=null && assignmentsEqual(assignments, previousAssignments)){
				break; //Converged - assignments stabilized between iterations
			}
			
			for(int i=0; i<k; i++){
				Centroid newCentroid=calculateCentroid(assignments.get(i)); //Update centroid to represent current cluster
				centroids.set(i, newCentroid);
			}
			
			previousAssignments=deepCopyAssignments(assignments); //Store for next iteration comparison
		}
		
		ArrayList<Cluster> result=new ArrayList<>(k); //Convert assignments to proper Cluster objects
		
		for(int i=0; i<k; i++){
			if(assignments.get(i).isEmpty()){return null;} //Final check for empty clusters
			
			Contig firstContig=assignments.get(i).get(0); //Use first contig's ID for cluster identification
			Cluster cluster=new Cluster(firstContig.id());
			
			if(debug && splitAttempts%500==0){
				System.err.println("CrystalChamber2 DEBUG: Creating cluster "+cluster.id()+" with "+assignments.get(i).size()+" contigs");
			}
			
			for(Contig contig : assignments.get(i)){
				cluster.add(contig); //Add contig and update its cluster pointer
				if(contig.cluster!=cluster){ //Verify pointer integrity for debugging
					System.err.println("ERROR: Contig "+contig.id()+" cluster pointer not updated! Points to "+
						(contig.cluster==null ? "null" : contig.cluster.id())+" but should point to "+cluster.id());
				}
			}
			result.add(cluster);
		}
		
		return result;
	}
	
	/**
	 * Initializes centroids using farthest-first strategy to maximize initial cluster separation.
	 * Selects first centroid randomly from available contigs, then iteratively chooses subsequent centroids
	 * to maximize the minimum distance (1 - similarity) to all previously selected centroids.
	 * This approach provides better initial separation than random selection and avoids the problematic
	 * K-means++ probabilistic weighting which doesn't work well with Oracle similarity functions.
	 * Each iteration examines all remaining contigs and selects the one most dissimilar to the closest
	 * existing centroid, creating well-separated initial cluster centers for effective convergence.
	 * @param contigs Available contigs for centroid selection and initialization
	 * @param k Number of centroids to initialize (must be <= contigs.size())
	 * @return ArrayList of k initialized centroids with maximum separation, or null if insufficient contigs
	 */
	private ArrayList<Centroid> initializeCentroids(ArrayList<Contig> contigs, int k){
		if(contigs.size()<k){return null;} //Cannot select more centroids than available contigs
		
		ArrayList<Centroid> centroids=new ArrayList<>(k);
		ArrayList<Contig> chosen=new ArrayList<>(k); //Track selected centroids to avoid duplicates
		
		Contig first=contigs.get(random.nextInt(contigs.size())); //Random first centroid for reproducible diversity
		chosen.add(first);
		centroids.add(new Centroid(first));
		
		for(int i=1; i<k; i++){ //Select remaining centroids using farthest-first heuristic
			Contig best=null;
			float maxMinDistance=-1;
			
			for(Contig candidate : contigs){
				if(chosen.contains(candidate)){continue;} //Skip already chosen centroids
				
				float minDistance=Float.MAX_VALUE; //Find minimum distance to any existing centroid
				for(Contig existing : chosen){
					float similarity=oracle.similarity(candidate, existing, 1.0f);
					float distance=1.0f-similarity; //Convert Oracle similarity to distance metric
					minDistance=Math.min(minDistance, distance);
				}
				
				if(minDistance>maxMinDistance){ //Select candidate that maximizes minimum distance
					maxMinDistance=minDistance;
					best=candidate;
				}
			}
			
			if(best==null){return null;} //Fallback if no suitable candidate found
			chosen.add(best);
			centroids.add(new Centroid(best));
		}
		
		return centroids;
	}
	
	/**
	 * Finds the centroid with highest similarity to the given contig for cluster assignment.
	 * Iterates through all available centroids and calculates Oracle similarity scores,
	 * returning the index of the centroid with maximum similarity. Uses "nearest" in the similarity sense
	 * (highest similarity score) rather than traditional distance minimization, since Oracle provides
	 * similarity rather than distance metrics. This assignment strategy ensures contigs are placed
	 * in clusters where they are most similar to the representative centroid.
	 * @param contig Contig to assign to the most similar centroid
	 * @param centroids Available centroids for similarity-based assignment
	 * @return Index of centroid with highest similarity to the input contig
	 */
	private int findNearestCentroid(Contig contig, ArrayList<Centroid> centroids){
		int best=0;
		float bestSimilarity=-1;
		
		for(int i=0; i<centroids.size(); i++){
			float similarity=centroids.get(i).similarityTo(contig, oracle); //Calculate similarity using Oracle
			if(similarity>bestSimilarity){ //Higher similarity = better assignment
				bestSimilarity=similarity;
				best=i;
			}
		}
		
		return best;
	}
	
	/**
	 * Calculates centroid representative for a group of contigs using size-based selection.
	 * Uses the largest contig in the group as the centroid representative, under the assumption
	 * that larger contigs provide more stable similarity calculations and better represent
	 * the cluster's characteristics. This is a practical heuristic since true feature averaging
	 * is complex for biological sequences and may not improve clustering effectiveness.
	 * For single-contig groups, returns the contig itself as the centroid representative.
	 * The size-based selection provides deterministic centroid updates and avoids computational
	 * complexity of feature space averaging while maintaining reasonable clustering behavior.
	 * @param contigs Group of contigs requiring centroid representation for clustering
	 * @return Centroid object with largest contig as representative, or null if group is empty
	 */
	private Centroid calculateCentroid(ArrayList<Contig> contigs){
		if(contigs.isEmpty()){return null;} //Cannot calculate centroid for empty cluster
		if(contigs.size()==1){return new Centroid(contigs.get(0));} //Single contig is its own centroid
		
		Contig largest=contigs.get(0); //Use largest contig as most representative
		for(Contig c : contigs){
			if(c.size()>largest.size()){largest=c;} //Size comparison for stability
		}
		
		return new Centroid(largest);
		//TODO: Could implement proper averaging of features for more sophisticated centroid calculation
	}
	
	/**
	 * Checks if two assignment sets are identical for convergence detection in clustering iterations.
	 * Performs deep comparison of all contig assignments across all clusters to determine if
	 * assignments have stabilized between iterations. This provides more reliable convergence
	 * detection than traditional centroid movement thresholds, since it directly measures
	 * the clustering objective (stable assignments) rather than intermediate values.
	 * Uses element-wise comparison assuming assignment lists maintain consistent ordering.
	 * @param a First assignment set from current iteration
	 * @param b Second assignment set from previous iteration  
	 * @return true if all assignments are identical (converged), false if any differences exist
	 */
	private boolean assignmentsEqual(ArrayList<ArrayList<Contig>> a, ArrayList<ArrayList<Contig>> b){
		if(a.size()!=b.size()) return false; //Different number of clusters
		
		for(int i=0; i<a.size(); i++){
			if(a.get(i).size()!=b.get(i).size()) return false; //Different cluster sizes
			for(int j=0; j<a.get(i).size(); j++){
				if(!a.get(i).get(j).equals(b.get(i).get(j))) return false; //Different contig assignments
			}
		}
		return true; //All assignments identical - clustering has converged
	}
	
	/**
	 * Creates deep copy of assignment sets for comparison in subsequent clustering iterations.
	 * Prevents accidental modification of stored assignments during convergence checking by creating
	 * independent copies of all cluster assignment lists. Essential for assignment stability detection
	 * since the original assignments will be cleared and rebuilt in each iteration.
	 * Only copies list structure and references, not the Contig objects themselves, which is sufficient
	 * for assignment comparison purposes and avoids unnecessary object duplication.
	 * @param original Original assignment set from current iteration to preserve
	 * @return Deep copy of assignment structure with independent ArrayList instances
	 */
	private ArrayList<ArrayList<Contig>> deepCopyAssignments(ArrayList<ArrayList<Contig>> original){
		ArrayList<ArrayList<Contig>> copy=new ArrayList<>();
		for(ArrayList<Contig> list : original){
			copy.add(new ArrayList<>(list)); //Create new ArrayList with same Contig references
		}
		return copy;
	}
	
	/**
	 * Represents a cluster centroid using a single representative contig for similarity calculations.
	 * Provides a lightweight abstraction for centroid-based clustering where centroids are actual
	 * contigs rather than computed feature vectors. This approach works well with Oracle similarity
	 * functions that operate on pairs of biological sequences, avoiding the complexity of
	 * feature space averaging while maintaining effective clustering behavior.
	 */
	private static class Centroid {
		/** Representative contig serving as the centroid for similarity calculations */
		final Contig representative;
		
		/**
		 * Creates centroid with specified representative contig as the cluster center.
		 * The representative contig will be used for all similarity calculations with this centroid.
		 * @param rep Representative contig that defines this centroid's position in similarity space
		 */
		Centroid(Contig rep){representative=rep;}
		
		/**
		 * Calculates Oracle similarity between this centroid and a target contig for assignment.
		 * Uses the representative contig as the centroid's position in similarity space,
		 * delegating to the Oracle's similarity function for biological sequence comparison.
		 * Higher similarity scores indicate better assignment matches for clustering.
		 * @param contig Target contig for similarity-based distance calculation
		 * @param oracle Oracle instance providing biological sequence similarity calculations
		 * @return Similarity score between centroid representative and target contig (higher = more similar)
		 */
		float similarityTo(Contig contig, Oracle oracle){
			return oracle.similarity(representative, contig, 1.0f);
		}
	}
	
	/** Oracle instance for biological sequence similarity calculations during clustering and validation */
	private final Oracle oracle;
	/** Maximum iterations allowed for centroid convergence before terminating clustering attempts */
	private final int maxIterations;
	/** Traditional convergence threshold for centroid movement (unused in assignment stability detection) */
	private final float convergenceThreshold;
	/** Minimum similarity difference required between clusters to justify split (prevents over-fragmentation) */
	private final float minSplitImprovement;
	/** Random number generator with fixed seed for reproducible centroid initialization across runs */
	private final Random random;
	
	/** Enable debugging output for split success rate analysis and similarity threshold monitoring */
	private final boolean debug;
	/** Total number of split attempts across all inputs for debugging and performance statistics */
	private int splitAttempts;
	/** Number of successful splits that passed quality thresholds for debugging and effectiveness analysis */
	private int successfulSplits;
	
	@Override
	ArrayList<structures.IntHashSet> refineToIntSets(Bin input) {
		ArrayList<Bin> refined = refine(input);
		if(refined == null) return null;
		
		ArrayList<structures.IntHashSet> result = new ArrayList<>();
		for(Bin bin : refined) {
			if(bin.isCluster()) {
				Cluster cluster = (Cluster) bin;
				structures.IntHashSet intSet = new structures.IntHashSet();
				for(Contig contig : cluster.contigs) {
					intSet.add(contig.id());
				}
				result.add(intSet);
			}
		}
		return result.isEmpty() ? null : result;
	}
}