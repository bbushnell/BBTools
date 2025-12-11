package bin;

import java.util.ArrayList;
import java.util.Random;

class CrystalChamber2 extends AbstractRefiner {
	
	/**
	 * Creates a CrystalChamber2 refiner with the specified Oracle for contig similarity calculations.
	 * Initializes recrystallization parameters for binary clustering, including convergence limits, split thresholds, and random seed.
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
	 * Performs binary splitting refinement on an input bin using a recrystallization algorithm.
	 * Uses centroid-based clustering with k=2 and Oracle similarity to decide whether a split is beneficial.
	 * @param input Input bin to refine through binary splitting
	 * @return List containing two refined bins, or null if splitting fails or is not beneficial
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
	 * Determines an optimal number of clusters based on simple size heuristics.
	 * Currently unused because the implementation always performs binary splits (k=2), but preserved for future flexibility.
	 * @param contigs List of contigs to analyze for potential k selection
	 * @return Recommended number of clusters based on heuristic rules
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
	 * Uses farthest-first centroid initialization, repeated assignment/centroid updates, and aborts if any cluster becomes empty.
	 * @param contigs List of contigs to partition into clusters using similarity-based assignments
	 * @param k Number of clusters to create (typically 2 for binary splitting)
	 * @return List of k clusters with contigs assigned, or null if clustering fails
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
	 * Initializes centroids using a farthest-first strategy to maximize initial cluster separation.
	 * Selects the first centroid at random, then iteratively chooses contigs that are far from existing centroids under the Oracle metric.
	 * @param contigs Available contigs for centroid selection and initialization
	 * @param k Number of centroids to initialize (must be <= contigs.size())
	 * @return List of k initialized centroids, or null if insufficient contigs
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
	 * Iterates through all centroids, scores similarity via the Oracle, and returns the best index.
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
	 * Calculates a centroid representative for a group of contigs using size-based selection.
	 * Uses the largest contig in the group as the representative to provide stable similarity calculations.
	 * @param contigs Group of contigs requiring centroid representation for clustering
	 * @return Centroid object with the largest contig as representative, or null if the group is empty
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
	 * Performs deep comparison of contig assignments across clusters to see if assignments have stabilized.
	 * @param a First assignment set from the current iteration
	 * @param b Second assignment set from the previous iteration
	 * @return true if all assignments are identical, false otherwise
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
	 * Creates a deep copy of assignment sets for comparison in subsequent clustering iterations.
	 * Copies the list structure and references so that convergence checks are not affected by later mutations.
	 * @param original Original assignment set from the current iteration
	 * @return Deep copy of the assignment structure with independent ArrayList instances
	 */
	private ArrayList<ArrayList<Contig>> deepCopyAssignments(ArrayList<ArrayList<Contig>> original){
		ArrayList<ArrayList<Contig>> copy=new ArrayList<>();
		for(ArrayList<Contig> list : original){
			copy.add(new ArrayList<>(list)); //Create new ArrayList with same Contig references
		}
		return copy;
	}
	
	/** Represents a cluster centroid using a single representative contig for similarity calculations.
	 * Provides a lightweight abstraction where centroids are contigs rather than averaged feature vectors. */
	private static class Centroid {
		final Contig representative;
		
		/** Creates a centroid with the specified representative contig as the cluster center.
		 * @param rep Representative contig used for similarity calculations */
		Centroid(Contig rep){representative=rep;}
		
		/**
		 * Calculates Oracle similarity between this centroid and a target contig for assignment.
		 * @param contig Target contig for similarity-based distance calculation
		 * @param oracle Oracle instance providing biological sequence similarity calculations
		 * @return Similarity score between centroid representative and target contig (higher = more similar)
		 */
		float similarityTo(Contig contig, Oracle oracle){
			return oracle.similarity(representative, contig, 1.0f);
		}
	}
	
	private final Oracle oracle;
	private final int maxIterations;
	private final float convergenceThreshold;
	private final float minSplitImprovement;
	private final Random random;
	
	private final boolean debug;
	private int splitAttempts;
	private int successfulSplits;
	
	/**
	 * Converts refined clusters to IntHashSet representations for compatibility with integer-based algorithms.
	 * Runs binary splitting refinement, then extracts contig IDs from each resulting cluster into IntHashSet collections.
	 * @param input Input cluster to refine and convert to integer set representation
	 * @return List of IntHashSets containing contig IDs for each refined cluster, or null if refinement fails
	 */
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