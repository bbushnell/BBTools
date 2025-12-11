package bin;

import java.util.ArrayList;
import java.util.Random;

import structures.IntHashSet;
import structures.LongHashMap;

class EvidenceRefiner extends AbstractRefiner {
    
    /**
     * Creates an EvidenceRefiner with the specified Oracle and default parameters.
     * Uses the default epsilon, minPoints, and minClusterSize values defined in EvidenceRefinerParams.
     * @param oracle_ Oracle instance for contig similarity evaluation
     */
    public EvidenceRefiner(Oracle oracle_) {
        this(oracle_, new AbstractRefiner.EvidenceRefinerParams());
    }
    
    /**
     * Creates an EvidenceRefiner with specified Oracle and custom parameters.
     * Initializes DBSCAN epsilon, minPoints, minClusterSize, and random seed from the provided params.
     * @param oracle_ Oracle instance for contig similarity evaluation
     * @param params EvidenceRefiner-specific parameters including epsilon, minPoints, and minClusterSize
     */
    public EvidenceRefiner(Oracle oracle_, AbstractRefiner.EvidenceRefinerParams params) {
        oracle = oracle_;
        epsilon = params.epsilon;
        minPoints = params.minPoints;
        minClusterSize = params.minClusterSize;
        random = new Random(params.seed);
        debug = true;
        splitAttempts = 0;
        successfulSplits = 0;
    }
    
    /**
     * Refines a bin using DBSCAN clustering and returns index-based clusters.
     * Applies density-based clustering to identify natural contig groupings from Oracle similarity metrics.
     * @param input The bin to refine (must be a Cluster with at least 4 contigs)
     * @return List of IntHashSets representing refined clusters, or null if refinement fails
     */
    @Override
    ArrayList<IntHashSet> refineToIntSets(Bin input) {
        if(input == null || input.numContigs() < 4) return null;
        if(!input.isCluster()) return null;
        
        Cluster cluster = (Cluster) input;
        ArrayList<Contig> contigs = new ArrayList<>(cluster.contigs);
        
        if(contigs.size() < 4) return null;
        
        splitAttempts++;
        
        // Apply DBSCAN clustering
        DBSCANResult result = performDBSCAN(contigs);
        if(result == null) {
            return null;
        }
        
        // If DBSCAN found no clusters (all noise), split into individual contigs
        if(result.clusters.size() == 0) {
            ArrayList<IntHashSet> individualClusters = new ArrayList<>();
            for(int i = 0; i < contigs.size(); i++) {
                IntHashSet individual = new IntHashSet();
                individual.add(i);
                individualClusters.add(individual);
            }
            return individualClusters;
        }
        
        // If DBSCAN found only 1 cluster, no improvement over original
        if(result.clusters.size() < 2) {
            return null;
        }
        
        // Validate cluster quality - check internal cohesion vs external separation
        if(!hasGoodSeparation(result.clusters, contigs)) {
            return null;
        }
        
        successfulSplits++;
        return result.clusters;
    }
    
    /**
     * Main refinement method that converts an input bin into refined Cluster objects.
     * Calls refineToIntSets for clustering, converts results to Bin format, and validates that the split is beneficial.
     * @param input The bin to refine
     * @return List of refined Bin objects, or null if refinement is not beneficial
     */
    @Override
    ArrayList<Bin> refine(Bin input) {
        if(input == null || input.numContigs() < 4) return null;
        if(!input.isCluster()) return null;
        
        Cluster cluster = (Cluster) input;
        ArrayList<Contig> contigs = new ArrayList<>(cluster.contigs);
        
        if(contigs.size() < 4) return null;
        
        splitAttempts++;
        
        // Get IntHashSet clusters using dual-mode method
        ArrayList<IntHashSet> clusters = refineToIntSets(input);
        if(clusters == null) {
            restoreOriginalCluster(cluster);
            return null;
        }
        
        // Convert to Cluster objects
        ArrayList<Bin> refinedBins = convertToCluster(clusters, contigs);
        if(!isSplitBeneficial(input, refinedBins)) {
            restoreOriginalCluster(cluster);
            if(debug && splitAttempts % 100 == 0) 
                System.err.println("EvidenceRefiner: Split not beneficial, attempts=" + splitAttempts);
            return null;
        }
        
        cluster.clear(); // Clean up original cluster
        
        if(debug && splitAttempts % 100 == 0) {
            System.err.println("EvidenceRefiner DEBUG: clusters=" + clusters.size() + 
                " attempts=" + splitAttempts + " successes=" + successfulSplits);
        }
        
        return refinedBins;
    }
    
    /**
     * Executes DBSCAN clustering algorithm on contigs using cached similarity calculations.
     * Identifies dense regions of similar contigs separated by sparse boundaries and classifies noise points.
     * @param contigs List of contigs to cluster
     * @return DBSCANResult containing identified clusters and noise points
     */
    private DBSCANResult performDBSCAN(ArrayList<Contig> contigs) {
        int n = contigs.size();
        ContigStatus[] status = new ContigStatus[n];
        for(int i = 0; i < n; i++) {
            status[i] = ContigStatus.UNVISITED;
        }
        
        // Initialize similarity cache for this cluster
        LongHashMap similarityCache = new LongHashMap();
        
        ArrayList<IntHashSet> clusters = new ArrayList<>();
        IntHashSet noise = new IntHashSet();
        
        // Process each unvisited point
        for(int i = 0; i < n; i++) {
            if(status[i] != ContigStatus.UNVISITED) continue;
            
            // Find neighbors within epsilon distance
            IntHashSet neighbors = findNeighbors(i, contigs, similarityCache);
            
            if(neighbors.size() < minPoints) {
                // Not enough neighbors - mark as noise for now
                status[i] = ContigStatus.NOISE;
            } else {
                // Start new cluster
                IntHashSet cluster = new IntHashSet();
                expandCluster(i, neighbors, cluster, status, contigs, similarityCache);
                if(cluster.size() >= minClusterSize) {
                    clusters.add(cluster);
                }
            }
        }
        
        // Collect remaining noise points
        for(int i = 0; i < n; i++) {
            if(status[i] == ContigStatus.NOISE) {
                noise.add(i);
            }
        }
        
        return new DBSCANResult(clusters, noise);
    }
    
    /**
     * Finds all contigs within epsilon distance of the specified contig.
     * Uses cached Oracle similarity calculations and converts similarity scores to distances using (1 - similarity).
     * @param index Index of the query contig
     * @param contigs List of all contigs in the cluster
     * @param similarityCache Cache storing pairwise similarity values
     * @return Set of neighbor indices within epsilon distance
     */
    private IntHashSet findNeighbors(int index, ArrayList<Contig> contigs, LongHashMap similarityCache) {
        IntHashSet neighbors = new IntHashSet();
        Contig query = contigs.get(index);
        
        for(int i = 0; i < contigs.size(); i++) {
            if(i == index) continue;
            
            // Create cache key: (minId << 32) | maxId to ensure consistent ordering
            int minId = Math.min(index, i);
            int maxId = Math.max(index, i);
            long cacheKey = (((long)minId) << 32) | ((long)maxId);
            
            float similarity;
            if(similarityCache.containsKey(cacheKey)) {
                // Retrieve cached similarity (stored as similarity * 1000000)
                int cachedValue = similarityCache.get(cacheKey);
                similarity = cachedValue / 1000000.0f;
            } else {
                // Calculate and cache similarity
                similarity = oracle.similarity(query, contigs.get(i), 1.0f);
                int cachedValue = (int)(similarity * 1000000);
                similarityCache.put(cacheKey, cachedValue);
            }
            
            if(similarity > 0) { // Only consider compatible pairs
                float distance = 1.0f - similarity; // Convert to distance like CrystalChamber
                if(distance <= epsilon) {
                    neighbors.add(i);
                }
            }
        }
        
        return neighbors;
    }
    
    /**
     * Recursively expands a cluster by adding density-connected points.
     * Implements the cluster expansion phase of DBSCAN by visiting neighbors and promoting new core points.
     * @param corePoint Index of the initial core point
     * @param neighbors Initial neighbors of the core point
     * @param cluster The cluster being expanded
     * @param status Array tracking processing status of each contig
     * @param contigs List of all contigs
     * @param similarityCache Cache for similarity calculations
     */
    private void expandCluster(int corePoint, IntHashSet neighbors, IntHashSet cluster, 
                             ContigStatus[] status, ArrayList<Contig> contigs, LongHashMap similarityCache) {
        
        cluster.add(corePoint);
        status[corePoint] = ContigStatus.CLUSTERED;
        
        // Process all neighbors
        ArrayList<Integer> toProcess = new ArrayList<>();
        int[] neighborArray = neighbors.toArray();
        for(int neighbor : neighborArray) {
            toProcess.add(neighbor);
        }
        
        int processed = 0;
        while(processed < toProcess.size()) {
            int current = toProcess.get(processed);
            processed++;
            
            if(status[current] == ContigStatus.NOISE) {
                // Border point - add to cluster
                status[current] = ContigStatus.CLUSTERED;
                cluster.add(current);
            } else if(status[current] == ContigStatus.UNVISITED) {
                status[current] = ContigStatus.CLUSTERED;
                cluster.add(current);
                
                // Check if this point is also a core point
                IntHashSet newNeighbors = findNeighbors(current, contigs, similarityCache);
                if(newNeighbors.size() >= minPoints) {
                    // Add new neighbors to processing queue
                    int[] newNeighborArray = newNeighbors.toArray();
                    for(int newNeighbor : newNeighborArray) {
                        if(!toProcess.contains(newNeighbor)) {
                            toProcess.add(newNeighbor);
                        }
                    }
                }
            }
        }
    }
    
    /**
     * Validates cluster quality by checking internal cohesion versus external separation.
     * Requires internal similarity to exceed average external similarity by at least 0.1.
     * @param clusters List of identified clusters
     * @param contigs List of all contigs
     * @return true if clusters show good separation, false otherwise
     */
    private boolean hasGoodSeparation(ArrayList<IntHashSet> clusters, ArrayList<Contig> contigs) {
        if(clusters.size() < 2) return false;
        
        // Calculate average internal similarity for each cluster
        float[] internalSimilarity = new float[clusters.size()];
        for(int c = 0; c < clusters.size(); c++) {
            internalSimilarity[c] = calculateInternalSimilarity(clusters.get(c), contigs);
        }
        
        // Calculate average external similarity between clusters
        float totalExternalSimilarity = 0.0f;
        int comparisons = 0;
        
        for(int i = 0; i < clusters.size(); i++) {
            for(int j = i + 1; j < clusters.size(); j++) {
                float external = calculateExternalSimilarity(clusters.get(i), clusters.get(j), contigs);
                totalExternalSimilarity += external;
                comparisons++;
            }
        }
        
        float avgExternalSimilarity = totalExternalSimilarity / comparisons;
        
        // Require internal similarity to be meaningfully higher than external
        for(float internal : internalSimilarity) {
            if(internal <= avgExternalSimilarity + 0.1f) {
                return false; // Poor separation
            }
        }
        
        return true;
    }
    
    /**
     * Calculates average pairwise similarity within a single cluster.
     * Uses Oracle similarity measurements for all contig pairs within the cluster.
     * @param cluster Set of contig indices in the cluster
     * @param contigs List of all contigs
     * @return Average internal similarity, or 1.0 for single-contig clusters
     */
    private float calculateInternalSimilarity(IntHashSet cluster, ArrayList<Contig> contigs) {
        if(cluster.size() < 2) return 1.0f;
        
        float totalSimilarity = 0.0f;
        int comparisons = 0;
        
        int[] indices = cluster.toArray();
        for(int i = 0; i < indices.length; i++) {
            for(int j = i + 1; j < indices.length; j++) {
                float similarity = oracle.similarity(contigs.get(indices[i]), contigs.get(indices[j]), 1.0f);
                totalSimilarity += similarity;
                comparisons++;
            }
        }
        
        return comparisons > 0 ? totalSimilarity / comparisons : 0.0f;
    }
    
    /**
     * Calculates average similarity between contigs in two different clusters.
     * Measures inter-cluster separation by computing all pairwise similarities between contigs from the two clusters.
     * @param cluster1 First cluster
     * @param cluster2 Second cluster
     * @param contigs List of all contigs
     * @return Average similarity between the two clusters
     */
    private float calculateExternalSimilarity(IntHashSet cluster1, IntHashSet cluster2, ArrayList<Contig> contigs) {
        float totalSimilarity = 0.0f;
        int comparisons = 0;
        
        int[] indices1 = cluster1.toArray();
        int[] indices2 = cluster2.toArray();
        for(int index1 : indices1) {
            for(int index2 : indices2) {
                float similarity = oracle.similarity(contigs.get(index1), contigs.get(index2), 1.0f);
                totalSimilarity += similarity;
                comparisons++;
            }
        }
        
        return comparisons > 0 ? totalSimilarity / comparisons : 0.0f;
    }
    
    /**
     * Converts index-based cluster sets into Cluster objects for system compatibility.
     * Creates new Cluster instances using the lowest-indexed contig ID as cluster identifier and adds member contigs.
     * @param clusterSets List of IntHashSets containing contig indices
     * @param contigs List of all contigs
     * @return List of Cluster objects corresponding to the input sets
     */
    private ArrayList<Bin> convertToCluster(ArrayList<IntHashSet> clusterSets, ArrayList<Contig> contigs) {
        ArrayList<Bin> result = new ArrayList<>();
        
        for(IntHashSet clusterSet : clusterSets) {
            if(clusterSet.size() == 0) continue;
            
            // Find first contig in cluster for ID
            int firstIndex = -1;
            int[] indices = clusterSet.toArray();
            for(int index : indices) {
                if(firstIndex == -1 || index < firstIndex) {
                    firstIndex = index;
                }
            }
            
            Cluster cluster = new Cluster(contigs.get(firstIndex).id());
            for(int index : indices) {
                cluster.add(contigs.get(index));
            }
            result.add(cluster);
        }
        
        return result;
    }
    
    /**
     * Restores contig-to-cluster references when refinement is rejected.
     * Prevents reference corruption by ensuring all contigs point back to the original cluster (similar to CrystalChamber fix).
     * @param originalCluster The original cluster to restore references to
     */
    private void restoreOriginalCluster(Cluster originalCluster) {
        for(Contig contig : originalCluster.contigs) {
            contig.cluster = originalCluster;
        }
    }
    
    /**
     * Status tracking for DBSCAN algorithm.
     */
    private enum ContigStatus {
        UNVISITED,
        CLUSTERED, 
        NOISE
    }
    
    /** Container for DBSCAN clustering results.
     * Separates identified clusters from noise points for downstream processing. */
    private static class DBSCANResult {
        final ArrayList<IntHashSet> clusters;
        final IntHashSet noise;
        
        DBSCANResult(ArrayList<IntHashSet> clusters, IntHashSet noise) {
            this.clusters = clusters;
            this.noise = noise;
        }
    }
    
    private final Oracle oracle;
    private final float epsilon;
    /**
     * DBSCAN minPoints parameter; minimum neighbors required for core point status
     */
    private final int minPoints;
    /** Minimum viable cluster size after DBSCAN processing */
    private final int minClusterSize;
    /** Random number generator for deterministic processing with fixed seed */
    private final Random random;
    
    /** Flag enabling debugging output for cluster analysis */
    private final boolean debug;
    /** Counter tracking total number of refinement attempts */
    private int splitAttempts;
    /** Counter tracking number of successful cluster splits */
    private int successfulSplits;
}