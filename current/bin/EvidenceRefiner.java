package bin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;
import structures.IntHashSet;
import structures.LongHashMap;

/**
 * Evidence-based bin refinement using DBSCAN-style density clustering.
 * Identifies dense regions of similar contigs separated by sparse boundaries,
 * automatically determining cluster count and identifying outliers.
 * 
 * "Like ecology - find dense populations separated by empty niches,
 * automatically discovering natural community boundaries." - EvidenceRefiner Philosophy
 * 
 * @author UMP45
 */
class EvidenceRefiner extends AbstractRefiner {
    
    /**
     * Creates an EvidenceRefiner with specified Oracle for similarity calculations.
     * @param oracle_ Oracle instance for contig similarity evaluation
     */
    public EvidenceRefiner(Oracle oracle_) {
        this(oracle_, new AbstractRefiner.EvidenceRefinerParams());
    }
    
    /**
     * Creates an EvidenceRefiner with specified Oracle and parameters.
     * @param oracle_ Oracle instance for contig similarity evaluation
     * @param params EvidenceRefiner-specific parameters
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
     * Performs DBSCAN clustering on contigs using Oracle similarity as distance metric.
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
     * Finds all contigs within epsilon similarity of the given contig using cached similarities.
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
     * Expands cluster by recursively adding density-connected points.
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
     * Validates that clusters have good internal cohesion and external separation.
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
     * Calculates average internal similarity within a cluster.
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
     * Calculates average similarity between two different clusters.
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
     * Converts IntHashSet clusters to Cluster objects for compatibility.
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
     * Restores original cluster references when refinement fails or is rejected.
     * Prevents cluster reference corruption similar to CrystalChamber fix.
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
    
    /**
     * Result container for DBSCAN clustering.
     */
    private static class DBSCANResult {
        final ArrayList<IntHashSet> clusters;
        final IntHashSet noise;
        
        DBSCANResult(ArrayList<IntHashSet> clusters, IntHashSet noise) {
            this.clusters = clusters;
            this.noise = noise;
        }
    }
    
    /** Oracle instance for contig similarity calculations */
    private final Oracle oracle;
    /** DBSCAN epsilon parameter - similarity threshold for neighborhood */
    private final float epsilon;
    /** DBSCAN minPoints parameter - minimum points for core region */
    private final int minPoints;
    /** Minimum cluster size to be considered viable */
    private final int minClusterSize;
    /** Random number generator for deterministic processing */
    private final Random random;
    
    /** Enable debugging output */
    private final boolean debug;
    /** Split attempt counter */
    private int splitAttempts;
    /** Successful split counter */
    private int successfulSplits;
}