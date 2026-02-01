package bin;

import java.util.ArrayList;
import java.util.HashMap;

import map.IntHashSet;
import shared.Random;

/**
 * Graph-based bin refinement using modularity maximization.
 * Constructs a contig similarity graph and applies community detection algorithms to find natural cluster boundaries that may be missed by centroid methods.
 * @author UMP45
 */
class GraphRefiner extends AbstractRefiner {
    
    /** Creates a GraphRefiner with the specified Oracle for similarity calculations using default GraphRefinerParams.
     * @param oracle_ Oracle instance for contig similarity evaluation */
    public GraphRefiner(Oracle oracle_) {
        this(oracle_, new AbstractRefiner.GraphRefinerParams());
    }
    
    /**
     * Creates a GraphRefiner with the specified Oracle and parameters and initializes the random number generator for deterministic behavior.
     * @param oracle_ Oracle instance for contig similarity evaluation
     * @param params GraphRefiner-specific parameters including edge weight threshold
     */
    public GraphRefiner(Oracle oracle_, AbstractRefiner.GraphRefinerParams params) {
        oracle = oracle_;
        minEdgeWeight = params.minEdgeWeight;
        maxIterations = params.maxIterations;
        random = shared.Shared.random(params.seed);
        debug = true;
        splitAttempts = 0;
        successfulSplits = 0;
    }
    
    /**
     * Refines a bin into communities represented as IntHashSet collections.
     * Builds a similarity graph, applies modularity-based community detection, and validates that new communities improve modularity over a single cluster.
     * @param input Bin to refine (must be Cluster with 4+ contigs)
     * @return List of communities as IntHashSet, or null if refinement unsuccessful
     */
    @Override
    ArrayList<IntHashSet> refineToIntSets(Bin input) {
        if(input == null || input.numContigs() < 4) return null;
        if(!input.isCluster()) return null;
        
        Cluster cluster = (Cluster) input;
        ArrayList<Contig> contigs = new ArrayList<>(cluster.contigs);
        
        if(contigs.size() < 4) return null;
        
        splitAttempts++;
        
        // Build similarity graph
        SimilarityGraph graph = buildSimilarityGraph(contigs);
        if(graph == null || graph.getEdgeCount() < 2) {
            return null;
        }
        
        // Apply modularity-based community detection
        ArrayList<IntHashSet> communities = detectCommunities(graph, contigs);
        if(communities == null || communities.size() < 2) {
            return null;
        }
        
        // Additional validation: check modularity improvement
        float originalModularity = calculateModularity(graph, createSingleCommunity(contigs));
        float newModularity = calculateModularity(graph, communities);
        if(newModularity <= originalModularity + 0.1f) {
            return null;
        }
        
        successfulSplits++;
        return communities;
    }
    
    /**
     * Refines a bin into multiple Cluster objects using graph-based community detection.
     * Converts IntHashSet communities to Cluster objects, checks that the split is beneficial, and restores the original cluster if refinement is not beneficial.
     * @param input Bin to refine (must be Cluster with 4+ contigs)
     * @return List of refined Cluster objects, or null if refinement unsuccessful
     */
    @Override
    ArrayList<Bin> refine(Bin input) {
        if(input == null || input.numContigs() < 4) return null;
        if(!input.isCluster()) return null;
        
        Cluster cluster = (Cluster) input;
        ArrayList<Contig> contigs = new ArrayList<>(cluster.contigs);
        
        // Get IntHashSet communities using dual-mode method
        ArrayList<IntHashSet> communities = refineToIntSets(input);
        if(communities == null) {
            restoreOriginalCluster(cluster);
            return null;
        }
        
        // Convert to Cluster objects
        ArrayList<Bin> result = convertToCluster(communities, contigs);
        if(!isSplitBeneficial(input, result)) {
            restoreOriginalCluster(cluster);
            if(debug && splitAttempts % 100 == 0) 
                System.err.println("GraphRefiner: Split not beneficial, attempts=" + splitAttempts);
            return null;
        }
        
        cluster.clear(); // Clean up original cluster
        
        if(debug && splitAttempts % 100 == 0) {
            System.err.println("GraphRefiner DEBUG: communities=" + communities.size() + 
                " attempts=" + splitAttempts + " successes=" + successfulSplits);
        }
        
        return result;
    }
    
    /**
     * Builds a weighted similarity graph from contigs.
     * Creates edges between contig pairs whose distance is below the configured threshold, converting similarity to distance for threshold comparison but keeping similarity as edge weight.
     * @param contigs List of contigs to build graph from
     * @return SimilarityGraph with edges above minimum weight threshold
     */
    private SimilarityGraph buildSimilarityGraph(ArrayList<Contig> contigs) {
        SimilarityGraph graph = new SimilarityGraph(contigs.size());
        
        // Create edges for all pairs above threshold
        for(int i = 0; i < contigs.size(); i++) {
            for(int j = i + 1; j < contigs.size(); j++) {
                float similarity = oracle.similarity(contigs.get(i), contigs.get(j), 1.0f);
                if(similarity > 0) { // Only process valid similarities
                    float distance = 1.0f - similarity; // Convert similarity to distance like CrystalChamber
                    float weight = similarity; // Keep similarity as edge weight for modularity calculation
                    if(distance <= minEdgeWeight) { // Use distance for threshold comparison
                        graph.addEdge(i, j, weight);
                    }
                }
            }
        }
        
        return graph;
    }
    
    /**
     * Applies Louvain-style modularity maximization for community detection.
     * Iteratively moves nodes to communities that maximize modularity gain, processing nodes in random order to avoid bias in community assignment.
     * @param graph Similarity graph of contigs
     * @param contigs Original contig list for size validation
     * @return List of detected communities as IntHashSet, or null if no improvement
     */
    private ArrayList<IntHashSet> detectCommunities(SimilarityGraph graph, ArrayList<Contig> contigs) {
        int n = contigs.size();
        int[] communities = new int[n];
        
        // Initialize: each node in its own community
        for(int i = 0; i < n; i++) {
            communities[i] = i;
        }
        
        boolean improved = true;
        int iterations = 0;
        
        while(improved && iterations < maxIterations) {
            improved = false;
            iterations++;
            
            // Process nodes in random order to avoid bias
            int[] nodeOrder = generateRandomOrder(n);
            
            for(int nodeIndex : nodeOrder) {
                int bestCommunity = communities[nodeIndex];
                float bestGain = 0.0f;
                
                // Try moving to neighbor communities
                for(int neighbor : graph.getNeighbors(nodeIndex)) {
                    int neighborCommunity = communities[neighbor];
                    if(neighborCommunity != communities[nodeIndex]) {
                        float gain = calculateModularityGain(graph, nodeIndex, neighborCommunity, communities);
                        if(gain > bestGain) {
                            bestGain = gain;
                            bestCommunity = neighborCommunity;
                        }
                    }
                }
                
                // Move node if beneficial
                if(bestCommunity != communities[nodeIndex]) {
                    communities[nodeIndex] = bestCommunity;
                    improved = true;
                }
            }
        }
        
        // Convert community array to IntHashSet list
        return groupIntoCommunities(communities);
    }
    
    /**
     * Calculates modularity gain from moving a node to a different community.
     * Compares internal edge weights in the current vs proposed community assignment.
     * @param graph Similarity graph containing edge weights
     * @param node Node index to potentially move
     * @param newCommunity Target community ID
     * @param communities Current community assignments for all nodes
     * @return Modularity gain (positive indicates beneficial move)
     */
    private float calculateModularityGain(SimilarityGraph graph, int node, int newCommunity, int[] communities) {
        float currentContribution = 0.0f;
        float newContribution = 0.0f;
        
        // Calculate change in internal edges
        for(int neighbor : graph.getNeighbors(node)) {
            float weight = graph.getWeight(node, neighbor);
            if(communities[neighbor] == communities[node]) {
                currentContribution += weight;
            }
            if(communities[neighbor] == newCommunity) {
                newContribution += weight;
            }
        }
        
        return newContribution - currentContribution;
    }
    
    /**
     * Groups nodes by community ID into IntHashSet collections and filters out singleton communities to ensure meaningful clustering.
     * @param communities Array mapping node indices to community IDs
     * @return List of communities with 2+ members, or null if insufficient communities
     */
    private ArrayList<IntHashSet> groupIntoCommunities(int[] communities) {
        HashMap<Integer, IntHashSet> communityMap = new HashMap<>();
        
        for(int i = 0; i < communities.length; i++) {
            int communityId = communities[i];
            if(!communityMap.containsKey(communityId)) {
                communityMap.put(communityId, new IntHashSet());
            }
            communityMap.get(communityId).add(i);
        }
        
        // Filter out singleton communities
        ArrayList<IntHashSet> result = new ArrayList<>();
        for(IntHashSet community : communityMap.values()) {
            if(community.size() > 1) {
                result.add(community);
            }
        }
        
        return result.size() >= 2 ? result : null;
    }
    
    /**
     * Calculates modularity score for a given community structure.
     * Uses a Newman-Girvan-style modularity formula comparing internal edges to expected internal edges; higher scores indicate better community structure.
     * @param graph Similarity graph with edge weights
     * @param communities Current community partition
     * @return Modularity score (higher is better, typically in the range -0.5 to 1.0)
     */
    private float calculateModularity(SimilarityGraph graph, ArrayList<IntHashSet> communities) {
        float modularity = 0.0f;
        float totalEdgeWeight = graph.getTotalWeight();
        
        if(totalEdgeWeight == 0) return 0.0f;
        
        for(IntHashSet community : communities) {
            float internalWeight = 0.0f;
            float totalDegree = 0.0f;
            
            // Calculate internal edges and total degree for this community
            int[] nodes = community.toArray();
            for(int node1 : nodes) {
                totalDegree += graph.getDegree(node1);
                for(int node2 : nodes) {
                    if(node1 < node2) {
                        internalWeight += graph.getWeight(node1, node2);
                    }
                }
            }
            
            // Modularity contribution: (internal edges) - (expected internal edges)
            float expectedInternal = (totalDegree * totalDegree) / (4.0f * totalEdgeWeight);
            modularity += (internalWeight / totalEdgeWeight) - (expectedInternal / totalEdgeWeight);
        }
        
        return modularity;
    }
    
    /**
     * Creates a single community containing all nodes for baseline modularity calculation.
     * Used to compare against multi-community partitions to validate improvement.
     * @param contigs List of contigs to include in the single community
     * @return Single-element list containing one community with all contig indices
     */
    private ArrayList<IntHashSet> createSingleCommunity(ArrayList<Contig> contigs) {
        IntHashSet singleCommunity = new IntHashSet();
        for(int i = 0; i < contigs.size(); i++) {
            singleCommunity.add(i);
        }
        ArrayList<IntHashSet> result = new ArrayList<>();
        result.add(singleCommunity);
        return result;
    }
    
    /**
     * Converts IntHashSet communities to Cluster objects for compatibility.
     * Uses the lowest-index contig from each community as cluster ID and populates Cluster objects with the corresponding contigs.
     * @param communities List of communities as IntHashSet of indices
     * @param contigs Original contig list for creating clusters
     * @return List of Cluster objects containing grouped contigs
     */
    private ArrayList<Bin> convertToCluster(ArrayList<IntHashSet> communities, ArrayList<Contig> contigs) {
        ArrayList<Bin> result = new ArrayList<>();
        
        for(IntHashSet community : communities) {
            if(community.size() == 0) continue;
            
            // Find first contig in community for cluster ID
            int firstIndex = -1;
            int[] indices = community.toArray();
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
     * Prevents cluster reference corruption by restoring contig cluster assignments to the original cluster.
     * @param originalCluster The original cluster to restore contig references to
     */
    private void restoreOriginalCluster(Cluster originalCluster) {
        for(Contig contig : originalCluster.contigs) {
            contig.cluster = originalCluster;
        }
    }
    
    /**
     * Generates a random permutation of indices 0..n-1 for unbiased processing.
     * Uses Fisher-Yates shuffle algorithm for uniform randomness.
     * @param n Number of elements to permute
     * @return Array containing randomized order of indices 0 through n-1
     */
    private int[] generateRandomOrder(int n) {
        int[] order = new int[n];
        for(int i = 0; i < n; i++) {
            order[i] = i;
        }
        
        // Fisher-Yates shuffle
        for(int i = n - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int temp = order[i];
            order[i] = order[j];
            order[j] = temp;
        }
        
        return order;
    }
    
    /** Internal representation of weighted similarity graph.
     * Stores contig similarity relationships as a weighted adjacency list optimized for modularity calculations and community detection algorithms. */
    private static class SimilarityGraph {
        private final ArrayList<ArrayList<Edge>> adjacencyList;
        private final int nodeCount;
        private float totalWeight;
        
        public SimilarityGraph(int nodeCount) {
            this.nodeCount = nodeCount;
            this.adjacencyList = new ArrayList<>(nodeCount);
            this.totalWeight = 0.0f;
            
            for(int i = 0; i < nodeCount; i++) {
                adjacencyList.add(new ArrayList<Edge>());
            }
        }
        
        public void addEdge(int u, int v, float weight) {
            adjacencyList.get(u).add(new Edge(v, weight));
            adjacencyList.get(v).add(new Edge(u, weight));
            totalWeight += weight;
        }
        
        public ArrayList<Integer> getNeighbors(int node) {
            ArrayList<Integer> neighbors = new ArrayList<>();
            for(Edge edge : adjacencyList.get(node)) {
                neighbors.add(edge.target);
            }
            return neighbors;
        }
        
        public float getWeight(int u, int v) {
            for(Edge edge : adjacencyList.get(u)) {
                if(edge.target == v) return edge.weight;
            }
            return 0.0f;
        }
        
        public float getDegree(int node) {
            float degree = 0.0f;
            for(Edge edge : adjacencyList.get(node)) {
                degree += edge.weight;
            }
            return degree;
        }
        
        public int getEdgeCount() {
            int count = 0;
            for(ArrayList<Edge> edges : adjacencyList) {
                count += edges.size();
            }
            return count / 2; // Each edge counted twice
        }
        
        public float getTotalWeight() {
            return totalWeight;
        }
        
        private static class Edge {
            final int target;
            final float weight;
            
            Edge(int target, float weight) {
                this.target = target;
                this.weight = weight;
            }
        }
    }
    
    private final Oracle oracle;
    private final float minEdgeWeight;
    private final int maxIterations;
    private final Random random;
    
    private final boolean debug;
    private int splitAttempts;
    private int successfulSplits;
}