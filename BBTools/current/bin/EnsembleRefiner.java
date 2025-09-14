package bin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import structures.IntHashSet;

/**
 * Ensemble-based bin refinement combining multiple clustering strategies.
 * Applies CrystalChamber, GraphRefiner, and EvidenceRefiner in parallel,
 * then uses consensus voting to make robust refinement decisions.
 * 
 * "Like a scientific peer review - get multiple expert opinions,
 * then trust only the decisions supported by consensus." - EnsembleRefiner Philosophy
 * 
 * @author UMP45
 */
class EnsembleRefiner extends AbstractRefiner {
    
    /**
     * Creates an EnsembleRefiner with specified Oracle for similarity calculations.
     * @param oracle_ Oracle instance for contig similarity evaluation
     */
    public EnsembleRefiner(Oracle oracle_) {
        this(oracle_, new AbstractRefiner.EnsembleRefinerParams());
    }
    
    /**
     * Creates an EnsembleRefiner with specified Oracle and parameters.
     * @param oracle_ Oracle instance for contig similarity evaluation
     * @param params EnsembleRefiner-specific parameters
     */
    public EnsembleRefiner(Oracle oracle_, AbstractRefiner.EnsembleRefinerParams params) {
        oracle = oracle_;
        
        // Create constituent refiners with different seeds for diversity
        crystalRefiner = new CrystalChamber(oracle, params.seed);
        graphRefiner = new GraphRefiner(oracle, params.graphParams);
        evidenceRefiner = new EvidenceRefiner(oracle, params.evidenceParams);
        
        consensusThreshold = params.consensusThreshold;
        minMethodsAgreeing = params.minMethodsAgreeing;
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
        
        // Apply all refinement methods
        ArrayList<RefinerResult> results = applyAllRefiners(input, contigs);
        if(results.isEmpty()) {
            return null;
        }
        
        // Build consensus from multiple results
        ConsensusResult consensus = buildConsensus(results, contigs);
        if(consensus == null || consensus.clusters.size() < 2) {
            return null;
        }
        
        // Additional validation: ensure consensus is strong
        if(consensus.confidence < consensusThreshold) {
            return null;
        }
        
        successfulSplits++;
        return consensus.clusters;
    }
    
    @Override
    ArrayList<Bin> refine(Bin input) {
        if(input == null || input.numContigs() < 4) return null;
        if(!input.isCluster()) return null;
        
        Cluster cluster = (Cluster) input;
        ArrayList<Contig> contigs = new ArrayList<>(cluster.contigs);
        
        if(contigs.size() < 4) return null;
        
        splitAttempts++;
        
        // Get IntHashSet consensus using dual-mode method
        ArrayList<IntHashSet> consensusClusters = refineToIntSets(input);
        if(consensusClusters == null) {
            restoreOriginalCluster(cluster);
            return null;
        }
        
        // Convert consensus to Cluster objects
        ArrayList<Bin> refinedBins = convertToCluster(consensusClusters, contigs);
        if(!isSplitBeneficial(input, refinedBins)) {
            restoreOriginalCluster(cluster);
            if(debug && splitAttempts % 100 == 0) 
                System.err.println("EnsembleRefiner: Consensus split not beneficial, attempts=" + splitAttempts);
            return null;
        }
        
        cluster.clear(); // Clean up original cluster
        
        if(debug && splitAttempts % 100 == 0) {
            System.err.println("EnsembleRefiner DEBUG: consensus_clusters=" + consensusClusters.size() + 
                " attempts=" + splitAttempts + " successes=" + successfulSplits);
        }
        
        return refinedBins;
    }
    
    /**
     * Applies all available refinement methods to the input cluster.
     */
    private ArrayList<RefinerResult> applyAllRefiners(Bin input, ArrayList<Contig> contigs) {
        ArrayList<RefinerResult> results = new ArrayList<>();
        
        // Try CrystalChamber
        try {
            ArrayList<Bin> crystalResult = crystalRefiner.refine(input);
            if(crystalResult != null && crystalResult.size() >= 2) {
                results.add(new RefinerResult("CRYSTAL", convertToIntHashSets(crystalResult, contigs), 1.0f));
            }
        } catch(Exception e) {
            // Silently continue if one method fails
        }
        
        // Try GraphRefiner
        try {
            ArrayList<Bin> graphResult = graphRefiner.refine(input);
            if(graphResult != null && graphResult.size() >= 2) {
                results.add(new RefinerResult("GRAPH", convertToIntHashSets(graphResult, contigs), 1.0f));
            }
        } catch(Exception e) {
            // Silently continue if one method fails
        }
        
        // Try EvidenceRefiner
        try {
            ArrayList<Bin> evidenceResult = evidenceRefiner.refine(input);
            if(evidenceResult != null && evidenceResult.size() >= 2) {
                results.add(new RefinerResult("EVIDENCE", convertToIntHashSets(evidenceResult, contigs), 1.0f));
            }
        } catch(Exception e) {
            // Silently continue if one method fails
        }
        
        return results;
    }
    
    /**
     * Converts Bin results to IntHashSet format for consensus building.
     */
    private ArrayList<IntHashSet> convertToIntHashSets(ArrayList<Bin> bins, ArrayList<Contig> contigs) {
        ArrayList<IntHashSet> result = new ArrayList<>();
        
        // Create mapping from contig ID to index
        Map<Integer, Integer> idToIndex = new HashMap<>();
        for(int i = 0; i < contigs.size(); i++) {
            idToIndex.put(contigs.get(i).id(), i);
        }
        
        for(Bin bin : bins) {
            if(!bin.isCluster()) continue;
            
            Cluster cluster = (Cluster) bin;
            IntHashSet indexSet = new IntHashSet();
            
            for(Contig contig : cluster.contigs) {
                Integer index = idToIndex.get(contig.id());
                if(index != null) {
                    indexSet.add(index);
                }
            }
            
            if(indexSet.size() > 0) {
                result.add(indexSet);
            }
        }
        
        return result;
    }
    
    /**
     * Builds consensus clustering from multiple refiner results.
     * Uses co-occurrence counting to identify robust cluster boundaries.
     */
    private ConsensusResult buildConsensus(ArrayList<RefinerResult> results, ArrayList<Contig> contigs) {
        if(results.size() < minMethodsAgreeing) return null;
        
        int n = contigs.size();
        
        // Build co-occurrence matrix
        float[][] coOccurrence = new float[n][n];
        for(RefinerResult result : results) {
            for(IntHashSet cluster : result.clusters) {
                // Add co-occurrence for all pairs in this cluster
                int[] indices = cluster.toArray();
                for(int i = 0; i < indices.length; i++) {
                    for(int j = i + 1; j < indices.length; j++) {
                        coOccurrence[indices[i]][indices[j]] += result.confidence;
                        coOccurrence[indices[j]][indices[i]] += result.confidence;
                    }
                }
            }
        }
        
        // Normalize by number of methods
        float maxPossibleScore = results.size();
        for(int i = 0; i < n; i++) {
            for(int j = i + 1; j < n; j++) {
                coOccurrence[i][j] /= maxPossibleScore;
                coOccurrence[j][i] /= maxPossibleScore;
            }
        }
        
        // Build consensus clusters using threshold-based connectivity
        boolean[] assigned = new boolean[n];
        ArrayList<IntHashSet> consensusClusters = new ArrayList<>();
        
        for(int seed = 0; seed < n; seed++) {
            if(assigned[seed]) continue;
            
            IntHashSet cluster = new IntHashSet();
            growConsensusCluster(seed, cluster, assigned, coOccurrence, consensusThreshold);
            
            if(cluster.size() >= 2) {
                consensusClusters.add(cluster);
            }
        }
        
        // Calculate overall consensus confidence
        float totalConfidence = 0.0f;
        int totalPairs = 0;
        
        for(IntHashSet cluster : consensusClusters) {
            int[] indices = cluster.toArray();
            for(int i = 0; i < indices.length; i++) {
                for(int j = i + 1; j < indices.length; j++) {
                    totalConfidence += coOccurrence[indices[i]][indices[j]];
                    totalPairs++;
                }
            }
        }
        
        float averageConfidence = totalPairs > 0 ? totalConfidence / totalPairs : 0.0f;
        
        return new ConsensusResult(consensusClusters, averageConfidence);
    }
    
    /**
     * Grows a consensus cluster using connectivity threshold.
     */
    private void growConsensusCluster(int seed, IntHashSet cluster, boolean[] assigned, 
                                    float[][] coOccurrence, float threshold) {
        
        ArrayList<Integer> toProcess = new ArrayList<>();
        toProcess.add(seed);
        int processed = 0;
        
        while(processed < toProcess.size()) {
            int current = toProcess.get(processed);
            processed++;
            
            if(assigned[current]) continue;
            
            assigned[current] = true;
            cluster.add(current);
            
            // Find strongly connected neighbors
            for(int neighbor = 0; neighbor < coOccurrence.length; neighbor++) {
                if(!assigned[neighbor] && coOccurrence[current][neighbor] >= threshold) {
                    if(!toProcess.contains(neighbor)) {
                        toProcess.add(neighbor);
                    }
                }
            }
        }
    }
    
    /**
     * Converts IntHashSet consensus clusters to Cluster objects.
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
     * Container for individual refiner results.
     */
    private static class RefinerResult {
        final String method;
        final ArrayList<IntHashSet> clusters;
        final float confidence;
        
        RefinerResult(String method, ArrayList<IntHashSet> clusters, float confidence) {
            this.method = method;
            this.clusters = clusters;
            this.confidence = confidence;
        }
    }
    
    /**
     * Container for consensus clustering result.
     */
    private static class ConsensusResult {
        final ArrayList<IntHashSet> clusters;
        final float confidence;
        
        ConsensusResult(ArrayList<IntHashSet> clusters, float confidence) {
            this.clusters = clusters;
            this.confidence = confidence;
        }
    }
    
    /** Oracle instance for contig similarity calculations */
    private final Oracle oracle;
    /** Individual refiner instances */
    private final CrystalChamber crystalRefiner;
    private final GraphRefiner graphRefiner;  
    private final EvidenceRefiner evidenceRefiner;
    /** Consensus threshold for co-occurrence decisions */
    private final float consensusThreshold;
    /** Minimum number of methods that must agree */
    private final int minMethodsAgreeing;
    /** Random number generator for deterministic consensus tie-breaking */
    private final Random random;
    
    /** Enable debugging output */
    private final boolean debug;
    /** Split attempt counter */
    private int splitAttempts;
    /** Successful split counter */
    private int successfulSplits;
}