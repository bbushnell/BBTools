package bin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import map.IntHashSet;
import shared.Random;

class EnsembleRefiner extends AbstractRefiner {
    
    public EnsembleRefiner(Oracle oracle_) {
        this(oracle_, new AbstractRefiner.EnsembleRefinerParams());
    }
    
    public EnsembleRefiner(Oracle oracle_, AbstractRefiner.EnsembleRefinerParams params) {
        oracle = oracle_;
        
        // Create constituent refiners with different seeds for diversity
        crystalRefiner = new CrystalChamber(oracle, params.seed);
        graphRefiner = new GraphRefiner(oracle, params.graphParams);
        evidenceRefiner = new EvidenceRefiner(oracle, params.evidenceParams);
        
        consensusThreshold = params.consensusThreshold;
        minMethodsAgreeing = params.minMethodsAgreeing;
        random = shared.Shared.random(params.seed);
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
    
    private void restoreOriginalCluster(Cluster originalCluster) {
        for(Contig contig : originalCluster.contigs) {
            contig.cluster = originalCluster;
        }
    }
    
    /** Container for individual refiner results with method identification and confidence.
     * Stores clustering output from a single refinement algorithm along with metadata. */
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
    
    private static class ConsensusResult {
        final ArrayList<IntHashSet> clusters;
        final float confidence;
        
        ConsensusResult(ArrayList<IntHashSet> clusters, float confidence) {
            this.clusters = clusters;
            this.confidence = confidence;
        }
    }
    
    private final Oracle oracle;
    private final CrystalChamber crystalRefiner;
    private final GraphRefiner graphRefiner;  
    private final EvidenceRefiner evidenceRefiner;
    private final float consensusThreshold;
    private final int minMethodsAgreeing;
    private final Random random;
    
    private final boolean debug;
    private int splitAttempts;
    private int successfulSplits;
}