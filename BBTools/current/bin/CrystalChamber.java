package bin;

import java.util.ArrayList;
import java.util.Random;

/**
 * Recrystallization-based bin refinement using centroid clustering.
 * Dissolves clusters and rebuilds them using iterative centroid assignment
 * to find natural partitions that greedy algorithms might miss.
 * 
 * "Like chemistry - sometimes you need to dissolve and re-precipitate 
 * to allow components to find their optimal, pure state." - Kairos
 * 
 * @author UMP45
 */
class CrystalChamber extends AbstractRefiner {
    
    public CrystalChamber(Oracle oracle_) {
        this(oracle_, 12345); // Default seed for backward compatibility
    }
    
    /**
     * Creates a CrystalChamber refiner with specified Oracle and configurable random seed.
     * @param oracle_ Oracle instance for contig similarity evaluation during clustering
     * @param seed Random seed for deterministic centroid initialization and reproducible results
     */
    public CrystalChamber(Oracle oracle_, long seed) {
        oracle=oracle_;
        maxIterations=50;
        convergenceThreshold=0.01f;
        minSplitImprovement=0.1f;
        random=new Random(seed); // Deterministic results with cluster-specific seed
    }
    
    @Override
    ArrayList<Bin> refine(Bin input) {
        if(input==null || input.numContigs()<4){return null;}
        if(!input.isCluster()){return null;}
        
        Cluster cluster=(Cluster) input;
        ArrayList<Contig> contigs=new ArrayList<>(cluster.contigs);
        
        if(contigs.size()<4){return null;} // Need minimum size for meaningful split
        
        // Attempt binary recrystallization
        ArrayList<Cluster> crystals=recrystallize(contigs, 2);
        
        if(crystals==null || crystals.size() != 2){return null;}
        
        // Validate split quality
        ArrayList<Bin> result=new ArrayList<Bin>(crystals);
        if(!isSplitBeneficial(input, result)){return null;}
        
        // Additional validation: would Oracle recommend merging these back?
        if(shouldMergeBack(crystals.get(0), crystals.get(1))){return null;}
        
        return result;
    }
    
    /**
     * Performs iterative centroid-based clustering to separate contigs.
     */
    private ArrayList<Cluster> recrystallize(ArrayList<Contig> contigs, int k) {
        if(contigs.size()<k){return null;}
        
        // Initialize centroids using most dissimilar contigs
        ArrayList<Centroid> centroids=initializeCentroids(contigs, k);
        if(centroids==null){return null;}
        
        ArrayList<ArrayList<Contig>> assignments=new ArrayList<>(k);
        for(int i=0; i<k; i++) {
            assignments.add(new ArrayList<Contig>());
        }
        
        // Iterative refinement
        for(int iter=0; iter<maxIterations; iter++) {
            // Clear previous assignments
            for(ArrayList<Contig> list : assignments) {
                list.clear();
            }
            
            // Assign each contig to nearest centroid
            for(Contig contig : contigs) {
                int bestCentroid=findNearestCentroid(contig, centroids);
                assignments.get(bestCentroid).add(contig);
            }
            
            // Check for empty clusters (bad initialization)
            boolean hasEmpty=false;
            for(ArrayList<Contig> list : assignments) {
                if(list.isEmpty()) {
                    hasEmpty=true;
                    break;
                }
            }
            if(hasEmpty){return null;}
            
            // Update centroids and check convergence
            boolean converged=true;
            for(int i=0; i<k; i++) {
                Centroid newCentroid=calculateCentroid(assignments.get(i));
                if(centroidDistance(centroids.get(i), newCentroid)>convergenceThreshold) {
                    converged=false;
                }
                centroids.set(i, newCentroid);
            }
            
            if(converged) {break;}
        }
        
        // Convert assignments to clusters
        ArrayList<Cluster> result=new ArrayList<>(k);
        for(int i=0; i<k; i++) {
            if(assignments.get(i).isEmpty()){return null;}
            
            Cluster cluster=new Cluster(BinObject.globalTime++);
            for(Contig contig : assignments.get(i)) {
                cluster.add(contig);
            }
            result.add(cluster);
        }
        
        return result;
    }
    
    /**
     * Initialize centroids by finding most dissimilar contigs.
     */
    private ArrayList<Centroid> initializeCentroids(ArrayList<Contig> contigs, int k) {
        if(contigs.size()<k){return null;}
        
        ArrayList<Centroid> centroids=new ArrayList<>(k);
        ArrayList<Contig> chosen=new ArrayList<>(k);
        
        // Choose first centroid randomly
        Contig first=contigs.get(random.nextInt(contigs.size()));
        chosen.add(first);
        centroids.add(new Centroid(first));
        
        // Choose remaining centroids to maximize dissimilarity
        for(int i=1; i<k; i++) {
            Contig best=null;
            float maxMinDistance=-1;
            
            for(Contig candidate : contigs) {
                if(chosen.contains(candidate)) {continue;}
                
                // Find minimum distance to existing centroids
                float minDistance=Float.MAX_VALUE;
                for(Contig existing : chosen) {
                    float similarity=oracle.similarity(candidate, existing, 1.0f);
                    float distance=1.0f - similarity; // Convert similarity to distance
                    minDistance=Math.min(minDistance, distance);
                }
                
                // Choose candidate with maximum minimum distance (furthest from all)
                if(minDistance>maxMinDistance) {
                    maxMinDistance=minDistance;
                    best=candidate;
                }
            }
            
            if(best==null){return null;}
            chosen.add(best);
            centroids.add(new Centroid(best));
        }
        
        return centroids;
    }
    
    /**
     * Find the centroid nearest to given contig.
     */
    private int findNearestCentroid(Contig contig, ArrayList<Centroid> centroids) {
        int best=0;
        float bestSimilarity=-1;
        
        for(int i=0; i<centroids.size(); i++) {
            float similarity=centroids.get(i).similarityTo(contig, oracle);
            if(similarity>bestSimilarity) {
                bestSimilarity=similarity;
                best=i;
            }
        }
        
        return best;
    }
    
    /**
     * Calculate centroid of a group of contigs.
     */
    private Centroid calculateCentroid(ArrayList<Contig> contigs) {
        if(contigs.isEmpty()){return null;}
        if(contigs.size()==1){return new Centroid(contigs.get(0));}
        
        // For now, use the largest contig as representative centroid
        // TODO: Could implement proper averaging of features
        Contig largest=contigs.get(0);
        for(Contig c : contigs) {
            if(c.size()>largest.size()) {largest=c;}
        }
        
        return new Centroid(largest);
    }
    
    /**
     * Calculate distance between two centroids.
     */
    private float centroidDistance(Centroid a, Centroid b) {
        if(a==null || b==null){return Float.MAX_VALUE;}
        float similarity=a.similarityTo(b.representative, oracle);
        return 1.0f - similarity;
    }
    
    /**
     * Test whether Oracle would recommend merging the split clusters back together.
     */
    private boolean shouldMergeBack(Cluster a, Cluster b) {
        float similarity=oracle.similarity(a, b, 1.0f);
        return similarity>minSplitImprovement; // If high similarity, don't split
    }
    
    /**
     * Inner class representing a cluster centroid.
     */
    private static class Centroid {
        final Contig representative;
        
        Centroid(Contig rep) {representative=rep;}
        
        float similarityTo(Contig contig, Oracle oracle) {
            return oracle.similarity(representative, contig, 1.0f);
        }
    }
    
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
    
    // Configuration parameters
    private final Oracle oracle;
    private final int maxIterations;
    private final float convergenceThreshold;
    private final float minSplitImprovement;
    private final Random random;
}
