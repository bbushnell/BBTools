package aligner;

import dna.AminoAcid;
import structures.IntList;
import structures.IntListHashMap;

/**
 * k-mer index for fast sequence alignment candidate identification.
 * Builds a hash map from k-mers to their positions in a reference sequence.
 * Designed for efficient lookup of potential alignment positions.
 * @author Brian Bushnell
 */
public class IntIndex {
    /** Length of k-mers used for indexing */
    final int k;
    /** Hash map storing k-mer values to lists of reference positions */
    final IntListHashMap map;  // kmer -> IntList of positions
    
    /**
     * Constructs an IntIndex by building k-mer position mappings for a reference.
     * Creates a hash map with capacity based on reference length and indexes
     * all valid k-mers found in the sequence.
     *
     * @param ref Reference sequence to index
     * @param k_ Length of k-mers to extract and index
     */
    public IntIndex(byte[] ref, int k_) {
        k = k_;
        map = new IntListHashMap(ref.length * 2);
        indexRef(ref);
    }
    
    /**
     * Indexes all k-mers in the reference sequence into the position map.
     * Uses rolling hash approach to efficiently compute k-mer values.
     * Skips ambiguous bases and resets k-mer computation when encountered.
     * Stores start position (i - k + 1) for each valid k-mer.
     *
     * @param ref Reference sequence bytes to index
     */
    private void indexRef(byte[] ref) {
        if(ref.length < k) return;
        
        int kmer = 0;
        int mask = (1 << (2*k)) - 1;  // For k<=15
        int len = 0;
        
        for(int i = 0; i < ref.length; i++) {
            byte b = ref[i];
            int x = AminoAcid.baseToNumber[b];
            
            if(x < 0) {
                len = 0;
                kmer = 0;
            } else {
                kmer = ((kmer << 2) | x) & mask;
                if(++len >= k) {
                    IntList list = map.get(kmer);
                    if(list == null) {
                        list = new IntList(4);
                        map.put(kmer, list);
                    }
                    list.add(i - k + 1);  // Start position
                }
            }
        }
    }
    
    /**
     * FunctionUnclear: Method appears to be a placeholder for finding alignment
     * candidates from a query sequence but is not implemented. Expected to extract
     * k-mers from query and lookup their positions in the index.
     *
     * @param query Query sequence to find candidates for
     * @param maxHits Maximum number of candidate hits to return
     * @return List of candidate alignment positions (currently empty)
     */
    public IntList getCandidates(byte[] query, int maxHits) {
        IntList candidates = new IntList();
        // Extract k-mers from query and lookup
        // Could use spaced seeds or multiple k values
        return candidates;
    }
}