package aligner;

import dna.AminoAcid;
import structures.IntListHashMap;
import structures.IntList;

public class IntIndex {
    final int k;
    final IntListHashMap map;  // kmer -> IntList of positions
    
    public IntIndex(byte[] ref, int k_) {
        k = k_;
        map = new IntListHashMap(ref.length * 2);
        indexRef(ref);
    }
    
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
    
    public IntList getCandidates(byte[] query, int maxHits) {
        IntList candidates = new IntList();
        // Extract k-mers from query and lookup
        // Could use spaced seeds or multiple k values
        return candidates;
    }
}