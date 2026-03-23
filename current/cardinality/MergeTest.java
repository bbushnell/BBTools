package cardinality;
public class MergeTest {
    public static void main(String[] args) {
        // Test merge at various thread-like splits
        long N=10000000;
        int splits=16;
        UltraDynamicLogLog6 ref=new UltraDynamicLogLog6(2048, 31, -1, 0);
        UltraDynamicLogLog6[] parts=new UltraDynamicLogLog6[splits];
        for(int s=0; s<splits; s++){parts[s]=new UltraDynamicLogLog6(2048, 31, -1, 0);}
        for(long i=0; i<N; i++){
            ref.add(i);
            parts[(int)(i%splits)].add(i);
        }
        // Merge all into parts[0]
        for(int s=1; s<splits; s++){parts[0].add(parts[s]);}
        CardinalityTracker.clampToAdded=false;
        System.out.printf("UDLL6 16-way merge: ref=%d merged=%d ratio=%.4f%n",
            ref.cardinality(), parts[0].cardinality(), (double)parts[0].cardinality()/ref.cardinality());

        // Same for ErtlULL
        ErtlULL eref=new ErtlULL(2048, 31, -1, 0);
        ErtlULL[] eparts=new ErtlULL[splits];
        for(int s=0; s<splits; s++){eparts[s]=new ErtlULL(2048, 31, -1, 0);}
        for(long i=0; i<N; i++){
            eref.add(i);
            eparts[(int)(i%splits)].add(i);
        }
        for(int s=1; s<splits; s++){eparts[0].add(eparts[s]);}
        System.out.printf("ErtlULL 16-way merge: ref=%d merged=%d ratio=%.4f%n",
            eref.cardinality(), eparts[0].cardinality(), (double)eparts[0].cardinality()/eref.cardinality());

        // Also test with master+worker pattern (like LogLogWrapper)
        UltraDynamicLogLog6 master=new UltraDynamicLogLog6(2048, 31, -1, 0);
        UltraDynamicLogLog6 worker=new UltraDynamicLogLog6(2048, 31, -1, 0);
        for(long i=0; i<N; i++){worker.add(i);}
        master.add(worker);
        System.out.printf("UDLL6 master+worker: ref=%d merged=%d ratio=%.4f (master mz=%d worker mz=%d)%n",
            ref.cardinality(), master.cardinality(), (double)master.cardinality()/ref.cardinality(),
            master.getMinZeros(), worker.getMinZeros());
    }
}
