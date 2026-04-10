package cardinality;

public class BCLLTraceTest {
    public static void main(String[] args) {
        int numWords = 512; // 3072 buckets
        BankedCeilingLogLog bcll = new BankedCeilingLogLog(numWords * 6, 31, 42, 0);
        
        System.err.println("Created BCLL: " + bcll.getNumWords() + " words, "
            + bcll.getModBuckets() + " buckets");
        
        // Add elements and trace
        long[] checkpoints = {10, 100, 1000, 5000, 10000, 50000, 100000};
        int ci = 0;
        
        for (long i = 1; i <= 100001 && ci < checkpoints.length; i++) {
            bcll.add(i);
            if (i == checkpoints[ci]) {
                double[] raw = bcll.rawEstimates();
                // Count tier distribution
                int[] tierHist = new int[16];
                for (int w = 0; w < bcll.getNumWords(); w++) {
                    int word = bcll.getWord(w);
                    int le = (word >>> 12) & 0xF;
                    int absCeiling = bcll.getGlobalCeiling() + le;
                    if (absCeiling < 16) tierHist[absCeiling]++;
                }
                System.err.printf("trueCard=%d  globalCeiling=%d  Mean=%.1f  LC=%.1f  Hybrid=%.1f  DLC=%.1f  HC=%.1f%n",
                    i, bcll.getGlobalCeiling(), raw[0], raw[5], raw[6], raw[11], raw[9]);
                System.err.print("  tier distribution: ");
                for (int t = 0; t < 16; t++) {
                    if (tierHist[t] > 0) System.err.printf("T%d=%d ", t, tierHist[t]);
                }
                System.err.println();
                ci++;
            }
        }
    }
}
