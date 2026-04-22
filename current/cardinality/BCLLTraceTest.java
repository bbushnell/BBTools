package cardinality;

/**
 * Integration test that traces BankedCeilingLogLog estimates at several
 * cardinality checkpoints (10 through 100000).
 * Prints globalCeiling, raw estimate variants, and tier distributions
 * at each checkpoint.
 *
 * @author Brian Bushnell
 */
public class BCLLTraceTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int numWords=512; // 3072 buckets
		final BankedCeilingLogLog bcll=new BankedCeilingLogLog(numWords*6, 31, 42, 0);

		System.err.println("Created BCLL: "+bcll.getNumWords()+" words, "
			+bcll.getModBuckets()+" buckets");

		// Add elements and trace
		final long[] checkpoints={10, 100, 1000, 5000, 10000, 50000, 100000};
		int ci=0;

		for(long i=1; i<=100001 && ci<checkpoints.length; i++){
			bcll.add(i);
			if(i==checkpoints[ci]){
				final double[] raw=bcll.rawEstimates();
				// Count tier distribution
				final int[] tierHist=new int[16];
				for(int w=0; w<bcll.getNumWords(); w++){
					final int word=bcll.getWord(w);
					final int le=(word>>>12)&0xF;
					final int absCeiling=bcll.getGlobalCeiling()+le;
					if(absCeiling<16){tierHist[absCeiling]++;}
				}
				System.err.printf("trueCard=%d  globalCeiling=%d  Mean=%.1f  LC=%.1f  Hybrid=%.1f  DLC=%.1f  HC=%.1f%n",
					i, bcll.getGlobalCeiling(), raw[0], raw[5], raw[6], raw[11], raw[9]);
				System.err.print("  tier distribution: ");
				for(int t=0; t<16; t++){
					if(tierHist[t]>0){System.err.printf("T%d=%d ", t, tierHist[t]);}
				}
				System.err.println();
				ci++;
			}
		}
	}
}
