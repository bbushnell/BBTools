package cardinality;

/**
 * Validation harness: dump CDLL5 per-tier hist-pattern distribution
 * and compare against MantissaCompare2 mode=ctll bits=2 expectations.
 *
 * Simulator steady-state (tier 7) for 2-bit hist:
 *   P(0)=0.098  P(1)=0.333  P(2)=0.056  P(3)=0.514
 * Low-tier saturated (tier 3-5):
 *   P(3) dominates (~0.94) — both below-tier obs have fired.
 * If CDLL5 produces wildly different shape, carry-shift is buggy.
 */
public class TestCDLL5Hist {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int buckets=2048;
		long N=4_000_000L;
		final long seed=42L;
		if(args.length>=1){N=Long.parseLong(args[0]);}
		if(args.length>=2){buckets=Integer.parseInt(args[1]);}

		final CompressedDynamicLogLog5 cdll=new CompressedDynamicLogLog5(buckets, 31, seed, 0);
		final java.util.Random r=new java.util.Random(seed);
		for(long i=0; i<N; i++){cdll.hashAndStore(r.nextLong());}

		final int globalNLZ=cdll.getMinZeros()-1;
		final int[][] counts=new int[20][4];
		int floorCount=0;
		final int[] floorHist=new int[4];
		for(int b=0; b<buckets; b++){
			final int reg=cdll.readBucket(b);
			final int tp=(reg>>>2)&0x7;
			final int h=reg&0x3;
			if(tp==0){floorCount++; floorHist[h]++;}
			else{
				final int absTier=tp+globalNLZ;
				if(absTier<counts.length){counts[absTier][h]++;}
			}
		}

		System.err.println("=== CDLL5 hist distribution ===");
		System.err.println("globalNLZ="+globalNLZ+" N="+N+" buckets="+buckets
			+" occupancy="+String.format("%.4f", cdll.occupancy()));
		System.err.println("tier  count  P(0)    P(1)    P(2)    P(3)");
		for(int t=0; t<counts.length; t++){
			final int tot=counts[t][0]+counts[t][1]+counts[t][2]+counts[t][3];
			if(tot==0){continue;}
			System.err.printf("%4d  %5d  %.4f  %.4f  %.4f  %.4f%n", t, tot,
				counts[t][0]/(double)tot, counts[t][1]/(double)tot,
				counts[t][2]/(double)tot, counts[t][3]/(double)tot);
		}
		if(floorCount>0){
			System.err.printf("floor(nlz=%d) count=%d  P(0)=%.4f P(1)=%.4f P(2)=%.4f P(3)=%.4f%n",
				globalNLZ, floorCount,
				floorHist[0]/(double)floorCount, floorHist[1]/(double)floorCount,
				floorHist[2]/(double)floorCount, floorHist[3]/(double)floorCount);
		}
	}
}
