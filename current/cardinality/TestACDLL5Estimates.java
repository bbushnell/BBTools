package cardinality;

/**
 * Head-to-head DLC vs LDLC comparison between CDLL5 and ACDLL5.
 * Feeds identical hashes, compares estimates at checkpoints.
 * Should be identical until overflow (ACDLL5 has 10 states vs 8).
 */
public class TestACDLL5Estimates{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int buckets=2048;
		final CompressedDynamicLogLog5 cdll5=new CompressedDynamicLogLog5(buckets, 31, 42, 0);
		final ArithmeticCompressedDynamicLogLog5 acdll5=new ArithmeticCompressedDynamicLogLog5(buckets, 31, 42, 0);

		// Checkpoints to compare estimates
		final long[] checkpoints={100, 500, 1000, 5000, 10000, 50000, 100000, 500000};

		System.out.printf("%-12s %-14s %-14s %-14s %-14s %-14s %-14s%n",
			"Adds", "C_DLC", "A_DLC", "C_LDLC", "A_LDLC", "C_Hybrid", "A_Hybrid");
		System.out.println("=".repeat(100));

		int cpIdx=0;
		for(long i=0; i<checkpoints[checkpoints.length-1]; i++){
			cdll5.hashAndStore(i);
			acdll5.hashAndStore(i);

			if(cpIdx<checkpoints.length && i+1==checkpoints[cpIdx]){
				// Force re-summarize
				cdll5.lastCardinality=-1;
				acdll5.lastCardinality=-1;

				// Get raw estimates via CardStats
				final double[] cRaw=cdll5.rawEstimates();
				final CardStats cStats=cdll5.consumeLastSummarized();
				final double[] aRaw=acdll5.rawEstimates();
				final CardStats aStats=acdll5.consumeLastSummarized();

				final double cDlc=cStats.dlcRaw();
				final double aDlc=aStats.dlcRaw();
				final double cLdlc=cStats.ldlc();
				final double aLdlc=aStats.ldlc();
				final double cHyb=cStats.hybridDLL();
				final double aHyb=aStats.hybridDLL();

				System.out.printf("%-12d %-14.4f %-14.4f %-14.4f %-14.4f %-14.4f %-14.4f%n",
					checkpoints[cpIdx], cDlc, aDlc, cLdlc, aLdlc, cHyb, aHyb);

				// Flag any divergence
				if(Math.abs(cDlc-aDlc)>0.01){
					System.out.println("  *** DLC DIVERGENCE ***");
				}
				if(Math.abs(cLdlc-aLdlc)>0.01){
					System.out.println("  *** LDLC DIVERGENCE ***");
				}

				cpIdx++;
			}
		}

		// Also compare bucket states at final checkpoint
		System.out.println("\nFinal minZeros: CDLL5="+cdll5.getMinZeros()+" ACDLL5="+acdll5.getMinZeros());
		int bucketMismatches=0;
		for(int b=0; b<buckets; b++){
			final int cReg=cdll5.readBucket(b);
			final int cExp=(cReg>>>2)&0x7;
			final int cHist=cReg&0x3;
			final int aExp=acdll5.readExponent(b);
			final int aHist=acdll5.readHistory(b);
			if(cExp!=aExp || cHist!=aHist){
				bucketMismatches++;
				if(bucketMismatches<=10){
					System.out.printf("  b%d: CDLL5(exp=%d,hist=%d) ACDLL5(exp=%d,hist=%d)%n",
						b, cExp, cHist, aExp, aHist);
				}
			}
		}
		System.out.println("Bucket mismatches: "+bucketMismatches+" / "+buckets);
	}

}
