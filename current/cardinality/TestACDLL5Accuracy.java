package cardinality;

/**
 * Head-to-head CDLL5 vs ACDLL5 accuracy comparison.
 * Tracks overflows, state changes, DLC, LDLC, Mean estimates at checkpoints.
 * Verifies ACDLL5 mean >= CDLL5 mean (more information retained).
 */
public class TestACDLL5Accuracy{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int buckets=2048;
		final CompressedDynamicLogLog5 cdll5=new CompressedDynamicLogLog5(buckets, 31, 42, 0);
		final ArithmeticCompressedDynamicLogLog5 acdll5=new ArithmeticCompressedDynamicLogLog5(buckets, 31, 42, 0);

		final long[] checkpoints={1000, 5000, 10000, 50000, 100000, 200000, 500000,
			1000000, 2000000, 5000000, 10000000};

		// Track cumulative state changes and overflows
		long cStateChanges=0, aStateChanges=0;

		// Previous bucket states for detecting state changes
		final int[] cPrevExp=new int[buckets];
		final int[] cPrevHist=new int[buckets];
		final int[] aPrevExp=new int[buckets];
		final int[] aPrevHist=new int[buckets];

		System.out.printf("%-12s  %-6s %-6s  %-8s %-8s  %-12s %-12s  %-12s %-12s  %-12s %-12s  %-12s %-12s  %-5s%n",
			"Adds", "C_chg", "A_chg", "C_ovfl", "A_ovfl",
			"C_DLC", "A_DLC", "C_LDLC", "A_LDLC",
			"C_Mean", "A_Mean", "C_Hybrid", "A_Hybrid", "A>=C?");
		System.out.println("=".repeat(170));

		int cpIdx=0;
		for(long i=0; i<checkpoints[checkpoints.length-1]; i++){
			// Snapshot before
			for(int b=0; b<buckets; b++){
				cPrevExp[b]=(cdll5.readBucket(b)>>>2)&0x7;
				cPrevHist[b]=cdll5.readBucket(b)&0x3;
				aPrevExp[b]=acdll5.readExponent(b);
				aPrevHist[b]=acdll5.readHistory(b);
			}

			cdll5.hashAndStore(i);
			acdll5.hashAndStore(i);

			// Detect state changes
			for(int b=0; b<buckets; b++){
				final int cExp=(cdll5.readBucket(b)>>>2)&0x7;
				final int cHist=cdll5.readBucket(b)&0x3;
				final int aExp=acdll5.readExponent(b);
				final int aHist=acdll5.readHistory(b);

				if(cExp!=cPrevExp[b] || cHist!=cPrevHist[b]){cStateChanges++;}
				if(aExp!=aPrevExp[b] || aHist!=aPrevHist[b]){aStateChanges++;}
			}

			if(cpIdx<checkpoints.length && i+1==checkpoints[cpIdx]){
				// Count overflow buckets: where ACDLL5 exp > 7 (CDLL5's max)
				int cAtMax=0, aAtMax=0;
				int bucketDivergences=0;
				for(int b=0; b<buckets; b++){
					final int cExp=(cdll5.readBucket(b)>>>2)&0x7;
					final int aExp=acdll5.readExponent(b);
					if(cExp==7){cAtMax++;}
					if(aExp==9){aAtMax++;}
					if(cExp!=aExp || (cdll5.readBucket(b)&0x3)!=acdll5.readHistory(b)){
						bucketDivergences++;
					}
				}

				// Get estimates
				cdll5.lastCardinality=-1;
				acdll5.lastCardinality=-1;
				final double[] cRaw=cdll5.rawEstimates();
				final CardStats cStats=cdll5.consumeLastSummarized();
				final double[] aRaw=acdll5.rawEstimates();
				final CardStats aStats=acdll5.consumeLastSummarized();

				final double cDlc=cStats.dlcRaw();
				final double aDlc=aStats.dlcRaw();
				final double cLdlc=cStats.ldlc();
				final double aLdlc=aStats.ldlc();
				final double cMean=cStats.meanCF();
				final double aMean=aStats.meanCF();
				final double cHyb=cStats.hybridDLL();
				final double aHyb=aStats.hybridDLL();

				final long trueCard=i+1;
				final boolean aMeanGeCMean=aMean>=cMean-0.01; // tolerance for float

				System.out.printf("%-12d  %-6d %-6d  %-8d %-8d  %12.2f %12.2f  %12.2f %12.2f  %12.2f %12.2f  %12.2f %12.2f  %-5s%n",
					trueCard, cStateChanges, aStateChanges, cAtMax, aAtMax,
					cDlc, aDlc, cLdlc, aLdlc, cMean, aMean, cHyb, aHyb,
					aMeanGeCMean ? "YES" : "NO!");

				// Print relative errors for key estimators
				System.out.printf("  errors:  DLC: C=%+.4f A=%+.4f  LDLC: C=%+.4f A=%+.4f  Mean: C=%+.4f A=%+.4f  divBuckets=%d minZ: C=%d A=%d%n",
					cDlc/trueCard-1, aDlc/trueCard-1,
					cLdlc/trueCard-1, aLdlc/trueCard-1,
					cMean/trueCard-1, aMean/trueCard-1,
					bucketDivergences, cdll5.getMinZeros(), acdll5.getMinZeros());

				cpIdx++;
			}
		}
	}

}
