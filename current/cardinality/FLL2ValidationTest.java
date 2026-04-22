package cardinality;

import rand.FastRandomXoshiro;

/**
 * Direct validation: creates a real FutureLogLog2 with 6 buckets (1 word),
 * feeds it random values, and checks estimate vs truth at each step.
 * Uses the exact same code path as the real estimator.
 * Reports per-tier signed and absolute error across many iterations.
 *
 * @author Brian Bushnell
 */
public class FLL2ValidationTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int iters=10000;
		int maxCard=10_000_000;

		for(final String arg : args){
			final String[] kv=arg.split("=", 2);
			if(kv.length!=2){continue;}
			if("iters".equals(kv[0])){iters=Integer.parseInt(kv[1]);}
			if("maxcard".equals(kv[0])){maxCard=Integer.parseInt(kv[1]);}
		}

		FutureLogLog2.loadCFTable();
		CardinalityTracker.clampToAdded=false;

		double sumSigned=0, sumAbs=0;
		long totalObs=0;

		// Per-tier tracking
		final double[] tierSigned=new double[16];
		final double[] tierAbs=new double[16];
		final long[] tierObs=new long[16];

		for(int iter=0; iter<iters; iter++){
			final long seed=iter+1;
			final FutureLogLog2 fll=new FutureLogLog2(6, 31, seed, 0);
			final FastRandomXoshiro rng=new FastRandomXoshiro(seed+1000000);

			for(long trueCard=1; trueCard<=maxCard; trueCard++){
				fll.add(rng.nextLong());

				final long est=fll.cardinality();
				final double err=(est-trueCard)/(double)trueCard;
				sumSigned+=err;
				sumAbs+=Math.abs(err);
				totalObs++;

				final int tier=Math.min(fll.getGlobalExp(), 15);
				tierSigned[tier]+=err;
				tierAbs[tier]+=Math.abs(err);
				tierObs[tier]++;

				if(fll.getGlobalExp()>12){break;}
			}
		}

		System.err.println("=== FLL2 Real Estimator Validation (1 word, 6 buckets) ===");
		System.err.println("iters="+iters+"  maxCard="+maxCard);
		System.err.printf("Overall: signedErr=%.6f  absErr=%.6f  observations=%d%n",
			sumSigned/totalObs, sumAbs/totalObs, totalObs);
		System.err.printf("%-6s %-14s %-14s %-12s%n", "Tier", "SignedErr", "AbsErr", "Observations");
		for(int t=0; t<16; t++){
			if(tierObs[t]==0){continue;}
			System.err.printf("%-6d %-14.6f %-14.6f %-12d%n",
				t, tierSigned[t]/tierObs[t], tierAbs[t]/tierObs[t], tierObs[t]);
		}
	}
}
