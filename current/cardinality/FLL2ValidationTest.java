package cardinality;

import rand.FastRandomXoshiro;

/**
 * Direct validation: create a real FutureLogLog2 with 6 buckets (1 word),
 * feed it random values, check estimate vs truth at each step.
 * This uses the EXACT same code path as the real estimator.
 */
public class FLL2ValidationTest {

	public static void main(String[] args) {
		int iters = 10000;
		int maxCard = 10_000_000;

		for (String arg : args) {
			String[] kv = arg.split("=", 2);
			if (kv.length != 2) continue;
			if ("iters".equals(kv[0])) iters = Integer.parseInt(kv[1]);
			if ("maxcard".equals(kv[0])) maxCard = Integer.parseInt(kv[1]);
		}

		FutureLogLog2.loadCFTable();
		CardinalityTracker.clampToAdded = false;

		double sumSigned = 0, sumAbs = 0;
		long totalObs = 0;

		// Per-tier tracking
		double[] tierSigned = new double[16];
		double[] tierAbs = new double[16];
		long[] tierObs = new long[16];

		for (int iter = 0; iter < iters; iter++) {
			long seed = iter + 1;
			FutureLogLog2 fll = new FutureLogLog2(6, 31, seed, 0);
			FastRandomXoshiro rng = new FastRandomXoshiro(seed + 1000000);

			for (long trueCard = 1; trueCard <= maxCard; trueCard++) {
				fll.add(rng.nextLong());

				long est = fll.cardinality();
				double err = (est - trueCard) / (double) trueCard;
				sumSigned += err;
				sumAbs += Math.abs(err);
				totalObs++;

				int tier = Math.min(fll.getGlobalExp(), 15);
				tierSigned[tier] += err;
				tierAbs[tier] += Math.abs(err);
				tierObs[tier]++;

				if (fll.getGlobalExp() > 12) break;
			}
		}

		System.err.println("=== FLL2 Real Estimator Validation (1 word, 6 buckets) ===");
		System.err.println("iters=" + iters + "  maxCard=" + maxCard);
		System.err.printf("Overall: signedErr=%.6f  absErr=%.6f  observations=%d%n",
			sumSigned / totalObs, sumAbs / totalObs, totalObs);
		System.err.printf("%-6s %-14s %-14s %-12s%n", "Tier", "SignedErr", "AbsErr", "Observations");
		for (int t = 0; t < 16; t++) {
			if (tierObs[t] == 0) continue;
			System.err.printf("%-6d %-14.6f %-14.6f %-12d%n",
				t, tierSigned[t] / tierObs[t], tierAbs[t] / tierObs[t], tierObs[t]);
		}
	}
}
