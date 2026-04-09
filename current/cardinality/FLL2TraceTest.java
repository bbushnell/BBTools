package cardinality;

import rand.FastRandomXoshiro;

/**
 * Trace a single FLL2 word: print state, CF lookup, and compare with
 * what the simulator would compute.
 */
public class FLL2TraceTest {

	public static void main(String[] args) {
		FutureLogLog2.loadCFTable();

		long seed = 42;
		FutureLogLog2 fll = new FutureLogLog2(6, 31, seed, 0);
		FastRandomXoshiro rng = new FastRandomXoshiro(seed + 1000000);

		System.err.printf("%-8s %-4s %-4s %-6s %-6s %-12s %-12s %-12s %-12s %-10s%n",
			"TrueC", "GE", "LE", "Hist", "Idx84", "TierAvg", "Mult", "CF_Est", "FLL2_Est", "Error%");

		for (long trueCard = 1; trueCard <= 500; trueCard++) {
			fll.add(rng.nextLong());

			int w = fll.getWord(0);
			int ge = fll.getGlobalExp();
			int le = (w >>> 12) & 0xF;
			int hist = w & 0xFFF;
			int absNLZ = ge + le;
			int tier = Math.min(absNLZ, FutureLogLog2.CF_TABLE_TIERS - 1);
			int eq = FutureLogLog2.idx84(hist);

			double ta = (absNLZ < FutureLogLog2.tierAvg.length) ?
				FutureLogLog2.tierAvg[absNLZ] :
				FutureLogLog2.tierAvg[FutureLogLog2.tierAvg.length-1] *
				Math.pow(2.0, absNLZ - FutureLogLog2.tierAvg.length + 1);
			double mult = FutureLogLog2.cfTable[tier][eq];
			double cfEst = (w == 0) ? 0 : ta * mult;
			long fllEst = fll.cardinality();
			double err = (fllEst - trueCard) / (double) trueCard;

			boolean print = trueCard <= 40 || trueCard % 50 == 0;
			if (print) {
				System.err.printf("%-8d %-4d %-4d 0x%03X  %-6d %-12.2f %-12.6f %-12.2f %-12d %-10.4f%n",
					trueCard, ge, le, hist, eq, ta, mult, cfEst, fllEst, err);
			}
		}
	}
}
