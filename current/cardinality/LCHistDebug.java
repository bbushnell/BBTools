package cardinality;

/**
 * Debug: 1 estimator, 8 buckets, print state after every add.
 */
public class LCHistDebug {

	public static void main(String[] args) {
		// Set up PLL16c with 8 buckets, 2-bit history
		ProtoLogLog16c.HISTORY_BITS = 2;
		ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_HISTORY);
		CorrectionFactor.sbsFile = "?cardinalityCorrectionLC2BitHist_8b.tsv.gz";
		CorrectionFactor.loadSbsTable();
		System.err.println("LC_2BIT_CF_TABLE: " +
			(CorrectionFactor.SBS_CF_TABLE != null ? "loaded, buckets=" + CorrectionFactor.sbsBuckets : "NULL"));

		int buckets = 8;
		ProtoLogLog16c pll = new ProtoLogLog16c(buckets, 31, -1, 0);

		// Print table row for filled=1..8 so we can cross-check
		System.out.println("=== LC History Table (rows 1-8) ===");
		String[] stateNames = {"0.00","1.00","1.10","2.00","2.01","2.10","2.11",
			"3+.00","3+.01","3+.10","3+.11"};
		System.out.print("Filled");
		for (String n : stateNames) System.out.print("\t" + n);
		System.out.println();
		float[][] table = CorrectionFactor.SBS_CF_TABLE;
		if (table != null) {
			for (int f = 1; f <= Math.min(table.length - 1, 8); f++) {
				System.out.print(f);
				for (int s = 0; s < 11; s++) {
					System.out.printf("\t%.4f", table[f][s]);
				}
				System.out.println();
			}
		}
		System.out.println();

		// Add elements and print state after each
		System.out.println("=== Simulation ===");
		System.out.println("Add#\tBucket\tRawNLZ\tAction\tFilled\tLC\tSBS\tBucket_States");

		java.util.Random rng = new java.util.Random(42);
		for (int i = 0; i < 30; i++) {
			long key = rng.nextLong();
			pll.add(key);

			// Get estimates
			// Access the internal state
			short[] maxArray = pll.maxArray;
			int filled = 0;
			StringBuilder states = new StringBuilder();
			for (int b = 0; b < buckets; b++) {
				int stored = maxArray[b] & 0xFFFF;
				if (stored > 0) {
					filled++;
					int nlzShift = 16 - 6; // 6 NLZ bits for PLL16c
					int absNlz = stored >>> nlzShift;
					int extra = stored & ((1 << nlzShift) - 1);
					int hshift = nlzShift - 2; // history shift (2 bits at top of extra)
					int hist = (extra >>> hshift) & 3;
					int rawNlz = absNlz - 1;
					int nlzBin = Math.min(rawNlz, 3);
					int si = CorrectionFactor.sbsStateIndex(nlzBin, hist);
					states.append(String.format("[%d:nlz%d.%02d si=%d]",
						b, rawNlz, Integer.parseInt(Integer.toBinaryString(hist)), si));
				} else {
					states.append(String.format("[%d:empty]", b));
				}
			}

			// Compute LC
			int V = buckets - filled;
			double lc = V > 0 ? buckets * Math.log((double) buckets / V) : Double.POSITIVE_INFINITY;

			// Compute SBS manually — direct lookup, B matches table
			double sbs = 0;
			if (table != null && filled > 0) {
				float[] row = table[Math.min(filled, table.length - 1)];

				for (int b = 0; b < buckets; b++) {
					int stored = maxArray[b] & 0xFFFF;
					if (stored > 0) {
						int nlzShift = 16 - 6;
						int absNlz = stored >>> nlzShift;
						int extra = stored & ((1 << nlzShift) - 1);
						int hshift = nlzShift - 2;
						int hist = (extra >>> hshift) & 3;
						int rawNlz = absNlz - 1;
						int nlzBin = Math.min(rawNlz, 3);
						int si = CorrectionFactor.sbsStateIndex(nlzBin, hist);
						if (si >= 0) sbs += row[si];
					}
				}
			}

			// What bucket did this key hash to?
			// We can't easily know, so just print the add number
			System.out.printf("%d\t-\t-\t-\t%d\t%.4f\t%.4f\t%s%n",
				i + 1, filled, lc, sbs, states);
		}
	}
}
