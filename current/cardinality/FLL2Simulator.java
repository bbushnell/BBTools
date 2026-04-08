package cardinality;

import java.util.Arrays;
import java.util.HashSet;

/**
 * FLL2 (FutureLogLog 2-bit) Simulator.
 *
 * Simulates a single 16-bit FLL2 word to build per-tier correction factor tables.
 * The word encodes one bank of 6 buckets sharing a 4-bit local exponent.
 *
 * Word layout: [15:12]=localExp [11:0]=6×2-bit future bitmap
 * Per-bucket state (2 bits):
 *   Bit 0 (LSB) = "seen hash with NLZ == absNLZ"     (floor hit)
 *   Bit 1 (MSB) = "seen hash with NLZ == absNLZ+1"   (future hit)
 *   All four states {00,01,10,11} are reachable.
 *
 * Promotion fires when all 6 LSBs are set: (word & LSB_MASK) == LSB_MASK.
 * Promotion: MSBs shift right to LSB positions, MSBs clear, localExp++.
 * Max cascade = 2 (all-{11} → all-{01} → all-{00}).
 *
 * Output: double[16][4096] multiplier table and double[16] tier averages.
 * Equivalent states (same bucket histogram) are combined before output.
 *
 * Run with: java -ea cardinality.FLL2Simulator [iters=N] [threads=N] [maxTier=N]
 *
 * @author Brian Bushnell, Chloe von Einzbern-Bushnell
 * @date April 8, 2026
 */
public class FLL2Simulator {

	// --- Word layout constants ---
	static final int BUCKETS_PER_WORD = 6;
	static final int NUM_TIERS        = 16;
	static final int HISTORY_BITS     = 12;
	static final int HISTORY_STATES   = 1 << HISTORY_BITS; // 4096
	/** Bits 0,2,4,6,8,10 — floor-level bits (LSB of each bucket) */
	static final int LSB_MASK = 0x555;
	/** Bits 1,3,5,7,9,11 — future-level bits (MSB of each bucket) */
	static final int MSB_MASK = 0xAAA;

	// --- Simulation parameters ---
	final int iters;
	final int threads;
	final int maxTier; // end trial when word's tier exceeds this (0-15)

	// --- Accumulated statistics [tier][history] ---
	long[][]   counts;
	double[][] sums;

	// --- Diagnostic counters ---
	long[] promotionHist; // [n] = number of adds causing exactly n promotions
	long overflowCount;   // adds where delta > 1 (ignored, IOT mode)
	long underflowCount;  // adds where delta < 0 (ignored)
	long totalAdds;

	public FLL2Simulator(int iters, int threads, int maxTier) {
		assert iters > 0      : "iters must be positive: " + iters;
		assert threads > 0    : "threads must be positive: " + threads;
		assert maxTier >= 0 && maxTier < NUM_TIERS : "maxTier out of range: " + maxTier;
		this.iters   = iters;
		this.threads = threads;
		this.maxTier = maxTier;
		this.counts  = new long[NUM_TIERS][HISTORY_STATES];
		this.sums    = new double[NUM_TIERS][HISTORY_STATES];
		this.promotionHist = new long[4]; // 0, 1, 2, 3+
	}

	// ---------------------------------------------------------------------------
	// Hash function (matches Tools.hash64shift in BBTools)
	// ---------------------------------------------------------------------------

	static long hash64shift(long key) {
		key = (~key) + (key << 21);
		key ^= (key >>> 24);
		key += (key << 3) + (key << 8);
		key ^= (key >>> 14);
		key += (key << 2) + (key << 4);
		key ^= (key >>> 28);
		key += (key << 31);
		return key;
	}

	// ---------------------------------------------------------------------------
	// Canonicalize: order-independent history state
	// Extract 6 two-bit bucket values, sort ascending, re-pack.
	// All permutations of the same bucket histogram map to the same canonical value.
	// ---------------------------------------------------------------------------

	static int canonicalize(int history) {
		int[] vals = new int[BUCKETS_PER_WORD];
		for (int i = 0; i < BUCKETS_PER_WORD; i++) {
			vals[i] = (history >>> (i * 2)) & 0x3;
		}
		Arrays.sort(vals); // ascending
		int canonical = 0;
		for (int i = BUCKETS_PER_WORD - 1; i >= 0; i--) {
			canonical = (canonical << 2) | vals[i];
		}
		return canonical;
	}

	// ---------------------------------------------------------------------------
	// Promote: shift MSBs to LSBs, increment localExp.
	// Cascades until no promotion needed or localExp==15.
	// Returns updated word; increments promHist[promotionCount].
	// ---------------------------------------------------------------------------

	static int promote(int word, long[] promHist) {
		int promotionCount = 0;
		while (true) {
			if ((word & LSB_MASK) != LSB_MASK) break; // not all LSBs set
			int localExp = (word >>> 12) & 0xF;
			if (localExp >= 15) break; // can't promote past tier 15
			// MSBs → LSBs, clear MSBs, increment localExp
			word = ((word & MSB_MASK) >>> 1) | ((localExp + 1) << 12);
			promotionCount++;
			assert promotionCount <= 2
				: "Cascade exceeded 2 — impossible with 2-bit future: count=" + promotionCount;
			assert (word & 0xFFFF) == word : "Word overflow after promotion: " + word;
		}
		if (promHist != null) {
			promHist[Math.min(promotionCount, promHist.length - 1)]++;
		}
		return word;
	}

	// ---------------------------------------------------------------------------
	// Invariant checker for a single word state.
	// ---------------------------------------------------------------------------

	static boolean repOK(int word) {
		int localExp = (word >>> 12) & 0xF;
		assert localExp >= 0 && localExp <= 15
			: "localExp out of range: " + localExp;
		assert (word & 0xFFFF) == word
			: "Word exceeds 16 bits: " + Integer.toHexString(word);
		// No promotable state should survive (unless maxed out)
		assert localExp >= 15 || (word & LSB_MASK) != LSB_MASK
			: "Promotable state should have been promoted: word=0x"
			+ Integer.toHexString(word) + " localExp=" + localExp;
		return true;
	}

	// ---------------------------------------------------------------------------
	// Core simulation
	// ---------------------------------------------------------------------------

	void simulate() throws InterruptedException {
		// Per-thread result storage
		long[][][]   threadCounts    = new long[threads][NUM_TIERS][HISTORY_STATES];
		double[][][] threadSums      = new double[threads][NUM_TIERS][HISTORY_STATES];
		long[][]     threadPromHist  = new long[threads][4];
		long[]       threadOverflow  = new long[threads];
		long[]       threadUnderflow = new long[threads];
		long[]       threadTotal     = new long[threads];

		int base = iters / threads;
		int rem  = iters % threads;

		Thread[] threadArr = new Thread[threads];
		for (int t = 0; t < threads; t++) {
			final int tid     = t;
			final int myIters = base + (tid < rem ? 1 : 0);

			threadArr[t] = new Thread(() -> {
				java.util.Random rng = new java.util.Random(tid * 6364136223846793005L + 1442695040888963407L);
				long[][] lCounts   = threadCounts[tid];
				double[][] lSums   = threadSums[tid];
				long[] lPromHist   = threadPromHist[tid];
				long lOverflow = 0, lUnderflow = 0, lTotal = 0;

				for (int iter = 0; iter < myIters; iter++) {
					int word = 0;
					long eeMask = -1L; // (-1L) >>> localExp; all hashes valid at tier 0

					for (long trueCard = 1; trueCard <= 10_000_000L; trueCard++) {
						long hash = hash64shift(rng.nextLong());
						lTotal++;

						// Fast skip: hashes above eeMask have NLZ < localExp (below floor).
						// Still counts as an added element (trueCard advances).
						if (Long.compareUnsigned(hash, eeMask) > 0) {
							lUnderflow++;
						} else {
							int hashNLZ = Long.numberOfLeadingZeros(hash);
							int tier    = (word >>> 12) & 0xF;
							int delta   = hashNLZ - tier; // >= 0 guaranteed by eeMask

							if (delta > 1) {
								lOverflow++; // IOT: ignore
							} else {
								// Set the appropriate bit: delta=0 → LSB, delta=1 → MSB
								int bucket   = (int)((hash & 0x7FFF_FFFFL) % BUCKETS_PER_WORD);
								int bitToSet = 1 << (delta + bucket * 2);
								if ((word & bitToSet) == 0) {
									word |= bitToSet;
									if ((word & LSB_MASK) == LSB_MASK) {
										word = promote(word, lPromHist);
										// Update eeMask after promotion
										eeMask = -1L >>> ((word >>> 12) & 0xF);
									} else {
										lPromHist[0]++;
									}
								} else {
									lPromHist[0]++;
								}
							}
							assert repOK(word);
						}

						// Record state after every element (including underflows)
						int curTier = (word >>> 12) & 0xF;
						int history = word & 0xFFF;
						lCounts[curTier][history]++;
						lSums[curTier][history] += trueCard;

						// End trial when word advances past maxTier
						if (curTier > maxTier) break;
					}
				}

				threadOverflow[tid]  = lOverflow;
				threadUnderflow[tid] = lUnderflow;
				threadTotal[tid]     = lTotal;
			});
			threadArr[t].start();
		}

		for (Thread t : threadArr) t.join();

		// Combine per-thread results
		for (int t = 0; t < threads; t++) {
			for (int tier = 0; tier < NUM_TIERS; tier++) {
				for (int h = 0; h < HISTORY_STATES; h++) {
					counts[tier][h] += threadCounts[t][tier][h];
					sums[tier][h]   += threadSums[t][tier][h];
				}
			}
			for (int i = 0; i < 4; i++) promotionHist[i] += threadPromHist[t][i];
			overflowCount  += threadOverflow[t];
			underflowCount += threadUnderflow[t];
			totalAdds      += threadTotal[t];
		}
	}

	// ---------------------------------------------------------------------------
	// Combine equivalent states (same bucket histogram → same multiplier)
	// ---------------------------------------------------------------------------

	void combineEquivalentStates() {
		long[][]   combined_counts = new long[NUM_TIERS][HISTORY_STATES];
		double[][] combined_sums   = new double[NUM_TIERS][HISTORY_STATES];

		// Accumulate into canonical forms
		for (int tier = 0; tier < NUM_TIERS; tier++) {
			for (int h = 0; h < HISTORY_STATES; h++) {
				if (counts[tier][h] > 0) {
					int canon = canonicalize(h);
					combined_counts[tier][canon] += counts[tier][h];
					combined_sums[tier][canon]   += sums[tier][h];
				}
			}
		}

		// Expand back: every state gets its canonical's combined stats
		for (int tier = 0; tier < NUM_TIERS; tier++) {
			for (int h = 0; h < HISTORY_STATES; h++) {
				int canon = canonicalize(h);
				counts[tier][h] = combined_counts[tier][canon];
				sums[tier][h]   = combined_sums[tier][canon];
			}
		}
	}

	// ---------------------------------------------------------------------------
	// Build output tables
	// ---------------------------------------------------------------------------

	/** Average cardinality for each tier, summed across all history states. */
	double[] buildTierAverages() {
		double[] tierAvg = new double[NUM_TIERS];
		for (int tier = 0; tier < NUM_TIERS; tier++) {
			long   totalCount = 0;
			double totalSum   = 0;
			for (int h = 0; h < HISTORY_STATES; h++) {
				totalCount += counts[tier][h];
				totalSum   += sums[tier][h];
			}
			tierAvg[tier] = totalCount > 0 ? totalSum / totalCount : 0;
		}
		return tierAvg;
	}

	/**
	 * Per-state multiplier: stateAvg / tierAvg.
	 * Gives a value close to 1.0 for "typical" states at each tier.
	 * States with no observations get multiplier 1.0 (neutral).
	 */
	double[][] buildMultipliers(double[] tierAvg) {
		double[][] mult = new double[NUM_TIERS][HISTORY_STATES];
		for (int tier = 0; tier < NUM_TIERS; tier++) {
			// Default to 1.0 for unobserved states
			Arrays.fill(mult[tier], 1.0);
			if (tierAvg[tier] <= 0) continue;
			for (int h = 0; h < HISTORY_STATES; h++) {
				if (counts[tier][h] > 0) {
					double stateAvg = sums[tier][h] / counts[tier][h];
					mult[tier][h] = stateAvg / tierAvg[tier];
				}
			}
		}
		return mult;
	}

	// ---------------------------------------------------------------------------
	// Print results
	// ---------------------------------------------------------------------------

	void printResults(double[] tierAvg, double[][] mult) {
		System.out.println("=== FLL2 Simulator Results ===");
		System.out.println("iters=" + iters + "  threads=" + threads + "  maxTier=" + maxTier);
		System.out.println("totalAdds=" + totalAdds);
		System.out.printf("overflowRate =%.4f  (delta>1, IOT-ignored)%n",
			totalAdds > 0 ? (double) overflowCount / totalAdds : 0);
		System.out.printf("underflowRate=%.4f  (delta<0, ignored)%n",
			totalAdds > 0 ? (double) underflowCount / totalAdds : 0);

		long totalProm = 0;
		for (long v : promotionHist) totalProm += v;
		System.out.printf("Promotions per add  0=%.4f  1=%.4f  2=%.4f  3+=%.4f%n",
			totalProm > 0 ? (double) promotionHist[0] / totalProm : 0,
			totalProm > 0 ? (double) promotionHist[1] / totalProm : 0,
			totalProm > 0 ? (double) promotionHist[2] / totalProm : 0,
			totalProm > 0 ? (double) promotionHist[3] / totalProm : 0);

		System.out.println("\n--- Per-Tier Average Cardinality ---");
		System.out.printf("%-6s %-16s %-10s %-12s%n", "Tier", "AvgCard", "Growth", "Observations");
		for (int tier = 0; tier <= maxTier; tier++) {
			long obs = 0;
			for (int h = 0; h < HISTORY_STATES; h++) obs += counts[tier][h];
			double growth = (tier > 0 && tierAvg[tier-1] > 0)
				? tierAvg[tier] / tierAvg[tier-1] : Double.NaN;
			System.out.printf("%-6d %-16.2f %-10.4f %-12d%n",
				tier, tierAvg[tier], growth, obs);
		}

		System.out.println("\n--- Distinct States Per Tier ---");
		for (int tier = 0; tier <= maxTier; tier++) {
			int rawDistinct = 0;
			HashSet<Integer> canonSet = new HashSet<>();
			for (int h = 0; h < HISTORY_STATES; h++) {
				if (counts[tier][h] > 0) {
					rawDistinct++;
					canonSet.add(canonicalize(h));
				}
			}
			System.out.printf("Tier %2d: %4d raw states,  %2d canonical equivalence classes%n",
				tier, rawDistinct, canonSet.size());
		}

		System.out.println("\n--- Sample Multipliers Tier 0 (canonical states only) ---");
		System.out.printf("%-16s %-30s %-12s %-12s%n", "History", "Buckets", "Multiplier", "Count");
		int printed = 0;
		for (int h = 0; h < HISTORY_STATES && printed < 30; h++) {
			// Only print canonical states (avoids duplicates)
			if (counts[0][h] > 0 && canonicalize(h) == h) {
				System.out.printf("0x%03X            %-30s %-12.6f %-12d%n",
					h, historyToString(h), mult[0][h], counts[0][h]);
				printed++;
			}
		}

		System.out.println("\n--- Min/Max Multiplier Per Tier ---");
		for (int tier = 0; tier <= maxTier; tier++) {
			double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;
			for (int h = 0; h < HISTORY_STATES; h++) {
				if (counts[tier][h] > 0) {
					min = Math.min(min, mult[tier][h]);
					max = Math.max(max, mult[tier][h]);
				}
			}
			System.out.printf("Tier %2d: mult range [%.4f, %.4f]%n", tier, min, max);
		}
	}

	/** Human-readable bucket state string, e.g. {11,01,00,00,01,10} */
	static String historyToString(int history) {
		StringBuilder sb = new StringBuilder("{");
		for (int i = BUCKETS_PER_WORD - 1; i >= 0; i--) {
			if (i < BUCKETS_PER_WORD - 1) sb.append(',');
			int v = (history >>> (i * 2)) & 0x3;
			sb.append((v >> 1) & 1).append(v & 1); // binary: "00","01","10","11"
		}
		sb.append('}');
		return sb.toString();
	}

	// ---------------------------------------------------------------------------
	// main
	// ---------------------------------------------------------------------------

	public static void main(String[] args) throws InterruptedException {
		int iters   = 10000;
		int threads = 8;
		int maxTier = 15;

		for (String arg : args) {
			String[] kv = arg.split("=", 2);
			if (kv.length != 2) continue;
			switch (kv[0]) {
				case "iters":   iters   = Integer.parseInt(kv[1]); break;
				case "threads": threads = Integer.parseInt(kv[1]); break;
				case "maxTier": maxTier = Integer.parseInt(kv[1]); break;
				default: System.err.println("Unknown arg: " + arg);
			}
		}

		System.out.println("FLL2Simulator: iters=" + iters
			+ " threads=" + threads + " maxTier=" + maxTier);

		long t0 = System.currentTimeMillis();
		FLL2Simulator sim = new FLL2Simulator(iters, threads, maxTier);
		sim.simulate();
		long t1 = System.currentTimeMillis();
		System.out.println("Simulation time: " + (t1 - t0) + " ms");

		sim.combineEquivalentStates();
		double[] tierAvg = sim.buildTierAverages();
		double[][] mult  = sim.buildMultipliers(tierAvg);
		sim.printResults(tierAvg, mult);
	}
}
