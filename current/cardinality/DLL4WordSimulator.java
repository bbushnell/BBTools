package cardinality;

import fileIO.ByteStreamWriter;
import rand.FastRandomXoshiro;
import structures.ByteBuilder;
import structures.IntHashMap;

/**
 * DLL4 Word-Based Simulator.
 *
 * Simulates a single 4-bucket DLL4 word to build per-tier state tables.
 * Uses a char[65536] remap built at init time: for each of the 65536 possible
 * raw 16-bit register states, canonicalize (subtract min, sort, pack) and map
 * to a compact index. Accumulates into small arrays indexed by these compact indices.
 *
 * Run: java -ea cardinality.DLL4WordSimulator [iters=N] [threads=N] [maxTier=N]
 *
 * @author Brian Bushnell, Chloe von Einzbern-Bushnell
 * @date April 9, 2026
 */
public class DLL4WordSimulator {

	static final int BUCKETS_PER_WORD = 4;

	// --- Remap: raw 16-bit state -> compact canonical index ---
	static final char[] REMAP = new char[65536];
	static final int NUM_CANONICAL;
	/** Inverse: canonical index -> canonical packed key (for printing). */
	static final int[] CANONICAL_KEYS;

	static {
		IntHashMap map = new IntHashMap(1024);
		for (int state = 0; state < 65536; state++) {
			int r0 = state & 0xF;
			int r1 = (state >> 4) & 0xF;
			int r2 = (state >> 8) & 0xF;
			int r3 = (state >> 12) & 0xF;
			int wMin = Math.min(Math.min(r0, r1), Math.min(r2, r3));
			int d0 = r0 - wMin, d1 = r1 - wMin, d2 = r2 - wMin, d3 = r3 - wMin;
			// Sort ascending (5-comparator network)
			if (d0 > d1) { int t = d0; d0 = d1; d1 = t; }
			if (d2 > d3) { int t = d2; d2 = d3; d3 = t; }
			if (d0 > d2) { int t = d0; d0 = d2; d2 = t; }
			if (d1 > d3) { int t = d1; d1 = d3; d3 = t; }
			if (d1 > d2) { int t = d1; d1 = d2; d2 = t; }
			int canonKey = d0 | (d1 << 4) | (d2 << 8) | (d3 << 12);
			assert d0 == 0 : "Canonical form must start with 0";
			if (!map.containsKey(canonKey)) {
				map.put(canonKey, map.size());
			}
			REMAP[state] = (char) map.get(canonKey);
		}
		NUM_CANONICAL = map.size();
		CANONICAL_KEYS = new int[NUM_CANONICAL];
		int[] keys = map.keys();
		for (int k : keys) {
			if (k < 0) continue; // IntHashMap sentinel
			CANONICAL_KEYS[map.get(k)] = k;
		}
	}

	static String keyToString(int key) {
		return "{" + (key & 0xF) + "," + ((key >> 4) & 0xF) + ","
			+ ((key >> 8) & 0xF) + "," + ((key >> 12) & 0xF) + "}";
	}

	// --- Simulation parameters ---
	final int iters;
	final int threads;
	final int maxTier;
	final int tiers;

	// Accumulated statistics [tier][canonicalIdx]
	long[][]   counts;
	double[][] sums;
	double[][] sumSq;

	long belowFloorCount;
	long totalAdds;
	long advanceCount;

	public DLL4WordSimulator(int iters, int threads, int maxTier) {
		this.iters   = iters;
		this.threads = threads;
		this.maxTier = maxTier;
		this.tiers   = maxTier + 2;
		this.counts  = new long[tiers][NUM_CANONICAL];
		this.sums    = new double[tiers][NUM_CANONICAL];
		this.sumSq   = new double[tiers][NUM_CANONICAL];
	}

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

	void simulate() throws InterruptedException {
		long[][][]   threadCounts = new long[threads][tiers][NUM_CANONICAL];
		double[][][] threadSums   = new double[threads][tiers][NUM_CANONICAL];
		double[][][] threadSumSq  = new double[threads][tiers][NUM_CANONICAL];
		long[] threadBelow   = new long[threads];
		long[] threadTotal   = new long[threads];
		long[] threadAdvance = new long[threads];

		int base = iters / threads;
		int rem  = iters % threads;

		Thread[] threadArr = new Thread[threads];
		for (int t = 0; t < threads; t++) {
			final int tid     = t;
			final int myIters = base + (tid < rem ? 1 : 0);

			threadArr[t] = new Thread(() -> {
				FastRandomXoshiro rng = new FastRandomXoshiro(tid + 1);
				long[][] lCounts   = threadCounts[tid];
				double[][] lSums   = threadSums[tid];
				double[][] lSumSq  = threadSumSq[tid];
				long lBelow = 0, lTotal = 0, lAdvance = 0;

				for (int iter = 0; iter < myIters; iter++) {
					int[] reg = new int[4];
					int globalExp = 0;
					int minZeroCount = 4;
					long eeMask = -1L;

					for (long trueCard = 1; trueCard <= 10_000_000L; trueCard++) {
						long hash = hash64shift(rng.nextLong());
						lTotal++;

						if (Long.compareUnsigned(hash, eeMask) > 0) {
							lBelow++;
						} else {
							int hashNlz = Long.numberOfLeadingZeros(hash);
							int bucket = (int)((hash & 0x7FFF_FFFFL) % BUCKETS_PER_WORD);
							int relNlz = hashNlz - globalExp;

							if (relNlz < 0) {
								lBelow++;
							} else {
								int newStored = Math.min(relNlz + 1, 15);
								int oldStored = reg[bucket];

								if (newStored > oldStored) {
									reg[bucket] = newStored;

									if (oldStored == 0) {
										minZeroCount--;
										while (minZeroCount == 0 && globalExp < 64) {
											globalExp++;
											eeMask >>>= 1;
											minZeroCount = 0;
											for (int i = 0; i < 4; i++) {
												reg[i] = Math.max(0, reg[i] - 1);
												if (reg[i] == 0) minZeroCount++;
											}
											lAdvance++;
										}
									}
								}
							}
						}

						// Record using remap
						int wMin = Math.min(Math.min(reg[0], reg[1]), Math.min(reg[2], reg[3]));
						int tier = globalExp + wMin;
						int rawKey = reg[0] | (reg[1] << 4) | (reg[2] << 8) | (reg[3] << 12);
						int idx = REMAP[rawKey];

						if (tier < tiers) {
							lCounts[tier][idx]++;
							lSums[tier][idx] += trueCard;
							lSumSq[tier][idx] += (double) trueCard * trueCard;
						}

						if (tier > maxTier) break;
					}
				}

				threadBelow[tid]  = lBelow;
				threadTotal[tid]  = lTotal;
				threadAdvance[tid] = lAdvance;
			});
			threadArr[t].start();
		}

		for (Thread t : threadArr) t.join();

		for (int t = 0; t < threads; t++) {
			for (int tier = 0; tier < tiers; tier++) {
				for (int s = 0; s < NUM_CANONICAL; s++) {
					counts[tier][s] += threadCounts[t][tier][s];
					sums[tier][s]   += threadSums[t][tier][s];
					sumSq[tier][s]  += threadSumSq[t][tier][s];
				}
			}
			belowFloorCount += threadBelow[t];
			totalAdds       += threadTotal[t];
			advanceCount    += threadAdvance[t];
		}
	}

	double[][] buildStateAvg() {
		double[][] avg = new double[tiers][NUM_CANONICAL];
		for (int tier = 0; tier < tiers; tier++) {
			for (int s = 0; s < NUM_CANONICAL; s++) {
				avg[tier][s] = counts[tier][s] > 0 ? sums[tier][s] / counts[tier][s] : 0;
			}
		}
		return avg;
	}

	double[][] buildCV() {
		double[][] cv = new double[tiers][NUM_CANONICAL];
		for (int tier = 0; tier < tiers; tier++) {
			for (int s = 0; s < NUM_CANONICAL; s++) {
				cv[tier][s] = 1.0;
				if (counts[tier][s] > 1) {
					double mean = sums[tier][s] / counts[tier][s];
					if (mean > 0) {
						double variance = sumSq[tier][s] / counts[tier][s] - mean * mean;
						if (variance > 0) {
							cv[tier][s] = Math.sqrt(variance) / mean;
						} else {
							cv[tier][s] = 0.001;
						}
					}
				}
			}
		}
		return cv;
	}

	void printResults(double[][] stateAvg, double[][] cv) {
		System.err.println("=== DLL4 Word Simulator Results ===");
		System.err.println("iters=" + iters + "  threads=" + threads + "  maxTier=" + maxTier);
		System.err.println("numCanonical=" + NUM_CANONICAL);
		System.err.println("totalAdds=" + totalAdds);
		System.err.printf("belowFloorRate=%.6f%n",
			totalAdds > 0 ? (double) belowFloorCount / totalAdds : 0);
		System.err.printf("advanceCount  =%d%n", advanceCount);
		System.err.println();

		System.err.printf("%-6s %-16s %-10s %-12s %-10s%n", "Tier", "AvgCard", "Growth", "Observations", "States");
		double prevAvg = 0;
		for (int tier = 0; tier < tiers; tier++) {
			long obs = 0;
			double totalSum = 0;
			int statesUsed = 0;
			for (int s = 0; s < NUM_CANONICAL; s++) {
				obs += counts[tier][s];
				totalSum += sums[tier][s];
				if (counts[tier][s] > 0) statesUsed++;
			}
			if (obs == 0) continue;
			double avg = totalSum / obs;
			double growth = (prevAvg > 0) ? avg / prevAvg : Double.NaN;
			System.err.printf("%-6d %-16.2f %-10.4f %-12d %-10d%n",
				tier, avg, growth, obs, statesUsed);
			prevAvg = avg;
		}
		System.err.println();

		// Sparse output: tier \t idx \t canonKey \t stateAvg \t count \t cv
		ByteStreamWriter bsw = new ByteStreamWriter("stdout", true, false, false);
		bsw.start();
		bsw.println("#tier\tidx\tcanonKey\tstateAvg\tcount\tcv");
		ByteBuilder bb = new ByteBuilder(256);
		for (int tier = 0; tier < tiers; tier++) {
			for (int idx = 0; idx < NUM_CANONICAL; idx++) {
				if (counts[tier][idx] > 0) {
					bb.clear();
					bb.append(tier).append('\t').append(idx).append('\t');
					bb.append(CANONICAL_KEYS[idx]).append('\t');
					bb.append(String.valueOf(stateAvg[tier][idx])).append('\t');
					bb.append(counts[tier][idx]).append('\t');
					bb.append(String.valueOf(cv[tier][idx]));
					bsw.println(bb);
				}
			}
		}
		bsw.poisonAndWait();
	}

	public static void main(String[] args) throws InterruptedException {
		int iters   = 10000;
		int threads = 8;
		int maxTier = 12;

		for (String arg : args) {
			String[] kv = arg.split("=", 2);
			if (kv.length != 2) continue;
			switch (kv[0]) {
				case "iters":   iters   = Integer.parseInt(kv[1]);  break;
				case "threads": threads = Integer.parseInt(kv[1]);  break;
				case "maxTier": maxTier = Integer.parseInt(kv[1]);  break;
				default: System.err.println("Unknown arg: " + arg);
			}
		}

		System.err.println("DLL4WordSimulator: iters=" + iters + " threads=" + threads
			+ " maxTier=" + maxTier + " numCanonical=" + NUM_CANONICAL);

		long t0 = System.currentTimeMillis();
		DLL4WordSimulator sim = new DLL4WordSimulator(iters, threads, maxTier);
		sim.simulate();
		long t1 = System.currentTimeMillis();
		System.err.println("Simulation time: " + (t1 - t0) + " ms");

		double[][] stateAvg = sim.buildStateAvg();
		double[][] cv = sim.buildCV();
		sim.printResults(stateAvg, cv);
	}
}
