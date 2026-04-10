package cardinality;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import rand.FastRandomXoshiro;
import structures.ByteBuilder;
import structures.IntHashMap;

/**
 * DLL4 Word-Based Simulator.
 *
 * Simulates DLL4 words to build per-tier state tables and/or validate
 * estimation accuracy against a loaded table.
 *
 * Recording modes:
 *   every - record state at every cardinality increment (original behavior)
 *   entry - record only at state transitions (new state at trueCard)
 *   both  - record at transitions (new state at trueCard, old state at trueCard-1)
 *
 * Error tracking:
 *   Load a state table with table=file and set words=N to simulate a full
 *   N-word DLL4. At each transition, computes the running estimate (sum of
 *   stateAvg across all words) and records error vs trueCard.
 *
 * Run: java -ea cardinality.DLL4WordSimulator [iters=N] [threads=N]
 *      [maxTier=N] [words=N] [mode=every|entry|both] [table=file]
 *
 * @author Brian Bushnell, Chloe
 * @date April 10, 2026
 */
public class DLL4WordSimulator {

	static final int BUCKETS_PER_WORD = 4;

	static final int MODE_EVERY = 0;
	static final int MODE_ENTRY = 1;
	static final int MODE_BOTH = 2;

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
			if (k < 0) continue;
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
	final int words;
	final int mode;

	// Table-building accumulators [tier][canonicalIdx]
	long[][]   counts;
	double[][] sums;
	double[][] sumSq;

	// Loaded state table for error tracking (null = table-building mode)
	double[][] loadedTable;

	// Error tracking accumulators (global)
	double logErrSum;
	double logSignedSum;
	double linErrNum;
	double linSignedNum;
	double linErrDen;
	long transitionCount;

	// Per-tier error tracking
	double[] tierLogErr;
	double[] tierLogSigned;
	double[] tierLinNum;
	double[] tierLinSigned;
	double[] tierLinDen;
	long[] tierTransitions;

	// Diagnostics
	long belowFloorCount;
	long totalAdds;
	long advanceCount;

	public DLL4WordSimulator(int iters, int threads, int maxTier, int words, int mode) {
		this.iters   = iters;
		this.threads = threads;
		this.maxTier = maxTier;
		this.tiers   = maxTier + 2;
		this.words   = words;
		this.mode    = mode;
		this.counts  = new long[tiers][NUM_CANONICAL];
		this.sums    = new double[tiers][NUM_CANONICAL];
		this.sumSq   = new double[tiers][NUM_CANONICAL];
		this.tierLogErr     = new double[tiers];
		this.tierLogSigned  = new double[tiers];
		this.tierLinNum     = new double[tiers];
		this.tierLinSigned  = new double[tiers];
		this.tierLinDen     = new double[tiers];
		this.tierTransitions = new long[tiers];
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

	void loadTable(String fname) {
		ByteFile bf = ByteFile.makeByteFile(fname, false);
		int maxLoadedTier = 0;
		// First pass: find max tier
		for (byte[] line = bf.nextLine(); line != null; line = bf.nextLine()) {
			if (line.length < 1 || line[0] == '#') continue;
			String s = new String(line);
			String[] parts = s.split("\t");
			int tier = Integer.parseInt(parts[0]);
			if (tier > maxLoadedTier) maxLoadedTier = tier;
		}
		bf.close();

		int tableTiers = maxLoadedTier + 1;
		double[][] table = new double[tableTiers][NUM_CANONICAL];

		bf = ByteFile.makeByteFile(fname, false);
		int entries = 0;
		for (byte[] line = bf.nextLine(); line != null; line = bf.nextLine()) {
			if (line.length < 1 || line[0] == '#') continue;
			String s = new String(line);
			String[] parts = s.split("\t");
			int tier = Integer.parseInt(parts[0]);
			int idx = Integer.parseInt(parts[1]);
			double avg = Double.parseDouble(parts[3]);
			if (idx < NUM_CANONICAL) {
				table[tier][idx] = avg;
				entries++;
			}
		}
		bf.close();
		this.loadedTable = table;
		System.err.println("Loaded state table: " + entries + " entries, " + tableTiers + " tiers from " + fname);
	}

	double getWordEstimate(int[] reg, int globalExp) {
		if (reg[0] == 0 && reg[1] == 0 && reg[2] == 0 && reg[3] == 0 && globalExp == 0) return 0;
		int wMin = Math.min(Math.min(reg[0], reg[1]), Math.min(reg[2], reg[3]));
		int tier = globalExp + wMin;
		int raw = reg[0] | (reg[1] << 4) | (reg[2] << 8) | (reg[3] << 12);
		int idx = REMAP[raw];
		if (tier < loadedTable.length) return loadedTable[tier][idx];
		// Beyond table: use last tier's value (best we have)
		return loadedTable[loadedTable.length - 1][idx];
	}

	double computeFullEstimate(int[][] regs, int globalExp, double[] wordContrib) {
		double est = 0;
		for (int w = 0; w < words; w++) {
			wordContrib[w] = getWordEstimate(regs[w], globalExp);
			est += wordContrib[w];
		}
		return est;
	}

	void simulate() throws InterruptedException {
		long[][][]   tCounts  = new long[threads][tiers][NUM_CANONICAL];
		double[][][] tSums    = new double[threads][tiers][NUM_CANONICAL];
		double[][][] tSumSq   = new double[threads][tiers][NUM_CANONICAL];
		long[]   tBelow       = new long[threads];
		long[]   tTotal       = new long[threads];
		long[]   tAdvance     = new long[threads];
		double[] tLogErr      = new double[threads];
		double[] tLogSigned   = new double[threads];
		double[] tLinNum      = new double[threads];
		double[] tLinSigned   = new double[threads];
		double[] tLinDen      = new double[threads];
		long[]   tTransitions = new long[threads];
		double[][] tTierLogErr    = new double[threads][tiers];
		double[][] tTierLogSigned = new double[threads][tiers];
		double[][] tTierLinNum    = new double[threads][tiers];
		double[][] tTierLinSigned = new double[threads][tiers];
		double[][] tTierLinDen    = new double[threads][tiers];
		long[][]   tTierTrans     = new long[threads][tiers];

		int base = iters / threads;
		int rem  = iters % threads;

		Thread[] threadArr = new Thread[threads];
		for (int t = 0; t < threads; t++) {
			final int tid     = t;
			final int myIters = base + (tid < rem ? 1 : 0);

			threadArr[t] = new Thread(() -> {
				FastRandomXoshiro rng = new FastRandomXoshiro(tid + 1);
				long[][]   lCounts = tCounts[tid];
				double[][] lSums   = tSums[tid];
				double[][] lSumSq  = tSumSq[tid];
				long lBelow = 0, lTotal = 0, lAdvance = 0;
				double lLogErr = 0, lLogSigned = 0, lLinNum = 0, lLinSigned = 0, lLinDen = 0;
				long lTransitions = 0;
				double[] lTierLogErr    = tTierLogErr[tid];
				double[] lTierLogSigned = tTierLogSigned[tid];
				double[] lTierLinNum    = tTierLinNum[tid];
				double[] lTierLinSigned = tTierLinSigned[tid];
				double[] lTierLinDen    = tTierLinDen[tid];
				long[]   lTierTrans     = tTierTrans[tid];

				final int totalBuckets = BUCKETS_PER_WORD * words;
				final boolean trackErrors = (loadedTable != null);
				final boolean singleWord = (words == 1);

				for (int iter = 0; iter < myIters; iter++) {
					int[][] regs = new int[words][4];
					int globalExp = 0;
					int totalZeros = totalBuckets;
					long eeMask = -1L;

					double[] wordContrib = trackErrors ? new double[words] : null;
					double runningEst = 0;

					for (long trueCard = 1; trueCard <= 10_000_000L; trueCard++) {
						long hash = hash64shift(rng.nextLong());
						lTotal++;

						// --- Determine if hash causes a register update ---
						boolean skip = false;
						boolean changed = false;
						int wordIdx = 0;
						int oldTier = -1, oldIdx = -1;
						int newTier = -1, newIdx = -1;

						if (Long.compareUnsigned(hash, eeMask) > 0) {
							skip = true;
						} else {
							int hashNlz = Long.numberOfLeadingZeros(hash);
							int bucket = (int)((hash & 0x7FFF_FFFFL) % totalBuckets);
							wordIdx = bucket / BUCKETS_PER_WORD;
							int localBucket = bucket % BUCKETS_PER_WORD;
							int relNlz = hashNlz - globalExp;

							if (relNlz < 0) {
								skip = true;
							} else {
								int newStored = Math.min(relNlz + 1, 15);
								int oldStored = regs[wordIdx][localBucket];

								if (newStored > oldStored) {
									// Save old state of affected word
									int[] wr = regs[wordIdx];
									int oldRaw = wr[0] | (wr[1] << 4) | (wr[2] << 8) | (wr[3] << 12);
									oldIdx = REMAP[oldRaw];
									int oldWMin = Math.min(Math.min(wr[0], wr[1]), Math.min(wr[2], wr[3]));
									oldTier = globalExp + oldWMin;

									// Update register
									regs[wordIdx][localBucket] = newStored;

									// Handle globalExp advance
									int oldGlobalExp = globalExp;
									if (oldStored == 0) {
										totalZeros--;
										while (totalZeros == 0 && globalExp < 64) {
											globalExp++;
											eeMask >>>= 1;
											totalZeros = 0;
											for (int w = 0; w < words; w++) {
												for (int b = 0; b < 4; b++) {
													regs[w][b] = Math.max(0, regs[w][b] - 1);
													if (regs[w][b] == 0) totalZeros++;
												}
											}
											lAdvance++;
										}
									}

									// Compute new state
									wr = regs[wordIdx];
									int newRaw = wr[0] | (wr[1] << 4) | (wr[2] << 8) | (wr[3] << 12);
									newIdx = REMAP[newRaw];
									int newWMin = Math.min(Math.min(wr[0], wr[1]), Math.min(wr[2], wr[3]));
									newTier = globalExp + newWMin;

									changed = (oldTier != newTier || oldIdx != newIdx);

									// Error tracking: update running estimate
									if (trackErrors && changed) {
										if (globalExp != oldGlobalExp) {
											runningEst = computeFullEstimate(regs, globalExp, wordContrib);
										} else {
											runningEst -= wordContrib[wordIdx];
											wordContrib[wordIdx] = getWordEstimate(wr, globalExp);
											runningEst += wordContrib[wordIdx];
										}
									}
								}
								// else: newStored <= oldStored, no change
							}
						}
						if (skip) lBelow++;

						// --- Recording ---
						// Per-word cardinality: each word sees ~1/words of total hashes
						double recordCard = (double) trueCard / words;

						if (mode == MODE_EVERY) {
							if (singleWord) {
								int[] wr = regs[0];
								int wm = Math.min(Math.min(wr[0], wr[1]), Math.min(wr[2], wr[3]));
								int tier = globalExp + wm;
								int raw = wr[0] | (wr[1] << 4) | (wr[2] << 8) | (wr[3] << 12);
								int idx = REMAP[raw];
								if (tier < tiers) {
									lCounts[tier][idx]++;
									lSums[tier][idx] += recordCard;
									lSumSq[tier][idx] += recordCard * recordCard;
								}
							}
						} else if (changed) {
							// Entry: new state at recordCard
							if (newTier < tiers) {
								lCounts[newTier][newIdx]++;
								lSums[newTier][newIdx] += recordCard;
								lSumSq[newTier][newIdx] += recordCard * recordCard;
							}
							if (mode == MODE_BOTH) {
								// Exit: old state at (trueCard-1)/words
								double exitCard = (double)(trueCard - 1) / words;
								if (exitCard > 0 && oldTier < tiers) {
									lCounts[oldTier][oldIdx]++;
									lSums[oldTier][oldIdx] += exitCard;
									lSumSq[oldTier][oldIdx] += exitCard * exitCard;
								}
							}
						}

						// --- Error tracking ---
						if (trackErrors && changed) {
							if (trueCard > 10 && runningEst > 0) {
								double diff = runningEst - trueCard;
								double relErr = Math.abs(diff) / trueCard;
								double signedRel = diff / trueCard;
								lLogErr += relErr;
								lLogSigned += signedRel;
								lLinNum += Math.abs(diff);
								lLinSigned += diff;
								lLinDen += trueCard;
								lTransitions++;
								// Per-tier: bin by tier of the transitioning word
								if (newTier >= 0 && newTier < tiers) {
									lTierLogErr[newTier] += relErr;
									lTierLogSigned[newTier] += signedRel;
									lTierLinNum[newTier] += Math.abs(diff);
									lTierLinSigned[newTier] += diff;
									lTierLinDen[newTier] += trueCard;
									lTierTrans[newTier]++;
								}
							}
						}

						// Termination: for single-word, check tier; for multi-word, rely on 10M limit
						if (singleWord && newTier > maxTier) break;
					}
				}

				tBelow[tid]       = lBelow;
				tTotal[tid]       = lTotal;
				tAdvance[tid]     = lAdvance;
				tLogErr[tid]      = lLogErr;
				tLogSigned[tid]   = lLogSigned;
				tLinNum[tid]      = lLinNum;
				tLinSigned[tid]   = lLinSigned;
				tLinDen[tid]      = lLinDen;
				tTransitions[tid] = lTransitions;
			});
			threadArr[t].start();
		}

		for (Thread th : threadArr) th.join();

		for (int t = 0; t < threads; t++) {
			for (int tier = 0; tier < tiers; tier++) {
				for (int s = 0; s < NUM_CANONICAL; s++) {
					counts[tier][s] += tCounts[t][tier][s];
					sums[tier][s]   += tSums[t][tier][s];
					sumSq[tier][s]  += tSumSq[t][tier][s];
				}
			}
			belowFloorCount += tBelow[t];
			totalAdds       += tTotal[t];
			advanceCount    += tAdvance[t];
			logErrSum       += tLogErr[t];
			logSignedSum    += tLogSigned[t];
			linErrNum       += tLinNum[t];
			linSignedNum    += tLinSigned[t];
			linErrDen       += tLinDen[t];
			transitionCount += tTransitions[t];
			for (int tier = 0; tier < tiers; tier++) {
				tierLogErr[tier]     += tTierLogErr[t][tier];
				tierLogSigned[tier]  += tTierLogSigned[t][tier];
				tierLinNum[tier]     += tTierLinNum[t][tier];
				tierLinSigned[tier]  += tTierLinSigned[t][tier];
				tierLinDen[tier]     += tTierLinDen[t][tier];
				tierTransitions[tier] += tTierTrans[t][tier];
			}
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

	/** Weighted average CV per tier (weighted by observation count). */
	double[] buildTierCV(double[][] cv) {
		double[] tierCV = new double[tiers];
		for (int tier = 0; tier < tiers; tier++) {
			double wSum = 0, wTot = 0;
			for (int s = 0; s < NUM_CANONICAL; s++) {
				if (counts[tier][s] > 0) {
					wSum += cv[tier][s] * counts[tier][s];
					wTot += counts[tier][s];
				}
			}
			tierCV[tier] = wTot > 0 ? wSum / wTot : 0;
		}
		return tierCV;
	}

	void printResults(double[][] stateAvg, double[][] cv) {
		System.err.println("=== DLL4 Word Simulator Results ===");
		System.err.println("iters=" + iters + "  threads=" + threads + "  maxTier=" + maxTier
			+ "  words=" + words + "  mode=" + modeName(mode));
		System.err.println("numCanonical=" + NUM_CANONICAL);
		System.err.println("totalAdds=" + totalAdds);
		System.err.printf("belowFloorRate=%.6f%n",
			totalAdds > 0 ? (double) belowFloorCount / totalAdds : 0);
		System.err.printf("advanceCount  =%d%n", advanceCount);

		if (loadedTable != null && transitionCount > 0) {
			System.err.println();
			System.err.println("=== Error Tracking ===");
			System.err.printf("transitions    =%d%n", transitionCount);
			System.err.printf("logAvgErr      =%.4f%% (mean |est-true|/true)%n", 100.0 * logErrSum / transitionCount);
			System.err.printf("logSignedErr   =%+.4f%% (mean (est-true)/true)%n", 100.0 * logSignedSum / transitionCount);
			System.err.printf("linAvgErr      =%.4f%% (sum|est-true| / sumTrue)%n", 100.0 * linErrNum / linErrDen);
			System.err.printf("linSignedErr   =%+.4f%% (sum(est-true) / sumTrue)%n", 100.0 * linSignedNum / linErrDen);
		}

		double[] tierCV = buildTierCV(cv);
		boolean hasErrors = (loadedTable != null && transitionCount > 0);
		System.err.println();
		if (hasErrors) {
			System.err.printf("%-6s %-14s %-10s %-10s %-6s %-8s %-10s %-10s%n",
				"Tier", "AvgCard", "Growth", "Obs", "St", "CV", "LogErr%", "Signed%");
		} else {
			System.err.printf("%-6s %-14s %-10s %-10s %-6s %-8s%n",
				"Tier", "AvgCard", "Growth", "Obs", "St", "CV");
		}
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
			if (hasErrors && tierTransitions[tier] > 0) {
				double tLogPct = 100.0 * tierLogErr[tier] / tierTransitions[tier];
				double tSignedPct = 100.0 * tierLogSigned[tier] / tierTransitions[tier];
				System.err.printf("%-6d %-14.2f %-10.4f %-10d %-6d %-8.4f %-10.2f %+-10.2f%n",
					tier, avg, growth, obs, statesUsed, tierCV[tier], tLogPct, tSignedPct);
			} else {
				System.err.printf("%-6d %-14.2f %-10.4f %-10d %-6d %-8.4f%n",
					tier, avg, growth, obs, statesUsed, tierCV[tier]);
			}
			prevAvg = avg;
		}

		// Sparse TSV output
		System.err.println();
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

	static String modeName(int m) {
		switch (m) {
			case MODE_EVERY: return "every";
			case MODE_ENTRY: return "entry";
			case MODE_BOTH:  return "both";
			default: return "unknown";
		}
	}

	static int parseMode(String s) {
		switch (s.toLowerCase()) {
			case "every": case "all": return MODE_EVERY;
			case "entry": case "in":  return MODE_ENTRY;
			case "both":  case "io":  return MODE_BOTH;
			default:
				System.err.println("Unknown mode: " + s + ", using 'every'");
				return MODE_EVERY;
		}
	}

	public static void main(String[] args) throws InterruptedException {
		int iters   = 10000;
		int threads = 8;
		int maxTier = 12;
		int words   = 1;
		int mode    = MODE_EVERY;
		String tableFile = null;

		for (String arg : args) {
			String[] kv = arg.split("=", 2);
			if (kv.length != 2) continue;
			switch (kv[0]) {
				case "iters":   iters     = Integer.parseInt(kv[1]);  break;
				case "threads": threads   = Integer.parseInt(kv[1]);  break;
				case "maxTier": maxTier   = Integer.parseInt(kv[1]);  break;
				case "words":   words     = Integer.parseInt(kv[1]);  break;
				case "mode":    mode      = parseMode(kv[1]);         break;
				case "table":   tableFile = kv[1];                    break;
				default: System.err.println("Unknown arg: " + arg);
			}
		}

		System.err.println("DLL4WordSimulator: iters=" + iters + " threads=" + threads
			+ " maxTier=" + maxTier + " words=" + words + " mode=" + modeName(mode)
			+ " numCanonical=" + NUM_CANONICAL
			+ (tableFile != null ? " table=" + tableFile : ""));

		long t0 = System.currentTimeMillis();
		DLL4WordSimulator sim = new DLL4WordSimulator(iters, threads, maxTier, words, mode);
		if (tableFile != null) {
			sim.loadTable(tableFile);
		}
		sim.simulate();
		long t1 = System.currentTimeMillis();
		System.err.println("Simulation time: " + (t1 - t0) + " ms");

		double[][] stateAvg = sim.buildStateAvg();
		double[][] cv = sim.buildCV();
		sim.printResults(stateAvg, cv);
	}
}
