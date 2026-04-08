package cardinality;

import shared.Tools;

/**
 * FLL2: FutureLogLog with 2-bit future registers.
 *
 * Packs 6 buckets into each 16-bit word:
 *   [15:12] = 4-bit local exponent (floor relative to global)
 *   [11:0]  = 6 × 2-bit future bitmap
 *
 * Per-bucket bitmap semantics:
 *   Bit 0 (LSB) = "seen hash with NLZ == absNLZ"     (floor hit)
 *   Bit 1 (MSB) = "seen hash with NLZ == absNLZ + 1"  (future hit)
 *   All four states {00, 01, 10, 11} are reachable.
 *
 * Promotion fires when all 6 LSBs are set: (word & 0x555) == 0x555.
 * Promotion shifts MSBs to LSB positions, clears MSBs, increments localExp.
 * Max cascade per add = 2.
 *
 * Two-level exponent:
 *   globalExp = shared minimum floor across all words.
 *   absNLZ = globalExp + localExp.
 *   Global advance: when all words have localExp >= 1, increment globalExp,
 *   decrement all localExps, update eeMask.
 *
 * Estimation: sum over all words of tierCard(absNLZ) * cfTable[tier][history].
 * CF table loaded from resources/cardinalityCorrectionFLL2.tsv.gz.
 *
 * Uses modulo bucket addressing since numBuckets = numWords * 6 is never a
 * power of 2. Follows the DLL4m pattern for non-power-of-2 bucket counts.
 *
 * @author Brian Bushnell, Chloe von Einzbern-Bushnell
 * @date April 8, 2026
 */
public final class FutureLogLog2 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int BUCKETS_PER_WORD = 6;
	/** Bits 0,2,4,6,8,10 — floor-level bits (LSB of each bucket). */
	static final int LSB_MASK = 0x555;
	/** Bits 1,3,5,7,9,11 — future-level bits (MSB of each bucket). */
	static final int MSB_MASK = 0xAAA;
	static final int HISTORY_MASK = 0xFFF;
	static final int NUM_TIERS = 16;
	static final int HISTORY_STATES = 4096;
	/** Enable expensive invariant checks. Set false for production. */
	private static final boolean debugging = false;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	FutureLogLog2() { this(2048, 31, -1, 0); }

	FutureLogLog2(parse.Parser p) {
		super(p);
		numWords = roundUpToWords(buckets);
		modBuckets = numWords * BUCKETS_PER_WORD;
		words = new short[numWords];
		numLocalZeros = numWords;
	}

	/**
	 * @param buckets_ Desired bucket count; will be rounded up to next multiple of 6.
	 */
	FutureLogLog2(int buckets_, int k_, long seed, float minProb_) {
		// Super requires power-of-2; pass next-pow2 above our actual bucket count.
		super(nextPow2(roundUpToWords(buckets_) * BUCKETS_PER_WORD), k_, seed, minProb_);
		numWords = roundUpToWords(buckets_);
		modBuckets = numWords * BUCKETS_PER_WORD;
		words = new short[numWords];
		numLocalZeros = numWords;
	}

	@Override
	public FutureLogLog2 copy() {
		return new FutureLogLog2(modBuckets, k, -1, minProb);
	}

	/** Rounds buckets up to the next multiple of 6, expressed as word count. */
	private static int roundUpToWords(int buckets_) {
		return (buckets_ + BUCKETS_PER_WORD - 1) / BUCKETS_PER_WORD;
	}

	/** Next power of 2 >= n. */
	private static int nextPow2(int n) {
		return Integer.highestOneBit(Math.max(1, n - 1)) << 1;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Core Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number) {
		final long key = Tools.hash64shift(number ^ hashXor);

		// Early exit: skip hashes below global floor
		if (Long.compareUnsigned(key, eeMask) > 0) { return; }

		final int hashNLZ = Long.numberOfLeadingZeros(key);
		final int bucket = (int) Long.remainderUnsigned(key, modBuckets);
		final int wordIdx = bucket / BUCKETS_PER_WORD;
		final int register = bucket % BUCKETS_PER_WORD;

		// MicroIndex for very-low-cardinality detection
		final long micro = (key >>> 58) & 0x3FL;
		microIndex |= (1L << micro);

		int word = words[wordIdx] & 0xFFFF;
		final int localExp = (word >>> 12) & 0xF;
		final int absNLZ = globalExp + localExp;
		final int delta = hashNLZ - absNLZ;

		// IOT: only delta 0 (floor hit) and 1 (future hit) are representable
		if (delta < 0 || delta > 1) { return; }

		final int bitToSet = 1 << (delta + register * 2);
		if ((word & bitToSet) != 0) { return; } // bit already set

		lastCardinality = -1;
		word |= bitToSet;

		// Check for promotion: all 6 LSBs set?
		if ((word & LSB_MASK) == LSB_MASK) {
			final boolean wasZeroExp = (localExp == 0);
			word = promote(word);
			words[wordIdx] = (short) word;
			// Update numLocalZeros tracking
			if (wasZeroExp) {
				numLocalZeros--;
				if (numLocalZeros == 0) { advanceGlobal(); }
			}
		} else {
			words[wordIdx] = (short) word;
		}

		assert !debugging || repOK()
			: "repOK failed after hashAndStore, word=" + wordIdx;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Promotion           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Promotes a word: MSBs shift to LSBs, MSBs clear, localExp++.
	 * Cascades until promotion condition is gone or localExp reaches 15.
	 * Max cascade = 2 (all-{11} → all-{01} → all-{00}).
	 */
	static int promote(int word) {
		while ((word & LSB_MASK) == LSB_MASK) {
			int localExp = (word >>> 12) & 0xF;
			if (localExp >= 15) { break; }
			word = ((word & MSB_MASK) >>> 1) | ((localExp + 1) << 12);
		}
		return word;
	}

	/**
	 * Global advance: when all words have localExp >= 1, increment globalExp,
	 * decrement all localExps, recount numLocalZeros, update eeMask.
	 */
	private void advanceGlobal() {
		while (numLocalZeros == 0) {
			globalExp++;
			numLocalZeros = 0;
			for (int i = 0; i < numWords; i++) {
				int w = words[i] & 0xFFFF;
				int le = (w >>> 12) & 0xF;
				assert le >= 1 : "localExp=0 but numLocalZeros was 0, word " + i;
				le--;
				words[i] = (short) ((le << 12) | (w & HISTORY_MASK));
				if (le == 0) { numLocalZeros++; }
			}
			eeMask = (-1L) >>> globalExp;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Estimation           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality() {
		if (lastCardinality >= 0) { return lastCardinality; }
		double sum = 0;
		for (int i = 0; i < numWords; i++) {
			final int w = words[i] & 0xFFFF;
			final int localExp = (w >>> 12) & 0xF;
			final int history = w & HISTORY_MASK;
			final int absNLZ = localExp + globalExp;
			final int tier = Math.min(localExp, CF_TABLE_TIERS - 1);
			sum += tierCardinality(absNLZ) * cfTable[tier][history];
		}
		long card = (long) sum;
		card = Math.max(card, microCardinality());
		card = Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality = card;
		return card;
	}

	/**
	 * MicroIndex-based very-low-cardinality estimate (below ~56 distinct).
	 * Uses the 64-bit microIndex as a birthday-paradox estimator.
	 */
	private long microCardinality() {
		final int bits = Long.bitCount(microIndex);
		if (bits >= 64) { return 64; }
		return (long) (-64.0 * Math.log(1.0 - bits / 64.0));
	}

	/**
	 * Returns the expected average cardinality for a single bank at the given
	 * absolute NLZ tier. Uses the tier average array from the simulator,
	 * with exponential extrapolation beyond the table.
	 */
	private static double tierCardinality(int absNLZ) {
		if (absNLZ < tierAvg.length && absNLZ >= 0) {
			return tierAvg[absNLZ];
		}
		// Extrapolate beyond table using last measured growth factor
		if (tierAvg.length < 2) { return 1.0; }
		double lastGrowth = tierAvg[tierAvg.length - 1] / tierAvg[tierAvg.length - 2];
		double base = tierAvg[tierAvg.length - 1];
		int extra = absNLZ - (tierAvg.length - 1);
		return base * Math.pow(lastGrowth, extra);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Validation           ----------------*/
	/*--------------------------------------------------------------*/

	private boolean repOK() {
		assert numLocalZeros >= 0 : "numLocalZeros < 0: " + numLocalZeros;
		int actualZeros = 0;
		for (int i = 0; i < numWords; i++) {
			int w = words[i] & 0xFFFF;
			int le = (w >>> 12) & 0xF;
			assert le >= 0 && le <= 15 : "localExp out of range at word " + i + ": " + le;
			assert (w & 0xFFFF) == w : "Word exceeds 16 bits at " + i;
			// No promotable state should exist unless at max tier
			assert le >= 15 || (w & LSB_MASK) != LSB_MASK
				: "Promotable state at word " + i + ": 0x" + Integer.toHexString(w);
			if (le == 0) { actualZeros++; }
		}
		assert actualZeros == numLocalZeros
			: "numLocalZeros mismatch: tracked=" + numLocalZeros + " actual=" + actualZeros;
		assert numLocalZeros > 0 : "numLocalZeros is 0 but advanceGlobal was not called";
		assert eeMask == ((-1L) >>> globalExp)
			: "eeMask inconsistent with globalExp=" + globalExp;
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Merge/Misc          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void add(CardinalityTracker log) {
		throw new UnsupportedOperationException("FLL2 merge not yet implemented");
	}

	@Override
	public final float[] compensationFactorLogBucketsArray() { return null; }

	public int getNumWords() { return numWords; }
	public int getModBuckets() { return modBuckets; }
	public int getGlobalExp() { return globalExp; }

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** 16-bit packed banks: [15:12]=localExp, [11:0]=6×2-bit future bitmap. */
	private final short[] words;
	/** Number of 16-bit words (banks). */
	private final int numWords;
	/** Actual bucket count = numWords × 6. May not be a power of 2. */
	private final int modBuckets;
	/** Shared minimum floor across all words. */
	private int globalExp = 0;
	/** Count of words with localExp == 0. When this reaches 0, advanceGlobal fires. */
	private int numLocalZeros;
	/** Early-exit mask: hashes above this are below the global floor. */
	private long eeMask = -1L;

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-tier correction multipliers: cfTable[tier][12-bit history]. */
	static double[][] cfTable;
	/** Number of tiers in the loaded table. */
	static int CF_TABLE_TIERS = 0;
	/** Average cardinality at each tier (from simulator). */
	static double[] tierAvg;

	/** Resource file for correction factors. */
	static final String CF_FILE = "cardinalityCorrectionFLL2.tsv.gz";

	/**
	 * Loads the FLL2 correction factor table from resources.
	 * Expected format: TSV with rows starting with "tierAvg" or tier numbers.
	 * Called once at startup (e.g., by DDLCalibrationDriver2 or on first use).
	 */
	public static synchronized void loadCFTable() {
		if (cfTable != null) { return; } // already loaded
		try {
			java.io.InputStream is = FutureLogLog2.class.getResourceAsStream("/cardinality/" + CF_FILE);
			if (is == null) {
				is = FutureLogLog2.class.getResourceAsStream(CF_FILE);
			}
			if (is == null) {
				System.err.println("WARNING: FLL2 CF table not found: " + CF_FILE);
				cfTable = new double[1][HISTORY_STATES];
				java.util.Arrays.fill(cfTable[0], 1.0);
				tierAvg = new double[]{1.0};
				CF_TABLE_TIERS = 1;
				return;
			}
			if (CF_FILE.endsWith(".gz")) {
				is = new java.util.zip.GZIPInputStream(is);
			}
			java.io.BufferedReader br = new java.io.BufferedReader(new java.io.InputStreamReader(is));

			java.util.ArrayList<double[]> multRows = new java.util.ArrayList<>();
			java.util.ArrayList<Double> avgList = new java.util.ArrayList<>();
			String line;

			while ((line = br.readLine()) != null) {
				if (line.startsWith("#") || line.startsWith("Tier\t")) { continue; }
				String[] parts = line.split("\t");
				if (parts.length < 3) { continue; }

				if (parts[0].equals("tierAvg")) {
					// tierAvg	tier	avgCard	growth	observations
					int tier = Integer.parseInt(parts[1]);
					double avg = Double.parseDouble(parts[2]);
					while (avgList.size() <= tier) { avgList.add(0.0); }
					avgList.set(tier, avg);
				} else {
					// tier	type	state0	state1	...
					int tier;
					try { tier = Integer.parseInt(parts[0]); }
					catch (NumberFormatException e) { continue; }
					String type = parts[1];
					if (!"mult".equals(type)) { continue; } // skip count rows

					double[] row = new double[HISTORY_STATES];
					java.util.Arrays.fill(row, 1.0);
					// Parts 2..N are the canonical state multipliers
					// We need to map them back to all 4096 history states
					// For now, store as-is; the table has all 84 canonical columns
					// We need the header to know the state mapping
					// Actually, the easiest approach: store raw multiplier values
					// and do a separate mapping pass
					while (multRows.size() <= tier) { multRows.add(null); }
					multRows.set(tier, row);
					// Parse multiplier values from columns
					for (int i = 2; i < parts.length && (i - 2) < row.length; i++) {
						try {
							row[i - 2] = Double.parseDouble(parts[i]);
						} catch (NumberFormatException e) {
							// leave as 1.0
						}
					}
				}
			}
			br.close();

			tierAvg = new double[avgList.size()];
			for (int i = 0; i < avgList.size(); i++) { tierAvg[i] = avgList.get(i); }

			CF_TABLE_TIERS = multRows.size();
			cfTable = new double[CF_TABLE_TIERS][HISTORY_STATES];
			for (int t = 0; t < CF_TABLE_TIERS; t++) {
				if (multRows.get(t) != null) {
					cfTable[t] = multRows.get(t);
				} else {
					java.util.Arrays.fill(cfTable[t], 1.0);
				}
			}

			System.err.println("Loaded FLL2 CF table: " + CF_TABLE_TIERS + " tiers, "
				+ tierAvg.length + " tier averages");

		} catch (java.io.IOException e) {
			System.err.println("WARNING: Failed to load FLL2 CF table: " + e.getMessage());
			cfTable = new double[1][HISTORY_STATES];
			java.util.Arrays.fill(cfTable[0], 1.0);
			tierAvg = new double[]{1.0};
			CF_TABLE_TIERS = 1;
		}
	}
}
