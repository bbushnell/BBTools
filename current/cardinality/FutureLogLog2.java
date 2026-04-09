package cardinality;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
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
	/** Number of equivalence classes for 6 buckets with 4 states = C(9,3). */
	static final int NUM_EQUIV = 84;
	/** Prefix sums for idx84 computation. */
	static final int[] OFFSET_A = {0, 28, 49, 64, 74, 80, 83};
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

		if (CLAMP_OVERFLOW) {
			// Overflow hashes (delta>1) set the MSB: proves NLZ >= absNLZ+1
			if (delta < 0) { return; }
		} else {
			// IOT: only delta 0 (floor hit) and 1 (future hit) are representable
			if (delta < 0 || delta > 1) { return; }
		}
		final int clampedDelta = Math.min(delta, 1);

		final int bitToSet = 1 << (clampedDelta + register * 2);
		if ((word & bitToSet) != 0) { return; } // bit already set

		lastCardinality = -1;
		lastEstimates = null;
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
	/*----------------    Equivalence Class Index   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts a 12-bit history state to its 84-slot equivalence class index.
	 * Bucket histogram (a,b,c,d) where a=count(00), b=count(01), c=count(10),
	 * d=count(11), a+b+c+d=6. Zero allocation, ~15 ops, 3 popcounts.
	 */
	static int idx84(int history) {
		int h1 = history >>> 1;
		int d = Integer.bitCount(history & h1 & LSB_MASK);
		int a = 6 - Integer.bitCount((history | h1) & LSB_MASK);
		int c = Integer.bitCount(~history & h1 & LSB_MASK);
		int b = 6 - a - c - d;
		int r = 6 - a;
		return OFFSET_A[a] + b * (r + 1) - b * (b - 1) / 2 + c;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Estimation           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality() {
		if (lastCardinality >= 0) { return lastCardinality; }
		double fllEst = rawCFEstimate();
		double microEst = microCardinalityDouble();
		long card = (long) blendEstimate(fllEst, microEst);
		card = Math.max(card, 0);
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

	/** Fraction of words with any bucket activity (localExp > 0 or history != 0). */
	public double occupancy() {
		int filled = 0;
		for (int i = 0; i < numWords; i++) {
			if ((words[i] & 0xFFFF) != 0) { filled++; }
		}
		return (double) filled / numWords;
	}

	/** Flat correction multiplier for state-table bias. Modifiable at runtime via fll2mult=. */
	public static double TERMINAL_CORRECTION = 0.883319848;
	/** When true, overflow hashes (delta>1) set the MSB (future bit) instead of being ignored. */
	public static boolean CLAMP_OVERFLOW = false;

	/** Raw CF-table estimate without clamp or microCardinality. */
	private double rawCFEstimate() {
		if (CF_TABLE_TIERS <= 0) { return 0; }
		double sum = 0;
		for (int i = 0; i < numWords; i++) {
			final int w = words[i] & 0xFFFF;
			if (w == 0) { continue; }
			final int localExp = (w >>> 12) & 0xF;
			final int history = w & HISTORY_MASK;
			final int absNLZ = localExp + globalExp;
			final int tier = Math.min(absNLZ, CF_TABLE_TIERS - 1);
			sum += tierCardinality(absNLZ) * cfTable[tier][idx84(history)];
		}
		sum *= TERMINAL_CORRECTION;
		// Apply cardinality-keyed CF iteratively if available
		if (cardCFKeys != null && cardCFKeys.length > 0) {
			double keyScale = (cardCFBuckets > 0 && modBuckets != cardCFBuckets)
				? (double) cardCFBuckets / modBuckets : 1.0;
			double cf = CorrectionFactor.getCF(sum, sum,
				cardCFValues, cardCFKeys, 5, 0.0001, keyScale);
			sum *= cf;
		}
		return sum;
	}

	@Override
	public double[] rawEstimates() {
		if (lastEstimates != null && lastCardinality >= 0) { return lastEstimates; }
		final double fllEst = rawCFEstimate();
		final double microEst = microCardinalityDouble();
		final double hybridEst = blendEstimate(fllEst, microEst);
		final int total = 17 + AbstractCardStats.NUM_DLC_TIERS + AbstractCardStats.NUM_EXTRA;
		final double[] r = new double[total];
		java.util.Arrays.fill(r, hybridEst);
		r[0] = fllEst;    // Mean = FLLPure (CF table always generated from this)
		r[5] = microEst;   // LC = microIndex estimate
		r[6] = hybridEst;  // Hybrid = blended FLL+micro
		lastEstimates = r;
		return r;
	}

	/** Floating-point microCardinality (LC on 64-bit microIndex),
	 *  floored by the total set bits across all word histories. */
	private double microCardinalityDouble() {
		final int bits = Long.bitCount(microIndex);
		double micro = (bits >= 64) ? 64.0 : -64.0 * Math.log(1.0 - bits / 64.0);
		// Floor: total set history bits across all words.
		// Each element sets at most 1 bit, so this is a lower bound.
		int totalSetBits = 0;
		for (int i = 0; i < numWords; i++) {
			totalSetBits += Integer.bitCount(words[i] & HISTORY_MASK);
		}
		return Math.max(micro, totalSetBits);
	}

	private static final double BLEND_LO = 8.0;
	private static final double BLEND_HI = 100.0;
	private static final double BLEND_LOG_LO = Math.log(BLEND_LO);
	private static final double BLEND_LOG_RANGE = Math.log(BLEND_HI) - Math.log(BLEND_LO);

	/**
	 * Blends micro and FLL estimates on an exponential scale.
	 * 100% micro below 8, 100% FLL above 80, log-linear blend between.
	 * Midpoint (50/50) at sqrt(8×80) ≈ 25.3.
	 *
	 * Bootstrap: uses microEst as the blend key when the microIndex isn't
	 * saturated (reliable below ~50), switches to fllEst when it is.
	 */
	private double blendEstimate(double fllEst, double microEst) {
		final int microBits = Long.bitCount(microIndex);
		// Micro is reliable when not near saturation (64 bits)
		double key = (microBits < 58) ? microEst : fllEst;
		key = Math.max(key, 1.0);
		if (key <= BLEND_LO) { return microEst; }
		if (key >= BLEND_HI) { return fllEst; }
		double frac = (Math.log(key) - BLEND_LOG_LO) / BLEND_LOG_RANGE;
		return microEst * (1.0 - frac) + fllEst * frac;
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
	/** Debug: return raw word value at index i. */
	public int getWord(int i) { return words[i] & 0xFFFF; }

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
	/** Cached rawEstimates array; invalidated when lastCardinality is set to -1. */
	private double[] lastEstimates;

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-tier correction multipliers: cfTable[tier][idx84]. */
	static double[][] cfTable;
	/** Number of tiers in the loaded table. */
	static int CF_TABLE_TIERS = 0;
	/** Average cardinality at each tier (from simulator). */
	static double[] tierAvg;

	/** Cardinality-keyed CF values (from calibration). */
	static float[] cardCFValues;
	/** Cardinality keys for binary search into cardCFValues. */
	static float[] cardCFKeys;
	/** Bucket count the cardinality CF table was built at. */
	static int cardCFBuckets = 0;

	/** Resource file for per-state correction factors (84-slot equivalence classes). */
	static final String STATE_TABLE_FILE = "stateTableFLL2.tsv.gz";
	/** Resource file for cardinality-keyed correction factors (from calibration). */
	static final String CF_FILE = "cardinalityCorrectionFLL2.tsv.gz";

	/** Static initializer: auto-load state table and cardinality CF on first class access.
	 *  MUST appear after all static field declarations to avoid initializer-order bugs. */
	static {
		loadCFTable();
		loadCardCFTable();
	}

	/**
	 * Loads the FLL2 correction factor table from resources.
	 * Format: tierAvg lines, then header, then per-tier count/mult rows.
	 * Each mult row has 84 columns (one per equivalence class).
	 * Called once at startup (e.g., by DDLCalibrationDriver2 or on first use).
	 */
	public static synchronized void loadCFTable() {
		if (cfTable != null) { return; }
		String path = Data.findPath("?" + STATE_TABLE_FILE);
		if (path == null) {
			System.err.println("WARNING: FLL2 CF table not found: " + CF_FILE);
			setFallbackCF();
			return;
		}
		FileFormat ff = FileFormat.testInput(path, null, false);
		if (ff == null) { setFallbackCF(); return; }
		ByteFile bf = ByteFile.makeByteFile(ff, 1);
		if (bf == null) { setFallbackCF(); return; }

		java.util.ArrayList<double[]> multRows = new java.util.ArrayList<>();
		java.util.ArrayList<Double> avgList = new java.util.ArrayList<>();

		for (byte[] line = bf.nextLine(); line != null; line = bf.nextLine()) {
			if (line.length == 0) { continue; }
			String s = new String(line).trim();
			if (s.startsWith("#") || s.startsWith("Tier\t")) { continue; }
			String[] parts = s.split("\t");
			if (parts.length < 3) { continue; }

			if (parts[0].equals("tierAvg")) {
				int tier = Integer.parseInt(parts[1]);
				double avg = Double.parseDouble(parts[2]);
				while (avgList.size() <= tier) { avgList.add(0.0); }
				avgList.set(tier, avg);
			} else {
				int tier;
				try { tier = Integer.parseInt(parts[0]); }
				catch (NumberFormatException e) { continue; }
				if (!"mult".equals(parts[1])) { continue; }

				double[] row = new double[NUM_EQUIV];
				java.util.Arrays.fill(row, 1.0);
				for (int i = 2; i < parts.length && (i - 2) < NUM_EQUIV; i++) {
					try { row[i - 2] = Double.parseDouble(parts[i]); }
					catch (NumberFormatException e) { /* leave as 1.0 */ }
				}
				while (multRows.size() <= tier) { multRows.add(null); }
				multRows.set(tier, row);
			}
		}
		bf.close();

		tierAvg = new double[avgList.size()];
		for (int i = 0; i < avgList.size(); i++) { tierAvg[i] = avgList.get(i); }

		CF_TABLE_TIERS = multRows.size();
		cfTable = new double[CF_TABLE_TIERS][NUM_EQUIV];
		for (int t = 0; t < CF_TABLE_TIERS; t++) {
			if (multRows.get(t) != null) {
				cfTable[t] = multRows.get(t);
			} else {
				java.util.Arrays.fill(cfTable[t], 1.0);
			}
		}

		System.err.println("Loaded FLL2 CF table: " + CF_TABLE_TIERS + " tiers, "
			+ tierAvg.length + " tier averages, " + NUM_EQUIV + " equiv classes");
	}

	/**
	 * Loads the cardinality-keyed CF table (V5 format from calibration).
	 * Columns: #TrueCard Mean_cf HMean_cf ... — we use Mean_cf (column 1).
	 */
	public static synchronized void loadCardCFTable() {
		if (cardCFKeys != null) { return; }
		String path = Data.findPath("?" + CF_FILE);
		if (path == null) {
			System.err.println("WARNING: FLL2 cardinality CF table not found: " + CF_FILE);
			return;
		}
		FileFormat ff = FileFormat.testInput(path, null, false);
		if (ff == null) { return; }
		ByteFile bf = ByteFile.makeByteFile(ff, 1);
		if (bf == null) { return; }

		structures.FloatList keys = new structures.FloatList();
		structures.FloatList vals = new structures.FloatList();

		for (byte[] line = bf.nextLine(); line != null; line = bf.nextLine()) {
			if (line.length == 0) { continue; }
			String s = new String(line).trim();
			if (s.startsWith("#Buckets")) {
				String[] parts = s.split("\t");
				if (parts.length >= 2) { cardCFBuckets = Integer.parseInt(parts[1].trim()); }
				continue;
			}
			if (s.startsWith("#") || !Character.isDigit(s.charAt(0))) { continue; }
			String[] parts = s.split("\t");
			if (parts.length < 2) { continue; }
			float key = Float.parseFloat(parts[0]);
			float cf = Float.parseFloat(parts[1]); // Mean_cf column
			keys.add(key);
			vals.add(cf);
		}
		bf.close();

		cardCFKeys = keys.toArray();
		cardCFValues = vals.toArray();
		System.err.println("Loaded FLL2 cardinality CF table: " + cardCFKeys.length
			+ " entries, buckets=" + cardCFBuckets);
	}

	private static void setFallbackCF() {
		cfTable = new double[1][NUM_EQUIV];
		java.util.Arrays.fill(cfTable[0], 1.0);
		tierAvg = new double[]{1.0};
		CF_TABLE_TIERS = 1;
	}
}
