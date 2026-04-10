package cardinality;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;

/**
 * BCLL: BankedCeilingLogLog — ceiling-based 2-bit register cardinality estimator.
 *
 * Same physical structure as FLL2 (16-bit words, 4-bit localExp, 6x2-bit
 * registers), but localExp represents a CEILING (max NLZ in word) instead
 * of a floor (min NLZ).  This makes updates idempotent: duplicates never
 * change state, because the ceiling only rises and old hashes can only
 * fall further below it.
 *
 * Word layout:
 *   [15:12] = 4-bit localExp  (offset from globalCeiling)
 *   [11:0]  = 6 x 2-bit registers
 *
 * Per-bucket register semantics:
 *   Bit 1 (MSB) = "seen hash at NLZ == absCeiling"
 *   Bit 0 (LSB) = "seen hash at NLZ == absCeiling - 1"
 *   State 00 = not observed at ceiling or ceiling-1
 *   State 01 = observed at ceiling-1 only
 *   State 10 = observed at ceiling only
 *   State 11 = observed at both levels
 *
 * Update rule (hashAndStore):
 *   delta = hashNLZ - absCeiling
 *   delta > 0  : raise ceiling. delta==1: MSBs to LSBs. delta>=2: wipe.
 *                Set MSB for arriving bucket.
 *   delta == 0 : set MSB (ceiling hit)
 *   delta == -1: set LSB (one below ceiling)
 *   delta <= -2: ignore (too far below ceiling)
 *   Overflow: if newLocalExp > 15, ignore hash entirely.
 *
 * Two-level exponent (identical mechanism to FLL2):
 *   globalCeiling = shared minimum ceiling across all words.
 *   absCeiling = globalCeiling + localExp.
 *   advanceGlobal fires when all words have localExp >= 1.
 *
 * Estimation: state-table CF (84 equivalence classes, same as FLL2),
 * plus DLC/HC/LDLC via variable-observation-count framework.
 *
 * @author Brian Bushnell, Chloe von Einzbern-Bushnell
 * @date April 9, 2026
 */
public final class BankedCeilingLogLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int BUCKETS_PER_WORD = 6;
	/** Bits 0,2,4,6,8,10 — LSB of each bucket (ceiling-1 level). */
	static final int LSB_MASK = 0x555;
	/** Bits 1,3,5,7,9,11 — MSB of each bucket (ceiling level). */
	static final int MSB_MASK = 0xAAA;
	static final int HISTORY_MASK = 0xFFF;
	static final int NUM_TIERS = 16;
	/** Number of equivalence classes for 6 buckets with 4 states = C(9,3). */
	static final int NUM_EQUIV = 84;
	/** Prefix sums for idx84 computation. */
	static final int[] OFFSET_A = {0, 28, 49, 64, 74, 80, 83};
	private static final boolean debugging = false;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	BankedCeilingLogLog() { this(2048, 31, -1, 0); }

	BankedCeilingLogLog(parse.Parser p) {
		super(p);
		numWords = roundUpToWords(buckets);
		modBuckets = numWords * BUCKETS_PER_WORD;
		words = new short[numWords];
		numLocalZeros = numWords;
	}

	BankedCeilingLogLog(int buckets_, int k_, long seed, float minProb_) {
		super(nextPow2(roundUpToWords(buckets_) * BUCKETS_PER_WORD), k_, seed, minProb_);
		numWords = roundUpToWords(buckets_);
		modBuckets = numWords * BUCKETS_PER_WORD;
		words = new short[numWords];
		numLocalZeros = numWords;
	}

	@Override
	public BankedCeilingLogLog copy() {
		return new BankedCeilingLogLog(modBuckets, k, -1, minProb);
	}

	private static int roundUpToWords(int buckets_) {
		return (buckets_ + BUCKETS_PER_WORD - 1) / BUCKETS_PER_WORD;
	}

	private static int nextPow2(int n) {
		return Integer.highestOneBit(Math.max(1, n - 1)) << 1;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Core Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number) {
		final long key = Tools.hash64shift(number ^ hashXor);

		// Early exit: skip hashes below globalCeiling - 1
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
		final int absCeiling = globalCeiling + localExp;
		final int delta = hashNLZ - absCeiling;

		if (delta <= -2) { return; } // too far below ceiling

		if (delta > 0) {
			// Hash exceeds ceiling — must raise it
			final int newLocalExp = localExp + delta;
			if (newLocalExp > 15) { return; } // overflow: ignore entirely
			final boolean wasZeroExp = (localExp == 0);
			if (delta == 1) {
				// Shift by 1: MSBs become LSBs (old ceiling becomes ceiling-1)
				word = ((word & MSB_MASK) >>> 1) | (newLocalExp << 12);
			} else {
				// delta >= 2: wipe all history, only this bucket survives
				word = (newLocalExp << 12);
			}
			// Set MSB for the arriving bucket
			word |= (1 << (1 + register * 2));
			lastCardinality = -1;
			lastEstimates = null;
			words[wordIdx] = (short) word;
			// Update numLocalZeros tracking
			if (wasZeroExp) {
				numLocalZeros--;
				if (numLocalZeros == 0) { advanceGlobal(); }
			}
		} else {
			// delta == 0: set MSB (ceiling hit)
			// delta == -1: set LSB (one below ceiling)
			final int bitToSet = (delta == 0)
				? (1 << (1 + register * 2))   // MSB
				: (1 << (register * 2));       // LSB
			if ((word & bitToSet) != 0) { return; } // already set
			lastCardinality = -1;
			lastEstimates = null;
			word |= bitToSet;
			words[wordIdx] = (short) word;
		}

		assert !debugging || repOK()
			: "repOK failed after hashAndStore, word=" + wordIdx;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Global Advance         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * When all words have localExp >= 1, increment globalCeiling,
	 * decrement all localExps, recount numLocalZeros, update eeMask.
	 * Identical trigger and mechanism to FLL2's advanceGlobal.
	 */
	private void advanceGlobal() {
		while (numLocalZeros == 0) {
			globalCeiling++;
			numLocalZeros = 0;
			for (int i = 0; i < numWords; i++) {
				int w = words[i] & 0xFFFF;
				int le = (w >>> 12) & 0xF;
				assert le >= 1 : "localExp=0 but numLocalZeros was 0, word " + i;
				le--;
				words[i] = (short) ((le << 12) | (w & HISTORY_MASK));
				if (le == 0) { numLocalZeros++; }
			}
			// Reject hashes with NLZ < globalCeiling - 1
			eeMask = (globalCeiling <= 1) ? -1L : (-1L) >>> (globalCeiling - 1);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Equivalence Class Index   ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts a 12-bit history to its 84-slot equivalence class index. */
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
		double cfEst = rawCFEstimate();
		double microEst = microCardinalityDouble();
		double hybridEst = blendEstimate(cfEst, microEst);
		long card = (long) hybridEst;
		card = Math.max(card, 0);
		card = Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality = card;
		return card;
	}

	/**
	 * Raw CF-table estimate with CV-weighted averaging.
	 *
	 * weight_i = (1/CV_i)^CV_POWER.  High-CV states (recently wiped)
	 * get downweighted; low-CV states (stable, informative) dominate.
	 * est = sum(value_i * weight_i) / sum(weight_i)
	 */
	public static double CV_POWER = 2.0;

	private double rawCFEstimate() {
		if (CF_TABLE_TIERS <= 0) { return 0; }
		double weightedSum = 0;
		double weightSum = 0;
		for (int i = 0; i < numWords; i++) {
			final int w = words[i] & 0xFFFF;
			if (w == 0) { continue; }
			final int localExp = (w >>> 12) & 0xF;
			final int history = w & HISTORY_MASK;
			final int absCeiling = localExp + globalCeiling;
			final int tier = Math.min(absCeiling, CF_TABLE_TIERS - 1);
			final int eq = idx84(history);
			final double value = tierCardinality(absCeiling) * cfTable[tier][eq];
			final double cv = (cvTable != null && tier < cvTable.length)
				? Math.max(cvTable[tier][eq], 0.001) : 1.0;
			final double weight = Math.pow(1.0 / cv, CV_POWER);
			weightedSum += value * weight;
			weightSum += weight;
		}
		if (weightSum <= 0) { return 0; }
		double est = (weightedSum / weightSum) * numWords;
		est *= TERMINAL_CORRECTION;
		if (cardCFKeys != null && cardCFKeys.length > 0) {
			double keyScale = (cardCFBuckets > 0 && modBuckets != cardCFBuckets)
				? (double) cardCFBuckets / modBuckets : 1.0;
			double cf = CorrectionFactor.getCF(est, est,
				cardCFValues, cardCFKeys, 5, 0.0001, keyScale);
			est *= cf;
		}
		return est;
	}

	/*--------------------------------------------------------------*/
	/*----------------       LC Estimate            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * LC estimate: occupancy-based, accurate at very low cardinality.
	 * Counts buckets with state 00 (empty) as unoccupied.
	 */
	private double lcEstimate() {
		int emptyBuckets = 0;
		for (int i = 0; i < numWords; i++) {
			final int w = words[i] & 0xFFFF;
			final int history = w & HISTORY_MASK;
			for (int r = 0; r < BUCKETS_PER_WORD; r++) {
				if (((history >>> (r * 2)) & 0x3) == 0) { emptyBuckets++; }
			}
		}
		if (emptyBuckets == 0) { return modBuckets; }
		return -modBuckets * Math.log((double) emptyBuckets / modBuckets);
	}

	/*--------------------------------------------------------------*/
	/*----------------       DLC Estimates          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Collect per-tier "filled at tier T" counts using standard DLC model.
	 *
	 * B = modBuckets for every tier (all buckets contribute).
	 * For each word at ceiling C, per bucket:
	 *   T > C:  bucket is empty at tier T (ceiling would be >= T otherwise)
	 *   T == C: bucket is filled if MSB set (seen NLZ == C == T)
	 *   T < C:  bucket is filled if nonZero (seen NLZ >= C-1 >= T)
	 *           00 buckets treated as empty (conservative; may have seen [T..C-2])
	 *
	 * DLC_T = 2^T * LC(modBuckets, filled_T) — standard DLC formula.
	 */
	private int[] dlcFilled;
	private int dlcTierCount;

	private void collectTierCounts() {
		int maxCeiling = globalCeiling;
		for (int i = 0; i < numWords; i++) {
			int c = globalCeiling + ((words[i] >>> 12) & 0xF);
			if (c > maxCeiling) { maxCeiling = c; }
		}
		dlcTierCount = maxCeiling + 1;
		if (dlcFilled == null || dlcFilled.length < dlcTierCount) {
			dlcFilled = new int[dlcTierCount];
		}
		java.util.Arrays.fill(dlcFilled, 0, dlcTierCount, 0);

		for (int i = 0; i < numWords; i++) {
			final int w = words[i] & 0xFFFF;
			final int C = globalCeiling + ((w >>> 12) & 0xF);
			final int history = w & HISTORY_MASK;

			int nonZero = 0, msbSet = 0;
			for (int r = 0; r < BUCKETS_PER_WORD; r++) {
				final int regBits = (history >>> (r * 2)) & 0x3;
				if (regBits != 0) { nonZero++; }
				if ((regBits & 2) != 0) { msbSet++; }
			}

			// Tier C: filled += msbSet (seen NLZ == C)
			if (C < dlcTierCount) {
				dlcFilled[C] += msbSet;
			}
			// Tiers < C: filled += nonZero (seen NLZ >= C-1 >= T for all T < C)
			for (int t = 0; t < C && t < dlcTierCount; t++) {
				dlcFilled[t] += nonZero;
			}
		}
	}

	/** Standard DLC at tier T: 2^T * LC(modBuckets, filled_T). */
	private double dlcAtTier(int filled, int tier) {
		final int V = modBuckets - filled;
		if (V <= 0 || V >= modBuckets) { return -1; }
		return (1L << tier) * (-modBuckets * Math.log((double) V / modBuckets));
	}

	/**
	 * DLC estimate: pick the tier with occupancy closest to target fraction.
	 */
	private double dlcEstimate() {
		collectTierCounts();
		double bestEst = 0;
		double bestDist = Double.MAX_VALUE;
		final double target = AbstractCardStats.DLC_TARGET_FRAC;

		for (int t = 0; t < dlcTierCount; t++) {
			final int V = modBuckets - dlcFilled[t];
			if (V <= 0 || V >= modBuckets) { continue; }
			final double frac = (double) V / modBuckets;
			final double dist = Math.abs(Math.log(frac) - Math.log(target));
			if (dist < bestDist) {
				bestDist = dist;
				bestEst = dlcAtTier(dlcFilled[t], t);
			}
		}
		return bestEst;
	}

	/**
	 * DLCBest: pick the tier with maximum information (B * binary entropy).
	 */
	private double dlcBestEstimate() {
		double bestEst = 0;
		double bestInfo = -1;
		for (int t = 0; t < dlcTierCount; t++) {
			final int V = modBuckets - dlcFilled[t];
			if (V <= 0 || V >= modBuckets) { continue; }
			final double frac = (double) V / modBuckets;
			final double info = -modBuckets * (frac * Math.log(frac) + (1 - frac) * Math.log(1 - frac));
			if (info > bestInfo) {
				bestInfo = info;
				bestEst = dlcAtTier(dlcFilled[t], t);
			}
		}
		return bestEst;
	}

	/*--------------------------------------------------------------*/
	/*----------------       HC Estimate            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * HC: information-weighted geometric mean of per-tier DLC estimates.
	 */
	private double hcEstimate() {
		double logSum = 0;
		double weightSum = 0;
		for (int t = 0; t < dlcTierCount; t++) {
			final int V = modBuckets - dlcFilled[t];
			if (V <= 0 || V >= modBuckets) { continue; }
			final double est = dlcAtTier(dlcFilled[t], t);
			if (est <= 0) { continue; }
			final double frac = (double) V / modBuckets;
			final double info = -modBuckets * (frac * Math.log(frac) + (1 - frac) * Math.log(1 - frac));
			logSum += info * Math.log(est);
			weightSum += info;
		}
		if (weightSum <= 0) { return 0; }
		return Math.exp(logSum / weightSum);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Micro Estimate         ----------------*/
	/*--------------------------------------------------------------*/

	/** MicroIndex-based very-low-cardinality estimate. */
	private double microCardinalityDouble() {
		final int bits = Long.bitCount(microIndex);
		double micro = (bits >= 64) ? 64.0 : -64.0 * Math.log(1.0 - bits / 64.0);
		int totalSetBits = 0;
		for (int i = 0; i < numWords; i++) {
			totalSetBits += Integer.bitCount(words[i] & HISTORY_MASK);
		}
		return Math.max(micro, totalSetBits);
	}

	/** Expected cardinality at a given tier (from simulator tierAvg). */
	private static double tierCardinality(int absCeiling) {
		if (tierAvg != null && absCeiling < tierAvg.length && absCeiling >= 0) {
			return tierAvg[absCeiling];
		}
		if (tierAvg == null || tierAvg.length < 2) { return 1.0; }
		double lastGrowth = tierAvg[tierAvg.length - 1] / tierAvg[tierAvg.length - 2];
		double base = tierAvg[tierAvg.length - 1];
		int extra = absCeiling - (tierAvg.length - 1);
		return base * Math.pow(lastGrowth, extra);
	}

	private static final double BLEND_LO = 8.0;
	private static final double BLEND_HI = 100.0;
	private static final double BLEND_LOG_LO = Math.log(BLEND_LO);
	private static final double BLEND_LOG_RANGE = Math.log(BLEND_HI) - Math.log(BLEND_LO);

	/** Blend micro and CF estimates (100% micro below 8, 100% CF above 100). */
	private double blendEstimate(double cfEst, double microEst) {
		final int microBits = Long.bitCount(microIndex);
		double key = (microBits < 58) ? microEst : cfEst;
		key = Math.max(key, 1.0);
		if (key <= BLEND_LO) { return microEst; }
		if (key >= BLEND_HI) { return cfEst; }
		double frac = (Math.log(key) - BLEND_LOG_LO) / BLEND_LOG_RANGE;
		return microEst * (1.0 - frac) + cfEst * frac;
	}

	@Override
	public double[] rawEstimates() {
		if (lastEstimates != null && lastCardinality >= 0) { return lastEstimates; }
		final double cfEst = rawCFEstimate();
		final double microEst = microCardinalityDouble();
		final double hybridEst = blendEstimate(cfEst, microEst);
		final double lc = lcEstimate();
		collectTierCounts();
		final double dlc = dlcEstimate();
		final double dlcBest = dlcBestEstimate();
		final double hc = hcEstimate();
		final int total = 17 + AbstractCardStats.NUM_DLC_TIERS + AbstractCardStats.NUM_EXTRA;
		final double[] r = new double[total];
		java.util.Arrays.fill(r, hybridEst);
		r[0] = cfEst;      // Mean = CF table estimate
		r[5] = lc;          // LC = occupancy-based
		r[6] = hybridEst;   // Hybrid = blended CF+micro
		r[10] = dlc;        // DLCPure
		r[11] = dlc;        // DLC
		r[13] = dlcBest;    // DLCBest
		r[9] = hc;          // LCmin slot = HC estimate
		lastEstimates = r;
		return r;
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
			if (le == 0) { actualZeros++; }
		}
		assert actualZeros == numLocalZeros
			: "numLocalZeros mismatch: tracked=" + numLocalZeros + " actual=" + actualZeros;
		assert numLocalZeros > 0 : "numLocalZeros is 0 but advanceGlobal was not called";
		long expectedMask = (globalCeiling <= 1) ? -1L : (-1L) >>> (globalCeiling - 1);
		assert eeMask == expectedMask
			: "eeMask inconsistent with globalCeiling=" + globalCeiling;
		return true;
	}

	public double occupancy() {
		int filled = 0;
		for (int i = 0; i < numWords; i++) {
			if ((words[i] & 0xFFFF) != 0) { filled++; }
		}
		return (double) filled / numWords;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Merge/Misc          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void add(CardinalityTracker log) {
		throw new UnsupportedOperationException("BCLL merge not yet implemented");
	}

	@Override
	public final float[] compensationFactorLogBucketsArray() { return null; }

	public int getNumWords() { return numWords; }
	public int getModBuckets() { return modBuckets; }
	public int getGlobalCeiling() { return globalCeiling; }
	public int getWord(int i) { return words[i] & 0xFFFF; }

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final short[] words;
	private final int numWords;
	private final int modBuckets;
	/** Shared minimum ceiling across all words. */
	private int globalCeiling = 0;
	/** Count of words with localExp == 0. */
	private int numLocalZeros;
	/** Early-exit mask: reject hashes with NLZ < globalCeiling - 1. */
	private long eeMask = -1L;
	private double[] lastEstimates;

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	/** Flat correction multiplier — will be calibrated separately from FLL2. */
	public static double TERMINAL_CORRECTION = 1.0;

	static double[][] cfTable;
	/** Per-state CV (coefficient of variation) for inverse-variance weighting. */
	static double[][] cvTable;
	static int CF_TABLE_TIERS = 0;
	static double[] tierAvg;

	static float[] cardCFValues;
	static float[] cardCFKeys;
	static int cardCFBuckets = 0;

	static final String STATE_TABLE_FILE = "stateTableBCLL.tsv.gz";
	static final String CF_FILE = "cardinalityCorrectionBCLL.tsv.gz";

	/** MUST appear after all static field declarations. */
	static {
		loadCFTable();
		loadCardCFTable();
	}

	public static synchronized void loadCFTable() {
		if (cfTable != null) { return; }
		String path = Data.findPath("?" + STATE_TABLE_FILE);
		if (path == null) {
			setFallbackCF();
			return;
		}
		FileFormat ff = FileFormat.testInput(path, null, false);
		if (ff == null) { setFallbackCF(); return; }
		ByteFile bf = ByteFile.makeByteFile(ff, 1);
		if (bf == null) { setFallbackCF(); return; }

		java.util.ArrayList<double[]> multRows = new java.util.ArrayList<>();
		java.util.ArrayList<double[]> cvRows = new java.util.ArrayList<>();
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
				final String rowType = parts[1];
				if (!"mult".equals(rowType) && !"cv".equals(rowType)) { continue; }

				double[] row = new double[NUM_EQUIV];
				double fill = "cv".equals(rowType) ? 1.0 : 1.0;
				java.util.Arrays.fill(row, fill);
				for (int i = 2; i < parts.length && (i - 2) < NUM_EQUIV; i++) {
					try { row[i - 2] = Double.parseDouble(parts[i]); }
					catch (NumberFormatException e) { /* leave as default */ }
				}
				java.util.ArrayList<double[]> target = "cv".equals(rowType) ? cvRows : multRows;
				while (target.size() <= tier) { target.add(null); }
				target.set(tier, row);
			}
		}
		bf.close();

		tierAvg = new double[avgList.size()];
		for (int i = 0; i < avgList.size(); i++) { tierAvg[i] = avgList.get(i); }

		CF_TABLE_TIERS = multRows.size();
		cfTable = new double[CF_TABLE_TIERS][NUM_EQUIV];
		for (int t = 0; t < CF_TABLE_TIERS; t++) {
			if (multRows.get(t) != null) { cfTable[t] = multRows.get(t); }
			else { java.util.Arrays.fill(cfTable[t], 1.0); }
		}
		cvTable = new double[cvRows.size()][NUM_EQUIV];
		for (int t = 0; t < cvRows.size(); t++) {
			if (cvRows.get(t) != null) { cvTable[t] = cvRows.get(t); }
			else { java.util.Arrays.fill(cvTable[t], 1.0); }
		}

		System.err.println("Loaded BCLL CF table: " + CF_TABLE_TIERS + " tiers, "
			+ tierAvg.length + " tier averages, " + NUM_EQUIV + " equiv classes"
			+ ", " + cvTable.length + " CV tiers");
	}

	public static synchronized void loadCardCFTable() {
		if (cardCFKeys != null) { return; }
		String path = Data.findPath("?" + CF_FILE);
		if (path == null) { return; }
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
			float cf = Float.parseFloat(parts[1]);
			keys.add(key);
			vals.add(cf);
		}
		bf.close();

		cardCFKeys = keys.toArray();
		cardCFValues = vals.toArray();
		System.err.println("Loaded BCLL cardinality CF table: " + cardCFKeys.length
			+ " entries, buckets=" + cardCFBuckets);
	}

	private static void setFallbackCF() {
		cfTable = new double[1][NUM_EQUIV];
		java.util.Arrays.fill(cfTable[0], 1.0);
		tierAvg = new double[]{1.0};
		CF_TABLE_TIERS = 1;
	}
}
