package cardinality;

/**
 * Static utility for sub-NLZ state lookups and per-state correction factors.
 * Single owner of ALL per-state correction data:
 * - Additive NLZ corrections for Mean/HMean (tierMult on difSum)
 * - LCHist correction table (per-state LC estimates)
 * - Terminal CF constants (empirical steady-state Mean bias)
 *
 * No instance state. The caller applies CF_OFFSET and constructs
 * tierMult = 2^(-(cf + offset)).
 *
 * State space for 2-bit history (11 reachable states):
 * NLZ=0: {00} = 1 state. NLZ=1: {00,10} = 2 states.
 * NLZ=2: {00,01,10,11} = 4 states. NLZ>=3: same as NLZ=2.
 * Total reachable: 3*2^h - 1 (11 for hbits=2).
 *
 * @author Brian Bushnell, Ady, Chloe
 * @date March 2026
 */
public final class StateTable {

	/*--------------------------------------------------------------*/
	/*----------------    History Corrections       ----------------*/
	/*--------------------------------------------------------------*/

	/** Steady-state per-state CFs for 1-bit history (2 states). */
	static final double[] CF_HISTORY_1={-1.37003787, +0.16864767};
	/** Steady-state per-state CFs for 2-bit history (4 states). */
	static final double[] CF_HISTORY_2={-2.50813368, -1.16962885, -1.64933633, +0.20806313};
	/** Steady-state per-state CFs for 3-bit history (8 states). */
	static final double[] CF_HISTORY_3={-3.58851053, -2.29626211, -2.72632606, -1.11858788, -2.86374753, -1.57703255, -2.14003646, +0.21427608};

	/** Per-tier CFs for 2-bit history (tiers 0-2; tier 3+ uses steady-state). */
	static final double[][] CF_HISTORY_2_TIERS={
		{+0.00000000,  0.00000000,  0.00000000,  0.00000000},
		{-1.84448832,  0.00000000, +0.20510426,  0.00000000},
		{-3.17437230, -1.29310079, -1.92165476, +0.22896780},
	};
	/** Per-tier CFs for 3-bit history (tiers 0-3; tier 4+ uses steady-state). */
	static final double[][] CF_HISTORY_3_TIERS={
		{0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000},
		{-1.84448832, 0.00000000, 0.00000000, 0.00000000, +0.20510426, 0.00000000, 0.00000000, 0.00000000},
		{-3.17437230, 0.00000000, -1.29310079, 0.00000000, -1.92165476, 0.00000000, +0.22896780, 0.00000000},
		{-4.34036999, -2.51859936, -3.06284399, -1.22271738, -3.25845403, -1.66944315, -2.36112478, +0.22352032},
	};

	/**
	 * Additive NLZ correction for a bucket with the given history state.
	 * For nlzBin below the tier table depth, returns per-tier correction.
	 * For nlzBin at or above, returns steady-state correction.
	 */
	static double historyOffset(int nlzBin, int hbits, int histPattern){
		final double[] steadyState;
		final double[][] tierTables;
		if(hbits==1){steadyState=CF_HISTORY_1; tierTables=null;}
		else if(hbits==2){steadyState=CF_HISTORY_2; tierTables=CF_HISTORY_2_TIERS;}
		else if(hbits==3){steadyState=CF_HISTORY_3; tierTables=CF_HISTORY_3_TIERS;}
		else{return 0;}
		if(tierTables!=null && nlzBin<tierTables.length){
			return tierTables[nlzBin][histPattern];
		}
		return steadyState[histPattern];
	}

	/*--------------------------------------------------------------*/
	/*----------------    Mantissa Corrections      ----------------*/
	/*--------------------------------------------------------------*/

	static final double[] CF_MANTISSA_2={-0.4928, -0.2633, -0.0342, +0.2640};
	static final double[] CF_ANDTISSA_2={-0.3295, +0.0017, +0.2241, +0.3588};
	static final double[] CF_NLZ2_2={-0.3562, -0.0342, +0.1648, +0.3382};

	/** Additive NLZ correction for a bucket with the given mantissa value. */
	static double mantissaOffset(int mbits, int mantissaPattern){
		if(mbits==2){return CF_MANTISSA_2[mantissaPattern];}
		return 0;
	}

	/** Additive NLZ correction for andtissa mode. */
	static double andtissaOffset(int mbits, int pattern){
		if(mbits==2){return CF_ANDTISSA_2[pattern];}
		return 0;
	}

	/** Additive NLZ correction for nlz2 mode. */
	static double nlz2Offset(int mbits, int pattern){
		if(mbits==2){return CF_NLZ2_2[pattern];}
		return 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Luck Corrections        ----------------*/
	/*--------------------------------------------------------------*/

	static final double[] CF_LUCK_1={+0.1686, -1.3682};
	static final double[] CF_LUCK_2={+0.1686, -1.1662, -2.3051, -3.6058};
	static final double[] CF_LUCK_3={+0.1686, -1.1662, -2.3051, -3.3914, -4.4281, -5.5035, -6.5513, -7.9929};
	static final double[][] CF_LUCK_2_TIERS={
		{ 0.0000,  0.0000,  0.0000, +0.0000},
		{+0.2051,  0.0000,  0.0000, -1.8472},
		{+0.1855, -1.3051,  0.0000, -3.1800},
	};

	/** Additive NLZ correction for a bucket with the given luck gap. */
	static double luckOffset(int nlzBin, int lbits, int luckPattern){
		final double[] steadyState;
		final double[][] tierTables;
		if(lbits==1){steadyState=CF_LUCK_1; tierTables=null;}
		else if(lbits==2){steadyState=CF_LUCK_2; tierTables=CF_LUCK_2_TIERS;}
		else if(lbits==3){steadyState=CF_LUCK_3; tierTables=null;}
		else{return 0;}
		if(tierTables!=null && nlzBin<tierTables.length){
			return tierTables[nlzBin][luckPattern];
		}
		return steadyState[luckPattern];
	}

	/*--------------------------------------------------------------*/
	/*----------------    Combined Corrections      ----------------*/
	/*--------------------------------------------------------------*/

	/** Combined history+mantissa CFs (2h+2m, 16 states). State = mant | (hist << mbits). */
	static final double[] CF_HISTMANT_2H2M={
		-2.5995, -2.5518, -2.4953, -2.4397,
		-1.3324, -1.2396, -1.1559, -1.0402,
		-1.7472, -1.6938, -1.6204, -1.5905,
		-0.2122, -0.0389, +0.1694, +0.4223};
	/** Combined history+mantissa CFs (1h+3m, 16 states). */
	static final double[] CF_HISTMANT_1H3M={
		-1.5812, -1.5239, -1.4818, -1.4209, -1.3730, -1.3330, -1.2667, -1.1835,
		-0.3085, -0.2349, -0.1372, -0.0452, +0.0691, +0.1802, +0.3047, +0.4605};
	/** Combined history+mantissa CFs (3h+1m, 16 states). */
	static final double[] CF_HISTMANT_3H1M={
		-3.6133, -3.5716,
		-2.3542, -2.2627,
		-2.7558, -2.7068,
		-1.2263, -1.0449,
		-2.8923, -2.8413,
		-1.6402, -1.5328,
		-2.1591, -2.1219,
		-0.0999, +0.3347};
	/** Per-tier CFs for 1h+3m (tiers 0-1). */
	static final double[][] CF_HISTMANT_1H3M_TIERS={
		{-0.8227, -0.6814, -0.5253, -0.3557, -0.1542, +0.0391, +0.2449, +0.4656,
		       0,       0,       0,       0,       0,       0,       0,       0},
		{-2.1487, -2.0697, -2.0043, -1.9300, -1.8541, -1.7696, -1.6953, -1.5896,
		 -0.3752, -0.2654, -0.1551, -0.0362, +0.0784, +0.2230, +0.3644, +0.5282},
	};

	/**
	 * Additive NLZ correction when both history and mantissa bits are present.
	 * Uses simulator-derived combined tables when available, falls back to
	 * historyOffset + mantissaOffset.
	 */
	static double combinedOffset(int nlzBin, int hbits, int mbits,
			int histPattern, int mantPattern){
		final int combinedIdx=mantPattern|(histPattern<<mbits);

		// Check for simulator-derived combined table
		final double[] combined;
		final double[][] combinedTiers;
		if(hbits==2 && mbits==2){combined=CF_HISTMANT_2H2M; combinedTiers=null;}
		else if(hbits==1 && mbits==3){combined=CF_HISTMANT_1H3M; combinedTiers=CF_HISTMANT_1H3M_TIERS;}
		else if(hbits==3 && mbits==1){combined=CF_HISTMANT_3H1M; combinedTiers=null;}
		else{
			// No combined table; fall back to sum
			return historyOffset(nlzBin, hbits, histPattern)+mantissaOffset(mbits, mantPattern);
		}

		// Use per-tier table if available and applicable
		if(combinedTiers!=null && nlzBin<combinedTiers.length){
			return combinedTiers[nlzBin][combinedIdx];
		}

		// Check for outer-product tier tables for configs that have
		// history tier tables but no dedicated combined tier table
		final double[][] htiers=(hbits==2?CF_HISTORY_2_TIERS:hbits==3?CF_HISTORY_3_TIERS:null);
		if(htiers!=null && nlzBin<htiers.length){
			return CF_MANTISSA_2[Math.min(mantPattern, CF_MANTISSA_2.length-1)]+htiers[nlzBin][histPattern];
		}

		// Steady-state combined
		return combined[combinedIdx];
	}

	/*--------------------------------------------------------------*/
	/*----------------     State Enumeration        ----------------*/
	/*--------------------------------------------------------------*/

	/** Precomputed cumulative state offsets per NLZ bin.
	 *  STATE_OFFSETS[hbits][bin] = number of reachable states before this bin.
	 *  hbits=0 is unused (no history). */
	static final int[][] STATE_OFFSETS={
		{},              // hbits=0: no states
		{0, 1},          // hbits=1: bin0=1state, bin1+=2states
		{0, 1, 3, 7},   // hbits=2: bin0=1, bin1=2, bin2=4, bin3+=4
		{0, 1, 3, 7, 15},// hbits=3: bin0=1, bin1=2, bin2=4, bin3=8, bin4+=8
	};

	/** Total reachable states for given sub-NLZ configuration. */
	static int stateCount(int hbits, int mbits, int lbits){
		int count=1;
		if(hbits>0){count*=3*(1<<hbits)-1;}
		if(mbits>0){count*=(1<<mbits);}
		if(lbits>0){count*=(1<<lbits);}
		return count;
	}

	/**
	 * Maps (absNlz, histPattern) to a 0-based index into the reachable
	 * state enumeration. For hbits=2: returns 0-10.
	 */
	static int stateIndex(int absNlz, int hbits, int histPattern){
		if(hbits<=0 || hbits>=STATE_OFFSETS.length){return 0;}
		final int[] offsets=STATE_OFFSETS[hbits];
		final int bin=Math.min(absNlz, offsets.length-1);
		return offsets[bin]+histPattern;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Global Tuning Parameters   ----------------*/
	/*--------------------------------------------------------------*/

	/** Additive offset applied to all per-state CFs: tierMult = 2^(-(cf+CF_OFFSET)).
	 *  Settable by calibration drivers. Default 0. */
	public static double CF_OFFSET=0;

	/*--------------------------------------------------------------*/
	/*----------------    Terminal CF Constants      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Empirically measured terminal Mean CF. At high cardinality all
	 * history bits saturate, so every bucket gets the same tierMult.
	 * This constant captures the convergent bias.
	 * CF_OFFSET is NOT included.
	 */
	static double terminalCF(int hbits, int mbits){
		if(hbits==2 && mbits==0){return 0.929224472;}
		// Other configurations: return 1.0 (no correction) until measured
		return 1.0;
	}

	/*--------------------------------------------------------------*/

	private StateTable(){throw new AssertionError("No instances");}

}
