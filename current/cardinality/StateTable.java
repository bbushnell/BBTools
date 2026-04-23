package cardinality;

/**
 * Static utility for sub-NLZ state lookups and per-state correction factors.
 * Single owner of ALL per-state correction data:
 * - Additive NLZ corrections for Mean/HMean (tierMult on difSum)
 * - SBS correction table (per-state LC estimates)
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

	/** Steady-state per-state CFs for 1-bit history, standard 2x tiers (2 states).
	 *  Regenerated 2026-04-12 from mc2_hist1_128m_32k.txt tier 11 LinCF. */
	static final double[] CF_HISTORY_1={-1.35708218, +0.16616068};
	/** Steady-state per-state CFs for 1-bit history, CTLL 2*sqrt(2) tiers (2 states).
	 *  Regenerated 2026-04-19 from MC2 mode=ctll bits=1 outer=2M inner=32k, tier 7 LinCF. */
	static final double[] CF_CTLL_1={-2.22885953, +0.06747291};
	/** Steady-state per-state CFs for 2-bit history (4 states).
	 *  Geometric mean averaging (all/geo model). See CardinalityGuide Section 14.
	 *  Regenerated 2026-04-12 from mc2_hist2_128m_32k.txt tier 11 GeoCF. */
	static final double[] CF_HISTORY_2={-2.45693840, -1.03241245, -1.48290584, +0.33548548};
	/** Steady-state per-state CFs for 3-bit history (8 states).
	 *  Regenerated 2026-04-12 from mc2_hist3_128m_32k.txt tier 11 LinCF. */
	static final double[] CF_HISTORY_3={-3.53486901, -2.26872423, -2.70578580, -1.10730867, -2.84418969, -1.56163360, -2.13998342, +0.21177784};

	/** Per-tier CFs for 1-bit history, standard 2x tiers (tiers 0-1; tier 2+ uses steady-state).
	 *  Regenerated 2026-04-12 from mc2_hist1 LinCF rows. */
	static final double[][] CF_HISTORY_1_TIERS={
		{+0.00000000, +0.00000000},
		{-1.90683416, +0.19744068},
	};
	/** Per-tier CFs for 1-bit history, CTLL 2*sqrt(2) tiers (tiers 0-2; tier 3+ uses steady-state).
	 *  Regenerated 2026-04-19 from MC2 mode=ctll bits=1 outer=2M inner=32k LinCF rows. */
	static final double[][] CF_CTLL_1_TIERS={
		{+0.00000000, +0.00000000},
		{-2.92159826, +0.07409516},
		{-2.43085655, +0.06963695},
	};
	/** Steady-state per-state CFs for 2-bit history, CTLL 2*sqrt(2) tiers (4 states).
	 *  Generated 2026-04-19 from MC2 mode=ctll bits=2 outer=2M inner=32k, tier 7 LinCF. */
	static final double[] CF_CTLL_2={-3.86427049, -2.11815111, -3.06401825, +0.07339904};
	/** Per-tier CFs for 2-bit history, CTLL 2*sqrt(2) tiers (tiers 0-2; tier 3+ uses steady-state).
	 *  Generated 2026-04-19 from MC2 mode=ctll bits=2 outer=2M inner=32k LinCF rows. */
	static final double[][] CF_CTLL_2_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-2.92159826, +0.00000000, +0.07409516, +0.00000000},
		{-4.76961496, -2.30133622, -3.48567704, +0.07572406},
	};
	/** Per-tier CFs for 2-bit history (tiers 0-2; tier 3+ uses steady-state).
	 *  Geometric mean averaging (all/geo model). Regenerated 2026-04-12 from mc2_hist2 GeoCF rows. */
	static final double[][] CF_HISTORY_2_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-1.64516934, +0.00000000, +0.32902894, +0.00000000},
		{-2.91470444, -1.11691781, -1.63200820, +0.37491432},
	};
	/** Per-tier CFs for 3-bit history (tiers 0-3; tier 4+ uses steady-state).
	 *  Regenerated 2026-04-12 from mc2_hist3 LinCF rows. */
	static final double[][] CF_HISTORY_3_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-1.90683416, +0.00000000, +0.00000000, +0.00000000, +0.19744068, +0.00000000, +0.00000000, +0.00000000},
		{-3.26667625, +0.00000000, -1.33649006, +0.00000000, -1.95068638, +0.00000000, +0.22284749, +0.00000000},
		{-4.43057341, -2.53750546, -3.10739120, -1.18796076, -3.29890601, -1.69584646, -2.37783659, +0.22051950},
	};

	/** Steady-state CFs for 3-bit history, ETLL expanded tiers (halfNlz/4, TIER_SCALE=2.0).
	 *  Generated 2026-04-21 from MC2 mode=etll bits=3 outer=32k inner=32k, tier 5 LinCF. */
	static final double[] CF_ETLL_3_LIN={-0.17376441, -0.08848045, -0.05907605, -0.01104290, -0.03955570, -0.01740046, -0.01647446, +0.22568811};
	static final double[] CF_ETLL_3_GEO={-1.70572680, -0.65826663, -0.77383387, +0.12971612, -0.83497522, +0.02286040, -0.11328315, +0.82508404};
	static final double[] CF_ETLL_3_HARM={-1.67788945, -0.37901570, -0.48891880, +0.49768142, -0.54745193, +0.39651958, +0.26607104, +1.22967348};
	static final double[] CF_ETLL_3_ENTRY_GEO={-1.83930155, -0.29355145, -0.35220248, +0.63639596, -0.37945578, +0.58750517, +0.50182279, +1.35479159};
	static final double[] CF_ETLL_3_ENTRY_LIN={-1.74017876, -0.54315205, -0.60724981, +0.27158162, -0.63799416, +0.21883980, +0.12691903, +0.93748722};
	static final double[] CF_ETLL_3_ENTRY_HARM={-1.66251479, +0.37036598, +0.32095537, +1.48341377, +0.29568706, +1.43707870, +1.35959655, +2.28083573};
	static final double[] CF_ETLL_3_BOTH_GEO={-1.38437572, -0.10733116, -0.16778306, +0.69275752, -0.20328467, +0.65660996, +0.57945412, +1.20607043};
	/** 0=all/lin, 1=all/geo, 2=all/harm, 3=entry/geo, 4=entry/lin, 5=entry/harm, 6=both/geo */
	static int ETLL_HSB_MODE=1;
	static double[] CF_ETLL_3_OVERRIDE=null;
	static double[] etll3Table(){
		if(CF_ETLL_3_OVERRIDE!=null){return CF_ETLL_3_OVERRIDE;}
		switch(ETLL_HSB_MODE){
			case 1: return CF_ETLL_3_GEO;
			case 2: return CF_ETLL_3_HARM;
			case 3: return CF_ETLL_3_ENTRY_GEO;
			case 4: return CF_ETLL_3_ENTRY_LIN;
			case 5: return CF_ETLL_3_ENTRY_HARM;
			case 6: return CF_ETLL_3_BOTH_GEO;
			default: return CF_ETLL_3_LIN;
		}
	}
	/** Per-tier CFs for 3-bit history, ETLL expanded tiers (tiers 0-3; tier 4+ uses steady-state).
	 *  Generated 2026-04-22 from MC2 mode=etll bits=3 outer=2M inner=32k GeoCF rows.
	 *  Correct expanded mapping: tier=halfNlz, TIER_SCALE=0.5. all/geo model. */
	static final double[][] CF_ETLL_3_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-0.76128382, +0.00000000, +0.00000000, +0.00000000, +0.61321318, +0.00000000, +0.00000000, +0.00000000},
		{-1.43993655, +0.00000000, -0.08369830, +0.00000000, -0.21651589, +0.00000000, +0.83980362, +0.00000000},
		{-2.06795280, -0.72201024, -0.84480485, +0.18037009, -0.92528004, +0.06492381, -0.08881447, +0.92675886},
	};

	static final double[] CF_ETLL_2={-0.24650375, -0.05903061, -0.05085254, +0.22081960};
	static final double[][] CF_ETLL_2_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-0.93583304, +0.00000000, +0.46929410, +0.00000000},
		{-1.70521485, -0.30657758, -0.46663286, +0.63089713},
	};
	static final double[] CF_ETLL_1={-0.27190525, +0.19555480};
	static final double[][] CF_ETLL_1_TIERS={
		{+0.00000000, +0.00000000},
		{-0.93989503, +0.46832201},
	};

	/** AUDLL32: 3-state collapsed history (10+11→state 2), standard tiers.
	 *  State 2 emitted as: binary 10 at tier 1 (only high bit valid),
	 *  binary 11 at tier 2+ (both bits valid).
	 *  Positions: [0]=00(no obs), [1]=01(obs -2), [2]=10(state 2 at tier 1),
	 *  [3]=11(state 2 at tier 2+).
	 *  Generated from MC2 mode=history bits=2 outer=64k inner=32k, collapsed via
	 *  CF=log2((P(10)*2^CF(10)+P(11)*2^CF(11))/(P(10)+P(11))). April 2026. */
	static final double[] CF_AUDLL32_2={-2.48357976, -1.15357362, +0.16555091, +0.16555091};
	static final double[][] CF_AUDLL32_2_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-1.92572967, +0.00000000, +0.19631649, +0.00000000},
		{-3.26384437, -1.34192589, +0.18142796, +0.18142796},
	};
	/** When true, historyOffset uses AUDLL32's collapsed 3-state tables instead of CF_HISTORY_2. */
	static boolean USE_AUDLL32_HSB=false;

	/**
	 * Additive NLZ correction for a bucket with the given history state.
	 * For nlzBin below the tier table depth, returns per-tier correction.
	 * For nlzBin at or above, returns steady-state correction.
	 */
	static double historyOffset(int nlzBin, int hbits, int histPattern){
		final double[] steadyState;
		final double[][] tierTables;
		if(USE_TTLL_HSB && hbits==4){
			steadyState=(CF_TTLL_4_OVERRIDE!=null ? CF_TTLL_4_OVERRIDE : CF_TTLL_4);
			tierTables=CF_TTLL_4_TIERS;
		}else if(USE_AUDLL32_HSB && hbits==2){
			steadyState=CF_AUDLL32_2; tierTables=CF_AUDLL32_2_TIERS;
		}else if(hbits==1 && AbstractCardStats.TIER_SCALE>1.0){
			steadyState=CF_CTLL_1; tierTables=CF_CTLL_1_TIERS;
		}else if(hbits==2 && AbstractCardStats.TIER_SCALE>1.0){
			steadyState=CF_CTLL_2; tierTables=CF_CTLL_2_TIERS;
		}else if(hbits==1 && AbstractCardStats.TIER_SCALE<=0.5){
			steadyState=CF_ETLL_1; tierTables=CF_ETLL_1_TIERS;
		}else if(hbits==2 && AbstractCardStats.TIER_SCALE<=0.5){
			steadyState=CF_ETLL_2; tierTables=CF_ETLL_2_TIERS;
		}else if(hbits==3 && AbstractCardStats.TIER_SCALE<=0.5){
			steadyState=etll3Table(); tierTables=CF_ETLL_3_TIERS;
		}else if(hbits==1){steadyState=CF_HISTORY_1; tierTables=CF_HISTORY_1_TIERS;}
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

	/** Regenerated 2026-04-12 from mc2_mantissa / mc2_andtissa / mc2_nlz2 tier 11 LinCF. */
	static final double[] CF_MANTISSA_2={-0.47966150, -0.27170868, -0.02938129, +0.25891253};
	static final double[] CF_ANDTISSA_2={-0.32168520, -0.00131953, +0.21278764, +0.36323342};
	static final double[] CF_NLZ2_2={-0.35721832, -0.02938129, +0.17131840, +0.32547686};

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

	/** Regenerated 2026-04-12 from mc2_luck1/2/3 tier 11 LinCF. */
	static final double[] CF_LUCK_1={+0.16616068, -1.35708218};
	static final double[] CF_LUCK_2={+0.16616068, -1.15851633, -2.26872423, -3.53486901};
	static final double[] CF_LUCK_3={+0.16616068, -1.15851633, -2.26872423, -3.31759640, -4.34098592, -5.35666575, -6.37171780, -7.62134559};
	static final double[][] CF_LUCK_2_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{+0.19744068, +0.00000000, +0.00000000, -1.90683416},
		{+0.18054181, -1.33649006, +0.00000000, -3.26667625},
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

	/** Combined history+mantissa CFs (2h+2m, 16 states). State = mant | (hist << mbits).
	 *  Regenerated 2026-04-12 from mc2_hm2h2m tier 11 LinCF. */
	static final double[] CF_HISTMANT_2H2M={
		-2.55732601, -2.50936461, -2.46139657, -2.41237047,
		-1.31475843, -1.23006318, -1.13878523, -1.04007904,
		-1.73542236, -1.67827942, -1.61748433, -1.55602545,
		-0.21215887, -0.03748874, +0.16876559, +0.41789590};
	/** Combined history+mantissa CFs (1h+3m, 16 states).
	 *  Regenerated 2026-04-12 from mc2_hm1h3m tier 11 LinCF. */
	static final double[] CF_HISTMANT_1H3M={
		-1.55575619, -1.50946173, -1.46135009, -1.41214737, -1.35975175, -1.30762743, -1.25121743, -1.19382385,
		-0.31476344, -0.23001984, -0.13899047, -0.04048309, +0.06608254, +0.18164435, +0.30998730, +0.44847093};
	/** Combined history+mantissa CFs (3h+1m, 16 states).
	 *  Regenerated 2026-04-12 from mc2_hm3h1m tier 11 LinCF. */
	static final double[] CF_HISTMANT_3H1M={
		-3.55948071, -3.51219274,
		-2.31526672, -2.23000908,
		-2.73511936, -2.68014689,
		-1.21230225, -1.03710071,
		-2.87065536, -2.82069927,
		-1.62583161, -1.51363378,
		-2.17482525, -2.11118531,
		-0.09849967, +0.33145235};
	/** Per-tier CFs for 1h+3m (tiers 0-1). Regenerated 2026-04-12 from mc2_hm1h3m LinCF rows. */
	static final double[][] CF_HISTMANT_1H3M_TIERS={
		{-0.90676987, -0.72547070, -0.54171649, -0.35407978, -0.16142129, +0.03894549, +0.24951120, +0.47378379,
		 +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-2.27608278, -2.18571147, -2.09544112, -2.00388602, -1.91227373, -1.82033943, -1.72689046, -1.63187553,
		 -0.39387537, -0.28798832, -0.17491036, -0.05437216, +0.07349372, +0.21304422, +0.36411701, +0.53212391},
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
		{0, 1, 3, 7, 15, 31},// hbits=4: bin0=1, bin1=2, bin2=4, bin3=8, bin4=16, bin5+=16
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
	/*----------------    TwinTail Corrections       ----------------*/
	/*--------------------------------------------------------------*/

	/** Steady-state per-state CFs for TTLL symmetric encoding (mode 0).
	 *  State = (h1<<2)|h0.  Both tails bit-filtered by histBit.
	 *  States 0,1,4,5 are unreachable (CF=0).
	 *  Regenerated 2026-04-12 from mc2_twintail tier 11 LinCF (128m trials). */
	static final double[] CF_TWINTAIL_SYM={
		+0.00000000, +0.00000000, -1.48012967, -0.63565089,
		+0.00000000, +0.00000000, -0.63674532, +0.16588297,
		-1.47986386, -0.63602662, -0.78159160, -0.05659339,
		-0.63578397, +0.16664983, -0.05689921, +0.77854334,
	};

	/** Per-tier CFs for TTLL symmetric encoding (tiers 0-2).
	 *  Regenerated 2026-04-12 from mc2_twintail LinCF rows (128m trials). */
	static final double[][] CF_TWINTAIL_SYM_TIERS={
		{+0.00000000, +0.00000000, -0.58478361, +0.00000000,
		 +0.00000000, +0.00000000, +0.00000000, +0.00000000,
		 -0.58451724, +0.00000000, +0.73703846, +0.00000000,
		 +0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{+0.00000000, +0.00000000, -2.12939990, -0.81293009,
		 +0.00000000, +0.00000000, -0.81316157, +0.19126012,
		 -2.12887589, -0.81291704, -1.01388945, -0.07608401,
		 -0.81298650, +0.19145997, -0.07598090, +0.90112914},
		{+0.00000000, +0.00000000, -1.74317881, -0.71717685,
		 +0.00000000, +0.00000000, -0.71702835, +0.17460671,
		 -1.74326858, -0.71696719, -0.88543067, -0.06957253,
		 -0.71721475, +0.17468025, -0.06872554, +0.83884921},
	};

	/** HSB table for TTLL 4-bit combined_h, used by historyOffset() for Mean+H.
	 *  16 states: (h1<<2)|h0.
	 *  Regenerated 2026-04-23 from TTLLSimulator tier 8 LinCF (1M iters, 128 threads).
	 *  LinCF beat GeoCF (6.70% vs 17.41% Mean+H WidthWt at 16 DDLs). */
	static final double[] CF_TTLL_4={
		-1.49192135, -0.64982726, -1.49483430, -0.65015221,
		-0.65268397, +0.16042073, -0.65340451, +0.15885110,
		-1.48904227, -0.64950669, -0.80067450, -0.07292658,
		-0.65196680, +0.16198683, -0.07020482, +0.78354314,
	};
	/** Per-tier HSB for TTLL 4-bit combined_h (tiers 0-2).
	 *  Regenerated 2026-04-23 from TTLLSimulator LinCF (1M iters, 128 threads). */
	static final double[][] CF_TTLL_4_TIERS={
		{-0.58451970, -0.99828173, -0.58371294, -0.58371294,
		 -0.99828173, +0.00000000, -0.58371294, +0.00000000,
		 -0.58530616, -0.58530616, +0.73738649, +0.73738649,
		 -0.58530616, +0.00000000, +0.73738649, +0.00000000},
		{-2.12722297, -0.81130651, -2.12715875, -0.81153296,
		 -0.80889054, +0.19196520, -0.80881554, +0.19217468,
		 -2.12728725, -0.81108105, -1.01257778, -0.07569150,
		 -0.80896571, +0.19175376, -0.07561258, +0.90129211},
		{-1.74574997, -0.71609843, -1.74623535, -0.72073856,
		 -0.71981759, +0.17556117, -0.72047585, +0.17661725,
		 -1.74526699, -0.71149179, -0.88063150, -0.05687208,
		 -0.71916583, +0.17450372, -0.06881533, +0.83955537},
	};
	/** When true, historyOffset uses TTLL's 16-state HSB tables for hbits==4. */
	static boolean USE_TTLL_HSB=false;
	/** Command-line override for CF_TTLL_4 steady-state table (16 values). Null = use compiled table. */
	static double[] CF_TTLL_4_OVERRIDE=null;

	/** Steady-state per-state CFs for TTLL master/slave encoding (mode 1).
	 *  h0=master (pure NLZ history), h1=slave (bit-filtered).
	 *  States 4,5,12,13 are unreachable (h1 LSB implies h0 MSB).
	 *  Regenerated 2026-04-12 from mc2_masterslave tier 11 LinCF (128m trials). */
	static final double[] CF_TWINTAIL_MS={
		-2.53243964, -1.26893998, -1.84145164, -0.56098203,
		+0.00000000, +0.00000000, -1.62385196, +0.00258584,
		-2.43577229, -1.08405746, -1.73657328, -0.31227071,
		+0.00000000, +0.00000000, -1.49880296, +0.43378515,
	};

	/** Per-tier CFs for TTLL master/slave encoding (tiers 0-2).
	 *  Regenerated 2026-04-12 from mc2_masterslave LinCF rows (128m trials). */
	static final double[][] CF_TWINTAIL_MS_TIERS={
		{-0.58478361, +0.00000000, +0.00000000, +0.00000000,
		 +0.00000000, +0.00000000, +0.00000000, +0.00000000,
		 +0.22241615, +0.00000000, +0.00000000, +0.00000000,
		 +0.00000000, +0.00000000, +0.00000000, +0.00000000},
		{-2.12939990, +0.00000000, -0.81293009, +0.00000000,
		 +0.00000000, +0.00000000, -0.05217232, +0.00000000,
		 -1.75990503, +0.00000000, -0.46519935, +0.00000000,
		 +0.00000000, +0.00000000, +0.48143329, +0.00000000},
		{-3.36643895, -1.47360186, -2.23516850, -0.63168742,
		 +0.00000000, +0.00000000, -1.93217830, -0.00252344,
		 -3.18475962, -1.24541261, -2.08772666, -0.35144443,
		 +0.00000000, +0.00000000, -1.76880942, +0.47012507},
	};

	/** Additive NLZ correction for a TTLL bucket with the given combined_h state.
	 *  Selects symmetric or master/slave tables based on TwinTailLogLog.ENCODING_MODE.
	 *  For absNlz < 3, returns per-tier correction.  For absNlz >= 3, returns steady-state. */
	static double ttllOffset(int absNlz, int combinedH){
		final double[] ss=(TwinTailLogLog.ENCODING_MODE==0) ? CF_TWINTAIL_SYM : CF_TWINTAIL_MS;
		final double[][] tiers=(TwinTailLogLog.ENCODING_MODE==0) ? CF_TWINTAIL_SYM_TIERS : CF_TWINTAIL_MS_TIERS;
		if(tiers!=null && absNlz<tiers.length){
			return tiers[absNlz][combinedH];
		}
		return ss[combinedH];
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
	/** Override: set to non-zero to use this value instead of the default terminalCF. */
	static double terminalCFOverride=0;

	static double terminalCF(int hbits, int mbits){
		if(terminalCFOverride!=0){return terminalCFOverride;}
		if(mbits==0){
			if(hbits==1){
				// CTLL compressed tiers: 2^(-CF_CTLL_1[1]) = 2^(-0.0663) ≈ 0.9551
				if(AbstractCardStats.TIER_SCALE>1.0){return Math.pow(2.0, -CF_CTLL_1[1]);}
				return 0.84237;
			}
			if(hbits==2){return 0.772605900;} // geo model terminal CF
			if(hbits==3){return 0.97751;}
		}
		return 1.0;
	}

	/*--------------------------------------------------------------*/

	private StateTable(){throw new AssertionError("No instances");}

}
