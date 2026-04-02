package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * UltraLogLog8: 6-bit NLZ + 2-bit sub-NLZ history in 8-bit registers.
 * <p>
 * Extends the LL6 (HLL) baseline by tracking whether elements with
 * NLZ = (max-1) and NLZ = (max-2) were ever observed. These 2 history
 * bits refine the per-bucket contribution to the harmonic sum using
 * theoretically derived correction factors, reducing estimator variance.
 * <p>
 * Two history modes are supported:
 * <ul>
 *   <li><b>Frozen</b> (FROZEN_HISTORY=true, default): History bits are locked
 *       at the moment of tier promotion and never updated afterward. This
 *       enables full early-exit optimization (elements below max are pure
 *       no-ops). Correction factors reflect the entry-moment distribution.</li>
 *   <li><b>Active</b> (FROZEN_HISTORY=false): History bits continue to be
 *       updated as elements with NLZ = max-1 or max-2 arrive. This captures
 *       more information but requires processing elements that don't update
 *       the max. Correction factors reflect the time-averaged distribution.</li>
 * </ul>
 * <p>
 * Register encoding (8 bits):
 * <ul>
 *   <li>Bits 7-2: nlzStored = absNlz + 1 (1-63, or 0 = empty)</li>
 *   <li>Bit 1: sawMinus2 (observed NLZ = max - 2)</li>
 *   <li>Bit 0: sawMinus1 (observed NLZ = max - 1)</li>
 * </ul>
 * <p>
 * History carry on promotion uses a shift register:
 * {@code newHist = ((oldHist | 0x4) >> k) & 0x3} where k is the jump distance.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class UltraLogLog8 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	UltraLogLog8(){
		this(2048, 31, -1, 0);
	}

	UltraLogLog8(Parser p){
		super(p);
		maxArray=new byte[buckets];
	}

	UltraLogLog8(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
	}

	@Override
	public UltraLogLog8 copy(){return new UltraLogLog8(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Per-bucket summarize, modeled on DDL8's approach.
	 * Computes hllSumFilled (integer NLZ) for standard estimators,
	 * and hllSumFilledM (state-corrected) for the refined estimator.
	 * Uses the active CF table selected by FROZEN_HISTORY toggle.
	 */
	private CardinalityStats summarize(){
		final double[] baseCF=FROZEN_HISTORY ? STATE_CF_FROZEN : STATE_CF_ACTIVE;
		// Build per-tier multiplier tables: tiers 0-2 have their own CFs and offsets,
		// tier 3+ uses the steady-state CFs with STATE_CF_OFFSET.
		final double[][] tierMult=new double[64][4];
		for(int t=0; t<64; t++){
			final double[] cf;
			final double offset;
			if(t==0){cf=TIER0_CF; offset=TIER0_CF_OFFSET;}
			else if(t==1){cf=TIER1_CF; offset=TIER1_CF_OFFSET;}
			else if(t==2){cf=TIER2_CF; offset=TIER2_CF_OFFSET;}
			else{cf=baseCF; offset=STATE_CF_OFFSET;}
			for(int s=0; s<4; s++){
				tierMult[t][s]=Math.pow(2.0, -(cf[s]*STATE_POWER+offset));
			}
		}
		final float[][] cfMatrix=FROZEN_HISTORY ? CF_MATRIX_FROZEN : CF_MATRIX_ACTIVE;
		final int cfBuckets=FROZEN_HISTORY ? CF_BUCKETS_FROZEN : CF_BUCKETS_ACTIVE;

		double difSum=0;
		double hllSumFilled=0;
		double gSum=0;
		int count=0;
		final int[] nlzCounts=new int[66];

		for(int i=0; i<buckets; i++){
			final int stored=maxArray[i]&0xFF;
			if(stored>0){
				final int absNlz=(stored>>2)-1;
				final int state=stored&0x3;
				if(absNlz>=0 && absNlz<65){nlzCounts[absNlz+1]++;}
				final double m=tierMult[Math.min(absNlz, 63)][state];
				final double base=Math.pow(2.0, -absNlz)*m;
				final double dif=(absNlz==0 ? (double)Long.MAX_VALUE : (absNlz<64 ? (double)(1L<<(63-absNlz)) : 1.0))*m;
				difSum+=dif;
				hllSumFilled+=base;
				gSum+=Math.log(Tools.max(1, dif));
				count++;
			}
		}
		nlzCounts[0]=buckets-count;
		lastRawNlz=nlzCounts;
		lastCorrNlz=nlzCounts;
		return new CardinalityStats(difSum, hllSumFilled, hllSumFilled,
		                            gSum, count, buckets, null, cfMatrix, cfBuckets,
		                            CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys, microIndex,
		                            nlzCounts, 0);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardinalityStats s=summarize();
		final double rawHyb=s.hybridDLL();
		long card=(long)(rawHyb);
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((UltraLogLog8)log);
	}

	/** Merges another UltraLogLog8 via per-bucket max with history carry. */
	public void add(UltraLogLog8 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				mergeRegister(i, log.maxArray[i]&0xFF);
			}
		}
	}

	/**
	 * Merges a single register value into this sketch.
	 * Takes the register with the higher NLZ. On ties, OR the history bits.
	 * When the other has higher NLZ, carry my history via shift register.
	 */
	private void mergeRegister(int bucket, int otherStored){
		final int myStored=maxArray[bucket]&0xFF;
		if(otherStored<=0){return;}
		if(myStored==0){maxArray[bucket]=(byte)otherStored; return;}

		final int myNlz=(myStored>>2);
		final int otherNlz=(otherStored>>2);

		if(myNlz>otherNlz){
			final int k=myNlz-otherNlz;
			final int otherHist=otherStored&0x3;
			final int shifted=((otherHist|0x4)>>k)&0x3;
			maxArray[bucket]=(byte)(myStored|shifted);
		}else if(otherNlz>myNlz){
			final int k=otherNlz-myNlz;
			final int myHist=myStored&0x3;
			final int shifted=((myHist|0x4)>>k)&0x3;
			maxArray[bucket]=(byte)(otherStored|shifted);
		}else{
			maxArray[bucket]=(byte)(myStored|otherStored);
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int newNlzStored=Math.min(nlz+1, 63); // absNlz+1, clamped
		final int oldStored=maxArray[bucket]&0xFF;
		final int oldNlzStored=(oldStored>0) ? (oldStored>>2) : 0;

		if(FROZEN_HISTORY){
			// Frozen mode: only promotions change the register.
			// Elements at or below current max are pure early exits.
			if(newNlzStored<=oldNlzStored){return;}
		}else{
			// Active mode: elements below max can set history bits.
			if(newNlzStored<oldNlzStored){
				final int diff=oldNlzStored-newNlzStored;
				if(diff==1){
					maxArray[bucket]=(byte)(oldStored|0x02); // set sawMinus1 (bit1)
				}else if(diff==2){
					maxArray[bucket]=(byte)(oldStored|0x01); // set sawMinus2 (bit0)
				}
				return;
			}
			if(newNlzStored==oldNlzStored){return;} // same NLZ, no change
		}

		// New NLZ is higher: promote with shift-register history carry
		lastCardinality=-1;
		if(oldStored==0){filledBuckets++;}

		final int k=newNlzStored-oldNlzStored;
		// Shift register carry: bit2 = "old max existed" (only if bucket was non-empty)
		final int oldHist=(oldStored>0) ? (oldStored&0x3) : 0;
		final int carry=(oldStored>0) ? 0x4 : 0; // don't set "old max existed" for empty buckets
		final int newHist=((oldHist|carry)>>k)&0x3;

		maxArray[bucket]=(byte)((newNlzStored<<2)|newHist);
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(s.hybridDLL());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * One byte per bucket: bits 7-2 = nlzStored (absNlz+1),
	 * bit 1 = sawMinus2, bit 0 = sawMinus1.
	 */
	final byte[] maxArray;
	private int filledBuckets=0;

	int[] lastRawNlz, lastCorrNlz;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Toggle between frozen and active history modes.
	 * Frozen: history locked at promotion, enables full early exit.
	 * Active: history updated by sub-max elements, more information but slower.
	 */
	public static boolean FROZEN_HISTORY=false;

	/**
	 * Power applied to state corrections. 1.0 = full correction,
	 * 0.5 = sqrt (half the additive CF in NLZ-space), 0.0 = no correction.
	 */
	public static double STATE_POWER=1.0;

	/**
	 * Additive offset applied to CF values before conversion to multipliers.
	 * Default -0.7929 was determined empirically by minimizing the coefficient
	 * of variation (std / |1+meanErr|) of the Hybrid estimator via quadratic
	 * fit across offsets -1.2 to -0.3, using 16384 estimators, 512 buckets,
	 * maxmult=20, cf=f. The minimum Hybrid CV was at offset -0.7929.
	 * Applied to the empirically derived STATE_CF_ACTIVE values.
	 */
	public static double STATE_CF_OFFSET=-0.7929;

	/*--------------------------------------------------------------*/
	/*----------------     State Correction Tables     ----------------*/
	/*--------------------------------------------------------------*/

	/*
	 * DERIVATION OF STATE CORRECTION FACTORS
	 *
	 * Each bucket stores 2 history bits tracking whether NLZ values
	 * max-1 and max-2 were observed (via shift register carry on promotion).
	 * The 4 states (00, 01, 10, 11) encode how "established" the current
	 * NLZ is: state 11 = saw both neighbors (well-established), state 00 =
	 * saw neither (likely outlier from a large jump).
	 *
	 * The additive NLZ corrections (STATE_CF_*) were derived theoretically
	 * via numerical integration of E[log2(u) | state, max=X] where
	 * u = n / 2^X (elements per bucket normalized by expected count at
	 * this NLZ level). Two density functions were used:
	 *
	 *   FROZEN:  P(u) ∝ e^(-u)
	 *            Entry-moment density — state locked at promotion instant.
	 *
	 *   ACTIVE:  P(u | max=X) ∝ (e^(-u/2) - e^(-u)) / u
	 *            Time-averaged density over the bucket's lifetime at tier X.
	 *
	 * The theoretical CF values were then refined empirically: a uniform
	 * additive offset of -1.0 was found to minimize the coefficient of
	 * variation (see STATE_CF_OFFSET documentation above).
	 *
	 * The effective CF at runtime is: baseCF[state] * STATE_POWER + STATE_CF_OFFSET
	 * The effective multiplier is: 2^(-effectiveCF)
	 *
	 * Tiers 0, 1, 2 use per-tier correction tables (STATE_CF_TIER0/1/2)
	 * because the steady-state derivation does not apply to tiers where
	 * NLZ-1 or NLZ-2 cannot exist (e.g., at NLZ=0, state is always 00).
	 *
	 * Bit convention (MSB = nearer neighbor):
	 *   bit1 (0x02) = sawMinus1 (saw NLZ = max-1, the nearer neighbor)
	 *   bit0 (0x01) = sawMinus2 (saw NLZ = max-2, the farther neighbor)
	 *   state = (sawMinus1 << 1) | sawMinus2
	 *
	 * Index order: [state 00, state 01, state 10, state 11]
	 *   00 = saw neither
	 *   01 = saw max-2 only (farther only; unusual — typically from k=2 jump)
	 *   10 = saw max-1 only (nearer only; common at entry from k=1 jump)
	 *   11 = saw both (well-established)
	 */

	/**
	 * Additive NLZ corrections for frozen history (entry-moment density).
	 * Strictly ascending: [00=neither, 01=far only, 10=near only, 11=both].
	 * Negative values reduce the effective NLZ (lowering the cardinality
	 * estimate for lucky outlier buckets); positive values increase it
	 * (raising the estimate for well-established buckets).
	 */
	/** Theoretical CF for frozen history. TODO: re-derive empirically. */
	static final double[] STATE_CF_FROZEN={-1.976212, -0.350372, -0.010539, +1.257391};

	/**
	 * Empirical additive NLZ corrections for active history.
	 * Derived from single-bucket simulation (6.5M trials × 32k hashes).
	 * For each (tier, state), measured avg cardinality at observation time;
	 * CF = log2(stateAvgCard / tierAvgCard). Stable across tiers 3-11.
	 *
	 * Index: [00=neither, 01=saw NLZ-2 only, 10=saw NLZ-1 only, 11=both].
	 * Note: CF[10] < CF[01] because seeing the near neighbor (NLZ-1, which
	 * is 2× rarer than NLZ-2) WITHOUT the far neighbor indicates a stronger
	 * outlier than seeing the far without the near.
	 *
	 * Per-tier corrections for tiers 0-2 (non-steady-state):
	 *   Tier 0: {0.00, N/A,   N/A,   N/A  }  (only state 00 possible)
	 *   Tier 1: {-1.86, N/A,  +0.21, N/A  }  (only states 00 and 10 possible)
	 *   Tier 2: {-3.19, -1.29, -1.95, +0.23}
	 *   Tier 3+: use STATE_CF_ACTIVE below
	 */
	static final double[] STATE_CF_ACTIVE={-2.50813368, -1.16962885, -1.64933633, +0.20806313};

	/**
	 * Multipliers applied to each bucket's harmonic sum contribution.
	 * Computed as 2^(-CF). Larger values inflate the contribution to the
	 * harmonic sum, which REDUCES the cardinality estimate (estimate ∝ 1/sum).
	 */
	static final double[] STATE_MULT_FROZEN={3.934587, 1.274889, 1.007332, 0.418300};

	/** Empirical multipliers for active history mode (2^(-CF)). */
	static final double[] STATE_MULT_ACTIVE={5.68883675, 2.24953818, 3.13689302, 0.86569868};

	/*
	 * PER-TIER CORRECTIONS FOR LOW TIERS
	 *
	 * At NLZ=0: only state 00 is possible (no NLZ=-1 or NLZ=-2 to observe).
	 * At NLZ=1: states 00 and 10 are possible (can see NLZ=-1... wait, NLZ=0
	 *   exists, so sawMinus1 can be set via shift carry from a previous tier,
	 *   but sawMinus2 requires NLZ=-1 which doesn't exist).
	 * At NLZ=2: all 4 states are possible but the distribution is non-steady-state
	 *   because the bucket has only been through 0-2 promotions.
	 *
	 * TODO: Derive per-tier corrections empirically from simulation data.
	 * For now, tiers 0-2 use multiplier 1.0 (no correction).
	 */
	/**
	 * Per-tier CFs for tiers 0-2 where steady-state assumptions don't hold.
	 * Empirically derived from single-bucket simulation (6.5M outer × 32k inner).
	 * Index: [00, 01, 10, 11]. N/A states use 0.0 (multiplier=1.0, no correction).
	 *
	 * Tier 0: only state 00 possible (bucket just filled, no history).
	 * Tier 1: only states 00 and 10 possible (can see NLZ=0 but not NLZ=-1).
	 * Tier 2: all 4 states possible but non-steady-state distribution.
	 */
	static final double[] TIER0_CF={0.00000000, 0.00000000, 0.00000000, 0.00000000};
	static final double[] TIER1_CF={-1.84448832, 0.00000000, +0.20510426, 0.00000000};
	static final double[] TIER2_CF={-3.17437230, -1.29310079, -1.92165476, +0.22896780};

	/** Per-tier offsets. Tier 0 = 0 (no correction meaningful).
	 *  Tier 1/2 = -0.25 (approximate; not yet optimized). */
	static final double TIER0_CF_OFFSET=0.0;
	static final double TIER1_CF_OFFSET=-0.25;
	static final double TIER2_CF_OFFSET=-0.25;

	/*--------------------------------------------------------------*/
	/*----------------       CF Table Loading       ----------------*/
	/*--------------------------------------------------------------*/

	/** CF table file for frozen history mode. */
	public static final String CF_FILE_FROZEN="?cardinalityCorrectionULL8-FRZ.tsv.gz";
	/** CF table file for active history mode. */
	public static final String CF_FILE_ACTIVE="?cardinalityCorrectionULL8-ACT.tsv.gz";
	// TODO: Replace CF_FILE_ACTIVE with "?cardinalityCorrectionULL8Active.tsv.gz"
	// once a separate table is generated for active mode.

	/** Frozen mode CF table. */
	private static int CF_BUCKETS_FROZEN=2048;
	private static float[][] CF_MATRIX_FROZEN=initializeCF_Frozen(CF_BUCKETS_FROZEN);

	/** Active mode CF table. */
	private static int CF_BUCKETS_ACTIVE=2048;
	private static float[][] CF_MATRIX_ACTIVE=initializeCF_Active(CF_BUCKETS_ACTIVE);

	public static float[][] initializeCF_Frozen(int buckets){
		CF_BUCKETS_FROZEN=buckets;
		return CF_MATRIX_FROZEN=CorrectionFactor.loadFile(CF_FILE_FROZEN, buckets);
	}

	public static float[][] initializeCF_Active(int buckets){
		CF_BUCKETS_ACTIVE=buckets;
		return CF_MATRIX_ACTIVE=CorrectionFactor.loadFile(CF_FILE_ACTIVE, buckets);
	}

	/** Reinitialize both CF tables for a new bucket count. */
	public static void initializeCF(int buckets){
		initializeCF_Frozen(buckets);
		initializeCF_Active(buckets);
	}

}
