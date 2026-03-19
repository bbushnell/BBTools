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
		final double[] mult=FROZEN_HISTORY ? STATE_MULT_FROZEN : STATE_MULT_ACTIVE;
		final float[][] cfMatrix=FROZEN_HISTORY ? CF_MATRIX_FROZEN : CF_MATRIX_ACTIVE;
		final int cfBuckets=FROZEN_HISTORY ? CF_BUCKETS_FROZEN : CF_BUCKETS_ACTIVE;

		double difSum=0;
		double hllSumFilled=0;
		double gSum=0;
		int count=0;
		final int[] nlzCounts=new int[64];

		for(int i=0; i<buckets; i++){
			final int stored=maxArray[i]&0xFF;
			if(stored>0){
				final int absNlz=(stored>>2)-1;
				final int state=stored&0x3;
				if(absNlz>=0 && absNlz<64){nlzCounts[absNlz]++;}
				final double m=(absNlz>=MIN_CORRECTION_TIER ? mult[state] : 1.0);
				final double base=Math.pow(2.0, -absNlz)*m;
				final double dif=(absNlz==0 ? (double)Long.MAX_VALUE : (absNlz<64 ? (double)(1L<<(63-absNlz)) : 1.0))*m;
				difSum+=dif;
				hllSumFilled+=base;
				gSum+=Math.log(Tools.max(1, dif));
				count++;
			}
		}
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
		long card=(long)s.hybridDLL();
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
		if(USE_MICRO){
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
					maxArray[bucket]=(byte)(oldStored|0x01); // set sawMinus1
				}else if(diff==2){
					maxArray[bucket]=(byte)(oldStored|0x02); // set sawMinus2
				}
				return;
			}
			if(newNlzStored==oldNlzStored){return;} // same NLZ, no change
		}

		// New NLZ is higher: promote with shift-register history carry
		lastCardinality=-1;
		if(oldStored==0){filledBuckets++;}

		final int k=newNlzStored-oldNlzStored;
		final int oldHist=(oldStored>0) ? (oldStored&0x3) : 0;
		final int newHist=((oldHist|0x4)>>k)&0x3;

		maxArray[bucket]=(byte)((newNlzStored<<2)|newHist);
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(Math.max(s.hybridDLL(), s.microCardinality()));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * One byte per bucket: bits 7-2 = nlzStored (absNlz+1),
	 * bit 1 = sawMinus2, bit 0 = sawMinus1.
	 */
	private final byte[] maxArray;
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
	 * Minimum absNlz tier for applying state corrections.
	 * Tiers below this lack meaningful history (e.g., absNlz=0
	 * is always state 00 since there's no NLZ=-1 to observe).
	 * Buckets at tiers 0, 1, 2 use multiplier 1.0 (no correction).
	 */
	static final int MIN_CORRECTION_TIER=3;

	/**
	 * State multipliers for frozen history mode.
	 * Derived from entry-moment density P(u) ∝ e^(-u).
	 * Index = (sawMinus2 << 1) | sawMinus1.
	 * CF values: state 00 = -1.976, 01 = -0.350, 10 = -0.011, 11 = +1.257
	 */
	static final double[] STATE_MULT_FROZEN={
		3.934587,  // state 00: saw neither; outlier NLZ, inflate contribution
		1.274889,  // state 01: saw max-1 only (k=2 jump, old max = new max-2)
		1.007332,  // state 10: saw max-2 only (k=1 jump, nearly neutral)
		0.418300   // state 11: saw both; well-established, reduce contribution
	};

	/**
	 * State multipliers for active history mode.
	 * Derived from lifetime density P(u|max=X) ∝ (e^(-u/2) - e^(-u)) / u.
	 * CF values: state 00 = -2.379, 01 = -0.729, 10 = -0.338, 11 = +1.133
	 */
	static final double[] STATE_MULT_ACTIVE={
		5.202266,  // state 00: saw neither; rare survivor, strong outlier
		1.657683,  // state 01: saw max-1 only
		1.263895,  // state 10: saw max-2 only
		0.455923   // state 11: saw both; well-established
	};

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
