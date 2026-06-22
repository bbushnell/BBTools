package cardinality;

/**
 * Arithmetic Variable LogLog (AVLL) — self-contained cardinality estimator.
 * Paper companion class with all estimation logic embedded.
 * <p>
 * Packs 11 variable-history registers into each 64-bit word using
 * arithmetic encoding with RADIX=56.  Supports 2-bit history tracking
 * for state-based sampling (SBS).
 * <p>
 * Primary estimator: HLDLC = 0.5 × LDLC + 0.5 × Hybrid+2.
 * Secondary estimator: VLDLC = 0.64 × LDLC + 0.36 × VWMean (disabled).
 * <p>
 * LDLC blends DLC-SBS (dynamic linear counting with state-based
 * sampling fallback) and HC (history counting) across the
 * cardinality range.  VWMean uses per-tier-state expected distinct
 * counts from a precomputed table with harmonic-mean aggregation.
 *
 * @author Brian Bushnell
 * @contributor Noire
 * @date May 2026
 */
public class ArithmeticVariableLogLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	ArithmeticVariableLogLog(){this(2048, 31, -1, 0);}

	/** Parser constructor for DDLCalibrationDriver2 integration. */
	ArithmeticVariableLogLog(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	/**
	 * Creates an AVLL sketch with the specified parameters.
	 * Buckets are rounded up to a multiple of 11 (buckets per word),
	 * then to the next power of two for hash distribution.
	 * @param buckets_ requested number of registers (typical: 2048)
	 * @param k_ k-mer length for sequence hashing
	 * @param seed hash seed (-1 for default)
	 * @param minProb_ minimum base quality probability
	 */
	ArithmeticVariableLogLog(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*11<=0 ? buckets_ :       // parent gets power-of-2 bucket count
			Integer.highestOneBit(roundToWords(buckets_)*11-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*11;
		modBuckets=rounded>0 ? rounded : buckets;            // actual register count
		packedLen=(modBuckets+10)/11;                         // number of 64-bit words
		packed=new long[packedLen];
		floorCount=modBuckets;                               // all registers start at floor
	}

	private static int roundToWords(int b){return Math.max(1, (b+10)/11);}

	@Override public ArithmeticVariableLogLog copy(){return new ArithmeticVariableLogLog(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}
	@Override public float terminalMeanCF(){return TERMINAL_MEAN_CF;}
	@Override public float terminalMeanPlusCF(){return TERMINAL_MEANH_CF;}
	@Override public float hldlcWeight(){return HLDLC_WEIGHT;}
	@Override public float[] compensationFactorLogBucketsArray(){return null;}
	@Override public int bitsPerWord(){return 64;}
	@Override public int bucketsPerWord(){return BPW;}

	/*--------------------------------------------------------------*/
	/*----------------       Register Constants      ----------------*/
	/*--------------------------------------------------------------*/

	/** Base for arithmetic encoding. */
	private static final int RADIX=56;
	/** Buckets per 64-bit word (56^11 < 2^64). */
	private static final int BPW=11;
	/** History tier limit: tiers 2..HTL have 4 states.
	 *  HTL=13 (vs AVDLL64's 14) trades one rare history tier for 4 extra
	 *  non-history exponent levels, improving idempotency (MAX_EXP 18 vs 15). */
	private static final int HTL=13;
	/** Total states with history (T0:1 + T1:2 + T2-13:48 = 51). */
	private static final int HIST_STATES=4*HTL-1;
	/** Maximum exponent (=18, vs AVDLL64's 15). */
	private static final int MAX_EXP=HTL+(RADIX-HIST_STATES);
	/** Maximum register value (=55). */
	private static final int MAX_REGISTER=RADIX-1;
	/** NLZ margin below globalNLZ for floor detection. */
	private static final int HISTORY_MARGIN=2;
	/** Hash word length in bits. */
	private static final int WORDLEN=64;
	/** Number of history bits per register. */
	private static final int HIST_BITS=2;

	/** Powers of RADIX: POW[k] = 56^k. Used for arithmetic encode/decode. */
	private static final long[] POW=computePow();
	private static long[] computePow(){
		long[] p=new long[BPW];
		p[0]=1;
		for(int i=1; i<BPW; i++){p[i]=p[i-1]*RADIX;}
		return p;
	}

	/*--------------------------------------------------------------*/
	/*----------------      State ↔ Exp/Hist Maps    ----------------*/
	/*--------------------------------------------------------------*/

	/** Maps register state (0-55) to exponent (leading-zero count).
	 *  T0: state 0 → exp 0.  T1: states 1-2 → exp 1.
	 *  T2+: states 3-54 → exp (s-3)/4+2.  Non-hist: state 55 → exp 15. */
	static int stateToExp(int state){
		if(state==0){return 0;}                         // empty above floor
		if(state<=2){return 1;}                         // tier 1: 2 states
		if(state<HIST_STATES){return (state-3)/4+2;}    // tiers 2-14: 4 states each
		return HTL+1+(state-HIST_STATES);               // above history: no hist bits
	}

	/** Maps register state (0-55) to 2-bit history pattern.
	 *  T0: 0.  T1: (s-1)<<1 (bit 1 valid, bit 0 unknown).
	 *  T2+: (s-3)&3 (both bits valid).  Non-hist: 0. */
	static int stateToHist(int state){
		if(state==0){return 0;}
		if(state<=2){return (state-1)<<1;}              // tier 1: MSB only
		if(state<HIST_STATES){return (state-3)&3;}      // tiers 2+: both bits
		return 0;                                        // above history tier limit
	}

	/** Encodes (exponent, history) pair into a register state value. */
	static int toState(int exp, int hist){
		if(exp==0){return 0;}
		if(exp==1){return 1+(hist>>1);}                 // only MSB of hist matters at T1
		if(exp<=HTL){return 4*exp-5+hist;}              // tiers 2-14
		return HIST_STATES+(exp-HTL-1);                 // above history: no hist
	}

	/*--------------------------------------------------------------*/
	/*----------------    Register Bitmap Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Unpacks register state into a bitmap where bit k represents NLZ=k.
	 *  History bits are placed at positions (highest-2) and (highest-1). */
	private static long unpackReg(int reg){
		if(reg==0){return 0;}
		final int nlzPart=stateToExp(reg);
		final int hist=stateToHist(reg);
		long bitmap=1L<<nlzPart;                         // leading-zero marker bit
		if(nlzPart>=2){bitmap|=(long)(hist&1)<<(nlzPart-2);}   // history bit 0
		if(nlzPart>=1){bitmap|=(long)((hist>>>1)&1)<<(nlzPart-1);} // history bit 1
		return bitmap;
	}

	/** Packs a bitmap back into a register state.
	 *  Extracts the highest set bit (exponent) and the two bits below it (history). */
	private static int packReg(long bitmap){
		if(bitmap==0){return 0;}
		final int highest=63-Long.numberOfLeadingZeros(bitmap);
		int hist;
		if(highest>=2){
			hist=(int)((bitmap>>>(highest-2))&3);        // both history bits below top
		}else if(highest==1){
			hist=((int)(bitmap&1))<<1;                   // only MSB available
		}else{
			hist=0;                                      // tier 0: no history
		}
		if(highest>MAX_EXP){                             // clamp to maximum exponent
			if(MAX_EXP<=HTL){
				hist=(MAX_EXP>=2) ? (int)((bitmap>>>(MAX_EXP-2))&3) : 0;
				return toState(MAX_EXP, hist);
			}
			return toState(MAX_EXP, 0);
		}
		if(highest>HTL){hist=0;}                         // above HTL: no history stored
		return toState(highest, hist);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads register i from the packed array using unsigned arithmetic.
	 *  Position within word: i%11.  Decode: (word / 56^pos) % 56. */
	int getReg(int i){
		final int wordIdx=i/BPW;
		final int pos=i%BPW;
		final long word=packed[wordIdx];
		return (int)(Long.remainderUnsigned(Long.divideUnsigned(word, POW[pos]), RADIX));
	}

	/** Writes register i.  Uses modular arithmetic: word += (val-old)*56^pos.
	 *  Works because Java's long multiplication is modular (mod 2^64). */
	private void setReg(int i, int val){
		final int wordIdx=i/BPW;
		final int pos=i%BPW;
		long word=packed[wordIdx];
		final int oldVal=(int)(Long.remainderUnsigned(Long.divideUnsigned(word, POW[pos]), RADIX));
		word=word+(long)(val-oldVal)*POW[pos];           // modular update
		packed[wordIdx]=word;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Hash Function        ----------------*/
	/*--------------------------------------------------------------*/

	/** Thomas Wang's 64-bit integer hash (shift-multiply variant).
	 *  Fast, good avalanche properties, invertible. */
	private static long hash64shift(long key){
		key=(~key)+(key<<21);    // key = (key<<21) - key - 1
		key=key^(key>>>24);
		key=(key+(key<<3))+(key<<8); // key * 265
		key=key^(key>>>14);
		key=(key+(key<<2))+(key<<4); // key * 21
		key=key^(key>>>28);
		key=key+(key<<31);
		return key;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Add Element          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Hashes and stores a number into the sketch.
	 * Called by CardinalityTracker.add(long) which also increments the added counter.
	 * @param number the value to hash and store
	 */
	@Override
	public void hashAndStore(long number){
		final long key=hash64shift(number^hashXor);      // hash with parent's XOR seed

		final int idx=(int)(Long.remainderUnsigned(key, modBuckets)); // bucket selection
		final int nlz=Long.numberOfLeadingZeros(key);    // count leading zeros
		final int relNlz=nlz-(globalNLZ+1);              // NLZ relative to current floor

		final int bitPos=relNlz+HISTORY_MARGIN;          // position in register bitmap
		if(bitPos<0 || bitPos>=WORDLEN){return;}         // out of range
		// reads++;

		final int oldReg=getReg(idx);
		int newReg;
		if(bitPos>HTL && bitPos<=MAX_EXP){
			newReg=toState(bitPos, 0);                   // above HTL: no history tracking
		}else{
			long hashPrefix=unpackReg(oldReg);           // decode current bitmap
			hashPrefix|=1L<<bitPos;                      // set the new bit
			newReg=packReg(hashPrefix);                  // re-encode
		}

		if(newReg<=oldReg){return;}                      // register only increases
		// writes++;

		setReg(idx, newReg);
		lastCardinality=-1;                              // invalidate cached result

		if(oldReg==0){filledBuckets++;}                  // track fill count

		final int oldNlzPart=stateToExp(oldReg);
		final int newNlzPart=stateToExp(newReg);
		if(oldNlzPart<=HISTORY_MARGIN && newNlzPart>HISTORY_MARGIN){
			floorCount--;                                // register left the floor zone
			if(floorCount<=0){                           // floor is saturated
				while(floorCount==0 && globalNLZ<WORDLEN){
					globalNLZ++;                         // advance the global NLZ
					floorCount=countAndDecrement();      // demote all registers
				}
			}
		}
	}

	/** Merges another AVLL sketch into this one (register-wise OR). */
	@Override
	public void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((ArithmeticVariableLogLog)log);
	}

	public void add(ArithmeticVariableLogLog log){
		added+=log.added;
		lastCardinality=-1;
		final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
		final int deltaA=newGlobalNLZ-globalNLZ;
		final int deltaB=newGlobalNLZ-log.globalNLZ;
		for(int i=0; i<modBuckets; i++){
			final int rA=getReg(i);
			final int rB=log.getReg(i);
			final long hpA=(rA==0) ? 0 : unpackReg(rA)>>>deltaA;
			final long hpB=(rB==0) ? 0 : unpackReg(rB)>>>deltaB;
			final long merged=hpA|hpB;
			if(merged==0){setReg(i, 0);}
			else{setReg(i, Math.min(packReg(merged), MAX_REGISTER));}
		}
		globalNLZ=newGlobalNLZ;
		filledBuckets=0; floorCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=getReg(i);
			if(reg>0){filledBuckets++;}
			if(reg==0 || stateToExp(reg)<=HISTORY_MARGIN){floorCount++;}
		}
		while(floorCount==0 && globalNLZ<WORDLEN){
			globalNLZ++;
			floorCount=countAndDecrement();
		}
	}

	/** Demotes all registers by one tier (called when globalNLZ advances).
	 *  Returns the new floor count (registers at or below HISTORY_MARGIN). */
	private int countAndDecrement(){
		int newFloorCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){newFloorCount++; continue;}       // already at zero
			final int exp=stateToExp(reg);
			int newReg;
			if(exp>HTL+1){newReg=reg-1;}                 // above history: simple decrement
			else if(exp==HTL+1){newReg=toState(exp-1, 0);} // entering history zone: clear hist
			else if(exp>=3){newReg=reg-4;}               // tiers 2+: shift down one tier (4 states)
			else if(exp==2){newReg=toState(1, stateToHist(reg));} // T2→T1: keep MSB of hist
			else if(exp==1){newReg=0;}                   // T1→T0: becomes empty
			else{newReg=0;}                              // already empty
			setReg(i, newReg);
			if(stateToExp(newReg)<=HISTORY_MARGIN){newFloorCount++;}
		}
		return newFloorCount;
	}

	/*--------------------------------------------------------------*/
	/*----------------     SBS Constants & Formula   ----------------*/
	/*--------------------------------------------------------------*/

	/** SBS state count for 2-bit history: 1+2+4+4 = 11 states. */
	private static final int SBS_STATES=11;

	/** Base values per SBS state: minimum distinct count for this history pattern. */
	private static final int[] SBS_BASE={1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3};

	/** Cubic polynomial coefficients: y = base + a*L + b*L² + c*L³ where L = ln(B/V).
	 *  B-independent (verified on B=2048 and B=256).  @author Brian Bushnell, Ady */
	private static final double[] SBS_A={
		0.22105100, 0.12252876, 0.36952900, 0.06133797, 0.30751800,
		0.18526204, 0.43143072, 0.67720200, 0.70123350, 0.64688227, 0.73638900};
	private static final double[] SBS_B={
		0.03969600, 0.00756061, 0.03269600, 0.00236917, 0.02818300,
		0.00912048, 0.03450029, 0.03911000, 0.01746612, 0.01820189, 0.01219500};
	private static final double[] SBS_C={
		-0.00321100, -0.00065621, -0.00192200, -0.00029760, -0.00176600,
		-0.00068916, -0.00183925, -0.00248900, -0.00112447, -0.00079415, -0.00085800};

	/** Expected distinct elements for one bucket with this SBS state and occupancy.
	 *  @param si   SBS state index (0-10)
	 *  @param filled  number of non-empty buckets
	 *  @param B    total bucket count */
	private static double sbsFormula(int si, int filled, int B){
		final double V=Math.max(0.5, B-filled);          // clamp to avoid singularity
		final double L=Math.log((double)B/V);            // occupancy parameter
		return SBS_BASE[si]+SBS_A[si]*L+SBS_B[si]*L*L+SBS_C[si]*L*L*L;
	}

	/** Maps (absNlz, histPattern) to an SBS state index (0-10) for 2-bit history.
	 *  NLZ bins: 0 → 1 state, 1 → 2 states, 2 → 4 states, 3+ → 4 states.
	 *  Returns -1 for structurally impossible states (nonzero invalid bits). */
	private static int sbsStateIndex(int nlzBin, int histPattern){
		final int bin=Math.min(nlzBin, HIST_BITS+1);     // clamp to 3
		final int validSlots=Math.min(bin, HIST_BITS);   // how many hist bits are valid
		final int invalidBits=HIST_BITS-validSlots;      // how many bits MUST be zero
		if(invalidBits>0 && (histPattern&((1<<invalidBits)-1))!=0){return -1;} // invalid state
		final int offset=(1<<Math.min(bin, HIST_BITS+1))-1; // cumulative state count below this bin
		final int idx=(invalidBits>0) ? (histPattern>>>invalidBits) : histPattern;
		return offset+idx;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Mean CF Formula         ----------------*/
	/*--------------------------------------------------------------*/

	/** 3-Sigmoid + 2-Gaussian CF formula.
	 *  CF = a0 + a1*S(lc,c1,w1) + a2*S(lc,c2,w2) + a3*S(lc,c3,w3)
	 *       + g1*G(lc,gc1,gw1) + g2*G(lc,gc2,gw2)
	 *  where lc = log2(card), S(x,c,w) = (1+tanh((x-c)/w))/2,
	 *  G(x,c,w) = exp(-((x-c)/w)²).
	 *  Terminal (card→∞) = p[0]+p[1]+p[4]+p[7].  @author Eru */
	private static double meanCfFormula(double card, double[] p){
		if(card<=0){return p[0]+p[1]+p[4]+p[7];}        // high-card terminal value
		final double lc=Math.log(card)*INV_LOG2;         // log2(card)
		final double s1=(1+Math.tanh((lc-p[2])/p[3]))*0.5;  // sigmoid 1
		final double s2=(1+Math.tanh((lc-p[5])/p[6]))*0.5;  // sigmoid 2
		final double s3=(1+Math.tanh((lc-p[8])/p[9]))*0.5;  // sigmoid 3
		final double d1=(lc-p[11])/p[12]; final double g1=Math.exp(-d1*d1); // gaussian 1
		final double d2=(lc-p[14])/p[15]; final double g2=Math.exp(-d2*d2); // gaussian 2
		return p[0]+p[1]*s1+p[4]*s2+p[7]*s3+p[10]*g1+p[13]*g2;
	}

	/** 1/ln(2), for converting ln to log2. */
	private static final double INV_LOG2=1.0/Math.log(2.0);

	/** Iterative CF evaluation using card/B as formula input.
	 *  Mirrors the fixed-point iteration in CorrectionFactor.getCF():
	 *  seed at dlcRaw/B, then refine with corrected estimate/B. */
	private double iterativeCF(double seedEst, double rawEst, double[] coeffs){
		final double invB=1.0/modBuckets;
		double cf=meanCfFormula(seedEst*invB, coeffs);
		if(seedEst>10.0*modBuckets){
			final double bestEst=cf*rawEst;
			if(bestEst>0){
				final double newCf=meanCfFormula(bestEst*invB, coeffs);
				if(Math.abs(newCf-cf)>1e-6){cf=newCf;}
			}
		}
		return cf;
	}

	/** Mean CF residual, x=log2(card/B). Bucket-independent.
	 *  Fitted from resource CF table (B=1408). R²=0.9999969, terminal=0.99991. */
	private static final double[] MCF_MEAN={
		-4.134681757262724,  4.999999992726675,-11.794568553347373,  1.166556972002922,
		 0.328887715628532,  1.845414705069709,  1.071814615231428,
		-0.194293984917707, -5.783510608691229,  2.244741892627562,
		 0.185448568786388, -3.798664343259118,  3.684574870998338,
		 0.070445731787030,  1.581614638086199,  1.391853026540059};

	/** MeanH CF residual, x=log2(card/B). Bucket-independent.
	 *  Fitted from resource CF table (B=1408). R²=0.9999969, terminal=0.99963. */
	private static final double[] MCF_MEANH={
		-3.691071292933851, -0.363637261057505,  2.567263769729832,  1.066355956682916,
		 4.999999754294653,-11.189952323115953,  0.785742127708774,
		 0.054343414450697, -8.335873985296551,  1.101742014940036,
		-0.143358513726359, -0.794735788647950,  2.755625941876753,
		-0.187952399977993,  0.812771424938221,  2.145566972133455};

	private static final float TERMINAL_MEAN_CF=0.720952f;
	private static final float TERMINAL_MEANH_CF=0.855563f;

	/*--------------------------------------------------------------*/
	/*----------------        HC CF Formula          ----------------*/
	/*--------------------------------------------------------------*/

	/** HC CF: exp+sin model.  CF(card) = T + 1/(A + B·exp) + S0·sin(S1·exp) + C0·cos(C1·exp)
	 *  where exp = log2(card).  Generic UDLL6-family coefficients. */
	private static final double HC_CF_T =   0.998484777537405;
	private static final double HC_CF_A = 562.709161557896778;
	private static final double HC_CF_B =   1.011821147415617;
	private static final double HC_CF_S0=   0.000344111534286;
	private static final double HC_CF_S1=  -0.000037079356918;
	private static final double HC_CF_C0=   0.000424944596363;
	private static final double HC_CF_C1=  -0.000035494763272;

	private static double hcCfFormula(double card){
		if(card<=0){return HC_CF_T;}                     // terminal
		final double lc=Math.log(card)*INV_LOG2;         // log2(card)
		final double base=HC_CF_T-HC_CF_A*Math.exp(-HC_CF_B*lc);
		final double angle=TWO_PI*lc;
		return base+(HC_CF_S0+HC_CF_S1*lc)*Math.sin(angle)
		           +(HC_CF_C0+HC_CF_C1*lc)*Math.cos(angle);
	}

	/** 2π, cached for HC CF formula. */
	private static final double TWO_PI=2.0*Math.PI;
	/** √(2/π), theoretical DLC tier error constant. */
	private static final double SQRT_2_OVER_PI=Math.sqrt(2.0/Math.PI);

	/*--------------------------------------------------------------*/
	/*----------------     VWMean State Table        ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-tier-state expected distinct counts for VWMean.
	 *  Indexed as VW_TABLE[tier][histPattern], where histPattern is 0-3.
	 *  Beyond VW_MAX_TIER, values scale geometrically by 2× per tier.
	 *  Only tiers 0-4 are needed; higher tiers are exactly 2× the previous
	 *  (verified: T25/T4 cluster comparison shows zero difference).
	 *  Source: perTierStateUDLL6.tsv.gz, entry_lin section, columns 4-7. */
	private static final int VW_MAX_TIER=4;

	private static final double[][] VW_TABLE={
		{1.000000,   0.000000,   0.000000,   0.000000},  // tier 0
		{1.000000,   0.000000,   2.733337,   0.000000},  // tier 1
		{1.000000,   2.688307,   2.245420,   5.837820},  // tier 2
		{2.000006,   5.376626,   4.490847,  11.675610},  // tier 3
		{4.000001,  10.753230,   8.981696,  23.351283},  // tier 4
	};

	/*--------------------------------------------------------------*/
	/*----------------       Blend Constants         ----------------*/
	/*--------------------------------------------------------------*/

	/** DLC blend lower bound: below this V/B fraction, pure DLC. */
	private static final float DLC_BLEND_LO=0.088f;
	/** DLC blend upper bound: above this V/B fraction, pure lcMin. */
	private static final float DLC_BLEND_HI=0.952f;
	/** Info-power exponent for DLC tier weighting. */
	private static final float DLC_INFO_POWER=4.5f;
	/** Info-power exponent for HC tier weighting. */
	private static final float HC_INFO_POWER=2.0f;

	/** DLC-SBS blend: below 2×B, pure SBS. */
	private static final float DLCSBS_BLEND_LO=2.0f;
	/** DLC-SBS blend: above 6×B, pure DLC. */
	private static final float DLCSBS_BLEND_HI=6.0f;
	/** Hybrid+2 blend lower bound (×B). */
	private static final float HYBRID2_BLEND_LO=1.0f;
	/** Hybrid+2 blend upper bound (×B). */
	private static final float HYBRID2_BLEND_HI=6.0f;
	/** LDLC HC blend lower bound (×B). */
	private static final float LDLC_B_LO=0.5f;
	/** LDLC HC blend upper bound (×B). */
	private static final float LDLC_B_HI=4.5f;
	/** Maximum HC weight in LDLC blend. */
	private static final float LDLC_HC_WEIGHT=0.50f;
	/** LDLC weight in HLDLC blend (remainder goes to Hybrid+2). */
	private static final float HLDLC_WEIGHT=0.5f;

	/** LDLC weight in VLDLC blend. */
	private static final double VLDLC_LDLC=0.64;
	/** VWMean weight in VLDLC blend. */
	private static final double VLDLC_VW=0.36;

	/** VW blend: below 1×B, pure SBS for VWMean. */
	private static final float VW_BLEND_LO=1.0f;
	/** VW blend: above 6×B, pure VWMean. */
	private static final float VW_BLEND_HI=6.0f;

	/** VWMean CF: 3S2G formula with coefficients fitted for AVDLL64.
	 *  R²=0.9999997, maxErr=0.00017, terminal=1.36014.
	 *  Independent variable: raw VWMean estimate (before CF).
	 *  @author Nowi */
	/** VWMean CF, x=log2(card/B). Bucket-independent.
	 *  Fitted from resource CF table (B=1408). R²=0.9999983, terminal=1.00004. */
	private static final double[] VWCF={
		 0.523022169365246, -1.580388109567368, -9.550748054346389,  0.165245313710792,
		-0.753249924328123,-10.057664860210931,  0.461936297367938,
		 2.810654431769284,-10.281395080571443,  0.116736516458500,
		-0.000898922132070,  0.603776533689649,  3.876667688716435,
		 0.000937826217617,  0.799259479034718,  3.338171497444599};

	/*--------------------------------------------------------------*/
	/*----------------     Estimation Pipeline      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the HLDLC cardinality estimate (primary estimator).
	 * HLDLC = 0.5 × LDLC + 0.5 × Hybrid+2.
	 * @return estimated number of distinct elements added
	 */
	@Override
	public long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		lastSummarized=summarize();                      // populate for consumeLastSummarized()
		long card=Math.max(0, Math.round(estimate()[1]));
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	/**
	 * Returns the HLDLC cardinality estimate (secondary estimator).
	 * HLDLC = 0.5 × LDLC + 0.5 × Hybrid+2.
	 * @return estimated number of distinct elements added
	 */
	public long cardinalityHLDLC(){
		return Math.max(0, Math.round(estimate()[1]));
	}

	/** Fraction of non-empty registers. */
	public double occupancy(){
		return (double)filledBuckets/modBuckets;
	}

	/*--------------------------------------------------------------*/
	/*----------------   CardStats Integration      ----------------*/
	/*--------------------------------------------------------------*/

	/** Summarizes registers into a CardStats for calibration pipeline compatibility.
	 *  Uses the same register→NLZ mapping as AVDLL64. */
	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[modBuckets];
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){
				if(globalNLZ>=0 && globalNLZ<64){
					nlzCounts[globalNLZ+1]++;
					filledCount++;
				}
			}else{
				final int nlzPart=stateToExp(reg);
				final int histPattern=stateToHist(reg);
				if(nlzPart>=HISTORY_MARGIN){
					final int absNlz=nlzPart-HISTORY_MARGIN+(globalNLZ+1);
					if(absNlz<64){
						nlzCounts[absNlz+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((absNlz+1)<<HIST_BITS)|histPattern);
					}
				}else{
					if(globalNLZ>=0 && globalNLZ<64){
						nlzCounts[globalNLZ+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((globalNLZ+1)<<HIST_BITS)|histPattern);
					}
				}
			}
		}
		nlzCounts[0]=modBuckets-filledCount;
		return new CardStats(packedBuckets, nlzCounts, 0, HIST_BITS, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
	}

	@Override
	public double[] rawEstimates(){
		lastSummarized=summarize();
		final double hybridEst=lastSummarized.hybridDLL();
		final double[] r=AbstractCardStats.buildLegacyArray(lastSummarized, hybridEst);
		r[AbstractCardStats.HC_IDX+1]=lastSummarized.vwMean();
		return r;
	}

	/** Returns and clears the last summarized CardStats (for DDLCalibrationDriver2). */
	public CardStats consumeLastSummarized(){
		final CardStats s=lastSummarized;
		lastSummarized=null;
		return s;
	}

	/** CF table file (shared with AVDLL64 — same estimator family).
	 *  Loaded lazily via CardinalityParser (cf=t), not at class init. */
	public static final String CF_FILE="?cardinalityCorrectionAVDLL64.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX;

	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	private CardStats lastSummarized;

	/*--------------------------------------------------------------*/
	/*----------------      Internal Estimation      ----------------*/
	/*--------------------------------------------------------------*/

	/** Computes both estimators.  Returns {VLDLC, HLDLC}.
	 *  Package-private so DDLCalibrationDriver2 can read the calibrated
	 *  internal estimates (the CardStats pipeline lacks AVLL's VWMean CF). */
	double[] estimate(){
		/* Step 1: Summarize registers into NLZ histogram and packed buckets */
		final int[] counts=new int[66];                  // counts[k+1] = # regs with absNlz=k
		final char[] packedBuckets=new char[modBuckets];
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){                                  // empty register
				if(globalNLZ>=0 && globalNLZ<64){
					counts[globalNLZ+1]++;               // contributes to floor tier
					filledCount++;
				}
			}else{
				final int nlzPart=stateToExp(reg);
				final int histPattern=stateToHist(reg);
				if(nlzPart>=HISTORY_MARGIN){              // above floor zone
					final int absNlz=nlzPart-HISTORY_MARGIN+(globalNLZ+1);
					if(absNlz<64){
						counts[absNlz+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((absNlz+1)<<HIST_BITS)|histPattern);
					}
				}else{                                    // at or below floor zone
					if(globalNLZ>=0 && globalNLZ<64){
						counts[globalNLZ+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((globalNLZ+1)<<HIST_BITS)|histPattern);
					}
				}
			}
		}
		counts[0]=modBuckets-filledCount;                // empty bucket count

		/* Step 2: Compute fundamental sums from NLZ histogram */
		double difSum=0, hllSumFilled=0;
		int filled=0;
		for(int ci=1; ci<counts.length; ci++){
			final int n=counts[ci];
			if(n==0){continue;}
			filled+=n;
			final int absNlz=ci-1;
			final double dif=(absNlz==0 ? (double)Long.MAX_VALUE : Math.pow(2.0, 63-absNlz));
			difSum+=n*dif;                               // sum of dif values for Mean
			hllSumFilled+=n*Math.pow(2.0, -absNlz);     // sum for HLL/HMean
		}
		final int V=modBuckets-filled;                   // empty buckets
		final double alpha_m=0.7213/(1.0+1.079/modBuckets); // HLL bias correction
		final float correction=(filled+modBuckets)/(float)(modBuckets+modBuckets); // empty-bias correction

		/* Step 3: LC — basic linear counting */
		final double lcRaw=(double)modBuckets*Math.log((double)modBuckets/Math.max(V, 0.5));

		/* Step 4: SBS — state-based sampling from per-bucket history states */
		final int histMask=(1<<HIST_BITS)-1;             // =3 for 2-bit history
		final byte[] sbsStates=new byte[modBuckets];
		int sbsContributing=0;
		for(int i=0; i<modBuckets; i++){
			final int val=packedBuckets[i];
			if(val==0){sbsStates[i]=-1; continue;}      // empty bucket
			final int absNlz=(val>>>HIST_BITS)-1;
			final int histPattern=val&histMask;
			final int si=sbsStateIndex(absNlz, histPattern);
			sbsStates[i]=(byte)si;
			if(si>=0){sbsContributing++;}
		}
		double sbsEst=0;
		if(sbsContributing>0){
			for(int i=0; i<modBuckets; i++){
				final int si=sbsStates[i];
				if(si>=0 && si<SBS_STATES){
					sbsEst+=sbsFormula(si, sbsContributing, modBuckets);
				}
			}
		}else{
			sbsEst=lcRaw;                               // fallback to plain LC
		}
		sbsEst=Math.max(sbsEst, filled);                // floor at filled count

		/* Step 5: DLC — dynamic linear counting (info-power weighted tier blend) */
		final double dlcPure=dlcPure(counts, V, modBuckets);
		final double lcMin=lcMin(counts, V, modBuckets);
		final double dlcRaw=dlcBlendWithLcMin(dlcPure, lcMin, V, modBuckets);

		/* Step 6: Mean estimate + CF correction */
		final int div=Math.max(filled, 1);
		final double meanVal=difSum/div;                 // average dif
		final double meanRaw=2*(Long.MAX_VALUE/Math.max(1.0, meanVal))*div*correction*TERMINAL_MEAN_CF;
		final double meanCF=meanRaw*iterativeCF(dlcRaw, meanRaw, MCF_MEAN);

		/* Step 7: History-corrected Mean estimate + MeanH CF */
		double corrDifSum=0;
		if(filled>0){
			final double invTermCF=1.0/TERMINAL_MEANH_CF;
			for(int i=0; i<modBuckets; i++){
				final int val=packedBuckets[i];
				if(val==0){continue;}
				final int absNlz=(val>>>HIST_BITS)-1;
				final int histPattern=val&histMask;
				final int nlzBin=Math.min(absNlz, HIST_BITS+1);
				final double dif=(absNlz==0 ? (double)Long.MAX_VALUE : Math.pow(2.0, 63-absNlz));
				final double hcf=historyOffset(nlzBin, histPattern);
				final double tierMult=Math.pow(2.0, -(hcf+CF_OFFSET))*invTermCF;
				corrDifSum+=dif*tierMult;
			}
		}
		final double corrMeanVal=corrDifSum/div;
		final double meanCorrRaw=2*(Long.MAX_VALUE/Math.max(1.0, corrMeanVal))*div*correction;
		final double meanCorrCF=meanCorrRaw*iterativeCF(dlcRaw, meanCorrRaw, MCF_MEANH);

		/* Step 8: dlcSbs — DLC blended with SBS */
		final double dlcSbs=dlcBlendWithSbs(dlcPure, sbsEst, dlcRaw, modBuckets);

		/* Step 9: HC — history counting (per-tier exact LC with info-power weighting) */
		double hcF=0;
		{
			final int[] nlzBucketCount=new int[64];      // per-NLZ bucket counts
			final int[][] nlzHbitSet=new int[HIST_BITS][64]; // per-bit history set counts
			for(int i=0; i<modBuckets; i++){
				final int val=packedBuckets[i];
				if(val==0){continue;}
				final int absNlz=(val>>>HIST_BITS)-1;
				final int histPattern=val&histMask;
				if(absNlz>=0 && absNlz<64){
					nlzBucketCount[absNlz]++;
					for(int d=0; d<HIST_BITS; d++){
						if((histPattern&(1<<(HIST_BITS-1-d)))!=0){
							nlzHbitSet[d][absNlz]++;     // count of set bits at depth d
						}
					}
				}
			}
			int maxNlzHC=0;
			for(int t=63; t>=0; t--){
				if(nlzBucketCount[t]>0){maxNlzHC=t; break;}
			}
			final int maxTierHC=Math.min(maxNlzHC+2, 63);
			double hcSumW=0, hcSumWLogE=0;
			for(int t=0; t<=maxTierHC; t++){
				int hcBeff=0, hcUnseen=0;                // effective buckets, unseen at this tier
				for(int d=0; d<HIST_BITS; d++){
					final int sourceTier=t+d+1;
					if(sourceTier<64){
						hcBeff+=nlzBucketCount[sourceTier];
						hcUnseen+=(nlzBucketCount[sourceTier]-nlzHbitSet[d][sourceTier]);
					}
				}
				if(hcBeff>=8 && hcUnseen>=1 && hcUnseen<hcBeff){ // enough data for estimation
					final double est=Math.pow(2.0, t+1)/(2.0-1.0) // tierRatio=2 when TIER_SCALE=1
						*(double)modBuckets*Math.log((double)hcBeff/hcUnseen);
					if(est>0 && !Double.isNaN(est)){
						final int hcOcc=hcBeff-hcUnseen;
						final double hcErr=SQRT_2_OVER_PI
							*Math.sqrt((double)hcOcc/((double)hcBeff*hcUnseen))
							/Math.log((double)hcBeff/hcUnseen);
						final double w=Math.pow(1.0/hcErr, HC_INFO_POWER);
						hcSumW+=w;
						hcSumWLogE+=w*Math.log(est);     // log-space weighted sum
					}
				}
			}
			final double hcRaw=(hcSumW>0 ? Math.exp(hcSumWLogE/hcSumW) : 0);
			hcF=(hcRaw>0 ? hcRaw*hcCfFormula(hcRaw) : 0); // HC_SCALE=1.0 for AVDLL64
		}

		/* Step 10: LDLC — DlcSbs blended with HC across [0.5B, 4.5B] */
		double ldlc;
		{
			final double B=(double)modBuckets;
			final double bLo=LDLC_B_LO*B, bHi=LDLC_B_HI*B;
			final boolean hcUsable=(hcF>0 && dlcSbs>0);
			if(dlcSbs<=bLo || !hcUsable){
				ldlc=dlcSbs;                             // low card: pure DLC-SBS
			}else if(dlcSbs<=bHi){
				final double t=(dlcSbs-bLo)/(bHi-bLo);  // linear ramp
				final double hcW=t*LDLC_HC_WEIGHT;       // HC weight grows from 0 to max
				ldlc=(1-hcW)*dlcSbs+hcW*hcF;
			}else{
				ldlc=(1-LDLC_HC_WEIGHT)*dlcSbs+LDLC_HC_WEIGHT*hcF; // max HC weight
			}
		}

		/* Step 11: Hybrid+2 — SBS → Mean+H blend using DLC as zone detector */
		double hybridPlus2;
		{
			final double hb0=HYBRID2_BLEND_LO*modBuckets; // =1×B
			final double hb1=HYBRID2_BLEND_HI*modBuckets; // =6×B
			if(dlcRaw<=hb0){
				hybridPlus2=sbsEst;                      // low card: pure SBS
			}else if(dlcRaw<hb1){
				final double t=Math.log(dlcRaw/hb0)/Math.log(hb1/hb0); // log interpolation
				hybridPlus2=(1-t)*sbsEst+t*meanCorrCF;   // SBS → Mean+H
			}else{
				hybridPlus2=meanCorrCF;                  // high card: pure Mean+H
			}
		}

		/* Step 12: VWMean — per-tier-state weighted harmonic mean */
		double vwMean=0;
		{
			double sumInvV=0;
			int used=0;
			for(int i=0; i<modBuckets; i++){
				final int val=packedBuckets[i];
				if(val==0){continue;}
				final int absNlz=(val>>>HIST_BITS)-1;
				final int histPattern=val&histMask;
				if(absNlz<0){continue;}
				final int tier=Math.min(absNlz, VW_MAX_TIER);
				if(tier>=VW_TABLE.length || histPattern>=VW_TABLE[tier].length){continue;}
				double v=VW_TABLE[tier][histPattern];    // expected distinct count
				if(v<=0){continue;}
				if(absNlz>VW_MAX_TIER){                  // geometric extrapolation beyond table
					v*=Math.pow(2.0, absNlz-VW_MAX_TIER);
				}
				used++;
				sumInvV+=1.0/v;                          // harmonic mean accumulator (w=1.0)
			}
			if(used>0 && sumInvV>0){
				vwMean=(double)used*used/sumInvV;         // used × harmonic_mean
				vwMean*=iterativeCF(dlcRaw, vwMean, VWCF); // CF correction keyed on dlcRaw/B
			}
		}

		/* Step 13: VW blend — SBS fallback at low cardinality */
		double vw=vwMean;
		if(vwMean>0){
			final double hb0=VW_BLEND_LO*modBuckets;    // =1×B
			final double hb1=VW_BLEND_HI*modBuckets;    // =6×B
			if(dlcRaw<=hb0){
				vw=sbsEst;                               // low card: pure SBS
			}else if(dlcRaw<hb1){
				final double t=Math.log(dlcRaw/hb0)/Math.log(hb1/hb0);
				vw=(1-t)*sbsEst+t*vwMean;               // blend SBS → VWMean
			}
		}

		/* Step 14: Compose final estimators */
		final double hldlc=HLDLC_WEIGHT*ldlc+(1-HLDLC_WEIGHT)*hybridPlus2;
		final double vldlc=VLDLC_LDLC*ldlc+VLDLC_VW*vw;

		return new double[]{vldlc, hldlc, meanCF, meanCorrCF, hybridPlus2, ldlc, vw, dlcRaw, hcF, sbsEst};
	}

	/*--------------------------------------------------------------*/
	/*----------------     DLC Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Tier-compensated LC: walks the NLZ histogram upward from V (empty buckets).
	 *  V_k = V + sum(counts[1..k]).  First k where V_k > 0 gives lcMin. */
	private static double lcMin(int[] counts, int V, int B){
		if(V>0){return (double)B*Math.log((double)B/Math.max(V, 0.5));} // tier 0 has empties
		int vk=V;                                        // starts at 0 when all filled
		for(int k=0; k<counts.length-1; k++){
			vk+=counts[k+1];                             // accumulate tier k buckets
			if(vk>0){                                     // first tier with "empties"
				return Math.pow(2.0, k+1)*(double)B*Math.log((double)B/Math.max(vk, 0.5));
			}
		}
		return 0;                                        // all tiers full (shouldn't happen)
	}

	/** DLC estimate at a single tier: 2^tier × B × ln(B / max(V_k, 0.5)). */
	private static double dlcFromVk(int tier, int vk, int B){
		return Math.pow(2.0, tier)*(double)B*Math.log((double)B/Math.max(vk, 0.5));
	}

	/** Info-power weighted DLC (all tiers blended in log-space).
	 *  Uses theoretical error model: err = √(2/π) × √(occ/(B×V)) / ln(B/V). */
	private static double dlcPure(int[] counts, int V, int B){
		final int minVK=Math.max(Math.max(1, 2), (int)(B*0.002)); // min V_k for participation
		final int maxVK=B-minVK;
		final double Bd=(double)B;
		int startTier=0;
		for(int k=0; k<counts.length-1; k++){if(counts[k+1]>0){startTier=k; break;}}
		int vk=V;
		double sumW=0, sumWLogE=0;
		for(int tier=startTier; tier<64; tier++){
			if(vk>=minVK && vk<=maxVK){
				final double est=dlcFromVk(tier, vk, B);
				final int occ=B-vk;
				final double err=SQRT_2_OVER_PI*Math.sqrt((double)occ/(Bd*vk))/Math.log(Bd/vk);
				final double w=Math.pow(1.0/err, DLC_INFO_POWER);
				sumW+=w; sumWLogE+=w*Math.log(est);      // log-space weighted sum
			}
			if(tier+1<counts.length && tier<63){vk+=counts[tier+1];}
			if(vk>=B){break;}                            // all buckets accounted for
		}
		if(sumW<=0){return Double.NaN;}
		return Math.exp(sumWLogE/sumW);                  // geometric weighted mean
	}

	/** Blend DLC with lcMin at low occupancy (V/B space). */
	private static double dlcBlendWithLcMin(double dlcPure, double lcMin, int V, int B){
		if(Double.isNaN(dlcPure) || V>=DLC_BLEND_HI*B){return lcMin;}
		if(V>DLC_BLEND_LO*B){
			final double t=(V-DLC_BLEND_LO*B)/((DLC_BLEND_HI-DLC_BLEND_LO)*B);
			return t*lcMin+(1-t)*dlcPure;                // linear blend in V/B space
		}
		return dlcPure;
	}

	/** Blend DLC with SBS in cardinality space (log interpolation). */
	private static double dlcBlendWithSbs(double dlcPure, double sbs,
			double dlcRaw, int B){
		if(Double.isNaN(dlcPure) || dlcRaw<=DLCSBS_BLEND_LO*B){return sbs;}
		if(dlcRaw>=DLCSBS_BLEND_HI*B){return dlcPure;}
		final double t=Math.log(dlcRaw/(DLCSBS_BLEND_LO*B))
			/Math.log((double)DLCSBS_BLEND_HI/DLCSBS_BLEND_LO);
		return (1-t)*sbs+t*dlcPure;
	}

	/*--------------------------------------------------------------*/
	/*----------------   History Offset Table         ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-history-state correction offset for Mean+H computation.
	 *  historyOffset(nlzBin, histPattern) returns the log2 adjustment
	 *  that accounts for the information content of the history pattern.
	 *  Geometric mean averaging (all/geo model).
	 *  Source: StateTable.CF_HISTORY_2 and CF_HISTORY_2_TIERS.
	 *  Pre-trained values embedded to avoid external StateTable dependency. */
	private static final double CF_OFFSET=0.0;

	/** Per-tier history offset tables.  Tier 0-2 are per-tier; tier 3+ uses steady state.
	 *  Regenerated 2026-04-12 from mc2_hist2 GeoCF rows. */
	private static final double[][] HIST_TIERS={
		{+0.00000000, +0.00000000, +0.00000000, +0.00000000},  // tier 0: no history
		{-1.64516934, +0.00000000, +0.32902894, +0.00000000},  // tier 1: only MSB valid
		{-2.91470444, -1.11691781, -1.63200820, +0.37491432},  // tier 2: both bits valid
	};

	/** Steady-state per-state CFs for nlzBin >= 3 (2-bit history, 4 states).
	 *  Regenerated 2026-04-12 from mc2_hist2_128m_32k.txt tier 11 GeoCF. */
	private static final double[] HIST_STEADY={-2.45693840, -1.03241245, -1.48290584, +0.33548548};

	private static double historyOffset(int nlzBin, int histPattern){
		if(nlzBin<HIST_TIERS.length){
			return HIST_TIERS[nlzBin][histPattern];       // per-tier correction
		}
		return HIST_STEADY[histPattern];                  // steady-state for deep tiers
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Arithmetic-encoded register array. */
	private final long[] packed;
	/** Number of 64-bit words in packed array. */
	private final int packedLen;
	/** Actual register count (multiple of 11). */
	private final int modBuckets;
	/** Global NLZ floor (advances as registers fill). */
	private int globalNLZ=-1;
	/** Registers at or below HISTORY_MARGIN. */
	private int floorCount;
	/** Non-empty register count. */
	private int filledBuckets=0;

	long reads=0, writes=0;
	public long registerReads(){return reads;}
	public long registerWrites(){return writes;}

}
