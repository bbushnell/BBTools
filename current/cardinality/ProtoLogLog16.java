package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * ProtoLogLog16: 16-bit register prototype for testing sub-NLZ information.
 * <p>
 * NOT for production use. Designed for paper experiments comparing different
 * 2-bit extraction methods (mantissa, andtissa, nlz2, history, luck) and
 * their combinations, with configurable NLZ bit width.
 * <p>
 * Register layout (16 bits per bucket, stored as short[]):
 * <pre>
 *   Bits 15..10: nlzStored = absNlz+1 (6 bits, 0=empty, max NLZ=62)
 *   Bits 9..0:   extra info (10 bits available for any combination)
 * </pre>
 * The extra bits are partitioned by mode:
 * <ul>
 *   <li>MANTISSA: inverted top N bits after leading 1-bit (2-8 bits)</li>
 *   <li>NLZ2: min(2^N-1, NLZ(key << (NLZ(key)+1))) (2-8 bits)</li>
 *   <li>HISTORY: shift-register sawMinus1..sawMinusN (2-8 bits)</li>
 *   <li>LUCK: min(2^N-1, maxNlz - secondMaxNlz - 1) (2-8 bits)</li>
 *   <li>Combined: e.g. 2-bit mantissa + 2-bit history = 4-bit extra</li>
 * </ul>
 * <p>
 * The NLZ bit width is configurable (3-6, default 6) for testing DLL3-DLL6
 * equivalents. Extra bits = 16 - nlzBits.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class ProtoLogLog16 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------          Modes              ----------------*/
	/*--------------------------------------------------------------*/

	/** Sub-NLZ extraction modes. Can be combined (e.g., MANTISSA + HISTORY). */
	public static final int MODE_NONE=0;
	public static final int MODE_MANTISSA=1;   // inverted top bits after leading 1
	public static final int MODE_NLZ2=2;       // NLZ of remainder after first NLZ+1 bits
	public static final int MODE_HISTORY=4;    // shift-register sawMinus1..sawMinusN
	public static final int MODE_LUCK=8;       // min(cap, maxNlz - secondMaxNlz - 1)
	public static final int MODE_ANDTISSA=16;  // mantissa of (hash & (hash<<2))

	/** Active mode bitmask. Set via command line, e.g., mode=mantissa+history */
	public static int MODE=MODE_MANTISSA;

	/** Bits allocated to each active mode. Must sum to <= (16 - nlzBits). */
	public static int MANTISSA_BITS=2;
	public static int NLZ2_BITS=2;
	public static int HISTORY_BITS=2;
	public static int LUCK_BITS=2;
	public static int ANDTISSA_BITS=2;

	/** NLZ bits (3-6). Default 6 = LL6 equivalent. */
	public static int NLZ_BITS=6;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ProtoLogLog16(){this(2048, 31, -1, 0);}

	ProtoLogLog16(Parser p){
		super(p);
		maxArray=new short[buckets];
		luckSecond=usesLuck() ? new byte[buckets] : null;
	}

	ProtoLogLog16(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new short[buckets];
		luckSecond=usesLuck() ? new byte[buckets] : null;
	}

	@Override
	public ProtoLogLog16 copy(){return new ProtoLogLog16(buckets, k, -1, minProb);}

	private static boolean usesLuck(){return (MODE & MODE_LUCK)!=0;}
	private static boolean usesHistory(){return (MODE & MODE_HISTORY)!=0;}

	/*--------------------------------------------------------------*/
	/*----------------         Bit Layout          ----------------*/
	/*--------------------------------------------------------------*/

	/** Compute bit layout from current mode settings. */
	private static int nlzShift(){return 16-NLZ_BITS;}
	private static int nlzMask(){return (1<<NLZ_BITS)-1;}
	private static int maxNlzStored(){return (1<<NLZ_BITS)-1;} // 0=empty

	/**
	 * Extra bits layout: each enabled mode gets its bits packed from MSB down.
	 * Returns the shift and mask for each mode's bits within the extra field.
	 */
	private static int extraBits(){return 16-NLZ_BITS;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Extract absNlz from stored value. */
	private static int getAbsNlz(int stored){
		return (stored>>>nlzShift())-1;
	}

	/** Extract the full extra field (all non-NLZ bits). */
	private static int getExtra(int stored){
		return stored & ((1<<nlzShift())-1);
	}

	/**
	 * Compute the extra bits for a given hash, primary NLZ, and current bucket state.
	 * Each mode contributes its bits; they are packed together.
	 */
	private int computeExtra(long key, int nlz, int bucket, int oldStored){
		int extra=0;
		int shift=0; // current bit position (from LSB)

		if((MODE & MODE_MANTISSA)!=0){
			final int mbits=MANTISSA_BITS;
			final int mshift=63-nlz-mbits; // skip leading zeros + leading 1
			final int val=(mshift<0) ? 0 : (int)((~(key>>>mshift))&((1<<mbits)-1));
			extra|=(val<<shift);
			shift+=mbits;
		}

		if((MODE & MODE_ANDTISSA)!=0){
			final int abits=ANDTISSA_BITS;
			final long anded=key&(key<<2);
			final int ashift=63-nlz-abits;
			final int val=(ashift<0) ? 0 : (int)((anded>>>ashift)&((1<<abits)-1));
			extra|=(val<<shift);
			shift+=abits;
		}

		if((MODE & MODE_NLZ2)!=0){
			final int nbits=NLZ2_BITS;
			final int cap=(1<<nbits)-1;
			final int consumed=nlz+1;
			int val;
			if(consumed>=64){val=cap;}
			else{val=Math.min(cap, Long.numberOfLeadingZeros(key<<consumed));}
			extra|=(val<<shift);
			shift+=nbits;
		}

		if((MODE & MODE_HISTORY)!=0){
			final int hbits=HISTORY_BITS;
			// History carry from old state
			final int oldNlzStored=(oldStored>0) ? (oldStored>>>nlzShift()) : 0;
			final int newNlzStored=Math.min(nlz+1, maxNlzStored());
			int hist=0;
			if(oldStored>0){
				if(newNlzStored>oldNlzStored){
					// Promotion: shift register carry
					final int k=newNlzStored-oldNlzStored;
					final int oldHist=(oldStored>>shift)&((1<<hbits)-1);
					final int carry=(1<<hbits); // "old max existed" bit
					hist=((oldHist|carry)>>k)&((1<<hbits)-1);
				}else{
					// Tie or lower: preserve existing history
					hist=(oldStored>>shift)&((1<<hbits)-1);
				}
			}
			extra|=(hist<<shift);
			shift+=hbits;
		}

		if((MODE & MODE_LUCK)!=0){
			final int lbits=LUCK_BITS;
			final int cap=(1<<lbits)-1;
			final int oldMax=(oldStored>0) ? getAbsNlz(oldStored) : -1;
			final int effectiveMax=Math.max(nlz, oldMax);
			final int sec=(luckSecond[bucket]&0xFF)-1;
			int luckVal;
			if(effectiveMax>=0 && sec>=0){
				// Gap 1 (e.g. 5&4) -> State 0 [Established]
				// Gap cap+1 -> State cap [Outlier]
				luckVal=Math.min(cap, Math.max(0, effectiveMax-sec-1));
			}else{
				luckVal=cap; // No second best = Outlier
			}
			extra|=(luckVal<<shift);
			shift+=lbits;
		}

		return extra;
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

		final int newNlzStored=Math.min(nlz+1, maxNlzStored());
		final int oldStored=maxArray[bucket]&0xFFFF;
		final int oldNlzStored=(oldStored>0) ? (oldStored>>>nlzShift()) : 0;

		// Update luck leaderboard BEFORE any register changes.
		// Tracks top-2 all-time NLZ values. Uses +1 encoding (0=unset).
		if(usesLuck()){
			final int curMax=(oldStored>0) ? getAbsNlz(oldStored) : -1;
			final int sec=(luckSecond[bucket]&0xFF)-1;
			if(nlz>curMax){
				// Promotion: old max becomes second
				luckSecond[bucket]=(byte)(curMax+1);
			}else if(nlz>sec){
				// Silver medalist: only update if strictly better
				luckSecond[bucket]=(byte)(nlz+1);
			}
		}

		if(newNlzStored>oldNlzStored){
			// Promotion: compute new extra bits
			lastCardinality=-1;
			if(oldStored==0){filledBuckets++;}
			final int extra=computeExtra(key, nlz, bucket, oldStored);
			maxArray[bucket]=(short)((newNlzStored<<nlzShift())|extra);
		}else if(newNlzStored==oldNlzStored){
			// Same NLZ: update extra bits.
			// Mantissa/nlz2/andtissa: max-update (higher = rarer = better).
			// Luck: min-update (smaller gap = more established = better).
			final int newExtra=computeExtra(key, nlz, bucket, oldStored);
			final int oldExtra=getExtra(oldStored);
			final boolean improved=usesLuck() ? (newExtra<oldExtra) : (newExtra>oldExtra);
			if(improved){
				lastCardinality=-1;
				maxArray[bucket]=(short)((newNlzStored<<nlzShift())|newExtra);
			}
			// History: update bits for sub-max NLZ elements (active mode)
			// Handled via computeExtra returning updated history
		}else{
			// Lower NLZ: update history bits if active
			if(usesHistory()){
				final int diff=oldNlzStored-newNlzStored;
				final int hshift=getHistoryShift();
				final int hbits=HISTORY_BITS;
				if(diff>=1 && diff<=hbits){
					// Set the corresponding history bit
					final int bit=1<<(hbits-diff); // bit for NLZ-(diff)
					final int newStored=oldStored|(bit<<hshift);
					if(newStored!=oldStored){
						maxArray[bucket]=(short)newStored;
						lastCardinality=-1;
					}
				}
			}
		}
	}

	/** Get the bit shift for history within the extra field. */
	private static int getHistoryShift(){
		int shift=0;
		if((MODE & MODE_MANTISSA)!=0){shift+=MANTISSA_BITS;}
		if((MODE & MODE_ANDTISSA)!=0){shift+=ANDTISSA_BITS;}
		if((MODE & MODE_NLZ2)!=0){shift+=NLZ2_BITS;}
		return shift;
	}

	private CardinalityStats summarize(){
		// Build runtime multiplier table from CFs + offset
		final int totalExtraBits=getActiveExtraBits();
		final int numStates=1<<totalExtraBits;
		final double[] mult;
		if(CORRECTION_TABLE!=null && CORRECTION_TABLE.length==numStates){
			mult=new double[numStates];
			for(int s=0; s<numStates; s++){
				// History/mantissa/etc: outlier has negative CF, needs 2^(-CF) = mult > 1
			//   (inflate sum → reduce estimate for lucky outlier NLZ)
			// Luck: outlier has negative CF, needs 2^(CF) = mult < 1
			//   (reduce contribution → de-emphasize lucky outlier)
			final double cf=CORRECTION_TABLE[s]+CF_OFFSET;
			mult[s]=Math.pow(2.0, usesLuck() ? cf : -cf);
			}
		}else{
			mult=null; // no correction
		}

		final int[] nlzCounts=new int[64];
		double difSum=0, hllSumFilled=0, gSum=0;
		int count=0;
		final int extraMask=(1<<nlzShift())-1;
		for(int i=0; i<buckets; i++){
			final int stored=maxArray[i]&0xFFFF;
			if(stored>0){
				final int absNlz=getAbsNlz(stored);
				final int extra=stored&extraMask;
				if(absNlz>=0 && absNlz<64){nlzCounts[absNlz]++;}
				// Infancy guard: don't apply luck corrections until bucket has 2+ values.
			// A bucket with only one hash has luckSecond==0 and its "gap" is undefined.
			final boolean canCorrect=!usesLuck() || (luckSecond!=null && (luckSecond[i]&0xFF)>0);
			final double m=(mult!=null && absNlz>=3 && extra<mult.length && canCorrect) ? mult[extra] : 1.0;
				final double base=Math.pow(2.0, -absNlz)*m;
				final double dif=((absNlz==0 ? (double)Long.MAX_VALUE : (absNlz<64 ? (double)(1L<<(63-absNlz)) : 1.0)))*m;
				difSum+=dif;
				hllSumFilled+=base;
				gSum+=Math.log(Tools.max(1, dif));
				count++;
			}
		}
		lastRawNlz=nlzCounts;
		return new CardinalityStats(difSum, hllSumFilled, hllSumFilled,
			gSum, count, buckets, null, CF_MATRIX, CF_BUCKETS,
			CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys, microIndex,
			nlzCounts, 0);
	}

	/** Returns total extra bits used by active modes. */
	private static int getActiveExtraBits(){
		int bits=0;
		if((MODE & MODE_MANTISSA)!=0) bits+=MANTISSA_BITS;
		if((MODE & MODE_ANDTISSA)!=0) bits+=ANDTISSA_BITS;
		if((MODE & MODE_NLZ2)!=0) bits+=NLZ2_BITS;
		if((MODE & MODE_HISTORY)!=0) bits+=HISTORY_BITS;
		if((MODE & MODE_LUCK)!=0) bits+=LUCK_BITS;
		return bits;
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
		throw new UnsupportedOperationException("ProtoLogLog16 merge not implemented");
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

	final short[] maxArray;
	/** Per-bucket second-highest NLZ, for luck mode. */
	private final byte[] luckSecond;
	private int filledBuckets=0;
	int[] lastRawNlz;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Per-state correction table (additive CFs). Set externally based on mode.
	 * Length must match 2^(active extra bits). null = no correction.
	 * Larger CF values INCREASE the effective NLZ, raising the cardinality estimate.
	 * Set via setCorrectionTable() or command line.
	 */
	public static double[] CORRECTION_TABLE=null;

	/** Additive offset to all correction table entries.
	 *  Optimized via CV sweep for each mode. */
	public static double CF_OFFSET=0.0;

	/**
	 * Preset empirical correction tables for 2-bit modes (tier 8 steady-state).
	 * Call setMode() to select both the mode and its correction table.
	 */
	// 1-bit empirical CFs (tier 8 steady-state, single-bucket simulation)
	public static final double[] CF_MANTISSA_1={-0.3427, +0.1495};
	public static final double[] CF_ANDTISSA_1={-0.2036, +0.2647};
	public static final double[] CF_NLZ2_1={-0.3427, +0.1495};
	public static final double[] CF_HISTORY_1={-1.3686, +0.1694};
	// Luck: raw gap. State 0=established(small gap), cap=outlier(big gap). Min-update.
	public static final double[] CF_LUCK_1={+0.4318, -0.5923};

	// 2-bit empirical CFs (tier 8 steady-state, single-bucket simulation)
	public static final double[] CF_MANTISSA_2={-0.4928, -0.2633, -0.0342, +0.2640};
	public static final double[] CF_ANDTISSA_2={-0.3295, +0.0017, +0.2241, +0.3588};
	public static final double[] CF_NLZ2_2={-0.3562, -0.0342, +0.1648, +0.3382};
	public static final double[] CF_HISTORY_2={-2.5208, -1.1662, -1.6570, +0.2078};
	public static final double[] CF_LUCK_2={+0.4318, -0.3681, -1.3730, -2.6120};

	// 3-bit empirical CFs
	public static final double[] CF_MANTISSA_3={-0.5254, -0.4335, -0.3224, -0.1984, -0.0785, +0.0176, +0.1590, +0.3308};
	public static final double[] CF_ANDTISSA_3={-0.3769, -0.1509, -0.0180, +0.0894, +0.1511, +0.2919, +0.2871, +0.4103};
	public static final double[] CF_NLZ2_3={-0.3427, -0.0258, +0.1590, +0.2718, +0.3745, +0.3318, +0.4442, +0.4378};
	public static final double[] CF_HISTORY_3={-3.5958, -2.2916, -2.7710, -1.1179, -2.8612, -1.5939, -2.2068, +0.2139};
	public static final double[] CF_LUCK_3={+0.4318, -0.3681, -1.3730, -2.3818, -3.4227, -4.4743, -5.5593, -6.7761};

	// 4-bit empirical CFs
	public static final double[] CF_HISTORY_4={-4.6751, -3.3747, -3.7614, -2.2344, -3.9271, -2.6884, -3.1876, -1.1089, -3.9923, -2.7859, -3.0713, -1.5813, -3.3114, -2.1801, -2.7362, +0.2145};
	// 4-bit luck: TODO regenerate with inverted encoding and sufficient iterations
	public static final double[] CF_LUCK_4=null;

	/** Set mode and auto-select the matching empirical correction table
	 *  based on the current XXXX_BITS setting for that mode. */
	public static void setMode(int mode){
		MODE=mode;
		if(mode==MODE_MANTISSA){
			final int b=MANTISSA_BITS;
			CORRECTION_TABLE=(b==1 ? CF_MANTISSA_1 : b==3 ? CF_MANTISSA_3 : CF_MANTISSA_2);
		}else if(mode==MODE_ANDTISSA){
			final int b=ANDTISSA_BITS;
			CORRECTION_TABLE=(b==1 ? CF_ANDTISSA_1 : b==3 ? CF_ANDTISSA_3 : CF_ANDTISSA_2);
		}else if(mode==MODE_NLZ2){
			final int b=NLZ2_BITS;
			CORRECTION_TABLE=(b==1 ? CF_NLZ2_1 : b==3 ? CF_NLZ2_3 : CF_NLZ2_2);
		}else if(mode==MODE_HISTORY){
			final int b=HISTORY_BITS;
			CORRECTION_TABLE=(b==1 ? CF_HISTORY_1 : b==3 ? CF_HISTORY_3 :
			                  b==4 ? CF_HISTORY_4 : CF_HISTORY_2);
		}else if(mode==MODE_LUCK){
			final int b=LUCK_BITS;
			CORRECTION_TABLE=(b==1 ? CF_LUCK_1 : b==3 ? CF_LUCK_3 :
			                  b==4 ? CF_LUCK_4 : CF_LUCK_2);
		}else{CORRECTION_TABLE=null;}
	}

	public static final String CF_FILE="?cardinalityCorrectionLL6.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}
}
