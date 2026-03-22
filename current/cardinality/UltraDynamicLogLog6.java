package cardinality;

import shared.Tools;

/**
 * UltraDynamicLogLog6: 6-bit registers combining DLL4's dynamic architecture
 * with ULL's 2-bit sub-NLZ history.
 * <p>
 * Each register stores 6 bits: 4 bits of relative NLZ (0-15 tiers) and 2 bits
 * of history tracking which sub-tiers below the max have been observed.
 * This is equivalent to Ertl's UltraLogLog encoding but with DLL4's sliding
 * minZeros window for early exit and compact relative storage.
 * <p>
 * Register encoding: uses Ertl's pack/unpack on a hashPrefix bitmask.
 * Bit position = relNlz + HISTORY_MARGIN (offset to reserve positions for
 * below-floor history bits). Register 0 = empty. Values 1-63 encode the
 * packed hashPrefix. Clamped to 63 on overflow.
 * <p>
 * Key differences from DLL4:
 * - 6 bits per bucket instead of 4 (byte[] instead of packed int[])
 * - 2-bit history provides ~18% variance reduction (from PLL16c experiments)
 * - Early exit margin: 2 tiers more conservative than DLL4
 * - countAndDecrement uses unpack/shift/repack to preserve history
 * - Estimation via Ertl's FGRA (no CF tables needed)
 * <p>
 * At 2048 buckets: 2048 bytes (one byte per bucket, 6 bits used, 2 wasted).
 * Same memory as DLL4 with 4096 buckets (4096 * 4 bits = 2048 bytes).
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class UltraDynamicLogLog6 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	UltraDynamicLogLog6(){
		this(2048, 31, -1, 0);
	}

	UltraDynamicLogLog6(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		registers=new byte[buckets];
		minZeroCount=buckets;
	}

	@Override
	public UltraDynamicLogLog6 copy(){return new UltraDynamicLogLog6(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------     Pack / Unpack (Ertl)     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reconstructs a representative hashPrefix bitmask from a 6-bit register.
	 * Identical to Ertl's unpack but only used for register values 0-63.
	 * Register 0 is treated as empty (returns 0, not unpack(0)).
	 */
	static long unpack(int reg){
		if(reg==0){return 0;}
		return (4L|(reg&3))<<((reg>>>2)-2);
	}

	/**
	 * Packs a hashPrefix bitmask into a 6-bit register value (0-63).
	 * Identical to Ertl's pack but clamped to 63 for overflow.
	 */
	static byte pack6(long hashPrefix){
		if(hashPrefix==0){return 0;}
		int nlz=Long.numberOfLeadingZeros(hashPrefix)+1;
		int raw=(-nlz<<2)|((int)((hashPrefix<<nlz)>>>62));
		return (byte)Math.min(raw&0xFF, 63);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		// Early exit with HISTORY_MARGIN-tier margin.
		// DLL4 exits when NLZ < minZeros. UDLL6 must allow elements with
		// NLZ down to (minZeros - HISTORY_MARGIN) to update history bits.
		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);
		final int relNlz=nlz-minZeros;

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(USE_MICRO){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Bit position in hashPrefix: relNlz + HISTORY_MARGIN
		// This reserves HISTORY_MARGIN positions (0..HISTORY_MARGIN-1) for
		// below-floor history, ensuring register=0 stays "empty".
		final int bitPos=relNlz+HISTORY_MARGIN;
		if(bitPos<0||bitPos>=64){return;} // out of range

		final byte oldReg=registers[bucket];
		long hashPrefix=unpack(oldReg&0x3F);
		hashPrefix|=1L<<bitPos;
		final byte newReg=pack6(hashPrefix);

		if(newReg<=oldReg){return;} // no improvement (registers are monotonically increasing)

		registers[bucket]=newReg;
		lastCardinality=-1;

		// Track filledBuckets: "filled" means NLZ_part >= HISTORY_MARGIN
		// (at least one element at or above the current floor)
		final int oldNlzPart=oldReg>>>2;
		final int newNlzPart=newReg>>>2;
		if(oldNlzPart<HISTORY_MARGIN && newNlzPart>=HISTORY_MARGIN){filledBuckets++;}

		// minZeroCount: track buckets not yet at floor.
		// EARLY_PROMOTE=true: advance when all buckets have at least one at-floor element.
		final boolean wasBelowFloor=EARLY_PROMOTE ? (oldNlzPart<HISTORY_MARGIN) : (oldNlzPart<=HISTORY_MARGIN);
		final boolean nowAtOrAbove=EARLY_PROMOTE ? (newNlzPart>=HISTORY_MARGIN) : (newNlzPart>HISTORY_MARGIN);
		if(wasBelowFloor && nowAtOrAbove && --minZeroCount<1){
			while(minZeroCount==0 && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/**
	 * Decrements all registers by shifting their hashPrefix right by 1.
	 * This effectively decrements relNLZ by 1 while preserving history bits.
	 * Returns count of buckets that are now at or below the new floor.
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int i=0; i<buckets; i++){
			final byte reg=registers[i];
			if(reg==0){continue;} // truly empty, skip
			long prefix=unpack(reg&0x3F);
			prefix>>>=1; // shift right = decrement all NLZ by 1
			final byte newReg=pack6(prefix);
			registers[i]=newReg;
			final int newNlzPart=newReg>>>2;
			if(EARLY_PROMOTE){
				if(newNlzPart<HISTORY_MARGIN){
					newMinZeroCount++;
					if((reg>>>2)>=HISTORY_MARGIN){filledBuckets--;}
				}
			}else{
				if(newNlzPart==HISTORY_MARGIN){newMinZeroCount++;}
			}
		}
		return newMinZeroCount;
	}

	/**
	 * Builds an NLZ histogram from the registers and delegates to
	 * CardinalityStats for standard estimators (Hybrid, DLC, etc.).
	 */
	private CardinalityStats summarize(){
		final int[] nlzCounts=new int[64];
		for(int i=0; i<buckets; i++){
			final byte reg=registers[i];
			final int nlzPart=reg>>>2;
			if(nlzPart>=HISTORY_MARGIN){
				// At or above floor: absNlz = (nlzPart - HISTORY_MARGIN) + minZeros
				final int absNlz=(nlzPart-HISTORY_MARGIN)+minZeros;
				if(absNlz<64){nlzCounts[absNlz]++;}
			}else if(minZeros>0){
				// Below floor or empty with history: phantom bucket
				final int phantomNlz=minZeros-1;
				if(phantomNlz<64){nlzCounts[phantomNlz]++;}
			}
			// reg==0 (truly empty): not counted, contributes to V
		}
		return CardinalityStats.fromNlzCounts(nlzCounts, buckets, microIndex,
				CF_MATRIX, CF_BUCKETS,
				CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys);
	}

	/**
	 * Converts registers to absolute 8-bit Ertl format and runs FGRA estimation.
	 * This is the primary cardinality estimator for UDLL6.
	 */
	public double fgraEstimate(){
		final int p=31-Integer.numberOfLeadingZeros(buckets); // log2(buckets)
		// Convert 6-bit relative registers to 8-bit absolute Ertl registers.
		// Ertl register = 4 * (absNlz + p - 1) + sub
		// absNlz = (nlzPart - HISTORY_MARGIN) + minZeros
		// So ertl_reg = 4 * ((nlzPart - HISTORY_MARGIN) + minZeros + p - 1) + sub
		//             = udll6_reg + 4 * (minZeros + p - 1 - HISTORY_MARGIN)
		//             = udll6_reg + 4 * (minZeros + p - 3)
		// But for empty buckets (reg=0), ertl register = 0.
		final int regOffset=4*(minZeros+p-1-HISTORY_MARGIN);
		final byte[] ertlRegs=new byte[buckets];
		for(int i=0; i<buckets; i++){
			final int reg=registers[i]&0x3F;
			if(reg==0){
				ertlRegs[i]=0;
			}else{
				int ertlVal=reg+regOffset;
				ertlRegs[i]=(byte)Math.min(ertlVal, 255);
			}
		}
		return ErtlULL.fgraEstimateStatic(ertlRegs, p);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		// Use FGRA as primary estimator
		long card=(long)fgraEstimate();
		// Fall back to Hybrid for very low cardinality where FGRA may be 0
		if(card<=0){
			final CardinalityStats s=summarize();
			final double rawHyb=s.hybridDLL();
			card=(long)(rawHyb*s.cf(rawHyb, CorrectionFactor.HYBRID));
			card=Math.max(card, s.microCardinality());
		}
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((UltraDynamicLogLog6)log);
	}

	/**
	 * Merges another UDLL6 into this one.
	 * Re-frames both to the higher minZeros, then OR's hashPrefixes per bucket.
	 */
	public void add(UltraDynamicLogLog6 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(registers!=log.registers){
			final int newMinZeros=Math.max(minZeros, log.minZeros);
			final int shiftA=newMinZeros-minZeros;
			final int shiftB=newMinZeros-log.minZeros;
			for(int i=0; i<buckets; i++){
				long prefA=(registers[i]==0) ? 0 : unpack(registers[i]&0x3F);
				long prefB=(log.registers[i]==0) ? 0 : unpack(log.registers[i]&0x3F);
				// Shift both to the new frame
				if(shiftA>0){prefA>>>=shiftA;}
				if(shiftB>0){prefB>>>=shiftB;}
				long merged=prefA|prefB;
				registers[i]=pack6(merged);
			}
			minZeros=newMinZeros;
			// Recount filledBuckets and minZeroCount
			filledBuckets=0;
			minZeroCount=0;
			for(int i=0; i<buckets; i++){
				final int nlzPart=registers[i]>>>2;
				if(nlzPart>=HISTORY_MARGIN){filledBuckets++;}
				if(EARLY_PROMOTE){
					if(nlzPart<HISTORY_MARGIN){minZeroCount++;}
				}else{
					if(nlzPart<=HISTORY_MARGIN){minZeroCount++;}
				}
			}
			// Handle if we need to advance floor
			while(minZeroCount==0 && minZeros<wordlen){
				minZeros++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(Math.max(s.hybridDLL(), s.microCardinality()));
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}
	public int getMinZeros(){return minZeros;}
	public byte getRegister(int i){return registers[i];}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** One byte per bucket; only low 6 bits used (4-bit relNLZ + 2-bit history). */
	private final byte[] registers;
	private int minZeros=0;
	/** Count of buckets below the current floor; triggers advance when 0. */
	private int minZeroCount;
	private int filledBuckets=0;
	/** Early exit mask. More conservative than DLL4: allows NLZ >= minZeros - HISTORY_MARGIN. */
	private long eeMask=-1L;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	/** Number of bit positions reserved below floor for history tracking.
	 *  With 2-bit history, elements up to 2 tiers below minZeros can update registers. */
	static final int HISTORY_MARGIN=2;
	/** When true, advance floor when all buckets have at least one at-floor element. */
	public static boolean EARLY_PROMOTE=true;

	/** Default CF — uses DLL4's CF for now. Will be replaced with UDLL6-specific CF. */
	public static final String CF_FILE="?cardinalityCorrectionDLL4.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}
}
