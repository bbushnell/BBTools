package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * UltraDynamicLogLog6: 6-bit registers combining DLL4's dynamic relative-NLZ
 * architecture with ULL's 2-bit sub-NLZ history and FGRA estimator.
 * <p>
 * Derived from ULLd (relative 8-bit) by clamping registers to [0,63].
 * Register encoding: 6 bits = 4-bit relative NLZ + 2-bit history.
 * 15 usable NLZ tiers (nlzPart 2..16, but clamped to 15 = 63>>>2).
 * The 2-bit history provides ~18% variance reduction via FGRA estimator.
 * <p>
 * At 2048 buckets with byte[]: 2048 bytes. Same memory as DLL4 at 4096 buckets.
 * Future: int-packed (5 per int) with modulo bucket count for exact 2KB.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class UltraDynamicLogLog6 extends CardinalityTracker {

	final byte[] registers;
	private int minZeros=0;
	/** Buckets with nlzPart <= HISTORY_MARGIN (relNlz <= 0) or empty. When 0, advance. */
	private int floorCount;
	private long eeMask=-1L;
	private int filledBuckets=0;

	UltraDynamicLogLog6(){this(2048, 31, -1, 0);}
	UltraDynamicLogLog6(Parser p){
		super(p);
		registers=new byte[buckets];
		floorCount=buckets;
	}
	UltraDynamicLogLog6(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		registers=new byte[buckets];
		floorCount=buckets;
	}
	@Override public UltraDynamicLogLog6 copy(){return new UltraDynamicLogLog6(buckets, k, -1, minProb);}

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);

		// Early exit: reject elements whose NLZ < minZeros - HISTORY_MARGIN.
		if(Long.compareUnsigned(key, eeMask)>0){return;}
		branch1++;

		final int idx=(int)(key&bucketMask);
		final int nlz=Long.numberOfLeadingZeros(key);
		final int relNlz=nlz-minZeros;

		// Bit position in hashPrefix: relNlz + HISTORY_MARGIN.
		final int bitPos=relNlz+HISTORY_MARGIN;
		if(bitPos<0||bitPos>=64){return;}

		final byte oldReg=registers[idx];
		long hashPrefix=ErtlULL.unpack(oldReg);
		hashPrefix|=1L<<bitPos;
		final int rawReg=ErtlULL.pack(hashPrefix)&0xFF;
		final byte newReg;
		if(rawReg>MAX_REGISTER){
			if(SATURATE_ON_OVERFLOW){
				newReg=(byte)MAX_REGISTER; // 63 = nlzPart=15, sub=11
			}else{
				// Clamp nlzPart to 15, extract actual sub-bits from hashPrefix
				// In hashPrefix, bit N = main, bit N-1 = sub MSB, bit N-2 = sub LSB.
				// For nlzPart=15: main at bit 15, sub at bits 14,13.
				final int sub=(int)((hashPrefix>>>13)&3);
				newReg=(byte)((15<<2)|sub);
			}
		}else{
			newReg=(byte)rawReg;
		}
		if((newReg&0xFF)<=(oldReg&0xFF)){return;}
		branch2++;

		registers[idx]=newReg;
		lastCardinality=-1;

		if(oldReg==0){filledBuckets++;}

		// Floor tracking: "at floor" = nlzPart <= HISTORY_MARGIN (relNlz <= 0) or empty.
		final int oldNlzPart=((oldReg&0xFF)>>>2);
		final int newNlzPart=((newReg&0xFF)>>>2);
		if(oldNlzPart<=HISTORY_MARGIN && newNlzPart>HISTORY_MARGIN){
			floorCount--;
			if(floorCount<=0){
				while(floorCount==0 && minZeros<wordlen){
					minZeros++;
					floorCount=countAndDecrement();
				}
				int exitThreshold=Math.max(0, minZeros-HISTORY_MARGIN);
				eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
			}
		}
	}

	/**
	 * Subtract 4 from each non-empty register (decrease nlzPart by 1,
	 * preserving sub-bits). Called only when floorCount==0, meaning all
	 * registers have nlzPart > HISTORY_MARGIN, so subtracting 4 never
	 * creates negative values.
	 *
	 * @return new floorCount (registers with nlzPart <= HISTORY_MARGIN after decrement)
	 */
	private int countAndDecrement(){
		int newFloorCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=registers[i]&0xFF;
			if(reg==0){newFloorCount++; continue;}
			final int newReg=reg-4;
			registers[i]=(byte)newReg;
			if((newReg>>>2)<=HISTORY_MARGIN){newFloorCount++;}
		}
		return newFloorCount;
	}

	/**
	 * Converts 6-bit relative registers to absolute Ertl format and runs FGRA.
	 * ertl_reg = relative_reg + 4*(minZeros + p - 1 - HISTORY_MARGIN)
	 */
	public double fgraEstimate(){
		final int p=bucketBits;
		final int regOffset=4*(minZeros+p-1-HISTORY_MARGIN);
		final byte[] ertlRegs=new byte[buckets];
		for(int i=0; i<buckets; i++){
			final int reg=registers[i]&0xFF;
			if(reg==0){
				ertlRegs[i]=0;
			}else{
				ertlRegs[i]=(byte)Math.min(reg+regOffset, 255);
			}
		}
		return ErtlULL.fgraEstimateStatic(ertlRegs, p);
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		long card=Math.max(0, Math.round(fgraEstimate()));
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
	 * Converts both to absolute hashPrefixes, ORs them, converts back.
	 * Per-register max is WRONG for ULL because sub-bits encode independent
	 * history that may differ between instances for the same bucket.
	 */
	public void add(UltraDynamicLogLog6 log){
		added+=log.added;
		lastCardinality=-1;
		// Re-frame both to newMinZeros by shifting hashPrefixes RIGHT by delta,
		// then OR to merge history bits, then repack.
		final int newMinZeros=Math.max(minZeros, log.minZeros);
		final int deltaA=newMinZeros-minZeros;  // >= 0
		final int deltaB=newMinZeros-log.minZeros;  // >= 0
		for(int i=0; i<buckets; i++){
			final int rA=registers[i]&0xFF;
			final int rB=log.registers[i]&0xFF;
			long hpA=(rA==0) ? 0 : ErtlULL.unpack((byte)rA)>>>deltaA;
			long hpB=(rB==0) ? 0 : ErtlULL.unpack((byte)rB)>>>deltaB;
			long merged=hpA|hpB;
			if(merged==0){
				registers[i]=0;
			}else{
				registers[i]=(byte)Math.min(ErtlULL.pack(merged)&0xFF, MAX_REGISTER);
			}
		}
		minZeros=newMinZeros;
		// Recount filledBuckets and floorCount
		filledBuckets=0;
		floorCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=registers[i]&0xFF;
			if(reg>0){filledBuckets++;}
			if(reg==0 || (reg>>>2)<=HISTORY_MARGIN){floorCount++;}
		}
		// Advance if possible
		while(floorCount==0 && minZeros<wordlen){
			minZeros++;
			floorCount=countAndDecrement();
		}
		int exitThreshold=Math.max(0, minZeros-HISTORY_MARGIN);
		eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
	}
	@Override public final float[] compensationFactorLogBucketsArray(){return null;}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}
	public int getMinZeros(){return minZeros;}
	public byte getRegister(int i){return registers[i];}

	@Override
	public double[] rawEstimates(){
		final int total=11+4+CardinalityStats.NUM_DLC_TIERS;
		final double[] r=new double[total];
		final double fgra=fgraEstimate();
		r[0]=fgra; r[1]=fgra; r[4]=fgra; r[6]=fgra; r[8]=fgra;
		return r;
	}

	private static final int wordlen=64;
	static final int HISTORY_MARGIN=2;
	/** Maximum 6-bit register value. nlzPart can be 0-15, sub 0-3. */
	static final int MAX_REGISTER=63;
	/** When true, overflow sets history to 11 (maximum). When false, preserves actual sub-bits. */
	public static boolean SATURATE_ON_OVERFLOW=false;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}
}
