package cardinality;

import shared.Tools;

/**
 * ULLd: ULLc converted to RELATIVE encoding with 8-bit registers.
 * Intermediate step between ULLc (absolute) and UDLL6 (relative, 6-bit).
 * <p>
 * Registers store relative NLZ (nlz - minZeros) using Ertl's pack/unpack
 * on a hashPrefix bitmask offset by HISTORY_MARGIN. When minZeros advances,
 * countAndDecrement subtracts 4 from each register (= decrease nlzPart by 1,
 * preserving sub-bits). This is equivalent to shifting the hashPrefix
 * bitmask right by 1, which correctly maintains the 2-bit history.
 * <p>
 * Must produce EXACT same FGRA estimates as ULLc at all cardinalities.
 * For FGRA, converts relative registers to absolute by adding
 * 4*(minZeros + p - 1 - HISTORY_MARGIN).
 */
public final class ULLd extends CardinalityTracker {

	final byte[] registers;
	private int minZeros=0;
	/** Buckets with nlzPart <= HISTORY_MARGIN (relNlz <= 0) or empty. When 0, advance. */
	private int floorCount;
	private long eeMask=-1L;
	private int filledBuckets=0;

	ULLd(){this(2048, 31, -1, 0);}
	ULLd(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		registers=new byte[buckets];
		floorCount=buckets;
	}
	@Override public ULLd copy(){return new ULLd(buckets, k, -1, minProb);}

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
		// Reserves positions 0..(HISTORY_MARGIN-1) for below-floor history.
		final int bitPos=relNlz+HISTORY_MARGIN;
		if(bitPos<0||bitPos>=64){return;}

		final byte oldReg=registers[idx];
		long hashPrefix=ErtlULL.unpack(oldReg);
		hashPrefix|=1L<<bitPos;
		final byte newReg=ErtlULL.pack(hashPrefix);
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
	 * preserving sub-bits). This is equivalent to shifting the hashPrefix
	 * bitmask right by 1 — the sub-bits maintain their correct meaning
	 * because they are relative to the register's max NLZ.
	 * <p>
	 * Called only when floorCount==0, meaning all registers have
	 * nlzPart > HISTORY_MARGIN, so subtracting 4 never creates negative values.
	 *
	 * @return new floorCount (registers with nlzPart <= HISTORY_MARGIN after decrement)
	 */
	private int countAndDecrement(){
		int newFloorCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=registers[i]&0xFF;
			if(reg==0){newFloorCount++; continue;} // shouldn't happen when floorCount was 0
			final int newReg=reg-4;
			registers[i]=(byte)newReg;
			if((newReg>>>2)<=HISTORY_MARGIN){newFloorCount++;}
		}
		return newFloorCount;
	}

	/**
	 * Converts relative registers to absolute Ertl format and runs FGRA.
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

	@Override public final void add(CardinalityTracker log){throw new UnsupportedOperationException();}
	@Override public final float[] compensationFactorLogBucketsArray(){return null;}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}
	public int getMinZeros(){return minZeros;}

	@Override
	public double[] rawEstimates(){
		final int total=11+6+CardinalityStats.NUM_DLC_TIERS;
		final double[] r=new double[total];
		final double fgra=fgraEstimate();
		r[0]=fgra; r[1]=fgra; r[4]=fgra; r[6]=fgra; r[8]=fgra;
		return r;
	}

	private static final int wordlen=64;
	static final int HISTORY_MARGIN=2;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}
}
