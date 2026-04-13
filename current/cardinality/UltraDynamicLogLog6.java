package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * Int-packed variant of UltraDynamicLogLog6.
 * Stores 5 six-bit registers per 32-bit int (30 bits used, 2 wasted).
 * At 2048 buckets: ceil(2048/5) = 410 ints = 1,640 bytes.
 * Compared to byte[]: 2,048 bytes → 1,640 bytes (20% savings).
 * <p>
 * Functionally identical to UltraDynamicLogLog6. Exists solely to
 * demonstrate the memory savings of int-packed 6-bit registers for
 * timing comparisons.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class UltraDynamicLogLog6 extends CardinalityTracker {

	final int[] packed;
	/** Number of ints in packed array. */
	private final int packedLen;
	private int minZeros=0;
	private int floorCount;
	private long eeMask=-1L;
	private int filledBuckets=0;
	/** Lazy-allocated per-bucket LC history state indices for lcHist(). */
	private byte[] sbsStates;

	public UltraDynamicLogLog6(){this(2048, 31, -1, 0);}
	UltraDynamicLogLog6(Parser p){
		super(p);
		packedLen=(buckets+4)/5;
		packed=new int[packedLen];
		floorCount=buckets;
	}
	UltraDynamicLogLog6(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		packedLen=(buckets+4)/5;
		packed=new int[packedLen];
		floorCount=buckets;
	}
	@Override public UltraDynamicLogLog6 copy(){return new UltraDynamicLogLog6(buckets, k, -1, minProb);}

	/** Read 6-bit register i from packed int array. */
	int getReg(int i){
		return (packed[i/5]>>>((i%5)*6))&0x3F;
	}

	/** Write 6-bit register i into packed int array. */
	private void setReg(int i, int val){
		final int idx=i/5;
		final int shift=(i%5)*6;
		packed[idx]=(packed[idx] & ~(0x3F<<shift)) | ((val&0x3F)<<shift);
	}

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);

		if(Long.compareUnsigned(key, eeMask)>0){return;}
		branch1++;

		final int idx=(int)(key&bucketMask);
		final int nlz=Long.numberOfLeadingZeros(key);
		final int relNlz=nlz-minZeros;

		final int bitPos=relNlz+HISTORY_MARGIN;
		if(bitPos<0||bitPos>=64){return;}

		final int oldReg=getReg(idx);
		long hashPrefix=ErtlULL.unpack((byte)oldReg);
		hashPrefix|=1L<<bitPos;
		final int rawReg=ErtlULL.pack(hashPrefix)&0xFF;
		final int newReg;
		if(rawReg>MAX_REGISTER){
			final int sub=(int)((hashPrefix>>>13)&3);
			newReg=(15<<2)|sub;
		}else{
			newReg=rawReg;
		}
		if(newReg<=oldReg){return;}
		branch2++;

		setReg(idx, newReg);
		lastCardinality=-1;

		if(oldReg==0){filledBuckets++;}

		final int oldNlzPart=(oldReg>>>2);
		final int newNlzPart=(newReg>>>2);
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

	private int countAndDecrement(){
		int newFloorCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg==0){newFloorCount++; continue;}
			final int newReg=reg-4;
			setReg(i, newReg);
			if((newReg>>>2)<=HISTORY_MARGIN){newFloorCount++;}
		}
		return newFloorCount;
	}

	public double fgraEstimate(){
		if(!CALC_FGRA){return 0;}
		final int p=bucketBits;
		final int regOffset=4*(minZeros+p-1-HISTORY_MARGIN);
		final byte[] ertlRegs=new byte[buckets];
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
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
		final CardStats s=summarize();
		long card=Math.max(0, Math.round(s.ldlc()));
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((UltraDynamicLogLog6)log);
	}

	public void add(UltraDynamicLogLog6 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		final int newMinZeros=Math.max(minZeros, log.minZeros);
		final int deltaA=newMinZeros-minZeros;
		final int deltaB=newMinZeros-log.minZeros;
		for(int i=0; i<buckets; i++){
			final int rA=getReg(i);
			final int rB=log.getReg(i);
			long hpA=(rA==0) ? 0 : ErtlULL.unpack((byte)rA)>>>deltaA;
			long hpB=(rB==0) ? 0 : ErtlULL.unpack((byte)rB)>>>deltaB;
			long merged=hpA|hpB;
			if(merged==0){
				setReg(i, 0);
			}else{
				setReg(i, Math.min(ErtlULL.pack(merged)&0xFF, MAX_REGISTER));
			}
		}
		minZeros=newMinZeros;
		filledBuckets=0;
		floorCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg>0){filledBuckets++;}
			if(reg==0 || (reg>>>2)<=HISTORY_MARGIN){floorCount++;}
		}
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

	@Override
	public double[] rawEstimates(){
		lastSummarized=summarize();
		final double hybridEst=lastSummarized.hybridDLL();
		return AbstractCardStats.buildLegacyArray(lastSummarized, hybridEst);
	}

	/** Cached CardStats from rawEstimates(), consumed by ldlcEstimate(). */
	private CardStats lastSummarized;

	/** Build a CardStats from the register array.
	 *  Counts-only for now; per-state tier multipliers and sbsStates will
	 *  move to CardStats phase 2a when StateTable is implemented. */
	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[buckets];
		final int phantomNlz=minZeros-1;
		int filledCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg==0){
				// Truly empty or post-decrement empty: phantom if minZeros>0
				if(minZeros>0 && phantomNlz>=0 && phantomNlz<64){
					nlzCounts[phantomNlz+1]++;
					filledCount++;
					// No history for empty phantoms
				}
			}else{
				final int nlzPart=reg>>>2;
				final int histPattern=reg&3;
				if(nlzPart>=HISTORY_MARGIN){
					// Normal filled bucket: absNlz = nlzPart - HISTORY_MARGIN + minZeros
					final int absNlz=nlzPart-HISTORY_MARGIN+minZeros;
					if(absNlz<64){
						nlzCounts[absNlz+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((absNlz+1)<<2)|histPattern);
					}
				}else{
					// Sub-margin bucket (after decrement): phantom at phantomNlz.
					// History bits may still be valid.
					if(minZeros>0 && phantomNlz>=0 && phantomNlz<64){
						nlzCounts[phantomNlz+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((phantomNlz+1)<<2)|histPattern);
					}
				}
			}
		}
		nlzCounts[0]=buckets-filledCount;
		return new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
	}

	/** Public accessor for register value. */
	public int getRegPublic(int i){return getReg(i);}
	/** Alias for getRegPublic. */
	public int getRegister(int i){return getReg(i);}

	/** Legacy array format for backward compatibility.
	 *  Returns {ldlc, dlcSbs, hc, lcMin, fgra, hll, meanH, hybridPlus2}. */
	public double[] udlcEstimate(){
		final CardStats s=(lastSummarized!=null) ? lastSummarized : summarize();
		return new double[]{s.ldlc(), s.dlcSbs(), s.hc(), s.lcMin(),
			fgraEstimate(), s.hllRaw(), s.meanHistCF(), s.hybridPlus2()};
	}

	/** Debug: print all register values. */
	public void printRegisters(){
		for(int i=0; i<buckets; i++){System.err.print(getReg(i)+" ");}
		System.err.println();
	}

	/** Memory used by the packed array in bytes. */
	public int packedBytes(){return packedLen*4;}

	/*--------------------------------------------------------------*/
	/*----------------     LDLC Estimation          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the cached CardStats from the last rawEstimates() call, then clears it. */
	public CardStats consumeLastSummarized(){
		final CardStats s=lastSummarized;
		lastSummarized=null;
		return s;
	}

	private static final int wordlen=64;
	static final int HISTORY_MARGIN=2;
	static final int MAX_REGISTER=63;

	/** When true, clamp register at MAX_REGISTER instead of wrapping. */
	public static boolean SATURATE_ON_OVERFLOW=true;

	/** When true, fgraEstimate() is calculated. Default false to avoid per-call overhead.
	 *  ErtlULL always calculates FGRA regardless of this flag. */
	public static boolean CALC_FGRA=false;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/** Default resource file for UDLL6 correction factors.
	 *  Placeholder: uses DLL4's table until a UDLL6-specific table is generated. */
	public static final String CF_FILE="?cardinalityCorrectionUDLL6.tsv.gz";
	/** Bucket count used to build CF_MATRIX (for interpolation). */
	private static int CF_BUCKETS=2048;
	/** Per-class correction factor matrix; null until initializeCF() is called. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	/** Loads the UDLL6 correction factor matrix from CF_FILE. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Asymptotic meanRaw/trueCard ratio, measured 128k ddls maxmult=4096 (Apr 13 2026). */
	@Override public float terminalMeanCF(){return 0.721123f;}

	/** Asymptotic Mean+H ratio (2-bit history), measured 128k ddls maxmult=4096 (Apr 13 2026).
	 *  Replaces the old StateTable.terminalCF(2,0)=0.7726 which was measured under
	 *  different conditions. */
	@Override public float terminalMeanPlusCF(){return 0.856020f;}
}
