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
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

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
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0);
	}

	/** Public accessor for register value. */
	public int getRegPublic(int i){return getReg(i);}
	/** Alias for getRegPublic. */
	public int getRegister(int i){return getReg(i);}

	/** Alias for ldlcEstimate() for backward compatibility. */
	public double[] udlcEstimate(){return ldlcEstimate();}

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

	/**
	 * LDLC (Layered Dynamic Linear Counting).
	 * Blends surface DLC (cumulative) with history-only per-tier LC.
	 * At low cardinality, falls back to lcMin (tier-compensated LC).
	 * Returns {ldlc, dlc, hc, lcMin, fgra}.
	 */
	public double[] ldlcEstimate(){
		final int[] absNlzArr=new int[buckets];
		final int[] histArr=new int[buckets];
		final int[] nlzCounts=new int[66];
		int maxNlz=0, emptyCount=0, nonEmpty=0, hasHistory=0;
		int histVirtualTotal=0, histVirtualFilled=0;
		final int hbits=2; // UDLL6 always has 2 history bits

		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg==0){
				absNlzArr[i]=-1;
				emptyCount++;
			}else{
				final int absNlz=(reg>>>2)-HISTORY_MARGIN+minZeros;
				absNlzArr[i]=absNlz;
				final int hist=reg&3;
				histArr[i]=hist;
				if(absNlz>maxNlz) maxNlz=absNlz;
				if(absNlz>=0 && absNlz<65) nlzCounts[absNlz+1]++;
				nonEmpty++;
				if(hist!=0) hasHistory++;
				// Virtual buckets: tier t has min(hbits, max(0, t-1)) valid slots
				final int validSlots=Math.min(hbits, Math.max(0, absNlz-1));
				histVirtualTotal+=validSlots;
				if(validSlots>0){
					final int validMask=((1<<validSlots)-1)<<(hbits-validSlots);
					histVirtualFilled+=Integer.bitCount(hist&validMask);
				}
			}
		}

		// History-based LC floor: filled buckets + set history bits in valid positions
		final int lcFloor=(CardinalityTracker.USE_HISTORY_FOR_LC ? nonEmpty+histVirtualFilled : 0);

		// lcMin: tier-compensated LC at lowest non-full tier
		// microIndex detects collisions: when microFilled > nonEmpty, some elements
		// landed in the same main bucket but different micro positions.
		final int microFilled=(int)Long.bitCount(microIndex);
		final int V=buckets-Math.max(nonEmpty, microFilled);
		double lcMin;
		if(V>0){
			lcMin=Math.max((double)buckets*Math.log((double)buckets/Math.max(V, 0.5)), lcFloor);
		}else{
			double lcMinTmp=0;
			int vk=0;
			for(int k=0; k+1<nlzCounts.length; k++){
				vk+=nlzCounts[k+1];
				if(vk>0){
					lcMinTmp=(1L<<(k+1))*(double)buckets*Math.log((double)buckets/Math.max(vk, 0.5));
					break;
				}
			}
			lcMin=Math.max(lcMinTmp, lcFloor);
		}

		// Get DLCHist and lcHist from CardStats
		final CardStats s=summarize();
		final double sbs=s.sbs();
		final double dlc=s.dlcSbs(); // DLC with SBS fallback

		// HC: history-only per-tier exact LC
		final int[] nlzBucketCount=new int[64];
		final int[] nlzH1Set=new int[64];
		final int[] nlzH0Set=new int[64];
		for(int i=0; i<buckets; i++){
			final int n=absNlzArr[i];
			if(n>=0 && n<64){
				nlzBucketCount[n]++;
				if((histArr[i]&2)!=0) nlzH1Set[n]++;
				if((histArr[i]&1)!=0) nlzH0Set[n]++;
			}
		}

		final int maxTier=Math.min(maxNlz+2, 63);
		double hcSumW=0, hcSumWLogE=0;
		for(int t=0; t<=maxTier; t++){
			int hcBeff=0, hcUnseen=0;
			if(t+1<64){
				hcBeff+=nlzBucketCount[t+1];
				hcUnseen+=(nlzBucketCount[t+1]-nlzH1Set[t+1]);
			}
			if(t+2<64){
				hcBeff+=nlzBucketCount[t+2];
				hcUnseen+=(nlzBucketCount[t+2]-nlzH0Set[t+2]);
			}
			if(hcBeff>=8 && hcUnseen>=1 && hcUnseen<hcBeff){
				final double est=(1L<<(t+1))*(double)buckets*Math.log((double)hcBeff/hcUnseen);
				if(est>0 && !Double.isNaN(est)){
					final int hcOcc=hcBeff-hcUnseen;
					final double hcErr=AbstractCardStats.SQRT_2_OVER_PI
						*Math.sqrt((double)hcOcc/((double)hcBeff*hcUnseen))
						/Math.log((double)hcBeff/hcUnseen);
					final double w=Math.pow(1.0/hcErr, AbstractCardStats.HC_INFO_POWER);
					hcSumW+=w;
					hcSumWLogE+=w*Math.log(est);
				}
			}
		}
		final double hc=(hcSumW>0 ? Math.exp(hcSumWLogE/hcSumW) : 0);

		// LDLC: double blend using DLC as zone detector.
		// Blend A: LCHist → DLC over [aLo*B, aHi*B]
		// Blend B: HC ramps in over [bLo*B, bHi*B] to maxHcWeight
		// Overlap: average both blends
		final double B=(double)buckets;
		final double aLo=AbstractCardStats.LDLC_A_LO*B, aHi=AbstractCardStats.LDLC_A_HI*B;
		final double bLo=AbstractCardStats.LDLC_B_LO*B, bHi=AbstractCardStats.LDLC_B_HI*B;
		final double maxHcWeight=CardinalityTracker.LDLC_HC_WEIGHT;
		final boolean hcUsable=(hc>0 && dlc>0);
		double ldlc;
		if(dlc<=aLo){
			ldlc=sbs;
		}else if(dlc<=bLo || !hcUsable){
			// Only Blend A active (or HC unavailable)
			final double tA=Math.min(1, (dlc-aLo)/(aHi-aLo));
			ldlc=(1-tA)*sbs+tA*dlc;
		}else if(dlc<=aHi){
			// Both blends active: average
			final double tA=(dlc-aLo)/(aHi-aLo);
			final double blendA=(1-tA)*sbs+tA*dlc;
			final double tB=Math.min(1, (dlc-bLo)/(bHi-bLo));
			final double hcW=tB*maxHcWeight;
			final double blendB=(1-hcW)*dlc+hcW*hc;
			ldlc=(blendA+blendB)*0.5;
		}else if(dlc<=bHi){
			// Only Blend B: HC ramping in
			final double tB=(dlc-bLo)/(bHi-bLo);
			final double hcW=tB*maxHcWeight;
			ldlc=(1-hcW)*dlc+hcW*hc;
		}else{
			// Final ratio
			ldlc=(1-maxHcWeight)*dlc+maxHcWeight*hc;
		}

		// HLL: standard HyperLogLog from max NLZ values
		double hllSum=0;
		for(int i=0; i<buckets; i++){
			final int reg=(absNlzArr[i]>=0) ? absNlzArr[i]+1 : 0;
			hllSum+=(reg==0 ? 1.0 : Math.pow(2, -reg));
		}
		final double alpha_m=0.7213/(1.0+1.079/buckets);
		double hll=alpha_m*(double)buckets*(double)buckets/hllSum;
		if(hll<=2.5*buckets && V>0){
			hll=(double)buckets*Math.log((double)buckets/(double)V);
		}

		final double fgra=fgraEstimate();
		final double meanH=s.meanHistCF();
		final double hybridPlus2=s.hybridPlus2();
		return new double[]{ldlc, dlc, hc, lcMin, fgra, hll, meanH, hybridPlus2};
	}

	private static final int wordlen=64;
	static final int HISTORY_MARGIN=2;
	static final int MAX_REGISTER=63;

	/** When true, clamp register at MAX_REGISTER instead of wrapping. */
	public static boolean SATURATE_ON_OVERFLOW=true;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/** Terminal CF for mean estimator with 2-bit history.
	 *  At high cardinality, the uncorrected mean converges to trueCard/MEAN_TERMINAL_CF.
	 *  Folded into tierMult during summarize() so the external CF table converges to ~1.0.
	 *  Will be refined once UDLL6-specific CF tables are generated. */
	public static double MEAN_TERMINAL_CF=0.929224472;

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
}
