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
public final class UltraDynamicLogLog6i extends CardinalityTracker {

	final int[] packed;
	/** Number of ints in packed array. */
	private final int packedLen;
	private int minZeros=0;
	private int floorCount;
	private long eeMask=-1L;
	private int filledBuckets=0;

	UltraDynamicLogLog6i(){this(2048, 31, -1, 0);}
	UltraDynamicLogLog6i(Parser p){
		super(p);
		packedLen=(buckets+4)/5;
		packed=new int[packedLen];
		floorCount=buckets;
	}
	UltraDynamicLogLog6i(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		packedLen=(buckets+4)/5;
		packed=new int[packedLen];
		floorCount=buckets;
	}
	@Override public UltraDynamicLogLog6i copy(){return new UltraDynamicLogLog6i(buckets, k, -1, minProb);}

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
		add((UltraDynamicLogLog6i)log);
	}

	public void add(UltraDynamicLogLog6i log){
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
		final int total=11+4+CardinalityStats.NUM_DLC_TIERS;
		final double[] r=new double[total];
		final double fgra=fgraEstimate();
		r[0]=fgra; r[1]=fgra; r[4]=fgra; r[6]=fgra; r[8]=fgra;
		return r;
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
		final int[] nlzCounts=new int[64];
		int maxNlz=0, emptyCount=0, nonEmpty=0, hasHistory=0;

		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg==0){
				absNlzArr[i]=-1;
				emptyCount++;
			}else{
				final int absNlz=(reg>>>2)-HISTORY_MARGIN+minZeros;
				absNlzArr[i]=absNlz;
				histArr[i]=reg&3;
				if(absNlz>maxNlz) maxNlz=absNlz;
				if(absNlz>=0 && absNlz<64) nlzCounts[absNlz]++;
				nonEmpty++;
				if((reg&3)!=0) hasHistory++;
			}
		}

		// lcMin: tier-compensated LC at lowest non-full tier
		// microIndex detects collisions: when microFilled > nonEmpty, some elements
		// landed in the same main bucket but different micro positions.
		final int microFilled=(int)Long.bitCount(microIndex);
		final int V=buckets-Math.max(nonEmpty, microFilled);
		double lcMin;
		if(V>0){
			lcMin=(double)buckets*Math.log((double)buckets/Math.max(V, 0.5));
		}else{
			lcMin=0;
			int vk=0;
			for(int k=0; k<nlzCounts.length; k++){
				vk+=nlzCounts[k];
				if(vk>0){
					lcMin=(1L<<(k+1))*(double)buckets*Math.log((double)buckets/Math.max(vk, 0.5));
					break;
				}
			}
		}

		// DLC: cumulative log-space blend, with lcMin fallback
		final double dlc=dlcWithLcMin(nlzCounts, V, buckets, lcMin);

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
					final double w=(double)hcBeff/buckets;
					hcSumW+=w;
					hcSumWLogE+=w*Math.log(est);
				}
			}
		}
		final double hc=(hcSumW>0 ? Math.exp(hcSumWLogE/hcSumW) : 0);

		// LDLC blend: 60% DLC + 40% HC, with adaptive ramp
		double ldlc;
		if(hc<=0 || dlc<=0){
			ldlc=dlc;
		}else{
			final double coverage=(nonEmpty>0 ? (double)hasHistory/nonEmpty : 0);
			final double maxHcWeight=0.40;
			final double ramp=Math.max(0, Math.min(1, (coverage-0.30)/0.50));
			final double hcWeight=maxHcWeight*ramp;
			ldlc=(1.0-hcWeight)*dlc+hcWeight*hc;
		}

		// MicroIndex floor: LC over 64 virtual buckets (collision resistance at very low C)
		final int microEmpty=Math.max(Math.min(63, 64-microFilled), 0);
		final double microEst=(microIndex==0 ? 0 :
			64*Math.log((double)64/Math.max(microEmpty, 0.5)));

		// HLL: standard HyperLogLog from max NLZ values
		// Register value = absNlz + 1 for filled buckets, 0 for empty
		// hllSum = Σ 2^{-register}
		double hllSum=0;
		for(int i=0; i<buckets; i++){
			final int reg=(absNlzArr[i]>=0) ? absNlzArr[i]+1 : 0;
			hllSum+=(reg==0 ? 1.0 : Math.pow(2, -reg));
		}
		final double alpha_m=0.7213/(1.0+1.079/buckets);
		double hll=alpha_m*(double)buckets*(double)buckets/hllSum;
		// LC correction at low cardinality (standard HLL range correction)
		if(hll<=2.5*buckets && V>0){
			hll=(double)buckets*Math.log((double)buckets/(double)V);
		}

		final double fgra=fgraEstimate();
		return new double[]{ldlc, dlc, hc, lcMin, fgra, hll};
	}

	/** DLC log-space blend with lcMin fallback at low occupancy. */
	private double dlcWithLcMin(int[] nlzCounts, int V, int B, double lcMin){
		if(V>=B/2){return lcMin;}
		final double target=B*0.25;
		final double alpha=9.0/B;
		int vk=V;
		double sumW=0, sumWLogE=0;
		for(int tier=0; tier<64; tier++){
			if(vk>=1 && vk<B){
				final double est=(1L<<tier)*(double)B*Math.log((double)B/Math.max(vk, 0.5));
				final double w=Math.exp(-alpha*Math.abs(vk-target));
				sumW+=w;
				sumWLogE+=w*Math.log(est);
			}
			if(tier<nlzCounts.length) vk+=nlzCounts[tier];
			if(vk>=B) break;
		}
		if(sumW<=0){return lcMin;}
		final double blendEst=Math.exp(sumWLogE/sumW);
		if(V>0.3*B){
			final double t=((double)V-0.3*B)/(0.2*B);
			return t*lcMin+(1-t)*blendEst;
		}
		return blendEst;
	}

	private static final int wordlen=64;
	static final int HISTORY_MARGIN=2;
	static final int MAX_REGISTER=63;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}
}
