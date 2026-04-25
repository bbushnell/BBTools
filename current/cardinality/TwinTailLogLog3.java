package cardinality;

import shared.Tools;

/**
 * TTLL3: TwinTailLogLog with 3-bit tails — ceiling-based dual-history
 * cardinality estimator with 10-bit registers packed 3 per 32-bit int.
 *
 * Each register is 10 bits:
 *   [9:6] = 4-bit localExp (offset from globalExp)
 *   [5:3] = h1: 3-bit tail for histBit==1 (MSB=tier, mid=tier-1, LSB=tier-2)
 *   [2:0] = h0: 3-bit tail for histBit==0 (same semantics)
 *
 * Packing: 3 registers per int (30 bits used, 2 unused).
 * Default: 1536 buckets = 512 ints = 2KB.
 *
 * Bucket selection: modulo (non-power-of-2).
 *   bucket = Long.remainderUnsigned(key, numBuckets)
 *   histBit = Long.divideUnsigned(key, numBuckets) &amp; 1
 *
 * Update rule (symmetric only):
 *   delta = hashNLZ - (globalExp + localExp)
 *   delta &lt; -2 : ignore
 *   delta in [-2, 0]: set bit (2+delta) of history[histBit]
 *   delta &gt; 0 : advance localExp; shift both tails right by min(delta,3);
 *                set MSB of history[histBit]
 *   Overflow (newLocalExp &gt; 15): cap at 15.
 *
 * combined_h = (h1 &lt;&lt; 3) | h0 = 6-bit, 64 possible states.
 *
 * @author Brian, Chloe
 * @date April 25, 2026
 */
public final class TwinTailLogLog3 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int NUM_TIERS=16;
	static final int TAIL_BITS=3;
	static final int TAIL_MASK=0x7;
	static final int EXP_SHIFT=2*TAIL_BITS; // 6
	static final int TAILS_MASK=(1<<EXP_SHIFT)-1; // 0x3F
	static final int REG_BITS=EXP_SHIFT+4; // 10
	static final int REG_MASK=(1<<REG_BITS)-1; // 0x3FF
	static final int NUM_COMBINED=1<<EXP_SHIFT; // 64

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	TwinTailLogLog3(){this(1536, 31, -1, 0);}

	TwinTailLogLog3(parse.Parser p){
		super(p);
		numBuckets=p.loglogbuckets;
		regs=new int[(numBuckets+2)/3];
		numFloorBuckets=numBuckets;
	}

	TwinTailLogLog3(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		numBuckets=buckets_;
		regs=new int[(numBuckets+2)/3];
		numFloorBuckets=numBuckets;
	}

	@Override
	public TwinTailLogLog3 copy(){return new TwinTailLogLog3(numBuckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------       Register Access        ----------------*/
	/*--------------------------------------------------------------*/

	private int getReg(int bucket){
		return (regs[bucket/3]>>>((bucket%3)*REG_BITS))&REG_MASK;
	}

	private void setReg(int bucket, int val){
		final int word=bucket/3;
		final int shift=(bucket%3)*REG_BITS;
		regs[word]=(regs[word]&~(REG_MASK<<shift))|((val&REG_MASK)<<shift);
	}

	private static int getRegStatic(int[] regs, int bucket){
		return (regs[bucket/3]>>>((bucket%3)*REG_BITS))&REG_MASK;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Core Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);

		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int hashNLZ=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(Long.remainderUnsigned(key, numBuckets));
		final int histBit=(int)(Long.divideUnsigned(key, numBuckets)&1);

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int reg=getReg(bucket);
		final int localExp=(reg>>>EXP_SHIFT)&0xF;
		final int absNLZ=globalExp+localExp;
		final int delta=hashNLZ-absNLZ;

		if(delta<-(TAIL_BITS-1)){return;}

		if(delta>0){
			final int newLocalExp=Math.min(localExp+delta, 15);
			final boolean wasZeroExp=(localExp==0);
			final int shiftAmt=Math.min(delta, TAIL_BITS);
			int h0=((reg&TAIL_MASK)>>>shiftAmt);
			int h1=(((reg>>>TAIL_BITS)&TAIL_MASK)>>>shiftAmt);
			if(histBit==0){h0|=(1<<(TAIL_BITS-1));}
			else          {h1|=(1<<(TAIL_BITS-1));}
			lastCardinality=-1;
			lastSummarized=null;
			setReg(bucket, (newLocalExp<<EXP_SHIFT)|(h1<<TAIL_BITS)|h0);
			if(wasZeroExp){
				numFloorBuckets--;
				if(numFloorBuckets==0){advanceGlobal();}
			}
		}else{
			final int bitPos=(TAIL_BITS-1)+delta;
			final int regBitPos=(histBit==0) ? bitPos : (bitPos+TAIL_BITS);
			final int bitToSet=1<<regBitPos;
			if((reg&bitToSet)!=0){return;}
			lastCardinality=-1;
			lastSummarized=null;
			setReg(bucket, reg|bitToSet);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Global Advance         ----------------*/
	/*--------------------------------------------------------------*/

	private void advanceGlobal(){
		while(numFloorBuckets==0){
			globalExp++;
			numFloorBuckets=0;
			for(int i=0; i<numBuckets; i++){
				int reg=getReg(i);
				int le=(reg>>>EXP_SHIFT)&0xF;
				assert le>=1 : "localExp=0 but numFloorBuckets was 0, bucket "+i;
				le--;
				setReg(i, (le<<EXP_SHIFT)|(reg&TAILS_MASK));
				if(le==0){numFloorBuckets++;}
			}
			eeMask=(globalExp<=(TAIL_BITS-1)) ? -1L : (-1L)>>>(globalExp-(TAIL_BITS-1));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Estimation           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		final double rawHyb=s.hybridDLL();
		long card=(long)rawHyb;
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public double[] rawEstimates(){
		lastSummarized=summarize();
		final double hybridEst=lastSummarized.hybridDLL();
		double[] r=AbstractCardStats.buildLegacyArray(lastSummarized, hybridEst);
		if(r.length<=MEAN16_RAW_IDX){
			r=java.util.Arrays.copyOf(r, MEAN16_RAW_IDX+1);
		}
		r[MEAN16_RAW_IDX]=ttllMeanEstimate();
		return r;
	}

	public static final int MEAN16_RAW_IDX=AbstractCardStats.HC_IDX+1;

	public double[] ldlcEstimate(){
		final CardStats s=(lastSummarized!=null) ? lastSummarized : summarize();
		lastSummarized=null;
		return new double[]{s.ldlc(), s.dlcRaw(), s.hc(), s.lcMin(),
				0, s.hllRaw(), s.meanHistCF(), s.hybridPlus2(),
				ttllMeanEstimate(), dualBucketLC()};
	}

	/*--------------------------------------------------------------*/
	/*----------------        Summarize             ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[numBuckets];
		int filledCount=0;

		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0 && globalExp==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;

			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				filledCount++;
				final int combinedH=reg&TAILS_MASK;
				packedBuckets[i]=(char)(((absNlz+1)<<EXP_SHIFT)|combinedH);
			}
		}
		nlzCounts[0]=numBuckets-filledCount;
		lastRawNlz=nlzCounts;

		final CardStats cs=new CardStats(packedBuckets, nlzCounts, 0, EXP_SHIFT, 0, 0,
				numBuckets, microIndex, added, null, 0, 0.0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), 1.0, 1);
		cs.overrideHC(hcTTLL3(regs, globalExp, numBuckets));
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Dual-Bucket LC/DLC       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * LC estimate treating each tail as an independent virtual bucket.
	 * 3-bit tails give finer tier resolution than TTLL's 2-bit version.
	 *
	 * For each tail (h0, h1):
	 *   bit 2 set (MSB): filled at absNlz
	 *   bit 2 unset, bit 1 set: filled at absNlz-1
	 *   bit 2 unset, bit 1 unset, bit 0 set: filled at absNlz-2
	 *   all zero, absNlz &lt; 3: empty virtual bucket
	 *   all zero, absNlz &gt;= 3: unknown, excluded
	 */
	public double dualBucketLC(){
		int virtualTotal=0;
		int virtualEmpty=0;

		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0 && globalExp==0){
				virtualTotal+=2;
				virtualEmpty+=2;
				continue;
			}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;

			for(int t=0; t<2; t++){
				final int tail=(t==0) ? (reg&TAIL_MASK) : ((reg>>>TAIL_BITS)&TAIL_MASK);
				if((tail&4)!=0){
					virtualTotal++;
				}else if((tail&2)!=0){
					virtualTotal++;
				}else if((tail&1)!=0){
					virtualTotal++;
				}else{
					if(absNlz<TAIL_BITS){
						virtualTotal++;
						virtualEmpty++;
					}
				}
			}
		}

		if(virtualTotal==0 || virtualEmpty==0){return 0;}
		return (double)virtualTotal*Math.log((double)virtualTotal/virtualEmpty);
	}

	/*--------------------------------------------------------------*/
	/*----------------    64-State Mean Estimate    ----------------*/
	/*--------------------------------------------------------------*/

	/** Placeholder: Phase 3 will implement with 64-state HSB tables. */
	public double ttllMeanEstimate(){
		return 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------     HC for 3-bit tails       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * HC estimate using 3-bit tail observations.
	 * Three observation types per tail:
	 *   MSB (bit 2): conditional probability model (same as TTLL)
	 *   Mid (bit 1): unconditional half-population LC at depth 1 (same as TTLL LSB)
	 *   LSB (bit 0): unconditional half-population LC at depth 2
	 *
	 * All combined via info-power weighted geometric mean.
	 */
	static double hcTTLL3(int[] regs, int globalExp, int numBuckets){
		final int[] nlzCount=new int[64];
		final int[] bothMSB=new int[64];
		final int[] h0MidSet=new int[64];
		final int[] h1MidSet=new int[64];
		final int[] h0LSBset=new int[64];
		final int[] h1LSBset=new int[64];

		for(int i=0; i<numBuckets; i++){
			final int reg=getRegStatic(regs, i);
			if(reg==0 && globalExp==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;
			if(absNlz<0 || absNlz>=64){continue;}

			final int h0=reg&TAIL_MASK;
			final int h1=(reg>>>TAIL_BITS)&TAIL_MASK;

			nlzCount[absNlz]++;
			if(((h0>>>2)&1)!=0 && ((h1>>>2)&1)!=0){bothMSB[absNlz]++;}
			if(((h0>>>1)&1)!=0){h0MidSet[absNlz]++;}
			if(((h1>>>1)&1)!=0){h1MidSet[absNlz]++;}
			if((h0&1)!=0){h0LSBset[absNlz]++;}
			if((h1&1)!=0){h1LSBset[absNlz]++;}
		}

		int maxNlz=0;
		for(int t=63; t>=0; t--){
			if(nlzCount[t]>0){maxNlz=t; break;}
		}
		final int maxTierHC=Math.min(maxNlz+3, 63);

		double hcSumW=0, hcSumWLogE=0;
		for(int t=0; t<=maxTierHC; t++){
			final double scale=2.0*Math.pow(2.0, t+1)*(double)numBuckets;

			// --- Depth-1: mid bits (h0_mid + h1_mid pooled) ---
			if(t+1<64 && nlzCount[t+1]>0){
				final int beff=2*nlzCount[t+1];
				final int unseen=(nlzCount[t+1]-h0MidSet[t+1])
					+(nlzCount[t+1]-h1MidSet[t+1]);
				if(beff>=8 && unseen>=1 && unseen<beff){
					final double est=scale*Math.log((double)beff/unseen);
					if(est>0 && !Double.isNaN(est)){
						final int occ=beff-unseen;
						final double err=AbstractCardStats.SQRT_2_OVER_PI
							*Math.sqrt((double)occ/((double)beff*unseen))
							/Math.log((double)beff/unseen);
						final double w=Math.pow(1.0/err, AbstractCardStats.HC_INFO_POWER);
						hcSumW+=w;
						hcSumWLogE+=w*Math.log(est);
					}
				}
			}

			// --- Depth-2: LSB bits (h0_LSB + h1_LSB pooled) ---
			if(t+2<64 && nlzCount[t+2]>0){
				final int beff=2*nlzCount[t+2];
				final int unseen=(nlzCount[t+2]-h0LSBset[t+2])
					+(nlzCount[t+2]-h1LSBset[t+2]);
				if(beff>=8 && unseen>=1 && unseen<beff){
					final double est=scale*Math.log((double)beff/unseen);
					if(est>0 && !Double.isNaN(est)){
						final int occ=beff-unseen;
						final double err=AbstractCardStats.SQRT_2_OVER_PI
							*Math.sqrt((double)occ/((double)beff*unseen))
							/Math.log((double)beff/unseen);
						final double w=Math.pow(1.0/err, AbstractCardStats.HC_INFO_POWER);
						hcSumW+=w;
						hcSumWLogE+=w*Math.log(est);
					}
				}
			}

			// --- Depth-0: MSB conditional model ---
			if(t<64 && nlzCount[t]>=8){
				final int beff=nlzCount[t];
				final int unseen=beff-bothMSB[t];
				if(unseen>=1 && unseen<beff){
					final double frac=(double)unseen/beff;
					final double q=frac/(2.0-frac);
					if(q>0 && q<1){
						final double est=scale*(-Math.log(q));
						if(est>0 && !Double.isNaN(est)){
							final double logInvQ=-Math.log(q);
							final double err=AbstractCardStats.SQRT_2_OVER_PI
								*Math.sqrt((1.0-q)/((double)beff*q))
								/logInvQ;
							final double w=Math.pow(1.0/err, AbstractCardStats.HC_INFO_POWER);
							hcSumW+=w;
							hcSumWLogE+=w*Math.log(est);
						}
					}
				}
			}
		}
		final double hcRaw=(hcSumW>0 ? Math.exp(hcSumWLogE/hcSumW) : 0);
		return (hcRaw>0 ? hcRaw*CorrectionFactor.hcCfFormula(hcRaw)*AbstractCardStats.HC_SCALE : 0);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Merge/Misc          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void add(CardinalityTracker log){
		throw new UnsupportedOperationException("TTLL3 merge not yet implemented");
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public int getGlobalExp(){return globalExp;}

	@Override
	public int actualBuckets(){return numBuckets;}

	public double occupancy(){
		int filled=0;
		for(int i=0; i<numBuckets; i++){
			if(getReg(i)!=0){filled++;}
		}
		return (double)filled/numBuckets;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] regs;
	private final int numBuckets;
	private int globalExp=0;
	private int numFloorBuckets;
	private long eeMask=-1L;
	private CardStats lastSummarized;
	int[] lastRawNlz;

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	/** Asymptotic meanRaw/trueCard ratio, measured 128k DDLs maxmult=8192 tmcf=1 (Apr 25 2026). */
	@Override public float terminalMeanCF(){return 1.387272f;}
	/** Asymptotic Mean+H ratio (6-bit history, 63-state SBS, entry/lin), measured 128k DDLs maxmult=8192 tmpcf=1 (Apr 25 2026). */
	@Override public float terminalMeanPlusCF(){return 1.051190f;}

	/** HLDLC blend: 0.70*LDLC + 0.30*Hybrid+2. Swept 0.36-0.80, 32k DDLs, hcweight=0.65 (Apr 25 2026). */
	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.70f;}

	/** Default HC/DLC blend weight inside LDLC. Swept 0.30-0.80, optimized for LDLC (Apr 25 2026). */
	static final double DEFAULT_HC_WEIGHT=0.65;

	/** MUST appear after all static field declarations. */
	static {
		StateTable.USE_TTLL3_HSB=true;
		if(!AbstractCardStats.HC_SCALE_EXPLICIT){AbstractCardStats.HC_SCALE=0.996029f;}
		if(OVERRIDE_LDLC_HC_WEIGHT<0){LDLC_HC_WEIGHT=DEFAULT_HC_WEIGHT;}
		CorrectionFactor.sbsFile="?sbsTTLL3_1536.tsv.gz";
		CorrectionFactor.SBS_CF_TABLE=null;
		CorrectionFactor.loadSbsTable();
	}
}
