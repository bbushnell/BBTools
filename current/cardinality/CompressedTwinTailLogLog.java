package cardinality;

import shared.Tools;

/**
 * CTTLL: Compressed-Tier TwinTailLogLog.
 * Fork of TTLL3 with compressed tiers and variable tail geometry.
 *
 * Register layout (variable width, packed in 64-bit longs):
 *   [REG_BITS-1 : TAIL_BITS] = 4-bit localExp (in compressed tier space)
 *   [TAIL_BITS-1 : 0] = NUM_TAILS tails, each HIST_LEN bits
 *
 * Tier compression: compressedTier = (2 * rawNLZ) / 3.
 * TIER_SCALE = 1.5 (average rawNLZ per compressed tier).
 *
 * Estimation: per-tier state table lookup + harmonic aggregation + CF correction.
 *
 * @author Brian, Chloe
 * @date April 2026
 */
public final class CompressedTwinTailLogLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------      Variable Configuration     ----------------*/
	/*--------------------------------------------------------------*/

	static int HIST_LEN=2;
	static int NUM_TAILS=2;

	static int TAIL_BITS=NUM_TAILS*HIST_LEN;
	static int EXP_SHIFT=TAIL_BITS;
	static int REG_BITS=4+TAIL_BITS;
	static int REGS_PER_WORD=64/REG_BITS;
	static long REG_MASK_L=(1L<<REG_BITS)-1;
	static int TAIL_MASK=(1<<HIST_LEN)-1;
	static int TAILS_MASK=(1<<TAIL_BITS)-1;
	static int NUM_COMBINED=1<<TAIL_BITS;

	static final int TIER_NUMER=2;
	static final int TIER_DENOM=3;
	static final double TIER_SCALE=(double)TIER_DENOM/TIER_NUMER;

	/** Mantissa threshold for log-uniform half-NLZ steps: (2-sqrt(2)) * 2^20. */
	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	static void reconfigure(){
		TAIL_BITS=NUM_TAILS*HIST_LEN;
		EXP_SHIFT=TAIL_BITS;
		REG_BITS=4+TAIL_BITS;
		REGS_PER_WORD=64/REG_BITS;
		REG_MASK_L=(1L<<REG_BITS)-1;
		TAIL_MASK=(1<<HIST_LEN)-1;
		TAILS_MASK=(1<<TAIL_BITS)-1;
		NUM_COMBINED=1<<TAIL_BITS;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	CompressedTwinTailLogLog(){this(defaultBuckets(), 31, -1, 0);}

	CompressedTwinTailLogLog(parse.Parser p){
		super(p);
		numBuckets=p.loglogbuckets;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	CompressedTwinTailLogLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		numBuckets=buckets_;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	@Override
	public CompressedTwinTailLogLog copy(){return new CompressedTwinTailLogLog(numBuckets, k, -1, minProb);}

	static int defaultBuckets(){return 256*REGS_PER_WORD;}

	/*--------------------------------------------------------------*/
	/*----------------       Register Access        ----------------*/
	/*--------------------------------------------------------------*/

	private int getReg(int bucket){
		return (int)((regs[bucket/REGS_PER_WORD]>>>((bucket%REGS_PER_WORD)*REG_BITS))&REG_MASK_L);
	}

	private void setReg(int bucket, int val){
		final int word=bucket/REGS_PER_WORD;
		final int shift=(bucket%REGS_PER_WORD)*REG_BITS;
		regs[word]=(regs[word]&~(REG_MASK_L<<shift))|(((long)val&REG_MASK_L)<<shift);
	}

	private static int getRegStatic(long[] regs, int bucket){
		return (int)((regs[bucket/REGS_PER_WORD]>>>((bucket%REGS_PER_WORD)*REG_BITS))&REG_MASK_L);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Core Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);
		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int rawNLZ=Long.numberOfLeadingZeros(key);
		final int mBits;
		if(rawNLZ>=43){
			mBits=(int)((key<<(rawNLZ-42))&0xFFFFF);
		}else{
			mBits=(int)((key>>>(42-rawNLZ))&0xFFFFF);
		}
		final int mantissa=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		final int compressedNLZ=(2*rawNLZ+mantissa)/3;
		final int bucket=(int)(Long.remainderUnsigned(key, numBuckets));
		final int tailIdx=(int)(Long.divideUnsigned(key, numBuckets)&(NUM_TAILS-1));

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int reg=getReg(bucket);
		final int localExp=(reg>>>EXP_SHIFT)&0xF;
		final int absCompTier=globalExp+localExp;
		final int delta=compressedNLZ-absCompTier;

		if(delta<-(HIST_LEN-1)){return;}

		if(delta>0){
			final int newLocalExp=Math.min(localExp+delta, 15);
			final boolean wasZeroExp=(localExp==0);
			final int shiftAmt=Math.min(delta, HIST_LEN);
			int combinedTails=0;
			for(int ti=0; ti<NUM_TAILS; ti++){
				int tail=((reg>>>(ti*HIST_LEN))&TAIL_MASK)>>>shiftAmt;
				if(ti==tailIdx){tail|=(1<<(HIST_LEN-1));}
				combinedTails|=(tail<<(ti*HIST_LEN));
			}
			lastCardinality=-1;
			lastSummarized=null;
			setReg(bucket, (newLocalExp<<EXP_SHIFT)|combinedTails);
			if(wasZeroExp){
				numFloorBuckets--;
				if(numFloorBuckets==0){advanceGlobal();}
			}
		}else{
			final int bitPos=(HIST_LEN-1)+delta;
			final int regBitPos=tailIdx*HIST_LEN+bitPos;
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
			updateEeMask();
		}
	}

	private void updateEeMask(){
		final int targetTier=globalExp-(HIST_LEN-1);
		if(targetTier<=0){
			eeMask=-1L;
		}else{
			// With mantissa: compressedTier = (2*rawNLZ+mantissa)/3.
			// Minimum rawNLZ where max compressed tier (mantissa=1) >= targetTier:
			// (2*rawNLZ+1)/3 >= targetTier → rawNLZ >= (3*targetTier-1)/2 = ceil = (3*targetTier)/2
			final int minRawNLZ=(3*targetTier)/2;
			eeMask=(minRawNLZ>=64) ? 0 : (-1L)>>>minRawNLZ;
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
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), TIER_SCALE, 1);
		cs.overrideHC(hcCTTLL(regs, globalExp, numBuckets));
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Dual-Bucket LC/DLC       ----------------*/
	/*--------------------------------------------------------------*/

	public double dualBucketLC(){
		int virtualTotal=0;
		int virtualEmpty=0;

		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0 && globalExp==0){
				virtualTotal+=NUM_TAILS;
				virtualEmpty+=NUM_TAILS;
				continue;
			}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;

			for(int ti=0; ti<NUM_TAILS; ti++){
				final int tail=(reg>>>(ti*HIST_LEN))&TAIL_MASK;
				if(tail!=0){
					virtualTotal++;
				}else if(absNlz<HIST_LEN){
					virtualTotal++;
					virtualEmpty++;
				}
			}
		}

		if(virtualTotal==0 || virtualEmpty==0){return 0;}
		return (double)virtualTotal*Math.log((double)virtualTotal/virtualEmpty);
	}

	/*--------------------------------------------------------------*/
	/*----------------   Per-Tier State Mean Estimate   ----------------*/
	/*--------------------------------------------------------------*/

	double ttllMeanEstimateRaw(){
		if(PER_TIER_STATE_TABLE==null){return 0;}
		double sumInv=0;
		int filled=0;
		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;
			final int ch=reg&TAILS_MASK;
			final int tier=Math.min(absNlz, PER_TIER_MAX_TIER);
			final double v=PER_TIER_STATE_TABLE[tier][ch];
			if(v>0){sumInv+=1.0/v; filled++;}
		}
		if(filled==0){return 0;}
		return (double)filled*filled/sumInv;
	}

	public double ttllMeanEstimate(){
		double est=ttllMeanEstimateRaw();
		if(est<=0){return 0;}
		if(CorrectionFactor.USE_CORRECTION && CorrectionFactor.v1Matrix!=null
				&& CorrectionFactor.MEAN16<CorrectionFactor.v1Matrix.length
				&& CorrectionFactor.v1Matrix[CorrectionFactor.MEAN16]!=null){
			final double keyScale=(CorrectionFactor.v1Buckets>0 && numBuckets!=CorrectionFactor.v1Buckets)
				? (double)CorrectionFactor.v1Buckets/numBuckets : 1.0;
			final double cf=CorrectionFactor.getCF(est, est,
				CorrectionFactor.v1Matrix[CorrectionFactor.MEAN16],
				CorrectionFactor.v1Keys, 5, 0.0001, keyScale);
			est*=cf;
		}
		return est;
	}

	/*--------------------------------------------------------------*/
	/*----------------     HC — Variable Depth       ----------------*/
	/*--------------------------------------------------------------*/

	static double hcCTTLL(long[] regs, int globalExp, int numBuckets){
		final int maxTier=64;
		final int[] nlzCount=new int[maxTier];
		final int[] allMSB=new int[maxTier];
		final int[][] depthSet=new int[HIST_LEN-1][maxTier];

		for(int i=0; i<numBuckets; i++){
			final int reg=getRegStatic(regs, i);
			if(reg==0 && globalExp==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;
			if(absNlz<0 || absNlz>=maxTier){continue;}

			nlzCount[absNlz]++;

			boolean allMsbSet=true;
			for(int ti=0; ti<NUM_TAILS; ti++){
				final int tail=(reg>>>(ti*HIST_LEN))&TAIL_MASK;
				if(((tail>>>(HIST_LEN-1))&1)==0){allMsbSet=false;}
				for(int d=1; d<HIST_LEN; d++){
					if(((tail>>>(HIST_LEN-1-d))&1)!=0){
						depthSet[d-1][absNlz]++;
					}
				}
			}
			if(allMsbSet){allMSB[absNlz]++;}
		}

		int maxNlz=0;
		for(int t=maxTier-1; t>=0; t--){
			if(nlzCount[t]>0){maxNlz=t; break;}
		}
		final int maxTierHC=Math.min(maxNlz+HIST_LEN+1, maxTier-1);
		final double tierRatio=Math.pow(2.0, TIER_SCALE);

		double hcSumW=0, hcSumWLogE=0;
		for(int t=0; t<=maxTierHC; t++){
			final double scale=2.0*Math.pow(2.0, (t+1)*TIER_SCALE)/(tierRatio-1.0)
				*(double)numBuckets;

			// Depth 1..HIST_LEN-1: unconditional LC, pooled across all tails
			for(int d=1; d<HIST_LEN; d++){
				final int srcTier=t+d;
				if(srcTier>=maxTier || nlzCount[srcTier]<=0){continue;}
				final int beff=NUM_TAILS*nlzCount[srcTier];
				final int unseen=NUM_TAILS*nlzCount[srcTier]-depthSet[d-1][srcTier];
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

			// Depth 0: MSB conditional model (all tails have MSB set)
			if(t<maxTier && nlzCount[t]>=8){
				final int beff=nlzCount[t];
				final int unseen=beff-allMSB[t];
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
		throw new UnsupportedOperationException("CTTLL merge not yet implemented");
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

	private final long[] regs;
	private final int numBuckets;
	private int globalExp=0;
	private int numFloorBuckets;
	private long eeMask=-1L;
	private CardStats lastSummarized;
	int[] lastRawNlz;

	/*--------------------------------------------------------------*/
	/*----------------     Per-Tier State Table      ----------------*/
	/*--------------------------------------------------------------*/

	static double[][] PER_TIER_STATE_TABLE;
	static int PER_TIER_MAX_TIER;

	static String PTIER_SECTION=System.getProperty("ptier.section", "all_lin");

	static void loadPerTierTable(){
		final String fname=dna.Data.findPath("?perTierStateCTTLL.tsv.gz");
		if(fname==null){return;}
		fileIO.ByteFile bf=fileIO.ByteFile.makeByteFile(fname, false);
		final double[][] table=new double[16][NUM_COMBINED];
		int maxTier=0;
		boolean inSection=false;
		final String sectionHeader="#"+PTIER_SECTION;
		for(byte[] raw=bf.nextLine(); raw!=null; raw=bf.nextLine()){
			if(raw.length==0){continue;}
			final String line=new String(raw).trim();
			if(line.startsWith(sectionHeader+"\t") || line.equals(sectionHeader)){
				inSection=true; continue;
			}
			if(inSection && line.startsWith("#Tier")){continue;}
			if(inSection && line.startsWith("#")){break;}
			if(!inSection){continue;}
			final String[] parts=line.split("\t");
			if(parts.length<2){continue;}
			final int tier;
			try{tier=Integer.parseInt(parts[0]);}catch(NumberFormatException e){continue;}
			if(tier>=16){continue;}
			for(int s=1; s<parts.length && (s-1)<NUM_COMBINED; s++){
				try{table[tier][s-1]=Double.parseDouble(parts[s]);}catch(NumberFormatException e){}
			}
			if(tier>maxTier){maxTier=tier;}
		}
		bf.close();
		PER_TIER_STATE_TABLE=table;
		PER_TIER_MAX_TIER=maxTier;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	static final double DEFAULT_HC_WEIGHT=0.50;

	@Override public float terminalMeanCF(){return 1.0f;}
	@Override public float terminalMeanPlusCF(){return 1.0f;}
	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.68f;}

	/** Derive HSB tables from per-tier state table for Mean+H. */
	static void computeHSB(){
		if(PER_TIER_STATE_TABLE==null || PER_TIER_MAX_TIER<3){return;}

		// Steady state from highest tier
		final int highTier=PER_TIER_MAX_TIER;
		StateTable.CF_CTTLL_4=hsbRow(PER_TIER_STATE_TABLE[highTier]);

		// Per-tier for first 3 tiers
		StateTable.CF_CTTLL_4_TIERS=new double[3][NUM_COMBINED];
		for(int t=0; t<3 && t<=PER_TIER_MAX_TIER; t++){
			StateTable.CF_CTTLL_4_TIERS[t]=hsbRow(PER_TIER_STATE_TABLE[t]);
		}
		StateTable.USE_CTTLL_HSB=true;
	}

	private static double[] hsbRow(double[] row){
		final double[] out=new double[NUM_COMBINED];
		double sum=0; int count=0;
		for(int s=0; s<NUM_COMBINED; s++){
			if(row[s]>0){sum+=row[s]; count++;}
		}
		if(count==0){return out;}
		final double avg=sum/count;
		final double LOG2=Math.log(2.0);
		for(int s=0; s<NUM_COMBINED; s++){
			if(row[s]>0 && avg>0){
				out[s]=Math.log(row[s]/avg)/LOG2;
			}
		}
		return out;
	}

	static {
		reconfigure();
		if(OVERRIDE_LDLC_HC_WEIGHT<0){LDLC_HC_WEIGHT=DEFAULT_HC_WEIGHT;}
		loadPerTierTable();
		computeHSB();
	}
}
