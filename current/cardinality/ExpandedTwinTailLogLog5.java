package cardinality;

import shared.Tools;

/**
 * ETTLL5: Expanded-Tier TwinTailLogLog with 5-bit tails.
 *
 * Register layout (15 bits, 4 per 64-bit long):
 *   [14:10] = 5-bit localExp (expanded tier, sqrt(2) spacing)
 *   [9:5]   = 5-bit tail 1 history
 *   [4:0]   = 5-bit tail 0 history
 *
 * Tier geometry: expanded (sqrt(2) ratio). expandedNLZ = 2*rawNLZ + mantissa.
 * TIER_SCALE = 0.5 (each tier covers half an NLZ).
 *
 * Condensed state space: 3^5 = 243 classes (per-depth symmetry across 2 tails).
 * Estimation: VWMean (variance-weighted linear mean) with per-tier condensed
 * state table, falling back to harmonic mean when variance data unavailable.
 *
 * @author Nowi, Chloe, Brian Bushnell
 * @date April 2026
 */
public final class ExpandedTwinTailLogLog5 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int HIST_LEN=5;
	static final int NUM_TAILS=2;
	static final int EXP_BITS=5;

	static final int TAIL_BITS=NUM_TAILS*HIST_LEN; // 10
	static final int EXP_SHIFT=TAIL_BITS; // 10
	static final int REG_BITS=EXP_BITS+TAIL_BITS; // 15
	static final int REGS_PER_WORD=64/REG_BITS; // 4
	static final long REG_MASK_L=(1L<<REG_BITS)-1;
	static final int EXP_MASK=(1<<EXP_BITS)-1; // 31
	static final int TAIL_MASK=(1<<HIST_LEN)-1; // 31
	static final int TAILS_MASK=(1<<TAIL_BITS)-1; // 1023

	static final int TIER_NUMER=2;
	static final int TIER_DENOM=1;
	static final double TIER_SCALE=(double)TIER_DENOM/TIER_NUMER; // 0.5

	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	static final int NUM_CONDENSED=243; // 3^5

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ExpandedTwinTailLogLog5(){this(defaultBuckets(), 31, -1, 0);}

	ExpandedTwinTailLogLog5(parse.Parser p){
		super(p);
		numBuckets=p.loglogbuckets;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	ExpandedTwinTailLogLog5(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		numBuckets=buckets_;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	@Override
	public ExpandedTwinTailLogLog5 copy(){return new ExpandedTwinTailLogLog5(numBuckets, k, -1, minProb);}

	static int defaultBuckets(){return 256*REGS_PER_WORD;} // 1024

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
		final int expandedNLZ=TIER_NUMER*rawNLZ+mantissa;
		final int bucket=(int)(Long.remainderUnsigned(key, numBuckets));
		final int tailIdx=(int)(Long.divideUnsigned(key, numBuckets)&(NUM_TAILS-1));

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int reg=getReg(bucket);
		final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
		final int absCompTier=globalExp+localExp;
		final int delta=expandedNLZ-absCompTier;

		if(delta<-(HIST_LEN-1)){return;}

		if(delta>0){
			final int newLocalExp=Math.min(localExp+delta, EXP_MASK);
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
				int le=(reg>>>EXP_SHIFT)&EXP_MASK;
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
			final int minRawNLZ=(TIER_DENOM*targetTier)/TIER_NUMER;
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
		r[MEAN16_RAW_IDX]=vwMeanEstimate();
		return r;
	}

	public static final int MEAN16_RAW_IDX=AbstractCardStats.HC_IDX+1;

	public double[] ldlcEstimate(){
		final CardStats s=(lastSummarized!=null) ? lastSummarized : summarize();
		lastSummarized=null;
		return new double[]{s.ldlc(), s.dlcRaw(), s.hc(), s.lcMin(),
				0, s.hllRaw(), s.meanHistCF(), s.hybridPlus2(),
				vwMeanEstimate(), dualBucketLC()};
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
			final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
			final int absNlz=globalExp+localExp;

			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				filledCount++;
				final int combinedH=reg&TAILS_MASK;
				packedBuckets[i]=(char)(((absNlz+1)<<EXP_SHIFT)|combinedH);
			}
		}
		nlzCounts[0]=numBuckets-filledCount;

		final CardStats cs=new CardStats(packedBuckets, nlzCounts, 0, EXP_SHIFT, 0, 0,
				numBuckets, microIndex, added, null, 0, 0.0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), TIER_SCALE, 1);
		cs.overrideHC(hcETTLL(regs, globalExp, numBuckets));
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Condensed State Class     ----------------*/
	/*--------------------------------------------------------------*/

	/** Maps combined 10-bit (h1:h0) to 3^5 condensed class.
	 *  Per-depth pair: h0_bit + h1_bit in {0,1,2}, encoded base-3. */
	static int condensedClass(int combinedH){
		int idx=0, pow3=1;
		for(int d=0; d<HIST_LEN; d++){
			final int b0=(combinedH>>>d)&1;
			final int b1=(combinedH>>>(d+HIST_LEN))&1;
			idx+=(b0+b1)*pow3;
			pow3*=3;
		}
		return idx;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Dual-Bucket LC/DLC       ----------------*/
	/*--------------------------------------------------------------*/

	public double dualBucketLC(){
		int virtualTotal=0, virtualEmpty=0;
		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0 && globalExp==0){
				virtualTotal+=NUM_TAILS;
				virtualEmpty+=NUM_TAILS;
				continue;
			}
			final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
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
	/*----------------       VWMean Estimate        ----------------*/
	/*--------------------------------------------------------------*/

	/** Variance-weighted linear mean using per-tier condensed state table. */
	double vwMeanEstimateRaw(){
		if(PER_TIER_STATE_TABLE==null){return 0;}
		double sumWV=0, sumW=0;
		int used=0;
		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){continue;}
			used++;
			final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
			final int absNlz=globalExp+localExp;
			final int ch=reg&TAILS_MASK;
			final int tier=Math.min(absNlz, PER_TIER_MAX_TIER);
			final int idx=condensedClass(ch);
			double v=PER_TIER_STATE_TABLE[tier][idx];
			if(v<=0){continue;}
			if(absNlz>PER_TIER_MAX_TIER){
				v*=Math.pow(TIER_GROWTH_RATIO, absNlz-PER_TIER_MAX_TIER);
			}
			double w=1.0;
			if(PER_TIER_VARIANCE_TABLE!=null && tier<PER_TIER_VARIANCE_TABLE.length){
				final double variance=PER_TIER_VARIANCE_TABLE[tier][idx];
				if(variance>0){w=1.0/Math.pow(variance, VW_POWER);}
			}
			sumWV+=w*v;
			sumW+=w;
		}
		if(sumW==0){return 0;}
		return used*sumWV/sumW;
	}

	/** Harmonic mean fallback when variance table unavailable. */
	double harmonicMeanEstimateRaw(){
		if(PER_TIER_STATE_TABLE==null){return 0;}
		double sumInv=0;
		int filled=0;
		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
			final int absNlz=globalExp+localExp;
			final int ch=reg&TAILS_MASK;
			final int tier=Math.min(absNlz, PER_TIER_MAX_TIER);
			final int idx=condensedClass(ch);
			double v=PER_TIER_STATE_TABLE[tier][idx];
			if(v<=0){continue;}
			if(absNlz>PER_TIER_MAX_TIER){
				v*=Math.pow(TIER_GROWTH_RATIO, absNlz-PER_TIER_MAX_TIER);
			}
			sumInv+=1.0/v; filled++;
		}
		if(filled==0){return 0;}
		return (double)filled*filled/sumInv;
	}

	public double vwMeanEstimate(){
		double est=(PER_TIER_VARIANCE_TABLE!=null) ? vwMeanEstimateRaw() : harmonicMeanEstimateRaw();
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

	static double hcETTLL(long[] regs, int globalExp, int numBuckets){
		final int maxTier=64;
		final int[] nlzCount=new int[maxTier];
		final int[] allMSB=new int[maxTier];
		final int[][] depthSet=new int[HIST_LEN][maxTier];

		for(int i=0; i<numBuckets; i++){
			final int reg=getRegStatic(regs, i);
			if(reg==0 && globalExp==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
			final int absNlz=globalExp+localExp;
			if(absNlz<0 || absNlz>=maxTier){continue;}

			nlzCount[absNlz]++;

			boolean allMsbSet=true;
			for(int ti=0; ti<NUM_TAILS; ti++){
				final int tail=(reg>>>(ti*HIST_LEN))&TAIL_MASK;
				for(int d=0; d<HIST_LEN; d++){
					if(((tail>>>(HIST_LEN-1-d))&1)!=0){
						depthSet[d][absNlz]++;
					}
				}
				if(((tail>>>(HIST_LEN-1))&1)==0){allMsbSet=false;}
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

			for(int d=1; d<HIST_LEN; d++){
				final int srcTier=t+d;
				if(srcTier>=maxTier || nlzCount[srcTier]<=0){continue;}
				final int beff=NUM_TAILS*nlzCount[srcTier];
				final int unseen=beff-depthSet[d][srcTier];
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
		throw new UnsupportedOperationException("ETTLL5 merge not yet implemented");
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

	public CardStats consumeLastSummarized(){
		final CardStats cs=lastSummarized;
		lastSummarized=null;
		return cs;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Per-Tier State Table      ----------------*/
	/*--------------------------------------------------------------*/

	static double[][] PER_TIER_STATE_TABLE;
	static double[][] PER_TIER_VARIANCE_TABLE;
	static int PER_TIER_MAX_TIER;
	static double VW_POWER=Double.parseDouble(System.getProperty("vw.power", "0.4"));
	static double TIER_GROWTH_RATIO=Math.sqrt(2.0);

	static String PTIER_SECTION=System.getProperty("ptier.section", "all_lin");

	static void loadPerTierTable(){
		final String fname=dna.Data.findPath(perTierFile());
		if(fname==null){return;}
		fileIO.ByteFile bf=fileIO.ByteFile.makeByteFile(fname, false);
		final int maxTiers=(1<<EXP_BITS)*2;
		final double[][] rawTable=new double[maxTiers][1<<TAIL_BITS];
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
			if(tier>=maxTiers){continue;}
			for(int s=1; s<parts.length && (s-1)<(1<<TAIL_BITS); s++){
				try{rawTable[tier][s-1]=Double.parseDouble(parts[s]);}catch(NumberFormatException e){}
			}
			if(tier>maxTier){maxTier=tier;}
		}
		bf.close();

		final double[][] condensed=new double[maxTier+1][NUM_CONDENSED];
		final int[][] counts=new int[maxTier+1][NUM_CONDENSED];
		for(int tier=0; tier<=maxTier; tier++){
			for(int s=0; s<(1<<TAIL_BITS); s++){
				if(rawTable[tier][s]>0){
					final int cls=condensedClass(s);
					condensed[tier][cls]+=rawTable[tier][s];
					counts[tier][cls]++;
				}
			}
			for(int c=0; c<NUM_CONDENSED; c++){
				if(counts[tier][c]>0){condensed[tier][c]/=counts[tier][c];}
			}
		}
		PER_TIER_STATE_TABLE=condensed;
		PER_TIER_MAX_TIER=maxTier;

		if(maxTier>=2){
			double sumA=0, sumB=0;
			int n=0;
			final int lo=Math.max(0, maxTier-6);
			for(int tier=lo+1; tier<=maxTier; tier++){
				for(int c=0; c<NUM_CONDENSED; c++){
					if(condensed[tier][c]>0 && condensed[tier-1][c]>0){
						sumB+=condensed[tier][c];
						sumA+=condensed[tier-1][c];
						n++;
					}
				}
			}
			TIER_GROWTH_RATIO=(n>0 && sumA>0) ? sumB/sumA : Math.sqrt(2.0);
		}

		loadVarianceTable(fname);
	}

	static void loadVarianceTable(String fname){
		String varSection;
		if(PTIER_SECTION.startsWith("all_")){varSection="var_all";}
		else if(PTIER_SECTION.startsWith("entry_")){varSection="var_entry";}
		else if(PTIER_SECTION.startsWith("entryexit_")){varSection="var_entryexit";}
		else if(PTIER_SECTION.startsWith("pctile_")){varSection="var_pctile";}
		else{varSection="var_all";}

		fileIO.ByteFile bf=fileIO.ByteFile.makeByteFile(fname, false);
		final int maxTiers=(1<<EXP_BITS)*2;
		final double[][] varTable=new double[maxTiers][NUM_CONDENSED];
		int maxTier=0;
		boolean inSection=false;
		final String sectionHeader="#"+varSection;
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
			if(tier>=maxTiers){continue;}
			for(int s=1; s<parts.length && (s-1)<NUM_CONDENSED; s++){
				try{varTable[tier][s-1]=Double.parseDouble(parts[s]);}catch(NumberFormatException e){}
			}
			if(tier>maxTier){maxTier=tier;}
		}
		bf.close();
		if(maxTier>0){PER_TIER_VARIANCE_TABLE=varTable;}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	public static final String SBS_FILE="?sbsETTLL5.tsv.gz";

	static String perTierFile(){return "?perTierStateETTLL5.tsv.gz";}

	@Override public float terminalMeanCF(){return 1.0f;}
	@Override public float terminalMeanPlusCF(){return 1.0f;}
	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.68f;}

	static {
		loadPerTierTable();
	}
}
