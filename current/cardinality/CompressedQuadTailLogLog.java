package cardinality;

import shared.Tools;

/**
 * CQuadTLL: Compressed-Tier QuadTailLogLog.
 * Fork of CTTLL with 4 tails and 2x NLZ compression.
 *
 * Register layout (12 bits per bucket, 5 per 64-bit word):
 *   [11:8] = 4-bit localExp (in compressed tier space)
 *   [7:0]  = 4 tails, each 2 bits
 *
 * Tier compression: compressedTier = rawNLZ / 2.
 * TIER_SCALE = 2.0 (average rawNLZ per compressed tier).
 *
 * @author Brian, Nahida
 * @date June 2026
 */
public final class CompressedQuadTailLogLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------      Variable Configuration     ----------------*/
	/*--------------------------------------------------------------*/

	static int HIST_LEN=2;
	static int NUM_TAILS=4;

	static int EXP_BITS=4;
	static int TAIL_BITS=NUM_TAILS*HIST_LEN;
	static int EXP_SHIFT=TAIL_BITS;
	static int REG_BITS=EXP_BITS+TAIL_BITS;
	static int REGS_PER_WORD=64/REG_BITS;
	static long REG_MASK_L=(1L<<REG_BITS)-1;
	static int EXP_MASK=(1<<EXP_BITS)-1;
	static int TAIL_MASK=(1<<HIST_LEN)-1;
	static int TAILS_MASK=(1<<TAIL_BITS)-1;
	static int NUM_COMBINED=1<<TAIL_BITS;

	static int TIER_NUMER=1;
	static int TIER_DENOM=2;
	static double TIER_SCALE=(double)TIER_DENOM/TIER_NUMER;

	static void reconfigure(){
		TAIL_BITS=NUM_TAILS*HIST_LEN;
		EXP_SHIFT=TAIL_BITS;
		REG_BITS=EXP_BITS+TAIL_BITS;
		REGS_PER_WORD=64/REG_BITS;
		REG_MASK_L=(1L<<REG_BITS)-1;
		EXP_MASK=(1<<EXP_BITS)-1;
		TAIL_MASK=(1<<HIST_LEN)-1;
		TAILS_MASK=(1<<TAIL_BITS)-1;
		NUM_COMBINED=1<<TAIL_BITS;
		TIER_SCALE=(double)TIER_DENOM/TIER_NUMER;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	CompressedQuadTailLogLog(){this(defaultBuckets(), 31, -1, 0);}

	CompressedQuadTailLogLog(parse.Parser p){
		super(p);
		numBuckets=p.loglogbuckets;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	CompressedQuadTailLogLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		numBuckets=buckets_;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	@Override
	public CompressedQuadTailLogLog copy(){return new CompressedQuadTailLogLog(numBuckets, k, -1, minProb);}

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
		final int compressedNLZ=rawNLZ/2;
		final int bucket=(int)(Long.remainderUnsigned(key, numBuckets));
		final int tailIdx=(int)(Long.divideUnsigned(key, numBuckets)&(NUM_TAILS-1));

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int reg=getReg(bucket);
		final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
		final int absCompTier=globalExp+localExp;
		final int delta=compressedNLZ-absCompTier;

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
		r[MEAN16_RAW_IDX]=ttllMeanEstimate();
		return r;
	}

	public static final int MEAN16_RAW_IDX=AbstractCardStats.HC_IDX+1;

	public double[] ldlcEstimate(){
		final CardStats s=(lastSummarized!=null) ? lastSummarized : summarize();
		lastSummarized=null;
		final double blc=bloomLCEstimate();
		return new double[]{s.ldlc(), s.dlcRaw(), s.hc(), s.lcMin(),
				0, s.hllRaw(), s.meanHistCF(), s.hybridPlus2(),
				ttllMeanEstimate(), (NUM_TAILS>=4 ? blc : dualBucketLC())};
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
		lastRawNlz=nlzCounts;

		final CardStats cs=new CardStats(packedBuckets, nlzCounts, 0, EXP_SHIFT, 0, 0,
				numBuckets, microIndex, added, null, 0, 0.0,
				Integer.MAX_VALUE, null, terminalMeanCF(), terminalMeanPlusCF(), TIER_SCALE, 1);
		double hcEst=hcCTTLL(regs, globalExp, numBuckets);
		if(CorrectionFactor.USE_CORRECTION && CorrectionFactor.cfTable!=null
				&& CorrectionFactor.HC<CorrectionFactor.cfTable.length
				&& CorrectionFactor.cfTable[CorrectionFactor.HC]!=null
				&& CorrectionFactor.cfKeys!=null){
			final double ks=(CorrectionFactor.cfTableBuckets>0 && numBuckets!=CorrectionFactor.cfTableBuckets)
				? (double)CorrectionFactor.cfTableBuckets/numBuckets : 1.0;
			hcEst*=CorrectionFactor.getCF(cs.dlcRaw(), hcEst,
				CorrectionFactor.cfTable[CorrectionFactor.HC],
				CorrectionFactor.cfKeys, 5, 0.0001, ks);
		}
		cs.overrideHC(hcEst);
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
	/*----------------   Per-Tier State Mean Estimate   ----------------*/
	/*--------------------------------------------------------------*/

	double ttllMeanEstimateRaw(){
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
			final int idx=USE_CONDENSED ? condensedClass(ch) : ch;
			final double v=PER_TIER_STATE_TABLE[tier][idx];
			if(v>0){sumInv+=1.0/v; filled++;}
		}
		if(filled==0){return 0;}
		return (double)filled*filled/sumInv;
	}

	public double ttllMeanEstimate(){
		double est=ttllMeanEstimateRaw();
		if(est<=0){return 0;}
		if(CorrectionFactor.USE_CORRECTION && CorrectionFactor.cfTable!=null
				&& CorrectionFactor.MEAN16<CorrectionFactor.cfTable.length
				&& CorrectionFactor.cfTable[CorrectionFactor.MEAN16]!=null){
			final double keyScale=(CorrectionFactor.cfTableBuckets>0 && numBuckets!=CorrectionFactor.cfTableBuckets)
				? (double)CorrectionFactor.cfTableBuckets/numBuckets : 1.0;
			final double cf=CorrectionFactor.getCF(est, est,
				CorrectionFactor.cfTable[CorrectionFactor.MEAN16],
				CorrectionFactor.cfKeys, 5, 0.0001, keyScale);
			est*=cf;
		}
		return est;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Bloom Filter LC Estimate    ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-bucket Bloom filter LC: treat each depth row of NUM_TAILS bits
	 *  as a tiny Bloom filter and do LC on each row independently. */
	public double bloomLCEstimate(){
		if(NUM_TAILS<4){return 0;}
		final double tierRatio=Math.pow(2.0, TIER_SCALE);
		double sumInv=0;
		int filled=0;

		for(int i=0; i<numBuckets; i++){
			final int reg=getReg(i);
			if(reg==0 && globalExp==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&EXP_MASK;
			final int absNlz=globalExp+localExp;

			double bucketEst=0;
			double bucketW=0;

			for(int d=0; d<HIST_LEN; d++){
				int setCount=0;
				for(int ti=0; ti<NUM_TAILS; ti++){
					final int tail=(reg>>>(ti*HIST_LEN))&TAIL_MASK;
					if(((tail>>>(HIST_LEN-1-d))&1)!=0){setCount++;}
				}
				final int srcTier=absNlz+d;
				final double scale=2.0*Math.pow(2.0, (srcTier+1)*TIER_SCALE)
					/(tierRatio-1.0)*(double)numBuckets;

				// MSB row: one tail is forced set by the routing, so effective n=NUM_TAILS-1
				final int n=(d==0) ? (NUM_TAILS-1) : NUM_TAILS;
				final int effectiveSet=(d==0) ? Math.max(setCount-1, 0) : setCount;
				final int unseen=n-effectiveSet;

				if(n<=0 || effectiveSet<=0){continue;}
				final double v=(unseen<=0) ? 0.5 : (double)unseen;
				final double lc=Math.log((double)n/v);
				final double est=scale*lc;
				if(est>0){
					final double occ=n-v;
					final double w=(occ>0 && v>0) ? occ/((double)n*v) : 0.01;
					bucketEst+=est*w;
					bucketW+=w;
				}
			}

			if(bucketW>0){
				final double perBucket=bucketEst/bucketW;
				if(perBucket>0){sumInv+=1.0/perBucket; filled++;}
			}
		}
		if(filled==0){return 0;}
		return (double)filled*filled/sumInv;
	}

	/*--------------------------------------------------------------*/
	/*----------------     HC — Variable Depth       ----------------*/
	/*--------------------------------------------------------------*/

	static double hcCTTLL(long[] regs, int globalExp, int numBuckets){
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

			if(NUM_TAILS>2){
				// Pool all depths (including MSB) across all tails as independent LC
				for(int d=0; d<HIST_LEN; d++){
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
			}else{
				// 2-tail: depth 1+ pooled LC, depth 0 conditional model
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
		}
		final double hcRaw=(hcSumW>0 ? Math.exp(hcSumWLogE/hcSumW) : 0);
		return (hcRaw>0 ? hcRaw*CorrectionFactor.hcCfFormula(hcRaw)*AbstractCardStats.HC_SCALE : 0);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Merge/Misc          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void add(CardinalityTracker log){
		throw new UnsupportedOperationException("CQuadTLL merge not yet implemented");
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
	static double[][] PER_TIER_FULL_TABLE;
	static int PER_TIER_MAX_TIER;
	static boolean USE_CONDENSED=false;
	static final int NUM_CLASSES=25;

	static int condensedClass(int combinedH){
		return Integer.bitCount(combinedH&0xAA)*5+Integer.bitCount(combinedH&0x55);
	}

	static String PTIER_SECTION=System.getProperty("ptier.section", "all_lin");

	static void loadPerTierTable(){
		final String fname=dna.Data.findPath(perTierFile());
		if(fname==null){return;}
		fileIO.ByteFile bf=fileIO.ByteFile.makeByteFile(fname, false);
		final int maxTiers=(1<<EXP_BITS)*2;
		final double[][] table=new double[maxTiers][NUM_COMBINED];
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
			for(int s=1; s<parts.length && (s-1)<NUM_COMBINED; s++){
				try{table[tier][s-1]=Double.parseDouble(parts[s]);}catch(NumberFormatException e){}
			}
			if(tier>maxTier){maxTier=tier;}
		}
		bf.close();

		if(NUM_TAILS>=4){
			final double[][] condensed=new double[maxTier+1][NUM_CLASSES];
			final int[][] counts=new int[maxTier+1][NUM_CLASSES];
			for(int tier=0; tier<=maxTier; tier++){
				for(int s=0; s<NUM_COMBINED; s++){
					if(table[tier][s]>0){
						final int cls=condensedClass(s);
						condensed[tier][cls]+=table[tier][s];
						counts[tier][cls]++;
					}
				}
				for(int c=0; c<NUM_CLASSES; c++){
					if(counts[tier][c]>0){condensed[tier][c]/=counts[tier][c];}
				}
			}
			PER_TIER_STATE_TABLE=condensed;
			PER_TIER_FULL_TABLE=table;
			USE_CONDENSED=true;
		}else{
			PER_TIER_STATE_TABLE=table;
			PER_TIER_FULL_TABLE=table;
			USE_CONDENSED=false;
		}
		PER_TIER_MAX_TIER=maxTier;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	static final double DEFAULT_HC_WEIGHT=0.50;

	public static final String SBS_FILE="?sbsCQuadTLL_1280.tsv.gz";

	static String sbsFile(){
		return SBS_FILE;
	}

	static String perTierFile(){
		return "?perTierStateCQuadTLL.tsv.gz";
	}

	/** Terminal CFs — placeholder 1.0f until measured for 2x compression. */
	@Override public float terminalMeanCF(){return 1.0f;}
	@Override public float terminalMeanPlusCF(){return 1.0f;}
	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.68f;}

	/** Derive HSB tables from per-tier state table for Mean+H. */
	static void computeHSB(){
		if(PER_TIER_FULL_TABLE==null || PER_TIER_MAX_TIER<3){return;}

		final int highTier=PER_TIER_MAX_TIER;
		StateTable.CF_CTTLL_4=hsbRow(PER_TIER_FULL_TABLE[highTier]);
		StateTable.CF_CTTLL_4_TIERS=new double[3][NUM_COMBINED];
		for(int t=0; t<3 && t<=PER_TIER_MAX_TIER; t++){
			StateTable.CF_CTTLL_4_TIERS[t]=hsbRow(PER_TIER_FULL_TABLE[t]);
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
