package cardinality;

import shared.Tools;

/**
 * TTLL4: TwinTailLogLog with 4-bit tails — ceiling-based dual-history
 * cardinality estimator with 12-bit registers packed 5 per 64-bit long.
 *
 * Each register is 12 bits:
 *   [11:8] = 4-bit localExp (offset from globalExp)
 *   [7:4]  = h1: 4-bit tail for histBit==1 (MSB=tier, bit2=tier-1, bit1=tier-2, LSB=tier-3)
 *   [3:0]  = h0: 4-bit tail for histBit==0 (same semantics)
 *
 * Packing: 5 registers per long (60 bits used, 4 unused).
 * Default: 1280 buckets = 256 longs = 2KB.
 *
 * Bucket selection: modulo (non-power-of-2).
 *   bucket = Long.remainderUnsigned(key, numBuckets)
 *   histBit = Long.divideUnsigned(key, numBuckets) &amp; 1
 *
 * Update rule (symmetric only):
 *   delta = hashNLZ - (globalExp + localExp)
 *   delta &lt; -3 : ignore
 *   delta in [-3, 0]: set bit (3+delta) of history[histBit]
 *   delta &gt; 0 : advance localExp; shift both tails right by min(delta,4);
 *               set MSB of history[histBit]
 *   Overflow (newLocalExp &gt; 15): cap at 15.
 *
 * combined_h = (h1 &lt;&lt; 4) | h0 = 8-bit, 256 possible states.
 *
 * @author Brian, Chloe
 * @date April 26, 2026
 */
public final class TwinTailLogLog4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int NUM_TIERS=16;
	static final int TAIL_BITS=4;
	static final int TAIL_MASK=0xF;
	static final int EXP_SHIFT=2*TAIL_BITS; // 8
	static final int TAILS_MASK=(1<<EXP_SHIFT)-1; // 0xFF
	static final int REG_BITS=EXP_SHIFT+4; // 12
	static final long REG_MASK_L=(1L<<REG_BITS)-1; // 0xFFF
	static final int NUM_COMBINED=1<<EXP_SHIFT; // 256
	static final int REGS_PER_WORD=5; // 5×12=60 bits per long

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	TwinTailLogLog4(){this(1280, 31, -1, 0);}

	TwinTailLogLog4(parse.Parser p){
		super(p);
		numBuckets=p.loglogbuckets;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	TwinTailLogLog4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		numBuckets=buckets_;
		regs=new long[(numBuckets+REGS_PER_WORD-1)/REGS_PER_WORD];
		numFloorBuckets=numBuckets;
	}

	@Override
	public TwinTailLogLog4 copy(){return new TwinTailLogLog4(numBuckets, k, -1, minProb);}

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
		cs.overrideHC(hcTTLL4(regs, globalExp, numBuckets));
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
				virtualTotal+=2;
				virtualEmpty+=2;
				continue;
			}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;

			for(int t=0; t<2; t++){
				final int tail=(t==0) ? (reg&TAIL_MASK) : ((reg>>>TAIL_BITS)&TAIL_MASK);
				if((tail&8)!=0){
					virtualTotal++;
				}else if((tail&4)!=0){
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
	/*----------------   256-State Mean Estimate    ----------------*/
	/*--------------------------------------------------------------*/

	public double ttllMeanEstimate(){
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
		double est=(double)filled*filled/sumInv;
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
	/*----------------     HC for 4-bit tails       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * HC estimate using 4-bit tail observations.
	 * Four observation types per tail:
	 *   MSB (bit 3): conditional probability model (same as TTLL)
	 *   Bit 2: unconditional half-population LC at depth 1
	 *   Bit 1: unconditional half-population LC at depth 2
	 *   LSB (bit 0): unconditional half-population LC at depth 3
	 */
	static double hcTTLL4(long[] regs, int globalExp, int numBuckets){
		final int[] nlzCount=new int[64];
		final int[] bothMSB=new int[64];
		final int[] h0Bit2Set=new int[64];
		final int[] h1Bit2Set=new int[64];
		final int[] h0Bit1Set=new int[64];
		final int[] h1Bit1Set=new int[64];
		final int[] h0Bit0Set=new int[64];
		final int[] h1Bit0Set=new int[64];

		for(int i=0; i<numBuckets; i++){
			final int reg=getRegStatic(regs, i);
			if(reg==0 && globalExp==0){continue;}
			final int localExp=(reg>>>EXP_SHIFT)&0xF;
			final int absNlz=globalExp+localExp;
			if(absNlz<0 || absNlz>=64){continue;}

			final int h0=reg&TAIL_MASK;
			final int h1=(reg>>>TAIL_BITS)&TAIL_MASK;

			nlzCount[absNlz]++;
			if(((h0>>>3)&1)!=0 && ((h1>>>3)&1)!=0){bothMSB[absNlz]++;}
			if(((h0>>>2)&1)!=0){h0Bit2Set[absNlz]++;}
			if(((h1>>>2)&1)!=0){h1Bit2Set[absNlz]++;}
			if(((h0>>>1)&1)!=0){h0Bit1Set[absNlz]++;}
			if(((h1>>>1)&1)!=0){h1Bit1Set[absNlz]++;}
			if((h0&1)!=0){h0Bit0Set[absNlz]++;}
			if((h1&1)!=0){h1Bit0Set[absNlz]++;}
		}

		int maxNlz=0;
		for(int t=63; t>=0; t--){
			if(nlzCount[t]>0){maxNlz=t; break;}
		}
		final int maxTierHC=Math.min(maxNlz+TAIL_BITS, 63);

		double hcSumW=0, hcSumWLogE=0;
		for(int t=0; t<=maxTierHC; t++){
			final double scale=2.0*Math.pow(2.0, t+1)*(double)numBuckets;

			// --- Depth-1: bit 2 (h0_bit2 + h1_bit2 pooled) ---
			if(t+1<64 && nlzCount[t+1]>0){
				final int beff=2*nlzCount[t+1];
				final int unseen=(nlzCount[t+1]-h0Bit2Set[t+1])
					+(nlzCount[t+1]-h1Bit2Set[t+1]);
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

			// --- Depth-2: bit 1 (h0_bit1 + h1_bit1 pooled) ---
			if(t+2<64 && nlzCount[t+2]>0){
				final int beff=2*nlzCount[t+2];
				final int unseen=(nlzCount[t+2]-h0Bit1Set[t+2])
					+(nlzCount[t+2]-h1Bit1Set[t+2]);
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

			// --- Depth-3: bit 0 / LSB (h0_bit0 + h1_bit0 pooled) ---
			if(t+3<64 && nlzCount[t+3]>0){
				final int beff=2*nlzCount[t+3];
				final int unseen=(nlzCount[t+3]-h0Bit0Set[t+3])
					+(nlzCount[t+3]-h1Bit0Set[t+3]);
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
		throw new UnsupportedOperationException("TTLL4 merge not yet implemented");
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
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	/** Placeholder — needs measurement with 128k DDLs maxmult=8192 tmcf=1. */
	@Override public float terminalMeanCF(){return 1.387f;}
	/** Placeholder — needs measurement with 128k DDLs maxmult=8192 tmpcf=1. */
	@Override public float terminalMeanPlusCF(){return 1.05f;}

	@Override public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.70f;}

	static final double DEFAULT_HC_WEIGHT=0.65;

	/*--------------------------------------------------------------*/
	/*----------------   Per-Tier State Table Est.   ----------------*/
	/*--------------------------------------------------------------*/

	static double[][] PER_TIER_STATE_TABLE;
	static int PER_TIER_MAX_TIER;

	static String PTIER_SECTION=System.getProperty("ptier.section", "all_lin");

	static void loadPerTierTable(){
		final String fname=dna.Data.findPath("?perTierStateTTLL4.tsv.gz");
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

	static {
		StateTable.USE_TTLL4_HSB=true;
		if(!AbstractCardStats.HC_SCALE_EXPLICIT){AbstractCardStats.HC_SCALE=0.996f;}
		if(OVERRIDE_LDLC_HC_WEIGHT<0){LDLC_HC_WEIGHT=DEFAULT_HC_WEIGHT;}
		loadPerTierTable();
	}
}
