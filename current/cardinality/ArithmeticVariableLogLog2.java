package cardinality;


/**
 * AVLL2 — cleanroom copy of AVDLL64, one change: HTL 14→13 (MAX_EXP 15→18).
 * RADIX=56, BPW=11, unsigned 64-bit. Uses CardStats + CF tables for estimation.
 *
 * @author Brian Bushnell
 * @date May 2026
 */
public final class ArithmeticVariableLogLog2 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ArithmeticVariableLogLog2(){this(2048, 31, -1, 0);}

	ArithmeticVariableLogLog2(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ArithmeticVariableLogLog2(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*11<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*11-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*11;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=(modBuckets+10)/11;
		packed=new long[packedLen];
		floorCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+10)/11);}

	@Override public ArithmeticVariableLogLog2 copy(){return new ArithmeticVariableLogLog2(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}
	@Override public int bitsPerWord(){return 64;}
	@Override public int bucketsPerWord(){return 11;}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	int getReg(int i){
		final int wordIdx=i/11;
		final int pos=i%11;
		final long word=packed[wordIdx];
		return (int)(Long.remainderUnsigned(Long.divideUnsigned(word, POW[pos]), RADIX));
	}

	private void setReg(int i, int val){
		final int wordIdx=i/11;
		final int pos=i%11;
		long word=packed[wordIdx];
		final int oldVal=(int)(Long.remainderUnsigned(Long.divideUnsigned(word, POW[pos]), RADIX));
		word=word+(long)(val-oldVal)*POW[pos];
		packed[wordIdx]=word;
	}

	/*--------------------------------------------------------------*/
	/*----------------      State Mapping Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	static int stateToExp(int state){
		if(state==0){return 0;}
		if(state<=2){return 1;}
		if(state<HIST_STATES){return (state-3)/4+2;}
		return HIST_TIER_LIMIT+1+(state-HIST_STATES);
	}

	static int stateToHist(int state){
		if(state==0){return 0;}
		if(state<=2){return (state-1)<<1;}
		if(state<HIST_STATES){return (state-3)&3;}
		return 0;
	}

	static int toState(int exp, int hist){
		if(exp==0){return 0;}
		if(exp==1){return 1+(hist>>1);}
		if(exp<=HIST_TIER_LIMIT){return 4*exp-5+hist;}
		return HIST_STATES+(exp-HIST_TIER_LIMIT-1);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Register Bitmap Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	private static long unpackReg(int reg){
		if(reg==0){return 0;}
		final int nlzPart=stateToExp(reg);
		final int hist=stateToHist(reg);
		long bitmap=1L<<nlzPart;
		if(nlzPart>=2){bitmap|=(long)(hist&1)<<(nlzPart-2);}
		if(nlzPart>=1){bitmap|=(long)((hist>>>1)&1)<<(nlzPart-1);}
		return bitmap;
	}

	private static int packReg(long bitmap){
		if(bitmap==0){return 0;}
		final int highest=63-Long.numberOfLeadingZeros(bitmap);
		int hist;
		if(highest>=2){
			hist=(int)((bitmap>>>(highest-2))&3);
		}else if(highest==1){
			hist=((int)(bitmap&1))<<1;
		}else{
			hist=0;
		}
		if(highest>MAX_EXP){
			if(MAX_EXP<=HIST_TIER_LIMIT){
				hist=(MAX_EXP>=2) ? (int)((bitmap>>>(MAX_EXP-2))&3) : 0;
				return toState(MAX_EXP, hist);
			}
			return toState(MAX_EXP, 0);
		}
		if(highest>HIST_TIER_LIMIT){hist=0;}
		return toState(highest, hist);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[modBuckets];
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){
				if(globalNLZ>=0 && globalNLZ<64){
					nlzCounts[globalNLZ+1]++;
					filledCount++;
				}
			}else{
				final int nlzPart=stateToExp(reg);
				final int histPattern=stateToHist(reg);
				if(nlzPart>=HISTORY_MARGIN){
					final int absNlz=nlzPart-HISTORY_MARGIN+(globalNLZ+1);
					if(absNlz<64){
						nlzCounts[absNlz+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((absNlz+1)<<2)|histPattern);
					}
				}else{
					if(globalNLZ>=0 && globalNLZ<64){
						nlzCounts[globalNLZ+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((globalNLZ+1)<<2)|histPattern);
					}
				}
			}
		}
		nlzCounts[0]=modBuckets-filledCount;
		return new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		lastSummarized=s;
		long card=Math.max(0, Math.round(s.ldlc()));
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((ArithmeticVariableLogLog2)log);
	}

	public void add(ArithmeticVariableLogLog2 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
		final int deltaA=newGlobalNLZ-globalNLZ;
		final int deltaB=newGlobalNLZ-log.globalNLZ;
		for(int i=0; i<modBuckets; i++){
			final int rA=getReg(i);
			final int rB=log.getReg(i);
			final long hpA=(rA==0) ? 0 : unpackReg(rA)>>>deltaA;
			final long hpB=(rB==0) ? 0 : unpackReg(rB)>>>deltaB;
			final long merged=hpA|hpB;
			if(merged==0){
				setReg(i, 0);
			}else{
				setReg(i, Math.min(packReg(merged), MAX_REGISTER));
			}
		}
		globalNLZ=newGlobalNLZ;
		filledBuckets=0;
		floorCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=getReg(i);
			if(reg>0){filledBuckets++;}
			if(reg==0 || stateToExp(reg)<=HISTORY_MARGIN){floorCount++;}
		}
		while(floorCount==0 && globalNLZ<wordlen){
			globalNLZ++;
			floorCount=countAndDecrement();
		}
		final int exitThreshold=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
	}

	@Override
	public final void hashAndStore(final long number){
		final long key=hash64shift(number^hashXor);

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		if(Long.compareUnsigned(key, eeMask)>0){return;}
		branch1++;

		final int idx=(int)(Long.remainderUnsigned(key, modBuckets));
		final int nlz=Long.numberOfLeadingZeros(key);
		final int relNlz=nlz-(globalNLZ+1);

		final int bitPos=relNlz+HISTORY_MARGIN;
		if(bitPos<0||bitPos>=64){return;}

		final int oldReg=getReg(idx);

		int newReg;
		if(bitPos>HIST_TIER_LIMIT && bitPos<=MAX_EXP){
			newReg=toState(bitPos, 0);
		}else{
			long hashPrefix=unpackReg(oldReg);
			hashPrefix|=1L<<bitPos;
			newReg=packReg(hashPrefix);
		}

		if(newReg<=oldReg){return;}
		branch2++;

		setReg(idx, newReg);
		lastCardinality=-1;

		if(oldReg==0){filledBuckets++;}

		final int oldNlzPart=stateToExp(oldReg);
		final int newNlzPart=stateToExp(newReg);
		if(oldNlzPart<=HISTORY_MARGIN && newNlzPart>HISTORY_MARGIN){
			floorCount--;
			if(floorCount<=0){
				while(floorCount==0 && globalNLZ<wordlen){
					globalNLZ++;
					floorCount=countAndDecrement();
				}
				final int exitThreshold=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
				eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
			}
		}
	}

	private int countAndDecrement(){
		int newFloorCount=0;
		for(int i=0; i<modBuckets; i++){
			final int reg=getReg(i);
			if(reg==0){newFloorCount++; continue;}
			final int exp=stateToExp(reg);
			int newReg;
			if(exp>HIST_TIER_LIMIT+1){
				newReg=reg-1;
			}else if(exp==HIST_TIER_LIMIT+1){
				newReg=toState(exp-1, 0);
			}else if(exp>=3){
				newReg=reg-4;
			}else if(exp==2){
				newReg=toState(1, stateToHist(reg));
			}else if(exp==1){
				newReg=0;
			}else{
				newReg=0;
			}
			setReg(i, newReg);
			if(stateToExp(newReg)<=HISTORY_MARGIN){newFloorCount++;}
		}
		return newFloorCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/modBuckets;}
	public int getMinZeros(){return globalNLZ+1;}

	@Override public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		lastSummarized=summarize();
		final double hybridEst=lastSummarized.hybridDLL();
		final double[] r=AbstractCardStats.buildLegacyArray(lastSummarized, hybridEst);
		r[AbstractCardStats.HC_IDX+1]=lastSummarized.vwMean();
		return r;
	}

	public int getRegPublic(int i){return getReg(i);}
	public int getRegister(int i){return getReg(i);}

	public double[] udlcEstimate(){
		final CardStats s=(lastSummarized!=null) ? lastSummarized : summarize();
		return new double[]{s.ldlc(), s.dlcSbs(), s.hc(), s.lcMin(),
			0, s.hllRaw(), s.meanHistCF(), s.hybridPlus2()};
	}

	public void printRegisters(){
		for(int i=0; i<modBuckets; i++){System.err.print(getReg(i)+" ");}
		System.err.println();
	}

	public int packedBytes(){return packedLen*8;}

	public CardStats consumeLastSummarized(){
		final CardStats s=lastSummarized;
		lastSummarized=null;
		return s;
	}

	public static void reconfigure(){
		HIST_STATES=4*HIST_TIER_LIMIT-1;
		MAX_EXP=HIST_TIER_LIMIT+(RADIX-HIST_STATES);
		MAX_REGISTER=RADIX-1;
		POW=computePow(RADIX);
	}

	private static long hash64shift(long key){
		key=(~key)+(key<<21);
		key=key^(key>>>24);
		key=(key+(key<<3))+(key<<8);
		key=key^(key>>>14);
		key=(key+(key<<2))+(key<<4);
		key=key^(key>>>28);
		key=key+(key<<31);
		return key;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final long[] packed;
	private final int packedLen;
	private final int modBuckets;
	private int globalNLZ=-1;
	private int floorCount;
	private long eeMask=-1L;
	private int filledBuckets=0;
	private CardStats lastSummarized;

	public long branch1=0, branch2=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	static final int HISTORY_MARGIN=2;

	public static boolean SATURATE_ON_OVERFLOW=true;

	public static int DEMOTION_MODE=1;
	public static int HIST_TIER_LIMIT=13;
	static int RADIX=56;
	static int HIST_STATES=4*HIST_TIER_LIMIT-1;
	public static int MAX_EXP=HIST_TIER_LIMIT+(RADIX-HIST_STATES);
	static int MAX_REGISTER=RADIX-1;

	public static final String SBS_FILE="?cardinalityCorrectionLC2BitHist.tsv.gz";

	public static final String CF_FILE="?cardinalityCorrectionAVDLL64.tsv.gz";

	/** Mean CF residual, x=log2(card/B). Fitted from resource table (B=1408).
	 *  R²=0.9999969, transition MaxErr=0.024%, terminal=0.99991. */
	public static final double[] MCF_MEAN_RESIDUAL={
		-4.134681757262724,  4.999999992726675,-11.794568553347373,  1.166556972002922,
		 0.328887715628532,  1.845414705069709,  1.071814615231428,
		-0.194293984917707, -5.783510608691229,  2.244741892627562,
		 0.185448568786388, -3.798664343259118,  3.684574870998338,
		 0.070445731787030,  1.581614638086199,  1.391853026540059};

	/** MeanH CF residual, x=log2(card/B). Fitted from resource table (B=1408).
	 *  R²=0.9999969, transition MaxErr=0.033%, terminal=0.99963. */
	public static final double[] MCF_MEANH_RESIDUAL={
		-3.691071292933851, -0.363637261057505,  2.567263769729832,  1.066355956682916,
		 4.999999754294653,-11.189952323115953,  0.785742127708774,
		 0.054343414450697, -8.335873985296551,  1.101742014940036,
		-0.143358513726359, -0.794735788647950,  2.755625941876753,
		-0.187952399977993,  0.812771424938221,  2.145566972133455};
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);

	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	@Override public float terminalMeanCF(){return 0.721124f;}
	@Override public float terminalMeanPlusCF(){return 0.855976f;}

	private static long[] POW=computePow(RADIX);

	private static long[] computePow(int radix){
		long[] pow=new long[12];
		pow[0]=1;
		for(int i=1; i<12; i++){pow[i]=pow[i-1]*radix;}
		return pow;
	}

}
