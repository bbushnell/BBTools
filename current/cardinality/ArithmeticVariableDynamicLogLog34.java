package cardinality;

import shared.Tools;

/**
 * Arithmetic-encoded variable-history DLL with 6 buckets per 34-bit word.
 * RADIX=50, 34 bits = 32-bit int + 2-bit extra. ~5.6% smaller than UDLL36.
 * T0=1 state, T1=2 states, T2+=4 states. HIST_STATES = 4*HTL - 1.
 *
 * @author Noire
 * @date May 2026
 */
public final class ArithmeticVariableDynamicLogLog34 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ArithmeticVariableDynamicLogLog34(){this(2048, 31, -1, 0);}

	ArithmeticVariableDynamicLogLog34(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ArithmeticVariableDynamicLogLog34(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*6<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*6-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*6;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=(modBuckets+5)/6;
		packed=new int[packedLen];
		extra=new byte[(packedLen+3)/4];
		floorCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+5)/6);}

	@Override public ArithmeticVariableDynamicLogLog34 copy(){return new ArithmeticVariableDynamicLogLog34(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}

	/*--------------------------------------------------------------*/
	/*----------------      Extra Array Access       ----------------*/
	/*--------------------------------------------------------------*/

	private int fetch2bits(int w){
		return (extra[w>>>2]>>>((w&3)<<1))&0x3;
	}

	private void store2bits(int w, int v){
		final int byteIdx=w>>>2;
		final int shift=(w&3)<<1;
		extra[byteIdx]=(byte)((extra[byteIdx]&~(0x3<<shift))|((v&0x3)<<shift));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	int getReg(int i){
		final int wordIdx=i/6;
		final int pos=i%6;
		final long word=(packed[wordIdx]&0xFFFFFFFFL)|((long)fetch2bits(wordIdx)<<32);
		return (int)((word/POW[pos])%RADIX);
	}

	private void setReg(int i, int val){
		final int wordIdx=i/6;
		final int pos=i%6;
		long word=(packed[wordIdx]&0xFFFFFFFFL)|((long)fetch2bits(wordIdx)<<32);
		final int oldVal=(int)((word/POW[pos])%RADIX);
		word=word+(long)(val-oldVal)*POW[pos];
		packed[wordIdx]=(int)word;
		store2bits(wordIdx, (int)(word>>>32)&0x3);
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
		add((ArithmeticVariableDynamicLogLog34)log);
	}

	public void add(ArithmeticVariableDynamicLogLog34 log){
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
		final long key=Tools.hash64shift(number^hashXor);

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
		return AbstractCardStats.buildLegacyArray(lastSummarized, hybridEst);
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

	public int packedBytes(){return packedLen*4+extra.length;}

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

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final int[] packed;
	private final byte[] extra;
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
	public static int HIST_TIER_LIMIT=12;
	static int RADIX=50;
	static int HIST_STATES=4*HIST_TIER_LIMIT-1;
	public static int MAX_EXP=HIST_TIER_LIMIT+(RADIX-HIST_STATES);
	static int MAX_REGISTER=RADIX-1;

	public static final String SBS_FILE="?cardinalityCorrectionAVDLL34_LC2BitHistSBS.tsv.gz";

	public static final String CF_FILE="?cardinalityCorrectionUDLL6.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);

	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	@Override public float terminalMeanCF(){return 1.386726f;}
	@Override public float terminalMeanPlusCF(){return 0.856020f;}

	private static long[] POW=computePow(RADIX);

	private static long[] computePow(int radix){
		long[] pow=new long[6];
		pow[0]=1;
		for(int i=1; i<6; i++){pow[i]=pow[i-1]*radix;}
		return pow;
	}

}
