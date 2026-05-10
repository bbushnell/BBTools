package cardinality;

import shared.Tools;

/**
 * 36-bit word variant of UltraDynamicLogLog6m.
 * Packs 6 six-bit registers per 36-bit logical word using a secondary
 * nibble array for the extra 4 bits beyond the 32-bit int.
 * ~6.6% more registers at equal memory vs UDLL6m.
 *
 * @author UMP45
 * @date April 2026
 */
public final class UltraDynamicLogLog36 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	UltraDynamicLogLog36(){this(2048, 31, -1, 0);}

	UltraDynamicLogLog36(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	UltraDynamicLogLog36(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*6<=0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*6-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*6;
		modBuckets=rounded>0 ? rounded : buckets;
		packedLen=(modBuckets+5)/6;
		packed=new int[packedLen];
		extra=new byte[(packedLen+1)/2];
		floorCount=modBuckets;
	}

	private static int roundToWords(int b){return Math.max(1, (b+5)/6);}

	@Override public UltraDynamicLogLog36 copy(){return new UltraDynamicLogLog36(modBuckets, k, -1, minProb);}
	@Override public int actualBuckets(){return modBuckets;}
	@Override public int bitsPerWord(){return 36;}
	@Override public int bucketsPerWord(){return 6;}

	/*--------------------------------------------------------------*/
	/*----------------      Extra Array Access       ----------------*/
	/*--------------------------------------------------------------*/

	private int fetch4bits(int w){
		return (extra[w>>>1]>>>((w&1)<<2))&0xF;
	}

	private void store4bits(int w, int v){
		final int byteIdx=w>>>1;
		final int shift=(w&1)<<2;
		extra[byteIdx]=(byte)((extra[byteIdx]&~(0xF<<shift))|((v&0xF)<<shift));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	int getReg(int i){
		final int wordIdx=i/6;
		final int pos=i%6;
		final long word=(packed[wordIdx]&0xFFFFFFFFL)|((long)fetch4bits(wordIdx)<<32);
		return (int)((word>>>(pos*6))&0x3F);
	}

	private void setReg(int i, int val){
		final int wordIdx=i/6;
		final int pos=i%6;
		final int shift=pos*6;
		long word=(packed[wordIdx]&0xFFFFFFFFL)|((long)fetch4bits(wordIdx)<<32);
		word=(word&~(0x3FL<<shift))|((long)(val&0x3F)<<shift);
		packed[wordIdx]=(int)word;
		store4bits(wordIdx, (int)(word>>>32)&0xF);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Register Bitmap Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Reconstruct 64-bit observation bitmap from 6-bit packed register. */
	private static long unpackReg(int reg){
		if(reg==0){return 0;}
		final int nlzPart=reg>>>2;
		final int hist=reg&3;
		long bitmap=1L<<nlzPart;
		if(nlzPart>=2){bitmap|=(long)(hist&1)<<(nlzPart-2);}
		if(nlzPart>=1){bitmap|=(long)((hist>>>1)&1)<<(nlzPart-1);}
		return bitmap;
	}

	/** Pack 64-bit observation bitmap back to 6-bit register, clamping at MAX_REGISTER. */
	private static int packReg(long bitmap){
		if(bitmap==0){return 0;}
		final int highest=63-Long.numberOfLeadingZeros(bitmap);
		final int hist;
		if(highest>=2){
			hist=(int)((bitmap>>>(highest-2))&3);
		}else if(highest==1){
			hist=((int)(bitmap&1))<<1;
		}else{
			hist=0;
		}
		final int raw=(highest<<2)|hist;
		if(raw>MAX_REGISTER){
			return (15<<2)|((int)((bitmap>>>13)&3));
		}
		return raw;
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
				final int nlzPart=reg>>>2;
				final int histPattern=reg&3;
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
		add((UltraDynamicLogLog36)log);
	}

	public void add(UltraDynamicLogLog36 log){
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
			if(reg==0 || (reg>>>2)<=HISTORY_MARGIN){floorCount++;}
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
		long hashPrefix=unpackReg(oldReg);
		hashPrefix|=1L<<bitPos;
		final int newReg=packReg(hashPrefix);
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
			final int newReg=reg-4;
			setReg(i, newReg);
			if((newReg>>>2)<=HISTORY_MARGIN){newFloorCount++;}
		}
		return newFloorCount;
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/modBuckets;}
	public int getMinZeros(){return globalNLZ+1;}

	@Override public final float[] compensationFactorLogBucketsArray(){return null;}

	public static final int MEAN16_RAW_IDX=AbstractCardStats.HC_IDX+1;

	@Override
	public double[] rawEstimates(){
		lastSummarized=summarize();
		final double hybridEst=lastSummarized.hybridDLL();
		double[] r=AbstractCardStats.buildLegacyArray(lastSummarized, hybridEst);
		if(VW_STATE_TABLE!=null){
			if(r.length<=MEAN16_RAW_IDX){
				r=java.util.Arrays.copyOf(r, MEAN16_RAW_IDX+1);
			}
			r[MEAN16_RAW_IDX]=vwMeanEstimate();
		}
		return r;
	}

	public int getRegPublic(int i){return getReg(i);}
	public int getRegister(int i){return getReg(i);}

	public static final int VWMEAN_IDX=8;

	public double[] udlcEstimate(){
		final CardStats s=(lastSummarized!=null) ? lastSummarized : summarize();
		return new double[]{s.ldlc(), s.dlcSbs(), s.hc(), s.lcMin(),
			0, s.hllRaw(), s.meanHistCF(), s.hybridPlus2(),
			vwMeanEstimate()};
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
	static final int MAX_REGISTER=63;

	public static boolean SATURATE_ON_OVERFLOW=true;

	public static final String CF_FILE=UltraDynamicLogLog6.CF_FILE;
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);

	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	@Override public float terminalMeanCF(){return 0.721255f;}
	@Override public float terminalMeanPlusCF(){return 0.856115f;}

	/*--------------------------------------------------------------*/
	/*----------------       VWMean Estimate        ----------------*/
	/*--------------------------------------------------------------*/

	static double[][] VW_STATE_TABLE;  // alias for CardStats.VW_STATE_TABLE
	static int VW_MAX_TIER;  // alias for CardStats.VW_MAX_TIER

	public double vwMeanEstimate(){
		if(lastSummarized!=null){return lastSummarized.vwMean();}
		if(VW_STATE_TABLE==null){return 0;}
		final CardStats cs=summarize();
		return cs.vwMean();
	}

	static void loadVWTable(){
		final String override=System.getProperty("vw.table");
		final String fname=(override!=null) ? override : dna.Data.findPath("?perTierStateUDLL6.tsv.gz");
		if(fname==null){return;}
		fileIO.ByteFile bf=fileIO.ByteFile.makeByteFile(fname, false);
		final int maxTiers=64;
		final int numStates=8;
		final double[][] table=new double[maxTiers][numStates];
		int maxTier=0;
		int stopTier=-1;
		boolean inSection=false;
		final String sectionHeader="#"+CardStats.VW_SECTION;
		for(byte[] raw=bf.nextLine(); raw!=null; raw=bf.nextLine()){
			if(raw.length==0){continue;}
			final String line=new String(raw).trim();
			if(line.startsWith(sectionHeader+"\t") || line.equals(sectionHeader)){
				inSection=true;
				if(line.contains("\t")){
					try{stopTier=Integer.parseInt(line.substring(line.indexOf('\t')+1).trim());}
					catch(NumberFormatException e){}
				}
				continue;
			}
			if(inSection && line.startsWith("#Tier")){continue;}
			if(inSection && line.startsWith("#")){break;}
			if(!inSection){continue;}
			final String[] parts=line.split("\t");
			if(parts.length<2){continue;}
			final int tier;
			try{tier=Integer.parseInt(parts[0]);}catch(NumberFormatException e){continue;}
			if(tier>=maxTiers){continue;}
			for(int s=1; s<parts.length && (s-1)<numStates; s++){
				try{table[tier][s-1]=Double.parseDouble(parts[s]);}catch(NumberFormatException e){}
			}
			if(tier>maxTier){maxTier=tier;}
		}
		bf.close();
		final int safeTier=(stopTier>0 ? stopTier-CardStats.VW_TAIL_BITS : maxTier);
		CardStats.VW_MAX_TIER=Math.min(maxTier, safeTier);
		CardStats.VW_STATE_TABLE=table;
		VW_STATE_TABLE=table;
		VW_MAX_TIER=CardStats.VW_MAX_TIER;

		if(CardStats.VW_MAX_TIER>=2){
			double sumA=0, sumB=0;
			for(int tier=Math.max(0, CardStats.VW_MAX_TIER-4)+1; tier<=CardStats.VW_MAX_TIER; tier++){
				for(int s=0; s<numStates; s++){
					if(table[tier][s]>0 && table[tier-1][s]>0){
						sumB+=table[tier][s];
						sumA+=table[tier-1][s];
					}
				}
			}
			if(sumA>0){CardStats.VW_TIER_GROWTH=sumB/sumA;}
		}

		CardStats.VW_EXTRAP=new double[64];
		CardStats.VW_EXTRAP[0]=1.0;
		for(int n=1; n<64; n++){
			int t=CardStats.VW_MAX_TIER+n;
			double g=CardStats.VW_GROWTH_L-CardStats.VW_GROWTH_A*Math.exp(-CardStats.VW_GROWTH_B*(t-CardStats.VW_GROWTH_C));
			CardStats.VW_EXTRAP[n]=CardStats.VW_EXTRAP[n-1]*g;
		}

		loadVWVarianceTable(fname);
	}

	static void loadVWVarianceTable(String fname){
		String varSection;
		if(CardStats.VW_SECTION.startsWith("all_")){varSection="var_all";}
		else if(CardStats.VW_SECTION.startsWith("entry_")){varSection="var_entry";}
		else if(CardStats.VW_SECTION.startsWith("entryexit_")){varSection="var_entryexit";}
		else if(CardStats.VW_SECTION.startsWith("pctile_")){varSection="var_pctile";}
		else{varSection="var_all";}

		fileIO.ByteFile bf=fileIO.ByteFile.makeByteFile(fname, false);
		final int maxTiers=64;
		final int numCondensed=8;
		final double[][] varTable=new double[maxTiers][numCondensed];
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
			for(int s=1; s<parts.length && (s-1)<numCondensed; s++){
				try{varTable[tier][s-1]=Double.parseDouble(parts[s]);}catch(NumberFormatException e){}
			}
			if(tier>maxTier){maxTier=tier;}
		}
		bf.close();
		if(maxTier>0){CardStats.VW_VARIANCE_TABLE=varTable;}
	}

	static {
		loadVWTable();
	}
}
