package cardinality;

import parse.Parser;
import shared.Tools;
import structures.LongList;

/**
 * BankedDynamicLogLog3: DLL3 variant with per-word banked shared exponents.
 * <p>
 * Like DLL3, packs 10 buckets into each int, 3 bits per bucket (30 bits).
 * The 2 unused bits (30-31) become a local bank exponent per word.
 * <p>
 * Absolute NLZ = globalNLZ + 1 + bankExponent + (stored - 1) = stored + globalNLZ + bank.
 * globalNLZ == -1 means nothing has been seen yet; >= 0 means all buckets
 * have absNlz >= globalNLZ.
 * <p>
 * Bank promotion: when a register would overflow (localRelNlz >= 7), if all
 * 10 registers in the word are >= 1, subtract 1 from all and increment bank.
 * This localizes overflow handling — only 10 registers shift, not the whole
 * structure. Extends effective per-word range from 7 to 10 tiers.
 * <p>
 * Global promotion absorption: when globalNLZ advances, words with
 * bankExponent > 0 just decrement their bank instead of touching registers.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class BankedDynamicLogLog3 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	BankedDynamicLogLog3(){
		this(2048, 31, -1, 0);
	}

	BankedDynamicLogLog3(Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	BankedDynamicLogLog3(int buckets_, int k_, long seed, float minProb_){
		super(roundToWords(buckets_)*10 <= 0 ? buckets_ :
			Integer.highestOneBit(roundToWords(buckets_)*10-1)<<1,
			k_, seed, minProb_);
		final int rounded=roundToWords(buckets_)*10;
		modBuckets=rounded>0 ? rounded : buckets;
		words=modBuckets/10;
		maxArray=new int[words];
		minZeroCount=modBuckets;
		xOverflow=modBuckets*Math.log(2.0*modBuckets)/256.0;
		overflowExpFactor=Math.exp(-xOverflow/modBuckets);
		storedOverflow=new int[64];
		promoteThreshold=(PROMOTE_FRAC>0 ? (int)(modBuckets*PROMOTE_FRAC) : 0);
		if(FAST_COUNT){nlzCounts=new int[66];}
	}

	/** Rounds bucket count to the nearest multiple of 10 (minimum 10). */
	private static int roundToWords(int b){return Math.max(1, (b+5)/10);}

	@Override
	public BankedDynamicLogLog3 copy(){return new BankedDynamicLogLog3(modBuckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Bucket Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads the 3-bit stored value for bucket i (0=empty, 1-7=relNlz+1). */
	private int readBucket(final int i){
		return (maxArray[i/10]>>>((i%10)*3))&0x7;
	}

	/** Writes a 3-bit stored value for bucket i. val must be in [0,7]. */
	private void writeBucket(final int i, final int val){
		final int wordIdx=i/10;
		final int shift=(i%10)*3;
		maxArray[wordIdx]=(maxArray[wordIdx]&~(0x7<<shift))|((val&0x7)<<shift);
	}

	/** Reads the 2-bit bank exponent for the word containing bucket i. */
	private int readBank(final int wordIdx){
		return (maxArray[wordIdx]>>>30)&0x3;
	}

	/** Writes the 2-bit bank exponent for a word. val must be in [0,3]. */
	private void writeBank(final int wordIdx, final int val){
		maxArray[wordIdx]=(maxArray[wordIdx]&0x3FFFFFFF)|((val&0x3)<<30);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		double rawHyb=s.hybridDLL();
		if(IGNORE_OVERFLOW && USE_IO_BIAS){rawHyb=ioBiasCorrect(rawHyb, modBuckets);}
		long card=(long)(rawHyb);
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((BankedDynamicLogLog3)log);
	}

	/**
	 * Merges another BDLL3 into this one.
	 * Bank exponents reset to 0 in the merged result for simplicity.
	 * Each bucket is converted to absolute NLZ then reframed.
	 */
	public void add(BankedDynamicLogLog3 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
			// Clear all bank exponents — merge resets to bank=0
			for(int w=0; w<maxArray.length; w++){
				writeBank(w, 0);
			}
			for(int i=0; i<modBuckets; i++){
				final int wordA=i/10, wordB=i/10;
				final int bankA=readBank(wordA);
				final int bankB=log.readBank(wordB);
				final int sA=readBucket(i);
				final int sB=log.readBucket(i);
				// Convert to absolute, then to new relative frame (bank=0)
				final int absA=sA+globalNLZ+bankA;
				final int absB=sB+log.globalNLZ+bankB;
				final int absMax=Math.max(absA, absB);
				final int newStored=(absMax<newGlobalNLZ+1 ? 0 : Math.min(absMax-newGlobalNLZ, 7));
				writeBucket(i, newStored);
			}
			globalNLZ=newGlobalNLZ;
			filledBuckets=0;
			minZeroCount=0;
			if(FAST_COUNT){java.util.Arrays.fill(nlzCounts, 0);}
			for(int i=0; i<modBuckets; i++){
				final int s=readBucket(i);
				if(s>0){
					filledBuckets++;
					if(FAST_COUNT){final int absNlz=s+globalNLZ; if(absNlz>=0 && absNlz<64){nlzCounts[absNlz+1]++;}}
				}else if(FAST_COUNT){
					final int absNlz=globalNLZ; // stored=0+globalNLZ+0 (bank cleared)
					if(absNlz>=0 && absNlz<64){nlzCounts[absNlz+1]++;}
				}
				if(s==0 || (!EARLY_PROMOTE && s==1)){minZeroCount++;}
			}
			while(minZeroCount<=promoteThreshold && globalNLZ<wordlen){
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		if(Long.compareUnsigned(key, eeMask)>0){return;}
		branch1++;
		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
		final int wordIdx=bucket/10;
		final int bank=readBank(wordIdx);
		int localRelNlz=nlz-globalNLZ-1-bank;

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		// Try bank promotion if overflow would occur
		if(localRelNlz>=7 && bank<3){
			if(canPromoteBank(wordIdx)){
				promoteBank(wordIdx);
				localRelNlz++; // bank increased by 1, so localRelNlz decreases by 1
				// Wait — localRelNlz = nlz - minZeros - bank.
				// After promotion: bank is now bank+1, so localRelNlz = nlz - minZeros - (bank+1) = old - 1
				localRelNlz=nlz-globalNLZ-1-readBank(wordIdx); // recompute to be safe
			}
		}

		if(localRelNlz+1>7){totalOverflows++;}
		if(IGNORE_OVERFLOW && localRelNlz+1>7){return;} // silently ignore overflow
		final int newStored=Math.min(localRelNlz+1, 7);
		if(newStored<1){return;} // below current floor
		// Read oldStored and currentBank AFTER potential bank promotion, so both
		// reflect the same post-promotion state of the word.
		final int curBank=readBank(wordIdx);
		final int oldStored=readBucket(bucket);

		if(newStored<=oldStored){return;}
		branch2++;
		lastCardinality=-1;

		writeBucket(bucket, newStored);
		if(oldStored==0){filledBuckets++;}

		if(FAST_COUNT){
			// Remove old slot from nlzCounts
			final int oldAbsNlz=oldStored+globalNLZ+curBank;
			if(oldAbsNlz>=0 && oldAbsNlz<64){
				nlzCounts[oldAbsNlz+1]--;
			}else{
				nlzCounts[0]--; // was truly empty (globalNLZ==-1 and stored==0)
			}
			// Add new slot
			final int newAbsNlz=newStored+globalNLZ+curBank;
			if(newAbsNlz<64){nlzCounts[newAbsNlz+1]++;}
		}

		// minZeroCount tracking.
		// A stored=0 slot is in minZeroCount only when curBank==0 (truly at the floor).
		// Banked-phantom stored=0 slots (curBank>0) are NOT tracked in minZeroCount —
		// they had data, bank promotion decremented them, and they are counted only
		// when countAndDecrement drops their bank to 0.  Decrementing minZeroCount
		// for a banked-phantom fill would drain the counter prematurely.
		// NOTE: curBank and oldStored are both read post-promotion, so they are
		// consistent — if bank promotion just fired, curBank reflects the new bank.
		final int oldRelNlz=(oldStored==0 ? 0 : oldStored-1);
		final boolean shouldDecrement=EARLY_PROMOTE ? (oldStored==0 && curBank==0) : (localRelNlz>oldRelNlz && oldRelNlz==0 && curBank==0);
		if(shouldDecrement && --minZeroCount<=promoteThreshold){
			while(minZeroCount<=promoteThreshold && globalNLZ<wordlen){
				if(USE_STORED_OVERFLOW){
					// Group stored=7 counts by bank: a stored=7 bucket in bank=b
					// overflows at tier 8+globalNLZ+b (= 7+minZeros+b).
					int[] topByBank=new int[4];
					for(int i=0; i<modBuckets; i++){
						if(readBucket(i)==7){topByBank[readBank(i/10)]++;}
					}
					for(int b=0; b<4; b++){
						int tier=8+globalNLZ+b;
						if(tier<64){storedOverflow[tier]+=topByBank[b]/2;}
					}
				}
				globalNLZ++;
				eeMask>>>=1;
				minZeroCount=countAndDecrement();
			}
		}
	}

	/** Returns true if all 10 registers in the word are non-empty (stored >= 1). */
	private boolean canPromoteBank(final int wordIdx){
		final int word=maxArray[wordIdx]&0x3FFFFFFF; // mask off bank bits
		for(int b=0; b<10; b++){
			if(((word>>>(b*3))&0x7)==0){return false;}
		}
		return true;
	}

	/** Subtracts 1 from all 10 registers in a word and increments its bank exponent. */
	private void promoteBank(final int wordIdx){
		int word=maxArray[wordIdx];
		final int oldBank=(word>>>30)&0x3;
		int result=0;
		for(int b=0; b<10; b++){
			final int shift=b*3;
			int stored=(word>>>shift)&0x7;
			stored--; // guaranteed >= 1 by canPromoteBank check
			result|=(stored<<shift);
		}
		// Write registers + new bank
		maxArray[wordIdx]=result|((oldBank+1)<<30);
	}

	/**
	 * Global tier promotion. For words with bank > 0, absorb by decrementing bank.
	 * For words with bank == 0, decrement registers as in DLL3.
	 */
	private int countAndDecrement(){
		int newMinZeroCount=0;
		for(int w=0; w<maxArray.length; w++){
			final int bank=readBank(w);
			if(bank>0){
				// Absorb global promotion: decrement bank, leave registers untouched
				final int newBank=bank-1;
				writeBank(w, newBank);
				// Only count stored=0 slots if the bank has dropped to 0.
				// While bank>0, stored=0 slots are banked phantoms (had data, bank
				// promotion decremented them) — they are NOT at the minimum floor
				// and must NOT be added to minZeroCount.  Only when bank reaches 0
				// do they become true floor-level entries and count toward promotion.
				if(newBank==0){
					for(int b=0; b<10; b++){
						final int stored=(maxArray[w]>>>(b*3))&0x7;
						if(EARLY_PROMOTE){
							if(stored==0){newMinZeroCount++;}
						}else{
							if(stored==0 || stored==1){newMinZeroCount++;}
						}
					}
				}
			}else{
				// bank==0: standard DLL3 decrement
				int word=maxArray[w];
				if((word&0x3FFFFFFF)==0){
					// All registers empty — count them toward minZeroCount
					final int regsInWord=Math.min(10, modBuckets-w*10);
					newMinZeroCount+=regsInWord;
					continue;
				}
				int result=0;
				for(int b=0; b<10; b++){
					final int shift=b*3;
					int stored=(word>>>shift)&0x7;
					if(stored>0){
						stored--;
						if(EARLY_PROMOTE){
							if(stored==0){newMinZeroCount++; filledBuckets--;}
						}else{
							if(stored==1){newMinZeroCount++;}
						}
					}
					result|=(stored<<shift);
				}
				maxArray[w]=result; // bank stays 0
			}
		}
		return newMinZeroCount;
	}

	public int filledBuckets(){return filledBuckets;}
	/** True occupancy: a bucket is "empty" only when stored==0 AND minZeros+bank==0.
	 *  Once minZeros>0 (or any word has bank>0), phantom buckets are informative and
	 *  count as occupied. The raw filledBuckets field lags because EARLY_PROMOTE
	 *  global promotion decrements it when stored drops 1→0, even though that bucket
	 *  becomes a valid phantom at the new floor. */
	public double occupancy(){
		if(globalNLZ>=0){return 1.0;}
		int empty=0;
		for(int w=0; w<maxArray.length; w++){
			if(readBank(w)>0){continue;}
			final int word=maxArray[w];
			for(int b=0; b<10; b++){
				if(((word>>>(b*3))&0x7)==0){empty++;}
			}
		}
		return 1.0-(double)empty/modBuckets;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Builds absolute NLZ histogram including bank exponents, delegates to CardStats.
	 */
	private CardStats summarize(){
		if(nlzCounts==null){nlzCounts=new int[66];}
		else{java.util.Arrays.fill(nlzCounts, 0);}
		int filledCount=0;
		for(int i=0; i<modBuckets; i++){
			final int wordIdx=i/10;
			final int bank=readBank(wordIdx);
			final int stored=readBucket(i);
			final int absNlz=stored+globalNLZ+bank;
			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				filledCount++;
			}
			// else: absNlz==-1 means stored=0+globalNLZ(=-1)+bank(=0) → truly empty
		}
		nlzCounts[0]=modBuckets-filledCount;
		assert nlzCounts[0]>=0 : "Negative empties: "+nlzCounts[0]+" filledCount="+filledCount+" modBuckets="+modBuckets;
		final boolean doDebug=DEBUG_ONCE;
		if(doDebug){
			DEBUG_ONCE=false;
			System.err.println("DEBUG summarize(): globalNLZ="+globalNLZ+" filledBuckets="+filledBuckets+" modBuckets="+modBuckets);
			System.err.println("  nlzCounts (nonzero):");
			for(int t=0; t<66; t++){if(nlzCounts[t]>0){System.err.println("    tier "+(t-1)+": "+nlzCounts[t]);}}
		}
		final int[] counts;
		if(CORRECT_OVERFLOW && globalNLZ>=0){
			final int lo=Math.max(7, globalNLZ+1);
			final int hi=Math.min(10+globalNLZ, 63); // 7+maxBank(3): bank=b overflows at 7+(globalNLZ+1)+b
			final int[] cumRaw=new int[64];
			cumRaw[63]=nlzCounts[64]; // absNlz=63 at index 64
			for(int t=62; t>=0; t--){cumRaw[t]=cumRaw[t+1]+nlzCounts[t+1];}
			final int[] corrCum=cumRaw.clone();
			for(int t=lo; t<=hi; t++){
				final double x=(USE_STORED_OVERFLOW && storedOverflow[t]>0)
					? storedOverflow[t]*OVERFLOW_SCALE : xOverflow*OVERFLOW_SCALE;
				final int addToCum=(int)Math.round(
					(modBuckets-cumRaw[t])*(1.0-Math.exp(-x/modBuckets)));
				corrCum[t]+=addToCum;
			}
			final int maxHi=Math.min(hi+1, 63);
			if(corrCum[hi]>0 && maxHi>hi){
				corrCum[maxHi]=corrCum[hi]/2;
			}
			counts=nlzCounts.clone();
			if(lo>0){counts[lo]=corrCum[lo-1]-corrCum[lo];}
			for(int t=lo; t<maxHi; t++){counts[t+1]=corrCum[t]-corrCum[t+1];}
			counts[maxHi+1]=corrCum[maxHi];
			int corrSum=0; for(int t=1; t<66; t++){corrSum+=counts[t];}
			counts[0]=modBuckets-corrSum;
		}else{
			counts=nlzCounts;
		}
		lastRawNlz=nlzCounts.clone();
		lastCorrNlz=(counts==nlzCounts) ? lastRawNlz : counts;
		return new CardStats(null, counts, 0, 0, 0, 0, modBuckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
	}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		double[] r=AbstractCardStats.buildLegacyArray(s, hybridEst);
		if(IGNORE_OVERFLOW && USE_IO_BIAS){
			final double mult=ioBiasMult(hybridEst, modBuckets);
			for(int i=0; i<r.length; i++){r[i]*=mult;}
		}
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 3-bit bucket array with 2-bit bank exponent in bits 30-31. */
	private final int[] maxArray;
	private final double xOverflow;
	/** Actual bucket count (multiple of 10, may be non-power-of-2). */
	private final int modBuckets;
	/** Number of 32-bit words (modBuckets / 10). */
	private final int words;
	private final double overflowExpFactor;
	private final int[] storedOverflow;
	private final int promoteThreshold;
	/** globalNLZ==-1 means nothing seen yet; >=0 means all buckets have absNlz >= globalNLZ. */
	private int globalNLZ=-1;
	private int minZeroCount;
	private int filledBuckets=0;
	private long eeMask=-1L;

	/** Compatibility accessor: returns globalNLZ+1 to match legacy minZeros convention. */
	public int getMinZeros(){return globalNLZ+1;}
	@Override public int actualBuckets(){return modBuckets;}
	public int[] getStoredOverflow(){return storedOverflow;}
	int[] lastRawNlz, lastCorrNlz;

	public long branch1=0, branch2=0;
	public long totalOverflows=0;
	private int peakOverflows=0;
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/** Counts buckets currently at stored=7, updates peakOverflows. */
	public int countCurrentOverflows(){
		int count=0;
		for(int i=0; i<modBuckets; i++){if(readBucket(i)==7){count++;}}
		peakOverflows=Math.max(peakOverflows, count);
		return count;
	}

	public void printOverflowStats(){
		final int current=countCurrentOverflows();
		System.err.println("BDLL3 overflow: total="+totalOverflows
			+" updates="+branch2
			+" overflowRate="+(branch2>0 ? String.format("%.4f", (double)totalOverflows/branch2) : "N/A")
			+" current="+current
			+" peak="+peakOverflows
			+" peakFrac="+String.format("%.4f", (double)peakOverflows/modBuckets));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	private static final int wordlen=64;
	public static final boolean FAST_COUNT=false;
	public static boolean EARLY_PROMOTE=true;
	public static float PROMOTE_FRAC=0.004f;
	public static boolean CORRECT_OVERFLOW=true;
	public static double OVERFLOW_SCALE=1.0;
	public static boolean USE_STORED_OVERFLOW=true;
	/** When true, hashes that would overflow the register are silently ignored.
	 *  Eliminates overflow corruption at the cost of delayed information.
	 *  The CF table absorbs the resulting systematic bias. */
	public static boolean IGNORE_OVERFLOW=false;
	/** When true, apply per-tier IO_BIAS correction in IO mode. */
	public static boolean USE_IO_BIAS=true;
	public static boolean DEBUG_ONCE=false;

	/** Per-tier bias constants for IO mode. Indexed by (tier-IO_TIER_MIN)/IO_TIER_STEP.
	 *  Can be overridden at runtime via iobfile= flag. */
	private static double[] IO_BIAS={
		-0.002696, -0.004979, -0.006794, -0.007470, -0.007553, -0.007582, -0.007567, -0.007543, // tiers -8..-4.5
		-0.007556, -0.007578, -0.007660, -0.007811, -0.008060, -0.008445, -0.009044, -0.009981, // tiers -4..-0.5
		-0.011518, -0.014230, -0.018743, -0.024282, -0.030522, -0.036442, -0.042419, -0.047468, // tiers 0..3.5
		-0.053413, -0.058115, -0.063512, -0.067936, -0.071727, -0.075495, -0.076534, -0.078819, // tiers 4..7.5
		-0.077188, -0.078890, -0.076174, -0.078573, -0.075200, -0.078367, -0.074412, -0.078283, // tiers 8..11.5
		-0.073747                                                                                 // tier 12
	};
	static double IO_TIER_STEP=0.5;
	static final int IO_TIER_MIN=-8, IO_TIER_MAX=12;

	public static void loadIOBias(String path){
		try{
			java.util.ArrayList<Double> vals=new java.util.ArrayList<>();
			java.io.BufferedReader br=new java.io.BufferedReader(new java.io.FileReader(path));
			for(String line=br.readLine(); line!=null; line=br.readLine()){
				line=line.trim();
				if(line.isEmpty() || line.startsWith("#")){continue;}
				String[] parts=line.split("\\s+");
				vals.add(Double.parseDouble(parts[0]));
			}
			br.close();
			final int fullCount=IO_TIER_MAX-IO_TIER_MIN+1; // 21
			final int halfCount=(IO_TIER_MAX-IO_TIER_MIN)*2+1; // 41
			if(vals.size()==fullCount){
				IO_TIER_STEP=1.0;
			}else if(vals.size()==halfCount){
				IO_TIER_STEP=0.5;
			}else{
				throw new RuntimeException("IO_BIAS file "+path+" has "+vals.size()+" values, expected "+fullCount+" or "+halfCount);
			}
			IO_BIAS=new double[vals.size()];
			for(int i=0; i<vals.size(); i++){IO_BIAS[i]=vals.get(i);}
		}catch(java.io.IOException e){throw new RuntimeException("Failed to load IO_BIAS from "+path, e);}
	}

	public static void setIOBias(double[] bias){IO_BIAS=bias;}

	/** Returns IO-bias-corrected estimate: raw / (1 + bias[tier]). */
	static double ioBiasCorrect(double raw, int buckets){
		return raw*ioBiasMult(raw, buckets);
	}

	/** Returns the multiplicative correction factor for IO mode. */
	static double ioBiasMult(double raw, int buckets){
		if(raw<=0){return 1;}
		double ftier=Math.log(raw/buckets)/Math.log(2);
		double fidx=(ftier-IO_TIER_MIN)/IO_TIER_STEP;
		int lo=(int)Math.floor(fidx);
		lo=Math.max(0, Math.min(IO_BIAS.length-1, lo));
		int hi=Math.min(lo+1, IO_BIAS.length-1);
		double frac=fidx-Math.floor(fidx);
		double bias=IO_BIAS[lo]+(IO_BIAS[hi]-IO_BIAS[lo])*frac;
		return 1.0/(1.0+bias);
	}

	/** Use DLL3's CF table initially — will need its own table eventually. */
	public static final String CF_FILE="?cardinalityCorrectionDLL3.tsv.gz";
	private static int CF_BUCKETS=8192;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	/** Asymptotic meanRaw/trueCard ratio at saturation (co=f, buckets=2048).
	 *  Measured 16k ddls, maxmult=4000, tmcf=1, cf=f (Apr 15 2026): 1.4394.
	 *  Prior value 0.691154 produced +108% Mean_err at saturation — confirmed wrong. */
	@Override public float terminalMeanCF(){return 0.694734f;}

}
