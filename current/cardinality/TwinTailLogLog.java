package cardinality;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;

/**
 * TTLL: TwinTailLogLog — ceiling-based dual-history cardinality estimator.
 *
 * Each bucket is one 8-bit byte:
 *   [7:4] = 4-bit localExp (offset from globalExp)
 *   [3:2] = history1 (MSB=hit at ceiling, LSB=hit at ceiling-1; for histBit==1)
 *   [1:0] = history0 (same semantics; for histBit==0)
 *
 * histBit = the hash bit immediately above the bucket-index bits.
 * Power-of-2 bucket count only.
 *
 * Update rule:
 *   delta = hashNLZ - (globalExp + localExp)
 *   delta &lt; -1 : ignore
 *   delta == -1: set LSB of history[histBit]
 *   delta == 0 : set MSB of history[histBit]
 *   delta &gt; 0 : advance localExp; shift both histories right by min(delta,2);
 *                set MSB of history[histBit]; fire advanceGlobal if needed
 *   Overflow (newLocalExp &gt; 15): ignore hash entirely.
 *
 * Two-level exponent: globalExp = shared minimum floor across all buckets.
 * advanceGlobal fires when all buckets have localExp &gt;= 1.
 *
 * Estimation: state-table CF (16 combined_h states per tier), CV-weighted,
 * plus LC/DLC/HC via variable-observation-count framework.
 *
 * @author Brian, Ady
 * @date April 10, 2026
 */
public final class TwinTailLogLog extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int NUM_TIERS   =16;
	static final int NUM_COMBINED=16;// combined_h: (h1<<2)|h0
	/** MSB mask for combined_h: bits 1 and 3. */
	static final int CH_MSB_MASK =0xA;// 0b1010

	/** Encoding mode.
	 *  0 = symmetric: both tails use histBit to select which tail updates.
	 *      h0 responds to histBit==0, h1 responds to histBit==1.
	 *      MSB = "saw current NLZ with my histBit", LSB = "saw NLZ-1 with my histBit".
	 *  1 = master/slave: h0 is pure NLZ history (unconditional, UDLL6-style).
	 *      h0 MSB = "saw NLZ-1", h0 LSB = "saw NLZ-2".
	 *      h1 MSB = "saw current NLZ with histBit=1", h1 LSB = "saw NLZ-1 with histBit=1". */
	public static int ENCODING_MODE=0;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	TwinTailLogLog(){this(2048, 31, -1, 0);}

	TwinTailLogLog(parse.Parser p){
		super(p);
		regs=new byte[buckets];
		numFloorBuckets=buckets;
	}

	TwinTailLogLog(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		regs=new byte[buckets];
		numFloorBuckets=buckets;
	}

	@Override
	public TwinTailLogLog copy(){return new TwinTailLogLog(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------         Core Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);

		// Early exit
		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int hashNLZ=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);
		final int histBit=(int)((key>>>bucketBits)&1);

		// MicroIndex for very-low-cardinality detection
		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int b=regs[bucket]&0xFF;
		final int localExp=(b>>>4)&0xF;
		final int absNLZ=globalExp+localExp;
		final int delta=hashNLZ-absNLZ;

		if(ENCODING_MODE==0){
			// Symmetric: both tails selected by histBit
			if(delta<-1){return;}
			if(delta>0){
				final int newLocalExp=localExp+delta;
				if(newLocalExp>15){return;}
				final boolean wasZeroExp=(localExp==0);
				final int shiftAmt=Math.min(delta, 2);
				final int oldH=b&0xF;
				final int h0=(oldH&0x3)>>>shiftAmt;
				final int h1=((oldH>>>2)&0x3)>>>shiftAmt;
				int newByte=(newLocalExp<<4)|(h1<<2)|h0;
				final int bitPos=(histBit==0) ? 1 : 3;
				newByte|=(1<<bitPos);
				lastCardinality=-1;
				lastSummarized=null;
				regs[bucket]=(byte)newByte;
				if(wasZeroExp){
					numFloorBuckets--;
					if(numFloorBuckets==0){advanceGlobal();}
				}
			}else{
				final int bitPos;
				if(delta==0){
					bitPos=(histBit==0) ? 1 : 3;
				}else{
					bitPos=(histBit==0) ? 0 : 2;
				}
				final int bitToSet=1<<bitPos;
				if((b&bitToSet)!=0){return;}
				lastCardinality=-1;
				lastSummarized=null;
				regs[bucket]=(byte)(b|bitToSet);
			}
		}else{
			// Master/slave: h0=unconditional NLZ history, h1=bit-filtered
			if(delta<-2){return;}
			if(delta>0){
				final int newLocalExp=localExp+delta;
				if(newLocalExp>15){return;}
				final boolean wasZeroExp=(localExp==0);
				final int oldH0=b&0x3;
				final int oldH1=(b>>>2)&0x3;
				final int carry=(b!=0) ? (1<<2) : 0;
				final int h0=((oldH0|carry)>>>delta)&0x3;
				int h1=(oldH1>>>Math.min(delta, 2))&0x3;
				if(histBit==1){ h1|=0x2; }
				lastCardinality=-1;
				lastSummarized=null;
				regs[bucket]=(byte)((newLocalExp<<4)|(h1<<2)|h0);
				if(wasZeroExp){
					numFloorBuckets--;
					if(numFloorBuckets==0){advanceGlobal();}
				}
			}else if(delta==0){
				if(histBit==0){return;}
				final int bitToSet=1<<3;
				if((b&bitToSet)!=0){return;}
				lastCardinality=-1;
				lastSummarized=null;
				regs[bucket]=(byte)(b|bitToSet);
			}else if(delta==-1){
				int newBits=(1<<1);
				if(histBit==1){ newBits|=(1<<2); }
				if((b&newBits)==newBits){return;}
				lastCardinality=-1;
				lastSummarized=null;
				regs[bucket]=(byte)(b|newBits);
			}else if(delta==-2){
				final int bitToSet=1<<0;
				if((b&bitToSet)!=0){return;}
				lastCardinality=-1;
				lastSummarized=null;
				regs[bucket]=(byte)(b|bitToSet);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Global Advance         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * When all buckets have localExp &gt;= 1, increment globalExp,
	 * decrement all localExps, recount numFloorBuckets, update eeMask.
	 */
	private void advanceGlobal(){
		while(numFloorBuckets==0){
			globalExp++;
			numFloorBuckets=0;
			for(int i=0; i<buckets; i++){
				int b=regs[i]&0xFF;
				int le=(b>>>4)&0xF;
				assert le>=1 : "localExp=0 but numFloorBuckets was 0, bucket "+i;
				le--;
				regs[i]=(byte)((le<<4)|(b&0xF));
				if(le==0){numFloorBuckets++;}
			}
			final int margin=(ENCODING_MODE==0) ? 1 : 2;
			eeMask=(globalExp<=margin) ? -1L : (-1L)>>>(globalExp-margin);
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

	/** Index of Mean16 in rawEstimates() output. */
	public static final int MEAN16_RAW_IDX=AbstractCardStats.HC_IDX+1;

	/**
	 * LDLC estimate with 1-bit history support.
	 * Reuses cached CardStats from rawEstimates() to avoid double-summarize.
	 * Returns {ldlc, dlc, hc, lcMin, fgra, hll, meanH, hybridPlus2, mean16, dualLC}.
	 */
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

	/**
	 * Builds a CardStats from the current register state.
	 * Extracts 1-bit effective history per bucket.
	 * Mode 0 (symmetric): OR of both tails' LSBs — "saw NLZ-1 with either histBit".
	 * Mode 1 (master/slave): h0 MSB — "saw NLZ-1" (unconditional, UDLL6-style).
	 *
	 * Packed bucket format: ((absNlz+1) &lt;&lt; 1) | histBit
	 */
	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[buckets];
		int filledCount=0;

		for(int i=0; i<buckets; i++){
			final int b=regs[i]&0xFF;
			if(b==0 && globalExp==0){
				continue;// truly empty bucket
			}
			final int localExp=(b>>>4)&0xF;
			final int absNlz=globalExp+localExp;

			if(absNlz>=0 && absNlz<64){
				nlzCounts[absNlz+1]++;
				filledCount++;
				final int histBit;
				if(ENCODING_MODE==0){
					// Symmetric: OR of both LSBs
					final int ch=b&0xF;
					histBit=((ch&1)|((ch>>>2)&1));
				}else{
					// Master/slave: h0 MSB (bit 1)
					histBit=(b>>>1)&1;
				}
				packedBuckets[i]=(char)(((absNlz+1)<<1)|histBit);
			}
		}
		nlzCounts[0]=buckets-filledCount;
		lastRawNlz=nlzCounts;

		return new CardStats(packedBuckets, nlzCounts, 0, 1, 0, 0,
				buckets, microIndex, added, null, 0, 0.0,
				terminalMeanCF(), terminalMeanPlusCF());
	}

	/*--------------------------------------------------------------*/
	/*----------------     Dual-Bucket LC/DLC       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * LC estimate treating each tail as an independent virtual bucket.
	 * Symmetric encoding (mode 0) only.
	 *
	 * Each physical bucket contributes 0, 1, or 2 virtual buckets:
	 *   For each tail (h0, h1):
	 *     tail == 10 or 11: filled virtual bucket at tier absNlz
	 *     tail == 01:       filled virtual bucket at tier absNlz-1
	 *     tail == 00:
	 *       absNlz < 2: empty virtual bucket (all possible sub-tiers are visible)
	 *       absNlz >= 2: unknown (observation may have shifted out), exclude
	 *
	 * Total virtual buckets varies per register but is up to 2*physical.
	 * LC = virtualBuckets * ln(virtualBuckets / emptyVirtual)
	 *
	 * @return dual-bucket LC estimate, or 0 if not applicable
	 */
	public double dualBucketLC(){
		if(ENCODING_MODE!=0){return 0;}
		int virtualTotal=0;
		int virtualEmpty=0;
		int[] virtualNlzCounts=new int[66]; // for potential DLC use

		for(int i=0; i<buckets; i++){
			final int b=regs[i]&0xFF;
			if(b==0 && globalExp==0){
				// Truly empty register: 2 empty virtual buckets
				virtualTotal+=2;
				virtualEmpty+=2;
				virtualNlzCounts[0]+=2;
				continue;
			}
			final int localExp=(b>>>4)&0xF;
			final int absNlz=globalExp+localExp;

			// Process each tail independently
			for(int t=0; t<2; t++){
				final int tail=(t==0) ? (b&0x3) : ((b>>>2)&0x3);
				if(tail>=2){
					// MSB set (10 or 11): filled at absNlz
					virtualTotal++;
					if(absNlz>=0 && absNlz<64){virtualNlzCounts[absNlz+1]++;}
				}else if(tail==1){
					// Only LSB set (01): filled at absNlz-1
					virtualTotal++;
					final int prevNlz=absNlz-1;
					if(prevNlz>=0 && prevNlz<64){virtualNlzCounts[prevNlz+1]++;}
				}else{
					// tail == 00: empty if absNlz < 2, else unknown (exclude)
					if(absNlz<2){
						virtualTotal++;
						virtualEmpty++;
						virtualNlzCounts[0]++;
					}
					// absNlz >= 2: don't count this tail at all
				}
			}
		}

		if(virtualTotal==0 || virtualEmpty==0){return 0;}
		return (double)virtualTotal*Math.log((double)virtualTotal/virtualEmpty);
	}

	/*--------------------------------------------------------------*/
	/*----------------    16-State Mean Estimate    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * 16-state Mean estimate using MantissaCompare2-derived log-space CFs.
	 * Same structure as CardStats Phase 2a (Mean+H), but with 16 states
	 * from the 4-bit combined_h instead of 2-4 from histBits.
	 *
	 * Per-bucket correction: tierMult = 2^(-(cf + CF_OFFSET))
	 * where cf = StateTable.ttllOffset(absNlz, combinedH).
	 * Systematic bias is corrected by the v5 CF table (Mean16_cf column).
	 */
	public double ttllMeanEstimate(){
		// Precompute tierMult for all (nlzBin, combinedH) to avoid per-bucket Math.pow
		final int maxNlzBin=4; // tiers 0,1,2 have per-tier CFs; 3+ use steady-state
		final double[][] tierMultTable=new double[maxNlzBin][16];
		for(int nb=0; nb<maxNlzBin; nb++){
			for(int ch=0; ch<16; ch++){
				final double cf=StateTable.ttllOffset(nb, ch);
				tierMultTable[nb][ch]=Math.pow(2.0, -(cf+StateTable.CF_OFFSET));
			}
		}

		double corrDifSum=0;
		int filled=0;
		for(int i=0; i<buckets; i++){
			final int b=regs[i]&0xFF;
			if(b==0 && globalExp==0){continue;}
			final int localExp=(b>>>4)&0xF;
			final int absNlz=globalExp+localExp;
			final int ch=b&0xF;
			final int nlzBin=Math.min(absNlz, maxNlzBin-1);

			final long dif=(absNlz==0 ? Long.MAX_VALUE : (absNlz<64 ? 1L<<(63-absNlz) : 1L));
			corrDifSum+=dif*tierMultTable[nlzBin][ch];
			filled++;
		}
		if(filled==0){return 0;}
		final double mean=corrDifSum/filled;
		final float correction=(filled+buckets)/(float)(buckets+buckets);
		double est=2*(Long.MAX_VALUE/Tools.max(1.0, mean))*filled*correction;

		// Apply v5 CF table correction (Mean16_cf column)
		if(CorrectionFactor.USE_CORRECTION && CorrectionFactor.v1Matrix!=null
				&& CorrectionFactor.MEAN16<CorrectionFactor.v1Matrix.length
				&& CorrectionFactor.v1Matrix[CorrectionFactor.MEAN16]!=null){
			final double keyScale=(CorrectionFactor.v1Buckets>0 && buckets!=CorrectionFactor.v1Buckets)
				? (double)CorrectionFactor.v1Buckets/buckets : 1.0;
			final double cf=CorrectionFactor.getCF(est, est,
				CorrectionFactor.v1Matrix[CorrectionFactor.MEAN16],
				CorrectionFactor.v1Keys, 5, 0.0001, keyScale);
			est*=cf;
		}
		return est;
	}

	/*--------------------------------------------------------------*/
	/*----------------     CF Table Estimate        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Raw CF-table estimate with CV-weighted averaging.
	 * Uses the 16-state TTLL-specific correction table.
	 */
	public static double CV_POWER=2.0;

	private double rawCFEstimate(){
		if(CF_TABLE_TIERS<=0){return 0;}
		double weightedSum=0;
		double weightSum=0;
		for(int i=0; i<buckets; i++){
			final int b=regs[i]&0xFF;
			if(b==0){continue;}
			final int localExp=(b>>>4)&0xF;
			final int ch=b&0xF;
			final int absNLZ=localExp+globalExp;
			final int tier=Math.min(absNLZ, CF_TABLE_TIERS-1);
			final double value=tierCardinality(absNLZ)*cfTable[tier][ch];
			final double cv=(cvTable!=null && tier<cvTable.length)
				? Math.max(cvTable[tier][ch], 0.001) : 1.0;
			final double weight=Math.pow(1.0/cv, CV_POWER);
			weightedSum+=value*weight;
			weightSum+=weight;
		}
		if(weightSum<=0){return 0;}
		double est=(weightedSum/weightSum)*buckets;
		est*=TERMINAL_CORRECTION;
		if(cardCFKeys!=null && cardCFKeys.length>0){
			double keyScale=(cardCFBuckets>0 && buckets!=cardCFBuckets)
				? (double)cardCFBuckets/buckets : 1.0;
			double cf=CorrectionFactor.getCF(est, est,
				cardCFValues, cardCFKeys, 5, 0.0001, keyScale);
			est*=cf;
		}
		return est;
	}

	/** Expected cardinality at a given tier (from simulator tierAvg). */
	private static double tierCardinality(int absNLZ){
		if(tierAvg!=null && absNLZ<tierAvg.length && absNLZ>=0){
			return tierAvg[absNLZ];
		}
		if(tierAvg==null || tierAvg.length<2){return 1.0;}
		double lastGrowth=tierAvg[tierAvg.length-1]/tierAvg[tierAvg.length-2];
		double base=tierAvg[tierAvg.length-1];
		int extra=absNLZ-(tierAvg.length-1);
		return base*Math.pow(lastGrowth, extra);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Merge/Misc          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void add(CardinalityTracker log){
		throw new UnsupportedOperationException("TTLL merge not yet implemented");
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public int getGlobalExp(){return globalExp;}

	public double occupancy(){
		int filled=0;
		for(int i=0; i<buckets; i++){
			if((regs[i]&0xFF)!=0){filled++;}
		}
		return (double)filled/buckets;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final byte[] regs;
	/** Shared minimum floor across all buckets. */
	private int globalExp=0;
	/** Count of buckets with localExp == 0. */
	private int numFloorBuckets;
	/** Early-exit mask: reject hashes with NLZ &lt; globalExp - 1. */
	private long eeMask=-1L;
	/** Cached CardStats from rawEstimates(), consumed by ldlcEstimate(). */
	private CardStats lastSummarized;
	int[] lastRawNlz;

	/*--------------------------------------------------------------*/
	/*----------------        Correction Factors    ----------------*/
	/*--------------------------------------------------------------*/

	public static double TERMINAL_CORRECTION=1.0;

	static double[][] cfTable;  // [tier][combined_h] multipliers
	static double[][] cvTable;  // [tier][combined_h] CV values
	/** Normalized inverse of cfTable for Mean formula correction.
	 *  meanCorrTable[t][ch] = (1/cfTable[t][ch]) / mean(1/cfTable[t][*]).
	 *  Ensures average correction per tier = 1.0, avoiding Jensen's bias. */
	static double[][] meanCorrTable;
	static int CF_TABLE_TIERS=0;
	static double[] tierAvg;

	static float[] cardCFValues;
	static float[] cardCFKeys;
	static int cardCFBuckets=0;

	static final String STATE_TABLE_FILE="stateTableTTLL.tsv.gz";
	static final String CF_FILE="cardinalityCorrectionTTLL.tsv.gz";

	/** MUST appear after all static field declarations. */
	static {
		loadCFTable();
		loadCardCFTable();
	}

	/**
	 * Loads the state table produced by TTLLSimulator.
	 * Reads sections: #TierAvg (computes blend avg), #StateBlendMult, #StateCV.
	 * The table format uses # for section/column headers, no # for data rows.
	 */
	public static synchronized void loadCFTable(){
		if(cfTable!=null){return;}
		String path=Data.findPath("?"+STATE_TABLE_FILE);
		if(path==null){setFallbackCF(); return;}
		FileFormat ff=FileFormat.testInput(path, null, false);
		if(ff==null){setFallbackCF(); return;}
		ByteFile bf=ByteFile.makeByteFile(ff, 1);
		if(bf==null){setFallbackCF(); return;}

		java.util.ArrayList<double[]> multRows=new java.util.ArrayList<>();
		java.util.ArrayList<double[]> cvRows  =new java.util.ArrayList<>();
		java.util.ArrayList<Double>   avgList =new java.util.ArrayList<>();

		// Section IDs for state machine
		final int NONE=0, TIERAVG=1, MULT=2, CV=3;
		int section=NONE;

		for(byte[] raw=bf.nextLine(); raw!=null; raw=bf.nextLine()){
			if(raw.length==0){continue;}
			final String s=new String(raw).trim();
			if(s.isEmpty()){continue;}

			if(s.startsWith("#")){
				// Check for section headers
				if(s.startsWith("#TierAvg")){section=TIERAVG; continue;}
				else if(s.startsWith("#StateBlendMult")){section=MULT; continue;}
				else if(s.startsWith("#StateCV")){section=CV; continue;}
				else if(s.startsWith("#Tier")){continue;}// column header, skip
				else{section=NONE; continue;}// separator or other header
			}

			// Data row (no #)
			if(section==NONE){continue;}
			final String[] parts=s.split("\t");
			if(parts.length<2){continue;}
			final int tier;
			try{tier=Integer.parseInt(parts[0]);}
			catch(NumberFormatException e){continue;}

			if(section==TIERAVG){
				// Columns: tier, lin, geo, harm, obs
				if(parts.length<4){continue;}
				final double geo =Double.parseDouble(parts[2]);
				final double harm=Double.parseDouble(parts[3]);
				final double blend=(2.0*harm+geo)/3.0;
				while(avgList.size()<=tier){avgList.add(0.0);}
				avgList.set(tier, blend);
			}else if(section==MULT){
				final double[] row=new double[NUM_COMBINED];
				java.util.Arrays.fill(row, 1.0);
				for(int i=1; i<parts.length && (i-1)<NUM_COMBINED; i++){
					try{row[i-1]=Double.parseDouble(parts[i]);}
					catch(NumberFormatException e){/* leave 1.0 */}
				}
				while(multRows.size()<=tier){multRows.add(null);}
				multRows.set(tier, row);
			}else if(section==CV){
				final double[] row=new double[NUM_COMBINED];
				java.util.Arrays.fill(row, 1.0);
				for(int i=1; i<parts.length && (i-1)<NUM_COMBINED; i++){
					try{row[i-1]=Double.parseDouble(parts[i]);}
					catch(NumberFormatException e){/* leave 1.0 */}
				}
				while(cvRows.size()<=tier){cvRows.add(null);}
				cvRows.set(tier, row);
			}
		}
		bf.close();

		tierAvg=new double[avgList.size()];
		for(int i=0; i<avgList.size(); i++){tierAvg[i]=avgList.get(i);}

		CF_TABLE_TIERS=multRows.size();
		cfTable=new double[CF_TABLE_TIERS][NUM_COMBINED];
		for(int t=0; t<CF_TABLE_TIERS; t++){
			if(multRows.get(t)!=null){cfTable[t]=multRows.get(t);}
			else{java.util.Arrays.fill(cfTable[t], 1.0);}
		}
		cvTable=new double[cvRows.size()][NUM_COMBINED];
		for(int t=0; t<cvRows.size(); t++){
			if(cvRows.get(t)!=null){cvTable[t]=cvRows.get(t);}
			else{java.util.Arrays.fill(cvTable[t], 1.0);}
		}

		// Build meanCorrTable: normalized 1/cfTable per tier
		meanCorrTable=new double[CF_TABLE_TIERS][NUM_COMBINED];
		for(int t=0; t<CF_TABLE_TIERS; t++){
			double invSum=0;
			for(int ch=0; ch<NUM_COMBINED; ch++){
				final double v=cfTable[t][ch];
				invSum+=1.0/(v>0.001 ? v : 0.001);
			}
			final double invMean=invSum/NUM_COMBINED;
			for(int ch=0; ch<NUM_COMBINED; ch++){
				final double v=cfTable[t][ch];
				meanCorrTable[t][ch]=(1.0/(v>0.001 ? v : 0.001))/invMean;
			}
		}

		System.err.println("Loaded TTLL CF table: "+CF_TABLE_TIERS+" tiers, "
			+tierAvg.length+" tier averages, "+NUM_COMBINED+" states"
			+", "+cvTable.length+" CV tiers");
	}

	public static synchronized void loadCardCFTable(){
		if(cardCFKeys!=null){return;}
		String path=Data.findPath("?"+CF_FILE);
		if(path==null){return;}
		FileFormat ff=FileFormat.testInput(path, null, false);
		if(ff==null){return;}
		ByteFile bf=ByteFile.makeByteFile(ff, 1);
		if(bf==null){return;}

		structures.FloatList keys=new structures.FloatList();
		structures.FloatList vals=new structures.FloatList();

		for(byte[] raw=bf.nextLine(); raw!=null; raw=bf.nextLine()){
			if(raw.length==0){continue;}
			String s=new String(raw).trim();
			if(s.startsWith("#Buckets")){
				String[] parts=s.split("\t");
				if(parts.length>=2){cardCFBuckets=Integer.parseInt(parts[1].trim());}
				continue;
			}
			if(s.startsWith("#") || !Character.isDigit(s.charAt(0))){continue;}
			String[] parts=s.split("\t");
			if(parts.length<2){continue;}
			float key=Float.parseFloat(parts[0]);
			float cf=Float.parseFloat(parts[1]);
			keys.add(key);
			vals.add(cf);
		}
		bf.close();

		cardCFKeys=keys.toArray();
		cardCFValues=vals.toArray();
		System.err.println("Loaded TTLL cardinality CF table: "+cardCFKeys.length
			+" entries, buckets="+cardCFBuckets);
	}

	private static void setFallbackCF(){
		cfTable=new double[1][NUM_COMBINED];
		java.util.Arrays.fill(cfTable[0], 1.0);
		meanCorrTable=new double[1][NUM_COMBINED];
		java.util.Arrays.fill(meanCorrTable[0], 1.0);
		tierAvg=new double[]{1.0};
		CF_TABLE_TIERS=1;
	}

	/** Stub: measure from preliminary CF table, then replace 1f with actual ratio. */
	@Override public float terminalMeanCF(){return 1f;}

	/** Stub: TTLL has 1-bit history. */
	@Override public float terminalMeanPlusCF(){return 1f;}
}
