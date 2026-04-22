package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * UltraDynamicLogLog6: int-packed 6-bit register variant with ErtlULL encoding and 2-bit history.
 * <p>
 * Packs 5 six-bit registers per 32-bit int (30 bits used, 2 wasted):
 * <ul>
 * <li>At 2048 buckets: ceil(2048/5) = 410 ints = 1,640 bytes vs 2,048 bytes for byte[]
 *     (20% memory savings).
 * </ul>
 * Register encoding (ErtlULL format, 6-bit):
 * <ul>
 * <li>Bits [5:2]: nlzPart — absolute NLZ level above (globalNLZ+1), offset by HISTORY_MARGIN.
 * <li>Bits [1:0]: histPattern — 2-bit history of sub-floor observations.
 * <li>Special: reg==0 means "at floor level" (absNlz = globalNLZ).
 * </ul>
 * Absolute NLZ = nlzPart - HISTORY_MARGIN + (globalNLZ+1).
 * HISTORY_MARGIN=2: registers at nlzPart in {0,1,2} are at or below floor level.
 * <p>
 * On observation: unpack old register to hash-prefix bits, set the bit at
 * relNlz+HISTORY_MARGIN, repack. New register replaces old if strictly larger.
 * <p>
 * Functionally identical to the byte-array UltraDynamicLogLog; exists to demonstrate
 * int-packed memory savings for timing comparisons.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class UltraDynamicLogLog6 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor: 2048 buckets, k=31. */
	public UltraDynamicLogLog6(){this(2048, 31, -1, 0);}

	/** Construct from parsed command-line arguments. */
	UltraDynamicLogLog6(Parser p){
		super(p);
		packedLen=(buckets+4)/5;
		packed=new int[packedLen];
		floorCount=buckets;
	}

	/**
	 * Full constructor.
	 * @param buckets_ Number of buckets
	 * @param k_ Hash prefix length
	 * @param seed Random seed (-1 for default)
	 * @param minProb_ Minimum probability threshold
	 */
	UltraDynamicLogLog6(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		packedLen=(buckets+4)/5;
		packed=new int[packedLen];
		floorCount=buckets;
	}

	/** Create an independent copy with a fresh seed. */
	@Override public UltraDynamicLogLog6 copy(){return new UltraDynamicLogLog6(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read 6-bit register i from packed int array. */
	int getReg(int i){
		return (packed[i/5]>>>((i%5)*6))&0x3F;
	}

	/** Write 6-bit register i into packed int array. */
	private void setReg(int i, int val){
		final int idx=i/5;
		final int shift=(i%5)*6;
		packed[idx]=(packed[idx] & ~(0x3F<<shift)) | ((val&0x3F)<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Build NLZ histogram and per-bucket packed state for CardStats.
	 * Converts each bucket's register to absolute NLZ via ErtlULL encoding;
	 * floor-level and sub-margin buckets map to absNlz = globalNLZ.
	 * @return CardStats with all estimator values computed
	 */
	private CardStats summarize(){
		final int[] nlzCounts=new int[66];
		final char[] packedBuckets=new char[buckets];
		int filledCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg==0){
				// Floor-level bucket: absNlz = globalNLZ
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
					// Sub-margin bucket: at floor level, history preserved
					if(globalNLZ>=0 && globalNLZ<64){
						nlzCounts[globalNLZ+1]++;
						filledCount++;
						packedBuckets[i]=(char)(((globalNLZ+1)<<2)|histPattern);
					}
				}
			}
		}
		nlzCounts[0]=buckets-filledCount;
		return new CardStats(packedBuckets, nlzCounts, 0, 2, 0, 0,
				buckets, microIndex, added, CF_MATRIX, CF_BUCKETS, 0,
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
		add((UltraDynamicLogLog6)log);
	}

	/**
	 * Merge another UDLL6 into this one via ErtlULL hash-prefix union.
	 * Shifts each side's hash-prefix right by the delta between their globalNLZ
	 * and the merged globalNLZ before ORing, preserving history alignment.
	 */
	public void add(UltraDynamicLogLog6 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		final int newGlobalNLZ=Math.max(globalNLZ, log.globalNLZ);
		final int deltaA=newGlobalNLZ-globalNLZ;
		final int deltaB=newGlobalNLZ-log.globalNLZ;
		for(int i=0; i<buckets; i++){
			final int rA=getReg(i);
			final int rB=log.getReg(i);
			long hpA=(rA==0) ? 0 : ErtlULL.unpack((byte)rA)>>>deltaA;
			long hpB=(rB==0) ? 0 : ErtlULL.unpack((byte)rB)>>>deltaB;
			long merged=hpA|hpB;
			if(merged==0){
				setReg(i, 0);
			}else{
				setReg(i, Math.min(ErtlULL.pack(merged)&0xFF, MAX_REGISTER));
			}
		}
		globalNLZ=newGlobalNLZ;
		filledBuckets=0;
		floorCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg>0){filledBuckets++;}
			if(reg==0 || (reg>>>2)<=HISTORY_MARGIN){floorCount++;}
		}
		while(floorCount==0 && globalNLZ<wordlen){
			globalNLZ++;
			floorCount=countAndDecrement();
		}
		int exitThreshold=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
		eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
	}

	/**
	 * Hash a value and store it in the appropriate bucket.
	 * Pipeline: hash → micro-cardinality sketch → early-exit mask →
	 * bucket index → NLZ → ErtlULL pack/unpack → register update → floor check.
	 */
	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);

		if(Long.compareUnsigned(key, eeMask)>0){return;}
		branch1++;

		final int idx=(int)(key&bucketMask);
		final int nlz=Long.numberOfLeadingZeros(key);
		final int relNlz=nlz-(globalNLZ+1);

		final int bitPos=relNlz+HISTORY_MARGIN;
		if(bitPos<0||bitPos>=64){return;}

		final int oldReg=getReg(idx);
		long hashPrefix=ErtlULL.unpack((byte)oldReg);
		hashPrefix|=1L<<bitPos;
		final int rawReg=ErtlULL.pack(hashPrefix)&0xFF;
		final int newReg;
		if(rawReg>MAX_REGISTER){
			final int sub=(int)((hashPrefix>>>13)&3); // 2-bit history from hash-prefix bits 13-14
			newReg=(15<<2)|sub;
		}else{
			newReg=rawReg;
		}
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
				int exitThreshold=Math.max(0, (globalNLZ+1)-HISTORY_MARGIN);
				eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
			}
		}
	}

	/**
	 * Decrement all registers after a global floor advance.
	 * Registers at reg=0 remain 0 (they are still floor-level).
	 * @return New floorCount (buckets at or below new floor level)
	 */
	private int countAndDecrement(){
		int newFloorCount=0;
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg==0){newFloorCount++; continue;}
			final int newReg=reg-4; // decrement nlzPart by 1 (bits [5:2]), history unchanged
			setReg(i, newReg);
			if((newReg>>>2)<=HISTORY_MARGIN){newFloorCount++;}
		}
		return newFloorCount;
	}

	/**
	 * Compute FGRA (Flajolet-Girard-Rémy-Adler) estimate from current register state.
	 * Converts int-packed 6-bit registers to byte[] ErtlULL format for ErtlULL.fgraEstimateStatic().
	 * Only computed when CALC_FGRA is true.
	 */
	public double fgraEstimate(){
		if(!CALC_FGRA){return 0;}
		final int p=bucketBits;
		final int regOffset=4*((globalNLZ+1)+p-1-HISTORY_MARGIN);
		final byte[] ertlRegs=new byte[buckets];
		for(int i=0; i<buckets; i++){
			final int reg=getReg(i);
			if(reg==0){
				ertlRegs[i]=0;
			}else{
				ertlRegs[i]=(byte)Math.min(reg+regOffset, 255);
			}
		}
		return ErtlULL.fgraEstimateStatic(ertlRegs, p);
	}

	/** Number of buckets with reg > 0. */
	public int filledBuckets(){return filledBuckets;}
	/** Fraction of buckets that are occupied. */
	public double occupancy(){return (double)filledBuckets/buckets;}
	/** Compatibility accessor: returns globalNLZ+1 to match legacy minZeros convention. */
	public int getMinZeros(){return globalNLZ+1;}

	/** Not used; CF correction handled via CF_MATRIX. */
	@Override public final float[] compensationFactorLogBucketsArray(){return null;}

	/**
	 * Compute all estimator values and return as a legacy-format array.
	 * Caches the CardStats for retrieval by consumeLastSummarized().
	 */
	@Override
	public double[] rawEstimates(){
		lastSummarized=summarize();
		final double hybridEst=lastSummarized.hybridDLL();
		return AbstractCardStats.buildLegacyArray(lastSummarized, hybridEst);
	}

	/** Public accessor for register value (for testing). */
	public int getRegPublic(int i){return getReg(i);}
	/** Alias for getRegPublic. */
	public int getRegister(int i){return getReg(i);}

	/** Legacy array format: {ldlc, dlcSbs, hc, lcMin, fgra, hll, meanH, hybridPlus2}. */
	public double[] udlcEstimate(){
		final CardStats s=(lastSummarized!=null) ? lastSummarized : summarize();
		return new double[]{s.ldlc(), s.dlcSbs(), s.hc(), s.lcMin(),
			fgraEstimate(), s.hllRaw(), s.meanHistCF(), s.hybridPlus2()};
	}

	/** Debug: print all register values to stderr. */
	public void printRegisters(){
		for(int i=0; i<buckets; i++){System.err.print(getReg(i)+" ");}
		System.err.println();
	}

	/** Memory used by the packed array in bytes. */
	public int packedBytes(){return packedLen*4;}

	/** Return and clear the cached CardStats. Used by calibration drivers. */
	public CardStats consumeLastSummarized(){
		final CardStats s=lastSummarized;
		lastSummarized=null;
		return s;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed storage: 5 × 6-bit registers per 32-bit int. */
	final int[] packed;
	/** Number of packed ints: ceil(buckets/5). */
	private final int packedLen;
	/** Global floor NLZ minus 1: -1 means nothing seen; >=0 means all buckets have absNlz >= globalNLZ+1. */
	private int globalNLZ=-1;
	/** Count of buckets at or below floor level (nlzPart <= HISTORY_MARGIN). */
	private int floorCount;
	/** Early-exit mask: filters hashes below current floor minus HISTORY_MARGIN. */
	private long eeMask=-1L;
	/** Count of buckets with reg > 0. */
	private int filledBuckets=0;
	/** Cached CardStats from last rawEstimates() call. */
	private CardStats lastSummarized;

	/** Branch counters for performance profiling. */
	public long branch1=0, branch2=0;
	/** Rate of hashes passing the eeMask filter (per added). */
	public double branch1Rate(){return branch1/(double)Math.max(1, added);}
	/** Rate of register updates among hashes that passed the eeMask filter. */
	public double branch2Rate(){return branch2/(double)Math.max(1, branch1);}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum NLZ value (64-bit hash). */
	private static final int wordlen=64;
	/** eeMask relaxation: accept hashes this many NLZ levels below floor for history. */
	static final int HISTORY_MARGIN=2;
	/** Maximum valid 6-bit register value (63). */
	static final int MAX_REGISTER=63;

	/** When true, clamp register at MAX_REGISTER instead of wrapping. */
	public static boolean SATURATE_ON_OVERFLOW=true;
	/** When true, fgraEstimate() is calculated. Default false to avoid per-call overhead. */
	public static boolean CALC_FGRA=false;

	/** Auto-loaded UDLL6 CF table resource path. */
	public static final String CF_FILE="?cardinalityCorrectionUDLL6.tsv.gz";
	/** Bucket count the CF_MATRIX was generated for. */
	private static int CF_BUCKETS=2048;
	/** Per-cardinality correction factor table, set by initializeCF or setCFMatrix. */
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);

	/** Load the CF table from the resource file, scaled to the given bucket count. */
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	/** Set the CF matrix directly (used by CardinalityParser after global load). */
	public static void setCFMatrix(float[][] matrix, int buckets){
		CF_MATRIX=matrix; CF_BUCKETS=buckets;
	}

	/** Asymptotic meanRaw/trueCard ratio, measured 128k ddls maxmult=4096 (Apr 13 2026). */
	@Override public float terminalMeanCF(){return 1.386726f;}

	/** Asymptotic Mean+H ratio (2-bit history), measured 128k ddls maxmult=4096 (Apr 13 2026). */
	@Override public float terminalMeanPlusCF(){return 0.856020f;}

}
