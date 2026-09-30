package align2;

import java.util.ArrayList;

import shared.Shared;
import stream.SiteScore;

/**
 * Shared index configuration plus per-worker search state for BBIndex variants.
 * The block index, analyzed key counts, chromosome bounds, and mode flags are
 * process-wide static state. A worker owns its AbstractIndex instance, scoring
 * scratch, aligner, and counters. Loading, analysis, and clearing shared state
 * must finish outside active worker searches; this is not a separate-reference
 * concurrent-mapper API. Counters measure operations, not CPU time.
 *
 * @author Brian Bushnell
 * @date Oct 15, 2013
 */
public abstract class AbstractIndex{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs an AbstractIndex with specified indexing parameters.
	 * Initializes key length, key space size, scoring parameters, and chromosome bounds.
	 *
	 * @param keylen Length of k-mers for indexing (typically 10-15 bases)
	 * @param kfilter Minimum number of contiguous matches required
	 * @param pointsMatch Points awarded per matching base
	 * @param minChrom_ Minimum chromosome number to process
	 * @param maxChrom_ Maximum chromosome number to process
	 * @param msa_ Multi-state aligner instance for alignment operations
	 */
	AbstractIndex(final int keylen, final int kfilter, final int pointsMatch,
			final int minChrom_, final int maxChrom_, final MSA msa_){
		assert(keylen>0 && keylen<16) : "AbstractMapper accepts k=1..15; larger shifts wrap the int keyspace: k="+keylen;
		KEYLEN=keylen;
		KEYSPACE=1<<(2*KEYLEN);
		BASE_KEY_HIT_SCORE=pointsMatch*KEYLEN;
		KFILTER=kfilter;
		msa=msa_;

		minChrom=minChrom_;
		maxChrom=maxChrom_;
		assert(minChrom==MINCHROM) : "Worker minimum must match the loaded shared index: "+minChrom+" versus "+MINCHROM;
		assert(maxChrom==MAXCHROM) : "Worker maximum must match the loaded shared index: "+maxChrom+" versus "+MAXCHROM;
		assert(minChrom<=maxChrom) : "Index chromosome range is empty or reversed: "+minChrom+".."+maxChrom;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the analyzed usable count for a key and its reverse complement.
	 * Analysis combines strands (palindromes once) and can zero filtered/clumpy keys;
	 * this is not an unfiltered genome-occurrence count. Without COUNTS, the legacy
	 * fallback queries only index[0], respecting Block's removed-list marker.
	 * @param key The k-mer encoded as an integer
	 * @return Cached analyzed count, or the single-block fallback count
	 */
	final int count(final int key){
//		assert(false);
		if(COUNTS!=null){return COUNTS[key];} //TODO: Benchmark speed and memory usage with counts=null.  Probably only works for single-block genomes.
//		assert(false);
		//TODO: Probable limitation - COUNTS=null does not aggregate multiple blocks
		// and assumes slot0 is populated. Do not enable this fallback generally.
		final Block b=index[0];
		final int rkey=KeyRing.reverseComplementKey(key, KEYLEN);
		final int a=b.length(key);
		return key==rkey ? a : a+b.length(rkey);
	}

	/** Adds nonnegative analyzed hit counts, saturating before any int narrowing.
	 * BBIndex.analyzeIndex uses this for cross-block and reverse-complement totals. */
	static int addCounts(final int a, final int b){
		assert(a>=0 && b>=0) : "Index hit-list lengths and accumulated counts must be nonnegative: "+a+", "+b;
		return (int)Math.min(Integer.MAX_VALUE, (long)a+b);
	}

	/** Tests inclusive intervals; touching at one coordinate counts as overlap. */
	static final boolean overlap(final int a1, final int b1, final int a2, final int b2){
		assert(a1<=b1 && a2<=b2) : "Inclusive intervals require ordered endpoints: "+a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}

	/**
	 * Tests whether the first inclusive interval is contained within the second.
	 *
	 * @param a1 Start position of first interval
	 * @param b1 End position of first interval
	 * @param a2 Start position of second interval
	 * @param b2 End position of second interval
	 * @return true if first interval is within second interval, false otherwise
	 */
	static final boolean isWithin(final int a1, final int b1, final int a2, final int b2){
		assert(a1<=b1 && a2<=b2) : "Inclusive intervals require ordered endpoints: "+a1+", "+b1+", "+a2+", "+b2;
		return a1>=a2 && b1<=b2;
	}


	/**
	 * Calculates a score term based on the distance between the leftmost and rightmost
	 * perfect matches sharing the same genomic location. Increases score with alignment span.
	 *
	 * @param locs Array of genomic locations for k-mer hits
	 * @param centerIndex Index of the leftmost perfect match in the arrays
	 * @param offsets Array of query offsets corresponding to the hits
	 * Only the offsets.length active entries are searched, even if locs has extra
	 * scratch capacity. centerIndex is assumed to be the first contributing match;
	 * this method finds the last equal location and does not search left of center.
	 * @return Query-offset span in bases; zero if no later hit shares the location
	 */
	final static int scoreY(final int[] locs, final int centerIndex, final int[] offsets){
		assert(centerIndex>=0 && centerIndex<offsets.length && locs.length>=offsets.length) :
				"quickScore must supply an anchor inside the active seed offsets: center="+centerIndex+
				", offsets="+offsets.length+", locations="+locs.length;
		final int center=locs[centerIndex];
//		int rightIndex=centerIndex;
//		for(int i=centerIndex; i<offsets.length; i++){
//			if(locs[i]==center){
//				rightIndex=i;
//			}
//		}

		int rightIndex=-1;
		for(int i=offsets.length-1; rightIndex<centerIndex; i--){
			if(locs[i]==center){
				rightIndex=i;
			}
		}

		//Assumed to not be necessary.
//		for(int i=0; i<centerIndex; i++){
//			if(locs[i]==center){
//				centerIndex=i;
//			}
//		}

		return offsets[rightIndex]-offsets[centerIndex];
	}

	/** Returns worker scratch for per-key probabilities; callers overwrite active entries. */
	abstract float[] keyProbArray();
	/** Returns strand-specific base-score scratch; may reuse nonzero prior contents. */
	abstract byte[] getBaseScoreArray(int len, int strand);
	/** Returns strand-specific key-score scratch; may reuse nonzero prior contents. */
	abstract int[] getKeyScoreArray(int len, int strand);

	abstract int maxScore(int[] offsets, byte[] baseScores, int[] keyScores, int readlen, boolean useQuality);
	/**
	 * Performs advanced site finding and scoring for read alignment.
	 * Main method for identifying potential alignment sites with detailed scoring.
	 *
	 * @param basesP Forward strand bases
	 * @param basesM Reverse strand bases
	 * @param qual Quality scores for the read
	 * @param baseScoresP Base-level scoring array for positive strand
	 * @param keyScoresP K-mer-level scoring array for positive strand
	 * @param offsets K-mer position offsets within the read
	 * @param id Unique identifier for the read
	 * @return List of potential alignment sites with scores
	 */
	public abstract ArrayList<SiteScore> findAdvanced(byte[] basesP, byte[] basesM, byte[] qual, byte[] baseScoresP, int[] keyScoresP, int[] offsets, long id);

	long callsToScore=0;
	long callsToExtendScore=0;
	long initialKeys=0;
	long initialKeyIterations=0;
	long initialKeys2=0;
	long initialKeyIterations2=0;
	long usedKeys=0;
	long usedKeyIterations=0;

	static final int HIT_HIST_LEN=40;
	final long[] hist_hits=new long[HIT_HIST_LEN+1];
	final long[] hist_hits_score=new long[HIT_HIST_LEN+1];
	final long[] hist_hits_extend=new long[HIT_HIST_LEN+1];

	final int minChrom;
	final int maxChrom;

	static int MINCHROM=1;
	static int MAXCHROM=Integer.MAX_VALUE;

	static final boolean SUBSUME_SAME_START_SITES=true; //Not recommended if slow alignment is disabled.
	static final boolean SUBSUME_SAME_STOP_SITES=true; //Not recommended if slow alignment is disabled.

	/**
	 * Whether to limit site subsumption to alignments within 2x length difference
	 */
	static final boolean LIMIT_SUBSUMPTION_LENGTH_TO_2X=true;

	/** Whether to merge overlapping alignment sites */
	static final boolean SUBSUME_OVERLAPPING_SITES=false;

	static final boolean SHRINK_BEFORE_WALK=true;

	/** Whether to use extended scoring calculation for higher accuracy */
	static final boolean USE_EXTENDED_SCORE=true; //Calculate score more slowly by extending keys

	/** Whether to use affine gap penalty scoring for alignment compatibility */
	static final boolean USE_AFFINE_SCORE=true && USE_EXTENDED_SCORE; //Calculate score even more slowly


	public static final boolean RETAIN_BEST_SCORES=true;
	public static final boolean RETAIN_BEST_QCUTOFF=true;

	public static boolean QUIT_AFTER_TWO_PERFECTS=true;
	static final boolean DYNAMICALLY_TRIM_LOW_SCORES=true;


	static final boolean REMOVE_CLUMPY=true; //Remove keys like AAAAAA or GCGCGC that self-overlap and thus occur in clumps


	/**
	 * Whether to perform second search pass with relaxed parameters if no hits found
	 */
	static final boolean DOUBLE_SEARCH_NO_HIT=false;
	/** Multiplier for genome exclusion fraction during second search pass */
	static final float DOUBLE_SEARCH_THRESH_MULT=0.25f; //Must be less than 1.

	static boolean PERFECTMODE=false;
	static boolean SEMIPERFECTMODE=false;

	static boolean REMOVE_FREQUENT_GENOME_FRACTION=true;//Default true; false is more accurate
	static boolean TRIM_BY_GREEDY=true;//default: true

	/**
	 * Whether to ignore longest hit lists during alignment walks for performance
	 */
	static final boolean TRIM_LONG_HIT_LISTS=false; //Increases speed with tiny loss of accuracy.  Default: true for clean or synthetic, false for noisy real data

	public static int MIN_APPROX_HITS_TO_KEEP=1; //Default 2 for skimmer, 1 otherwise, min 1; lower is more accurate


	public static final boolean TRIM_BY_TOTAL_SITE_COUNT=false; //default: false

	/** Maximum analyzed hit-list count admitted by the ordinary key filter. */
	static int MAX_USABLE_LENGTH=Integer.MAX_VALUE;
	/** More permissive hit-list count limit for relaxed filtering. */
	static int MAX_USABLE_LENGTH2=Integer.MAX_VALUE;


	/** Releases shared index/count/histogram arrays; does not reset configuration or
	 * worker counters. Call only after workers stop referencing the shared index. */
	public static void clear(){
		index=null;
		lengthHistogram=null;
		COUNTS=null;
	}

	static Block[] index;
	static int[] lengthHistogram=null;
	static int[] COUNTS=null;

	final int KEYLEN; //default 12, suggested 10 ~ 13, max 15; bigger is faster but uses more RAM
	final int KEYSPACE;
	/** Minimum number of contiguous matches required for site consideration */
	final int KFILTER;
	final MSA msa;
	final int BASE_KEY_HIT_SCORE;


	boolean verbose=false;
	static boolean verbose2=false;

	static boolean SLOW=false;
	static boolean VSLOW=false;

	static int NUM_CHROM_BITS=3;
	static int CHROMS_PER_BLOCK=(1<<(NUM_CHROM_BITS));

	static final int MINGAP=Shared.MINGAP;
	static final int MINGAP2=(MINGAP+128); //Depends on read length...

	static boolean USE_CAMELWALK=false;

	static final boolean ADD_LIST_SIZE_BONUS=false;
	static final byte[] LIST_SIZE_BONUS=new byte[100];

	public static boolean GENERATE_KEY_SCORES_FROM_QUALITY=true; //True: Much faster and more accurate.
	public static boolean GENERATE_BASE_SCORES_FROM_QUALITY=true; //True: Faster, and at least as accurate.

	/**
	 * Calculates a scoring bonus based on array length.
	 * Used for adjusting alignment scores based on hit list sizes.
	 * @param array Array to calculate bonus for
	 * @return Bonus score based on array length
	 */
	static final int calcListSizeBonus(final int[] array){
		if(array==null || array.length>LIST_SIZE_BONUS.length-1){return 0;}
		return LIST_SIZE_BONUS[array.length];
	}

	/**
	 * Calculates a scoring bonus based on list size.
	 * Used for adjusting alignment scores based on hit list sizes.
	 * @param size Size of the list
	 * @return Bonus score based on list size
	 */
	static final int calcListSizeBonus(final int size){
		//TODO: Possible bug [align2/AbstractIndex#001] - no lower-bound guard: a negative `size` indexes
		//LIST_SIZE_BONUS[neg] -> AIOOBE. The int[] overload above is structurally safe (array.length>=0). LOW/latent:
		//all 15 callers are gated by ADD_LIST_SIZE_BONUS (compile-time false) AND pass hit-list sizes (>=0). Brian's call.
		if(size>LIST_SIZE_BONUS.length-1){return 0;}
		return LIST_SIZE_BONUS[size];
	}

	static{
		final int len=LIST_SIZE_BONUS.length;
//		for(int i=1; i<len; i++){
//			int x=(int)((len/(Math.sqrt(i)))/5)-1;
//			LIST_SIZE_BONUS[i]=(byte)(x/2);
//		}
		LIST_SIZE_BONUS[0]=3;
		LIST_SIZE_BONUS[1]=2;
		LIST_SIZE_BONUS[2]=1;
		LIST_SIZE_BONUS[len-1]=0;
//		System.err.println(Arrays.toString(LIST_SIZE_BONUS));
	}

}
