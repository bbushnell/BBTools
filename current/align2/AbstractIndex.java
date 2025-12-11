package align2;

import java.util.ArrayList;

import shared.Shared;
import stream.SiteScore;

/**
 * Abstract base class for sequence indexing in BBTools alignment framework.
 * Provides common functionality for k-mer based indexing with configurable key lengths and filtering.
 * Supports genomic interval operations, quality filtering thresholds, and performance monitoring.
 * Forms the foundation for BBIndex variants used in sequence alignment and mapping.
 *
 * @author Brian Bushnell
 * @date Oct 15, 2013
 */
public abstract class AbstractIndex {
	
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
	AbstractIndex(int keylen, int kfilter, int pointsMatch, int minChrom_, int maxChrom_, MSA msa_){
		KEYLEN=keylen;
		KEYSPACE=1<<(2*KEYLEN);
		BASE_KEY_HIT_SCORE=pointsMatch*KEYLEN;
		KFILTER=kfilter;
		msa=msa_;

		minChrom=minChrom_;
		maxChrom=maxChrom_;
		assert(minChrom==MINCHROM);
		assert(maxChrom==MAXCHROM);
		assert(minChrom<=maxChrom);
	}
	
	/**
	 * Returns the total count of occurrences for a k-mer key and its reverse complement.
	 * Uses precomputed counts array when available, otherwise queries the index block.
	 * @param key The k-mer encoded as an integer
	 * @return Total count of key and reverse complement occurrences
	 */
	final int count(int key){
//		assert(false);
		if(COUNTS!=null){return COUNTS[key];} //TODO: Benchmark speed and memory usage with counts=null.  Probably only works for single-block genomes.
//		assert(false);
		final Block b=index[0];
		final int rkey=KeyRing.reverseComplementKey(key, KEYLEN);
		int a=b.length(key);
		return key==rkey ? a : a+b.length(rkey);
	}
	
	static final boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	/**
	 * Tests whether the first interval is completely contained within the second interval.
	 *
	 * @param a1 Start position of first interval
	 * @param b1 End position of first interval
	 * @param a2 Start position of second interval
	 * @param b2 End position of second interval
	 * @return true if first interval is within second interval, false otherwise
	 */
	static final boolean isWithin(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a1>=a2 && b1<=b2;
	}
	
	
	/**
	 * Calculates a score term based on the distance between the leftmost and rightmost
	 * perfect matches sharing the same genomic location. Increases score with alignment span.
	 *
	 * @param locs Array of genomic locations for k-mer hits
	 * @param centerIndex Index of the leftmost perfect match in the arrays
	 * @param offsets Array of query offsets corresponding to the hits
	 * @return Distance between leftmost and rightmost matches at the center location
	 */
	final static int scoreY(int[] locs, int centerIndex, int offsets[]){
		int center=locs[centerIndex];
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
	
	abstract float[] keyProbArray();
	abstract byte[] getBaseScoreArray(int len, int strand);
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
	
	/** Maximum sequence length that can be effectively processed */
	static int MAX_USABLE_LENGTH=Integer.MAX_VALUE;
	static int MAX_USABLE_LENGTH2=Integer.MAX_VALUE;

	
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
	static final int calcListSizeBonus(int[] array){
		if(array==null || array.length>LIST_SIZE_BONUS.length-1){return 0;}
		return LIST_SIZE_BONUS[array.length];
	}
	
	/**
	 * Calculates a scoring bonus based on list size.
	 * Used for adjusting alignment scores based on hit list sizes.
	 * @param size Size of the list
	 * @return Bonus score based on list size
	 */
	static final int calcListSizeBonus(int size){
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
