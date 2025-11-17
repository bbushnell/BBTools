package align2;

import java.util.ArrayList;

import shared.Shared;
import stream.SiteScore;

/**
 * @author Brian Bushnell
 * @date Oct 15, 2013
 *
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
	
	/**
	 * Tests whether two genomic intervals overlap.
	 *
	 * @param a1 Start position of first interval
	 * @param b1 End position of first interval
	 * @param a2 Start position of second interval
	 * @param b2 End position of second interval
	 * @return true if intervals overlap, false otherwise
	 */
	static final boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	/** Is (a1, b1) within (a2, b2) ? */
	static final boolean isWithin(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a1>=a2 && b1<=b2;
	}
	
	
	/** Generates a term that increases score with how far apart the two farthest perfect matches are.
	 * Assumes that the centerIndex corresponds to the leftmost perfect match. */
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
	
	/** Returns array of k-mer occurrence probabilities for scoring calculations */
	abstract float[] keyProbArray();
	/**
	 * Generates base-level scoring array for alignment evaluation.
	 * @param len Length of sequence to score
	 * @param strand Strand orientation (0 for forward, 1 for reverse)
	 * @return Array of base-level scores for alignment
	 */
	abstract byte[] getBaseScoreArray(int len, int strand);
	/**
	 * Generates k-mer-level scoring array for alignment evaluation.
	 * @param len Length of sequence to score
	 * @param strand Strand orientation (0 for forward, 1 for reverse)
	 * @return Array of k-mer-level scores for alignment
	 */
	abstract int[] getKeyScoreArray(int len, int strand);
	
	/**
	 * Calculates maximum possible alignment score for a read given scoring arrays.
	 *
	 * @param offsets Array of k-mer positions in the read
	 * @param baseScores Base-level scoring array
	 * @param keyScores K-mer-level scoring array
	 * @param readlen Length of the read sequence
	 * @param useQuality Whether to incorporate quality scores in calculation
	 * @return Maximum achievable alignment score
	 */
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
	
	/** Counter for number of calls to scoring methods */
	long callsToScore=0;
	/** Counter for number of calls to extended scoring methods */
	long callsToExtendScore=0;
	/** Counter for initial k-mer keys processed */
	long initialKeys=0;
	/** Counter for iterations during initial k-mer processing */
	long initialKeyIterations=0;
	/** Secondary counter for initial k-mer keys processed */
	long initialKeys2=0;
	/** Secondary counter for iterations during initial k-mer processing */
	long initialKeyIterations2=0;
	/** Counter for k-mer keys actually used in alignment */
	long usedKeys=0;
	/** Counter for iterations during used k-mer processing */
	long usedKeyIterations=0;
	
	/** Length of hit histogram arrays for performance tracking */
	static final int HIT_HIST_LEN=40;
	/** Histogram of hit counts for performance analysis */
	final long[] hist_hits=new long[HIT_HIST_LEN+1];
	/** Histogram of hit counts during scoring operations */
	final long[] hist_hits_score=new long[HIT_HIST_LEN+1];
	/** Histogram of hit counts during extended scoring operations */
	final long[] hist_hits_extend=new long[HIT_HIST_LEN+1];
	
	/** Minimum chromosome number to process during alignment */
	final int minChrom;
	/** Maximum chromosome number to process during alignment */
	final int maxChrom;
	
	/** Global minimum chromosome number for all index instances */
	static int MINCHROM=1;
	/** Global maximum chromosome number for all index instances */
	static int MAXCHROM=Integer.MAX_VALUE;

	/** Whether to merge alignment sites with identical start positions */
	static final boolean SUBSUME_SAME_START_SITES=true; //Not recommended if slow alignment is disabled.
	/** Whether to merge alignment sites with identical stop positions */
	static final boolean SUBSUME_SAME_STOP_SITES=true; //Not recommended if slow alignment is disabled.
	
	/**
	 * True: Slightly slower.<br>
	 * False: Faster, but may mask detection of some ambiguously mapping reads.
	 */
	static final boolean LIMIT_SUBSUMPTION_LENGTH_TO_2X=true;
	
	/** Not recommended if slow alignment is disabled.  Can conceal sites that should be marked as amiguous. */
	static final boolean SUBSUME_OVERLAPPING_SITES=false;
	
	/** Whether to shrink hit lists before performing alignment walks */
	static final boolean SHRINK_BEFORE_WALK=true;

	/** More accurate but uses chromosome arrays while mapping */
	static final boolean USE_EXTENDED_SCORE=true; //Calculate score more slowly by extending keys
	
	/** Even more accurate but even slower than normal extended score calculation.
	 * Scores are compatible with slow-aligned scores. */
	static final boolean USE_AFFINE_SCORE=true && USE_EXTENDED_SCORE; //Calculate score even more slowly

	
	/** Whether to retain best alignment scores during processing */
	public static final boolean RETAIN_BEST_SCORES=true;
	/** Whether to retain best quality cutoff values during processing */
	public static final boolean RETAIN_BEST_QCUTOFF=true;
	
	/** Whether to stop searching after finding two perfect alignment matches */
	public static boolean QUIT_AFTER_TWO_PERFECTS=true;
	/** Whether to dynamically remove low-scoring alignments during processing */
	static final boolean DYNAMICALLY_TRIM_LOW_SCORES=true;

	
	/**
	 * Whether to remove self-overlapping repetitive k-mers like AAAAAA or GCGCGC
	 */
	static final boolean REMOVE_CLUMPY=true; //Remove keys like AAAAAA or GCGCGC that self-overlap and thus occur in clumps
	
	
	/** If no hits are found, search again with slower parameters (less of genome excluded) */
	static final boolean DOUBLE_SEARCH_NO_HIT=false;
	/** Only this fraction of the originally removed genome fraction (FRACTION_GENOME_TO_EXCLUDE)
	 * is removed for the second pass */
	static final float DOUBLE_SEARCH_THRESH_MULT=0.25f; //Must be less than 1.
	
	/** Whether alignment operates in perfect match mode only */
	static boolean PERFECTMODE=false;
	/** Whether alignment operates in semi-perfect match mode */
	static boolean SEMIPERFECTMODE=false;
	
	/** Whether to remove frequently occurring genomic regions from consideration */
	static boolean REMOVE_FREQUENT_GENOME_FRACTION=true;//Default true; false is more accurate
	/** Whether to use greedy algorithm for trimming alignment candidates */
	static boolean TRIM_BY_GREEDY=true;//default: true
	
	/** Ignore longest site list(s) when doing a slow walk. */
	static final boolean TRIM_LONG_HIT_LISTS=false; //Increases speed with tiny loss of accuracy.  Default: true for clean or synthetic, false for noisy real data
	
	/** Minimum number of approximate hits required to retain an alignment site */
	public static int MIN_APPROX_HITS_TO_KEEP=1; //Default 2 for skimmer, 1 otherwise, min 1; lower is more accurate
	
	
	/** Whether to trim alignment candidates based on total site count */
	public static final boolean TRIM_BY_TOTAL_SITE_COUNT=false; //default: false
	/** Length histogram index of maximum average hit list length to use.
	 * The max number of sites to search is calculated by (#keys)*(lengthHistogram[chrom][MAX_AVERAGE_SITES_TO_SEARCH]).
	 * Then, while the actual number of sites exceeds this, the longest hit list should be removed.
	 */
	
	static int MAX_USABLE_LENGTH=Integer.MAX_VALUE;
	/** Secondary maximum sequence length threshold for processing */
	static int MAX_USABLE_LENGTH2=Integer.MAX_VALUE;

	
	/** Clears static data structures to free memory */
	public static void clear(){
		index=null;
		lengthHistogram=null;
		COUNTS=null;
	}
	
	/** Array of index blocks containing k-mer position data */
	static Block[] index;
	/** Histogram of sequence lengths for statistical analysis */
	static int[] lengthHistogram=null;
	/** Precomputed counts of k-mer occurrences for fast lookup */
	static int[] COUNTS=null;
	
	/** Length of k-mers used for indexing */
	final int KEYLEN; //default 12, suggested 10 ~ 13, max 15; bigger is faster but uses more RAM
	/** Total number of possible k-mer keys (4^KEYLEN) */
	final int KEYSPACE;
	/** Site must have at least this many contiguous matches */
	final int KFILTER;
	/** Multi-state aligner instance for performing detailed alignments */
	final MSA msa;
	/** Base score value for a k-mer hit match */
	final int BASE_KEY_HIT_SCORE;
	
	
	/** Whether to enable verbose debug output for this index instance */
	boolean verbose=false;
	/** Global verbose debug flag for all index instances */
	static boolean verbose2=false;

	/** Whether to enable slower but more thorough alignment modes */
	static boolean SLOW=false;
	/** Whether to enable very slow but most thorough alignment modes */
	static boolean VSLOW=false;
	
	/** Number of bits used to encode chromosome information */
	static int NUM_CHROM_BITS=3;
	/** Number of chromosomes that can be stored per index block */
	static int CHROMS_PER_BLOCK=(1<<(NUM_CHROM_BITS));

	/** Minimum gap size for alignment processing */
	static final int MINGAP=Shared.MINGAP;
	/** Secondary minimum gap size threshold adjusted for read length */
	static final int MINGAP2=(MINGAP+128); //Depends on read length...
	
	/** Whether to use camel walk algorithm for alignment extension */
	static boolean USE_CAMELWALK=false;
	
	/** Whether to add scoring bonuses based on hit list sizes */
	static final boolean ADD_LIST_SIZE_BONUS=false;
	/** Array of bonus scores indexed by list size */
	static final byte[] LIST_SIZE_BONUS=new byte[100];
	
	/**
	 * Whether to generate k-mer scores from quality values for improved accuracy
	 */
	public static boolean GENERATE_KEY_SCORES_FROM_QUALITY=true; //True: Much faster and more accurate.
	/** Whether to generate base scores from quality values for improved accuracy */
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
