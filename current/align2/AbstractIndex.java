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
	
	abstract float[] keyProbArray();
	abstract byte[] getBaseScoreArray(int len, int strand);
	abstract int[] getKeyScoreArray(int len, int strand);
	
	abstract int maxScore(int[] offsets, byte[] baseScores, int[] keyScores, int readlen, boolean useQuality);
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
	 * True: Slightly slower.<br>
	 * False: Faster, but may mask detection of some ambiguously mapping reads.
	 */
	static final boolean LIMIT_SUBSUMPTION_LENGTH_TO_2X=true;
	
	/** Not recommended if slow alignment is disabled.  Can conceal sites that should be marked as amiguous. */
	static final boolean SUBSUME_OVERLAPPING_SITES=false;
	
	static final boolean SHRINK_BEFORE_WALK=true;

	/** More accurate but uses chromosome arrays while mapping */
	static final boolean USE_EXTENDED_SCORE=true; //Calculate score more slowly by extending keys
	
	/** Even more accurate but even slower than normal extended score calculation.
	 * Scores are compatible with slow-aligned scores. */
	static final boolean USE_AFFINE_SCORE=true && USE_EXTENDED_SCORE; //Calculate score even more slowly

	
	public static final boolean RETAIN_BEST_SCORES=true;
	public static final boolean RETAIN_BEST_QCUTOFF=true;
	
	public static boolean QUIT_AFTER_TWO_PERFECTS=true;
	static final boolean DYNAMICALLY_TRIM_LOW_SCORES=true;

	
	static final boolean REMOVE_CLUMPY=true; //Remove keys like AAAAAA or GCGCGC that self-overlap and thus occur in clumps
	
	
	/** If no hits are found, search again with slower parameters (less of genome excluded) */
	static final boolean DOUBLE_SEARCH_NO_HIT=false;
	/** Only this fraction of the originally removed genome fraction (FRACTION_GENOME_TO_EXCLUDE)
	 * is removed for the second pass */
	static final float DOUBLE_SEARCH_THRESH_MULT=0.25f; //Must be less than 1.
	
	static boolean PERFECTMODE=false;
	static boolean SEMIPERFECTMODE=false;
	
	static boolean REMOVE_FREQUENT_GENOME_FRACTION=true;//Default true; false is more accurate
	static boolean TRIM_BY_GREEDY=true;//default: true
	
	/** Ignore longest site list(s) when doing a slow walk. */
	static final boolean TRIM_LONG_HIT_LISTS=false; //Increases speed with tiny loss of accuracy.  Default: true for clean or synthetic, false for noisy real data
	
	public static int MIN_APPROX_HITS_TO_KEEP=1; //Default 2 for skimmer, 1 otherwise, min 1; lower is more accurate
	
	
	public static final boolean TRIM_BY_TOTAL_SITE_COUNT=false; //default: false
	/** Length histogram index of maximum average hit list length to use.
	 * The max number of sites to search is calculated by (#keys)*(lengthHistogram[chrom][MAX_AVERAGE_SITES_TO_SEARCH]).
	 * Then, while the actual number of sites exceeds this, the longest hit list should be removed.
	 */
	
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
	/** Site must have at least this many contiguous matches */
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
	
	static final int calcListSizeBonus(int[] array){
		if(array==null || array.length>LIST_SIZE_BONUS.length-1){return 0;}
		return LIST_SIZE_BONUS[array.length];
	}
	
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
