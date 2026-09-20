package align2;

import java.util.ArrayList;
import java.util.Arrays;

import bloom.BloomFilter;
import dna.AminoAcid;
import dna.Data;
import jgi.CoveragePileup;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/**
 * Based on MapTestThread11f
 * 
 * @author Brian Bushnell
 * @date Dec 22, 2012
 *
 */
public final class BBMapThread extends AbstractMapThread{
	private final int hybridTipSearchCeiling;
	@Override
	protected int tipDeletionSearchRange(){
		if((hybridPair && index.maxIndel()==16000 && index.maxIndel2()==32000) ||
				(hybridMaxIndelConfig!=null && index.maxIndel()==hybridMaxIndelConfig.highPrimary &&
				index.maxIndel2()==hybridMaxIndelConfig.highSum)){
			return Tools.min(hybridTipSearchCeiling, index.maxIndel());
		}
		return super.tipDeletionSearchRange();
	}
	
	/** Number of columns in alignment matrices */
	static final int ALIGN_COLUMNS=BBIndex.ALIGN_COLUMNS;
	/** Number of rows in alignment matrices */
	static final int ALIGN_ROWS=601;
	
	

	/** Don't trim for local alignments unless at least this many bases will be clipped */
	private final int LOCAL_ALIGN_TIP_LENGTH=1;
	/** Range is 0-1; a lower number makes trimming more aggressive */
	private final float LOCAL_ALIGN_MATCH_POINT_RATIO=1f;
	
	/** Ratio of the points for a match of a single base needed to declare unambiguous.  1 SNP is currently about 2.57 */
	public final float CLEARZONE_RATIOP=1.6f; //default 1.3f, which makes read ambiguous if there is 1 N in an alternate site.
	/** Clearzone ratio for single-edit distance alignments */
	public final float CLEARZONE_RATIO1=2.0f;
	/** Clearzone ratio for moderate quality alignments */
	public final float CLEARZONE_RATIO1b=2.6f;
	/** Clearzone ratio for lower quality alignments */
	public final float CLEARZONE_RATIO1c=4.6f;
	/** Clearzone ratio for ambiguity penalties */
	public final float CLEARZONE_RATIO3=8.0f;
	/** Max allowed number of sites within 1 edit (excluding primary site) */
	public final int CLEARZONE_LIMIT1e=40; //Needs to be redone to assign a quality penalty rather than simply marking as ambiguous
	/** Clearzone threshold for perfect matches */
	public final int CLEARZONEP;
	/** Clearzone threshold for single-edit distance matches */
	public final int CLEARZONE1;
	/** Clearzone threshold for moderate quality matches */
	public final int CLEARZONE1b;
	/** Clearzone threshold for lower quality matches */
	public final int CLEARZONE1c;
	//public final int CLEARZONE1e;
	/** Clearzone threshold for ambiguity penalties */
	public final int CLEARZONE3;
	/** Inverse of CLEARZONE3 for efficient calculations */
	public final float INV_CLEARZONE3;
	/** Flat ratio component for clearzone 1b cutoff calculations */
	public final float CLEARZONE1b_CUTOFF_FLAT_RATIO=12;//3f;
	/** Flat cutoff value for clearzone 1b threshold */
	public final float CLEARZONE1b_CUTOFF_FLAT;
	/** Scale factor for clearzone 1b cutoff calculations */
	public final float CLEARZONE1b_CUTOFF_SCALE=0.97f;
	/** Flat ratio component for clearzone 1c cutoff calculations */
	public final float CLEARZONE1c_CUTOFF_FLAT_RATIO=26;//7f;
	/** Flat cutoff value for clearzone 1c threshold */
	public final float CLEARZONE1c_CUTOFF_FLAT;
	/** Scale factor for clearzone 1c cutoff calculations */
	public final float CLEARZONE1c_CUTOFF_SCALE=0.92f;
	
	/** BBMap index used for k-mer lookups and alignment */
	public final BBIndex index;
	
	
	/** Minimum number of sites to retain after trimming for single-end reads */
	private final int MIN_TRIM_SITES_TO_RETAIN_SINGLE=3;
	/** Minimum number of sites to retain after trimming for paired-end reads */
	private final int MIN_TRIM_SITES_TO_RETAIN_PAIRED=2;
	
	/**
	 * Warning method indicating EXPECTED_SITES is not valid for this class.
	 * Prints a deprecation warning when called.
	 * @param x Parameter value (ignored)
	 */
	public static void setExpectedSites(int x){
		System.err.println("Warning: EXPECTED_SITES is not valid for "+(new Object() { }.getClass().getEnclosingClass().getName()));
	}
	
	@Override
	public final int ALIGN_COLUMNS(){return ALIGN_COLUMNS;}
	@Override
	public final int ALIGN_ROWS(){return ALIGN_ROWS;}
	@Override
	public final int maxReadLength(){return ALIGN_ROWS-1;}
	@Override
	final AbstractIndex index(){return index;}
	@Override
	final int CLEARZONE1(){return CLEARZONE1;}

	/**
	 * Constructs a BBMapThread with comprehensive alignment and processing parameters.
	 * Initializes clearzone thresholds based on alignment mode and creates the BBIndex.
	 *
	 * @param cris_ Input stream for reading sequences
	 * @param keylen_ K-mer length for indexing
	 * @param pileup_ Coverage pileup tracker
	 * @param SMITH_WATERMAN_ Enable Smith-Waterman alignment
	 * @param THRESH_ Score threshold for alignments
	 * @param minChrom_ Minimum chromosome number to process
	 * @param maxChrom_ Maximum chromosome number to process
	 * @param keyDensity_ Target k-mer density
	 * @param maxKeyDensity_ Maximum k-mer density allowed
	 * @param minKeyDensity_ Minimum k-mer density required
	 * @param maxDesiredKeys_ Maximum number of keys desired
	 * @param REMOVE_DUPLICATE_BEST_ALIGNMENTS_ Remove duplicate top alignments
	 * @param SAVE_AMBIGUOUS_XY_ Save ambiguous alignment coordinates
	 * @param MINIMUM_ALIGNMENT_SCORE_RATIO_ Minimum score ratio for keeping alignments
	 * @param TRIM_LIST_ Enable trimming of site lists
	 * @param MAKE_MATCH_STRING_ Generate CIGAR/match strings
	 * @param QUICK_MATCH_STRINGS_ Use fast match string generation
	 * @param outStream_ Main output stream
	 * @param outStreamMapped_ Output stream for mapped reads
	 * @param outStreamUnmapped_ Output stream for unmapped reads
	 * @param outStreamBlack_ Output stream for blacklisted reads
	 * @param SLOW_ALIGN_PADDING_ Padding for slow alignment
	 * @param SLOW_RESCUE_PADDING_ Padding for rescue operations
	 * @param DONT_OUTPUT_UNMAPPED_READS_ Skip outputting unmapped reads
	 * @param DONT_OUTPUT_BLACKLISTED_READS_ Skip outputting blacklisted reads
	 * @param MAX_SITESCORES_TO_PRINT_ Maximum secondary alignments to output
	 * @param PRINT_SECONDARY_ALIGNMENTS_ Include secondary alignments in output
	 * @param REQUIRE_CORRECT_STRANDS_PAIRS_ Enforce proper strand orientation for pairs
	 * @param SAME_STRAND_PAIRS_ Allow same-strand pairs
	 * @param KILL_BAD_PAIRS_ Remove improperly paired reads
	 * @param RCOMP_MATE_ Reverse complement mate reads
	 * @param PERFECTMODE_ Enable perfect match mode
	 * @param SEMIPERFECTMODE_ Enable semi-perfect match mode
	 * @param FORBID_SELF_MAPPING_ Prevent reads from mapping to themselves
	 * @param TIP_DELETION_SEARCH_RANGE_ Range for searching tip deletions
	 * @param AMBIGUOUS_RANDOM_ Handle ambiguous reads randomly
	 * @param AMBIGUOUS_ALL_ Report all ambiguous alignments
	 * @param KFILTER_ K-mer filtering threshold
	 * @param IDFILTER_ Identity filtering threshold
	 * @param TRIM_LEFT_ Enable left-end trimming
	 * @param TRIM_RIGHT_ Enable right-end trimming
	 * @param UNTRIM_ Enable untrimming
	 * @param TRIM_QUAL_ Quality score for trimming
	 * @param TRIM_MIN_LEN_ Minimum length after trimming
	 * @param LOCAL_ALIGN_ Enable local alignment
	 * @param RESCUE_ Enable mate rescue
	 * @param STRICT_MAX_INDEL_ Enforce strict indel limits
	 * @param MSA_TYPE_ Type of multiple sequence aligner to use
	 * @param bloomFilter_ Bloom filter for contamination removal
	 */
	public BBMapThread(ConcurrentReadInputStream cris_, int keylen_,
			CoveragePileup pileup_, boolean SMITH_WATERMAN_, int THRESH_, int minChrom_,
			int maxChrom_, float keyDensity_, float maxKeyDensity_, float minKeyDensity_, int maxDesiredKeys_,
			boolean REMOVE_DUPLICATE_BEST_ALIGNMENTS_, boolean SAVE_AMBIGUOUS_XY_,
			float MINIMUM_ALIGNMENT_SCORE_RATIO_, boolean TRIM_LIST_, boolean MAKE_MATCH_STRING_, boolean QUICK_MATCH_STRINGS_,
			ConcurrentReadOutputStream outStream_, ConcurrentReadOutputStream outStreamMapped_, ConcurrentReadOutputStream outStreamUnmapped_, ConcurrentReadOutputStream outStreamBlack_,
			int SLOW_ALIGN_PADDING_, int SLOW_RESCUE_PADDING_, boolean DONT_OUTPUT_UNMAPPED_READS_, boolean DONT_OUTPUT_BLACKLISTED_READS_,
			int MAX_SITESCORES_TO_PRINT_, boolean PRINT_SECONDARY_ALIGNMENTS_,
			boolean REQUIRE_CORRECT_STRANDS_PAIRS_, boolean SAME_STRAND_PAIRS_, boolean KILL_BAD_PAIRS_, boolean RCOMP_MATE_,
			boolean PERFECTMODE_, boolean SEMIPERFECTMODE_, boolean FORBID_SELF_MAPPING_, int TIP_DELETION_SEARCH_RANGE_,
			boolean AMBIGUOUS_RANDOM_, boolean AMBIGUOUS_ALL_, int KFILTER_, float IDFILTER_, boolean TRIM_LEFT_, boolean TRIM_RIGHT_, boolean UNTRIM_, float TRIM_QUAL_, int TRIM_MIN_LEN_,
			boolean LOCAL_ALIGN_, boolean RESCUE_, boolean STRICT_MAX_INDEL_, String MSA_TYPE_, BloomFilter bloomFilter_){
		this(cris_, keylen_, pileup_, SMITH_WATERMAN_, THRESH_, minChrom_, maxChrom_, keyDensity_, maxKeyDensity_, minKeyDensity_, maxDesiredKeys_, REMOVE_DUPLICATE_BEST_ALIGNMENTS_, SAVE_AMBIGUOUS_XY_, MINIMUM_ALIGNMENT_SCORE_RATIO_, TRIM_LIST_, MAKE_MATCH_STRING_, QUICK_MATCH_STRINGS_, outStream_, outStreamMapped_, outStreamUnmapped_, outStreamBlack_, SLOW_ALIGN_PADDING_, SLOW_RESCUE_PADDING_, DONT_OUTPUT_UNMAPPED_READS_, DONT_OUTPUT_BLACKLISTED_READS_, MAX_SITESCORES_TO_PRINT_, PRINT_SECONDARY_ALIGNMENTS_, REQUIRE_CORRECT_STRANDS_PAIRS_, SAME_STRAND_PAIRS_, KILL_BAD_PAIRS_, RCOMP_MATE_, PERFECTMODE_, SEMIPERFECTMODE_, FORBID_SELF_MAPPING_, TIP_DELETION_SEARCH_RANGE_, AMBIGUOUS_RANDOM_, AMBIGUOUS_ALL_, KFILTER_, IDFILTER_, TRIM_LEFT_, TRIM_RIGHT_, UNTRIM_, TRIM_QUAL_, TRIM_MIN_LEN_, LOCAL_ALIGN_, RESCUE_, STRICT_MAX_INDEL_, MSA_TYPE_, bloomFilter_, TIP_DELETION_SEARCH_RANGE_);
	}

	public BBMapThread(ConcurrentReadInputStream cris_, int keylen_,
			CoveragePileup pileup_, boolean SMITH_WATERMAN_, int THRESH_, int minChrom_,
			int maxChrom_, float keyDensity_, float maxKeyDensity_, float minKeyDensity_, int maxDesiredKeys_,
			boolean REMOVE_DUPLICATE_BEST_ALIGNMENTS_, boolean SAVE_AMBIGUOUS_XY_,
			float MINIMUM_ALIGNMENT_SCORE_RATIO_, boolean TRIM_LIST_, boolean MAKE_MATCH_STRING_, boolean QUICK_MATCH_STRINGS_,
			ConcurrentReadOutputStream outStream_, ConcurrentReadOutputStream outStreamMapped_, ConcurrentReadOutputStream outStreamUnmapped_, ConcurrentReadOutputStream outStreamBlack_,
			int SLOW_ALIGN_PADDING_, int SLOW_RESCUE_PADDING_, boolean DONT_OUTPUT_UNMAPPED_READS_, boolean DONT_OUTPUT_BLACKLISTED_READS_,
			int MAX_SITESCORES_TO_PRINT_, boolean PRINT_SECONDARY_ALIGNMENTS_,
			boolean REQUIRE_CORRECT_STRANDS_PAIRS_, boolean SAME_STRAND_PAIRS_, boolean KILL_BAD_PAIRS_, boolean RCOMP_MATE_,
			boolean PERFECTMODE_, boolean SEMIPERFECTMODE_, boolean FORBID_SELF_MAPPING_, int TIP_DELETION_SEARCH_RANGE_,
			boolean AMBIGUOUS_RANDOM_, boolean AMBIGUOUS_ALL_, int KFILTER_, float IDFILTER_, boolean TRIM_LEFT_, boolean TRIM_RIGHT_, boolean UNTRIM_, float TRIM_QUAL_, int TRIM_MIN_LEN_,
			boolean LOCAL_ALIGN_, boolean RESCUE_, boolean STRICT_MAX_INDEL_, String MSA_TYPE_, BloomFilter bloomFilter_, final int hybridTipSearchCeiling_){
		this(cris_, keylen_, pileup_, SMITH_WATERMAN_, THRESH_, minChrom_, maxChrom_, keyDensity_, maxKeyDensity_, minKeyDensity_, maxDesiredKeys_, REMOVE_DUPLICATE_BEST_ALIGNMENTS_, SAVE_AMBIGUOUS_XY_, MINIMUM_ALIGNMENT_SCORE_RATIO_, TRIM_LIST_, MAKE_MATCH_STRING_, QUICK_MATCH_STRINGS_, outStream_, outStreamMapped_, outStreamUnmapped_, outStreamBlack_, SLOW_ALIGN_PADDING_, SLOW_RESCUE_PADDING_, DONT_OUTPUT_UNMAPPED_READS_, DONT_OUTPUT_BLACKLISTED_READS_, MAX_SITESCORES_TO_PRINT_, PRINT_SECONDARY_ALIGNMENTS_, REQUIRE_CORRECT_STRANDS_PAIRS_, SAME_STRAND_PAIRS_, KILL_BAD_PAIRS_, RCOMP_MATE_, PERFECTMODE_, SEMIPERFECTMODE_, FORBID_SELF_MAPPING_, TIP_DELETION_SEARCH_RANGE_, AMBIGUOUS_RANDOM_, AMBIGUOUS_ALL_, KFILTER_, IDFILTER_, TRIM_LEFT_, TRIM_RIGHT_, UNTRIM_, TRIM_QUAL_, TRIM_MIN_LEN_, LOCAL_ALIGN_, RESCUE_, STRICT_MAX_INDEL_, MSA_TYPE_, bloomFilter_, hybridTipSearchCeiling_, false);
	}

	public BBMapThread(ConcurrentReadInputStream cris_, int keylen_,
			CoveragePileup pileup_, boolean SMITH_WATERMAN_, int THRESH_, int minChrom_,
			int maxChrom_, float keyDensity_, float maxKeyDensity_, float minKeyDensity_, int maxDesiredKeys_,
			boolean REMOVE_DUPLICATE_BEST_ALIGNMENTS_, boolean SAVE_AMBIGUOUS_XY_,
			float MINIMUM_ALIGNMENT_SCORE_RATIO_, boolean TRIM_LIST_, boolean MAKE_MATCH_STRING_, boolean QUICK_MATCH_STRINGS_,
			ConcurrentReadOutputStream outStream_, ConcurrentReadOutputStream outStreamMapped_, ConcurrentReadOutputStream outStreamUnmapped_, ConcurrentReadOutputStream outStreamBlack_,
			int SLOW_ALIGN_PADDING_, int SLOW_RESCUE_PADDING_, boolean DONT_OUTPUT_UNMAPPED_READS_, boolean DONT_OUTPUT_BLACKLISTED_READS_,
			int MAX_SITESCORES_TO_PRINT_, boolean PRINT_SECONDARY_ALIGNMENTS_,
			boolean REQUIRE_CORRECT_STRANDS_PAIRS_, boolean SAME_STRAND_PAIRS_, boolean KILL_BAD_PAIRS_, boolean RCOMP_MATE_,
			boolean PERFECTMODE_, boolean SEMIPERFECTMODE_, boolean FORBID_SELF_MAPPING_, int TIP_DELETION_SEARCH_RANGE_,
			boolean AMBIGUOUS_RANDOM_, boolean AMBIGUOUS_ALL_, int KFILTER_, float IDFILTER_, boolean TRIM_LEFT_, boolean TRIM_RIGHT_, boolean UNTRIM_, float TRIM_QUAL_, int TRIM_MIN_LEN_,
			boolean LOCAL_ALIGN_, boolean RESCUE_, boolean STRICT_MAX_INDEL_, String MSA_TYPE_, BloomFilter bloomFilter_, final int hybridTipSearchCeiling_, final boolean hybridPair_){
		this(cris_, cris_.paired(), keylen_, pileup_, SMITH_WATERMAN_, THRESH_, minChrom_, maxChrom_,
				keyDensity_, maxKeyDensity_, minKeyDensity_, maxDesiredKeys_, REMOVE_DUPLICATE_BEST_ALIGNMENTS_, SAVE_AMBIGUOUS_XY_,
				MINIMUM_ALIGNMENT_SCORE_RATIO_, TRIM_LIST_, MAKE_MATCH_STRING_, QUICK_MATCH_STRINGS_,
				outStream_, outStreamMapped_, outStreamUnmapped_, outStreamBlack_,
				SLOW_ALIGN_PADDING_, SLOW_RESCUE_PADDING_, DONT_OUTPUT_UNMAPPED_READS_, DONT_OUTPUT_BLACKLISTED_READS_,
				MAX_SITESCORES_TO_PRINT_, PRINT_SECONDARY_ALIGNMENTS_, REQUIRE_CORRECT_STRANDS_PAIRS_, SAME_STRAND_PAIRS_,
				KILL_BAD_PAIRS_, RCOMP_MATE_, PERFECTMODE_, SEMIPERFECTMODE_, FORBID_SELF_MAPPING_, TIP_DELETION_SEARCH_RANGE_,
				AMBIGUOUS_RANDOM_, AMBIGUOUS_ALL_, KFILTER_, IDFILTER_, TRIM_LEFT_, TRIM_RIGHT_, UNTRIM_, TRIM_QUAL_, TRIM_MIN_LEN_,
				LOCAL_ALIGN_, RESCUE_, STRICT_MAX_INDEL_, MSA_TYPE_, bloomFilter_, hybridTipSearchCeiling_, hybridPair_, null);
	}

	/** Mapping engine constructor for callers that own input and output externally. */
	BBMapThread(boolean paired_, int keylen_,
			CoveragePileup pileup_, boolean SMITH_WATERMAN_, int THRESH_, int minChrom_,
			int maxChrom_, float keyDensity_, float maxKeyDensity_, float minKeyDensity_, int maxDesiredKeys_,
			boolean REMOVE_DUPLICATE_BEST_ALIGNMENTS_, boolean SAVE_AMBIGUOUS_XY_,
			float MINIMUM_ALIGNMENT_SCORE_RATIO_, boolean TRIM_LIST_, boolean MAKE_MATCH_STRING_, boolean QUICK_MATCH_STRINGS_,
			ConcurrentReadOutputStream outStream_, ConcurrentReadOutputStream outStreamMapped_, ConcurrentReadOutputStream outStreamUnmapped_, ConcurrentReadOutputStream outStreamBlack_,
			int SLOW_ALIGN_PADDING_, int SLOW_RESCUE_PADDING_, boolean DONT_OUTPUT_UNMAPPED_READS_, boolean DONT_OUTPUT_BLACKLISTED_READS_,
			int MAX_SITESCORES_TO_PRINT_, boolean PRINT_SECONDARY_ALIGNMENTS_,
			boolean REQUIRE_CORRECT_STRANDS_PAIRS_, boolean SAME_STRAND_PAIRS_, boolean KILL_BAD_PAIRS_, boolean RCOMP_MATE_,
			boolean PERFECTMODE_, boolean SEMIPERFECTMODE_, boolean FORBID_SELF_MAPPING_, int TIP_DELETION_SEARCH_RANGE_,
			boolean AMBIGUOUS_RANDOM_, boolean AMBIGUOUS_ALL_, int KFILTER_, float IDFILTER_, boolean TRIM_LEFT_, boolean TRIM_RIGHT_, boolean UNTRIM_, float TRIM_QUAL_, int TRIM_MIN_LEN_,
			boolean LOCAL_ALIGN_, boolean RESCUE_, boolean STRICT_MAX_INDEL_, String MSA_TYPE_, BloomFilter bloomFilter_, final int hybridTipSearchCeiling_, final boolean hybridPair_, final HybridMaxIndelConfig hybridMaxIndelConfig_){
		this(null, paired_, keylen_, pileup_, SMITH_WATERMAN_, THRESH_, minChrom_, maxChrom_,
				keyDensity_, maxKeyDensity_, minKeyDensity_, maxDesiredKeys_, REMOVE_DUPLICATE_BEST_ALIGNMENTS_, SAVE_AMBIGUOUS_XY_,
				MINIMUM_ALIGNMENT_SCORE_RATIO_, TRIM_LIST_, MAKE_MATCH_STRING_, QUICK_MATCH_STRINGS_,
				outStream_, outStreamMapped_, outStreamUnmapped_, outStreamBlack_,
				SLOW_ALIGN_PADDING_, SLOW_RESCUE_PADDING_, DONT_OUTPUT_UNMAPPED_READS_, DONT_OUTPUT_BLACKLISTED_READS_,
				MAX_SITESCORES_TO_PRINT_, PRINT_SECONDARY_ALIGNMENTS_, REQUIRE_CORRECT_STRANDS_PAIRS_, SAME_STRAND_PAIRS_,
				KILL_BAD_PAIRS_, RCOMP_MATE_, PERFECTMODE_, SEMIPERFECTMODE_, FORBID_SELF_MAPPING_, TIP_DELETION_SEARCH_RANGE_,
				AMBIGUOUS_RANDOM_, AMBIGUOUS_ALL_, KFILTER_, IDFILTER_, TRIM_LEFT_, TRIM_RIGHT_, UNTRIM_, TRIM_QUAL_, TRIM_MIN_LEN_,
				LOCAL_ALIGN_, RESCUE_, STRICT_MAX_INDEL_, MSA_TYPE_, bloomFilter_, hybridTipSearchCeiling_, hybridPair_, hybridMaxIndelConfig_);
	}

	private BBMapThread(ConcurrentReadInputStream cris_, boolean paired_, int keylen_,
			CoveragePileup pileup_, boolean SMITH_WATERMAN_, int THRESH_, int minChrom_,
			int maxChrom_, float keyDensity_, float maxKeyDensity_, float minKeyDensity_, int maxDesiredKeys_,
			boolean REMOVE_DUPLICATE_BEST_ALIGNMENTS_, boolean SAVE_AMBIGUOUS_XY_,
			float MINIMUM_ALIGNMENT_SCORE_RATIO_, boolean TRIM_LIST_, boolean MAKE_MATCH_STRING_, boolean QUICK_MATCH_STRINGS_,
			ConcurrentReadOutputStream outStream_, ConcurrentReadOutputStream outStreamMapped_, ConcurrentReadOutputStream outStreamUnmapped_, ConcurrentReadOutputStream outStreamBlack_,
			int SLOW_ALIGN_PADDING_, int SLOW_RESCUE_PADDING_, boolean DONT_OUTPUT_UNMAPPED_READS_, boolean DONT_OUTPUT_BLACKLISTED_READS_,
			int MAX_SITESCORES_TO_PRINT_, boolean PRINT_SECONDARY_ALIGNMENTS_,
			boolean REQUIRE_CORRECT_STRANDS_PAIRS_, boolean SAME_STRAND_PAIRS_, boolean KILL_BAD_PAIRS_, boolean RCOMP_MATE_,
			boolean PERFECTMODE_, boolean SEMIPERFECTMODE_, boolean FORBID_SELF_MAPPING_, int TIP_DELETION_SEARCH_RANGE_,
			boolean AMBIGUOUS_RANDOM_, boolean AMBIGUOUS_ALL_, int KFILTER_, float IDFILTER_, boolean TRIM_LEFT_, boolean TRIM_RIGHT_, boolean UNTRIM_, float TRIM_QUAL_, int TRIM_MIN_LEN_,
			boolean LOCAL_ALIGN_, boolean RESCUE_, boolean STRICT_MAX_INDEL_, String MSA_TYPE_, BloomFilter bloomFilter_, final int hybridTipSearchCeiling_, final boolean hybridPair_, final HybridMaxIndelConfig hybridMaxIndelConfig_){
		
		super(cris_, paired_,
				outStream_, outStreamMapped_, outStreamUnmapped_, outStreamBlack_,
				pileup_, SMITH_WATERMAN_, LOCAL_ALIGN_, REMOVE_DUPLICATE_BEST_ALIGNMENTS_,
				AMBIGUOUS_RANDOM_, AMBIGUOUS_ALL_, TRIM_LEFT_, TRIM_RIGHT_, UNTRIM_, TRIM_QUAL_, TRIM_MIN_LEN_, THRESH_,
				minChrom_, maxChrom_, KFILTER_, IDFILTER_, KILL_BAD_PAIRS_, SAVE_AMBIGUOUS_XY_,
				REQUIRE_CORRECT_STRANDS_PAIRS_,
				SAME_STRAND_PAIRS_, RESCUE_, STRICT_MAX_INDEL_, SLOW_ALIGN_PADDING_, SLOW_RESCUE_PADDING_,
				MSA_TYPE_, keylen_, PERFECTMODE_, SEMIPERFECTMODE_, FORBID_SELF_MAPPING_, RCOMP_MATE_,
				MAKE_MATCH_STRING_, DONT_OUTPUT_UNMAPPED_READS_, DONT_OUTPUT_BLACKLISTED_READS_, PRINT_SECONDARY_ALIGNMENTS_,
				QUICK_MATCH_STRINGS_, MAX_SITESCORES_TO_PRINT_, MINIMUM_ALIGNMENT_SCORE_RATIO_,
				keyDensity_, maxKeyDensity_, minKeyDensity_, maxDesiredKeys_,
				BBIndex.MIN_APPROX_HITS_TO_KEEP, BBIndex.USE_EXTENDED_SCORE,
				BBIndex.BASE_HIT_SCORE, BBIndex.USE_AFFINE_SCORE, BBIndex.MAX_INDEL, TRIM_LIST_, TIP_DELETION_SEARCH_RANGE_, bloomFilter_);
		assert(hybridTipSearchCeiling_>=TIP_DELETION_SEARCH_RANGE_);
		this.hybridTipSearchCeiling=hybridTipSearchCeiling_;
		this.hybridPair=hybridPair_;
		this.hybridMaxIndelConfig=hybridMaxIndelConfig_;
		if(hybridPair && hybridMaxIndelConfig!=null){
			throw new IllegalArgumentException("Legacy hybridpair and general hybridmaxindel are mutually exclusive");
		}
		quantumOnly=Boolean.getBoolean("bbmap3.quantumOnly");
		quantumTieredShadow=Boolean.getBoolean("bbmap3.quantumTieredShadow");
		quantumTieredMutate=Boolean.getBoolean("bbmap3.quantumTieredMutate");
		quantumHybrid=Boolean.getBoolean("bbmap3.quantumHybrid");
		final boolean tieredEnabled=quantumTieredShadow || quantumTieredMutate || quantumOnly || quantumHybrid;
		quantumTieredRanker=(tieredEnabled ? new QuantumRanker() : null);
		quantumTieredStats=(quantumTieredShadow || quantumTieredMutate ?
				new TieredQuantumStats() : null);
		quantumOnlyStats=(quantumOnly ? new QuantumOnlyStats() : null);
		quantumHybridStats=(quantumHybrid ? new QuantumHybridStats() : null);
		assert(!quantumOnly || msa==null) : "Quantum-only workers must not allocate an MSA";
		assert(!(quantumOnly && PSEUDO_ONLY)) : "Quantum-only and pseudoalignment modes are exclusive";
		assert(!PSEUDO_ONLY || msa==null) : "Pseudoalignment workers must not allocate an MSA";
		final boolean pairedHybrid=hybridPair || (paired_ && hybridMaxIndelConfig!=null);
		secondAttempt=pairedHybrid ? new PairAttemptAccounting() : null;
		entryState=pairedHybrid ? new PairSearchState() : null;
		firstState=pairedHybrid ? new PairSearchState() : null;
		singleSecondAttempt=hybridMaxIndelConfig==null ? null : new SingleAttemptAccounting();
		singleEntryState=hybridMaxIndelConfig==null ? null : new ReadSearchState();
		singleFirstState=hybridMaxIndelConfig==null ? null : new ReadSearchState();
		hybridMaxIndelStats=hybridMaxIndelConfig==null ? null : new HybridMaxIndelStats();
		if(hybridMaxIndelConfig!=null){hybridMatchCache=new HybridMatchCache(hybridMaxIndelStats);}
		else if(hybridPair && Boolean.getBoolean("bbmap3.hybridMatchReuse")){hybridMatchCache=new HybridMatchCache(null);}
		
		assert(SLOW_ALIGN_PADDING>=0);
		assert(!(RCOMP_MATE/* || FORBID_SELF_MAPPING*/)) : "RCOMP_MATE: TODO";
		
		if(SLOW_ALIGN || MAKE_MATCH_STRING){
//			msa=MSA.makeMSA(ALIGN_ROWS, ALIGN_COLUMNS, MSA_TYPE);
//			POINTS_MATCH=msa.POINTS_MATCH();
//			POINTS_MATCH2=msa.POINTS_MATCH2();
			CLEARZONE1=(int)(CLEARZONE_RATIO1*POINTS_MATCH2);
			CLEARZONE1b=(int)(CLEARZONE_RATIO1b*POINTS_MATCH2);
			CLEARZONE1c=(int)(CLEARZONE_RATIO1c*POINTS_MATCH2);
			CLEARZONEP=(int)(CLEARZONE_RATIOP*POINTS_MATCH2);
			CLEARZONE3=PENALIZE_AMBIG ? (int)(CLEARZONE_RATIO3*POINTS_MATCH2) : 0;
//			CLEARZONE1e=(int)(2*POINTS_MATCH2-POINTS_MATCH-msa.POINTS_SUB())+1;
		}else{
//			POINTS_MATCH=70;
//			POINTS_MATCH2=100;
//			msa=null;
			CLEARZONE1=0;
			CLEARZONE1b=0;
			CLEARZONE1c=0;
			CLEARZONEP=0;
			CLEARZONE3=0;
//			CLEARZONE1e=0;
		}
		
		CLEARZONE1b_CUTOFF_FLAT=CLEARZONE1b_CUTOFF_FLAT_RATIO*POINTS_MATCH2;
		CLEARZONE1c_CUTOFF_FLAT=CLEARZONE1c_CUTOFF_FLAT_RATIO*POINTS_MATCH2;
		INV_CLEARZONE3=(CLEARZONE3==0 ? 0 : 1f/CLEARZONE3);
		
		index=new BBIndex(KEYLEN, minChrom, maxChrom, KFILTER, msa);
	}
	
	
	@Override
	public int trimList(ArrayList<SiteScore> list, boolean retainPaired, int maxScore, boolean specialCasePerfect, int minSitesToRetain, int maxSitesToRetain){
		if(list==null || list.size()==0){return -99999;}
		if(list.size()==1){return list.get(0).score;}
		
		final int highestScore;
		if(USE_AFFINE_SCORE){
			
			highestScore=Tools.trimSiteList(list, .6f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
			if(highestScore==maxScore && specialCasePerfect){
				Tools.trimSiteList(list, .94f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				if(list.size()>8){Tools.trimSiteList(list, .99f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
				return highestScore;
			}
			
			final int mstr2=(minSitesToRetain<=1 ? 1 : minSitesToRetain+1);
			
//			if(list.size()>6){Tools.trimSiteList(list, .65f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
////			//			System.out.print(", "+list.size());
//			if(list.size()>10){Tools.trimSiteList(list, .7f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
////			//			System.out.print(", "+list.size());
//			if(list.size()>14){Tools.trimSiteList(list, .75f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
////			//			System.out.print(", "+list.size());
//			if(list.size()>18){Tools.trimSiteList(list, .8f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
////			////			System.out.print(", "+list.size());
//			if(list.size()>22){Tools.trimSiteList(list, .85f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
////			////			System.out.print(", "+list.size());
//			if(list.size()>26){Tools.trimSiteList(list, .9f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
////			//			System.out.print(", "+list.size());
//			if(list.size()>34){Tools.trimSiteList(list, .95f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
//			if(list.size()>42){Tools.trimSiteList(list, .98f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
//			if(list.size()>50){Tools.trimSiteList(list, .99f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
////			if(list.size()>64){Tools.trimSiteList(list, .995f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
			
			if(list.size()>4){Tools.trimSiteList(list, .65f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
			if(list.size()>8){Tools.trimSiteList(list, .7f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
			if(list.size()>12){Tools.trimSiteList(list, .75f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
			if(list.size()>16){Tools.trimSiteList(list, .8f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			////			System.out.print(", "+list.size());
			if(list.size()>20){Tools.trimSiteList(list, .85f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			////			System.out.print(", "+list.size());
			if(list.size()>24){Tools.trimSiteList(list, .9f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
			if(list.size()>32){Tools.trimSiteList(list, .95f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>40){Tools.trimSiteList(list, .97f, retainPaired, true, mstr2, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>48){Tools.trimSiteList(list, .99f, retainPaired, true, mstr2, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
//			if(list.size()>64){Tools.trimSiteList(list, .995f, retainPaired, true, mstr2, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			

		}else if(USE_EXTENDED_SCORE){
			highestScore=Tools.trimSiteList(list, .75f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
			
			if(list.size()>8){Tools.trimSiteList(list, .8f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>16){Tools.trimSiteList(list, .85f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>24){Tools.trimSiteList(list, .90f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>36){Tools.trimSiteList(list, .92f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			if(list.size()>40){Tools.trimSiteList(list, .94f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			if(list.size()>48){Tools.trimSiteList(list, .96f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			if(list.size()>56){Tools.trimSiteList(list, .97f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>64){Tools.trimSiteList(list, .98f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>80){Tools.trimSiteList(list, .99f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			

		}else{
			//			System.out.print("\n\nSize:\t"+list.size());


			highestScore=Tools.trimSiteList(list, .6f, retainPaired, true, minSitesToRetain, maxSitesToRetain);

			if(list.size()>12){Tools.trimSiteList(list, .65f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>16){Tools.trimSiteList(list, .7f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>24){Tools.trimSiteList(list, .74f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>28){Tools.trimSiteList(list, .8f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>32){Tools.trimSiteList(list, .85f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>48){Tools.trimSiteList(list, .90f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			//			if(list.size()>40){Tools.trimSiteList(list, .95f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			//			if(list.size()>48){Tools.trimSiteList(list, .96f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			//			if(list.size()>56){Tools.trimSiteList(list, .97f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
		}
		
		return highestScore;
	}
	
	
	@Override
	public void scoreSlow(final ArrayList<SiteScore> list, final byte[] basesP, final byte[] basesM,
			final int maxSwScore, final int maxImperfectSwScore){
		if(quantumHybrid && list.size()>1){
			int next=0;
			for(int i=0; i<list.size(); i++){
				final SiteScore ss=list.get(i);
				if(ss.gaps==null){
					if(i>next){list.remove(i);list.add(next, ss);}
					next++;
				}
			}
		}
		int minMsaLimit;
		if(PAIRED){
			minMsaLimit=-CLEARZONE1e+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE*maxSwScore);
		}else{
			minMsaLimit=-CLEARZONE1e+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore);
		}
		assert(Read.CHECKSITES(list, basesP, basesM, -1));
		
		int minMatch=Tools.max(-300, minMsaLimit-CLEARZONE3); //Score must exceed this to generate quick match string
		if(verbose){
			System.err.println("Slow-scoring.  maxSwScore="+maxSwScore+", maxImperfectSwScore="+maxImperfectSwScore+", minMsaLimit="+minMsaLimit+", minMatch="+minMatch);
		}
		int bestOrdinaryScore=Integer.MIN_VALUE;
		SiteScore bestOrdinarySite=null;
		for(int i=0; i<list.size(); i++){
			final SiteScore ss=list.get(i);
			assert(ss.lengthsAgree());
			final byte[] bases=(ss.strand==Shared.PLUS ? basesP : basesM);
			final boolean compressed=(ss.gaps!=null);
			final int scoreAwareRelax=(maxSwScore-maxImperfectSwScore)*2/3;
			final boolean scoreAwareCandidate=quantumHybrid && compressed && bestOrdinarySite!=null &&
					bestOrdinaryScore>=maxImperfectSwScore-scoreAwareRelax &&
					ss.hits<bestOrdinarySite.hits && ss.quickScore<bestOrdinarySite.quickScore;
			final boolean applyScoreAware=scoreAwareCandidate && quantumHybridStats.adaptiveOpportunity();
			if(quantumHybrid && compressed && bestOrdinaryScore>=maxImperfectSwScore){
				quantumHybridStats.compressedSkipped();
				ss.match=null;ss.setSlowScore(0);ss.setScore(0);
				ss.perfect=ss.semiperfect=false;
				continue;
			}else if(applyScoreAware){
				quantumHybridStats.scoreAwareSkipped();
				ss.match=null;ss.setSlowScore(0);ss.setScore(0);
				ss.perfect=ss.semiperfect=false;
				continue;
			}else if(quantumHybrid && compressed){
				quantumHybridStats.compressedEvaluated();
			}
			
			if(SEMIPERFECTMODE){
				assert(ss.stop-ss.start==bases.length-1);
				assert(ss.semiperfect);
			}
			
			if(verbose){System.err.println("\nSlow-scoring "+ss);}
			if(ss.stop-ss.start!=bases.length-1){
				assert(ss.stop-ss.start>bases.length-1) : bases.length+", "+ss.toText();
				assert(!ss.semiperfect) : "\n"+bases.length+", "+ss.toText()+", "+ss.perfect+", "+ss.semiperfect+", "+maxSwScore+"\n"+new String(basesP)+"\n";
				ss.setSlowScore(0);
				ss.semiperfect=false;
				ss.perfect=false;
			}
			
			final int swscoreNoIndel=ss.slowScore;
			int[] swscoreArray=null;
			
			boolean clipped=true, setLimits=false;
			if(swscoreNoIndel<maxImperfectSwScore && !ss.semiperfect){
				if(verbose && ss.stop-ss.start>4000){
					System.err.println(ss.toText());
					System.err.println(list.size());
					System.err.println();
				}
				
				int expectedLen=GapTools.calcGrefLen(ss);
				if(verbose){System.err.println("expectedLen="+expectedLen);}
				if(expectedLen>=EXPECTED_LEN_LIMIT){
					 //TODO: Alternately, I could kill the site.
					ss.setStop(ss.start+Tools.min(basesP.length+40, EXPECTED_LEN_LIMIT));
					if(verbose){System.err.println("expectedLen="+expectedLen+"; ss="+ss);}
				}
				
				int pad=SLOW_ALIGN_PADDING;
				final int minscore=Tools.max(swscoreNoIndel, minMsaLimit);
				final int minscore2=Tools.max(swscoreNoIndel-MSA.MIN_SCORE_ADJUST, minMsaLimit);
				if(verbose){System.err.println("Sent to msa with start="+ss.start+", stop="+ss.stop+", pad="+pad+", limit="+minscore+", gaps="+GapTools.toString(ss.gaps));}
				boolean bypass=false;
				QuantumRanker.Result tieredResult=null;
				int tieredScore=Integer.MIN_VALUE;
				final int projectedLengthDelta=ss.stop-ss.start+1-bases.length;
				if(quantumOnly){
					if(ss.gaps!=null){
						quantumOnlyStats.gapArray();
						ss.setSlowScore(0);
					}else{
						final byte[] ref=Data.getChromosome(ss.chrom).array;
						final int refStart=Math.max(0, ss.start-pad);
						final int refStop=Math.min(ref.length-1, ss.stop+pad);
						tieredResult=quantumTieredRanker.align(bases, ref, refStart, refStop,
								ss.start, Shared.SIMD, true, Integer.MAX_VALUE);
						if(tieredResult.supported && tieredResult.match!=null){
							tieredScore=quantumScaledScore(tieredResult.score);
						}
						final boolean accepted=tieredResult.supported && !tieredResult.uncertain &&
								tieredResult.match!=null && tieredScore>=minscore;
						quantumOnlyStats.observe(tieredResult, tieredScore, minscore);
						if(accepted){
							bypass=true;
							ss.match=tieredResult.match;
							ss.setSlowScore(tieredScore);
							ss.setLimits(tieredResult.rStart, tieredResult.rStop);
							setLimits=true;
						}else{
							ss.match=null;
							ss.setSlowScore(0);
						}
					}
				}else if(quantumHybrid && ss.gaps==null && Math.abs(projectedLengthDelta)<=1){
					quantumHybridStats.ordinaryAttempt();
					final byte[] ref=Data.getChromosome(ss.chrom).array;
					final int refStart=Math.max(0, ss.start-pad);
					final int refStop=Math.min(ref.length-1, ss.stop+pad);
					final int editBudget=quantumEditBudget(maxSwScore, minscore, bases.length);
					tieredResult=quantumTieredRanker.align(bases, ref, refStart, refStop,
							ss.start, Shared.SIMD, true, editBudget);
					final boolean selectedTraceUncertain=tieredResult.match!=null && matchContainsN(tieredResult.match);
					if(tieredResult.supported && tieredResult.match!=null){tieredScore=msa.score(tieredResult.match);}
					final int indelBases=tieredResult.insertions+tieredResult.deletions;
					if(tieredResult.supported && tieredResult.match!=null && !selectedTraceUncertain &&
							indelBases<=1 && tieredScore>=minscore){
						bypass=true;ss.match=tieredResult.match;ss.setSlowScore(tieredScore);
						ss.setLimits(tieredResult.rStart, tieredResult.rStop);setLimits=true;
						quantumHybridStats.accepted();
					}else if(!tieredResult.supported){quantumHybridStats.fallbackUnsupported();}
					else if(tieredResult.match==null){quantumHybridStats.fallbackMissing();}
					else if(selectedTraceUncertain){quantumHybridStats.fallbackUncertain();}
					else if(indelBases>1){quantumHybridStats.longTraceFallback();}
					else{quantumHybridStats.fallbackThreshold();}
				}else if(quantumHybrid && ss.gaps==null){
					quantumHybridStats.geometryFallback();
				}else if((quantumTieredShadow || quantumTieredMutate) && ss.gaps==null &&
						Math.abs(projectedLengthDelta)==1){
					final byte[] ref=Data.getChromosome(ss.chrom).array;
					final int refStart=Math.max(0, ss.start-pad);
					final int refStop=Math.min(ref.length-1, ss.stop+pad);
					final int editBudget=quantumEditBudget(maxSwScore, minscore, bases.length);
					tieredResult=quantumTieredRanker.align(bases, ref, refStart, refStop,
							ss.start, Shared.SIMD, true, editBudget);
					if(tieredResult.supported && tieredResult.match!=null){
						tieredScore=msa.score(tieredResult.match);
					}
				}
				if(!quantumOnly && quantumTieredMutate && tieredResult!=null && tieredResult.supported &&
						!tieredResult.uncertain && tieredResult.insertions+tieredResult.deletions<=1 &&
						tieredScore>=minscore){
					bypass=true;
					ss.match=tieredResult.match;
					ss.setSlowScore(tieredScore);
					ss.setLimits(tieredResult.rStart, tieredResult.rStop);
					setLimits=true;
				}
				if(!quantumOnly && !bypass){swscoreArray=msa.fillAndScoreLimited(bases, ss, pad, minscore);}
				if(verbose){System.err.println("Received "+Arrays.toString(swscoreArray));}
				
				if(!quantumOnly && swscoreArray!=null && swscoreArray.length>6 && (swscoreArray[3]+swscoreArray[4]+expectedLen<EXPECTED_LEN_LIMIT)){
					int[] oldArray=swscoreArray.clone();
					assert(swscoreArray.length==8);
					int extraPadLeft=swscoreArray[6];
					int extraPadRight=swscoreArray[7];
					
					if(verbose){
						System.err.println("msa returned "+Arrays.toString(swscoreArray)+", re-running.");
						System.err.println("Added extra padding: "+ss.toText()+", "+Arrays.toString(oldArray));
					}
					
					ss.setLimits(ss.start-extraPadLeft, ss.stop+extraPadRight);
					pad=SLOW_ALIGN_PADDING+EXTRA_PADDING;
					if(verbose){System.err.println("Sent to msa with start="+ss.start+", stop="+ss.stop+", pad="+pad+", limit="+minscore+", gaps="+GapTools.toString(ss.gaps));}
					swscoreArray=msa.fillAndScoreLimited(bases, ss, pad, minscore);
					
					if(verbose){System.err.println("Result of extra padding: "+ss.toText()+", "+Arrays.toString(swscoreArray));}
					if(swscoreArray==null || swscoreArray[0]<oldArray[0]){
						if(verbose){
							System.err.println("Result was inferior.");
						}
						swscoreArray=oldArray;
					}
				}
				if((quantumTieredShadow || quantumTieredMutate) && tieredResult!=null){
					final int legacyScore=(swscoreArray==null ? Integer.MIN_VALUE : swscoreArray[0]);
					final int legacyStart=(swscoreArray==null ? Integer.MIN_VALUE : swscoreArray[1]);
					final int legacyStop=(swscoreArray==null ? Integer.MIN_VALUE : swscoreArray[2]);
					quantumTieredStats.observe(tieredResult, tieredScore, legacyScore,
							legacyStart, legacyStop, minscore);
				}
				assert(ss.lengthsAgree());
				if(verbose){
					System.err.println(QUICK_MATCH_STRINGS+", "+(swscoreArray==null ? "null" : (swscoreArray.length+", "+swscoreArray[0]+" >=? "+minscore)));
					System.err.println("start="+ss.start+", stop="+ss.stop+", len="+ss.mappedLength());
				}
				if(bypass){
					assert(ss.match!=null && matchQueryLength(ss.match)==bases.length) :
							"Accepted one-base Quantum result requires a reusable match string";
				}else if(QUICK_MATCH_STRINGS && swscoreArray!=null && swscoreArray.length==6 && swscoreArray[0]>=minscore2 && (PRINT_SECONDARY_ALIGNMENTS || (USE_SS_MATCH_FOR_PRIMARY && swscoreArray[0]>minMatch))){
					if(verbose){System.err.println("Generating match string.");}
					assert(swscoreArray.length==6) : swscoreArray.length;
					assert(swscoreArray[0]>=minscore2) : "\n"+Arrays.toString(swscoreArray)+"\n"+minscore+"\n"+minMatch;
					ss.match=msa.traceback(bases, Data.getChromosome(ss.chrom).array, ss.start-pad, ss.stop+pad, swscoreArray[3], swscoreArray[4], swscoreArray[5], ss.gaps!=null);
					if(ss.match!=null){
						assert(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore); //123
						ss.setLimits(swscoreArray[1], swscoreArray[2]);
						setLimits=true;
						// Early traceback skips the later !setLimits score transfer. Save DP score before clipping may rescore.
						ss.setSlowScore(swscoreArray[0]);
						assert(ss.lengthsAgree());
						clipped=ss.fixXY(bases, true, msa);
						assert(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore); //123
						clipped=ss.clipTipIndels(bases, basesM, 4, 10, msa) || clipped;
						assert(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore); //123
						assert(ss.lengthsAgree());
					}
				}else{
					ss.match=null;
				}
			}
			if(quantumOnly){
				assert(msa==null) : "Quantum-only scoring must not construct an MSA";
			}else if(swscoreArray!=null && !setLimits){
				if(verbose){System.err.println("msa returned "+Arrays.toString(swscoreArray));}
				ss.setSlowScore(swscoreArray[0]);
				ss.setLimits(swscoreArray[1], swscoreArray[2]);
				assert(ss.lengthsAgree());
			}else{
				assert(swscoreNoIndel<=maxSwScore) : swscoreNoIndel+", "+maxImperfectSwScore+", "+maxSwScore+", "+new String(basesP);
				assert(clipped || swscoreNoIndel==-1 || msa.scoreNoIndels(bases, ss.chrom, ss.start)==swscoreNoIndel) :
					setLimits+", "+clipped+", "+(swscoreArray==null)+", "+
					swscoreNoIndel+" != "+msa.scoreNoIndels(bases, ss.chrom, ss.start)+"\n"+
					ss.toText()+"\n"+(ss.stop-ss.start)+", "+bases.length; //Slow
			}
			assert(ss.lengthsAgree());
			ss.setScore(ss.slowScore);
			if(quantumHybrid && !compressed && ss.slowScore>bestOrdinaryScore){
				bestOrdinaryScore=ss.slowScore;bestOrdinarySite=ss;
			}
			minMatch=Tools.max(minMatch, ss.slowScore);
			minMsaLimit=Tools.max(minMsaLimit, ss.slowScore-CLEARZONE3);
			assert(ss.slowScore<=maxSwScore);
			assert(!(ss.perfect && ss.slowScore<maxSwScore));
			ss.perfect=(ss.slowScore==maxSwScore);
			if(ss.perfect){ss.semiperfect=true;}
			else if(!ss.semiperfect){ss.setPerfect(bases);}
			
			if(verbose){System.err.println(" -> "+ss);}
		}
	}

	TieredQuantumStats quantumTieredStats(){return quantumTieredStats;}
	QuantumOnlyStats quantumOnlyStats(){return quantumOnlyStats;}
	QuantumHybridStats quantumHybridStats(){return quantumHybridStats;}
	private static int matchQueryLength(final byte[] match){
		int length=0;
		for(byte op : match){if(op!='D'){length++;}}
		return length;
	}
	private static boolean matchContainsN(final byte[] match){
		if(match==null){return false;}
		for(byte op : match){if(op=='N'){return true;}}
		return false;
	}
	private int quantumEditBudget(final int maxSwScore, final int minimumScore,
			final int readLength){
		if(quantumOnly){
			final int deficit=Math.max(0, maxSwScore-minimumScore);
			return Math.min(readLength, deficit/QUANTUM_EDIT_COST);
		}
		final int approximateEditCost=(POINTS_MATCH-msa.POINTS_SUB()+1)/2;
		assert(approximateEditCost>0) : "Half a substitution penalty must be positive";
		final int deficit=Math.max(0, maxSwScore-minimumScore);
		return Math.min(readLength, deficit/approximateEditCost);
	}

	private static int quantumMaxScore(final int readLength){
		return QUANTUM_SCORE_SCALE*readLength+QUANTUM_SCORE_OFFSET;
	}

	private static int quantumScaledScore(final int simpleScore){
		return QUANTUM_SCORE_SCALE*simpleScore+QUANTUM_SCORE_OFFSET;
	}

	private int scoreNoIndelsQuantum(final Read read, final byte[] basesP,
			final byte[] basesM, final int maxSwScore, final int maxImperfectSwScore){
		if(read.numSites()==0){return 0;}
		int nearPerfect=0;
		boolean forceQuantum=false;
		for(SiteScore ss : read.sites){
			final byte[] bases=(ss.strand==Shared.PLUS ? basesP : basesM);
			final byte[] ref=Data.getChromosome(ss.chrom).array;
			final int oldScore=ss.score;
			int start=ss.start;
			int score=quantumGaplessScore(bases, ref, start, null);
			if(score<oldScore && oldScore>=maxImperfectSwScore &&
					ss.stop-ss.start+1!=bases.length){
				final int alternateStart=ss.stop-bases.length+1;
				final int alternateScore=quantumGaplessScore(bases, ref, alternateStart, null);
				if(alternateScore>score){score=alternateScore;start=alternateStart;}
			}
			ss.setStart(start);
			ss.setSlowScore(score);
			ss.setScore(score);
			if(score>=maxImperfectSwScore){
				nearPerfect++;
				ss.setStop(start+bases.length-1);
				ss.gaps=null;
				final byte[] match=new byte[bases.length];
				final int check=quantumGaplessScore(bases, ref, start, match);
				assert(check==score) : "Quantum gapless score changed while tracing: "+check+" != "+score;
				ss.match=match;
				if(score==maxSwScore){ss.perfect=ss.semiperfect=true;}
				else{ss.setPerfect(bases);}
			}else if(oldScore>=maxImperfectSwScore || PRINT_SECONDARY_ALIGNMENTS){
				forceQuantum=true;
			}
		}
		return forceQuantum ? -nearPerfect : nearPerfect;
	}

	private static int quantumGaplessScore(final byte[] query, final byte[] ref,
			final int refStart, final byte[] match){
		if(refStart<0 || refStart+query.length>ref.length){return Integer.MIN_VALUE/4;}
		int simpleScore=0;
		for(int q=0, r=refStart; q<query.length; q++, r++){
			final byte qb=query[q], rb=ref[r];
			if(qb=='N' || rb=='N'){
				if(match!=null){match[q]='N';}
			}else if(qb==rb){
				simpleScore++;
				if(match!=null){match[q]='m';}
			}else{
				simpleScore--;
				if(match!=null){match[q]='S';}
			}
		}
		return quantumScaledScore(simpleScore);
	}

	private void removeQuantumRejectedSites(final Read read){
		if(!quantumOnly || read.sites==null){return;}
		int removed=0;
		for(int i=read.sites.size()-1; i>=0; i--){
			final SiteScore ss=read.sites.get(i);
			if(ss.match==null || ss.slowScore<=0){
				read.sites.set(i, null);
				removed++;
			}
		}
		if(removed>0){Tools.condenseStrict(read.sites);}
		if(read.sites.isEmpty()){read.clearMapping();}
	}
	
	
	@Override
	public void processRead(final Read r, final byte[] basesM){
		if(idmodulo>1 && r.numericID%idmodulo!=1){return;}
		if(hybridMaxIndelConfig!=null){
			processReadHybridMaxIndel(r,basesM);
			return;
		}
		singleAttempt.reset();
		searchAndScoreRead(r,basesM,singleAttempt);
		singleAttempt.commit(this);
		if(!singleAttempt.rejected){finishRead(r,singleAttempt);}
	}

	/** First functional route: retry only a nonrejected single read left unmapped by the low attempt. */
	private void processReadHybridMaxIndel(final Read r, final byte[] basesM){
		final HybridMaxIndelConfig config=hybridMaxIndelConfig;
		if(index.maxIndel()!=config.lowPrimary || index.maxIndel2()!=config.lowSum){
			throw new IllegalStateException("Selective max-indel worker must enter at low bounds "+config+
				"; got "+index.maxIndel()+"/"+index.maxIndel2());
		}
		singleEntryState.capture(r);
		singleAttempt.reset();
		SingleAttemptAccounting chosen=singleAttempt;
		try{
			searchAndScoreRead(r,basesM,singleAttempt);
			singleFirstState.capture(r);
			if(!singleAttempt.rejected && !r.mapped()){
				hybridMaxIndelStats.retryAttempted();
				hybridMaxIndelStats.singleNoAccepted();
				singleEntryState.restore();
				r.sites=null;
				index.setRuntimeIndelLimits(config.highPrimary,config.highSum);
				singleSecondAttempt.reset();
				searchAndScoreRead(r,basesM,singleSecondAttempt);
				if(!singleSecondAttempt.rejected && r.mapped() && retryMapqAccepted(r,config.retryMinMapq)){
					singleAttempt.discard();
					chosen=singleSecondAttempt;
					hybridMaxIndelStats.wideSelected();
				}else{
					if(!singleSecondAttempt.rejected && r.mapped()){hybridMaxIndelStats.mapqRejected();}
					singleSecondAttempt.discard();
					singleFirstState.restore();
					hybridMaxIndelStats.lowRestored();
				}
			}
			chosen.commit(this);
			if(!chosen.rejected){finishRead(r,chosen);}
		}finally{
			index.setRuntimeIndelLimits(config.lowPrimary,config.lowSum);
			singleEntryState.clear();singleFirstState.clear();
		}
	}

	private static boolean retryMapqAccepted(final Read r,final int minimum){
		return minimum<=0 || !r.mapped() || SamLine.toMapq(r,null)>=minimum;
	}

	private static boolean retryPairMapqAccepted(final Read r,final int minimum){
		if(minimum<=0){return true;}
		final Read r2=r.mate;
		return (r.mapped() || r2.mapped()) && retryMapqAccepted(r,minimum) && retryMapqAccepted(r2,minimum);
	}

	/** One complete single-read search attempt; final statistics commit separately. */
	private void searchAndScoreRead(final Read r, final byte[] basesM, final SingleAttemptAccounting attempt){
		assert(!attempt.committed) : "Single-read attempt accounting must be reset before search";
		final byte[] basesP=r.bases;
		
//		System.err.print(" rd#"+r.numericID+" ");
//		if(r.numericID==25967){
//			verbose=true;
//			msa.verbose=true;
//			GapTools.verbose=true;
//			index.verbose=true;
//			tcr.verbose=true;
//		}
		
		if(verbose){System.err.println("\nProcessing "+r);}
		
		final int maxPossibleQuickScore=attempt.quick=quickMap(r, basesM);
		if(verbose){System.err.println("\nQuick Map: \t"+r.sites);}
		
		if(maxPossibleQuickScore<0){
			r.sites=null;
			attempt.rejected=true;
			attempt.rejectedBases=basesP.length;
			r.setDiscarded(true);
			return;
		}
		attempt.initialSiteSum=r.numSites();
		if(verbose){System.err.println("\ninitialSiteSum1: "+(initialSiteSum1+attempt.initialSiteSum));}
		
		int maxSwScore=0;
		int maxImperfectSwScore=0;

		if(PSEUDO_ONLY){
			maxSwScore=Tools.max(1, maxPossibleQuickScore);
			maxImperfectSwScore=maxSwScore;
		}else if(SLOW_ALIGN || USE_AFFINE_SCORE){
			maxSwScore=(quantumOnly ? quantumMaxScore(r.length()) : msa.maxQuality(r.length()));
			maxImperfectSwScore=(quantumOnly ? maxSwScore-QUANTUM_EDIT_COST :
				msa.maxImperfectScore(r.length()));
		}
		attempt.max=maxSwScore;
		attempt.imperfect=maxImperfectSwScore;
		
		if(TRIM_LIST && r.numSites()>1){
			if(MIN_TRIM_SITES_TO_RETAIN_SINGLE>1){Shared.sort(r.sites);}
			int highestQuickScore=trimList(r.sites, false, maxSwScore, true, MIN_TRIM_SITES_TO_RETAIN_SINGLE, MAX_TRIM_SITES_TO_RETAIN);
		}
		attempt.postTrimSiteSum=r.numSites();
		if(verbose){System.err.println("\nAfter trim: \t"+r.sites);}
		
		assert(Read.CHECKSITES(r, basesM));
		
		
		if(SLOW_ALIGN && !PSEUDO_ONLY && r.numSites()>0){
			
			int numNearPerfectScores=(quantumOnly ?
					scoreNoIndelsQuantum(r, basesP, basesM, maxSwScore, maxImperfectSwScore) :
					scoreNoIndels(r, basesP, basesM, maxSwScore, maxImperfectSwScore));

			Shared.sort(r.sites); //Puts higher scores first to better trigger the early exit based on perfect scores
			assert(Read.CHECKSITES(r, basesM));
			
//			int numPerfectScores=0;
//			if(numNearPerfectScores>0){
//				for(SiteScore ss : r.list){
//					if(ss.perfect){numPerfectScores++;}
//					else{break;}
//				}
//			}
			
			if(verbose){
				System.err.println("\nAfter scoreNoIndels: \t"+r.sites);
			}
			
			if(numNearPerfectScores<1){
				if(!quantumOnly && FIND_TIP_DELETIONS){
					findTipDeletions(r, basesP, basesM, maxSwScore, maxImperfectSwScore);
				}
			}
			
			if(verbose){
				System.err.println("\nAfter findTipDeletions: \t"+r.sites);
			}
			
			//TODO: This causes problems with perfect matches that are mapped to areas longer than the read length
			//***Above note should be resolved now, but needs to be verified.
			
			if(numNearPerfectScores<1){
				scoreSlow(r.sites, basesP, basesM, maxSwScore, maxImperfectSwScore);
			}
			
			if(STRICT_MAX_INDEL){
				int removed=removeLongIndels(r.sites, index.maxIndel());
				if(r.numSites()==0){r.clearMapping();}
			}
			removeQuantumRejectedSites(r);
			
			if(verbose){System.err.println("\nAfter scoreSlow: \t"+r.sites);}
			assert(Read.CHECKSITES(r, basesM, false));
		}


		if(r.numSites()>0){
			attempt.mapped=1;
			try {
				Tools.mergeDuplicateSites(r.sites, true, true);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				throw new RuntimeException("\n\n"+r.toText(false)+"\n\n");
			}
			Shared.sort(r.sites);
		}
		
		if(r.numSites()>1){
			SiteScore ss1=r.topSite();
			SiteScore ss2=r.sites.get(1);
			//Ensure no duplicates
			assert(ss1.chrom!=ss2.chrom || ss1.strand!=ss2.strand || ss1.start!=ss2.start || ss1.stop!=ss2.stop) : r.toText(false);
		}
		assert(Read.CHECKSITES(r, basesM));
		
		if(r.numSites()>=1){
			assert(r.topSite().score==r.topSite().slowScore) : r.topSite();
		}
		
		if(!PSEUDO_ONLY && (SLOW_ALIGN || USE_AFFINE_SCORE)){r.setPerfectFlag(maxSwScore);}
		
		if(r.numSites()>1){
			
			final int clearzone;
			final int score=r.topSite().score;
			if(r.perfect() || (PSEUDO_ONLY && score>=maxSwScore)){clearzone=CLEARZONEP;}
			else{
				assert(score<maxSwScore);
				final float cz1blimit=(maxSwScore*CLEARZONE1b_CUTOFF_SCALE-CLEARZONE1b_CUTOFF_FLAT);
				final float cz1climit=(maxSwScore*CLEARZONE1c_CUTOFF_SCALE-CLEARZONE1c_CUTOFF_FLAT);
				if(score>cz1blimit){
//					clearzone=CLEARZONE1;
					clearzone=(int)(((maxSwScore-score)*CLEARZONE1b+(score-cz1blimit)*CLEARZONE1)/(maxSwScore-cz1blimit));
				}else if(score>cz1climit){
//					clearzone=CLEARZONE1b;
					clearzone=(int)(((cz1blimit-score)*CLEARZONE1c+(score-cz1climit)*CLEARZONE1b)/(cz1blimit-cz1climit));
				}else{
					clearzone=CLEARZONE1c;
				}
//				assert(false) : x+", "+cz1blimit+", "+cz1climit+", "+CLEARZONE1b_CUTOFF_FLAT+", "+clearzone;
			}
			
			
//			final int clearzone=r.perfect() ? CLEARZONEP :
//				r.list.get(0).score>=(int)(maxSwScore*CLEARZONE1b_CUTOFF) ? CLEARZONE1 :
//					(r.list.get(0).score>=(int)(maxSwScore*CLEARZONE1c_CUTOFF) ? (CLEARZONE1b_CUTOFF-)CLEARZONE1b : CLEARZONE1c);
			int numBestSites1=Tools.countTopScores(r.sites, clearzone);
			if(numBestSites1>1){
				//Ambiguous alignment
				assert(r.sites.size()>1);
				boolean b=processAmbiguous(r.sites, true, AMBIGUOUS_TOSS, clearzone, SAVE_AMBIGUOUS_XY); //Never gets executed anymore, so always returns true
				r.setAmbiguous(b);
			}else{
				final int lim=(r.perfect() ? (int)(4f*CLEARZONE_LIMIT1e) : score+CLEARZONE1e>=maxSwScore ? 2*CLEARZONE_LIMIT1e : CLEARZONE_LIMIT1e)+1;
				if(r.sites.size()>lim && clearzone<CLEARZONE1e){
					numBestSites1=Tools.countTopScores(r.sites, CLEARZONE1e);
					if(numBestSites1>lim){
						boolean b=processAmbiguous(r.sites, true, AMBIGUOUS_TOSS, clearzone, SAVE_AMBIGUOUS_XY);
						r.setAmbiguous(b);
					}
				}
			}
		}
		
		if(verbose){System.err.println("A: "+r);}
		
		if((SLOW_ALIGN || USE_AFFINE_SCORE) && r.numSites()>0){
			int lim=(int)(maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO);
			if(r.topSite().score<lim){r.sites=null;}
			else{Tools.removeLowQualitySitesUnpaired(r.sites, Tools.min(lim, Tools.max(1, lim-CLEARZONE3)));}
		}
		if(r.numSites()==0){r.sites=null;r.mapScore=0;}
		r.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
		if(PSEUDO_ONLY && r.mapped()){
			assert(r.match!=null) : "Mapped pseudoalignment must retain its polycrystalline trace";
			r.setShortMatch(true);
		}
		assert(Read.CHECKSITES(r, basesM));
		
		if(verbose){System.err.println("B: "+r);}
		
		//Unimportant anomaly due to ambiguous reads that later have low quality sites removed and become unmapped.
//		assert(!r.mapped() || new SamLine(r, 0).toRead(true).ambiguous()==r.ambiguous()) : "\n"+r+"\n\n"+new SamLine(r, 0)+"\n\n"+new SamLine(r, 0).toRead(true)+"\n\n"+
//		"ambi="+ambi+", r.ambiguous()="+r.ambiguous()+", new SamLine(r, 0).toRead(true).ambiguous()="+new SamLine(r, 0).toRead(true).ambiguous()+"\n\n"+
//		"r.mapped="+r.mapped()+", sl.mapped()="+new SamLine(r, 0).mapped()+", sl.toRead(true).mapped()="+new SamLine(r, 0).toRead(true).mapped();
//		assert(r.ambiguous()==ambi) : r;
		
		assert(r.gaps==null || r.gaps[0]==r.start && r.gaps[r.gaps.length-1]==r.stop);
		assert(r.sites==null || r.mapScore>0) : r.sites+", "+r.mapScore+"\n"+r;
		
		if(r.numSites()>1){
			assert(r.topSite().score==r.topSite().slowScore) : "\n"+r.toText(false)+"\n";
			assert(r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
		}
		
		if(verbose){System.err.println("C: "+r);}
		
		//***$
		if(MAKE_MATCH_STRING && r.numSites()>0){
			if(USE_SS_MATCH_FOR_PRIMARY && r.topSite().match!=null){
				r.match=r.topSite().match;
			}else{
				if(r.sites.size()>1){
					assert(r.topSite().score>=r.sites.get(1).score) : "\n"+r.topSite().toText()+"\t<\t"+r.sites.get(1).toText()+"\n"+r.toText(false)+"\n";
				}
				int mapScore=r.mapScore;

				assert(r.mate!=null || r.numSites()==0 || r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";

				if(verbose){System.err.println("D: "+r);}
				
				{
					boolean firstIter=true;
					do{//
						if(!firstIter){
							Shared.sort(r.sites);
							r.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
						}
						genMatchString(r, basesP, basesM, maxImperfectSwScore, maxSwScore, true, true);
						assert(r.mate!=null || r.numSites()==0 || r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
//						TODO: Fix this; it should never happen.
//						if(mapScore>r.mapScore){
//							System.err.println("genMatchString reduced mapping score: "+mapScore+" -> "+r.mapScore+" in read "+r.numericID);
//						}
						if(STRICT_MAX_INDEL && hasLongIndel(r.match, index.maxIndel())){
							SiteScore ss=r.topSite();
							r.mapScore=Tools.min(ss.score, -9999);
							ss.setScore(r.mapScore);
							ss.setSlowPairedScore(ss.score, ss.score);
						}
						r.topSite().setScore(r.topSite().slowScore);
						firstIter=false;
					}while(r.sites.size()>1 && r.topSite().score<r.sites.get(1).score);
				}
				
				if(r.numSites()>1){
					assert(r.topSite().score==r.topSite().slowScore) : "\n"+r.toText(false)+"\n";
					assert(r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
				}

				if(verbose){System.err.println("E: "+r);}
			}
		}
		
		if(r.numSites()>1){
			assert(r.topSite().score==r.topSite().slowScore) : "\n"+r.toText(false)+"\n";
			assert(r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
			removeDuplicateBestSites(r);
		}
		if(r.numSites()>0){r.topSite().match=r.match;}
		
		
		if(r.sites!=null && r.mapScore<=0){//This came from BBMapThreadPacBio; not sure if needed for other modes
			if(!STRICT_MAX_INDEL && !Shared.anomaly){
				System.err.println("Note: Read "+r.id+" failed cigar string generation and will be marked as unmapped.\t"+(r.match==null)+"\t"+r.mapScore+"\t"+r.topSite()+"\t"+new String(r.bases));
				if(MSA.bandwidth>0 || MSA.bandwidthRatio>0 || MSA.flatMode){Shared.anomaly=true;}
			}
			r.mapScore=0;
			r.setMapped(false);
			r.sites=null;
		}
		
		
		
		//This block is to prevent an assertion from firing.  Generally caused by alignment being lost during match generation.
		//TODO: Fix cause.
		if(r.mapScore>0 && r.sites==null){
			if(!Shared.anomaly){System.err.println("Anomaly: mapScore>0 and list==null.\n"+r+"\n");}
			Shared.anomaly=true;
			r.clearMapping();
		}else if(r.mapScore<=0 && r.sites!=null){
			if(BANDWIDTH<1){
				if(!Shared.anomaly){System.err.println("Anomaly1: mapScore<=0 and list!=null.\n"+r+"\n");}
				Shared.anomaly=true;
			}
			r.clearMapping();
		}
		assert(r.sites==null || r.mapScore>0) :
			(PSEUDO_ONLY ? "Pseudoalignment retained sites with nonpositive quick score: "+r :
			(quantumOnly ? "Quantum-only mapping retained sites with nonpositive score: "+r :
			"\nmapScore = "+r.mapScore+"\nread = "+r.toText(false)+"\nscore thresh = "+(-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore))+"\n"+
			"msa unlimited return = "+Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 10), 0))+"\n"+
			"msa limited return = "+Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 10), (-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore))))+"\n\n"+
			"msa vert limit: "+msa.showVertLimit()+"\n\nmsa horz limit: "+msa.showHorizLimit()+"\n\n"));
		
//		assert(r.list==null || r.mapScore>0) : r.mapScore+"\n"+r.list==null ? "null" : r.list.toString();
		
		if((CLEARZONE3>CLEARZONE1 || CLEARZONE3>CLEARZONEP) && r.sites!=null && !r.ambiguous()){
			
			assert(r.mapScore>0);
			float cz3v2=(CLEARZONE3*Tools.min(1.25f, (maxSwScore/(float)r.mapScore)));
			
//			boolean changed=applyClearzone3(r, CLEARZONE3, INV_CLEARZONE3);
			boolean changed=applyClearzone3(r, (int)cz3v2, 1/cz3v2);
			if(changed){
				int minScore=(int)(maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO);
				if(r.mapScore<minScore){
					assert(!r.ambiguous());
					r.setAmbiguous(true);
				}
			}
		}
		
		if(r.ambiguous() && AMBIGUOUS_TOSS){r.sites=null; r.clearSite(); r.setMapped(false);}
		
		if(r.mapped() && r.numSites()>1 && PRINT_SECONDARY_ALIGNMENTS){
			ensureMatchStringsOnSiteScores(r, basesM, maxImperfectSwScore, maxSwScore);
			assert(Read.CHECKSITES(r, basesM));
		}
		
		assert(checkTopSite(r));
		if(!quantumOnly && !PSEUDO_ONLY && r.mapped() && (LOCAL_ALIGN || r.containsXYC())){
			msa.toLocalAlignment(r, r.topSite(), basesM, r.containsXYC() ? 1 : LOCAL_ALIGN_TIP_LENGTH, LOCAL_ALIGN_MATCH_POINT_RATIO);
			assert(Read.CHECKSITES(r, basesM));
		}
		
		if(r.numSites()==0 || (!r.ambiguous() && r.mapScore<maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO)){
			r.clearMapping();
		}

//		assert(false) : "\n\n"+r.sites+"\n\n"+r.toSam()+"\n\n"+r+"\n\n";
		
		postFilterRead(r, basesM, maxImperfectSwScore, maxSwScore);
		if(MAKE_MATCH_STRING){ensureMatchStringOnPrimary(r, basesM, maxImperfectSwScore, maxSwScore);}
		
		if(PENALIZE_AMBIG){
			int penalty=calcTipScorePenalty(r, maxSwScore, 7);
			applyScorePenalty(r, penalty);
		}
		
//		if(r.ambiguous() && r.sites!=null){
//			r.setAmbiguous(false);
//			r.mapScore/=3;
//			for(SiteScore ss : r.sites){
//				ss.slowScore/=3;
//				ss.score/=3;
//			}
//		}
		
//		//Penalize quality score of long deletions, as they are less likely to be accurate.  Did not seem to be useful.
//		if(r.mapped() && !r.ambiguous()){
//			int delta=absdif(r.start+r.length()-1, r.stop);
//			if(delta>100){
//				float penalty2=0.0004f*(500f*delta)/(500f+delta);
//				r.mapScore=(int)(r.mapScore*(1-penalty2));
//			}
//		}
		
	}

	/** Finalize the selected single-read attempt after its accounting is committed. */
	private void finishRead(final Read r, final SingleAttemptAccounting attempt){
		assert(attempt.committed && !attempt.rejected) :
			"Finalize only one selected, nonrejected single-read attempt";
		if(CALC_STATISTICS){calcStatistics1(r,attempt.max,attempt.quick);}
	}
	
	
	/** Returns number of perfect pairs */
	@Override
	public int pairSiteScoresInitial(Read r, Read r2, boolean trim){
		
		if(r.numSites()<1 || r2.numSites()<1){return 0;}
		
		SiteScore.PCOMP.sort(r.sites);
		SiteScore.PCOMP.sort(r2.sites);
		
		for(SiteScore ss : r.sites){ss.setPairedScore(0);}
		for(SiteScore ss : r2.sites){ss.setPairedScore(0);}
		
//		ArrayList<SiteScorePair> pairs=new ArrayList<SiteScorePair>(Tools.min(8, Tools.min(r.list.size(), r2.list.size())));

		int maxPairedScore1=-1;
		int maxPairedScore2=-1;
		
		
//		for(SiteScore ss : r.list){
//			System.out.println(ss.toText());
//		}
		
//		int i=0, j=0;
		final int ilimit=r.sites.size()-1;
		final int jlimit=r2.sites.size()-1;
		final int maxReadLen=Tools.max(r.length(), r2.length());
		
//		final int outerDistLimit=MIN_PAIR_DIST+r.length()+r2.length();
		final int outerDistLimit=(Tools.max(r.length(), r2.length())*(OUTER_DIST_MULT))/OUTER_DIST_DIV;//-(SLOW_ALIGN ? 100 : 0);
		final int innerDistLimit=MAX_PAIR_DIST;//+(FIND_TIP_DELETIONS ? TIP_DELETION_SEARCH_RANGE : 0);
		final int expectedFragLength=AVERAGE_PAIR_DIST+r.length()+r2.length();
		
		int numPerfectPairs=0;
		
		for(int i=0, j=0; i<=ilimit && j<=jlimit; i++){
			SiteScore ss1=r.sites.get(i);
			SiteScore ss2=r2.sites.get(j);
			
			while(j<jlimit && (ss2.chrom<ss1.chrom || (ss2.chrom==ss1.chrom && ss1.start-ss2.stop>innerDistLimit))){
				j++;
				ss2=r2.sites.get(j);
			}

			for(int k=j; k<=jlimit; k++){
				ss2=r2.sites.get(k);

				if(ss2.chrom>ss1.chrom){break;}
				if(ss2.start-ss1.stop>innerDistLimit){break;}

//				int dist=0;
//
//				if(ss1.start<=ss2.start){
//					dist=ss2.start-ss1.stop;
//				}else if(ss1.start>ss2.start){
//					dist=ss1.start-ss2.stop;
//				}
				
				
//				int innerdist=0;
//				int outerdist=0;
//
//				if(ss1.start<=ss2.start){
//					innerdist=ss2.start-ss1.stop;
//					outerdist=ss2.stop-ss1.start;
//				}else if(ss1.start>ss2.start){
//					innerdist=ss1.start-ss2.stop;
//					outerdist=ss1.stop-ss2.start;
//				}
				
				final int innerdist, outerdist;
				//assert(!SAME_STRAND_PAIRS) : "TODO";
				
				if(REQUIRE_CORRECT_STRANDS_PAIRS){
					if(ss1.strand!=ss2.strand){
						if(ss1.strand==Shared.PLUS){
							innerdist=ss2.start-ss1.stop;
							outerdist=ss2.stop-ss1.start;
						}else{
							innerdist=ss1.start-ss2.stop;
							outerdist=ss1.stop-ss2.start;
						}
					}else{
						if(ss1.start<=ss2.start){
							innerdist=ss2.start-ss1.stop;
							outerdist=ss2.stop-ss1.start;
						}else{
							innerdist=ss1.start-ss2.stop;
							outerdist=ss1.stop-ss2.start;
						}
					}
				}else{
					if(ss1.start<=ss2.start){
						innerdist=ss2.start-ss1.stop;
						outerdist=ss2.stop-ss1.start;
					}else{
						innerdist=ss1.start-ss2.stop;
						outerdist=ss1.stop-ss2.start;
					}
				}
				
				assert(outerdist>=innerdist);

				if(outerdist>=outerDistLimit && innerdist<=innerDistLimit){
					
					boolean strandOK=((ss1.strand==ss2.strand)==SAME_STRAND_PAIRS);

					if(strandOK || !REQUIRE_CORRECT_STRANDS_PAIRS){
						
						boolean paired1=false, paired2=false;
						
						int deviation=absdif(AVERAGE_PAIR_DIST, innerdist);

						final int pairedScore1;
						final int pairedScore2;
						if(strandOK){
//							pairedScore1=ss1.score+ss2.score/2;
//							pairedScore2=ss2.score+ss1.score/2;
							
							pairedScore1=ss1.score+1+Tools.max(1, ss2.score/2-(((deviation)*ss2.score)/(32*expectedFragLength+100)));
							pairedScore2=ss2.score+1+Tools.max(1, ss1.score/2-(((deviation)*ss1.score)/(32*expectedFragLength+100)));
						}else{//e.g. a junction
							pairedScore1=ss1.score+Tools.max(0, ss2.score/16);
							pairedScore2=ss2.score+Tools.max(0, ss1.score/16);
						}

						if(pairedScore1>ss1.pairedScore){
							paired1=true;
							ss1.setPairedScore(Tools.max(ss1.pairedScore, pairedScore1));
							maxPairedScore1=Tools.max(ss1.score, maxPairedScore1);
							//						System.out.println("Paired "+ss1.toText()+" with "+ss2.toText());
						}else{
							//						System.out.println(ss1.toText()+" already paired.");
						}
						if(pairedScore2>ss2.pairedScore){
							paired2=true;
							ss2.setPairedScore(Tools.max(ss2.pairedScore, pairedScore2));
							maxPairedScore2=Tools.max(ss2.score, maxPairedScore2);
						}
						
						if(paired1 && paired2 && outerdist>=maxReadLen && deviation<=expectedFragLength && ss1.perfect && ss2.perfect){
							numPerfectPairs++; //Lower bound.  Some perfect pairs may be the same.
						}
						
//						ss1.pairedScore=Tools.max(ss1.pairedScore, pairedScore1);
//						ss2.pairedScore=Tools.max(ss2.pairedScore, pairedScore2);
//						maxPairedScore1=Tools.max(ss1.score, maxPairedScore1);
//						maxPairedScore2=Tools.max(ss2.score, maxPairedScore2);
					}
				}
			}
			
		}
		
		
		
		for(SiteScore ss : r.sites){
			if(ss.pairedScore>ss.score){ss.setScore(ss.pairedScore);}
			else{assert(ss.pairedScore==0);}
//			ss.score=ss.pairedScore=Tools.max(ss.pairedScore, ss.score);
		}
		for(SiteScore ss : r2.sites){
			if(ss.pairedScore>ss.score){ss.setScore(ss.pairedScore);}
			else{assert(ss.pairedScore==0);}
//			ss.score=ss.pairedScore=Tools.max(ss.pairedScore, ss.score);
		}
		
		if(trim){
			if(numPerfectPairs>0){
//				System.out.print(".");
				Tools.trimSitesBelowCutoff(r.sites, (int)(maxPairedScore1*.94f), false, true, 1, MAX_TRIM_SITES_TO_RETAIN);
				Tools.trimSitesBelowCutoff(r2.sites, (int)(maxPairedScore2*.94f), false, true, 1, MAX_TRIM_SITES_TO_RETAIN);
			}else{
				if(r.sites.size()>4){
					Tools.trimSitesBelowCutoff(r.sites, (int)(maxPairedScore1*.9f), true, true, 1, MAX_TRIM_SITES_TO_RETAIN);
				}
				if(r2.sites.size()>4){
					Tools.trimSitesBelowCutoff(r2.sites, (int)(maxPairedScore2*.9f), true, true, 1, MAX_TRIM_SITES_TO_RETAIN);
				}
			}
		}
		
//		if(pairs.isEmpty()){return null;}
//
//		ArrayList<SiteScore> temp=new ArrayList<SiteScore>(Tools.max(r.list.size(), r2.list.size()));
//
//		for(SiteScore ss : r.list){
//			if(ss.score>maxPairedScore1){temp.add(ss);}
//		}
//		for(SiteScorePair ssp : pairs){
//			temp.add(ssp.a);
//		}
//		r.list.clear();
//		r.list.addAll(temp);
//
//		for(SiteScore ss : r2.list){
//			if(ss.score>maxPairedScore2){temp.add(ss);}
//		}
//		for(SiteScorePair ssp : pairs){
//			temp.add(ssp.b);
//		}
//		r2.list.clear();
//		r2.list.addAll(temp);
//
//		return pairs;
		
		return numPerfectPairs;
	}
	
	
	// Per-worker configuration supplied by the mapper; legacy constructors default off.
	private final boolean hybridPair;
	private final PairAttemptAccounting secondAttempt;
	private final PairSearchState entryState;
	private final PairSearchState firstState;
	@Override
	public void processReadPair(final Read r, final byte[] basesM1, final byte[] basesM2){
		if(idmodulo>1 && r.numericID%idmodulo!=1){return;}
		assert(r.mate!=null) : "Paired search requires a mate before either attempt";
		readsUsed1++;readsUsed2++;
		final PairAttemptAccounting first=pairAttempt;first.reset();
		if(hybridMatchCache!=null){hybridMatchCache.clear();}
		final boolean generalHybrid=(hybridMaxIndelConfig!=null);
		if(!hybridPair && !generalHybrid){
			searchAndScorePair(r,basesM1,basesM2,first);
			first.commit(this);
			if(!first.rejected){finishReadPair(r,basesM1,basesM2,first);}
			return;
		}
		final int lowPrimary=generalHybrid ? hybridMaxIndelConfig.lowPrimary : 50;
		final int lowSum=generalHybrid ? hybridMaxIndelConfig.lowSum : 100;
		final int highPrimary=generalHybrid ? hybridMaxIndelConfig.highPrimary : 16000;
		final int highSum=generalHybrid ? hybridMaxIndelConfig.highSum : 32000;
		if(index.maxIndel()!=lowPrimary || index.maxIndel2()!=lowSum ||
				QUICK_MATCH_STRINGS || STRICT_MAX_INDEL || PRINT_SECONDARY_ALIGNMENTS || BBIndex.PERFECTMODE || BBIndex.SEMIPERFECTMODE){
			throw new IllegalArgumentException("Paired hybrid worker entered with incompatible options or bounds; expected "+
				lowPrimary+"/"+lowSum+", got "+index.maxIndel()+"/"+index.maxIndel2());
		}
		final long ac0=initialSiteSum1,ac1=initialSiteSum2,ac2=postTrimSiteSum1,ac3=postTrimSiteSum2,ac4=postRescueSiteSum1,ac5=postRescueSiteSum2,ac6=mapped1,ac7=mapped2,ac8=lowQualityReadsDiscarded1,ac9=lowQualityReadsDiscarded2,ac10=lowQualityBasesDiscarded1,ac11=lowQualityBasesDiscarded2,ac12=readsUsed1,ac13=readsUsed2;
		entryState.capture(r);
		try{
			searchAndScorePair(r,basesM1,basesM2,first);firstState.capture(r);
			final HybridPairPolicy.Observation observation=first.rejected ? null : HybridPairPolicy.observe(this,r,basesM1,basesM2,first.imperfect1,first.max1,first.imperfect2,first.max2,false);
			PairAttemptAccounting chosen=first;
			final boolean route=observation!=null && (generalHybrid ?
				(observation.flags&HybridPairPolicy.HALF_ERRORS)!=0 : observation.route());
			if(route){
				if(generalHybrid){
					hybridMaxIndelStats.retryAttempted();
					hybridMaxIndelStats.pairReasons(observation.flags);
				}
				if(hybridMatchCache!=null){hybridMatchCache.clear();}
				entryState.restore();r.sites=null;r.mate.sites=null;
				index.setRuntimeIndelLimits(highPrimary,highSum);secondAttempt.reset();
				searchAndScorePair(r,basesM1,basesM2,secondAttempt);
				final String geometry=secondAttempt.rejected ? "REJECTED" : HybridPairPolicy.observe(this,r,basesM1,basesM2,secondAttempt.imperfect1,secondAttempt.max1,secondAttempt.imperfect2,secondAttempt.max2,true).geometry;
				final boolean mapqAccepted=!generalHybrid || retryPairMapqAccepted(r,hybridMaxIndelConfig.retryMinMapq);
				if(geometry.equals("NO_CONFLICT") && mapqAccepted){
					first.discard();chosen=secondAttempt;
					if(generalHybrid){hybridMaxIndelStats.wideSelected();hybridMaxIndelStats.pairOutcome(observation.flags,true);}
				}
				else{
					if(generalHybrid && geometry.equals("NO_CONFLICT") && !mapqAccepted){hybridMaxIndelStats.mapqRejected();}
					if(hybridMatchCache!=null){hybridMatchCache.clear();}
					secondAttempt.discard();firstState.restore();index.setRuntimeIndelLimits(lowPrimary,lowSum);
					if(generalHybrid){hybridMaxIndelStats.lowRestored();hybridMaxIndelStats.pairOutcome(observation.flags,false);}
				}
			}
			chosen.commit(this);
			assert(initialSiteSum1-ac0==chosen.initialSiteSum1 && initialSiteSum2-ac1==chosen.initialSiteSum2) : "Only selected attempt may commit initialSiteSum";
			assert(postTrimSiteSum1-ac2==chosen.postTrimSiteSum1 && postTrimSiteSum2-ac3==chosen.postTrimSiteSum2) : "Only selected attempt may commit postTrimSiteSum";
			assert(postRescueSiteSum1-ac4==chosen.postRescueSiteSum1 && postRescueSiteSum2-ac5==chosen.postRescueSiteSum2) : "Only selected attempt may commit postRescueSiteSum";
			assert(mapped1-ac6==chosen.mapped1 && mapped2-ac7==chosen.mapped2) : "Only selected attempt may commit mapped counts";
			assert(lowQualityReadsDiscarded1-ac8==(chosen.rejected ? 1 : 0) && lowQualityReadsDiscarded2-ac9==(chosen.rejected ? 1 : 0)) : "Rejected-read counts belong only to selected attempt";
			assert(lowQualityBasesDiscarded1-ac10==(chosen.rejected ? chosen.rejectedBases1 : 0) && lowQualityBasesDiscarded2-ac11==(chosen.rejected ? chosen.rejectedBases2 : 0)) : "Rejected-base counts belong only to selected attempt";
			assert(readsUsed1==ac12 && readsUsed2==ac13) : "Retry must not count input reads again";
			if(!chosen.rejected){
				if(hybridMatchCache!=null){hybridMatchCache.replaying=true;}
				try{finishReadPair(r,basesM1,basesM2,chosen);}
				finally{if(hybridMatchCache!=null){hybridMatchCache.replaying=false;}}
			}
		}finally{
			try{index.setRuntimeIndelLimits(lowPrimary,lowSum);}finally{
				if(hybridMatchCache!=null){hybridMatchCache.clear();}
				entryState.clear();firstState.clear();
			}
		}
	}

	private void finishReadPair(final Read r, final byte[] basesM1, final byte[] basesM2, final PairAttemptAccounting attempt){
		assert(!attempt.rejected && attempt.committed) : "Finalize only the chosen committed attempt";
		final Read r2=r.mate;
		final byte[] basesP1=r.bases, basesP2=r2.bases;
		final int len1=(basesP1==null ? 0 : basesP1.length), len2=(basesP2==null ? 0 : basesP2.length);
		final int maxPossibleQuickScore1=attempt.quick1, maxPossibleQuickScore2=attempt.quick2;
		final int maxSwScore1=attempt.max1, maxSwScore2=attempt.max2;
		final int maxImperfectSwScore1=attempt.imperfect1, maxImperfectSwScore2=attempt.imperfect2;

		if(!PSEUDO_ONLY && (SLOW_ALIGN || USE_AFFINE_SCORE)){
			r.setPerfectFlag(maxSwScore1);
			r2.setPerfectFlag(maxSwScore2);
//			assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		}
		

		if(r.numSites()>1){
			final int clearzone=(r.perfect() || (PSEUDO_ONLY && r.topSite().score>=maxSwScore1)) ? CLEARZONEP :
				r.topSite().score>=(int)(maxSwScore1*CLEARZONE1b_CUTOFF_SCALE-CLEARZONE1b_CUTOFF_FLAT) ? CLEARZONE1 :
					(r.topSite().score>=(int)(maxSwScore1*CLEARZONE1c_CUTOFF_SCALE-CLEARZONE1c_CUTOFF_FLAT) ? CLEARZONE1b : CLEARZONE1c);
			int numBestSites1=Tools.countTopScores(r.sites, clearzone);
			if(numBestSites1>1){
				//Ambiguous alignment
				assert(r.sites.size()>1);
				
				boolean b=processAmbiguous(r.sites, true, AMBIGUOUS_TOSS, clearzone, SAVE_AMBIGUOUS_XY);
				r.setAmbiguous(b);
			}
//			assert(Read.CHECKSITES(r, basesM1));
		}

		if(r2.numSites()>1){
			final int clearzone=(r2.perfect() || (PSEUDO_ONLY && r2.topSite().score>=maxSwScore2)) ? CLEARZONEP :
				r2.topSite().score>=(int)(maxSwScore2*CLEARZONE1b_CUTOFF_SCALE-CLEARZONE1b_CUTOFF_FLAT) ? CLEARZONE1 :
					(r2.topSite().score>=(int)(maxSwScore2*CLEARZONE1c_CUTOFF_SCALE-CLEARZONE1c_CUTOFF_FLAT) ? CLEARZONE1b : CLEARZONE1c);
			int numBestSites2=Tools.countTopScores(r2.sites, clearzone);
			if(numBestSites2>1){
				//Ambiguous alignment
				assert(r2.sites.size()>1);
				
				boolean b=processAmbiguous(r2.sites, false, AMBIGUOUS_TOSS, clearzone, SAVE_AMBIGUOUS_XY);
				r2.setAmbiguous(b);
			}
//			assert(Read.CHECKSITES(r2, basesM2));
		}
		if(verbose){System.err.println("\nAfter ambiguous removal:\nRead1:\t"+r+"\nRead2:\t"+r2);}
		
		if(r.numSites()>0 && r2.numSites()>0){
			SiteScore ss1=r.topSite();
			SiteScore ss2=r2.topSite();
			if(canPair(ss1, ss2, len1, len2, REQUIRE_CORRECT_STRANDS_PAIRS, SAME_STRAND_PAIRS, MAX_PAIR_DIST)){
				assert(SLOW_ALIGN ? ss1.pairedScore>ss1.slowScore : ss1.pairedScore>ss1.quickScore) :
					"\n"+ss1.toText()+"\n"+ss2.toText()+"\n"+r.toText(false)+"\n"+r2.toText(false)+"\n\n"+
						r.mapped()+", "+r.paired()+", "+r.strand()+", "+r.ambiguous()+"\n\n"+r2.mapped()+", "+r2.paired()+", "+r2.strand()+", "+r2.ambiguous()+"\n\n";
				assert(SLOW_ALIGN ? ss2.pairedScore>ss2.slowScore : ss2.pairedScore>ss2.quickScore) :
					"\n"+ss1.toText()+"\n"+ss2.toText()+"\n"+r.toText(false)+"\n"+r2.toText(false)+"\n\n";
				r.setPaired(true);
				r.mate.setPaired(true);
			}
		}

		if(r.numSites()==0){r.sites=null;r.mapScore=0;}
		if(r2.numSites()==0){r2.sites=null;r2.mapScore=0;}
		
		r.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
		r2.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
		if(PSEUDO_ONLY && r.mapped()){
			assert(r.match!=null) : "Mapped pseudoalignment must retain its polycrystalline trace";
			r.setShortMatch(true);
		}
		if(PSEUDO_ONLY && r2.mapped()){
			assert(r2.match!=null) : "Mapped pseudoalignment must retain its polycrystalline trace";
			r2.setShortMatch(true);
		}
		if(KILL_BAD_PAIRS){
			if(r.isBadPair(REQUIRE_CORRECT_STRANDS_PAIRS, SAME_STRAND_PAIRS, MAX_PAIR_DIST)){
				int x=r.mapScore/len1;
				int y=r2.mapScore/len2;
				if(x>=y){
					r2.clearAnswers(false);
				}else{
					r.clearAnswers(false);
				}
			}
		}
		if(verbose){System.err.println("\nAfter bad pair removal:\nRead1:\t"+r+"\nRead2:\t"+r2);}
		
		assert(r.sites==null || r.mapScore>0) : r.mapScore+"\n"+r.toText(false)+"\n\n"+r2.toText(false)+"\n";
		assert(r2.sites==null || r2.mapScore>0) : r2.mapScore+"\n"+r.toText(false)+"\n\n"+r2.toText(false)+"\n";
		if(MAKE_MATCH_STRING){
			if(r.numSites()>0){
				if(USE_SS_MATCH_FOR_PRIMARY && r.topSite().match!=null){
					r.match=r.topSite().match;
				}else{
					genMatchString(r, basesP1, basesM1, maxImperfectSwScore1, maxSwScore1, false, false);
					
					if(STRICT_MAX_INDEL && r.mapped()){
						if(hasLongIndel(r.match, index.maxIndel())){
							r.clearMapping();
							r2.setPaired(false);
						}
					}
				}
//				assert(Read.CHECKSITES(r, basesM1));
			}
			if(r2.numSites()>0){
				if(USE_SS_MATCH_FOR_PRIMARY && r2.topSite().match!=null){
					r2.match=r2.topSite().match;
				}else{
					genMatchString(r2, basesP2, basesM2, maxImperfectSwScore2, maxSwScore2, false, false);
					if(STRICT_MAX_INDEL && r2.mapped()){
						if(hasLongIndel(r2.match, index.maxIndel())){
							r2.clearMapping();
							r.setPaired(false);
						}
					}
				}
//				assert(Read.CHECKSITES(r2, basesM2));
			}
		}
		
		assert(checkTopSite(r)); // TODO remove this
		if(verbose){
			System.err.println("\nFinal:\nRead1:\t"+r+"\nRead2:\t"+r2);
			if(r.match!=null && r.shortmatch()){r.toLongMatchString(false);}
			if(r2.match!=null && r2.shortmatch()){r2.toLongMatchString(false);}
		}
		
		//Block to prevent assertion from firing.  Generally caused by alignment being lost during match generation.  TODO: Fix cause.
		if(r.mapScore>0 && r.sites==null){
			if(!Shared.anomaly){System.err.println("Anomaly: mapScore>0 and list==null.\n"+r+"\n");}
			Shared.anomaly=true;
			r.clearMapping();
			r2.setPaired(false);
		}else if(r.mapScore<=0 && r.sites!=null){
			if(!STRICT_MAX_INDEL && !Shared.anomaly){System.err.println("Anomaly2: mapScore<=0 and list!=null.\n"+r+"\n");}
			Shared.anomaly=true;
			r.clearMapping();
			r2.setPaired(false);
		}
		assert(checkTopSite(r)); // TODO remove this
		//Block to prevent assertion from firing.  Generally caused by alignment being lost during match generation.  TODO: Fix cause.
		if(r2.mapScore>0 && r2.sites==null){
			if(!Shared.anomaly){System.err.println("Anomaly: mapScore>0 and list==null.\n"+r+"\n");}
			Shared.anomaly=true;
			r2.clearMapping();
			r.setPaired(false);
		}else if(r2.mapScore<=0 && r2.sites!=null){
			if(!STRICT_MAX_INDEL && !Shared.anomaly){System.err.println("Anomaly3: mapScore<=0 and list!=null.\n"+r+"\n");}
			Shared.anomaly=true;
			r2.clearMapping();
			r.setPaired(false);
		}
		
		assert(r.sites==null || r.mapScore>0) :
			(PSEUDO_ONLY ? "Pseudoalignment paired read retained sites with nonpositive quick score: "+r :
			(quantumOnly ? "Quantum-only paired read retained sites with nonpositive score: "+r :
			r.mapScore+"\t"+r.sites+"\n"+(-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore1))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), 0))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), (-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore1))))+"\n\n"+
			msa.showVertLimit()+"\n\n"+msa.showHorizLimit()+"\n\n"+r+"\n\n"+r2+"\n\n"));
		assert(r2.sites==null || r2.mapScore>0) :
			(PSEUDO_ONLY ? "Pseudoalignment paired read retained sites with nonpositive quick score: "+r2 :
			(quantumOnly ? "Quantum-only paired read retained sites with nonpositive score: "+r2 :
			r2.mapScore+"\t"+r2.sites+"\n"+(-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore2))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r2.strand()==Shared.PLUS ? r2.bases :
			AminoAcid.reverseComplementBases(r2.bases), r2.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), 0))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r2.strand()==Shared.PLUS ? r2.bases :
			AminoAcid.reverseComplementBases(r2.bases), r2.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), (-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore2))))+"\n\n"+
			msa.showVertLimit()+"\n\n"+msa.showHorizLimit()+"\n\n"+r+"\n\n"+r2+"\n\n"));
		
		assert(!r.mapped() || !MAKE_MATCH_STRING || r.match!=null) : "Note that sometimes, VERY RARELY, match string generation fails.";
		assert(checkTopSite(r)); // TODO remove this
		removeDuplicateBestSites(r);
		removeDuplicateBestSites(r2);
		
		if(DYNAMIC_INSERT_LENGTH && numMated>1000 && r.paired()){
			AVERAGE_PAIR_DIST=(int)(innerLengthSum*1f/numMated);
		}
		assert(checkTopSite(r)); // TODO remove this
		if(r.ambiguous() && AMBIGUOUS_TOSS){
			if(r.sites!=null){r.sites=null;}
			r.clearSite();
			r.setMapped(false);
			r.setPaired(false);
			r2.setPaired(false);
		}else if(r.mapped() && r.numSites()>1 && PRINT_SECONDARY_ALIGNMENTS){
			ensureMatchStringsOnSiteScores(r, basesM1, maxImperfectSwScore1, maxSwScore1);
			assert(Read.CHECKSITES(r, basesM1));
		}
		if(r2.ambiguous() && AMBIGUOUS_TOSS){
			if(r2.sites!=null){r2.sites=null;}
			r2.clearSite();
			r2.setMapped(false);
			r.setPaired(false);
			r2.setPaired(false);
		}else if(r2.mapped() && r2.numSites()>1 && PRINT_SECONDARY_ALIGNMENTS){
			ensureMatchStringsOnSiteScores(r2, basesM2, maxImperfectSwScore2, maxSwScore2);
			assert(Read.CHECKSITES(r2, basesM2));
		}
//		assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		
		assert(checkTopSite(r));
		if(!quantumOnly && !PSEUDO_ONLY && r.mapped() && (LOCAL_ALIGN || r.containsXYC())){
			final SiteScore ss=r.topSite();
			ss.match=r.match;
			msa.toLocalAlignment(r, ss, basesM1, r.containsXYC() ? 1 : LOCAL_ALIGN_TIP_LENGTH, LOCAL_ALIGN_MATCH_POINT_RATIO);
//			System.err.println("\n\n*********\n\n"+r+"\n\n*********\n\n");
//			assert(Read.CHECKSITES(r, basesM1)); //TODO: This can fail; see bug#0001
		}
//		assert(false) : r.mapped()+", "+LOCAL_ALIGN+", "+r.containsXYC()+", "+new String(r.match);
		
		assert(checkTopSite(r2));
		if(!quantumOnly && !PSEUDO_ONLY && r2.mapped() && (LOCAL_ALIGN || r2.containsXYC())){
			final SiteScore ss=r2.topSite();
			ss.match=r2.match;
			msa.toLocalAlignment(r2, ss, basesM2, r2.containsXYC() ? 1 : LOCAL_ALIGN_TIP_LENGTH, LOCAL_ALIGN_MATCH_POINT_RATIO);
//			assert(Read.CHECKSITES(r2, basesM2)); //TODO: This can fail; see bug#0001
		}
		
		postFilterRead(r, basesM1, maxImperfectSwScore1, maxSwScore1);
		postFilterRead(r2, basesM2, maxImperfectSwScore2, maxSwScore2);
		if(MAKE_MATCH_STRING){
			ensureMatchStringOnPrimary(r, basesM1, maxImperfectSwScore1, maxSwScore1);
			ensureMatchStringOnPrimary(r2, basesM2, maxImperfectSwScore2, maxSwScore2);
		}
		
		if(CALC_STATISTICS){
			calcStatistics1(r, maxSwScore1, maxPossibleQuickScore1);
			calcStatistics2(r2, maxSwScore2, maxPossibleQuickScore2);
		}
	}
	

	/** One search attempt; mutates both Reads. This is not a rollback transaction. */
	private void searchAndScorePair(final Read r, final byte[] basesM1, final byte[] basesM2, final PairAttemptAccounting attempt){
		assert(r.mate!=null && !attempt.committed) : "Attempt accounting must be reset before paired search";
		final Read r2=r.mate;
		final byte[] basesP1=r.bases, basesP2=r2.bases;
		final int len1=(basesP1==null ? 0 : basesP1.length), len2=(basesP2==null ? 0 : basesP2.length);
		final int maxPossibleQuickScore1=attempt.quick1=quickMap(r, basesM1);
		final int maxPossibleQuickScore2=attempt.quick2=quickMap(r2, basesM2);
		
		if(verbose){
			System.err.println("\nAfter quick map:\nRead1:\t"+r+"\nRead2:\t"+r.mate);
		}
		
		if(maxPossibleQuickScore1<0 && maxPossibleQuickScore2<0){
			r.sites=null;
			r2.sites=null;
			attempt.rejected=true;
			attempt.rejectedBases1=len1;
			r.setDiscarded(true);
			attempt.rejected=true;
			attempt.rejectedBases2=len2;
			r2.setDiscarded(true);
			return;
		}
		
		//Not really needed due to subsumption
//		Tools.mergeDuplicateSites(r.list);
//		Tools.mergeDuplicateSites(r2.list);
		
		attempt.initialSiteSum1=r.numSites();
		attempt.initialSiteSum2=r2.numSites();
		
		//TODO: Fix this.  This is a workaround for an assertion error counting the number of reads used.
		//Discards need to be tracked separately for each end.
//		if(maxPossibleQuickScore2<0){lowQualityReadsDiscarded--;}
		
		final int maxSwScore1=attempt.max1=(PSEUDO_ONLY ? Tools.max(1, maxPossibleQuickScore1) :
				(quantumOnly ? quantumMaxScore(len1) : msa.maxQuality(len1)));
		final int maxImperfectSwScore1=attempt.imperfect1=(PSEUDO_ONLY ? maxSwScore1 :
				(quantumOnly ? maxSwScore1-QUANTUM_EDIT_COST : msa.maxImperfectScore(len1)));
		final int maxSwScore2=attempt.max2=(PSEUDO_ONLY ? Tools.max(1, maxPossibleQuickScore2) :
				(quantumOnly ? quantumMaxScore(len2) : msa.maxQuality(len2)));
		final int maxImperfectSwScore2=attempt.imperfect2=(PSEUDO_ONLY ? maxSwScore2 :
				(quantumOnly ? maxSwScore2-QUANTUM_EDIT_COST : msa.maxImperfectScore(len2)));
		
		pairSiteScoresInitial(r, r2, TRIM_LIST);
		if(verbose){System.err.println("\nAfter initial pair:\nRead1:\t"+r+"\nRead2:\t"+r2);}
		
		if(TRIM_LIST){

			if(MIN_TRIM_SITES_TO_RETAIN_PAIRED>1){
				if(r.numSites()>MIN_TRIM_SITES_TO_RETAIN_PAIRED){Shared.sort(r.sites);}
				if(r2.numSites()>MIN_TRIM_SITES_TO_RETAIN_PAIRED){Shared.sort(r2.sites);}
			}
			
			trimList(r.sites, true, maxSwScore1, false, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
			trimList(r2.sites, true, maxSwScore2, false, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
		}
		attempt.postTrimSiteSum1=r.numSites();
		attempt.postTrimSiteSum2=r2.numSites();
		
		{//Reset score to non-paired score
			if(r.sites!=null){
				for(SiteScore ss : r.sites){
					assert(ss.slowScore<=ss.quickScore);
					ss.setScore(ss.quickScore);
				}
			}
			if(r2.sites!=null){
				for(SiteScore ss : r2.sites){
					assert(ss.slowScore<=ss.quickScore);
					ss.setScore(ss.quickScore);
				}
			}
		}
		
		if(verbose){System.err.println("\nAfter trim:\nRead1:\t"+r.sites+"\nRead2:\t"+r2.sites);}
		
//		assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		
		if(SLOW_ALIGN && !PSEUDO_ONLY){
			
			if(r.numSites()>0){
				
				int numNearPerfectScores1=(quantumOnly ?
						scoreNoIndelsQuantum(r, basesP1, basesM1, maxSwScore1, maxImperfectSwScore1) :
						scoreNoIndels(r, basesP1, basesM1, maxSwScore1, maxImperfectSwScore1));
				Shared.sort(r.sites); //Puts higher scores first to better trigger the early exit based on perfect scores
				
				if(numNearPerfectScores1<1){
					if(!quantumOnly && FIND_TIP_DELETIONS){
						findTipDeletions(r, basesP1, basesM1, maxSwScore1, maxImperfectSwScore1);
					}
				}
				
				//TODO:
				//Note scoreSlow can be skipped under this circumstance:
				//When rescue is disabled, numNearPerfectScores>0, and there are no paired sites.
				scoreSlow(r.sites, basesP1, basesM1, maxSwScore1, maxImperfectSwScore1);
				if(STRICT_MAX_INDEL){
					int removed=removeLongIndels(r.sites, index.maxIndel());
					if(r.numSites()==0){r.clearMapping();}
				}
				Tools.mergeDuplicateSites(r.sites, true, true);
				removeQuantumRejectedSites(r);
			}
			
			if(r2.numSites()>0){
				int numNearPerfectScores2=(quantumOnly ?
						scoreNoIndelsQuantum(r2, basesP2, basesM2, maxSwScore2, maxImperfectSwScore2) :
						scoreNoIndels(r2, basesP2, basesM2, maxSwScore2, maxImperfectSwScore2));
				Shared.sort(r2.sites); //Puts higher scores first to better trigger the early exit based on perfect scores
				
				if(numNearPerfectScores2<1){
					if(!quantumOnly && FIND_TIP_DELETIONS){
						findTipDeletions(r2, basesP2, basesM2, maxSwScore2, maxImperfectSwScore2);
					}
				}
				
				scoreSlow(r2.sites, basesP2, basesM2, maxSwScore2, maxImperfectSwScore2);
				if(STRICT_MAX_INDEL){
					int removed=removeLongIndels(r2.sites, index.maxIndel());
					if(r2.numSites()<1){r2.clearMapping();}
				}
				Tools.mergeDuplicateSites(r2.sites, true, true);
				removeQuantumRejectedSites(r2);
			}
			
			
			if(verbose){System.err.println("\nAfter slow align:\nRead1:\t"+r+"\nRead2:\t"+r2);}
			assert(Read.CHECKSITES(r, basesM1, false) && Read.CHECKSITES(r2, basesM2, false));
			
			if(DO_RESCUE){
				int unpaired1=0;
				int unpaired2=0;
				if(r.sites!=null){
					for(SiteScore ss : r.sites){
						assert(ss.pairedScore<1 || ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) :
							"\n"+ss.toText()+"\n"+r.toText(false)+"\n\n"+r.toFastq()+"\n"+r2.toFastq()+"\napd="+AVERAGE_PAIR_DIST+"\n";
						if(ss.pairedScore==0){unpaired1++;}
					}
				}
				if(r2.sites!=null){
					for(SiteScore ss : r2.sites){
						assert(ss.pairedScore<1 || ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) :
							"\n"+ss.toText()+"\n"+r2.toText(false)+"\n\n"+r.toFastq()+"\n"+r2.toFastq()+"\napd="+AVERAGE_PAIR_DIST+"\n";
						if(ss.pairedScore==0){unpaired2++;}
					}
				}
				
				if(unpaired1>0 && r.numSites()>0){
					Shared.sort(r.sites);
					Tools.removeLowQualitySitesPaired(r.sites, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
					rescue(r, r2, basesP2, basesM2, Tools.min(MAX_PAIR_DIST, 2*AVERAGE_PAIR_DIST+100));
					Tools.mergeDuplicateSites(r2.sites, true, true);
				}
				
				if(unpaired2>0 && r2.numSites()>0){
					Shared.sort(r2.sites);
					Tools.removeLowQualitySitesPaired(r2.sites, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
					rescue(r2, r, basesP1, basesM1, Tools.min(MAX_PAIR_DIST, 2*AVERAGE_PAIR_DIST+100));
					Tools.mergeDuplicateSites(r.sites, true, true);
				}

				attempt.postRescueSiteSum1=r.numSites();
				attempt.postRescueSiteSum2=r2.numSites();
				
//				if(r.list!=null){Shared.sort(r.list);}
//				if(r2.list!=null){Shared.sort(r2.list);}
//
//				Tools.removeLowQualitySites(r.list, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
//				Tools.removeLowQualitySites(r2.list, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
				
				if(verbose){System.err.println("\nAfter rescue:\nRead1:\t"+r+"\nRead2:\t"+r2);}
				assert(Read.CHECKSITES(r, basesM1, false) && Read.CHECKSITES(r2, basesM2, false));
			}
		}else{
			Tools.mergeDuplicateSites(r.sites, true, false);
			Tools.mergeDuplicateSites(r2.sites, true, false);
			if(verbose){System.err.println("\nAfter merge:\nRead1:\t"+r+"\nRead2:\t"+r2);}
			assert(Read.CHECKSITES(r, basesM1, false) && Read.CHECKSITES(r2, basesM2, false));
		}
		
		if(r.numSites()>1){Shared.sort(r.sites);}
		if(r2.numSites()>1){Shared.sort(r2.sites);}
		assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		
		if(false){//This block is optional, but increases correctness by a tiny bit. (or maybe not!)
			if(SLOW_ALIGN || USE_AFFINE_SCORE){
				Tools.removeLowQualitySitesPaired(r.sites, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED);
				Tools.removeLowQualitySitesPaired(r2.sites, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED);
			}

			pairSiteScoresFinal(r, r2, false, false, MAX_PAIR_DIST, AVERAGE_PAIR_DIST, SAME_STRAND_PAIRS, REQUIRE_CORRECT_STRANDS_PAIRS, MAX_TRIM_SITES_TO_RETAIN);
			
			if(r.numSites()>1){Shared.sort(r.sites);}
			if(r2.numSites()>1){Shared.sort(r2.sites);}
		}
		
		if(SLOW_ALIGN || USE_AFFINE_SCORE){
			Tools.removeLowQualitySitesPaired(r.sites, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED);
			Tools.removeLowQualitySitesPaired(r2.sites, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED);
		}
		
		pairSiteScoresFinal(r, r2, true, true, MAX_PAIR_DIST, AVERAGE_PAIR_DIST, SAME_STRAND_PAIRS, REQUIRE_CORRECT_STRANDS_PAIRS, MAX_TRIM_SITES_TO_RETAIN);
		if(verbose){System.err.println("\nAfter final pairing:\nRead1:\t"+r+"\nRead2:\t"+r2);}
		
		if(r.numSites()>0){
			attempt.mapped1=1;
			Shared.sort(r.sites);
		}
		if(r2.numSites()>0){
			attempt.mapped2=1;
			Shared.sort(r2.sites);
		}
		assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		
	}

	// One reusable record per worker. No per-read allocation; currently one attempt only.
	private final PairAttemptAccounting pairAttempt=new PairAttemptAccounting();
	private final boolean quantumOnly;
	private final boolean quantumTieredShadow;
	private final boolean quantumTieredMutate;
	private final boolean quantumHybrid;
	private final QuantumRanker quantumTieredRanker;
	private final TieredQuantumStats quantumTieredStats;
	private final QuantumOnlyStats quantumOnlyStats;
	private final QuantumHybridStats quantumHybridStats;
	private final SingleAttemptAccounting singleAttempt=new SingleAttemptAccounting();
	private final HybridMaxIndelConfig hybridMaxIndelConfig;
	private final SingleAttemptAccounting singleSecondAttempt;
	private final ReadSearchState singleEntryState;
	private final ReadSearchState singleFirstState;
	private final HybridMaxIndelStats hybridMaxIndelStats;
	private static final int QUANTUM_SCORE_SCALE=100;
	private static final int QUANTUM_SCORE_OFFSET=-30;
	private static final int QUANTUM_EDIT_COST=2*QUANTUM_SCORE_SCALE;
	private static final class SingleAttemptAccounting {
		int quick,max,imperfect;
		int initialSiteSum,postTrimSiteSum,mapped,rejectedBases;
		boolean rejected,committed=true;
		void reset(){
			assert(committed) : "The previous single-read attempt must commit before storage reuse";
			quick=max=imperfect=initialSiteSum=postTrimSiteSum=mapped=rejectedBases=0;
			rejected=false;committed=false;
		}
		void discard(){
			assert(!committed) : "Only an uncommitted single-read attempt may be discarded";
			committed=true;
		}
		void commit(BBMapThread owner){
			assert(!committed) : "Single-read output and counters must commit exactly once";
			owner.readsUsed1++;
			owner.initialSiteSum1+=initialSiteSum;
			owner.postTrimSiteSum1+=postTrimSiteSum;
			owner.mapped1+=mapped;
			if(rejected){
				owner.lowQualityReadsDiscarded1++;
				owner.lowQualityBasesDiscarded1+=rejectedBases;
			}
			committed=true;
		}
	}
	HybridMaxIndelStats hybridMaxIndelStats(){return hybridMaxIndelStats;}
	private static final class PairAttemptAccounting {
		int quick1, quick2, max1, max2, imperfect1, imperfect2;
		int initialSiteSum1, initialSiteSum2, postTrimSiteSum1, postTrimSiteSum2;
		int postRescueSiteSum1, postRescueSiteSum2, mapped1, mapped2;
		int rejectedBases1, rejectedBases2;
		boolean rejected, committed=true;
		void reset(){
			assert(committed) : "The previous attempt must be committed before reusing its accounting storage";
			quick1=quick2=max1=max2=imperfect1=imperfect2=0;
			initialSiteSum1=initialSiteSum2=postTrimSiteSum1=postTrimSiteSum2=0;
			postRescueSiteSum1=postRescueSiteSum2=mapped1=mapped2=0;
			rejectedBases1=rejectedBases2=0;rejected=false;committed=false;
		}
		void discard(){
			assert(!committed) : "Only an uncommitted attempt may be discarded";
			committed=true;
		}
		void commit(BBMapThread owner){
			assert(!committed) : "Paired output and intermediate-site statistics must commit exactly once";
			owner.initialSiteSum1+=initialSiteSum1;owner.initialSiteSum2+=initialSiteSum2;
			owner.postTrimSiteSum1+=postTrimSiteSum1;owner.postTrimSiteSum2+=postTrimSiteSum2;
			owner.postRescueSiteSum1+=postRescueSiteSum1;owner.postRescueSiteSum2+=postRescueSiteSum2;
			owner.mapped1+=mapped1;owner.mapped2+=mapped2;
			if(rejected){
				owner.lowQualityReadsDiscarded1++;owner.lowQualityReadsDiscarded2++;
				owner.lowQualityBasesDiscarded1+=rejectedBases1;owner.lowQualityBasesDiscarded2+=rejectedBases2;
			}
			committed=true;
		}
	}
}
