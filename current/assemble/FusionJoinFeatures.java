package assemble;

import java.util.ArrayList;

/**
 * Versioned, reference-free features for a proposed exact-overlap fusion.
 * Input regions are sorted read-kmer depths: left, overlap, right, boundary,
 * and the extra subset spanning both boundaries. No truth or identity is accepted.
 * Flank minima/maxima are pointwise statistics, not a designated physical flank.
 * One instance supplies reusable scratch and output buffers for one worker;
 * concurrent calls on that instance are unsupported. Input arrays are borrowed
 * for the call and neither retained nor modified. The 107 columns are ordered
 * by {@link #names()}; changing their meaning or order requires a new schema.
 * @author Fischl
 */
public final class FusionJoinFeatures{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Fills a borrowed vector, overwritten by the next call. Coordinates and depths
	 * refer to the retained joined window; trims describe discarded terminal bases.
	 * Deciles use floor(q*(n-1)). Missing regions have zero values and a presence bit.
	 * Flank minima/maxima are taken separately for each summary column, so those
	 * columns need not all describe the same flank. Enrichment needs both flanks
	 * and overlap observations; its columns are zero when any of these is absent.
	 *
	 * @param depths Five arrays with sorted, nonnegative populated prefixes:
	 * left flank, overlap, right flank, crossing either boundary, crossing both
	 * @param sizes Populated lengths; the first four sum to words, and the fifth
	 * is a subset of the fourth rather than another part of that sum
	 * @param k Positive length of the counted words in bases
	 * @param overlap Positive length of the shared overlap in bases
	 * @param leftBases Retained window bases before the overlap
	 * @param rightBases Retained window bases after the overlap
	 * @param sourceTrim Discarded terminal bases from the source contig
	 * @param destTrim Discarded terminal bases from the destination contig
	 * @param words Positive number of word positions in the retained window
	 * @param missing Number of those positions with zero observed depth
	 * @param tested Number of directed continuation comparisons near the overlap
	 * @param half Comparisons with a positive alternative at least half the intended depth
	 * @param dominated Comparisons with a positive alternative at least the intended depth;
	 * no greater than half, which is no greater than tested
	 * @param maxAltRatio Largest alternative depth divided by max(1, intended depth),
	 * finite and nonnegative; the producer supplies zero when no comparison exists
	 * @return This instance's reusable 107-element vector; copy it to retain a result
	 */
	public float[] fill(final int[][] depths, final int[] sizes, final int k, final int overlap,
			final int leftBases, final int rightBases, final int sourceTrim, final int destTrim,
			final int words, final int missing, final int tested, final int half,
			final int dominated, final double maxAltRatio){
		assert(validRegions(depths, sizes)) : "Fusion quantiles require five sorted, nonnegative populated regions.";
		assert(k>0 && overlap>0 && leftBases>=0 && rightBases>=0 && sourceTrim>=0 && destTrim>=0) :
				"Fusion geometry must describe retained nonnegative flanks and discarded trims.";
		assert(words>0 && missing>=0 && missing<=words && tested>=0 && dominated>=0 &&
				dominated<=half && half<=tested && maxAltRatio>=0 && Double.isFinite(maxAltRatio)) :
				"Fusion fractions require real word/branch denominators and nested competitor counts.";
		assert(sizes[0]+sizes[1]+sizes[2]+sizes[3]==words && sizes[4]<=sizes[3]) :
				"Pure regions and boundary words partition the window; span is a boundary subset.";
		int p=0;
		vector[p++]=logDepth(k);
		vector[p++]=logDepth(overlap);
		vector[p++]=logOnePlus(overlap/(double)k);
		vector[p++]=logOnePlus(Math.min(leftBases, rightBases));
		vector[p++]=logOnePlus(Math.max(leftBases, rightBases));
		vector[p++]=logOnePlus(Math.min(sourceTrim, destTrim));
		vector[p++]=logOnePlus(Math.max(sourceTrim, destTrim));
		vector[p++]=logOnePlus(words);
		vector[p++]=missing/(float)words;
		vector[p++]=tested>0 ? 1 : 0;
		vector[p++]=logOnePlus(tested);
		vector[p++]=logOnePlus(maxAltRatio);
		vector[p++]=tested>0 ? half/(float)tested : 0;
		vector[p++]=tested>0 ? dominated/(float)tested : 0;
		summarize(depths[0], sizes[0], left);
		summarize(depths[2], sizes[2], right);
		for(int i=0; i<REGION_WIDTH; i++){vector[p++]=Math.min(left[i], right[i]);}
		for(int i=0; i<REGION_WIDTH; i++){vector[p++]=Math.max(left[i], right[i]);}
		for(int region : CENTRAL_REGIONS){
			summarize(depths[region], sizes[region], center);
			System.arraycopy(center, 0, vector, p, REGION_WIDTH);
			p+=REGION_WIDTH;
		}
		// Differences of scaled log-depths expose enrichment without unbounded ratios.
		final boolean available=sizes[0]>0 && sizes[1]>0 && sizes[2]>0;
		vector[p++]=available ? 1 : 0;
		for(int q=0; q<=10; q++){
			final float anchor=sizes[1]>0 ? logDepth(quantile(depths[1], sizes[1], q)) : 0;
			vector[p++]=available ? anchor-Math.min(left[3+q], right[3+q]) : 0;
			vector[p++]=available ? anchor-Math.max(left[3+q], right[3+q]) : 0;
		}
		assert(p==NAMES.length) : "Feature values and schema names must have identical width: "+p+" vs "+NAMES.length;
		for(float value : vector){assert(Float.isFinite(value)) : "Model inputs may not contain NaN or infinity.";}
		return vector;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Fills all fourteen output columns from a sorted prefix without changing it:
	 * presence, log position count, zero fraction, then eleven log depth deciles.
	 * An empty prefix produces fourteen zeros. */
	private static void summarize(final int[] values, final int n, final float[] result){
		assert(n>=0 && n<=values.length && result.length==REGION_WIDTH) : "Region summary buffers must match the schema.";
		result[0]=n>0 ? 1 : 0;
		result[1]=logOnePlus(n);
		int zeros=0;
		while(zeros<n && values[zeros]==0){zeros++;}
		result[2]=n>0 ? zeros/(float)n : 0;
		for(int q=0; q<=10; q++){result[3+q]=n>0 ? logDepth(quantile(values, n, q)) : 0;}
	}

	/** Returns the observed depth at floor(decile*(n-1)/10) in a nonempty sorted
	 * prefix. Endpoints select its minimum/maximum; there is no interpolation. */
	private static int quantile(final int[] values, final int n, final int decile){
		assert(n>0 && n<=values.length && decile>=0 && decile<=10) : "A requested decile needs at least one observed word.";
		return values[(int)((long)decile*(n-1)/10)];
	}

	/** Returns one eighth of log2(max(1, depth)), matching the PacBio depth transform.
	 * Zero and one both encode as zero; presence and zero-fraction columns distinguish them. */
	private static float logDepth(final int depth){return (float)(LOG_SCALE*Math.log(Math.max(1, depth)));}

	/** Returns one eighth of log2(1+value) for nonnegative geometry/counts/ratios.
	 * True zero remains zero and large values are not clipped. */
	private static float logOnePlus(final double value){return (float)(LOG_SCALE*Math.log1p(value));}

	/** Checks the five arrays and their nonnegative sorted prefixes. Unused suffixes
	 * are ignored; geometry and cross-region count constraints are checked by fill. */
	private static boolean validRegions(final int[][] depths, final int[] sizes){
		if(depths==null || sizes==null || depths.length!=5 || sizes.length!=5){return false;}
		for(int r=0; r<5; r++){
			if(depths[r]==null || sizes[r]<0 || sizes[r]>depths[r].length){return false;}
			for(int i=0; i<sizes[r]; i++){
				if(depths[r][i]<0 || (i>0 && depths[r][i]<depths[r][i-1])){return false;}
			}
		}
		return true;
	}

	/** Returns a fresh copy of all 107 column names in vector order, protecting the
	 * shared schema from caller mutations. */
	public static String[] names(){return NAMES.clone();}

	/** Builds named columns once; the schema version changes if semantics/order change. */
	private static String[] makeNames(){
		final ArrayList<String> names=new ArrayList<String>();
		for(String name : new String[]{"k_log", "overlap_log", "overlap_k_log1p", "flank_bases_min_log1p",
				"flank_bases_max_log1p", "trim_min_log1p", "trim_max_log1p", "words_log1p", "missing_fraction",
				"branch_present", "branch_tests_log1p", "max_alt_ratio_log1p", "alt_half_fraction", "alt_expected_fraction"}){
			names.add(name);
		}
		for(String region : new String[]{"flank_min", "flank_max", "overlap", "boundary", "span"}){
			names.add(region+"_present");
			names.add(region+"_positions_log1p");
			names.add(region+"_zero_fraction");
			for(int q=0; q<=100; q+=10){names.add(region+"_p"+q+"_log");}
		}
		names.add("enrichment_present");
		for(int q=0; q<=100; q+=10){
			names.add("overlap_minus_flank_min_p"+q);
			names.add("overlap_minus_flank_max_p"+q);
		}
		assert(names.size()==107) : "Schema v1 is a fixed 107-feature contract.";
		return names.toArray(new String[0]);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Independent flank summaries, combined by columnwise minimum and maximum. */
	private final float[] left=new float[REGION_WIDTH], right=new float[REGION_WIDTH];
	/** Scratch reused in turn for overlap, boundary and spanning-word summaries. */
	private final float[] center=new float[REGION_WIDTH];
	/** Borrowed result storage; each fill overwrites every column. */
	private final float[] vector=new float[NAMES.length];

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Schema identifier accompanying sidecar columns and trained model inputs. */
	public static final String VERSION="fusion_join_v1";
	/** Presence, position count, zero fraction and eleven deciles per region. */
	private static final int REGION_WIDTH=14;
	/** Converts natural logarithms into one-eighth-scaled base-two logarithms. */
	private static final double LOG_SCALE=.125/Math.log(2);
	/** Producer indices for overlap, either-boundary and both-boundaries regions. */
	private static final int[] CENTRAL_REGIONS={1, 3, 4};
	/** Immutable-by-convention schema order; public access returns a clone. */
	private static final String[] NAMES=makeNames();
}
