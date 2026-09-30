package align2;

/**
 * Ordered 42-field raw feature contract for neural-MAPQ extraction and export.
 *
 * The historical pilot name identifies raw rows, not the trained network's
 * input schema. NeuralMapqFeatureTransform selects and scales 37 inputs under
 * its own schema name. Paired extraction embeds two of these raw blocks.
 * Field order is a data-format contract shared by exporters and transforms.
 *
 * @author Collei
 */
public final class NeuralMapqFeatureSchema{

	private NeuralMapqFeatureSchema(){}

	/** Returns a defensive copy so a diagnostic caller cannot mutate the contract. */
	public static String[] names(){return NAMES.clone();}

	/** Returns the raw field name at an index in [0, WIDTH). */
	public static String name(final int index){
		assert(index>=0 && index<NAMES.length) : "Neural MAPQ feature index outside "+
				SCHEMA_NAME+" width "+NAMES.length+": "+index;
		return NAMES[index];
	}

	/** Returns the case-sensitive field index, or -1 for null or unknown names. */
	public static int indexOf(final String name){
		if(name==null){return -1;}
		for(int i=0; i<NAMES.length; i++){
			if(NAMES[i].equals(name)){return i;}
		}
		return -1;
	}

	/** Checks raw width and finiteness; semantic ranges belong to the transform.
	 * Does not validate flags, read length, ordering, or biological consistency. */
	public static void validateVector(final float[] vector){
		if(vector==null || vector.length!=WIDTH){
			throw new IllegalArgumentException("Neural MAPQ schema "+SCHEMA_NAME+
					" requires "+WIDTH+" features; observed "+
					(vector==null ? "null" : vector.length));
		}
		for(int i=0; i<vector.length; i++){
			final float value=vector[i];
			if(Float.isNaN(value) || Float.isInfinite(value)){
				throw new IllegalArgumentException("Neural MAPQ feature "+name(i)+
						" is not finite: "+value);
			}
		}
	}

	public static final String SCHEMA_NAME="bbmaps_mapq_single_pilot_v1";
	/** Loose-oracle endpoint tolerance in reference bases, not an extraction feature. */
	public static final int ORACLE_TOLERANCE=20;
	/** Upper bound for exported/calibrated MAPQ; accepted length caps can be lower. */
	public static final int MAX_MAPQ=50;

	private static final String[] NAMES={
		"read_length",
		"rescued",
		"perfect",
		"semiperfect",
		"ambiguous",
		"map_score",
		"map_score_per_base",
		"top_quick_score",
		"top_slow_score",
		"top_final_score",
		"top_score_per_base",
		"retained_site_count",
		"equal_top_count",
		"second_missing",
		"third_missing",
		"second_score",
		"third_score",
		"top_minus_second",
		"top_minus_third",
		"second_to_top_ratio",
		"third_to_top_ratio",
		"top_seed_hits",
		"primary_identity",
		"substitution_events",
		"substituted_bases",
		"insertion_events",
		"inserted_bases",
		"deletion_events",
		"deleted_bases",
		"total_edit_bases",
		"longest_indel",
		"gap_count",
		"clipped_bases",
		"alignment_n_bases",
		"mean_quality",
		"minimum_quality",
		"expected_errors",
		"read_n_fraction",
		"read_gc_fraction",
		"homopolymer_fraction",
		"monomer_entropy",
		"dimer_entropy"
	};

	public static final int WIDTH=NAMES.length;
}
