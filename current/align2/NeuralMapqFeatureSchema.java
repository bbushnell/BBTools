package align2;

/**
 * Ordered single-end feature contract for the neural-MAPQ exporter pilot.
 *
 * This pilot name must not be bundled with a trained model.  Freeze a new
 * non-pilot schema only after observed ranges and ablations establish the final
 * field set.
 *
 * @author Collei
 */
public final class NeuralMapqFeatureSchema {

	private NeuralMapqFeatureSchema(){}

	/** Returns a defensive copy so a diagnostic caller cannot mutate the contract. */
	public static String[] names(){return NAMES.clone();}

	public static String name(final int index){
		assert(index>=0 && index<NAMES.length) : "Neural MAPQ feature index outside "+
				SCHEMA_NAME+" width "+NAMES.length+": "+index;
		return NAMES[index];
	}

	public static int indexOf(final String name){
		if(name==null){return -1;}
		for(int i=0; i<NAMES.length; i++){
			if(NAMES[i].equals(name)){return i;}
		}
		return -1;
	}

	/** Fails loudly before a vector can be written or passed to a network. */
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
	public static final int ORACLE_TOLERANCE=20;
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
