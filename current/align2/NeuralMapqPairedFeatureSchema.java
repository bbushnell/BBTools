package align2;

/**
 * Ordered 108-field raw contract for a mapped anchor and its mate.
 * Two 42-field single-end blocks precede 24 pair fields. This historical pilot
 * schema is distinct from the 94 transformed pilot inputs and the accepted
 * 51-input NeuralMapqPairedFeatureTransformV1 network contract.
 * @author Collei
 */
public final class NeuralMapqPairedFeatureSchema{

	private NeuralMapqPairedFeatureSchema(){}

	/** Returns a defensive copy of the names in serialization order. */
	public static String[] names(){return NAMES.clone();}

	/** Checks width and finiteness only; flag and feature ranges belong to transforms. */
	public static void validateVector(final float[] vector){
		if(vector==null || vector.length!=WIDTH){
			throw new IllegalArgumentException("Paired neural MAPQ schema "+SCHEMA_NAME+
					" requires "+WIDTH+" features; observed "+
					(vector==null ? "null" : vector.length));
		}
		for(int i=0; i<vector.length; i++){
			if(!Float.isFinite(vector[i])){
				throw new IllegalArgumentException("Paired neural MAPQ feature "+NAMES[i]+
						" is not finite: "+vector[i]);
			}
		}
	}

	/** Prefixes the shared raw end names without changing their order. */
	private static String[] makeNames(){
		final String[] single=NeuralMapqFeatureSchema.names();
		final String[] names=new String[single.length*2+PAIR_NAMES.length];
		int i=0;
		for(final String name : single){names[i++]="anchor_"+name;}
		for(final String name : single){names[i++]="mate_"+name;}
		for(final String name : PAIR_NAMES){names[i++]=name;}
		return names;
	}

	public static final String SCHEMA_NAME="bbmaps_mapq_paired_pilot_v1";
	public static final int ORACLE_TOLERANCE=20;
	private static final String[] PAIR_NAMES={
		"anchor_pairnum","mate_mapped","anchor_paired","mate_paired",
		"same_chromosome","same_scaffold","same_strand","expected_orientation",
		"anchor_insert_valid","observed_insert_missing","observed_insert","stored_insert",
		"average_inner_distance","expected_fragment_length","inner_distance",
		"signed_inner_deviation","absolute_inner_deviation","anchor_paired_score",
		"mate_paired_score","anchor_pair_score_gain","mate_pair_score_gain",
		"pair_score_sum","combined_read_length","both_perfect"
	};
	private static final String[] NAMES=makeNames();
	public static final int WIDTH=NAMES.length;
	static{if(WIDTH!=108){throw new AssertionError("Paired pilot width differs: "+WIDTH);}}
}
