package align2;

/**
 * Deterministic transform from the raw neural-MAPQ pilot vector to the values
 * presented to a network.  Training and inference must call this same method.
 *
 * Scores are expressed relative to the read's maximum nominal match score,
 * query-consuming counts are divided by read length, and reference-spanning
 * edit lengths use a log ratio so a 40 kbp deletion cannot dominate the input
 * scale.  RegressionTrainer may still standardize these values and fold that
 * standardization into the first layer.
 *
 * @author Collei
 */
public final class NeuralMapqFeatureTransform {

	private NeuralMapqFeatureTransform(){}

	public static void transform(final float[] raw, final float[] out){
		NeuralMapqFeatureSchema.validateVector(raw);
		if(out==null || out.length!=WIDTH){
			throw new IllegalArgumentException("Neural MAPQ transform "+SCHEMA_NAME+
					" requires output width "+WIDTH+"; observed "+
					(out==null ? "null" : out.length));
		}
		final float length=raw[0];
		if(!(length>=1) || Float.isInfinite(length)){
			throw new IllegalArgumentException("Neural MAPQ read length must be finite and positive: "+length);
		}
		final float invLength=1f/length;
		final float invMaximumScore=0.01f*invLength;
		final float invLogLength=(float)(1.0/Math.log1p(length));

		out[0]=(float)Math.log1p(length);
		copy(raw,out,1,4); // boolean mapping flags
		out[5]=raw[5]*invMaximumScore;
		out[6]=raw[6]*0.01f;
		out[7]=raw[7]*invMaximumScore;
		out[8]=raw[8]*invMaximumScore;
		out[9]=raw[9]*invMaximumScore;
		out[10]=raw[10];
		out[11]=(float)Math.log1p(nonnegative(raw[11],11));
		out[12]=(float)Math.log1p(nonnegative(raw[12],12));
		copy(raw,out,13,14); // missing-competitor flags
		out[15]=raw[15]*invMaximumScore;
		out[16]=raw[16]*invMaximumScore;
		out[17]=raw[17]*invMaximumScore;
		out[18]=raw[18]*invMaximumScore;
		copy(raw,out,19,20); // competitor score ratios may legitimately be negative
		out[21]=nonnegative(raw[21],21)*invLength;
		out[22]=raw[22];
		out[23]=nonnegative(raw[23],23)*invLength;
		out[24]=nonnegative(raw[24],24)*invLength;
		out[25]=nonnegative(raw[25],25)*invLength;
		out[26]=nonnegative(raw[26],26)*invLength;
		out[27]=nonnegative(raw[27],27)*invLength;
		out[28]=logReadRatio(raw[28],invLogLength,28);
		out[29]=logReadRatio(raw[29],invLogLength,29);
		out[30]=logReadRatio(raw[30],invLogLength,30);
		out[31]=nonnegative(raw[31],31)*invLength;
		out[32]=nonnegative(raw[32],32)*invLength;
		out[33]=nonnegative(raw[33],33)*invLength;
		out[34]=raw[34]*0.02f;
		out[35]=raw[35]*0.02f;
		out[36]=nonnegative(raw[36],36)*invLength;
		validate(out);
	}

	private static void copy(final float[] raw, final float[] out, final int from, final int to){
		for(int i=from; i<=to; i++){out[i]=raw[i];}
	}

	private static float logReadRatio(final float value, final float invLogLength, final int index){
		return (float)Math.log1p(nonnegative(value,index))*invLogLength;
	}

	private static float nonnegative(final float value, final int index){
		if(value<0){
			throw new IllegalArgumentException("Neural MAPQ raw feature "+
					NeuralMapqFeatureSchema.name(index)+" must be nonnegative: "+value);
		}
		return value;
	}

	public static void validate(final float[] vector){
		if(vector==null || vector.length!=WIDTH){
			throw new IllegalArgumentException("Neural MAPQ transformed vector has wrong width");
		}
		for(int i=0; i<vector.length; i++){
			if(Float.isNaN(vector[i]) || Float.isInfinite(vector[i])){
				throw new IllegalArgumentException("Neural MAPQ transformed feature "+name(i)+
						" is not finite: "+vector[i]);
			}
		}
	}

	public static String name(final int index){
		if(index<0 || index>=WIDTH){throw new IndexOutOfBoundsException("feature "+index);}
		return NAMES[index];
	}

	public static String[] names(){return NAMES.clone();}

	/** Candidate shippable single-end schema after validation-only feature ablation. */
	public static final String SCHEMA_NAME="bbmaps_mapq_single_v1";
	public static final int WIDTH=37;

	private static final String[] NAMES={
		"log1p_read_length","rescued","perfect","semiperfect","ambiguous",
		"map_score_fraction","map_score_per_base_scaled","top_quick_score_fraction",
		"top_slow_score_fraction","top_final_score_fraction","top_score_per_base",
		"log1p_retained_site_count","log1p_equal_top_count","second_missing","third_missing",
		"second_score_fraction","third_score_fraction","top_minus_second_fraction",
		"top_minus_third_fraction","second_to_top_ratio","third_to_top_ratio",
		"top_seed_hits_per_base","primary_identity","substitution_events_per_base",
		"substituted_bases_per_base","insertion_events_per_base","inserted_bases_per_base",
		"deletion_events_per_base","deleted_bases_log_read_ratio","total_edit_bases_log_read_ratio",
		"longest_indel_log_read_ratio","gap_count_per_base","clipped_bases_per_base",
		"alignment_n_bases_per_base","mean_quality_scaled","minimum_quality_scaled",
		"expected_errors_per_base"
	};

	static {if(NAMES.length!=WIDTH){throw new AssertionError("Transform name width differs: "+NAMES.length);}}
}
