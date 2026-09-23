package align2;

/** Frozen51 paired MAPQ transform: anchor37 plus selected pair state and score inputs. */
public final class NeuralMapqPairedFeatureTransformV1 {
	private NeuralMapqPairedFeatureTransformV1(){}

	public static void transform(final float[] raw,final float[] out,final Scratch scratch){
		NeuralMapqPairedFeatureSchema.validateVector(raw);
		if(out==null||out.length!=WIDTH||scratch==null){throw new IllegalArgumentException("Paired V1 output/scratch differs from "+SCHEMA_NAME);}
		System.arraycopy(raw,0,scratch.anchorRaw,0,RAW_END_WIDTH);NeuralMapqFeatureTransform.transform(scratch.anchorRaw,scratch.anchorTransformed);System.arraycopy(scratch.anchorTransformed,0,out,0,END_WIDTH);
		final float anchorLength=positive(scratch.anchorRaw[0],"anchor_read_length");
		final float mateLength=(raw[PAIR_RAW+1]>0?positive(raw[RAW_END_WIDTH],"mate_read_length"):1);
		final float combinedLength=positive(raw[PAIR_RAW+22],"combined_read_length");
		final float invAnchorScore=0.01f/anchorLength,invMateScore=0.01f/mateLength,invPairScore=0.01f/combinedLength;
		int i=END_WIDTH;
		out[i++]=binary(raw[PAIR_RAW+1],"mate_mapped");out[i++]=binary(raw[PAIR_RAW+2],"anchor_paired");out[i++]=binary(raw[PAIR_RAW+3],"mate_paired");
		out[i++]=binary(raw[PAIR_RAW+4],"same_chromosome");out[i++]=binary(raw[PAIR_RAW+5],"same_scaffold");out[i++]=binary(raw[PAIR_RAW+6],"same_strand");out[i++]=binary(raw[PAIR_RAW+7],"expected_orientation");
		out[i++]=nonnegative(raw[PAIR_RAW+17],"anchor_paired_score")*invAnchorScore;out[i++]=nonnegative(raw[PAIR_RAW+18],"mate_paired_score")*invMateScore;
		out[i++]=nonnegative(raw[PAIR_RAW+19],"anchor_pair_score_gain")*invAnchorScore;out[i++]=nonnegative(raw[PAIR_RAW+20],"mate_pair_score_gain")*invMateScore;
		out[i++]=finite(raw[PAIR_RAW+21],"pair_score_sum")*invPairScore;out[i++]=(float)Math.log1p(combinedLength);out[i++]=binary(raw[PAIR_RAW+23],"both_perfect");
		assert(i==WIDTH) : "Frozen37+14 feature contract must fill width51; wrote "+i;validate(out);
	}

	public static void project(final float[] pilot,final float[] out){
		if(pilot==null||pilot.length!=PILOT_WIDTH||out==null||out.length!=WIDTH){throw new IllegalArgumentException("Paired V1 projection width differs");}
		for(int i=0;i<pilot.length;i++)if(!Float.isFinite(pilot[i])){throw new IllegalArgumentException("Nonfinite pilot feature "+i+": "+pilot[i]);}
		System.arraycopy(pilot,0,out,0,END_WIDTH);int next=END_WIDTH;for(int i=74;i<=80;i++){out[next++]=pilot[i];}for(int i=87;i<=93;i++){out[next++]=pilot[i];}
		assert(next==WIDTH) : "Projection retains indices0-36,74-80,87-93; wrote "+next;validate(out);
	}
	public static String[] names(){return NAMES.clone();}
	public static void validate(final float[] vector){if(vector==null||vector.length!=WIDTH){throw new IllegalArgumentException("Paired V1 width differs");}for(int i=0;i<vector.length;i++)if(!Float.isFinite(vector[i])){throw new IllegalArgumentException("Nonfinite paired V1 feature "+NAMES[i]+": "+vector[i]);}}
	private static float binary(final float value,final String name){if(value!=0&&value!=1){throw new IllegalArgumentException(name+" must be binary: "+value);}return value;}
	private static float positive(final float value,final String name){if(!(value>0)||!Float.isFinite(value)){throw new IllegalArgumentException(name+" must be positive: "+value);}return value;}
	private static float nonnegative(final float value,final String name){if(value<0||!Float.isFinite(value)){throw new IllegalArgumentException(name+" must be nonnegative: "+value);}return value;}
	private static float finite(final float value,final String name){if(!Float.isFinite(value)){throw new IllegalArgumentException(name+" must be finite: "+value);}return value;}
	private static String[] makeNames(){final String[] anchor=NeuralMapqFeatureTransform.names(),names=new String[WIDTH];int i=0;for(final String name:anchor){names[i++]="anchor_"+name;}for(final String name:PAIR_NAMES){names[i++]=name;}return names;}

	public static final class Scratch {final float[] anchorRaw=new float[RAW_END_WIDTH],anchorTransformed=new float[END_WIDTH];}
	public static final String SCHEMA_NAME="bbmaps_mapq_paired_v1";
	public static final int WIDTH=51,PILOT_WIDTH=94;
	private static final int RAW_END_WIDTH=NeuralMapqFeatureSchema.WIDTH,END_WIDTH=NeuralMapqFeatureTransform.WIDTH,PAIR_RAW=RAW_END_WIDTH*2;
	private static final String[] PAIR_NAMES={"mate_mapped","anchor_paired","mate_paired","same_chromosome","same_scaffold","same_strand","expected_orientation","anchor_paired_score_fraction","mate_paired_score_fraction","anchor_pair_score_gain_fraction","mate_pair_score_gain_fraction","pair_score_sum_fraction","log1p_combined_read_length","both_perfect"};
	private static final String[] NAMES=makeNames();
	static {if(END_WIDTH!=37||WIDTH!=END_WIDTH+PAIR_NAMES.length||NAMES.length!=WIDTH){throw new AssertionError("Paired V1 schema width differs");}}
}
