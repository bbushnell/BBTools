package ml;

import java.util.Arrays;

import structures.ByteBuilder;

/** Stores one input/goal row, its current network result, and cached error
 * statistics used for training prioritization and confusion-matrix counts.
 *
 * @author Brian Bushnell
 * @date June 3, 2025 */
public class Sample implements Comparable<Sample>{

	/** Constructs a Sample with input features, target outputs, weight, and ID.
	 * Initializes the result array and caches positivity from goal[0]>=0.5.
	 *
	 * @param in_ input feature vector
	 * @param out_ target output values
	 * @param weight_ sample-specific multiplier for training, must be >0
	 * @param id_ unique identifier for this sample */
	public Sample(float[] in_, float[] out_, float weight_, int id_){
		this(in_, out_, weight_, id_, 0.5f);
	}

	/** Constructs a Sample using the supplied positivity cutoff. */
	public Sample(float[] in_, float[] out_, float weight_, int id_, float positiveCutoff){
		in=in_;
		goal=out_;
		result=new float[goal.length];
		weight=weight_;
		id=id_;
		positive=(goal[0]>=positiveCutoff);
		assert(weight>0) : weight;
	}

	@Override
	public int compareTo(Sample o){
		final float a=pivot, b=o.pivot;
		return a>b ? -1 : b>a ? 1 : id-o.id;
	}

	/**
	 * Validates that the cached pivot value matches the calculated pivot.
	 * Used for debugging pivot calculation consistency.
	 * @return true if cached pivot equals calculated pivot
	 */
	public boolean checkPivot(){
		return pivot==calcPivot();
	}

	/** Recalculates the cached pivot while holding this sample's monitor. */
	synchronized void setPivot(){
		pivot=calcPivot();
	}

	/**
	 * Calculates the pivot value used for sample prioritization during training.
	 * Combines error magnitude with excess error penalty and epoch decay.
	 * Positive samples that are predicted too high get excess penalty.
	 * @return Calculated pivot value for sample sorting
	 */
	synchronized float calcPivot(){
		final float v=result[0];
		final boolean positiveError=v>goal[0];
		final boolean excess=(positiveError==positive);
		final float mult=(excess ? excessPivotMult*0.5f : 0.5f);
		return (errorMagnitude+weightedErrorMagnitude)*mult-epoch*EPOCH_MULT;
//		return (errorMagnitude+weightedErrorMagnitude)*0.5f-epoch*EPOCH_MULT;
	}

	/**
	 * Creates a detailed string representation showing sample classification metrics.
	 * Includes sample ID, positive/negative indicator, confusion matrix type
	 * (TP/TN/FP/FN), epoch, goal, result, error magnitude, pivot, and input features.
	 * @return Formatted string with comprehensive sample information
	 */
	public String toString(){
//		String s="S%d\t%s\t%s\tep=%d\tg=%4f\tr=%4f\tem=%6f\tev=%.6f\tpv=%.6f";
		String s="S%d\t%s\t%s\tep=%d\tg=%4f\tr=%4f\tem=%6f\tpv=%.6f";

		boolean gol=(goal[0]>=0.5f);
		boolean pred=(result[0]>=0.5f);
		String type=(gol && pred) ? "TP" : (!gol && !pred) ? "TN" :
			(!gol && pred) ? "FP" : (gol && !pred) ? "FN" : "??";
		String sign=(positive ? "+" : "-");

//		s=String.format(s, id, sign, type, epoch, goal[0], result[0], errorMagnitude, errorValue, calcPivot());
		s=String.format(s, id, sign, type, epoch, goal[0], result[0], errorMagnitude, calcPivot());
		return s+"\t"+Arrays.toString(in);
	}

	/** Converts sample to tab-separated byte representation.
	 * @return ByteBuilder containing input and goal values */
	public ByteBuilder toBytes(){
		return toBytes(new ByteBuilder());
	}

	/**
	 * Appends sample data to provided ByteBuilder as tab-separated values.
	 * Requests six fractional digits for input features and goals; exact integers
	 * are emitted without a fractional part by ByteBuilder.append(double, int).
	 * @param bb ByteBuilder to append data to
	 * @return The provided ByteBuilder with sample data appended
	 */
	public ByteBuilder toBytes(ByteBuilder bb){
		for(float f : in){bb.append(f, 6).tab();}
		for(float f : goal){bb.append(f, 6).tab();}
		bb.trimLast(1);
		bb.nl();
		return bb;
	}

//	synchronized boolean positive() {
//		return goal[0]>=0.5f;
//	}

	/**
	 * Calculates error metrics for this sample using current prediction results.
	 * Stores raw half-squared-error magnitude and class-weighted magnitude.
	 * @param weightMult multiplier for class-weighted error calculation */
	public void calcError(float weightMult){
		double error=0;
		for(int i=0; i<result.length; i++){
			float r=result[i];
			float g=goal[i];
			float e=calcError(g, r);
			assert(e>=0);
			error+=e;
		}
		errorMagnitude=(float)error;
		assert(error>=0);
		weightedErrorMagnitude=Cell.toWeightedError(error, result[0], goal[0], weightMult);
		assert(weightedErrorMagnitude>=0);
	}

	/** Gets the current training epoch for this sample */
	public synchronized int epoch(){return epoch;}
	/** Gets the last thread ID that processed this sample */
	public synchronized int lastTID(){return lastTID;}
	/** Sets the training epoch for this sample.
	 * @param x Epoch number (cast to int) */
	public synchronized void setEpoch(long x){
		epoch=(int)x;
	}

	/** Sets the last thread ID that processed this sample.
	 * @param x Thread ID */
	public synchronized void setLastTID(int x){
		lastTID=x;
	}

	/** Calculates half squared error between goal and prediction values.
	 *
	 * @param goal target value
	 * @param pred predicted value
	 * @return 0.5*(goal-pred)^2 */
	public static final float calcError(float goal, float pred){
		float e=goal-pred;
		return 0.5f*e*e;
	}

	/** Whether this sample met the positivity cutoff used for materialization */
	public final boolean positive;
	/**
	 * Raw error magnitude calculated from squared differences between goal and result
	 */
	public float errorMagnitude=1;
	/** Class-weighted error magnitude using Cell.toWeightedError calculation */
	public float weightedErrorMagnitude=1;
//	float errorValue=1;//Unused, commented for efficiency
	/** Training epoch when this sample was last processed */
	private int epoch=-1;
	/** Thread ID that last processed this sample */
	private int lastTID=-1;
	/** Cached pivot value used for sample prioritization during training */
	float pivot=0;

	/** Input feature vector for this training sample */
	final float[] in;
	/** Target output values for this training sample */
	public final float[] goal;
	/** Neural network prediction results for this sample */
	public final float[] result;//Can't be volatile
	/** Sample-specific multiplier applied to CellNet.processSample's weightMult before
	 * weighted error and backpropagation; an incoming zero remains zero. */
	final float weight;//TODO
	/** Unique identifier for this sample */
	final int id;

	//0.2f is good for binary classifiers
	/**
	 * Multiplier for excess error penalty in pivot calculation (0.2f works well for binary classifiers)
	 */
	public static float excessPivotMult=0.2f;
	/** Epoch decay multiplier for pivot calculation (1/256f) */
	public static final float EPOCH_MULT=1/256f;
}
