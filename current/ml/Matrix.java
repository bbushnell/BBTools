package ml;

import java.util.ArrayList;

import structures.ByteBuilder;

/**
 * Mutable holder for shared input, output, and weight rows loaded for neural-network training.
 * Output range fields cache detected or assigned statistics after in-place transforms.
 */
public class Matrix{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Detects output statistics, then applies optional binarization and affine scaling. */
	public void initializeRange(){
		detectRange();
		if(convertTo01){
			convertToZeroOne(outputMidpoint);
		}
		if(setTargetOutputRangeMin || setTargetOutputRangeMax){
			adjustRange();
		}
	}

	/*--------------------------------------------------------------*/

	/** Computes output statistics over all cells; asserts that the observed minimum is below the maximum. */
	void detectRange(){
		outputMin=Float.MAX_VALUE;
		outputMax=-Float.MAX_VALUE;
		double sum=0;
		long count=0;
		for(float[] line : outputs){
			for(float f : line){
				outputMin=Math.min(f, outputMin);
				outputMax=Math.max(f, outputMax);
				sum+=f;
				count++;
			}
		}
		assert(outputMin<outputMax) : outputMin+", "+outputMax;
		outputMean=(float)(sum/count);
		outputRange=outputMax-outputMin;
		outputMidpoint=outputMin+outputRange*0.5f;
	}

	/** Binarizes every output cell against {@code cutoff} and assigns nominal 0/1 statistics. */
	void convertToZeroOne(final float cutoff){
		double sum=0;
		long count=0;
		//Repeated row aliases would still be revisited here; DataLoader.balance now copies added targets (ML047).
		for(float[] line : outputs){
			for(int j=0; j<line.length; j++){
				final float f=line[j]<cutoff ? 0 : 1;
				line[j]=f;
				sum+=f;
				count++;
			}
		}
		outputMin=0;
		outputMax=1;
		outputMean=(float)(sum/count);
		outputRange=outputMax-outputMin;
		outputMidpoint=outputMin+outputRange*0.5f;
	}

	/**
	 * Affine-translates outputs into the requested target range and assigns that target
	 * as the cached bounds; omitted endpoints are first filled from the current cache.
	 */
	void adjustRange(){
		assert(setTargetOutputRangeMin || setTargetOutputRangeMax) : "Must set minoutput or maxoutput";
		assert(outputMin<outputMax) : outputMin+", "+outputMax;
		if(!setTargetOutputRangeMin){targetOutputRangeMin=outputMin;}
		if(!setTargetOutputRangeMax){targetOutputRangeMax=outputMax;}
		if(targetOutputRangeMin==outputMin && targetOutputRangeMax==outputMax){return;}//Nothing to do

		final float range2=targetOutputRangeMax-targetOutputRangeMin;
		final float mult=range2/outputRange;

		double sum=0;
		long count=0;
		//Transforms each entry in place. DataLoader.balance supplies independent target rows (ML047);
		//other callers that supply repeated row aliases would still transform those arrays repeatedly.
		for(float[] line : outputs){
			for(int i=0; i<line.length; i++){
				float f=((line[i]-outputMin)*mult)+targetOutputRangeMin;
				line[i]=f;
				sum+=f;
				count++;
			}
		}
		outputMin=targetOutputRangeMin;
		outputMax=targetOutputRangeMax;
		outputMean=(float)(sum/count);
		outputRange=outputMax-outputMin;
		outputMidpoint=outputMin+outputRange*0.5f;
	}

	/*--------------------------------------------------------------*/

	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		bb.append(columns.toString()).nl();
		bb.append("lines="+inputs.length).nl();
		bb.append("inputs="+numInputs).nl();
		bb.append("outputs="+numOutputs).nl();
		bb.append("mean="+outputMean).nl();
		bb.append("midpoint="+outputMidpoint).nl();
		bb.append("range="+outputRange).nl();
		bb.append("inputs="+numInputs).nl();
		return bb.toString();
	}

	int numInputs(){return numInputs;}

	int numOutputs(){return numOutputs;}

	public float outputMidpoint(){return outputMidpoint;}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Optional column names parsed from a ## header; shared by both exclusive and inclusive split views. */
	public ArrayList<String> columns;

	/** Header dimensions: input count, output count, and optional weight flag. */
	int[] dims;

	public int numInputs;

	public int numOutputs;

	int numPositive=0;

	int numNegative=0;

	public int validLines=0;

	int invalidLines=0;

	/** Cached lower output bound after the last range transform. */
	private float outputMin;

	/** Cached upper output bound after the last range transform. */
	private float outputMax;

	/** Cached mean across output cells after the last transform. */
	private float outputMean;

	/** Cached midpoint between outputMin and outputMax. */
	private float outputMidpoint;

	/** Cached outputMax-outputMin span. */
	private float outputRange;

	/** Shared row groups in input, output, weight order. */
	public float[][][] data;

	public float[][] inputs;

	public float[][] outputs;

	public float[][] weights;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether outputs should be binarized around their observed midpoint. */
	public static boolean convertTo01=false;

	/** Lower endpoint for optional affine output scaling; may be filled from observed outputMin. */
	public static float targetOutputRangeMin=0;

	/** Upper endpoint for optional affine output scaling; may be filled from observed outputMax. */
	public static float targetOutputRangeMax=0;

	public static boolean setTargetOutputRangeMin=false;

	public static boolean setTargetOutputRangeMax=false;
}
