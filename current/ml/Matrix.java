package ml;

import java.util.ArrayList;
import structures.ByteBuilder;

/**
 * Represents a data matrix for machine learning operations.
 * Handles input/output data transformation, range detection, and format conversion.
 * Supports automatic range adjustment and binary classification conversion.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 * @version 1.0
 */
public class Matrix {
    
    /*--------------------------------------------------------------*/
    /*----------------           Methods            ----------------*/
    /*--------------------------------------------------------------*/
    
    /**
     * Initializes the matrix by detecting data ranges and applying transformations.
     * Performs range detection, optional binary conversion, and range adjustment.
     */
    void initializeRange() {
        detectRange();
        if(convertTo01) {
            convertToZeroOne(outputMidpoint);
        }
        if(setTargetOutputRangeMin || setTargetOutputRangeMax) {
            adjustRange();
        }
    }

    /*--------------------------------------------------------------*/
    
    /**
     * Analyzes output data to determine min, max, mean, and midpoint values.
     * Calculates statistical properties of the output data for normalization
     * and classification threshold determination.
     */
    void detectRange() {
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
    
    /**
     * Converts continuous output values to binary (0 or 1) classification.
     * Values below the cutoff become 0, values at or above become 1.
     * Updates all statistical measures to reflect the binary conversion.
     * 
     * @param cutoff Threshold value for binary conversion
     */
    void convertToZeroOne(final float cutoff) {
        double sum=0;
        long count=0;
        for(float[] line : outputs){
            for(int j=0; j<line.length; j++){
                final float f=line[j]<cutoff ? 0 : 1;
                line[j]=f;
                sum+=f;
            }
        }
        outputMin=0;
        outputMax=1;
        outputMean=(float)(sum/count);
        outputRange=outputMax-outputMin;
        outputMidpoint=outputMin+outputRange*0.5f;
    }
    
    /**
     * Adjusts output values to fit within specified target range.
     * Linearly scales all output values to match the desired min/max bounds.
     * Recalculates all statistical measures after adjustment.
     */
    void adjustRange() {
        assert(setTargetOutputRangeMin || setTargetOutputRangeMax) : "Must set minoutput or maxoutput";
        assert(outputMin<outputMax) : outputMin+", "+outputMax;
        if(!setTargetOutputRangeMin){targetOutputRangeMin=outputMin;}
        if(!setTargetOutputRangeMax){targetOutputRangeMax=outputMax;}
        if(targetOutputRangeMin==outputMin && targetOutputRangeMax==outputMax) {return;}//Nothing to do
        
        final float range2=targetOutputRangeMax-targetOutputRangeMin;
        assert(range2!=outputRange);
        final float mult=range2/outputRange;
        
        double sum=0;
        long count=0;
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
    
    /**
     * Generates a string representation of the matrix properties.
     * Includes column headers, dimensions, and statistical measures.
     * 
     * @return Formatted string describing the matrix
     */
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

    /**
     * Returns the number of input features.
     * 
     * @return Number of input dimensions
     */
    int numInputs() {return numInputs;}
    
    /**
     * Returns the number of output dimensions.
     * 
     * @return Number of output values
     */
    int numOutputs() {return numOutputs;}
    
    /**
     * Returns the midpoint value of the output range.
     * Used as the default threshold for binary classification.
     * 
     * @return Output range midpoint
     */
    public float outputMidpoint() {return outputMidpoint;}
    
    /*--------------------------------------------------------------*/
    /*----------------           Fields             ----------------*/
    /*--------------------------------------------------------------*/
    
    /** Column headers for data interpretation */
    ArrayList<String> columns;
    
    /** Dimensional structure of the data */
    int[] dims;
    
    /** Number of input features per sample */
    int numInputs;
    
    /** Number of output values per sample */
    int numOutputs;
    
    /** Count of positive examples in the dataset */
    int numPositive=0;
    
    /** Count of negative examples in the dataset */
    int numNegative=0;
    
    /** Number of successfully parsed data lines */
    int validLines=0;
    
    /** Number of unparseable or invalid data lines */
    int invalidLines=0;
    
    /** Minimum value in the output data */
    private float outputMin;
    
    /** Maximum value in the output data */
    private float outputMax;
    
    /** Mean of all output values */
    private float outputMean;
    
    /** Midpoint of the output range (used for classification threshold) */
    private float outputMidpoint;
    
    /** Range of output values (max - min) */
    private float outputRange;
    
    /** Complete data structure: [inputs, outputs, weights] */
    float[][][] data;
    
    /** Input feature vectors for each sample */
    float[][] inputs;
    
    /** Output target values for each sample */
    float[][] outputs;
    
    /** Sample weights for training emphasis */
    float[][] weights;
    
    /*--------------------------------------------------------------*/
    /*----------------        Static Fields         ----------------*/
    /*--------------------------------------------------------------*/
    
    /** Flag to enable binary conversion of outputs */
    static boolean convertTo01=false;
    
    /** Target minimum value for range adjustment */
    static float targetOutputRangeMin=0;
    
    /** Target maximum value for range adjustment */
    static float targetOutputRangeMax=0;
    
    /** Flag indicating whether to set custom minimum range */
    static boolean setTargetOutputRangeMin=false;
    
    /** Flag indicating whether to set custom maximum range */
    static boolean setTargetOutputRangeMax=false;
}