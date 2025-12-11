package ml;

import java.util.ArrayList;

import structures.ByteBuilder;

/** No JavaDoc generation required - class already fully documented */
public class Matrix {
    
    /*--------------------------------------------------------------*/
    /*----------------           Methods            ----------------*/
    /*--------------------------------------------------------------*/
    
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

    int numInputs() {return numInputs;}
    
    int numOutputs() {return numOutputs;}
    
    public float outputMidpoint() {return outputMidpoint;}
    
    /*--------------------------------------------------------------*/
    /*----------------           Fields             ----------------*/
    /*--------------------------------------------------------------*/
    
    ArrayList<String> columns;
    
    int[] dims;
    
    int numInputs;
    
    int numOutputs;
    
    int numPositive=0;
    
    int numNegative=0;
    
    int validLines=0;
    
    int invalidLines=0;
    
    private float outputMin;
    
    private float outputMax;
    
    private float outputMean;
    
    private float outputMidpoint;
    
    private float outputRange;
    
    float[][][] data;
    
    float[][] inputs;
    
    float[][] outputs;
    
    float[][] weights;
    
    /*--------------------------------------------------------------*/
    /*----------------        Static Fields         ----------------*/
    /*--------------------------------------------------------------*/
    
    static boolean convertTo01=false;
    
    static float targetOutputRangeMin=0;
    
    static float targetOutputRangeMax=0;
    
    static boolean setTargetOutputRangeMin=false;
    
    static boolean setTargetOutputRangeMax=false;
}