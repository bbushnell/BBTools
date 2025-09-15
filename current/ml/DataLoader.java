package ml;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.LineParser1;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.IntList;

/**
 * Loads and preprocesses machine learning datasets from files.
 * Handles data parsing, shuffling, class balancing, and train/test splitting.
 * Supports multiple file formats and automatic data type inference.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 * @version 1.0
 */
public class DataLoader {
    
    /**
     * Constructs a DataLoader for the specified file(s).
     * Accepts either a single file path or comma-separated multiple files.
     * 
     * @param fname_ File path or comma-separated list of file paths
     */
    private DataLoader(String fname_){
        fnames=new File(fname_).exists() ? 
                new String[] {fname_} : Tools.commaPattern.split(fname_);
    }

    /**
     * Splits a matrix into training and validation sets.
     * Distributes samples between sets while maintaining class balance.
     * 
     * @param m Source matrix to split
     * @param fraction Fraction of samples to assign to validation set
     * @param maxLines1 Maximum lines allowed in validation set
     * @param exclusive If true, creates independent sets; if false, training set contains all data
     * @return Array of two matrices: [training, validation]
     */
    public static Matrix[] split(Matrix m, float fraction, int maxLines1, boolean exclusive) {
        if(m.inputs.length == 0) {
            throw new RuntimeException("Cannot split empty matrix");
        }
        
        if(!exclusive) {fraction=1;}//For convenience...
        
        Matrix[] array=new Matrix[] {new Matrix(), new Matrix()};
        @SuppressWarnings("unchecked")
        ArrayList<float[]>[] inputs=new ArrayList[] {new ArrayList<float[]>(), new ArrayList<float[]>()};
        @SuppressWarnings("unchecked")
        ArrayList<float[]>[] outputs=new ArrayList[] {new ArrayList<float[]>(), new ArrayList<float[]>()};
        @SuppressWarnings("unchecked")
        ArrayList<float[]>[] weights=new ArrayList[] {new ArrayList<float[]>(), new ArrayList<float[]>()};
        
        for(Matrix n : array) {
            n.dims=m.dims;
            n.numInputs=m.numInputs;
            n.numOutputs=m.numOutputs;
            n.columns=m.columns;
        }
        
        Random randy=new Random(0);
        for(int i=0; i<m.inputs.length; i++) {
            float[] in=m.inputs[i], out=m.outputs[i], wt=m.weights[i];
            int positive=(out[0]>=0.5f ? 1 : 0);
            int negative=positive^1;
            int pick=(randy.nextFloat()>=fraction || inputs[1].size()>=maxLines1) ? 0 : 1;
            inputs[pick].add(in);
            outputs[pick].add(out);
            weights[pick].add(wt);
            array[pick].numNegative+=negative;
            array[pick].numPositive+=positive;
            array[pick].validLines++;
        }
        
        // Check for empty datasets after splitting
        if(inputs[0].isEmpty() || inputs[1].isEmpty()) {
            throw new RuntimeException("Split would create empty dataset - try different fraction or more data");
        }
        
        for(int set=0; set<array.length; set++){
            Matrix n=array[set];
            n.inputs=new float[inputs[set].size()][];
            n.outputs=new float[outputs[set].size()][];
            n.weights=new float[weights[set].size()][];
            for(int i=0; i<n.inputs.length; i++) {
                n.inputs[i]=inputs[set].get(i);
                n.outputs[i]=outputs[set].get(i);
                n.weights[i]=weights[set].get(i);
            }
            n.data=new float[][][] {n.inputs, n.outputs, n.weights};
            if(exclusive || set>0) {n.detectRange();}
        }
        
        if(!exclusive) {array[0]=m;}
        
        return array;
    }
    
    /**
     * Loads dataset from file(s) and creates SampleSets for training.
     * Handles file parsing, optional shuffling, train/test splitting, and class balancing.
     * 
     * @param fname File path(s) to load data from
     * @param maxLines0 Maximum lines to load initially
     * @param shuffleRaw Whether to shuffle data after loading
     * @param splitFraction Fraction for train/test split (0 for no split)
     * @param maxLines1 Maximum lines in validation set
     * @param exclusive Whether split sets should be independent
     * @param balance Class balancing factor (0 for no balancing, 1 for perfect balance)
     * @return Array of SampleSets: [training] or [training, validation]
     */
    public static SampleSet[] load(String fname, int maxLines0, boolean shuffleRaw, 
            float splitFraction, int maxLines1, boolean exclusive, float balance) {
        DataLoader dl=new DataLoader(fname);
        dl.load(maxLines0, shuffleRaw || splitFraction>0, balance);
        
        lastValidLines=dl.validLines;
        lastInvalidLines=dl.invalidLines;
        
        if(splitFraction<=0) {
            return new SampleSet[] {new SampleSet(dl.matrix)};
        }else {
            Matrix[] array=split(dl.matrix, splitFraction, maxLines1, exclusive);
            return new SampleSet[] {new SampleSet(array[0]), new SampleSet(array[1])};
        }
    }
    
    /**
     * Internal method to load and preprocess data from files.
     * Parses headers, processes data lines, applies shuffling and balancing.
     * 
     * @param maxLines Maximum number of lines to process
     * @param shuffleRaw Whether to shuffle the raw data
     * @param balance Class balancing multiplier
     */
    private void load(final int maxLines, final boolean shuffleRaw, final float balance) {
        matrix=new Matrix();
        ArrayList<float[]> inputList=new ArrayList<float[]>();
        ArrayList<float[]> outputList=new ArrayList<float[]>();
        ArrayList<float[]> weightList=new ArrayList<float[]>();
        byte[] s=null;
        final int max=(shuffleRaw ? Shared.MAX_ARRAY_LEN : maxLines);
        
        for(String f : fnames) {
            FileFormat ff=FileFormat.testInput(f, FileFormat.TEXT, null, true, false);
            ByteFile bf=ByteFile.makeByteFile(ff);
            for(s=bf.nextLine(); s!=null && validLines<max; s=bf.nextLine()){
                if(s.length>0) {
                    if(s[0]=='#') {//Header processing
                        if(Tools.startsWith(s, "#dims")) {
                            matrix.dims=parseIntArray(s, delimiter, true);
                            matrix.numInputs=matrix.dims[0];
                            matrix.numOutputs=matrix.dims[1];
                            weighted=(matrix.dims.length>2 && matrix.dims[2]==1);
                            assert(matrix.dims.length>1) : matrix.dims.length+", "+Arrays.toString(matrix.dims)+", '"+new String(s)+"'";
                        }else if(Tools.startsWith(s, "#inputs")) {
                            matrix.numInputs=parseInt(s);
                        }else if(Tools.startsWith(s, "#outputs")) {
                            matrix.numOutputs=parseInt(s);
                        }else if(Tools.startsWith(s, "##")) {
                            matrix.columns=new ArrayList<String>(Arrays.asList(new String(s).split("\t")));
                            matrix.columns.set(0, matrix.columns.get(0).substring(2));//Trim ##
                        }else {
                            //comment - ignore
                        }
                    }else {
                        // Auto-infer dimensions if not specified in header
                        if(matrix.numInputs==0) {
                            int terms=Tools.split(s, 0, (byte)'\t').size();
                            matrix.numOutputs=1;
                            matrix.numInputs=terms-(matrix.numOutputs+(weighted ? 1 : 0));
                            System.err.println("Inferring "+matrix.numInputs+" inputs, "+matrix.numOutputs+" output, "+(weighted ? 1 : 0)+" weights.");
                        }
                        assert(matrix.numInputs>0 & matrix.numOutputs>0) : 
                            "Number of inputs and outputs must be in data file header, e.g. '#inputs 5'";
                        float[] inputs=new float[matrix.numInputs];
                        float[] outputs=new float[matrix.numOutputs];
                        float[] weights=new float[] {1};
                        
                        boolean valid=parseDataLine(s, inputs, outputs, weights);
                        if(valid) {
                            inputList.add(inputs);
                            outputList.add(outputs);
                            weightList.add(weights);
                            validLines++;
                        }else {
                            invalidLines++;
                        }
                    }
                }
            }
            bf.close();
            if(validLines>=max) {break;}
        }
        
        if(shuffleRaw) {shuffle(inputList, outputList, weightList, maxLines);}
        if(balance>0) {balance(inputList, outputList, weightList, balance);}
        
        matrix.inputs=new float[inputList.size()][];
        matrix.outputs=new float[outputList.size()][];
        matrix.weights=new float[weightList.size()][];
        for(int i=0; i<matrix.inputs.length; i++) {
            matrix.inputs[i]=inputList.get(i);
            matrix.outputs[i]=outputList.get(i);
            matrix.weights[i]=weightList.get(i);
        }
        matrix.data=new float[][][] {matrix.inputs, matrix.outputs, matrix.weights};
        matrix.initializeRange();
    }
    
    /**
     * Randomly shuffles the dataset using a reproducible seed.
     * Maintains correspondence between inputs, outputs, and weights.
     * 
     * @param inputList List of input vectors
     * @param outputList List of output vectors  
     * @param weightList List of weight vectors
     * @param maxLines Maximum number of samples to retain after shuffling
     */
    private static void shuffle(ArrayList<float[]> inputList, ArrayList<float[]> outputList, ArrayList<float[]> weightList, int maxLines) {
        final int size=inputList.size();
        ArrayList<Triple> list=new ArrayList<Triple>(inputList.size());
        for(int i=0; i<inputList.size(); i++){
            Triple p=new Triple(inputList.get(i), outputList.get(i), weightList.get(i));
            list.add(p);
        }
        Random randy=new Random(SampleSet.shuffleSeed);
        Collections.shuffle(list, randy);
        inputList.clear();
        outputList.clear();
        weightList.clear();
        for(int i=0, lim=Tools.min(maxLines, list.size()); i<lim; i++) {
            Triple p=list.get(i);
            inputList.add(p.in);
            outputList.add(p.out);
            weightList.add(p.w);
        }
        assert(inputList.size()==size);
        assert(outputList.size()==size);
        assert(weightList.size()==size);
    }
    
    /**
     * Balances class representation by cloning underrepresented samples.
     * Intelligently handles edge cases including single-class datasets.
     * 
     * @param inputList List of input vectors to balance
     * @param outputList List of output vectors to balance
     * @param weightList List of weight vectors to balance
     * @param mult Target balance multiplier (1.0 = perfect balance)
     */
    private static void balance(ArrayList<float[]> inputList, ArrayList<float[]> outputList, ArrayList<float[]> weightList, float mult) {
        int pos=0, neg=0;
        assert(mult>0 && mult<=1);
        for(float[] out : outputList) {
            if(out[0]>=0.5f) {pos++;} else {neg++;}
        }
        final int target=(int)(Tools.max(pos, neg)*mult);
        if(pos>=target && neg>=target) {return;}
        assert(pos>=target || neg>=target);
        if(pos<1 || neg<1) {
            throw new RuntimeException("Can't balance with zero examples: pos="+pos+", neg="+neg);
        }
        for(int i=0; i<inputList.size() && (pos<target || neg<target); i++) {
            float[] in=inputList.get(i), out=outputList.get(i), weight=weightList.get(i);
            if(out[0]>=0.5 && pos<target) {
                pos++;
                inputList.add(in);
                outputList.add(out);
                weightList.add(weight);
            }else if(out[0]<0.5f && neg<target) {
                neg++;
                inputList.add(in);
                outputList.add(out);
                weightList.add(weight);
            }
        }
        shuffle(inputList, outputList, weightList, inputList.size());
    }
    
    /**
     * Parses a single data line into input, output, and weight arrays.
     * Handles tab-separated values with optional weight column.
     * 
     * @param line Raw byte array containing the data line
     * @param inputs Array to populate with input values
     * @param outputs Array to populate with output values
     * @param weights Array to populate with weight values
     * @return true if parsing was successful, false otherwise
     */
    boolean parseDataLine(byte[] line, float[] inputs, float[] outputs, float[] weights) {
        lp.set(line);
        int pos=0;
        for(int i=0; i<inputs.length; i++) {
            inputs[i]=lp.parseFloat(pos);
            pos++;
        }
        if(weighted) {
            weights[0]=lp.parseFloat(pos);
            pos++;
        }else {
            weights[0]=1;
        }
        for(int i=0; i<outputs.length; i++) {
            outputs[i]=lp.parseFloat(pos);
            pos++;
        }
        assert(pos==lp.terms()) : "\nExtra characters for line '"+new String(line)+
            "'; numInputs="+matrix.numInputs+", numOutputs="+matrix.numOutputs+
            ", "+pos+", "+line.length+"\n"+Arrays.toString(inputs)+"\n"+Arrays.toString(outputs)+"\n";
        return true;
    }
    
    /*--------------------------------------------------------------*/
    /*----------------           Parsing            ----------------*/
    /*--------------------------------------------------------------*/
    
    /**
     * Parses an integer value from a delimited line.
     * 
     * @param line Byte array containing the line to parse
     * @return Parsed integer value
     */
    private static int parseInt(byte[] line){
        int idx=Tools.indexOf(line, delimiter);
        return Parse.parseInt(line, idx+1, line.length);
    }
    
    /**
     * Parses an array of integers from a delimited line.
     * Optionally skips a title field at the beginning.
     * 
     * @param line Byte array containing the line to parse
     * @param delimiter Character used to separate fields
     * @param parseTitle Whether to skip the first field as a title
     * @return Array of parsed integer values
     */
    public static int[] parseIntArray(final byte[] line, final byte delimiter, boolean parseTitle){
        int a=0, b=0;
        IntList list=new IntList(3);
        
        if(parseTitle) {
            while(b<line.length && line[b]!=delimiter){b++;}
            assert(b>a) : "Missing Title: "+new String(line);
            b++;
            a=b;
        }
        
        while(a<line.length) {
            while(b<line.length && line[b]!=delimiter){b++;}
            assert(b>a) : "Missing element "+list.size+": '"+new String(line)+"'";
            int x=Parse.parseInt(line, a, b);
            list.add(x);
            b++;
            a=b;
        }
        return list.toArray();
    }
    
    /*--------------------------------------------------------------*/
    
    /**
     * Helper class to maintain correspondence between inputs, outputs, and weights during shuffling.
     */
    private static class Triple{
        /**
         * Creates a Triple containing input, output, and weight vectors.
         * 
         * @param in_ Input feature vector
         * @param out_ Output target vector
         * @param w_ Weight vector
         */
        Triple(float[] in_, float[] out_, float[] w_){
            in=in_;
            out=out_;
            w=w_;
        }
        
        /** Input feature vector */
        final float[] in;
        
        /** Output target vector */
        final float[] out;
        
        /** Sample weight vector */
        final float[] w;
    }
    
    /*--------------------------------------------------------------*/
    /*----------------           Fields             ----------------*/
    /*--------------------------------------------------------------*/
    
    /** Line parser for handling delimited data */
    LineParser1 lp=new LineParser1(delimiter);
    
    /** Array of file names to process */
    String fnames[];
    
    /** Loaded data matrix */
    Matrix matrix;
    
    /** Current parsing position */
    int pos=0;
    
    /** Count of successfully parsed lines */
    long validLines=0;
    
    /** Count of invalid or unparseable lines */
    long invalidLines=0;
    
    /** Valid lines from most recent load operation */
    static long lastValidLines=0;
    
    /** Invalid lines from most recent load operation */
    static long lastInvalidLines=0;
    
    /** Flag indicating whether weight columns are present */
    static boolean weighted=false;
    
    /** Field delimiter character */
    public static final byte delimiter='\t';
    
    /** Flag to ignore malformed lines instead of failing */
    public static boolean IGNORE_BAD_LINES=false;
}