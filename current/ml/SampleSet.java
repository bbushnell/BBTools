package ml;

import java.util.ArrayList;
import java.util.Random;

import shared.Shared;
import shared.Tools;

/**
 * Manages a set of machine learning samples for training and evaluation.
 * 
 * Provides functionality for:
 * - Sample generation
 * - Subset creation
 * - Performance metric calculations
 * - ROC curve generation
 * 
 * @author Brian Bushnell
 * @contributor Nepgear
 * @version 1.0
 */
public class SampleSet implements Cloneable {
	
	/**
	 * Constructs a SampleSet from a given matrix of input data.
	 * 
	 * @param m Matrix containing input, output, and weight data
	 * @throws AssertionError if matrix weights are null
	 */
	SampleSet(Matrix m) {
		matrix=m;
		assert(matrix.weights!=null);
	}

	/** Generates samples using the maximum possible number of inputs. */
	void makeSamples() {makeSamples(Integer.MAX_VALUE);}

	/**
	 * Generates samples from the input matrix, limited to a maximum number.
	 * Categorizes samples as positive or negative based on output midpoint.
	 * 
	 * @param max Maximum number of samples to generate
	 */
	void makeSamples(final int max) {
		assert(samples==null);
		numPositive=numNegative=0;
		samples=new Sample[Tools.min(max, matrix.inputs.length)];
		
		for(int i=0; i<samples.length; i++) {
			Sample s=new Sample(matrix.inputs[i], matrix.outputs[i], matrix.weights[i][0], i);
			samples[i]=s;
		}
		
		final float midpoint=matrix.outputMidpoint(); //TODO: Could use some cutoff for this instead, taken as a parameter
		for(Sample s : samples){
			if(s.goal[0]>=midpoint) {numPositive++;}
			else {numNegative++;}
		}
		
		samplesSortedByResult=samples.clone();
		shuffle();
	}
	
	/**
	 * Divides samples into multiple subsets for cross-validation.
	 * 
	 * @param numSets Number of subsets to create
	 */
	void makeSubsets(int numSets) {
		subsets=new Subset[numSets];
		@SuppressWarnings("unchecked")
		ArrayList<Sample>[] lists=new ArrayList[numSets];
		for(int i=0; i<lists.length; i++) {
			lists[i]=new ArrayList<Sample>();
		}
		for(int i=0, j=0; i<samples.length; i++, j++) {
			j=j%numSets;
			lists[j].add(samples[i]);
		}
		for(int i=0; i<lists.length; i++) {
			subsets[i]=new Subset(lists[i]);
		}
		currentSubset=0;
	}
	
	/**
	 * Retrieves the current subset for a given training epoch.
	 * 
	 * @param epoch Current training epoch
	 * @return Current subset of samples
	 */
	Subset currentSubset(int epoch) {
		if(epoch>=nextSubsetEpoch){
			advanceSubset();
			nextSubsetEpoch=epoch+subsetInterval;
		}
		return subsets[currentSubset];
	}
	
	/** Advances to the next subset, potentially shuffling samples. */
	private void advanceSubset() {
		currentSubset=(currentSubset+1)%subsets.length;
		assert(currentSubset>=0);
		if(currentSubset==0 && shuffle && subsets.length>1) {
			shuffle();
			makeSubsets(subsets.length);
//			System.err.println("Shuffled");
		}
	}
	
	/** Randomly shuffles samples using a reproducible seed. */
	private void shuffle() {
		final long seed=numShuffles^Long.rotateLeft(shuffleSeed, 17);
//		System.err.println("Shuffled ("+seed+")");
//		new Exception().printStackTrace();
//		if(true) {return;}
		Random randy=new Random(seed);
		for(int i=0; i<samples.length; i++) {
			int idx=randy.nextInt(samples.length);
			Sample s=samples[idx];
			samples[idx]=samples[i];
			samples[i]=s;
		}
//		System.err.println(samples[0].positive);
		numShuffles++;
	}
	
	/** Sorts samples by their predicted value. */
	public void sortByValue() {
		Shared.sort(samplesSortedByResult, SampleValueComparator.COMPARATOR);
//		assert(checkSort());
	}
	
	/**
	 * Verifies that samples are correctly sorted by value.
	 * 
	 * @return true if samples are sorted, false otherwise
	 */
	public boolean checkSort() {
		for(int i=1; i<samplesSortedByResult.length; i++) {
			if(samplesSortedByResult[i-1].result[0]>samplesSortedByResult[i].result[0]) {
				assert(false) : i+"\n"+samplesSortedByResult[0].result[0]
						+"\n"+samplesSortedByResult[1].result[0]
								+"\n"+samplesSortedByResult[2].result[0]
										+"\n"+samplesSortedByResult[3].result[0];
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Calculates False Negative Rate for a given cutoff value.
	 * 
	 * @param cutoff Threshold for classification
	 * @return False Negative Rate
	 */
	public double calcFNRFromCutoff(final double cutoff) {
		//Should be sorted
		int fn=0, tn=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		for(int i=0; i<samplesSortedByResult.length; i++) {
			Sample s=samplesSortedByResult[i];
			if(s.result[0]>=cutoff) {
//				System.err.println(i+", "+s.result[0]+", "+cutoff);
				break;
			}
			if(s.positive) {fn++;}
			else {tn++;}
		}
//		assert(false) : samplesSortedByResult.length+", "+fn+", "+tn+", "+
//		samplesSortedByResult[0].result[0]+", "+samplesSortedByResult[samplesSortedByResult.length-1].result[0];
		return fn*invSamples;
	}
	
	/**
	 * Calculates cutoff value based on crossover point.
	 * 
	 * @param fpMult False positive multiplier
	 * @return Calculated cutoff value
	 */
	public float calcCutoffFromCrossover(double fpMult) {
		//Should be sorted
		int pos=0, neg=numNegative, i=0;
		for(; i<samplesSortedByResult.length; i++) {
			Sample s=samplesSortedByResult[i];
			if(s.positive) {
				pos++;
			}else{
				neg--;
			}
			if(pos*fpMult>=neg) {break;}
		}
		if(i==0) {return samplesSortedByResult[i].result[0];}
		if(pos*fpMult==neg) {return samplesSortedByResult[i].result[0]-0.00001f;}
		Sample a=samplesSortedByResult[i-1];
		Sample b=samplesSortedByResult[i];
		return 0.5f*(a.result[0]+b.result[0]);
	}
	
	/**
	 * Calculates False Positive Rate for a given cutoff.
	 * 
	 * @param cutoff Threshold for classification
	 * @return False Positive Rate
	 */
	public double calcFPRFromCutoff(final double cutoff) {
		//Should be sorted
		int fp=0, tp=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		for(int i=samplesSortedByResult.length-1; i>=0; i--) {
			Sample s=samplesSortedByResult[i];
			if(s.result[0]<cutoff) {
//				System.err.println(i+", "+s.result[0]+", "+cutoff);
				break;
			}
			if(s.positive) {tp++;}
			else {fp++;}
		}
//		assert(false) : samplesSortedByResult.length+", "+fp+", "+tp+", "+
//		samplesSortedByResult[0].result[0]+", "+samplesSortedByResult[samplesSortedByResult.length-1].result[0];
		return fp*invSamples;
	}
	
	public double calcCutoffFromFPR(final double fpr) {
		//Should be sorted
		int fp=0, tp=0;
//		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fpr*samplesSortedByResult.length);
		float lastCutoff=1.0f, prev=1.0f;
		for(int i=samplesSortedByResult.length-1; i>=0 && fp<=target; i--) {
			Sample s=samplesSortedByResult[i];
			if(s.positive) {tp++;}
			else {fp++;}
			prev=lastCutoff;
			lastCutoff=s.result[0];
		}
		
		return (lastCutoff+prev)*0.5;
	}
	
	public double calcCutoffFromFNR(final double fnr) {
		//Should be sorted
		int fn=0, tn=0;
//		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fnr*samplesSortedByResult.length);
		float lastCutoff=0.0f, prev=0.0f;
		for(int i=0; i<samplesSortedByResult.length && fn<=target; i++) {
			Sample s=samplesSortedByResult[i];
			if(s.positive) {fn++;}
			else {tn++;}
			prev=lastCutoff;
			lastCutoff=s.result[0];
		}
		
		return (lastCutoff+prev)*0.5;
	}
	
	public double calcFNRFromFPR(final double fpr) {
		//Should be sorted
		int fp=0, tp=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fpr*samplesSortedByResult.length);
		for(int i=samplesSortedByResult.length-1; i>=0 && fp<=target; i--) {
			Sample s=samplesSortedByResult[i];
			if(s.positive) {tp++;}
			else {fp++;}
		}
		
		int tn=numNegative-fp;
		int fn=numPositive-tp;
		return fn*invSamples;
	}
	
	public double calcFPRFromFNR(final double fnr) {
		//Should be sorted
		int fn=0, tn=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fnr*samplesSortedByResult.length);
		for(int i=0; i<samplesSortedByResult.length && fn<=target; i++) {
			Sample s=samplesSortedByResult[i];
			if(s.positive) {fn++;}
			else {tn++;}
		}
		
		int fp=numNegative-tn;
		int tp=numPositive-fn;
		return fp*invSamples;
	}
	
	// Additional methods continue here...

	/**
	 * Generates Receiver Operating Characteristic (ROC) curve.
	 * 
	 * @param bins Number of bins for ROC curve
	 * @return ROC curve data points
	 */
	public float[] calcROC(int bins) {
		bins=Tools.max(bins, 2);
//		final double invSamples=1.0/samplesSortedByResult.length;
		final double invPositive=1.0/numPositive;
		double binlen=(numNegative-1)/(double)(bins-1);
		float[] roc=new float[bins];
		
		//Should be sorted
		int pcount=0, ncount=0;
		int bin=1;
		int target=(int)Math.round(binlen*bin);
		for(int i=0; i<samplesSortedByResult.length; i++) {
			Sample s=samplesSortedByResult[i];
			if(s.positive) {pcount++;}
			else {ncount++;}
			while(ncount>=target) {
				roc[bin]=(float)(1-pcount*invPositive);
//				System.err.println("i="+i+", bin="+bin+", target="+target+", val="+roc[bin]+", pred="+s.result[0]);
				bin++;
				target=(int)Math.round(binlen*bin);
			}
		}
		roc[0]=1;
//		roc[roc.length-1]=0;//This is not strictly necessary; doesn't need to start at 0/0.
		Tools.reverseInPlace(roc);
		return roc;
	}
	
	/** Creates a copy of the current sample set. */
	SampleSet copy() {return copy(Integer.MAX_VALUE, 1f);}
	
	/**
	 * Creates a copy of the sample set with optional size and subset limitations.
	 * 
	 * @param maxSamples Maximum number of samples to copy
	 * @param subsetSizeFraction Fraction of subset size to preserve
	 * @return Copied sample set
	 */
	SampleSet copy(int maxSamples, float subsetSizeFraction) {
		SampleSet copy=null;
		try {
			copy = (SampleSet) this.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		copy.samples=null;
		copy.samplesSortedByResult=null;
		copy.subsets=null;
		copy.numShuffles=0;
		copy.currentSubset=0;
		copy.reset();
		if(samples!=null) {
			copy.makeSamples(maxSamples);
//			System.err.println("First="+copy.samples[0].id+", "+copy.samples[0].positive);
			if(subsets!=null && subsetSizeFraction>0) {
				int numSubsets=subsets.length;
				if(maxSamples<0.5f*samples.length) {
					numSubsets=Tools.max(1, (int)Math.ceil(numSubsets*(maxSamples/(subsetSizeFraction*samples.length))));
				}
				copy.makeSubsets(numSubsets);
			}
		}
		return copy;
	}
	
	/** Resets the sample set to its initial state. */
	void reset() {
		currentSubset=0;
		numShuffles=0;
		nextSubsetEpoch=subsetInterval;
		if(subsets!=null) {
			for(Subset ss : subsets) {
				ss.reset();
			}
		}
	}

	/**
	 * Returns the maximum size of a subset.
	 * 
	 * @return Maximum subset size
	 */
	public int maxSubsetSize() {
		return subsets==null || subsets.length<1 || subsets[0].samples==null ? 0 : subsets[0].samples.length;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final int numInputs() {return matrix.numInputs();}
	final int numOutputs() {return matrix.numOutputs();}
	final float outputMidpoint() {return matrix.outputMidpoint();}

	final Matrix matrix;
	int numPositive=0, numNegative=0;
	Sample[] samples;
	Sample[] samplesSortedByResult;
	private Subset[] subsets;
	
	private int currentSubset=0;
	private int numShuffles=0;
	static long shuffleSeed=0;
	static boolean shuffle=true;
	
	int nextSubsetEpoch=subsetInterval;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	static int subsetInterval=64;
	
}
