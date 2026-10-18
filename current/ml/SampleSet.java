package ml;

import java.util.ArrayList;
import shared.Random;

import shared.Shared;
import shared.Tools;

/**
 * Matrix-backed owner for lazily materialized {@link Sample} rows, randomized
 * training subsets, and result-sorted threshold helpers.
 * Positive/negative counts follow {@link Sample#positive} after any requested
 * row truncation.
 */
public class SampleSet implements Cloneable{
	
	/** Wraps {@code m} without materializing {@link Sample} rows. */
	public SampleSet(Matrix m){
		matrix=m;
		assert(matrix.weights!=null);
	}

	/** Materializes all matrix rows as Samples and shuffles them for training. */
	public void makeSamples(){makeSamples(Integer.MAX_VALUE);}

	/** Materializes all matrix rows with the supplied positivity cutoff. */
	public void makeSamples(final float positiveCutoff_){makeSamples(Integer.MAX_VALUE, positiveCutoff_);}

	/** Materializes at most {@code max} matrix rows before counting classes and shuffling. */
	public void makeSamples(final int max){makeSamples(max, 0.5f);}

	/** Materializes at most {@code max} matrix rows with the supplied positivity cutoff. */
	void makeSamples(final int max, final float positiveCutoff_){
		assert(samples==null);
		positiveCutoff=positiveCutoff_;
		numPositive=numNegative=0;
		samples=new Sample[Tools.min(max, matrix.inputs.length)];
		
		for(int i=0; i<samples.length; i++){
			Sample s=new Sample(matrix.inputs[i], matrix.outputs[i], matrix.weights[i][0], i, positiveCutoff);
			samples[i]=s;
			if(s.positive){numPositive++;}
			else{numNegative++;}
		}
		
		samplesSortedByResult=samples.clone();
		shuffle();
	}
	
	/** Partitions the current shuffled sample order round-robin into {@code numSets}. */
	void makeSubsets(int numSets){
		subsets=new Subset[numSets];
		@SuppressWarnings("unchecked")
		ArrayList<Sample>[] lists=new ArrayList[numSets];
		for(int i=0; i<lists.length; i++){
			lists[i]=new ArrayList<Sample>();
		}
		for(int i=0, j=0; i<samples.length; i++, j++){
			j=j%numSets;
			lists[j].add(samples[i]);
		}
		for(int i=0; i<lists.length; i++){
			subsets[i]=new Subset(lists[i]);
		}
		currentSubset=0;
	}
	
	/** Returns the active subset, advancing one partition when {@code epoch} reaches the next boundary. */
	Subset currentSubset(int epoch){
		if(epoch>=nextSubsetEpoch){
			advanceSubset();
			nextSubsetEpoch=epoch+subsetInterval;
		}
		return subsets[currentSubset];
	}
	
	/** Advances one subset and reshuffles/repartitions after wrapping. */
	private void advanceSubset(){
		currentSubset=(currentSubset+1)%subsets.length;
		assert(currentSubset>=0);
		if(currentSubset==0 && shuffle && subsets.length>1){
			shuffle();
			makeSubsets(subsets.length);
//			System.err.println("Shuffled");
		}
	}
	
	/** Deterministically shuffles {@link #samples} with a per-call seed. */
	private void shuffle(){
		final long seed=numShuffles^Long.rotateLeft(shuffleSeed, 17);
//		System.err.println("Shuffled ("+seed+")");
//		new Exception().printStackTrace();
//		if(true){return;}
		Random randy=shared.Shared.random(seed);
		for(int i=0; i<samples.length; i++){
			int idx=randy.nextInt(samples.length);
			Sample s=samples[idx];
			samples[idx]=samples[i];
			samples[i]=s;
		}
//		System.err.println(samples[0].positive);
		numShuffles++;
	}
	
	/** Sorts {@link #samplesSortedByResult} by ascending primary result. */
	public void sortByValue(){
		Shared.sort(samplesSortedByResult, SampleValueComparator.COMPARATOR);
//		assert(checkSort());
	}
	
	/** Checks the ordering established by {@link #sortByValue()}. */
	public boolean checkSort(){
		for(int i=1; i<samplesSortedByResult.length; i++){
			if(samplesSortedByResult[i-1].result[0]>samplesSortedByResult[i].result[0]){
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
	 * Returns false negatives per total sample from the current
	 * {@link #sortByValue()} ordering.
	 */
	public double calcFNRFromCutoff(final double cutoff){
		//Should be sorted
		int fn=0, tn=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		for(int i=0; i<samplesSortedByResult.length; i++){
			Sample s=samplesSortedByResult[i];
			if(s.result[0]>=cutoff){
//				System.err.println(i+", "+s.result[0]+", "+cutoff);
				break;
			}
			if(s.positive){fn++;}
			else{tn++;}
		}
//		assert(false) : samplesSortedByResult.length+", "+fn+", "+tn+", "+
//		samplesSortedByResult[0].result[0]+", "+samplesSortedByResult[samplesSortedByResult.length-1].result[0];
		return fn*invSamples;
	}
	
	/** Returns the {@link #sortByValue()} cutoff where positives scaled by {@code fpMult} meet remaining negatives. */
	public float calcCutoffFromCrossover(double fpMult){
		//Should be sorted
		int pos=0, neg=numNegative, i=0;
		for(; i<samplesSortedByResult.length; i++){
			Sample s=samplesSortedByResult[i];
			if(s.positive){
				pos++;
			}else{
				neg--;
			}
			if(pos*fpMult>=neg){break;}
		}
		if(i==0){return samplesSortedByResult[i].result[0];}
		if(pos*fpMult==neg){return samplesSortedByResult[i].result[0]-0.00001f;}
		Sample a=samplesSortedByResult[i-1];
		Sample b=samplesSortedByResult[i];
		return 0.5f*(a.result[0]+b.result[0]);
	}
	
	/**
	 * Returns false positives per total sample from the current
	 * {@link #sortByValue()} ordering.
	 */
	public double calcFPRFromCutoff(final double cutoff){
		//Should be sorted
		int fp=0, tp=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		for(int i=samplesSortedByResult.length-1; i>=0; i--){
			Sample s=samplesSortedByResult[i];
			if(s.result[0]<cutoff){
//				System.err.println(i+", "+s.result[0]+", "+cutoff);
				break;
			}
			if(s.positive){tp++;}
			else{fp++;}
		}
//		assert(false) : samplesSortedByResult.length+", "+fp+", "+tp+", "+
//		samplesSortedByResult[0].result[0]+", "+samplesSortedByResult[samplesSortedByResult.length-1].result[0];
		return fp*invSamples;
	}
	
	/** Returns a cutoff targeting the requested false-positive rate after {@link #sortByValue()}. */
	public double calcCutoffFromFPR(final double fpr){
		//Should be sorted
		int fp=0, tp=0;
//		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fpr*samplesSortedByResult.length);
		float lastCutoff=1.0f, prev=1.0f;
		for(int i=samplesSortedByResult.length-1; i>=0 && fp<=target; i--){
			Sample s=samplesSortedByResult[i];
			if(s.positive){tp++;}
			else{fp++;}
			prev=lastCutoff;
			lastCutoff=s.result[0];
		}
		
		return (lastCutoff+prev)*0.5;
	}
	
	/** Returns a cutoff targeting the requested false-negative rate after {@link #sortByValue()}. */
	public double calcCutoffFromFNR(final double fnr){
		//Should be sorted
		int fn=0, tn=0;
//		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fnr*samplesSortedByResult.length);
		float lastCutoff=0.0f, prev=0.0f;
		for(int i=0; i<samplesSortedByResult.length && fn<=target; i++){
			Sample s=samplesSortedByResult[i];
			if(s.positive){fn++;}
			else{tn++;}
			prev=lastCutoff;
			lastCutoff=s.result[0];
		}
		
		return (lastCutoff+prev)*0.5;
	}
	
	/** Returns false negatives per total sample at a false-positive target after {@link #sortByValue()}. */
	public double calcFNRFromFPR(final double fpr){
		//Should be sorted
		int fp=0, tp=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fpr*samplesSortedByResult.length);
		for(int i=samplesSortedByResult.length-1; i>=0 && fp<=target; i--){
			Sample s=samplesSortedByResult[i];
			if(s.positive){tp++;}
			else{fp++;}
		}
		
		int tn=numNegative-fp;
		int fn=numPositive-tp;
		return fn*invSamples;
	}
	
	/** Returns false positives per total sample at a false-negative target after {@link #sortByValue()}. */
	public double calcFPRFromFNR(final double fnr){
		//Should be sorted
		int fn=0, tn=0;
		final double invSamples=1.0/samplesSortedByResult.length;
		final int target=(int)Math.floor(fnr*samplesSortedByResult.length);
		for(int i=0; i<samplesSortedByResult.length && fn<=target; i++){
			Sample s=samplesSortedByResult[i];
			if(s.positive){fn++;}
			else{tn++;}
		}
		
		int fp=numNegative-tn;
		int tp=numPositive-fn;
		return fp*invSamples;
	}
	
	/** Builds a coarse ROC curve from the current {@link #sortByValue()} ordering. */
	public float[] calcROC(int bins){
		bins=Tools.max(bins, 2);
//		final double invSamples=1.0/samplesSortedByResult.length;
		final double invPositive=1.0/numPositive;
		double binlen=(numNegative-1)/(double)(bins-1);
		float[] roc=new float[bins];
		
		//Should be sorted
		int pcount=0, ncount=0;
		int bin=1;
		int target=(int)Math.round(binlen*bin);
		for(int i=0; i<samplesSortedByResult.length; i++){
			Sample s=samplesSortedByResult[i];
			if(s.positive){pcount++;}
			else{ncount++;}
			//TODO: Probable bug (ML048) - bin may advance past roc.length when small negative counts repeat rounded targets.
			while(ncount>=target){
				assert(bin<roc.length) : "calcROC output-bin overrun; repeated rounded "+
					"targets can advance bin past length (ML048): bin="+bin+
					", length="+roc.length+", ncount="+ncount+", target="+target;
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
	
	/**
	 * Shallow-copies this set and rematerializes all Sample wrappers only if
	 * this set already has materialized samples.
	 */
	SampleSet copy(){return copy(Integer.MAX_VALUE, 1f);}
	
	/**
	 * Shallow-copies this set and, only if the original samples already exist,
	 * rematerializes up to {@code maxSamples} rows and rebuilds subsets only when
	 * present and {@code subsetSizeFraction}>0.
	 */
	public SampleSet copy(int maxSamples, float subsetSizeFraction){
		SampleSet copy=null;
		try{
			copy=(SampleSet)this.clone();
		}catch(CloneNotSupportedException e){
			e.printStackTrace();
		}
		copy.samples=null;
		copy.samplesSortedByResult=null;
		copy.subsets=null;
		copy.numShuffles=0;
		copy.currentSubset=0;
		copy.reset();
		if(samples!=null){
			copy.makeSamples(maxSamples, positiveCutoff);
//			System.err.println("First="+copy.samples[0].id+", "+copy.samples[0].positive);
			if(subsets!=null && subsetSizeFraction>0){
				int numSubsets=subsets.length;
				if(maxSamples<0.5f*samples.length){
					numSubsets=Tools.max(1, (int)Math.ceil(numSubsets*(maxSamples/(subsetSizeFraction*samples.length))));
				}
				copy.makeSubsets(numSubsets);
			}
		}
		return copy;
	}
	
	/** Restarts subset iteration and clears per-subset progress. */
	void reset(){
		currentSubset=0;
		numShuffles=0;
		nextSubsetEpoch=subsetInterval;
		if(subsets!=null){
			for(Subset ss : subsets){
				ss.reset();
			}
		}
	}

	/** Returns the first subset length, or zero before subset construction. */
	public int maxSubsetSize(){
		return subsets==null || subsets.length<1 || subsets[0].samples==null ? 0 : subsets[0].samples.length;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final int numInputs(){return matrix.numInputs();}
	final int numOutputs(){return matrix.numOutputs();}
	public final float outputMidpoint(){return matrix.outputMidpoint();}

	/** Matrix backing store; materialized Samples borrow its row arrays. */
	public final Matrix matrix;
	/** Count of positive samples in the dataset; negative count is tracked separately. */
	public int numPositive=0, numNegative=0;
	/** Materialized rows in the current randomized training order. */
	public Sample[] samples;
	/** Shallow copy of {@link #samples} retained for threshold sorting. */
	Sample[] samplesSortedByResult;
	/** Round-robin partitions of {@link #samples}, rebuilt after full reshuffles. */
	private Subset[] subsets;
	/** Last materialization cutoff, reused by sample copies. */
	private float positiveCutoff=0.5f;
	
	/** Active index into {@link #subsets}. */
	private int currentSubset=0;
	/** Number of completed shuffle calls, used to advance the deterministic seed. */
	private int numShuffles=0;
	/** Base seed mixed with {@link #numShuffles}. */
	static long shuffleSeed=0;
	/** Enables reshuffling when subset rotation wraps around. */
	static boolean shuffle=true;
	
	/** Epoch at which the next call to {@link #currentSubset(int)} advances. */
	int nextSubsetEpoch=subsetInterval;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of epochs between subset rotations. */
	static int subsetInterval=64;
	
}
