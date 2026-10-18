package test.ml;

import ml.Matrix;
import ml.Sample;
import ml.SampleSet;
import ml.Subset;

import java.util.ArrayList;

/**
 * Historical pre-ML012 -ea baseline for adjusted-output Sample positivity
 * inconsistencies. Cases mutate Matrix globals, so execute them serially;
 * run(name) lets an external harness isolate cases that current code has fixed,
 * while no-arg main stops at the first failure.
 * @author Sucrose
 */
public class SampleSetMidpointConsistencyTest{

	public static void main(String[] args){
		if(args.length==1){
			run(args[0]);
			return;
		}
		assert(args.length==0) : "This fixture accepts only one optional case name";
		run("cached_classification");
		run("fnr_symptom");
		run("default01_control");
		run("subset_bucket_consistency");
	}

	/** Runs one selectable case so a baseline failure cannot mask later rows. */
	private static void run(String name){
		final boolean oldConvertTo01=Matrix.convertTo01;
		final float oldTargetMin=Matrix.targetOutputRangeMin;
		final float oldTargetMax=Matrix.targetOutputRangeMax;
		final boolean oldSetMin=Matrix.setTargetOutputRangeMin;
		final boolean oldSetMax=Matrix.setTargetOutputRangeMax;
		try{
			if(name.equals("cached_classification")){
				checkCachedClassification();
			}else if(name.equals("fnr_symptom")){
				checkFNR();
			}else if(name.equals("default01_control")){
				checkDefault01();
			}else if(name.equals("subset_bucket_consistency")){
				checkSubsetBuckets();
			}else{
				throw new IllegalArgumentException("Unknown fixture mode: "+name);
			}
			System.err.println("SampleSetMidpointConsistencyTest PASS "+name);
		}finally{
			Matrix.convertTo01=oldConvertTo01;
			Matrix.targetOutputRangeMin=oldTargetMin;
			Matrix.targetOutputRangeMax=oldTargetMax;
			Matrix.setTargetOutputRangeMin=oldSetMin;
			Matrix.setTargetOutputRangeMax=oldSetMax;
		}
	}

	/** Captures the split between midpoint totals and cached 0.5-based Samples. */
	private static void checkCachedClassification(){
		final Matrix m=makeMatrix(true);
		final SampleSet set=makeSampleSet(m);
		final float midpoint=m.outputMidpoint();
		final int cached=countPositive(set.samples);
		System.out.println("cached_classification\tmidpoint="+midpoint+
			"\tnumPositive="+set.numPositive+"\tnumNegative="+set.numNegative+
			"\tcachedPositive="+cached);
		assert(midpoint==3f) : "Adjusted 2..4 outputs should have midpoint 3, not "+midpoint;
		assert(set.numPositive==1 && set.numNegative==1) :
			"Midpoint counts should split 2/4 outputs around 3: "+
			set.numPositive+", "+set.numNegative;
		assert(cached==2) :
			"Sample caches are currently 0.5-based, so remapped 2/4 goals both cache positive";
	}

	/** Captures the negative FNR produced by mixed class definitions. */
	private static void checkFNR(){
		final SampleSet set=makeSampleSet(makeMatrix(true));
		copyGoalToResult(set);
		set.sortByValue();
		final double fnr=set.calcFNRFromFPR(0);
		System.out.println("fnr_symptom\tfnrFromFpr0="+fnr);
		assert(fnr==-0.5) :
			"Baseline adjusted outputs should expose impossible FNR -0.5, not "+fnr;
	}

	/** Default 0/1 targets use midpoint 0.5 and remain internally consistent. */
	private static void checkDefault01(){
		final Matrix m=makeMatrix(false);
		final SampleSet set=makeSampleSet(m);
		final float midpoint=m.outputMidpoint();
		final int cached=countPositive(set.samples);
		copyGoalToResult(set);
		set.sortByValue();
		final double fnr=set.calcFNRFromFPR(0);
		System.out.println("default01_control\tmidpoint="+midpoint+
			"\tnumPositive="+set.numPositive+"\tnumNegative="+set.numNegative+
			"\tcachedPositive="+cached+"\tfnrFromFpr0="+fnr);
		assert(midpoint==0.5f) : "Default 0/1 outputs should have midpoint 0.5";
		assert(set.numPositive==1 && set.numNegative==1) :
			"Default midpoint counts should split 0/1 outputs: "+
			set.numPositive+", "+set.numNegative;
		assert(cached==1) : "Default 0/1 Sample caches should keep one positive";
		assert(fnr==0) : "Default 0/1 outputs should have zero FNR at FPR 0, not "+fnr;
	}

	/** Historical expectation: Subset rederived adjusted classes from goal[0]>=0.5. */
	private static void checkSubsetBuckets(){
		final SampleSet set=makeSampleSet(makeMatrix(true));
		final ArrayList<Sample> list=new ArrayList<Sample>();
		for(Sample s : set.samples){list.add(s);}
		final Subset subset=new Subset(list);
		System.out.println("subset_bucket_consistency\tnumPositive="+set.numPositive+
			"\tnumNegative="+set.numNegative+"\tsubsetPositive="+
			subset.positive.length+"\tsubsetNegative="+subset.negative.length);
		assert(set.numPositive==1 && set.numNegative==1) :
			"SampleSet midpoint counts should split adjusted 2/4 outputs";
		assert(subset.positive.length==2 && subset.negative.length==0) :
			"Subset currently buckets adjusted 2/4 goals with a 0.5 cutoff: "+
			subset.positive.length+", "+subset.negative.length;
	}

	/** Builds and materializes the minimal two-row matrix. */
	private static SampleSet makeSampleSet(final Matrix m){
		final SampleSet set=new SampleSet(m);
		set.makeSamples();
		return set;
	}

	/** Constructs a two-row matrix, optionally remapping outputs to 2..4. */
	private static Matrix makeMatrix(final boolean adjusted){
		Matrix.convertTo01=false;
		Matrix.targetOutputRangeMin=0;
		Matrix.targetOutputRangeMax=0;
		Matrix.setTargetOutputRangeMin=false;
		Matrix.setTargetOutputRangeMax=false;
		if(adjusted){
			Matrix.targetOutputRangeMin=2;
			Matrix.targetOutputRangeMax=4;
			Matrix.setTargetOutputRangeMin=true;
			Matrix.setTargetOutputRangeMax=true;
		}
		final Matrix m=new Matrix();
		m.columns=new ArrayList<String>();
		m.columns.add("input");
		m.columns.add("target");
		m.numInputs=1;
		m.numOutputs=1;
		m.inputs=new float[][]{{0}, {1}};
		m.outputs=new float[][]{{0}, {1}};
		m.weights=new float[][]{{1}, {1}};
		m.data=new float[][][]{m.inputs, m.outputs, m.weights};
		m.validLines=m.inputs.length;
		m.initializeRange();
		return m;
	}

	/** Counts the Sample.positive flags cached by the legacy makeSamples path. */
	private static int countPositive(final Sample[] samples){
		int count=0;
		for(Sample s : samples){
			if(s.positive){count++;}
		}
		return count;
	}

	/** Copies each goal into the result slot used by threshold-ranking checks. */
	private static void copyGoalToResult(final SampleSet set){
		for(Sample s : set.samples){
			s.result[0]=s.goal[0];
		}
	}
}
