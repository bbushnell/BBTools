package test.ml;

import ml.Cell;
import ml.Matrix;
import ml.Sample;
import ml.SampleSet;
import ml.Subset;
import ml.Trainer;

import java.util.ArrayList;

/**
 * Selectable -ea expectations for SampleSet training/evaluation cutoffs and
 * materialized positive/negative caches. Cases mutate Matrix, Trainer and Cell
 * globals, so execute them serially; run(name) lets an external harness isolate
 * cases, while no-arg main stops at the first failure.
 * @author Sucrose
 */
public class SampleSetCutoffInitializationTest{

	public static void main(String[] args){
		if(args.length==1){
			run(args[0]);
			return;
		}
		assert(args.length==0) : "This fixture accepts only one optional case name";
		run("default_false");
		run("midpoint_true");
		run("explicit_cutoffs");
		run("validation_range");
		run("data_null");
		run("copy_subset");
		run("legacy_entrypoint");
		run("explicit_training_only");
		run("explicit_eval_only");
	}

	/** Runs one selectable case so a failure cannot mask later rows. */
	private static void run(String name){
		final boolean oldConvertTo01=Matrix.convertTo01;
		final float oldTargetMin=Matrix.targetOutputRangeMin;
		final float oldTargetMax=Matrix.targetOutputRangeMax;
		final boolean oldSetMin=Matrix.setTargetOutputRangeMin;
		final boolean oldSetMax=Matrix.setTargetOutputRangeMax;
		final float oldEvalCutoff=Trainer.cutoffForEvaluation;
		final boolean oldSetEval=Trainer.setCutoffForEvaluation;
		final float oldTrainingCutoff=Cell.cutoffForTraining;
		final boolean oldSetTraining=Cell.setCutoffForTraining;
		final boolean oldUseMidpoint=Cell.useMidpoint;
		try{
			resetGlobals();
			if(name.equals("default_false")){
				checkDefaultFalse();
			}else if(name.equals("midpoint_true")){
				checkMidpointTrue();
			}else if(name.equals("explicit_cutoffs")){
				checkExplicitCutoffs();
			}else if(name.equals("validation_range")){
				checkValidationRange();
			}else if(name.equals("data_null")){
				checkDataNull();
			}else if(name.equals("copy_subset")){
				checkCopySubset();
			}else if(name.equals("legacy_entrypoint")){
				checkLegacyEntrypoint();
			}else if(name.equals("explicit_training_only")){
				checkExplicitTrainingOnly();
			}else if(name.equals("explicit_eval_only")){
				checkExplicitEvalOnly();
			}else{
				throw new IllegalArgumentException("Unknown fixture mode: "+name);
			}
			System.err.println("SampleSetCutoffInitializationTest PASS "+name);
		}finally{
			Matrix.convertTo01=oldConvertTo01;
			Matrix.targetOutputRangeMin=oldTargetMin;
			Matrix.targetOutputRangeMax=oldTargetMax;
			Matrix.setTargetOutputRangeMin=oldSetMin;
			Matrix.setTargetOutputRangeMax=oldSetMax;
			Trainer.cutoffForEvaluation=oldEvalCutoff;
			Trainer.setCutoffForEvaluation=oldSetEval;
			Cell.cutoffForTraining=oldTrainingCutoff;
			Cell.setCutoffForTraining=oldSetTraining;
			Cell.useMidpoint=oldUseMidpoint;
		}
	}

	/** With usemidpoint=f, adjusted 2/4 rows preserve legacy 0.5 classing. */
	private static void checkDefaultFalse(){
		final SampleSet data=new SampleSet(makeMatrix(2, 4, 0, 1));
		Trainer.initializeSampleCutoffs(data, null);
		System.out.println("default_false\tcutoff="+Cell.cutoffForTraining+
			"\tnumPositive="+data.numPositive+"\tnumNegative="+data.numNegative+
			"\tcachedPositive="+countPositive(data.samples));
		assert(Cell.cutoffForTraining==0.5f) : Cell.cutoffForTraining;
		assert(Trainer.cutoffForEvaluation==0.5f) : Trainer.cutoffForEvaluation;
		assertCounts(data, 2, 0);
	}

	/** With usemidpoint=t, adjusted 2/4 training rows use midpoint 3. */
	private static void checkMidpointTrue(){
		Cell.useMidpoint=true;
		final SampleSet data=new SampleSet(makeMatrix(2, 4, 0, 1));
		Trainer.initializeSampleCutoffs(data, null);
		copyGoalToResult(data);
		data.sortByValue();
		final double fnr=data.calcFNRFromFPR(0);
		System.out.println("midpoint_true\tcutoff="+Cell.cutoffForTraining+
			"\tnumPositive="+data.numPositive+"\tnumNegative="+data.numNegative+
			"\tcachedPositive="+countPositive(data.samples)+"\tfnrFromFpr0="+fnr);
		assert(Cell.cutoffForTraining==3f) : Cell.cutoffForTraining;
		assert(Trainer.cutoffForEvaluation==3f) : Trainer.cutoffForEvaluation;
		assertCounts(data, 1, 1);
		assert(fnr==0) : fnr;
	}

	/** Explicit training and evaluation cutoffs override usemidpoint=t. */
	private static void checkExplicitCutoffs(){
		Cell.useMidpoint=true;
		Cell.cutoffForTraining=2.5f;
		Cell.setCutoffForTraining=true;
		Trainer.cutoffForEvaluation=3.5f;
		Trainer.setCutoffForEvaluation=true;
		final SampleSet data=new SampleSet(makeMatrix(2, 4, 0, 0.5f, 1));
		final SampleSet validate=new SampleSet(makeMatrix(2, 4, 0, 0.5f, 1));
		Trainer.initializeSampleCutoffs(data, validate);
		System.out.println("explicit_cutoffs\ttrainingCutoff="+Cell.cutoffForTraining+
			"\tevalCutoff="+Trainer.cutoffForEvaluation+"\ttrainPositive="+
			data.numPositive+"\tvalidatePositive="+validate.numPositive);
		assert(Cell.cutoffForTraining==2.5f) : Cell.cutoffForTraining;
		assert(Trainer.cutoffForEvaluation==3.5f) : Trainer.cutoffForEvaluation;
		assertCounts(data, 2, 1);
		assertCounts(validate, 1, 2);
	}

	/** Validation uses the training-derived global midpoint, not its own range. */
	private static void checkValidationRange(){
		Cell.useMidpoint=true;
		final SampleSet data=new SampleSet(makeMatrix(2, 4, 0, 1));
		final SampleSet validate=new SampleSet(makeMatrix(10, 12, 0, 1));
		Trainer.initializeSampleCutoffs(data, validate);
		System.out.println("validation_range\ttrainingMidpoint="+data.outputMidpoint()+
			"\tvalidationMidpoint="+validate.outputMidpoint()+"\tevalCutoff="+
			Trainer.cutoffForEvaluation+"\tvalidatePositive="+validate.numPositive);
		assert(data.outputMidpoint()==3f) : data.outputMidpoint();
		assert(validate.outputMidpoint()==11f) : validate.outputMidpoint();
		assert(Trainer.cutoffForEvaluation==3f) : Trainer.cutoffForEvaluation;
		assertCounts(data, 1, 1);
		assertCounts(validate, 2, 0);
	}

	/** Without training data, implicit midpoint is not derived from validation. */
	private static void checkDataNull(){
		Cell.useMidpoint=true;
		final SampleSet validate=new SampleSet(makeMatrix(10, 12, 0, 1));
		Trainer.initializeSampleCutoffs(null, validate);
		System.out.println("data_null\tevalCutoff="+Trainer.cutoffForEvaluation+
			"\tvalidationMidpoint="+validate.outputMidpoint()+"\tvalidatePositive="+
			validate.numPositive);
		assert(Trainer.cutoffForEvaluation==0.5f) : Trainer.cutoffForEvaluation;
		assert(Cell.cutoffForTraining==0.5f) : Cell.cutoffForTraining;
		assertCounts(validate, 2, 0);
	}

	/** Copies and Subsets preserve the materialized positivity cutoff. */
	private static void checkCopySubset(){
		final SampleSet data=new SampleSet(makeMatrix(2, 4, 0, 1));
		data.makeSamples(3f);
		final SampleSet copy=data.copy(2, 0);
		final ArrayList<Sample> list=new ArrayList<Sample>();
		for(Sample s : copy.samples){list.add(s);}
		final Subset subset=new Subset(list);
		System.out.println("copy_subset\tcopyPositive="+copy.numPositive+
			"\tcopyNegative="+copy.numNegative+"\tsubsetPositive="+
			subset.positive.length+"\tsubsetNegative="+subset.negative.length);
		assertCounts(copy, 1, 1);
		assert(subset.positive.length==1 && subset.negative.length==1) :
			subset.positive.length+", "+subset.negative.length;
	}

	/** Legacy makeSamples entrypoints keep the fixed 0.5 cutoff. */
	private static void checkLegacyEntrypoint(){
		final SampleSet all=new SampleSet(makeMatrix(2, 4, 0, 1));
		all.makeSamples();
		final SampleSet capped=new SampleSet(makeMatrix(2, 4, 0, 1));
		capped.makeSamples(1);
		System.out.println("legacy_entrypoint\tallPositive="+all.numPositive+
			"\tallNegative="+all.numNegative+"\tcappedPositive="+
			capped.numPositive+"\tcappedNegative="+capped.numNegative);
		assertCounts(all, 2, 0);
		assertCounts(capped, 1, 0);
	}

	/** Explicit training cutoff overrides midpoint while implicit eval uses midpoint. */
	private static void checkExplicitTrainingOnly(){
		Cell.useMidpoint=true;
		Cell.cutoffForTraining=2.5f;
		Cell.setCutoffForTraining=true;
		final SampleSet data=new SampleSet(makeMatrix(2, 4, 0, 0.5f, 1));
		final SampleSet validate=new SampleSet(makeMatrix(2, 4, 0, 0.5f, 1));
		Trainer.initializeSampleCutoffs(data, validate);
		System.out.println("explicit_training_only\ttrainingCutoff="+
			Cell.cutoffForTraining+"\tevalCutoff="+Trainer.cutoffForEvaluation+
			"\ttrainPositive="+data.numPositive+"\tvalidatePositive="+
			validate.numPositive);
		assert(Cell.cutoffForTraining==2.5f) : Cell.cutoffForTraining;
		assert(Trainer.cutoffForEvaluation==3f) : Trainer.cutoffForEvaluation;
		assertCounts(data, 2, 1);
		assertCounts(validate, 2, 1);
	}

	/** Implicit training uses midpoint while explicit eval cutoff overrides it. */
	private static void checkExplicitEvalOnly(){
		Cell.useMidpoint=true;
		Trainer.cutoffForEvaluation=3.5f;
		Trainer.setCutoffForEvaluation=true;
		final SampleSet data=new SampleSet(makeMatrix(2, 4, 0, 0.5f, 1));
		final SampleSet validate=new SampleSet(makeMatrix(2, 4, 0, 0.5f, 1));
		Trainer.initializeSampleCutoffs(data, validate);
		System.out.println("explicit_eval_only\ttrainingCutoff="+
			Cell.cutoffForTraining+"\tevalCutoff="+Trainer.cutoffForEvaluation+
			"\ttrainPositive="+data.numPositive+"\tvalidatePositive="+
			validate.numPositive);
		assert(Cell.cutoffForTraining==3f) : Cell.cutoffForTraining;
		assert(Trainer.cutoffForEvaluation==3.5f) : Trainer.cutoffForEvaluation;
		assertCounts(data, 2, 1);
		assertCounts(validate, 1, 2);
	}

	/** Restores default numeric cutoffs and their explicit-setting flags. */
	private static void resetGlobals(){
		Cell.cutoffForTraining=0.5f;
		Cell.setCutoffForTraining=false;
		Cell.useMidpoint=false;
		Trainer.cutoffForEvaluation=0.5f;
		Trainer.setCutoffForEvaluation=false;
	}

	/** Constructs a matrix whose raw goals are remapped into the target range. */
	private static Matrix makeMatrix(final float targetMin, final float targetMax,
			final float... rawGoals){
		Matrix.convertTo01=false;
		Matrix.targetOutputRangeMin=targetMin;
		Matrix.targetOutputRangeMax=targetMax;
		Matrix.setTargetOutputRangeMin=true;
		Matrix.setTargetOutputRangeMax=true;
		final Matrix m=new Matrix();
		m.columns=new ArrayList<String>();
		m.columns.add("input");
		m.columns.add("target");
		m.numInputs=1;
		m.numOutputs=1;
		m.inputs=new float[rawGoals.length][1];
		m.outputs=new float[rawGoals.length][1];
		m.weights=new float[rawGoals.length][1];
		for(int i=0; i<rawGoals.length; i++){
			m.inputs[i][0]=i;
			m.outputs[i][0]=rawGoals[i];
			m.weights[i][0]=1;
		}
		m.data=new float[][][]{m.inputs, m.outputs, m.weights};
		m.validLines=m.inputs.length;
		m.initializeRange();
		return m;
	}

	/** Counts the cached Sample.positive flags materialized by makeSamples. */
	private static int countPositive(final Sample[] samples){
		int count=0;
		for(Sample s : samples){
			if(s.positive){count++;}
		}
		return count;
	}

	/** Copies each expected output into result so ranking utilities can score the fixture. */
	private static void copyGoalToResult(final SampleSet set){
		for(Sample s : set.samples){
			s.result[0]=s.goal[0];
		}
	}

	/** Verifies stored class counts and the per-sample cached-positive total agree. */
	private static void assertCounts(final SampleSet set, final int positives, final int negatives){
		assert(set.numPositive==positives && set.numNegative==negatives) :
			set.numPositive+", "+set.numNegative;
		assert(countPositive(set.samples)==positives) : countPositive(set.samples);
	}
}
