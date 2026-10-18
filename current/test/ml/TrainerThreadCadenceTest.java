package test.ml;

import ml.CellNet;
import ml.JobData;
import ml.Matrix;
import ml.Sample;
import ml.SampleSet;
import ml.Subset;
import ml.Trainer;
import ml.TrainerThread;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

/**
 * SelectTrainingSubset cadence fixture that invokes the real TrainerThread method
 * without starting workers.
 *
 * @author Sucrose
 */
public final class TrainerThreadCadenceTest{

	/** Runs the retained cadence case and two sorting bypass controls. */
	public static void main(String[] args) throws Exception{
		if(args.length==2){
			run(args[0], Paths.get(args[1]));
			return;
		}
		assert(args.length==0) :
			"This fixture accepts either no args or a case name and artifact directory";
		run("retained-64", Paths.get("."));
		run("sort-disabled", Paths.get("."));
		run("sort-in-thread", Paths.get("."));
	}

	/** Runs one selectable case. */
	private static void run(final String name, final Path artifacts) throws Exception{
		if(name.equals("retained-64")){
			checkCadence(name, artifacts, true, false, expectedRetained64());
		}else if(name.equals("sort-disabled")){
			checkCadence(name, artifacts, false, false, expectedShuffleOnly());
		}else if(name.equals("sort-in-thread")){
			checkCadence(name, artifacts, true, true, expectedShuffleOnly());
		}else{
			throw new IllegalArgumentException("Unknown cadence case: "+name);
		}
		System.err.println("TrainerThreadCadenceTest PASS "+name);
	}

	/** Checks the recorded shuffle and sort calls for updates 1..128. */
	private static void checkCadence(final String name, final Path artifacts,
			final boolean sort, final boolean sortInThread, final String expected) throws Exception{
		final RecordingSubset subset=new RecordingSubset(samples());
		final TrainerHarness harness=trainer(subset, artifacts, sort, sortInThread);
		for(int epoch=1; epoch<=128; epoch++){
			subset.epoch=epoch;
			CURRENT_EPOCH.setInt(harness.trainer, epoch);
			SELECT_TRAINING_SUBSET.invoke(harness.trainer);
		}
		final String actual=subset.log.toString();
		assert(actual.equals(expected)) : "Cadence drifted for "+name+
			"\nExpected:\n"+expected+"Actual:\n"+actual;
	}

	/** Expected retained cadence when full sort uses the low 6 epoch bits. */
	private static String expectedRetained64(){
		final StringBuilder sb=new StringBuilder();
		for(int epoch=1; epoch<=128; epoch++){
			final int mod8=epoch&7;
			if(mod8==Trainer.SHUFFLEMOD){
				append(sb, epoch, "shuffle");
			}else if(epoch==5 || epoch==69){
				append(sb, epoch, "sortSamples");
			}else if(mod8==5){
				append(sb, epoch, "sortSamples2:1");
			}else if(mod8==1){
				append(sb, epoch, "sortSamples2:6x");
			}else if(mod8==3 || mod8==7){
				append(sb, epoch, "sortSamples2:3x");
			}
		}
		return sb.toString();
	}

	/** Expected control cadence when TrainerThread bypasses explicit sorting. */
	private static String expectedShuffleOnly(){
		final StringBuilder sb=new StringBuilder();
		for(int epoch=1; epoch<=128; epoch++){
			if((epoch&7)==Trainer.SHUFFLEMOD){append(sb, epoch, "shuffle");}
		}
		return sb.toString();
	}

	/** Appends one expected action row. */
	private static void append(final StringBuilder sb, final int epoch, final String action){
		sb.append(epoch).append('\t').append(action).append('\n');
	}

	/** Builds a TrainerThread around the recording subset and disables worker launch. */
	private static TrainerHarness trainer(final RecordingSubset subset, final Path artifacts,
			final boolean sort, final boolean sortInThread) throws Exception{
		final Trainer parent=parent(subset, artifacts, sort, sortInThread);
		final TrainerThread trainer=new TrainerThread(parent, network());
		DATA.set(trainer, parent.data);
		return new TrainerHarness(trainer);
	}

	/** Builds a real Trainer with one injected training subset. */
	private static Trainer parent(final RecordingSubset subset, final Path artifacts,
			final boolean sort, final boolean sortInThread) throws Exception{
		final Trainer parent=new Trainer(new String[]{
			"in="+inputPath(artifacts), "dims=1,1", "epochs=128", "t=1",
			"launchinthread=f", "setlock=f", "ordered=t", "shuffle2=t", "ow=t"});
		parent.jobsPerEpoch=JOBS_PER_EPOCH;
		parent.workerQueue=new ArrayBlockingQueue<JobData>(JOBS_PER_EPOCH);
		parent.data=sampleSet(subset);
		parent.launchInThread=false;
		parent.useSetLock=false;
		parent.orderedJobs=true;
		parent.sort=sort;
		parent.sortInThread=sortInThread;
		parent.shuffleSubset=true;
		parent.allowMTSort=true;
		parent.training=true;
		parent.maxEpochs=128;
		return parent;
	}

	/** Writes a tiny existing input beside the runner's other case artifacts. */
	private static Path inputPath(final Path artifacts) throws Exception{
		assert(Files.isDirectory(artifacts)) : "Missing artifact directory: "+artifacts;
		final Path input=artifacts.resolve("input.tsv");
		Files.write(input, "0\t0\n1\t1\n".getBytes(StandardCharsets.UTF_8));
		return input;
	}

	/** Builds a matrix-backed SampleSet and injects the recording Subset. */
	private static SampleSet sampleSet(final RecordingSubset subset) throws Exception{
		final Matrix matrix=new Matrix();
		matrix.numInputs=1;
		matrix.numOutputs=1;
		matrix.inputs=new float[subset.samples.length][1];
		matrix.outputs=new float[subset.samples.length][1];
		matrix.weights=new float[subset.samples.length][1];
		for(int i=0; i<subset.samples.length; i++){
			matrix.inputs[i][0]=i;
			matrix.outputs[i][0]=(i&1);
			matrix.weights[i][0]=1;
		}
		matrix.data=new float[][][]{matrix.inputs, matrix.outputs, matrix.weights};
		final SampleSet set=new SampleSet(matrix);
		set.makeSamples();
		SUBSETS.set(set, new Subset[]{subset});
		return set;
	}

	/** Builds stable samples for the recording subset. */
	private static ArrayList<Sample> samples(){
		final ArrayList<Sample> list=new ArrayList<Sample>();
		for(int i=0; i<8; i++){
			list.add(new Sample(new float[]{i}, new float[]{i&1}, 1, i, 0.5f));
		}
		return list;
	}

	/** Builds the one-input, one-output net copied into TrainerThread. */
	private static CellNet network(){
		return Trainer.randomNetwork(DIMS, 1L, 1f, 1f, Long.MAX_VALUE, 1,
			new ArrayList<String>(), true);
	}

	/** Looks up one private no-arg method. */
	private static Method method(final Class<?> clazz, final String name){
		try{
			final Method method=clazz.getDeclaredMethod(name);
			method.setAccessible(true);
			return method;
		}catch(NoSuchMethodException e){throw new AssertionError(e);}
	}

	/** Looks up one field on the supplied class. */
	private static Field field(final Class<?> clazz, final String name){
		try{
			final Field field=clazz.getDeclaredField(name);
			field.setAccessible(true);
			return field;
		}catch(NoSuchFieldException e){throw new AssertionError(e);}
	}

	/** Recording Subset that prevents sample reordering while preserving call shape. */
	private static final class RecordingSubset extends Subset{
		RecordingSubset(final ArrayList<Sample> samples){super(samples);}

		@Override
		protected void shuffle(){append(log, epoch, "shuffle");}

		@Override
		protected void sortSamples(final float fraction, final boolean allowMultithreadedSort){
			assert(fraction==1f) : "Unexpected full-sort fraction "+fraction;
			assert(allowMultithreadedSort) : "Expected multithreaded full sort";
			append(log, epoch, "sortSamples");
		}

		@Override
		protected void sortSamples2(final float fraction, final int useSamples,
				final boolean allowMultithreadedSort, final structures.FloatList pivots){
			assert(useSamples>0) : "Expected positive epoch sample count";
			assert(allowMultithreadedSort) : "Expected multithreaded partial sort";
			assert(pivots!=null) : "Expected reusable pivot buffer";
			append(log, epoch, action(fraction));
		}

		/** Labels the partial-sort branch selected by TrainerThread. */
		private static String action(final float fraction){
			if(fraction==1f){return "sortSamples2:1";}
			if(fraction>0.47f && fraction<0.49f){return "sortSamples2:6x";}
			if(fraction>0.23f && fraction<0.25f){return "sortSamples2:3x";}
			throw new AssertionError("Unexpected partial-sort fraction "+fraction);
		}

		final StringBuilder log=new StringBuilder();
		int epoch;
	}

	/** Real constructed TrainerThread with reflected data already injected. */
	private static final class TrainerHarness{
		TrainerHarness(final TrainerThread trainer_){trainer=trainer_;}
		final TrainerThread trainer;
	}

	private static final int JOBS_PER_EPOCH=4;
	private static final int[] DIMS={1, 1};
	private static final Method SELECT_TRAINING_SUBSET=
		method(TrainerThread.class, "selectTrainingSubset");
	private static final Field CURRENT_EPOCH=field(TrainerThread.class, "currentEpoch");
	private static final Field DATA=field(TrainerThread.class, "data");
	private static final Field SUBSETS=field(SampleSet.class, "subsets");

	private TrainerThreadCadenceTest(){}
}
