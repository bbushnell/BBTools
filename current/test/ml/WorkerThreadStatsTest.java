package test.ml;

import ml.JobData;
import ml.Sample;
import ml.WorkerThread;

import java.lang.reflect.Field;
import java.util.concurrent.ArrayBlockingQueue;

/**
 * Selectable fixture for WorkerThread's per-output confusion counters and
 * error accumulation without starting a worker thread.
 * @author Sucrose
 */
public class WorkerThreadStatsTest{

	public static void main(String[] args) throws Exception{
		if(args.length==1){
			run(args[0]);
			return;
		}
		assert(args.length==0) : "This fixture accepts only one optional case name";
		run("one_output_control");
		run("two_output_reproducer");
		run("four_output_all_counts");
	}

	/** Runs one selectable case; separate invocations let an external harness isolate failures. */
	private static void run(final String name) throws Exception{
		if(name.equals("one_output_control")){
			checkOneOutputControl();
		}else if(name.equals("two_output_reproducer")){
			checkTwoOutputReproducer();
		}else if(name.equals("four_output_all_counts")){
			checkFourOutputAllCounts();
		}else{
			throw new IllegalArgumentException("Unknown fixture mode: "+name);
		}
		System.err.println("WorkerThreadStatsTest PASS "+name);
	}

	/** Preserves the legacy one-output TP path while checking error sums. */
	private static void checkOneOutputControl() throws Exception{
		final WorkerThread worker=makeWorker();
		add(worker, sample(new float[]{1}, new float[]{1}, 0.25f, 0.75f));
		checkStats(worker, 1, 0, 0, 0, 0.25, 0.75);
	}

	/** Checks correct per-output counts for the row that exposed duplicated output-0 counting. */
	private static void checkTwoOutputReproducer() throws Exception{
		final WorkerThread worker=makeWorker();
		add(worker, sample(new float[]{1, 0}, new float[]{1, 1}, 1.5f, 2.5f));
		checkStats(worker, 1, 0, 1, 0, 1.5, 2.5);
	}

	/** Covers TP, TN, FP and FN with nonidentical predictions. */
	private static void checkFourOutputAllCounts() throws Exception{
		final WorkerThread worker=makeWorker();
		add(worker, sample(new float[]{1, 0, 0, 1}, new float[]{1, 0, 1, 0}, 3.25f, 4.5f));
		checkStats(worker, 1, 1, 1, 1, 3.25, 4.5);
	}

	/** Constructs a worker without starting its thread. */
	private static WorkerThread makeWorker(){
		return new WorkerThread(0, new ArrayBlockingQueue<JobData>(1), 0.5f);
	}

	/** Adds one sample through WorkerThread.addToStats(). */
	private static void add(final WorkerThread worker, final Sample sample){
		worker.addToStats(sample);
	}

	/** Builds a sample with precomputed prediction and error values. */
	private static Sample sample(final float[] goal, final float[] result,
			final float error, final float weightedError){
		final float[] in=new float[]{0};
		final Sample s=new Sample(in, goal, 1, 0, 0.5f);
		for(int i=0; i<result.length; i++){s.result[i]=result[i];}
		s.errorMagnitude=error;
		s.weightedErrorMagnitude=weightedError;
		return s;
	}

	/** Checks all counters and both error accumulators. */
	private static void checkStats(final WorkerThread worker, final int tp, final int tn,
			final int fp, final int fn, final double error, final double weightedError)
			throws Exception{
		check(getInt(worker, "tpSum")==tp, "tpSum", tp, getInt(worker, "tpSum"));
		check(getInt(worker, "tnSum")==tn, "tnSum", tn, getInt(worker, "tnSum"));
		check(getInt(worker, "fpSum")==fp, "fpSum", fp, getInt(worker, "fpSum"));
		check(getInt(worker, "fnSum")==fn, "fnSum", fn, getInt(worker, "fnSum"));
		check(getDouble(worker, "errorSum")==error,
			"errorSum", error, getDouble(worker, "errorSum"));
		check(getDouble(worker, "weightedErrorSum")==weightedError,
			"weightedErrorSum", weightedError, getDouble(worker, "weightedErrorSum"));
	}

	/** Throws a diagnostic assertion failure for a mismatched value. */
	private static void check(final boolean pass, final String name,
			final Object expected, final Object observed){
		if(!pass){throw new AssertionError(name+" expected "+expected+" observed "+observed);}
	}

	/** Reads a private integer counter from WorkerThread. */
	private static int getInt(final WorkerThread worker, final String name) throws Exception{
		final Field field=WorkerThread.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getInt(worker);
	}

	/** Reads a private double accumulator from WorkerThread. */
	private static double getDouble(final WorkerThread worker, final String name) throws Exception{
		final Field field=WorkerThread.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getDouble(worker);
	}
}
