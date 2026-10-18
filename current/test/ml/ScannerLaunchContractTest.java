package test.ml;

import ml.CellNet;
import ml.JobData;
import ml.Matrix;
import ml.Sample;
import ml.SampleSet;
import ml.ScannerThread;
import ml.Seed;
import ml.Trainer;
import ml.TrainerThread;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

/**
 * Direct ScannerThread and TrainerThread launch-contract checks.
 *
 * <p>These cases use real Trainer, SampleSet, ScannerThread and TrainerThread
 * constructors. They invoke private launch methods directly and inspect only
 * queue shapes and lock ownership; they do not start workers or gather results.
 *
 * @author Sucrose
 */
public final class ScannerLaunchContractTest{

	/** Runs one selectable launch-contract case. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one case name");}
		final String name=args[0];
		final Path inputPath=inputPath(name);
		if(name.equals("scanner-list-direct-scarce")){
			checkScannerListDirect(inputPath, 3);
		}else if(name.equals("scanner-list-direct-zero")){
			checkScannerListDirect(inputPath, 0);
		}else if(name.equals("scanner-list-async-scarce")){
			checkScannerListAsync(inputPath, 3);
		}else if(name.equals("trainer-list-direct")){
			checkTrainerListDirect(inputPath, 4);
		}else if(name.equals("trainer-list-direct-scarce")){
			checkTrainerListDirectScarce(inputPath, 3);
		}else if(name.equals("trainer-list-async")){
			checkTrainerListAsync(inputPath, 4);
		}else if(name.equals("trainer-list-async-scarce")){
			checkTrainerListAsyncScarce(inputPath, 3);
		}else if(name.equals("scanner-set-direct")){
			checkScannerSetDirect(inputPath, 5);
		}else if(name.equals("trainer-set-direct")){
			checkTrainerSetDirect(inputPath, 5);
		}else if(name.equals("scanner-set-async")){
			checkScannerSetAsync(inputPath, 5);
		}else if(name.equals("trainer-set-async")){
			checkTrainerSetAsync(inputPath, 5);
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
		System.err.println("ScannerLaunchContractTest PASS "+name);
	}

	/** Checks ScannerThread direct list launch returns exactly the emitted job count. */
	private static void checkScannerListDirect(final Path inputPath, final int sampleCount)
			throws Exception{
		final ScannerHarness h=scanner(false, false, inputPath, sampleCount);
		final int jobs=invokeScannerLaunch(h.scanner, h.net, h.samples, sampleCount);
		assertEquals(sampleCount, jobs, "Scanner direct list return");
		assertEquals(sampleCount, workerQueue(h.scanner).size(), "Scanner direct worker jobs");
		assertListJobs(workerQueue(h.scanner), sampleCount, h.samples, false);
	}

	/** Checks ScannerThread async list launch delegates one aggregate before conversion. */
	private static void checkScannerListAsync(final Path inputPath, final int sampleCount)
			throws Exception{
		final ScannerHarness h=scanner(true, false, inputPath, sampleCount);
		final int jobs=invokeScannerLaunch(h.scanner, h.net, h.samples, sampleCount);
		final ArrayBlockingQueue<JobData> launchQueue=launchQueue(h.scanner);
		assertEquals(1, launchQueue.size(), "Scanner async aggregate count before conversion");
		assertEquals(0, workerQueue(h.scanner).size(), "Scanner async worker count before conversion");
		assertEquals(sampleCount, jobs, "Scanner async list return");
		final int converted=invokeScannerInner(h.scanner, only(launchQueue, h.samples, sampleCount));
		assertEquals(sampleCount, converted, "Scanner async manual conversion return");
		assertListJobs(workerQueue(h.scanner), sampleCount, h.samples, false);
	}

	/** Checks TrainerThread direct list launch preserves jobsPerEpoch. */
	private static void checkTrainerListDirect(final Path inputPath, final int sampleCount)
			throws Exception{
		final TrainerHarness h=trainer(false, false, inputPath, sampleCount);
		final int jobs=invokeTrainerLaunch(h.trainer, h.net, h.samples, sampleCount);
		assertEquals(JOBS_PER_EPOCH, jobs, "Trainer direct list return");
		assertEquals(JOBS_PER_EPOCH, workerQueue(h.trainer).size(), "Trainer direct worker jobs");
		assertListJobs(workerQueue(h.trainer), JOBS_PER_EPOCH, h.samples, true);
	}

	/** Checks TrainerThread direct scarce list launch returns exactly emitted jobs. */
	private static void checkTrainerListDirectScarce(final Path inputPath, final int sampleCount)
			throws Exception{
		final TrainerHarness h=trainer(false, false, inputPath, sampleCount);
		final int jobs=invokeTrainerLaunch(h.trainer, h.net, h.samples, sampleCount);
		assertEquals(sampleCount, jobs, "Trainer direct scarce list return");
		assertEquals(sampleCount, workerQueue(h.trainer).size(), "Trainer direct scarce worker jobs");
		assertListJobs(workerQueue(h.trainer), sampleCount, h.samples, true);
	}

	/** Checks TrainerThread async list launch delegates one aggregate before conversion. */
	private static void checkTrainerListAsync(final Path inputPath, final int sampleCount)
			throws Exception{
		final TrainerHarness h=trainer(true, false, inputPath, sampleCount);
		final int jobs=invokeTrainerLaunch(h.trainer, h.net, h.samples, sampleCount);
		final ArrayBlockingQueue<JobData> launchQueue=launchQueue(h.trainer);
		assertEquals(1, launchQueue.size(), "Trainer async aggregate count before conversion");
		assertEquals(0, workerQueue(h.trainer).size(), "Trainer async worker count before conversion");
		assertEquals(JOBS_PER_EPOCH, jobs, "Trainer async list return");
		final int converted=invokeTrainerInner(h.trainer, only(launchQueue, h.samples, sampleCount));
		assertEquals(JOBS_PER_EPOCH, converted, "Trainer async manual conversion return");
		assertListJobs(workerQueue(h.trainer), JOBS_PER_EPOCH, h.samples, true);
	}

	/** Checks TrainerThread async scarce list launch delegates an exact aggregate. */
	private static void checkTrainerListAsyncScarce(final Path inputPath, final int sampleCount)
			throws Exception{
		final TrainerHarness h=trainer(true, false, inputPath, sampleCount);
		final int jobs=invokeTrainerLaunch(h.trainer, h.net, h.samples, sampleCount);
		final ArrayBlockingQueue<JobData> launchQueue=launchQueue(h.trainer);
		assertEquals(1, launchQueue.size(), "Trainer async scarce aggregate count");
		assertEquals(0, workerQueue(h.trainer).size(), "Trainer async scarce worker count");
		assertEquals(sampleCount, jobs, "Trainer async scarce list return");
		final int converted=invokeTrainerInner(h.trainer, only(launchQueue, h.samples, sampleCount));
		assertEquals(sampleCount, converted, "Trainer async scarce manual conversion return");
		assertListJobs(workerQueue(h.trainer), sampleCount, h.samples, true);
	}

	/** Checks ScannerThread direct set launch converts under the caller's write lock. */
	private static void checkScannerSetDirect(final Path inputPath, final int sampleCount)
			throws Exception{
		final ScannerHarness h=scanner(false, true, inputPath, sampleCount);
		checkSetLaunch(h.scanner, h.net, h.samples, false);
	}

	/** Checks TrainerThread direct set launch converts under the caller's write lock. */
	private static void checkTrainerSetDirect(final Path inputPath, final int sampleCount)
			throws Exception{
		final TrainerHarness h=trainer(false, true, inputPath, sampleCount);
		checkSetLaunch(h.trainer, h.net, h.samples, false);
	}

	/** Checks ScannerThread async set launch converts on the caller thread, not LaunchThread. */
	private static void checkScannerSetAsync(final Path inputPath, final int sampleCount)
			throws Exception{
		final ScannerHarness h=scanner(true, true, inputPath, sampleCount);
		checkSetLaunch(h.scanner, h.net, h.samples, true);
	}

	/** Checks TrainerThread async set launch converts on the caller thread, not LaunchThread. */
	private static void checkTrainerSetAsync(final Path inputPath, final int sampleCount)
			throws Exception{
		final TrainerHarness h=trainer(true, true, inputPath, sampleCount);
		checkSetLaunch(h.trainer, h.net, h.samples, true);
	}

	/** Invokes one set-backed launch while the owner thread holds the write lock. */
	private static void checkSetLaunch(final Object target, final CellNet net, final Sample[] samples,
			final boolean async) throws Exception{
		final ReentrantReadWriteLock setLock=setLock(target);
		setLock.writeLock().lock();
		try{
			final int jobs=(target instanceof ScannerThread ?
					invokeScannerLaunch((ScannerThread)target, net, samples, samples.length) :
					invokeTrainerLaunch((TrainerThread)target, net, samples, samples.length));
			final ArrayBlockingQueue<JobData> launchQueue=launchQueue(target);
			assertEquals(JOBS_PER_EPOCH, jobs, "Set-backed return count");
			if(async){assertEquals(0, launchQueue.size(), "Async set-backed aggregate count");}
			else{assertNull(launchQueue, "Direct set-backed launch queue");}
			assertEquals(0, setLock.getWriteHoldCount(), "Write lock holds after set-backed launch");
			assertEquals(1, setLock.getReadHoldCount(), "Read lock holds after set-backed launch");
			assertSetJobs(workerQueue(target), samples, setLock);
		}finally{
			while(setLock.getReadHoldCount()>0){setLock.readLock().unlock();}
			while(setLock.getWriteHoldCount()>0){setLock.writeLock().unlock();}
		}
	}

	/** Builds a real ScannerThread around a hand-built training set. */
	private static ScannerHarness scanner(final boolean launchInThread, final boolean useSetLock,
			final Path inputPath, final int sampleCount) throws Exception{
		final Trainer parent=parent(launchInThread, useSetLock, inputPath, sampleCount);
		final CellNet net=network();
		final ArrayBlockingQueue<ArrayList<Seed>> returnQueue=
				new ArrayBlockingQueue<ArrayList<Seed>>(1);
		final ScannerThread scanner=new ScannerThread(parent, DIMS, null, null, 1, 1, 1,
				sampleCount, 1L, returnQueue);
		return new ScannerHarness(scanner, net, parent.data.samples);
	}

	/** Builds a real TrainerThread around a hand-built training set. */
	private static TrainerHarness trainer(final boolean launchInThread, final boolean useSetLock,
			final Path inputPath, final int sampleCount) throws Exception{
		final Trainer parent=parent(launchInThread, useSetLock, inputPath, sampleCount);
		final CellNet net=network();
		final TrainerThread trainer=new TrainerThread(parent, net);
		return new TrainerHarness(trainer, net, parent.data.samples);
	}

	/** Builds a real Trainer with explicit retained input, thread, epoch and launch settings. */
	private static Trainer parent(final boolean launchInThread, final boolean useSetLock,
			final Path inputPath, final int sampleCount) throws Exception{
		final Trainer parent=new Trainer(new String[]{
				"in="+inputPath, "dims=1,1", "epochs=1", "t=1",
				"launchinthread="+launchInThread, "setlock="+useSetLock,
				"ordered=t", "shuffle2=f", "ow=t"});
		parent.jobsPerEpoch=JOBS_PER_EPOCH;
		parent.workerQueue=new ArrayBlockingQueue<JobData>(JOBS_PER_EPOCH);
		parent.data=sampleSet(sampleCount);
		parent.launchInThread=launchInThread;
		parent.useSetLock=useSetLock;
		parent.orderedJobs=true;
		parent.sort=false;
		parent.sortInThread=false;
		parent.shuffleSubset=false;
		parent.training=true;
		parent.maxEpochs=1;
		return parent;
	}

	/** Returns a tiny existing input path retained for Trainer's constructor checks. */
	private static Path inputPath(final String name) throws Exception{
		final Path root=Files.createTempDirectory("scanner-launch-contract-"+name+"-");
		final Path input=root.resolve("input.tsv");
		Files.write(input, "0\t0\n1\t1\n".getBytes(StandardCharsets.UTF_8));
		return input;
	}

	/** Builds a matrix-backed SampleSet with one input, one output and one weight. */
	private static SampleSet sampleSet(final int count){
		final Matrix matrix=new Matrix();
		matrix.numInputs=1;
		matrix.numOutputs=1;
		matrix.inputs=new float[count][1];
		matrix.outputs=new float[count][1];
		matrix.weights=new float[count][1];
		for(int i=0; i<count; i++){
			matrix.inputs[i][0]=i;
			matrix.outputs[i][0]=(i&1);
			matrix.weights[i][0]=1;
		}
		matrix.data=new float[][][]{matrix.inputs, matrix.outputs, matrix.weights};
		final SampleSet set=new SampleSet(matrix);
		set.makeSamples();
		return set;
	}

	/** Builds the one-input, one-output net copied into worker jobs. */
	private static CellNet network(){
		return Trainer.randomNetwork(DIMS, 1L, 1f, 1f, Long.MAX_VALUE, 1,
				new ArrayList<String>(), true);
	}

	/** Invokes ScannerThread.launchJobs. */
	private static int invokeScannerLaunch(final ScannerThread scanner, final CellNet net,
			final Sample[] samples, final int sampleCount) throws Exception{
		return ((Integer)SCANNER_LAUNCH.invoke(scanner, net, samples, sampleCount,
				false, 1f, false)).intValue();
	}

	/** Invokes TrainerThread.launchJobs. */
	private static int invokeTrainerLaunch(final TrainerThread trainer, final CellNet net,
			final Sample[] samples, final int sampleCount) throws Exception{
		return ((Integer)TRAINER_LAUNCH.invoke(trainer, net, samples, sampleCount,
				false, 1f, false)).intValue();
	}

	/** Manually converts one ScannerThread aggregate launch request. */
	private static int invokeScannerInner(final ScannerThread scanner, final JobData job)
			throws Exception{
		return ((Integer)SCANNER_INNER.invoke(scanner, job.immutableNet, job.set, job.maxSamples,
				job.epoch, job.alpha, job.backprop, job.weightMult, job.sort)).intValue();
	}

	/** Manually converts one TrainerThread aggregate launch request. */
	private static int invokeTrainerInner(final TrainerThread trainer, final JobData job)
			throws Exception{
		return ((Integer)TRAINER_INNER.invoke(trainer, job.immutableNet, job.set, job.maxSamples,
				job.epoch, job.alpha, job.backprop, job.weightMult, job.sort)).intValue();
	}

	/** Requires an aggregate queue to contain exactly one job and then drains it. */
	private static JobData only(final ArrayBlockingQueue<JobData> queue, final Sample[] samples,
			final int sampleCount){
		assertEquals(1, queue.size(), "Aggregate queue size");
		final JobData job=queue.poll();
		assertNotNull(job, "Aggregate job");
		assertEquals(0, queue.size(), "Aggregate queue after drain");
		assertNull(job.list, "Aggregate list source");
		assertSame(samples, job.set, "Aggregate set source");
		assertEquals(sampleCount, job.maxSamples, "Aggregate maxSamples");
		assertEquals(0, job.jid, "Aggregate jid");
		assertEquals(0, job.jobsPerEpoch, "Aggregate jobsPerEpoch");
		assertNull(job.setLock, "List aggregate setLock");
		return job;
	}

	/** Drains and verifies list-backed jobs while avoiding unsupported empty jobs. */
	private static void assertListJobs(final ArrayBlockingQueue<JobData> queue,
			final int expected, final Sample[] samples, final boolean trainer){
		for(int i=0; i<expected; i++){
			final JobData job=queue.poll();
			assertNotNull(job, "List job");
			assertNotNull(job.list, "List job samples");
			assertEquals(i, job.jid, "List jid");
			assertNull(job.set, "List job set");
			assertNull(job.setLock, "List job setLock");
			final int expectedSize=(samples.length-i+JOBS_PER_EPOCH-1)/JOBS_PER_EPOCH;
			assertEquals(expectedSize, job.list.size(), "List job sample count");
			assertEquals(expectedSize, job.maxSamples, "List maxSamples");
			for(int sampleIndex=i, listIndex=0; sampleIndex<samples.length;
					sampleIndex+=JOBS_PER_EPOCH, listIndex++){
				assertSame(samples[sampleIndex], job.list.get(listIndex), "List job sample");
			}
			assertEquals(trainer ? JOBS_PER_EPOCH : 0, job.jobsPerEpoch, "List jobsPerEpoch");
		}
		assertEquals(0, queue.size(), "List worker queue after drain");
	}

	/** Drains and verifies set-backed jobs. */
	private static void assertSetJobs(final ArrayBlockingQueue<JobData> queue,
			final Sample[] samples, final ReentrantReadWriteLock setLock){
		assertEquals(JOBS_PER_EPOCH, queue.size(), "Set worker jobs");
		for(int i=0; i<JOBS_PER_EPOCH; i++){
			final JobData job=queue.poll();
			assertNotNull(job, "Set job");
			assertEquals(i, job.jid, "Set jid");
			assertEquals(JOBS_PER_EPOCH, job.jobsPerEpoch, "Set jobsPerEpoch");
			assertEquals(samples.length, job.maxSamples, "Set maxSamples");
			assertNull(job.list, "Set job list");
			assertSame(samples, job.set, "Set job samples");
			assertSame(setLock, job.setLock, "Set job lock");
		}
		assertEquals(0, queue.size(), "Set worker queue after drain");
	}

	/** Returns the workerQueue field from either launch owner. */
	@SuppressWarnings("unchecked")
	private static ArrayBlockingQueue<JobData> workerQueue(final Object owner) throws Exception{
		return (ArrayBlockingQueue<JobData>)field(owner, "workerQueue").get(owner);
	}

	/** Returns the launchQueue field from either launch owner. */
	@SuppressWarnings("unchecked")
	private static ArrayBlockingQueue<JobData> launchQueue(final Object owner) throws Exception{
		return (ArrayBlockingQueue<JobData>)field(owner, "launchQueue").get(owner);
	}

	/** Returns the setLock field from either launch owner. */
	private static ReentrantReadWriteLock setLock(final Object owner) throws Exception{
		return (ReentrantReadWriteLock)field(owner, "setLock").get(owner);
	}

	/** Looks up one private launch method. */
	private static Method method(final Class<?> clazz, final String name){
		try{
			final Method method=clazz.getDeclaredMethod(name, CellNet.class, Sample[].class,
					int.class, boolean.class, float.class, boolean.class);
			method.setAccessible(true);
			return method;
		}catch(NoSuchMethodException e){throw new AssertionError(e);}
	}

	/** Looks up one private launch-inner method. */
	private static Method innerMethod(final Class<?> clazz, final String name){
		try{
			final Method method=clazz.getDeclaredMethod(name, CellNet.class, Sample[].class,
					int.class, int.class, double.class, boolean.class, float.class, boolean.class);
			method.setAccessible(true);
			return method;
		}catch(NoSuchMethodException e){throw new AssertionError(e);}
	}

	/** Looks up one private field. */
	private static Field field(final Object owner, final String name) throws Exception{
		final Field field=owner.getClass().getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	private static void assertEquals(final int expected, final int observed, final String label){
		if(expected!=observed){throw new AssertionError(label+" expected "+expected+
				" but observed "+observed);}
	}

	private static void assertSame(final Object expected, final Object observed, final String label){
		if(expected!=observed){throw new AssertionError(label+" expected same object");}
	}

	private static void assertNull(final Object observed, final String label){
		if(observed!=null){throw new AssertionError(label+" expected null but observed "+observed);}
	}

	private static void assertNotNull(final Object observed, final String label){
		if(observed==null){throw new AssertionError(label+" expected nonnull");}
	}

	/** Real constructed ScannerThread with its launch inputs. */
	private static final class ScannerHarness{
		ScannerHarness(final ScannerThread scanner_, final CellNet net_, final Sample[] samples_){
			scanner=scanner_;
			net=net_;
			samples=samples_;
		}
		final ScannerThread scanner;
		final CellNet net;
		final Sample[] samples;
	}

	/** Real constructed TrainerThread with its launch inputs. */
	private static final class TrainerHarness{
		TrainerHarness(final TrainerThread trainer_, final CellNet net_, final Sample[] samples_){
			trainer=trainer_;
			net=net_;
			samples=samples_;
		}
		final TrainerThread trainer;
		final CellNet net;
		final Sample[] samples;
	}

	private static final int JOBS_PER_EPOCH=4;
	private static final int[] DIMS={1, 1};
	private static final Method SCANNER_LAUNCH=method(ScannerThread.class, "launchJobs");
	private static final Method TRAINER_LAUNCH=method(TrainerThread.class, "launchJobs");
	private static final Method SCANNER_INNER=innerMethod(ScannerThread.class, "launchJobsInner");
	private static final Method TRAINER_INNER=innerMethod(TrainerThread.class, "launchJobsInner");

	private ScannerLaunchContractTest(){}

}
