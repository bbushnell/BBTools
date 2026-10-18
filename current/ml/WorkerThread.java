package ml;

import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Tools;

/** Consumes {@link JobData} from a shared queue and returns per-job training or
 * validation results to that job's result queue. */
public class WorkerThread extends Thread implements Comparable<WorkerThread>{

	/** Constructs a worker bound to a job queue and evaluation cutoff.
	 * @param tid_ deterministic worker id used for ordering and reporting
	 * @param wq queue supplying work and poison jobs
	 * @param cutoffForEvaluation_ threshold used for per-output confusion counts */
	public WorkerThread(final int tid_, final ArrayBlockingQueue<JobData> wq,	/*final Object LOCK_, */final float cutoffForEvaluation_){
		tid=tid_;
		jobQueue=wq;
//		LOCK=LOCK_;
		cutoffForEvaluation=cutoffForEvaluation_;
	}

	/** Processes jobs until a poison job is received. */
	@Override
	public void run(){
		//Do anything necessary prior to processing

		//Process the samples
		processInner();

		//Do anything necessary after processing

		//Indicate successful exit status
		success=true;
		//			System.err.println("Worker Done.");
	}

	/** Takes jobs from the queue, prepares their network state, processes their
	 * samples, and sends one {@link JobResults} for each non-poison job. */
	void processInner(){

		while(true){
			//				synchronized(LOCK){
			//				synchronized(this){
			tprof.reset();
			//					if(net!=null){net.clear();} //Not a good place since master could read it
			clearStats();
			tprof.log();//0
			//				}
			//				}
			
			
			final JobData job=getJob();
			tprof.log();//1: 47855878
			if(job==JobData.POISON){break;}
			
			synchronized(this){
				prepareForWork(job);
				assert(net!=null);
				tprof.log();//2: 209349

				if(job.sort){sortSamples(job);}
				tprof.log();//3: 0/1000000?
				
				synchronized(net){
					int processed=processSamples(job.weightMult);
					tprof.log();//4: 7417842/8404496
					sendResults(processed, job);
					if(job.setLock!=null){job.setLock.readLock().unlock();}
					this.job=null;
					tprof.log();//4: 45718/39325
				}
			}
		}
		//			System.err.println("Worker "+tid+" finished");
	}
	

//	SITT:
//	W	464	47855878	209349	?	7417842	45718	0
//	
//	SITF:
//	W	473	44546439	212696	?	8404496	39325	0

	/** Blocks until a job is available, reporting interruptions and retrying.
	 * @return the next queued job, including poison jobs */
	JobData getJob(){
		JobData job=null;
		while(job==null){
			try{
				job=jobQueue.take();//Could process any network here, with the same dimensions
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			//					System.err.println("Thread "+tid+" got "+token.seed+"; "+(token==poison)+", size="+workerQueue.size());
		}
		//System.err.println("W"+tid+" got permission for epoch e="+epoch+", ce="+currentEpoch);
		return job;
	}

	/** Copies or refreshes the job network and records sample-processing state.
	 * Jobs with setLock take their read lock here; the normal processing path
	 * releases that lock after sending results.
	 * @param job next non-poison job to prepare */
	void prepareForWork(JobData job){
		epoch=job.epoch;
		backprop=job.backprop;
		maxSamples=job.maxSamples;
		samples=job.list;
		this.job=job;
		if(job.setLock!=null){job.setLock.readLock().lock();}
		assert(samples==null || maxSamples<=samples.size());
		if(maxSamples<1){
			assert(false) : job.maxSamples+", "+samples.size()+", "+epoch;
			return;
		}
		
		
		if(job.doCopy){//TODO: If job size is zero, a null net can be returned
//			assert(Trainer.copyNetInWorkerThread);//Not currently true for scanner threads
			assert(job.mutableNet==null);
			net=job.immutableNet.copy(false);//Works, but the reason is uncertain
			net.transpose();
		}else{
//			synchronized(net){//Does not work; a new net is needed
//				assert(false);
//				net.setFrom(job.net, false);
//				net.clear();
//			}
			net=job.mutableNet;
			synchronized(net){
				assert(Trainer.setNetInWorkerThread==(job.immutableNet!=null)) : Trainer.setNetInWorkerThread;
				if(job.immutableNet!=null){
					net.setFrom(job.immutableNet, false);
				}
				net.clear();
				net.transpose();
			}
		}
		assert(net!=null);
		net.simdFF=job.simdFF; //Per-job SIMD feed-forward flag
		//Input-noise RNG (rand.FastRandomXoshiro) keyed off the NETWORK seed (net.seed), decorrelated per
		//epoch+job so parallel jobs don't draw identical noise; created only when noise is enabled.
		net.noiseRandy=(CellNet.inputNoise>0 ? shared.Shared.random(net.seed^(epoch*0x9E3779B97F4A7C15L)^job.jid) : null);
		//			}
		//			net.setFrom(job.net, false);
	}

	/** Processes each sample assigned to the current job.
	 * @param weightMult multiplier passed to the network's weighted error update
	 * @return number of synchronized sample visits */
	int processSamples(float weightMult){
		//			System.err.println("processEpoch()");

		int samplesProcessed=0;
		//System.err.println("W"+tid+" starts processing samples for "+epoch);

//		for(int pos=tid; pos<maxSampleT && pos<samples.length; pos+=threads){
//			final Sample s=samples[pos];
//			synchronized(s){//Syncing here on LOCK solves nondeterminism...
//				processSample(s, backprop);
//				samplesProcessed++;
//			}
//		}
//
//		for(final Sample s : samples){
//			synchronized(s){//Syncing here on LOCK solves nondeterminism...
//				processSample(s, backprop);
//				samplesProcessed++;
//			}
//		}

		if(job.set!=null){
			for(int i=job.jid; i<maxSamples; i+=job.jobsPerEpoch){
				Sample s=job.set[i];
				synchronized(s){
					processSample(s, backprop, weightMult);
					samplesProcessed++;
				}
			}
		}else{
			for(int i=0; i<maxSamples; i++){
				Sample s=samples.get(i);
				synchronized(s){
					processSample(s, backprop, weightMult);
					samplesProcessed++;
				}
			}
		}
		return samplesProcessed;
	}

	/** Builds a result from this job's processed network and accumulated
	 * statistics, nulls this worker's net, and enqueues the result.
	 * @param samplesProcessed number of samples processed for this job
	 * @param job finished job receiving the result */
	void sendResults(int samplesProcessed, JobData job){
		assert(maxSamples==samplesProcessed || samples==null) : maxSamples+", "+samplesProcessed+", "+(samples==null ? job.set.length : samples.size());
		JobResults jr=new JobResults(maxSamples>0 ? net : null, epoch, samplesProcessed, tid, job.jid,
				errorSum, weightedErrorSum, tpSum, tnSum, fpSum, fnSum);
		net=null;//This is necessary
		job.jobResultsQueue.add(jr);
	}

	/** Marks the sample's epoch and last-worker id, processes it through the
	 * job network, then updates its counters and pivot.
	 *
	 * @param sample training sample to process
	 * @param backprop whether to run network backpropagation
	 * @param weightMult weight multiplier for this job */
	void processSample(final Sample sample, boolean backprop, float weightMult){
		sample.setEpoch(Tools.max(epoch, sample.epoch()));
		sample.setLastTID(tid);
		
		net.processSample(sample, backprop, weightMult);
		
		addToStats(sample);
		sample.setPivot();
		
		linesProcessedT++;
	}
	
	/** Alternates negative and positive samples after sorting each group by pivot.
	 * @param job current job, reported if a pivot check fails */
	void sortSamples(JobData job){
//		if(true){return;}
		//			positive.clear();
		//			negative.clear();
		assert(positive.size()==0);
		assert(negative.size()==0);
		
		for(Sample s : samples){
//			if(s.epoch()<2){s.setPivot();}
			s.setPivot();
			assert(s.checkPivot()) : job;
			if(s.positive){positive.add(s);}
			else {negative.add(s);}
		}
		
		Collections.sort(positive);
		Collections.sort(negative);
//		assert(false) : positive.size()+", "+negative.size();

		final int limit=samples.size();
		samples.clear();
		int ppos=0, npos=0;
		while(samples.size()<limit){
			if(npos<negative.size()){
				samples.add(negative.get(npos));
				npos++;
			}
			if(ppos<positive.size()){
				samples.add(positive.get(ppos));
				ppos++;
			}
		}
		assert(limit==samples.size());
		assert(limit==positive.size()+negative.size());
		assert(ppos==positive.size());
		assert(npos==negative.size());
		positive.clear();
		negative.clear();//Avoids dangling references
	}
	
	/** Counts each output's thresholded prediction and accumulates this
	 * sample's error totals once.
	 * @param s processed sample with current result and error values */
	public void addToStats(Sample s){
		for(int i=0; i<s.result.length; i++){
//			boolean goal=(s.goal[i]>=Trainer.booleanCutoffGoal);
			boolean goal=(s.goal[i]>=cutoffForEvaluation);
			boolean pred=(s.result[i]>=cutoffForEvaluation);
			tpSum+=(goal && pred) ? 1 : 0;
			tnSum+=(!goal && !pred) ? 1 : 0;
			fpSum+=(!goal && pred) ? 1 : 0;
			fnSum+=(goal && !pred) ? 1 : 0;
		}
		assert(s.errorMagnitude>=0 && s.weightedErrorMagnitude>=0) : s;
		errorSum+=s.errorMagnitude;
		weightedErrorSum+=s.weightedErrorMagnitude;
	}

	/** Clears all per-job error and confusion accumulators. */
	synchronized private void clearStats(){
		errorSum=0;
		weightedErrorSum=0;
		tpSum=tnSum=fpSum=fnSum=0;
	}

	/** Compares worker threads by thread id for deterministic ordering. */
	@Override
	public int compareTo(WorkerThread o){
		return tid-o.tid;
	}
	
	/** Current job being processed, retained for set-backed iteration metadata. */
	private JobData job;
	
	/** Queue shared with the master thread for work dispatch. */
	private final ArrayBlockingQueue<JobData> jobQueue;
	/** Threshold that divides positive from negative goals and predictions. */
	private final float cutoffForEvaluation;

	/** Sum of raw sample errors for the current job. */
	private double errorSum=0;
	/** Sum of weighted sample errors for the current job. */
	private double weightedErrorSum=0;
	/** Current job's true-positive, true-negative, false-positive, and
	 * false-negative output counts. */
	private int tpSum=0, tnSum=0, fpSum=0, fnSum=0;

	//		private Sample[] samples;
	/** List-backed job samples; null for set-backed jobs. */
	private ArrayList<Sample> samples;
	/** Scratch list for positive samples while sorting. */
	private final ArrayList<Sample> positive=new ArrayList<Sample>();
	/** Scratch list for negative samples while sorting. */
	private final ArrayList<Sample> negative=new ArrayList<Sample>();
	/** List sample count or exclusive upper bound for strided set jobs. */
	private int maxSamples=0;
	/** Timing checkpoints for historical worker profiling. */
	final Profiler tprof=new Profiler("W", 7);

	/** Number of samples processed by this worker. */
	protected long linesProcessedT=0;

	/** Reserved for symmetry with other BBTools workers. */
	protected long linesOutT=0;

	/** Reserved per-thread error flag. */
	protected boolean errorStateT=false;

	/** True after the worker exits through the poison-job path. */
	boolean success=false;

	/** Epoch attached to the current job. */
	private int epoch=0;
	/** Whether current samples run network backpropagation. */
	private boolean backprop;

	/** Stable worker id. */
	final int tid;

	/** Current job network, nulled before its result is enqueued. */
	private CellNet net;
	//		private final CellNet net;

//	@Deprecated
//	private final Object LOCK; //Only for testing synchronization
}
