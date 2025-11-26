package ml;

import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ArrayBlockingQueue;

import shared.Tools;

/**
 * Manages parallel neural network processing tasks, enabling distributed
 * computation across multiple threads with fine-grained job processing.
 * Orchestrates thread-safe sample processing, network computation, and
 * result collection during distributed machine learning training.
 *
 * @author Brian Bushnell
 * @date 2013
 */
class WorkerThread extends Thread implements Comparable<WorkerThread> {
	//Constructor
	/**
	 * Constructs a worker thread with specified configuration parameters.
	 * @param tid_ Thread identifier for tracking and sorting
	 * @param wq Job queue for receiving work assignments
	 * @param cutoffForEvaluation_ Threshold for binary classification evaluation
	 */
	WorkerThread(final int tid_, final ArrayBlockingQueue<JobData> wq,	/*final Object LOCK_, */final float cutoffForEvaluation_){
		tid=tid_;
		jobQueue=wq;
//		LOCK=LOCK_;
		cutoffForEvaluation=cutoffForEvaluation_;
	}

	//Called by start()
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

	/** Iterate through the lines */
	void processInner(){

		while(true) {
			//				synchronized(LOCK) {
			//				synchronized(this) {
			tprof.reset();
			//					if(net!=null) {net.clear();} //Not a good place since master could read it
			clearStats();
			tprof.log();//0
			//				}
			//				}
			
			
			final JobData job=getJob();
			tprof.log();//1: 47855878
			if(job==JobData.POISON) {break;}
			
			synchronized(this) {
				prepareForWork(job);
				assert(net!=null);
				tprof.log();//2: 209349

				if(job.sort){sortSamples(job);}
				tprof.log();//3: 0/1000000?
				
				synchronized(net) {
					int processed=processSamples(job.weightMult);
					tprof.log();//4: 7417842/8404496
					sendResults(processed, job);
					if(job.setLock!=null) {job.setLock.readLock().unlock();}
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

	/**
	 * Retrieves next job from the blocking queue, handling interruption and
	 * blocking until work becomes available.
	 * @return Next job to process, or JobData.POISON for thread termination
	 */
	JobData getJob() {
		JobData job=null;
		while(job==null) {
			try {
				job=jobQueue.take();//Could process any network here, with the same dimensions
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			//					System.err.println("Thread "+tid+" got "+token.seed+"; "+(token==poison)+", size="+workerQueue.size());
		}
		//System.err.println("W"+tid+" got permission for epoch e="+epoch+", ce="+currentEpoch);
		return job;
	}

	/**
	 * Configures thread state for job execution including network setup,
	 * sample assignment, and synchronization. Handles both copied and shared
	 * network configurations based on job requirements.
	 * @param job Job configuration containing network, samples, and parameters
	 */
	void prepareForWork(JobData job) {
		epoch=job.epoch;
		backprop=job.backprop;
		maxSamples=job.maxSamples;
		samples=job.list;
		this.job=job;
		if(job.setLock!=null) {job.setLock.readLock().lock();}
		assert(samples==null || maxSamples<=samples.size());
		if(maxSamples<1) {
			assert(false) : job.maxSamples+", "+samples.size()+", "+epoch;
			return;
		}
		
		
		if(job.doCopy) {//TODO: If job size is zero, a null net can be returned
//			assert(Trainer.copyNetInWorkerThread);//Not currently true for scanner threads
			assert(job.mutableNet==null);
			net=job.immutableNet.copy(false);//Works, but the reason is uncertain
			net.transpose();
		}else{
//			synchronized(net) {//Does not work; a new net is needed
//				assert(false);
//				net.setFrom(job.net, false);
//				net.clear();
//			}
			net=job.mutableNet;
			synchronized(net) {
				assert(Trainer.setNetInWorkerThread==(job.immutableNet!=null)) : Trainer.setNetInWorkerThread;
				if(job.immutableNet!=null) {
					net.setFrom(job.immutableNet, false);
				}
				net.clear();
				net.transpose();
			}
		}
		assert(net!=null);
		//			}
		//			net.setFrom(job.net, false);
	}

	/**
	 * Processes assigned samples through the neural network with thread-safe
	 * synchronization. Handles both array-based and list-based sample sets
	 * with configurable weight multiplication.
	 *
	 * @param weightMult Multiplier for sample weights during processing
	 * @return Number of samples successfully processed
	 */
	int processSamples(float weightMult) {
		//			System.err.println("processEpoch()");

		int samplesProcessed=0;
		//System.err.println("W"+tid+" starts processing samples for "+epoch);

//		for(int pos=tid; pos<maxSampleT && pos<samples.length; pos+=threads) {
//			final Sample s=samples[pos];
//			synchronized(s) {//Syncing here on LOCK solves nondeterminism...
//				processSample(s, backprop);
//				samplesProcessed++;
//			}
//		}
//
//		for(final Sample s : samples) {
//			synchronized(s) {//Syncing here on LOCK solves nondeterminism...
//				processSample(s, backprop);
//				samplesProcessed++;
//			}
//		}

		if(job.set!=null) {
			for(int i=job.jid; i<maxSamples; i+=job.jobsPerEpoch) {
				Sample s=job.set[i];
				synchronized(s) {
					processSample(s, backprop, weightMult);
					samplesProcessed++;
				}
			}
		}else {
			for(int i=0; i<maxSamples; i++) {
				Sample s=samples.get(i);
				synchronized(s) {
					processSample(s, backprop, weightMult);
					samplesProcessed++;
				}
			}
		}
		return samplesProcessed;
	}

	/**
	 * Packages processing results and performance metrics into JobResults
	 * and queues for collection by coordinator thread.
	 * @param samplesProcessed Number of samples processed in this job
	 * @param job Original job configuration for result correlation
	 */
	void sendResults(int samplesProcessed, JobData job) {
		assert(maxSamples==samplesProcessed || samples==null) : maxSamples+", "+samplesProcessed+", "+(samples==null ? job.set.length : samples.size());
		JobResults jr=new JobResults(maxSamples>0 ? net : null, epoch, samplesProcessed, tid, job.jid,
				errorSum, weightedErrorSum, tpSum, tnSum, fpSum, fnSum);
		net=null;//This is necessary
		job.jobResultsQueue.add(jr);
	}

	/**
	 * Process a sample.
	 * @param line sample number
	 */
	void processSample(final Sample sample, boolean backprop, float weightMult){
		sample.setEpoch(Tools.max(epoch, sample.epoch()));
		sample.setLastTID(tid);
		
		net.processSample(sample, backprop, weightMult);
		
		addToStats(sample);
		sample.setPivot();
		
		linesProcessedT++;
	}
	
	/**
	 * Sorts samples by class (positive/negative) and error magnitude, then
	 * interleaves them for balanced processing. Maintains separate positive
	 * and negative sample lists during sorting.
	 * @param job Job containing samples to sort
	 */
	void sortSamples(JobData job) {
//		if(true) {return;}
		//			positive.clear();
		//			negative.clear();
		assert(positive.size()==0);
		assert(negative.size()==0);
		
		for(Sample s : samples) {
//			if(s.epoch()<2) {s.setPivot();}
			s.setPivot();
			assert(s.checkPivot()) : job;
			if(s.positive) {positive.add(s);}
			else {negative.add(s);}
		}
		
		Collections.sort(positive);
		Collections.sort(negative);
//		assert(false) : positive.size()+", "+negative.size();

		final int limit=samples.size();
		samples.clear();
		int ppos=0, npos=0;
		while(samples.size()<limit) {
			if(npos<negative.size()) {
				samples.add(negative.get(npos));
				npos++;
			}
			if(ppos<positive.size()) {
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
	
	/**
	 * Accumulates classification performance metrics from processed sample
	 * including true/false positives/negatives and error magnitudes.
	 * @param s Sample containing prediction results and target goals
	 */
	void addToStats(Sample s) {
		for(int i=0; i<s.result.length; i++){
//			boolean goal=(s.goal[0]>=Trainer.booleanCutoffGoal);
			boolean goal=(s.goal[0]>=cutoffForEvaluation);
			boolean pred=(s.result[0]>=cutoffForEvaluation);
			tpSum+=(goal && pred) ? 1 : 0;
			tnSum+=(!goal && !pred) ? 1 : 0;
			fpSum+=(!goal && pred) ? 1 : 0;
			fnSum+=(goal && !pred) ? 1 : 0;
		}
		assert(s.errorMagnitude>=0 && s.weightedErrorMagnitude>=0) : s;
		errorSum+=s.errorMagnitude;
		weightedErrorSum+=s.weightedErrorMagnitude;
	}

	/** Resets accumulated performance statistics to zero for new job processing */
	synchronized private void clearStats() {
		errorSum=0;
		weightedErrorSum=0;
		tpSum=tnSum=fpSum=fnSum=0;
	}

	@Override
	public int compareTo(WorkerThread o) {
		return tid-o.tid;
	}
	
	/** Current job being processed by this thread */
	private JobData job;
	
	/** Blocking queue for receiving job assignments from coordinator */
	private final ArrayBlockingQueue<JobData> jobQueue;
	/** Threshold value for binary classification evaluation metrics */
	private final float cutoffForEvaluation;

	/** Accumulated raw error magnitude across processed samples */
	private double errorSum=0;
	/** Accumulated weighted error magnitude across processed samples */
	private double weightedErrorSum=0;
	/** Count of false negative predictions */
	/** Count of false positive predictions */
	/** Count of true negative predictions */
	/** Count of true positive predictions */
	private int tpSum=0, tnSum=0, fpSum=0, fnSum=0;

	//		private Sample[] samples;
	/** List of samples assigned to this thread for processing */
	private ArrayList<Sample> samples;
	/** Temporary storage for positive samples during sorting operations */
	private final ArrayList<Sample> positive=new ArrayList<Sample>();
	/** Temporary storage for negative samples during sorting operations */
	private final ArrayList<Sample> negative=new ArrayList<Sample>();
	/** Maximum number of samples to process in current job */
	private int maxSamples=0;
	/** Performance profiler for timing different processing phases */
	final Profiler tprof=new Profiler("W", 7);

	/** Number of reads processed by this thread */
	protected long linesProcessedT=0;

	/** Number of reads retained by this thread */
	protected long linesOutT=0;

	/** Flag indicating whether this thread encountered an error state */
	protected boolean errorStateT=false;

	/** True only if this thread has completed successfully */
	boolean success=false;

	/** Current training epoch number for this job */
	private int epoch=0;
	/** Whether to perform backpropagation during sample processing */
	private boolean backprop;

	/** Thread ID */
	final int tid;

	/** Neural network instance used by this thread for sample processing */
	private CellNet net;
	//		private final CellNet net;

//	@Deprecated
//	private final Object LOCK; //Only for testing synchronization
}