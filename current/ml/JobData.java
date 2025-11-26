package ml;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

/**
 * Encapsulates configuration and state for parallel ML job processing with thread safety.
 * Thread-safe container for neural network training jobs, managing sample subsets,
 * network references, and job parameters during distributed learning.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class JobData {
	
	/**
	 * Constructs a JobData instance with complete training configuration parameters.
	 * Initializes job settings including network references, queue connections, and
	 * processing parameters.
	 *
	 * @param net_ The immutable neural network to use for this job
	 * @param jrq_ Queue for submitting job results after processing
	 * @param epoch_ Current training epoch number
	 * @param maxSamples_ Maximum number of samples to process in this job
	 * @param alpha_ Learning rate for gradient updates
	 * @param backprop_ Whether to perform backpropagation during training
	 * @param weightMult_ Weight multiplication factor for regularization
	 * @param sort_ Whether to sort samples before processing
	 * @param doCopy_ Whether to copy the network for thread safety
	 * @param list_ List of samples specific to this job
	 * @param set_ Complete set of all available samples
	 * @param setLock_ Thread synchronization lock for the sample set
	 * @param jid_ Unique job identifier (-1 for non-tracked jobs)
	 * @param jpe_ Number of jobs per epoch for coordination
	 */
	JobData(final CellNet net_, final ArrayBlockingQueue<JobResults> jrq_, 
			final int epoch_, final int maxSamples_, final double alpha_,
			boolean backprop_, final float weightMult_, final boolean sort_, final boolean doCopy_,
			ArrayList<Sample> list_, Sample[] set_, ReentrantReadWriteLock setLock_, final int jid_, final int jpe_){

//		if(jid_>-1) {
//			System.err.println("Called  new JobData("+(net_==null ? "null" : ""+net_.dims)+", "+
//					jrq_+", "+epoch_+", "+maxSamples_+", "+alpha_+", "+backprop_+", "+
//					sort_+", "+(list_==null ? "null" : ""+list_.size())+", "+set_+", "+jid_+")");
//
//			System.err.println("JD: maxSamples_="+maxSamples_+", list_="+
//					(list_==null ? "null" : ""+list_.size()));
//		}
		immutableNet=net_;
		jobResultsQueue=jrq_;
		
		epoch=epoch_;
		maxSamples=maxSamples_;
		alpha=alpha_;
		backprop=backprop_;
		weightMult=weightMult_;
		sort=sort_;
		doCopy=doCopy_;
		list=list_;
		set=set_;
		setLock=setLock_;
		jid=jid_;
		jobsPerEpoch=jpe_;
		assert(jid_==-1 || maxSamples>0) : 
			maxSamples+", "+(list_==null ? "null" : ""+list_.size())+", "+epoch+", "+jid_;
	}
	
	/** Returns a string representation of key job parameters for debugging.
	 * @return String containing epoch, max samples, list size, backprop, and sort settings */
	public String toString() {
		return "jD: epoch="+epoch+", max="+maxSamples+", len="+list.size()+", back="+backprop+", sort="+sort;
	}
	
	/** Immutable reference to the neural network used for this job */
	final CellNet immutableNet;
	/** Optional mutable network copy for thread-safe modifications */
	CellNet mutableNet;//Optional and mutable
	/** Thread-safe queue for submitting completed job results */
	final ArrayBlockingQueue<JobResults> jobResultsQueue;

	/** Current training epoch number for this job */
	final int epoch;
	/** Maximum number of samples to process in this job */
	final int maxSamples;
	/** Learning rate for gradient descent updates */
	final double alpha;
	/** Whether to perform backpropagation during sample processing */
	final boolean backprop;
	/** Weight multiplication factor applied during training for regularization */
	final float weightMult;
	/** Whether to sort samples by error or other criteria before processing */
	final boolean sort;
	/** Whether to create a mutable copy of the network for thread safety */
	final boolean doCopy;
	/** Unique job identifier, -1 for non-tracked jobs */
	final int jid;
	/** Total number of jobs per epoch for coordination and progress tracking */
	final int jobsPerEpoch;

	/** List of samples assigned specifically to this job for processing */
	final ArrayList<Sample> list;//Only samples for this job
	/** Complete array of all available samples in the dataset */
	final Sample[] set;//All samples
	/** Thread synchronization lock for safe access to the complete sample set */
	final ReentrantReadWriteLock setLock;
	
	/**
	 * Sentinel job instance used to signal thread termination in producer-consumer patterns
	 */
	static final JobData POISON=new JobData(null, null, -1, -1, -1, false, 0.0f, false, false, null, null, null, -1, -1);
}
