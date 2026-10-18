package ml;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

/**
 * Encapsulates one parallel ML job submitted by TrainerThread or ScannerThread.
 * A normal job carries either a per-job sample list or a shared sample array plus
 * strided job coordinates; POISON leaves both sample sources null. Callers that
 * set doCopy false attach mutableNet before the job runs.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class JobData{
	
	/**
	 * Constructs a worker job with its network source, sample source, and result queue.
	 *
	 * @param net_ source network copied by the worker or loaded into mutableNet;
	 *              null when the caller attaches a ready mutableNet after construction
	 * @param jrq_ Queue for submitting job results after processing
	 * @param epoch_ Current training epoch number
	 * @param maxSamples_ maximum list samples to process, or exclusive set index limit
	 * @param alpha_ Learning rate for gradient updates
	 * @param backprop_ Whether to perform backpropagation during training
	 * @param weightMult_ Weighting multiplier forwarded to CellNet.processSample
	 * @param sort_ Whether to sort samples before processing
	 * @param doCopy_ Whether the worker should copy net_ instead of using mutableNet
	 * @param list_ samples specific to this job; null for shared-array jobs and POISON
	 * @param set_ complete array of all available samples; null for list jobs and POISON
	 * @param setLock_ optional lock acquired by WorkerThread before sample processing
	 * @param jid_ unique job identifier; -1 for non-tracked jobs and POISON
	 * @param jpe_ number of jobs per epoch for coordination
	 */
	public JobData(final CellNet net_, final ArrayBlockingQueue<JobResults> jrq_,
			final int epoch_, final int maxSamples_, final double alpha_,
			boolean backprop_, final float weightMult_, final boolean sort_, final boolean doCopy_,
			ArrayList<Sample> list_, Sample[] set_, ReentrantReadWriteLock setLock_, final int jid_, final int jpe_){

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
	
	/** Returns key job parameters for debugging. The len field is the list
	 * size, shared set length, or zero when no sample source is attached.
	 * @return epoch, max samples, sample-source length, backprop flag, and sort flag */
	public String toString(){
		final int len=(list!=null ? list.size() : set==null ? 0 : set.length);
		return "jD: epoch="+epoch+", max="+maxSamples+", len="+len+", back="+backprop+", sort="+sort;
	}
	
	/** Source network copied or loaded by workers; may be null with a caller-supplied mutableNet. */
	public final CellNet immutableNet;
	/** Caller-supplied mutable network used when doCopy is false. */
	CellNet mutableNet;//Optional and mutable
	/** Per-job feed-forward SIMD toggle stamped by TrainerThread after construction. Only
	 * consulted inside the Shared.SIMD &amp;&amp; Shared.SIMD_FEED_FORWARD gate, so a batch can
	 * choose scalar versus SIMD feed-forward for the SIMD_FF robustness experiment. */
	boolean simdFF=false;//Optional and mutable
	/** Queue for submitting completed job results; null only for POISON. */
	final ArrayBlockingQueue<JobResults> jobResultsQueue;

	/** Current training epoch number for this job */
	public final int epoch;
	/** Maximum list samples to process, or exclusive index limit for the shared set. */
	public final int maxSamples;
	/** Learning rate for gradient descent updates */
	public final double alpha;
	/** Whether to perform backpropagation during sample processing */
	public final boolean backprop;
	/** Weighting multiplier forwarded to CellNet.processSample */
	public final float weightMult;
	/** Whether to sort samples by error or other criteria before processing */
	public final boolean sort;
	/** Whether the worker should copy immutableNet instead of using mutableNet */
	final boolean doCopy;
	/** Unique job identifier, -1 for non-tracked jobs */
	public final int jid;
	/** Total number of jobs per epoch for coordination and progress tracking */
	public final int jobsPerEpoch;

	/** Samples assigned specifically to this job; null for shared-array jobs and POISON. */
	public final ArrayList<Sample> list;//Only samples for this job
	/** Complete array of all available samples; null for list jobs and POISON. */
	public final Sample[] set;//All samples
	/** Optional lock acquired by WorkerThread before processing. */
	public final ReentrantReadWriteLock setLock;
	
	/** Sentinel job instance used to signal thread termination in producer-consumer patterns. */
	public static final JobData POISON=new JobData(null, null, -1, -1, -1, false, 0.0f, false, false, null, null, null, -1, -1);
}
