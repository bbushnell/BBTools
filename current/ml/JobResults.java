package ml;

/**
 * Encapsulates computational results and performance metrics for individual
 * machine learning job executions during neural network training.
 * Stores comprehensive job execution details including error metrics,
 * classification outcomes, and network reference for distributed training scenarios.
 *
 * @author Brian Bushnell
 */
public class JobResults implements Comparable<JobResults>{
	
	/**
	 * Constructs a JobResults instance with comprehensive training metrics and validation.
	 * Validates that error sums are non-negative unless epoch is -1 (poison instance).
	 *
	 * @param net_ Neural network reference for this job
	 * @param epoch_ Training epoch number
	 * @param numProcessed_ Number of samples processed in this job
	 * @param tid_ Thread ID that executed this job
	 * @param jid_ Job ID for ordering and identification
	 * @param errorSum_ Sum of raw error values across all samples
	 * @param weightedErrorSum_ Sum of weighted error values across all samples
	 * @param tpSum_ True positive count from classification
	 * @param tnSum_ True negative count from classification
	 * @param fpSum_ False positive count from classification
	 * @param fnSum_ False negative count from classification
	 */
	JobResults(final CellNet net_, final int epoch_, final int numProcessed_, int tid_, int jid_,
			final double errorSum_, final double weightedErrorSum_,
			final int tpSum_, final int tnSum_, final int fpSum_, final int fnSum_){
		net=net_;
		
		epoch=epoch_;
		numProcessed=numProcessed_;
		tid=tid_;
		jid=jid_;
		
		errorSum=errorSum_;
		weightedErrorSum=weightedErrorSum_;
		tpSum=tpSum_;
		tnSum=tnSum_;
		fpSum=fpSum_;
		fnSum=fnSum_;
		assert(errorSum>=0 || epoch==-1) : this;
		assert(weightedErrorSum>=0 || epoch==-1) : this;
	}
	
	/**
	 * Returns a compact string representation of job results for debugging.
	 * Format: "jR: e=[epoch], jid=[jobId], num=[processed], err=[errorSum],
	 * wer=[weightedErrorSum], fn=[falseNegatives], fp=[falsePositives],
	 * tn=[trueNegatives], tp=[truePositives]"
	 *
	 * @return Compact string summary of all job metrics
	 */
	public String toString() {
		return "jR: e="+epoch+", jid="+jid+", num="+numProcessed+", err="+errorSum+", wer="+weightedErrorSum+", fn="+fnSum+", fp="+fpSum+", tn="+tnSum+", tp="+tpSum;
	}

	@Override
	public int compareTo(JobResults o) {
		return epoch==o.epoch ? jid-o.jid : epoch-o.epoch;
	}
	
	/** Neural network reference associated with this job execution */
	final CellNet net;

	/** Training epoch number when this job was executed */
	final int epoch;
	/** Number of training samples processed during this job */
	final int numProcessed;
	/** Thread ID of the worker thread that executed this job */
	final int tid;
	/** Job ID for identification and ordering within an epoch */
	final int jid;
	
	/** Sum of raw error values across all samples processed in this job */
	final double errorSum;
	/** Sum of weighted error values across all samples processed in this job */
	final double weightedErrorSum;
	/** Total false negatives counted during this job. */
	/** Total false positives counted during this job. */
	/** Total true negatives counted during this job. */
	/** Total true positives counted during this job. */
	final int tpSum, tnSum, fpSum, fnSum;
	
	/** Sentinel instance used for thread termination in producer-consumer queues */
	static final JobResults POISON=new JobResults(null, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
}
