package ml;

/** Worker result envelope containing a processed network copy plus per-job counters. */
public class JobResults implements Comparable<JobResults>{
	
	/** Constructs a result for one WorkerThread job.
	 * @param net_ processed network copy, or null when no network should be accumulated
	 * @param epoch_ epoch that produced this result
	 * @param numProcessed_ number of samples processed by the job
	 * @param tid_ worker thread ID
	 * @param jid_ job ID within the epoch
	 * @param errorSum_ sum of raw sample error magnitudes
	 * @param weightedErrorSum_ sum of weighted sample error magnitudes
	 * @param tpSum_ true-positive count
	 * @param tnSum_ true-negative count
	 * @param fpSum_ false-positive count
	 * @param fnSum_ false-negative count */
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
	
	public String toString() {
		return "jR: e="+epoch+", jid="+jid+", num="+numProcessed+", err="+errorSum+", wer="+weightedErrorSum+", fn="+fnSum+", fp="+fpSum+", tn="+tnSum+", tp="+tpSum;
	}

	/**
	 * Compares JobResults for sorting by epoch then job ID.
	 * Primary sort by epoch number, secondary sort by job ID for deterministic ordering.
	 * @param o JobResults to compare against
	 * @return Negative if this comes before o, positive if after, zero if equal
	 */
	@Override
	public int compareTo(JobResults o) {
		return epoch==o.epoch ? jid-o.jid : epoch-o.epoch;
	}
	
	final CellNet net;

	final int epoch;
	final int numProcessed;
	final int tid;
	final int jid;
	
	final double errorSum;
	final double weightedErrorSum;
	/** Confusion-matrix counts from this job, in TP/TN/FP/FN order. */
	final int tpSum, tnSum, fpSum, fnSum;
	
	/** Sentinel result for future queue consumers that need an explicit terminator. */
	static final JobResults POISON=new JobResults(null, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
}
