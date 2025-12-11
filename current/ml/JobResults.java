package ml;

public class JobResults implements Comparable<JobResults>{
	
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
	/** Total false negatives counted during this job. */
	/** Total false positives counted during this job. */
	/** Total true negatives counted during this job. */
	final int tpSum, tnSum, fpSum, fnSum;
	
	static final JobResults POISON=new JobResults(null, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
}
