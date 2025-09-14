package ml;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

public class JobData {
	
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
	
	public String toString() {
		return "jD: epoch="+epoch+", max="+maxSamples+", len="+list.size()+", back="+backprop+", sort="+sort;
	}
	
	final CellNet immutableNet;
	CellNet mutableNet;//Optional and mutable
	final ArrayBlockingQueue<JobResults> jobResultsQueue;

	final int epoch;
	final int maxSamples;
	final double alpha;
	final boolean backprop;
	final float weightMult;
	final boolean sort;
	final boolean doCopy;
	final int jid;
	final int jobsPerEpoch;

	final ArrayList<Sample> list;//Only samples for this job
	final Sample[] set;//All samples
	final ReentrantReadWriteLock setLock;
	
	static final JobData POISON=new JobData(null, null, -1, -1, -1, false, 0.0f, false, false, null, null, null, -1, -1);
}
