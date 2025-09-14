package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import shared.Shared;
import shared.Tools;

/** 
 * The goal of this thread is to generate lots of networks,
 * evaluate them, and return the best candidates.
 * @author BBushnell
 *
 */
public class ScannerThread extends Thread {
	
	public ScannerThread(Trainer parent_, final int[] dims0_, final int[] minDims_, final int[] maxDims_,
			final int seedsToEvaluate_, final int seedsToReturn_, final int epochs_, 
			final int maxSamples_, final long netSeed0_, final ArrayBlockingQueue<ArrayList<Seed>> returnQueue_) {
//		System.err.println("Made scanner ste="+seedsToEvaluate_+
//				", str="+seedsToReturn_+", e="+epochs_+", s="+netSeed0_);
		parent=parent_;
		dims0=dims0_;
		minDims=minDims_;
		maxDims=maxDims_;
		seedsToEvaluate=seedsToEvaluate_;
		seedsToReturn=seedsToReturn_;
		returnQueue=returnQueue_;
		maxEpochs=epochs_;
		maxSamples=maxSamples_;
		netSeed0=netSeed0_;
		randy=Shared.threadLocalRandom(netSeed0);
		heap=new PriorityQueue<CellNet>(seedsToReturn+1);
		
		jobsPerEpoch=parent.jobsPerEpoch;
		
		orderedJobs=parent.orderedJobs;
		launchInThread=parent.launchInThread;
		jobResultsQueue=new ArrayBlockingQueue<JobResults>(jobsPerEpoch);
		workerQueue=parent.workerQueue;
		launchQueue=(launchInThread ? new ArrayBlockingQueue<JobData>(2) : null);
		
		training=parent.training;
		
		targetError=parent.targetError;
		targetFPR=parent.targetFPR;
		targetFNR=parent.targetFNR;
		crossoverFpMult=parent.crossoverFpMult;

		sortAll=parent.sortAll;
		sort=parent.sort;
		sortInThread=parent.sortInThread;
		shuffleSubset=parent.shuffleSubset;
		
		fractionPerEpoch=parent.fractionPerEpoch0;
		
		alpha=parent.alphaZero;
		setLock=parent.useSetLock ? new ReentrantReadWriteLock() : null;
	}
	
	@Override
	public void run() {
		if(setLock!=null) {
//			System.err.println("Initial Lock");
			setLock.writeLock().lock();
		}
		synchronized(this) {
			//I really only need part of the data for this.
			data=(parent.data.copy(maxSamples, 0.25f));
			validateSet=data;

			if(launchInThread) {new LaunchThread().start();}

			for(int i=0; i<seedsToEvaluate; i++) {
				final long seed=(i==0 ? netSeed0 : randy.nextLong()&Long.MAX_VALUE);
				net0=parent.randomNetwork(seed);
				data.reset();
				runEpochs();
				//			if(heap.size()<seedsToReturn){
				//				heap.offer(net0);
				//			}else{
				//				//Only offer if better than worst
				//			}
				assert(validateThisEpoch) : currentEpoch;
				heap.offer(net0);
				if(heap.size()>seedsToReturn) {
					CellNet evict=heap.poll();
//					System.err.println(evict.header());
					assert(evict.compareTo(net0)<=0) : evict+", "+net0;
				}
				net0=null;
			}

			ArrayList<Seed> list=new ArrayList<Seed>(seedsToReturn);
			while(!heap.isEmpty()){
				CellNet net=heap.poll();
//				System.err.println(net.header());
				
				list.add(new Seed(net.seed, /*net.annealSeed,*/ net.pivot()));
			}

			try {
				returnQueue.put(list);
			} catch (InterruptedException e) {
				//In this case it looks like it fails.  But there's no return so who knows.
				e.printStackTrace();
			}

			if(launchInThread) {launchQueue.add(JobData.POISON);}
		}
		if(setLock!=null) {
//			System.err.println("Final Unlock");
			setLock.writeLock().unlock();
		}
	}

	private int runEpochs() {
		mprof.reset();
		currentEpoch=0;
		while(currentEpoch<maxEpochs) {
			mprof.reset();
			currentEpoch++;
			
			if(training && currentEpoch<=maxEpochs) {
				assert(jobResultsQueue.size()==0);
				assert(parent.networksPerCycle>1 || workerQueue.size()==0);
				runTrainingInterval();
				assert(jobResultsQueue.size()==0);
				assert(parent.networksPerCycle>1 || workerQueue.size()==0);
			}
			
			mprof.log();//?: 11078/10597
			//System.err.println("M finished epoch "+currentEpoch);
		}
		
		validateThisEpoch=true;
		if(validateThisEpoch){
			assert(jobResultsQueue.size()==0);
			assert(parent.networksPerCycle>1 || workerQueue.size()==0);
			runTestingInterval(validateSet.samples);
			assert(jobResultsQueue.size()==0);
			assert(parent.networksPerCycle>1 || workerQueue.size()==0);
//			assert(false) : validateSet.samples.length;
			reviseCutoff();
		}
		mprof.log();//12: 
		
		return currentEpoch;
	}
	
	private void runTrainingInterval() {
		float weightMult=0.5f;//TODO: Test high and low.  Or maybe 0.5f.
		assert(training);
		clearStats();
		net0.clear();
		mprof.log();//0: 614 / 564

		selectTrainingSubset();
		mprof.log();//1: 197101 / 9032

		assert(samplesThisEpoch>0) : samplesThisEpoch+", "+currentSamples.length;

		assert(jobResultsQueue.size()==0);
		assert(parent.networksPerCycle>1 || workerQueue.size()==0);

		int jobs=launchJobs(net0, currentSamples, samplesThisEpoch, training, weightMult, sort); //Takes longer with sortInThread (or higher fpe) because more samples are sent
		mprof.log();//2: 90239 / 140357
		
		gatherResults(net0, jobResultsQueue, training, jobs);
		lock();
		mprof.log();//3: 561312/661228
		//System.err.println("M done waiting for threads.");
		
		synchronized(net0) {
			//System.err.println("M checking epochs.");
			assert(jobResultsQueue.size()==0);
			assert(parent.networksPerCycle>1 || workerQueue.size()==0);

			//System.err.println("M gathering.");
			mergeStats(samplesThisEpoch);
			//		errorRate=weightedErrorRate;
			mprof.log();//4: 154/143

			net0.applyChanges(samplesThisEpoch, (float)alpha);
			mprof.log();//5: 2356/2134
		}
	}
	
	private void runTestingInterval(Sample[] set) {
//		synchronized(LOCK) {
			clearStats();
			net0.clear();
			int jobs=launchJobs(net0, set, set.length, false, 1.0f, false);
			mprof.log();//8: 1738/1709
//		}
		
		gatherResults(null, jobResultsQueue, false, jobs);
		lock();
		//System.err.println("M done waiting for threads.");

//		synchronized(LOCK) {
			mprof.log();//9: 23373/23901
			//System.err.println("M checking epochs.");
			assert(jobResultsQueue.size()==0);
			assert(parent.networksPerCycle>1 || workerQueue.size()==0);

			//System.err.println("M gathering.");
			mergeStats(set.length);
			//		errorRate=weightedErrorRate;
			mprof.log();//10: 0
//		}
	}
	
	void lock() {
//		System.err.println("Lock");
		if(setLock!=null) {
			setLock.readLock().unlock();
			setLock.writeLock().lock();
		}
	}
	
	void unlock() {
//		System.err.println("Unlock");
		if(setLock!=null) {
			setLock.writeLock().unlock();
			setLock.readLock().lock();
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private class LaunchThread extends Thread{
		
		//Called by start()
		@Override
		public void run(){
			for(JobData job=getJob(); job!=JobData.POISON; job=getJob()) {
				launchJobsInner(job.immutableNet, job.set, job.maxSamples, job.epoch, job.alpha,  
						job.backprop, job.weightMult, job.sort);
			}
		}
		
		JobData getJob() {
			JobData job=null;
			while(job==null) {
				try {
					job=launchQueue.take();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			return job;
		}
			
	}
	
	private int launchJobs(CellNet net0, Sample[] set, int numSamples, boolean backprop, 
			float weightMult, boolean sort) {
		if(launchInThread) {
			JobData job=new JobData(net0, jobResultsQueue, currentEpoch, numSamples, alpha, 
					backprop, weightMult, sort, true, null, set, setLock, 0, 0);
			launchQueue.add(job);
			return jobsPerEpoch;
		}else {
			return launchJobsInner(net0, set, numSamples, currentEpoch, alpha, backprop, weightMult, sort);
		}
	}
	
	//Does not seem faster...
	private int launchJobsInner(CellNet net0, Sample[] set, int numSamples_, int epoch, double alpha, 
			boolean backprop, float weightMult, boolean sort) {
		if(setLock!=null) {return launchJobs_SetLock(net0, set, numSamples_, epoch, alpha, backprop, weightMult, sort);}
		//TODO: Eliminate this method if the above works
		
		
		//Note:  This is a little confusing because you may want to send more samples (samplesToSend)
		//than you actually want processed (numSamples) if sorting is being done per-thread.
		final int numSamples=Tools.min(numSamples_, set.length);
		final boolean sortFlag=sort && sortInThread && numSamples<set.length;
//		final float samplesPerThread=numSamples/(float)threads;
//		final int minSamplesPerThread=(int)samplesPerThread, maxSamplesPerThread=(int)Math.ceil(samplesPerThread);
		
		//Do not enable, or else different threadcounts will give different results
		//This is just for testing
		//Actually, this breaks the everything whether in sortInThread mode or not 
//		if(threads*maxSamplesPerThread<=set.length) {numSamples=threads*maxSamplesPerThread;}
		
		final int samplesToSend=(sortFlag ? set.length : numSamples);//This is the only difference
		final int listLen=(samplesToSend+jobsPerEpoch-1)/jobsPerEpoch;
		
//		int[] sizes=new int[threads];
		//This sizes[] is necessary or else it goes crazy in sortinthread mode
		//for reasons I do not understand
//		System.err.println("epoch="+epoch+", sts="+samplesToSend+", training="+training);

		int sent=0;
		int jobs=0;
		final CellNet net=net0.copy(false);
		for(int jid=0; jid<jobsPerEpoch && jid<numSamples; jid++){
			ArrayList<Sample> list=new ArrayList<Sample>(listLen);
			int idx=jid;
			while(idx<numSamples) {
				list.add(set[idx]);
				idx+=jobsPerEpoch;
			}
			final int toProcess=list.size();
			if(sortFlag) {
				while(idx<samplesToSend) {
					list.add(set[idx]);
					idx+=jobsPerEpoch;
				}
			}
			assert(toProcess>0 && toProcess<=list.size()) : toProcess+", "+list.size()+
				", e="+epoch+", j="+jid+", ns="+numSamples;
//			System.err.println("st1: toProcess="+toProcess+", size="+list.size());
//			
//			System.err.println("Calling new JobData("+net.dims+", "+jobResultsQueue.size()+", "+epoch+", "+
//					toProcess+", "+alpha+", "+backprop+", "+sort+", "+list.size()+", "+null+", "+jid+")");
			
			final JobData job=new JobData(net, jobResultsQueue, epoch, toProcess, 
					alpha, backprop, weightMult, sort, true, list, null, null, jid, 0);
//			System.err.println("st2: toProcess="+toProcess+", size="+list.size());
			jobs++;
//			System.err.println(job);
			sent+=list.size();
			workerQueue.add(job);
		}
		
		
		assert((sent==numSamples && !sort) || (sent==samplesToSend && sort)) : 
			"sort="+sort+", sent="+sent+", samples="+numSamples+", samples_="+numSamples_+", "+
			"toSend="+samplesToSend+", setlen="+set.length+", jobs="+jobsPerEpoch+", listlen="+listLen;
		
		if(sortFlag & shuffleSubset && ((epoch&7)==3)) {
			currentSubset.shuffle();
		}
		return jobs;
	}
	
	private int launchJobs_SetLock(CellNet net0, Sample[] set, int numSamples_, int epoch, double alpha, 
			boolean backprop, float weightMult, boolean sort) {
		
		final int numSamples=Tools.min(numSamples_, set.length);

//		assert(!sort);
		int sent=0;
		int jobs=0;
		
		final CellNet immutableNet=Trainer.copyNetInWorkerThread ||  Trainer.setNetInWorkerThread ? net0.copy(false) : null;
//		final CellNet immutableNet=Trainer.copyNetInWorkerThread ||  Trainer.setNetInWorkerThread ? net00.setFrom(net0, false) : null;
		unlock();
		for(int jid=0; jid<jobsPerEpoch; jid++){
			
			final int toProcess=(numSamples-jid+jobsPerEpoch-1)/jobsPerEpoch;//I think this is right
			final JobData job;
//			if(Trainer.copyNetInWorkerThread){
				job=new JobData(immutableNet, jobResultsQueue, epoch, numSamples, alpha, 
						backprop, weightMult, false, true, null, set, setLock, jid, jobsPerEpoch);
//			}else{
//				if(Trainer.setNetInWorkerThread) {
//					job=new JobData(Trainer.setNetInWorkerThread ? immutableNet : null, jobResultsQueue, epoch, numSamples, alpha, 
//							backprop, weightMult, false, false, null, set, setLock, jid, jobsPerEpoch);
//				}else{
//					synchronized(net0) {
//						synchronized(subnets[jid]) {
//							subnets[jid].setFrom(net0, false);
//						}
//					}
//					job=new JobData(null, jobResultsQueue, epoch, numSamples, alpha, 
//							backprop, weightMult, false, false, null, set, setLock, jid, jobsPerEpoch);
//				}
//				job.mutableNet=subnets[jid];
//			}
			jobs++;
			sent+=toProcess;
			workerQueue.add(job);
		}
		
		assert(sent==numSamples && jobs==jobsPerEpoch) : 
			"sort="+sort+", sent="+sent+", samples="+numSamples+", samples_="+numSamples_+", "+
			", setlen="+set.length+", jobs="+jobsPerEpoch;
		return jobs;
	}
	
	private void gatherResults(final CellNet net0, final ArrayBlockingQueue<JobResults> mq, 
			final boolean accumulate, final int numJobs) {
		if(orderedJobs) {
			gatherResultsOrdered(net0, mq, accumulate, numJobs);
		}else {
			gatherResultsDisordered(net0, mq, accumulate, numJobs);
		}
	}
	
	private void gatherResultsDisordered(final CellNet net0, final ArrayBlockingQueue<JobResults> mq, 
			final boolean accumulate, final int numJobs) {
		//System.err.println("M waiting for threads.");
		for(int i=0; i<numJobs; i++) {
			JobResults job=null;
			while(job==null) {
				try {job=mq.take();} 
				catch (InterruptedException e){e.printStackTrace();}
			}
			assert(job.epoch==currentEpoch) : job.epoch+", "+currentEpoch+", "+job.tid;
			gatherStats(job);
			if(accumulate && job.net!=null) {net0.accumulate(job.net);}
			else {assert(!accumulate || jobsPerEpoch>samplesThisEpoch);}
		}
	}
	
	private void gatherResultsOrdered(final CellNet net0, final ArrayBlockingQueue<JobResults> mq, 
			final boolean accumulate, final int numJobs) {
		JobResults[] results=new JobResults[numJobs];
		//System.err.println("M waiting for threads.");
		
		int next=0;
		for(int i=0; i<numJobs; i++) {
			{//Get a job
				JobResults job=null;
				while(job==null) {
					try {job=mq.take();} 
					catch (InterruptedException e){e.printStackTrace();}
				}
				assert(job.epoch==currentEpoch) : job.epoch+", "+currentEpoch+", "+job.tid;
				results[job.jid]=job;
			}
			
			//Process as many consecutive jobs as are available 
			while(next<numJobs && results[next]!=null){
				final JobResults job=results[next];
				gatherStats(job);
				if(accumulate && job.net!=null) {net0.accumulate(job.net);}
				else {assert(!accumulate || jobsPerEpoch>samplesThisEpoch);}
				next++;
			}
		}
	}
	
	private void selectTrainingSubset() {
		currentSubset=data.currentSubset(currentEpoch);
		currentSamples=currentSubset.samples;
		assert(currentSamples!=null && currentSamples.length>0) : currentSamples+", "+currentEpoch+", "+data.samples.length;
		if(currentEpoch<currentSubset.nextFullPassEpoch) {
			if(currentEpoch>=2) {
				samplesThisEpoch=(int)Tools.min(currentSamples.length, Tools.max(4, jobsPerEpoch, currentSamples.length*fractionPerEpoch));
				assert(samplesThisEpoch>0) : samplesThisEpoch+", "+currentSamples.length+", "+jobsPerEpoch+", "+fractionPerEpoch;
			}else{
				samplesThisEpoch=currentSamples.length;
				assert(samplesThisEpoch>0) : samplesThisEpoch+", "+currentSamples.length+", "+jobsPerEpoch+", "+fractionPerEpoch;
			}
		}else {
			samplesThisEpoch=currentSamples.length;
			currentSubset.nextFullPassEpoch=currentEpoch+SampleSet.subsetInterval;
//			System.err.println("B:"+currentSubset.nextFullPassEpoch+","+samplesThisEpoch);
			assert(samplesThisEpoch>0) : samplesThisEpoch+", "+currentSamples.length+", "+jobsPerEpoch+", "+fractionPerEpoch;
			
		}
//		System.err.println("samplesThisEpoch="+samplesThisEpoch+", ce="+currentEpoch+", dsl="+
//			data.samples.length+", csl="+currentSamples.length+", jpe="+jobsPerEpoch+", fpe="+fractionPerEpoch);
		
		if(!sort || samplesThisEpoch>=currentSamples.length || sortInThread) {return;}
		if(sortAll || (currentEpoch&3)==1) {
			currentSubset.sortSamples(1f, false);
		}else if((currentEpoch&3)==3){
			currentSubset.sortSamples(fractionPerEpoch*3, false);
		}
	}
	
	private void clearStats() {
		rawErrorSum=0;
		weightedErrorSum=0;
		tpSum=tnSum=fpSum=fnSum=0;
	}
	
	private void gatherStats(JobResults job) {
		rawErrorSum+=job.errorSum;
		weightedErrorSum+=job.weightedErrorSum;
		tpSum+=job.tpSum;
		tnSum+=job.tnSum;
		fpSum+=job.fpSum;
		fnSum+=job.fnSum;
	}
	
	private void mergeStats(int samples) {
		final double invSamples=1.0/samples;
		final double invOutputs=1.0/data.matrix.numOutputs;
//		final double e1=net0.errorSum*invSamples;
//		final double we1=net0.weightedErrorSum*invSamples;
		
//		assert(false) : "TP="+tpSum+", TN="+tnSum+", FP="+fpSum+", FN="+fnSum+"; sum="+(tpSum+tnSum+fpSum+fnSum);
		fpRate=fpSum*invSamples*invOutputs;
		fnRate=fnSum*invSamples*invOutputs;
		tpRate=tpSum*invSamples*invOutputs;
		tnRate=tnSum*invSamples*invOutputs;
		final double e3=rawErrorSum*invSamples;
		final double we3=weightedErrorSum*invSamples;
		
//		assert!Double.isNaN(e3) : invSamples;

		rawErrorRate=e3;//Tools.max(e1, e3);//, e2);
		weightedErrorRate=we3;//Tools.max(we1, we3);//, e2);
		setNetStats(net0);
	}
	
	void reviseCutoff() {
		
		SampleSet set=validateSet;
		if(validateThisEpoch){
			if(crossoverFpMult>0) {
				set.sortByValue();
				lastCutoffForTarget=set.calcCutoffFromCrossover(crossoverFpMult);
				fpRate=set.calcFPRFromCutoff(lastCutoffForTarget);
				fnRate=set.calcFNRFromCutoff(lastCutoffForTarget);
			}else if(targetFPR>=0) {
				set.sortByValue();
				fpRate=targetFPR;
				fnRate=set.calcFNRFromFPR(targetFPR);
				lastCutoffForTarget=set.calcCutoffFromFPR(fpRate);
			}else if(targetFNR>=0) {
				set.sortByValue();
				fnRate=targetFNR;
				fpRate=set.calcFPRFromFNR(targetFNR);
				lastCutoffForTarget=set.calcCutoffFromFNR(fnRate);//TODO: Test this function
			}else{
				lastCutoffForTarget=Trainer.cutoffForEvaluation;
//				fpRate=validateSet.calcFPRFromCutoff(lastCutoff);
//				fnRate=validateSet.calcFNRFromCutoff(lastCutoff);
			}
			tpRate=(set.numPositive/(double)set.samples.length)-fnRate;
			tnRate=(set.numNegative/(double)set.samples.length)-fpRate;
		}else {
			lastCutoffForTarget=Trainer.cutoffForEvaluation;
		}
		setNetStats(net0);
	}
	
	private void setNetStats(CellNet net) {
		net.fpRate=(float) fpRate;
		net.fnRate=(float) fnRate;
		net.tpRate=(float) tpRate;
		net.tnRate=(float) tnRate;
		net.errorRate=(float) rawErrorRate;
		net.weightedErrorRate=(float) weightedErrorRate;
		net.alpha=(float) alpha;
		net.annealStrength=0;
		net.epoch=currentEpoch;
		net.setCutoff((float) lastCutoffForTarget);
	}
	
	/*--------------------------------------------------------------*/
	
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final Trainer parent;
	private final int[] dims0;
	private final int[] minDims;
	private final int[] maxDims;
	private final int seedsToEvaluate;
	private final int seedsToReturn;
	private final long netSeed0;
	private final PriorityQueue<CellNet> heap;
	final ArrayBlockingQueue<ArrayList<Seed>> returnQueue;
	
	private CellNet net0;//Network being evaluated

	/*--------------------------------------------------------------*/

	private SampleSet data;
	private SampleSet validateSet;
	
	private Subset currentSubset;
	private Sample[] currentSamples;
	
	private final ReentrantReadWriteLock setLock;
	
//	private CellNet bestNetwork;//Best observed network

	/*--------------------------------------------------------------*/
	
	private final int jobsPerEpoch;
	
	private final boolean orderedJobs; //Without ordered, very very slight nondeterminism.
	private final ArrayBlockingQueue<JobResults> jobResultsQueue;
	private final ArrayBlockingQueue<JobData> workerQueue;
	private final ArrayBlockingQueue<JobData> launchQueue;
	final Profiler mprof=new Profiler("M", 13);
	
	private final boolean training;
	
	/*--------------------------------------------------------------*/
	
	final int maxEpochs;
	final int maxSamples;
	final float targetError;
	final float targetFPR;
	final float targetFNR;
	final float crossoverFpMult;
	
	/*--------------------------------------------------------------*/

	final boolean sortAll;
	final boolean sort;
	final boolean sortInThread;
	final boolean shuffleSubset; //Only if sortInThread is true
	final boolean launchInThread;
	
	/*--------------------------------------------------------------*/
	
	final double alpha;
	
	private final float fractionPerEpoch;
	
	/*--------------------------------------------------------------*/
	
	float bestErrorRate=999;
	float bestFNR=999;

	double rawErrorSum=0;
	double weightedErrorSum=0;
	long tpSum=0, tnSum=0, fpSum=0, fnSum=0;

	double rawErrorRate=999f;
	double weightedErrorRate=999f;
	
	double fpRate=0, fnRate=0, tpRate, tnRate, crossover;
	double lastCutoffForTarget=1.0f;
	
	/*--------------------------------------------------------------*/
	
	private int samplesThisEpoch=-1;
	private boolean validateThisEpoch=false;
	private int currentEpoch=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final Random randy;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
