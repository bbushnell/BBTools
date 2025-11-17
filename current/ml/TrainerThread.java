package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import shared.Tools;
import structures.FloatList;

/**
 * Worker thread that executes neural network training epochs for machine learning tasks.
 * Manages training batches, validation, annealing, and statistical collection across epochs.
 * Handles parallel execution of jobs across multiple worker threads with result synchronization.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class TrainerThread extends Thread {
	
	/**
	 * Constructs a TrainerThread with configuration copied from parent trainer.
	 * Initializes network copies, queues, and all training parameters from parent.
	 * Creates subnet copies for parallel worker thread processing.
	 *
	 * @param parent_ The parent Trainer that manages this thread
	 * @param net0_ The base neural network to train
	 */
	public TrainerThread(Trainer parent_, CellNet net0_) {
		parent=parent_;
		net0=net0_;
		net00=net0.copy(false);
		randyAnneal=new Random(net00.seed);
		
		jobsPerEpoch=parent.jobsPerEpoch;
		
		orderedJobs=parent.orderedJobs;
		launchInThread=parent.launchInThread;
		jobResultsQueue=new ArrayBlockingQueue<JobResults>(jobsPerEpoch);
		workerQueue=parent.workerQueue;
		launchQueue=(launchInThread ? new ArrayBlockingQueue<JobData>(2) : null);
		
		training=parent.training;
		
		maxEpochs=parent.maxEpochs;
		targetError=parent.targetError;
		targetFPR=parent.targetFPR;
		targetFNR=parent.targetFNR;
		crossoverFpMult=parent.crossoverFpMult;

		sortAll=parent.sortAll;
		sort=parent.sort;
		sortInThread=parent.sortInThread;
		shuffleSubset=parent.shuffleSubset;
		allowMTSort=parent.allowMTSort;
		
		alphaZero=parent.alphaZero;
		alphaMult=parent.alphaMult;
		alphaMult2=parent.alphaMult2;
		peakAlphaEpoch=parent.peakAlphaEpoch;
		
		alphaIncrease=parent.alphaIncrease;
		alphaEpochs=parent.alphaEpochs;
		alphaDropoff=parent.alphaDropoff;
		
		annealStrength0=parent.annealStrength0;
		annealProb=parent.annealProb;
		annealMult2=parent.annealMult2;
		annealDropoff0=parent.annealDropoff0;
		
		minAnnealEpoch=parent.minAnnealEpoch;
		maxAnnealEpoch=parent.maxAnnealEpoch;

		fractionPerEpoch0=parent.fractionPerEpoch0;
		fractionPerEpoch2=parent.fractionPerEpoch2;
		fpeStart=parent.fpeStart;
		
		positiveTriage=parent.positiveTriage;
		negativeTriage=parent.negativeTriage;
		startTriage=parent.startTriage;
		
		printStatus=parent.printStatus;
		printInterval=parent.printInterval;
		
		dumpRate=parent.dumpRate;
		dumpEpoch=parent.dumpEpoch;
		minWeightEpoch=parent.minWeightEpoch;
		minWeightEpochInverse=1f/minWeightEpoch;
		
		subnets=new CellNet[jobsPerEpoch];
		for(int i=0; i<subnets.length; i++){
			subnets[i]=net0.copy(false);
		}
		
		setLock=parent.useSetLock ? new ReentrantReadWriteLock() : null;
		flist=new FloatList(Tools.max(256, parent.data.maxSubsetSize()));
	}
	
	@Override
	public void run() {
		if(setLock!=null) {
//			System.err.println("Initial Lock");
			setLock.writeLock().lock();
		}
		data=(parent.networksPerCycle==1 ? parent.data : parent.data.copy());
		validateSet=(parent.validateSet==null ? null : 
			(parent.networksPerCycle==1 ? parent.validateSet : parent.validateSet.copy()));
		
		alpha=alphaZero;
		annealStrength=annealStrength0;
		int annealEpochs=Tools.min(maxEpochs, maxAnnealEpoch)-minAnnealEpoch;
		annealDropoff=1.0/Math.pow(annealMult2, 1.0/annealEpochs);
//		fractionPerEpoch=fractionPerEpoch0;
		
		if(launchInThread) {new LaunchThread().start();}
		
		runEpochs();
		
		parent.networkQueue.add(net0);///TODO: Ensure net's stats are correct
		if(launchInThread) {launchQueue.add(JobData.POISON);}
		if(setLock!=null) {
//			System.err.println("Final Unlock");
			setLock.writeLock().unlock();
		}
	}

	/**
	 * Executes training epochs until convergence criteria are met.
	 * Each epoch includes training interval, validation interval, and statistics calculation.
	 * Handles periodic dumping of samples and status printing based on configuration.
	 * @return Number of epochs completed
	 */
	private int runEpochs() {
		while(currentEpoch<maxEpochs && (bestErrorRate>targetError)) {
			mprof.reset();
			currentEpoch++;
			
			if(currentEpoch==dumpEpoch && dumpRate>0) {
				dump(data);
			}
			
			if(training) {
				assert(jobResultsQueue.size()==0);
				assert(parent.networksPerCycle>1 || workerQueue.size()==0);
				runTrainingInterval();
				assert(jobResultsQueue.size()==0);
				assert(parent.networksPerCycle>1 || workerQueue.size()==0);
			}
			
			final boolean print=handlePrintInterval();
			validateThisEpoch=print; //TODO: can make this more frequent, esp. when not printing
			
			boolean validated=false;
			if(validateThisEpoch | print){
				assert(jobResultsQueue.size()==0);
				assert(parent.networksPerCycle>1 || workerQueue.size()==0);
				runTestingInterval(validateSet.samples);
				validated=true;
				assert(jobResultsQueue.size()==0);
				assert(parent.networksPerCycle>1 || workerQueue.size()==0);
//				assert(false) : validateSet.samples.length;
			}
			
			if(validated) {
				calcNetStats(!training || Trainer.setCutoffForEvaluation);
				parent.compareWithBest(net0.copy(false));
			}
			
			mprof.log();//11: 11078/10597
			//System.err.println("M finished epoch "+currentEpoch);
		}
		return currentEpoch;
	}
	
	/**
	 * Legacy method for sample dumping during training.
	 * Retains high-error samples and discards low-error samples to focus training.
	 * Sorts positive and negative samples separately by error rate before selective retention.
	 *
	 * @param data The sample set to filter and reduce
	 * @deprecated Use dump(SampleSet) instead
	 */
	@Deprecated
	private void dump_old(SampleSet data){
//		System.err.println("Before: samples="+data.samples.length);
//		System.err.println("Before: subsets="+parent.subsets);
//		System.err.println("Before: fpeMult="+fpeMult);
//		System.err.println("Before: subslen="+data.currentSubset(currentEpoch).samples.length);
//		System.err.println("Before: fpe=    "+calcFractionPerEpoch());
//		System.err.println("Before: active ="+(int)(data.currentSubset(currentEpoch).samples.length*calcFractionPerEpoch()));
		
		final float retainFraction=(1-dumpRate);//Fraction completely retained
		final float retainFraction2=(parent.partialDumpFraction<1 ? 1-parent.partialDumpFraction*dumpRate : retainFraction);//Total fraction retained, including the partials
		runTestingInterval(data.samples);
		ArrayList<Sample> plist=new ArrayList<Sample>();
		ArrayList<Sample> nlist=new ArrayList<Sample>();
		
//		System.err.println("len="+data.samples.length+", ret="+retainFraction+", +ret2="+retainFraction2);
		
		for(Sample s : data.samples) {
			s.setPivot();//Necessary because the assertion failed once.  Usually works, though.
			assert(s.checkPivot()) : s.pivot+", "+s.calcPivot(); //TODO: Slow
			if(s.positive) {plist.add(s);}
			else {nlist.add(s);}
		}
//		System.err.println(plist.get(0).pivot+", "+plist.get(0).id);
//		System.err.println(nlist.get(0).pivot+", "+nlist.get(0).id);
		Collections.sort(plist);
		Collections.sort(nlist);
//		System.err.println(plist.get(0).pivot+", "+plist.get(0).id);
//		System.err.println(nlist.get(0).pivot+", "+nlist.get(0).id);
		final int pcount=(int)Math.ceil(plist.size()*retainFraction);
		final int ncount=(int)Math.ceil(nlist.size()*retainFraction);
		ArrayList<Sample> list=new ArrayList<Sample>(parent.partialDumpFraction<1 ? 
				data.samples.length : pcount+ncount);
		for(int i=0; i<pcount; i++){
			list.add(plist.get(i));
			assert(i==0 || plist.get(i).pivot<=plist.get(i-1).pivot);
		}
//		System.err.println("len="+list.size());
		
		//TODO:
		//It would be optimal to ensure the widest diversity of retained vectors, rather than discarding randomly.
		//It would also be prudent to retain relatively more high-error samples and fewer low-error samples. 
		if(parent.partialDumpFraction<1) {
			float x=0;
			final float y=1-parent.partialDumpFraction;
			for(int i=pcount; i<plist.size(); i+=1, x+=y){//Retain only some elements for the low-error samples
				if(x>=1){
					list.add(plist.get(i));
					x-=1;
				}
			}
		}
//		System.err.println("len="+list.size());
		for(int i=0; i<ncount; i++){
			list.add(nlist.get(i));
			assert(i==0 || nlist.get(i).pivot<=nlist.get(i-1).pivot);
		}
//		System.err.println("len="+list.size());
		if(parent.partialDumpFraction<1) {
			float x=0;
			final float y=1-parent.partialDumpFraction;
			for(int i=ncount; i<nlist.size(); i+=1, x+=y){//Retain only some element for the low-error samples
				if(x>=1){
					list.add(nlist.get(i));
					x-=1;
				}
			}
		}
//		System.err.println("len="+list.size());
		Collections.shuffle(list, new Random(SampleSet.shuffleSeed+1));
		data.samples=list.toArray(new Sample[0]);
		data.samplesSortedByResult=data.samples.clone();
		data.numPositive=pcount;
		data.numNegative=ncount;
		
		
		int subsets;
		boolean shrinkSubsets=parent.shrinkSubsets;
		//Note: shrinkSubsets is hard-coded as TRUE because it works better
		if(!shrinkSubsets){//This method of reducing subsets did not improve speed much.
			final int setsize=parent.setsize;
			subsets=(int)Math.ceil(parent.subsets*retainFraction2);
			if(setsize>0) {
				assert(setsize>=100) : "Setsize should be at least 100";
				subsets=Tools.max(1, data.samples.length/setsize);
				//			System.err.println("Data was organized into "+subsets+(subsets==1 ? " set." : " sets."));
			}
			subsets=Tools.mid(1, subsets, data.samples.length);
			fpeMult=1.0f;
		}else{//This method makes subsets smaller for less sorting, but also does not improve speed much (~10%).  Messes up convergence though.
			subsets=parent.subsets;
			fpeMult=1f/retainFraction2;
		}
		data.makeSubsets(subsets);
		
//		System.err.println("retainFraction2="+retainFraction2);
//		System.err.println("After:  samples="+data.samples.length);
//		System.err.println("After:  subsets="+subsets);
//		System.err.println("After:  fpeMult="+fpeMult);
//		System.err.println("After:  subslen="+data.currentSubset(currentEpoch).samples.length);
//		System.err.println("After:  fpe=    "+calcFractionPerEpoch());
//		System.err.println("After:  active ="+(int)(data.currentSubset(currentEpoch).samples.length*calcFractionPerEpoch()));

	}
	
	/**
	 * Filters sample set by retaining high-error samples and discarding low-error ones.
	 * Maintains minimum retention counts for both positive and negative samples.
	 * Reorganizes subset structure and adjusts training fraction multipliers accordingly.
	 * @param data The sample set to filter and reduce
	 */
	private void dump(SampleSet data){
		final float retainFraction=(1-dumpRate);//Fraction completely retained
		final float retainFraction2=(parent.partialDumpFraction<1 ? 1-parent.partialDumpFraction*dumpRate : retainFraction);//Total fraction retained, including the partials
		runTestingInterval(data.samples);
		ArrayList<Sample> plist=new ArrayList<Sample>();
		ArrayList<Sample> nlist=new ArrayList<Sample>();
		
//		System.err.println("len="+data.samples.length+", ret="+retainFraction+", +ret2="+retainFraction2);
		
		for(Sample s : data.samples) {
			s.setPivot();//Necessary because the assertion failed once.  Usually works, though.
			assert(s.checkPivot()) : s.pivot+", "+s.calcPivot(); //TODO: Slow
			if(s.positive) {plist.add(s);}
			else {nlist.add(s);}
		}
//		
		
		final int minRetainCount=(int)Math.ceil(Tools.max(plist.size(), nlist.size())*retainFraction);
		final int pcount=Tools.mid(minRetainCount, (int)Math.ceil(plist.size()*retainFraction), plist.size());
		final int ncount=Tools.mid(minRetainCount, (int)Math.ceil(nlist.size()*retainFraction), nlist.size());
		ArrayList<Sample> list=new ArrayList<Sample>(parent.partialDumpFraction<1 ? 
				data.samples.length : pcount+ncount);
//		System.err.println("samples="+data.samples.length+
//				", pcount="+pcount+"/"+plist.size()+", ncount="+ncount+"/"+nlist.size());

		dumpList(plist, list, pcount);
		dumpList(nlist, list, ncount);
		
//		System.err.println("len="+list.size());
		assert(data.samples.length>=list.size());
//		System.err.println(data.samples.length+", "+list.size()+", "+pcount+", "+ncount);
		final float sampleRatio=data.samples.length/(float)Tools.max(1, list.size());
		Collections.shuffle(list, new Random(SampleSet.shuffleSeed+1));
		data.samples=list.toArray(new Sample[0]);
//		System.err.println("samples="+data.samples.length);
		data.samplesSortedByResult=data.samples.clone();
		data.numPositive=pcount;
		data.numNegative=ncount;
		
		
		int subsets;
		boolean shrinkSubsets=parent.shrinkSubsets;
		//Note: shrinkSubsets is hard-coded as TRUE because it works better
		if(!shrinkSubsets){//This method of reducing subsets did not improve speed much.
			final int setsize=parent.setsize;
			subsets=(int)Math.ceil(parent.subsets*retainFraction2);
			if(setsize>0) {
				assert(setsize>=100) : "Setsize should be at least 100";
				subsets=Tools.max(1, data.samples.length/setsize);
				//			System.err.println("Data was organized into "+subsets+(subsets==1 ? " set." : " sets."));
			}
			subsets=Tools.mid(1, subsets, data.samples.length);
			fpeMult=1.0f;
		}else{//This method makes subsets smaller for less sorting, but also does not improve speed much (~10%).  Messes up convergence though.
			subsets=parent.subsets;
			fpeMult=sampleRatio;//1f/retainFraction2;
			assert(fpeMult>=1) : fpeMult;
//			System.err.println(fpeMult);
//			System.err.println(1f/retainFraction2);
		}
		data.makeSubsets(subsets);
		
//		System.err.println("retainFraction2="+retainFraction2);
//		System.err.println("After:  samples="+data.samples.length);
//		System.err.println("After:  subsets="+subsets);
//		System.err.println("After:  fpeMult="+fpeMult);
//		System.err.println("After:  subslen="+data.currentSubset(currentEpoch).samples.length);
//		System.err.println("After:  fpe=    "+calcFractionPerEpoch());
//		System.err.println("After:  active ="+(int)(data.currentSubset(currentEpoch).samples.length*calcFractionPerEpoch()));

	}
	
	/**
	 * Transfers samples from input list to output list with selective retention.
	 * Sorts input list by error rate and retains specified count plus partial fraction.
	 * Uses fractional retention for low-error samples beyond the retain count.
	 *
	 * @param inList Source list of samples sorted by error rate
	 * @param outList Destination list to receive retained samples
	 * @param retainCount Number of highest-error samples to retain completely
	 */
	private void dumpList(ArrayList<Sample> inList, ArrayList<Sample> outList, int retainCount) {
		Collections.sort(inList);
		for(int i=0; i<retainCount; i++){
			outList.add(inList.get(i));
			assert(i==0 || inList.get(i).pivot<=inList.get(i-1).pivot);
		}
		
		if(parent.partialDumpFraction<1) {
			float x=0;
			final float y=1-parent.partialDumpFraction;
			for(int i=retainCount; i<inList.size(); i+=1, x+=y){//Retain only some elements for the low-error samples
				if(x>=1){
					outList.add(inList.get(i));
					x-=1;
				}
			}
		}
	}
	
	/**
	 * Executes one training interval including job distribution and result collection.
	 * Clears network state, selects training subset, launches parallel jobs,
	 * gathers results, applies network changes, and performs annealing and alpha adjustment.
	 */
	private void runTrainingInterval() {
//		synchronized(LOCK) {
		
		assert(training);
		clearStats();
		net0.clear();
		mprof.log();//0: 614 / 564

		selectTrainingSubset();
		mprof.log();//1: 197101 / 9032

		assert(samplesThisEpoch>0) : samplesThisEpoch+", "+currentSamples.length;

		assert(jobResultsQueue.size()==0);
		assert(parent.networksPerCycle>1 || workerQueue.size()==0);
		final float weightMult=weightMult();
		
		int jobs=launchJobs(net0, currentSamples, samplesThisEpoch, training, weightMult, sort); 
		//Takes longer with sortInThread (or higher fpe) because more samples are sent
		mprof.log();//2: 90239 / 140357
//		}
		
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
			anneal();
			mprof.log();//6: 2635/2623
			adjustAlpha();
			triage();
			mprof.log();//7: 4998/5798
		}
	}
	
	/**
	 * Calculates weight multiplier for early training epochs.
	 * Returns 1.0 after minimum weight epoch, otherwise returns square root ramp-up factor.
	 * @return Weight multiplier for current epoch
	 */
	private final float weightMult() {
		if(currentEpoch>=minWeightEpoch){return 1.0f;}
		return (float)Math.sqrt((currentEpoch+1)*minWeightEpochInverse);
	}
	
	/**
	 * Executes validation/testing interval on provided sample set.
	 * Launches jobs without backpropagation, gathers results, and merges statistics.
	 * Does not modify network weights, only evaluates current performance.
	 * @param set Sample array to evaluate (typically validation set)
	 */
	private void runTestingInterval(Sample[] set) {
		final int vlines=Tools.min(parent.maxLinesV, set.length);
//		synchronized(LOCK) {
			clearStats();
			net0.clear();
			int jobs=launchJobs(net0, set, vlines, false, 1f, false);
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
			mergeStats(vlines);
			//		errorRate=weightedErrorRate;
			mprof.log();//10: 0
//		}
	}
	
	/** Acquires write lock by releasing read lock and acquiring write lock.
	 * Used to transition from shared access to exclusive access on sample sets. */
	void lock() {
//		System.err.println("Lock");
		if(setLock!=null) {
			setLock.readLock().unlock();
			setLock.writeLock().lock();
		}
	}
	
	/** Releases write lock and acquires read lock for shared access.
	 * Used to transition from exclusive access back to shared access on sample sets. */
	void unlock() {
//		System.err.println("Unlock");
		if(setLock!=null) {
			setLock.writeLock().unlock();
			setLock.readLock().lock();
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Applies simulated annealing to network weights based on current epoch.
	 * Randomly perturbs weights with decreasing strength over epochs.
	 * Annealing stops after maximum anneal epoch or when strength drops too low.
	 */
	private void anneal() {
		if(currentEpoch>maxAnnealEpoch) {annealStrength=0;}
		else if(currentEpoch>=minAnnealEpoch && annealStrength>0 && 
				annealProb>0 && currentEpoch+1<maxEpochs) {
			if(annealProb>randyAnneal.nextFloat()){
				net0.anneal((float)annealStrength, randyAnneal);
				//annealDropoff=annealDropoff*0.999f;
			}
			annealStrength=annealStrength*annealDropoff;

			if(annealDropoff0==annealDropoff && annealStrength*40<annealStrength0) {
				annealDropoff=(1-(1-annealDropoff)*0.25f);//Slow anneal dropoff
			}
		}
	}
	
	/** Adjusts learning rate (alpha) based on current epoch.
	 * Increases alpha until peak epoch, then applies exponential decay. */
	private void adjustAlpha() {
		if(currentEpoch<=peakAlphaEpoch){
			alpha+=alphaIncrease;
		}else {
			alpha*=alphaDropoff;
		}
	}
	
	/**
	 * Removes easy samples from training set to focus on difficult cases.
	 * Only activates after start triage epoch and processes all samples in current epoch.
	 * Applies separate triage rates for positive and negative samples.
	 */
	private void triage() {//Do this AFTER processing the epoch
		if(currentEpoch>=startTriage && samplesThisEpoch==currentSamples.length) {
			currentSubset.triage(currentEpoch, startTriage, positiveTriage, negativeTriage);
		}
	}
	
	/**
	 * Helper thread that launches training jobs asynchronously.
	 * Processes job requests from launch queue until poison pill received.
	 * Enables overlapping of job setup with job execution for better parallelism.
	 */
	private class LaunchThread extends Thread{
		
		//Called by start()
		@Override
		public void run(){
			for(JobData job=getJob(); job!=JobData.POISON; job=getJob()) {
				launchJobsInner(job.immutableNet, job.set, job.maxSamples, job.epoch, job.alpha, 
						job.backprop, job.weightMult, job.sort);
			}
		}
		
		/**
		 * Retrieves next job from launch queue, blocking until available.
		 * Handles interruption exceptions and continues waiting.
		 * @return Next job data from queue, or poison pill for termination
		 */
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
	
	/**
	 * Launches training or testing jobs, either in separate thread or directly.
	 * Delegates to launch thread if configured, otherwise calls inner launch method directly.
	 *
	 * @param net0 Base network for job processing
	 * @param set Sample array to process
	 * @param numSamples Number of samples to process from set
	 * @param backprop Whether to perform backpropagation (training) or just evaluation
	 * @param weightMult Weight multiplier for this epoch
	 * @param sort Whether to sort samples during processing
	 * @return Number of jobs launched
	 */
	private int launchJobs(CellNet net0, Sample[] set, int numSamples, 
			boolean backprop, float weightMult, boolean sort) {
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
	/**
	 * Core job launching logic that distributes samples across worker threads.
	 * Creates job data structures with sample sublists and submits to worker queue.
	 * Handles network copying strategy and sorting configuration per job.
	 *
	 * @param net0 Base network for job processing
	 * @param set Complete sample array
	 * @param numSamples_ Maximum samples to process
	 * @param epoch Current epoch number
	 * @param alpha Learning rate for this job
	 * @param backprop Whether jobs should perform backpropagation
	 * @param weightMult Weight multiplier for gradients
	 * @param sort Whether to enable per-thread sorting
	 * @return Number of jobs launched
	 */
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
		
		final int samplesToSend=(sortFlag ? set.length : numSamples);//This is the only difference
		final int listLen=(samplesToSend+jobsPerEpoch-1)/jobsPerEpoch;

		int sent=0;
		int jobs=0;
		final CellNet immutableNet=Trainer.copyNetInWorkerThread ||  Trainer.setNetInWorkerThread ? net0.copy(false) : null;
		for(int jid=0; jid<jobsPerEpoch; jid++){
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
//			final JobData job=new JobData(immutableNet, jobResultsQueue, epoch, toProcess, alpha, 
//					backprop, weightMult, sort, true, list, null, jid);
			
			final JobData job;
			if(Trainer.copyNetInWorkerThread){
				job=new JobData(immutableNet, jobResultsQueue, epoch, toProcess, alpha, 
						backprop, weightMult, sort, true, list, null, setLock, jid, jobsPerEpoch);
			}else{
				if(Trainer.setNetInWorkerThread) {
					job=new JobData(Trainer.setNetInWorkerThread ? immutableNet : null, jobResultsQueue, epoch, toProcess, alpha, 
							backprop, weightMult, sort, false, list, null, setLock, jid, jobsPerEpoch);
				}else{
					synchronized(net0) {
						synchronized(subnets[jid]) {
							subnets[jid].setFrom(net0, false);
						}
					}
					job=new JobData(null, jobResultsQueue, epoch, toProcess, alpha, 
							backprop, weightMult, sort, false, list, null, setLock, jid, jobsPerEpoch);
				}
				job.mutableNet=subnets[jid];
			}
			jobs++;
//			System.err.println(job);
			sent+=list.size();
			workerQueue.add(job);
		}
		
		
		assert((sent==numSamples && !sort) || (sent==samplesToSend && sort)) : 
			"sort="+sort+", sent="+sent+", samples="+numSamples+", samples_="+numSamples_+", "+
			"toSend="+samplesToSend+", setlen="+set.length+", jobs="+jobsPerEpoch+", listlen="+listLen;
		
//		if(sortFlag & shuffleSubset && ((epoch&7)==3)) {
//			currentSubset.shuffle();
//		}
		return jobs;
	}
	
	/**
	 * Alternative job launching method that uses set locking for thread safety.
	 * Distributes samples by job ID rather than creating explicit sublists.
	 * Unlocks before job submission to allow concurrent access to sample set.
	 *
	 * @param net0 Base network for job processing
	 * @param set Complete sample array shared across jobs
	 * @param numSamples_ Maximum samples to process
	 * @param epoch Current epoch number
	 * @param alpha Learning rate for this job
	 * @param backprop Whether jobs should perform backpropagation
	 * @param weightMult Weight multiplier for gradients
	 * @param sort Whether to enable sorting (disabled in this method)
	 * @return Number of jobs launched
	 */
	private int launchJobs_SetLock(CellNet net0, Sample[] set, int numSamples_, int epoch, double alpha, 
			boolean backprop, float weightMult, boolean sort) {
		
		final int numSamples=Tools.min(numSamples_, set.length);

//		assert(!sort);
		int sent=0;
		int jobs=0;
		
		//This does not seem to change anything...
//		final CellNet immutableNet=Trainer.copyNetInWorkerThread ||  Trainer.setNetInWorkerThread ? net0.copy(false) : null;
		final CellNet immutableNet=Trainer.copyNetInWorkerThread ||  Trainer.setNetInWorkerThread ? net00.setFrom(net0, false) : null;
		unlock();
		for(int jid=0; jid<jobsPerEpoch; jid++){
			
			final int toProcess=(numSamples-jid+jobsPerEpoch-1)/jobsPerEpoch;//I think this is right
			final JobData job;
			if(Trainer.copyNetInWorkerThread){
				job=new JobData(immutableNet, jobResultsQueue, epoch, numSamples, alpha, 
						backprop, weightMult, false, true, null, set, setLock, jid, jobsPerEpoch);
			}else{
				if(Trainer.setNetInWorkerThread) {
					job=new JobData(Trainer.setNetInWorkerThread ? immutableNet : null, jobResultsQueue, epoch, numSamples, alpha, 
							backprop, weightMult, false, false, null, set, setLock, jid, jobsPerEpoch);
				}else{
					synchronized(net0) {
						synchronized(subnets[jid]) {
							subnets[jid].setFrom(net0, false);
						}
					}
					job=new JobData(null, jobResultsQueue, epoch, numSamples, alpha, 
							backprop, weightMult, false, false, null, set, setLock, jid, jobsPerEpoch);
				}
				job.mutableNet=subnets[jid];
			}
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
			//Can be moved outside of the loop to ensure read-write exclusion on net0, if needed
			while(next<numJobs && results[next]!=null){//This loop actually seems to take very little time
				final JobResults job=results[next];
				gatherStats(job);
				if(accumulate && job.net!=null) {net0.accumulate(job.net);}
				else {assert(!accumulate || jobsPerEpoch>samplesThisEpoch);}
				next++;
			}
		}
	}
	
	/**
	 * Determines if status should be printed for current epoch.
	 * Uses exponentially increasing intervals early in training, then fixed intervals.
	 * Always prints on final epoch regardless of interval settings.
	 * @return true if status should be printed this epoch
	 */
	private boolean handlePrintInterval() {
		boolean print=!training || currentEpoch==maxEpochs;
		if(/*printStatus && */currentEpoch>=nextPrintEpoch) {
			print=true;
			if(currentEpoch<printInterval) {
				nextPrintEpoch=nextPrintEpoch*4;
				if(nextPrintEpoch>printInterval) {
					nextPrintEpoch=printInterval;
				}
			}
			if(currentEpoch>=nextPrintEpoch) {
				nextPrintEpoch+=printInterval;
			}
			nextPrintEpoch=Tools.min(nextPrintEpoch, maxEpochs);
		}
		return print;
	}
	
	/**
	 * Calculates fraction of training samples to use in current epoch.
	 * Uses linear interpolation from initial fraction to final fraction.
	 * Applies fraction multiplier to account for sample dumping effects.
	 * @return Fraction of samples to use this epoch (0.0 to 1.0+)
	 */
	float calcFractionPerEpoch() {
		if(currentEpoch<fpeStart){
			float f=fractionPerEpoch0+(1-(currentEpoch/(float)fpeStart))*(1-fractionPerEpoch0);
			assert(f<=1 && f>=fractionPerEpoch0) : f+", "+fractionPerEpoch0+", "+currentEpoch+", "+fpeStart;
			return f;
		}
		
		final int start=Tools.mid(fpeStart, 0, maxEpochs);
		final int fpeEpochs=maxEpochs-start;
		final int epochsSinceStart=currentEpoch-start;
		final float incr=(fractionPerEpoch2-fractionPerEpoch0)/fpeEpochs;
		float fractionPerEpoch=(fpeEpochs<1 ? fractionPerEpoch0 : 
			Tools.min(1, (fractionPerEpoch0+incr*(epochsSinceStart))));
		assert(Tools.mid(fractionPerEpoch2, fractionPerEpoch0, fractionPerEpoch)==fractionPerEpoch) : 
			"start="+start+", current="+currentEpoch+", fpeEpochs="+fpeEpochs+", epochsSinceStart="+
			epochsSinceStart+", incr="+incr+",\n fractionPerEpoch0="+fractionPerEpoch0+
			", fractionPerEpoch2="+fractionPerEpoch2+
			",\n fractionPerEpoch="+fractionPerEpoch;
		fractionPerEpoch*=fpeMult;
		return fractionPerEpoch;
	}
	
	//TODO: Use calcFractionPerEpoch instead of recalculating
	/**
	 * Calculates exact number of samples to process in current epoch.
	 * Applies fraction per epoch calculation to subset size.
	 * Ensures minimum of 4 samples or jobs per epoch, whichever is larger.
	 *
	 * @param currentSubset The subset being processed this epoch
	 * @return Number of samples to process from the subset
	 */
	int calcSamplesThisEpoch(Subset currentSubset) {
		final int len=currentSubset.samples.length;
		if(currentEpoch>=currentSubset.nextFullPassEpoch) {//This should really go outside the function.
			currentSubset.nextFullPassEpoch=currentEpoch+SampleSet.subsetInterval;
			return len;
		}
		if(currentEpoch<fpeStart){
			float f=fractionPerEpoch0+(1-(currentEpoch/(float)fpeStart))*(1-fractionPerEpoch0);
			assert(f<=1 && f>=fractionPerEpoch0) : f+", "+fractionPerEpoch0+", "+currentEpoch+", "+fpeStart;
			return (int)Tools.min(len, Tools.max(4, jobsPerEpoch, len*f));
		}
		
		final int start=Tools.mid(fpeStart, 0, maxEpochs);
		final int fpeEpochs=maxEpochs-start;
		final int epochsSinceStart=currentEpoch-start;
		final float incr=(fractionPerEpoch2-fractionPerEpoch0)/fpeEpochs;
		float fractionPerEpoch=(fpeEpochs<1 ? fractionPerEpoch0 : 
			Tools.min(1, (fractionPerEpoch0+incr*(epochsSinceStart))));
		assert(Tools.mid(fractionPerEpoch2, fractionPerEpoch0, fractionPerEpoch)==fractionPerEpoch) : 
			"start="+start+", current="+currentEpoch+", fpeEpochs="+fpeEpochs+", epochsSinceStart="+
			epochsSinceStart+", incr="+incr+",\n fractionPerEpoch0="+fractionPerEpoch0+
			", fractionPerEpoch2="+fractionPerEpoch2+
			",\n fractionPerEpoch="+fractionPerEpoch;
		fractionPerEpoch*=fpeMult;
		final int ste=(int)Tools.min(currentSamples.length, Tools.max(4, jobsPerEpoch, len*fractionPerEpoch));
		return ste;
	}
	
	/**
	 * Selects and prepares training subset for current epoch.
	 * Applies shuffling, sorting, or partial sorting based on epoch and configuration.
	 * Uses different sorting strategies on different epoch intervals for variety.
	 */
	private void selectTrainingSubset() {
		currentSubset=data.currentSubset(currentEpoch);
		currentSamples=currentSubset.samples;
		samplesThisEpoch=calcSamplesThisEpoch(currentSubset);
		assert(setLock==null || setLock.writeLock().isHeldByCurrentThread());
		final int mod8=currentEpoch&7, mod64=currentEpoch&7;
		if(shuffleSubset && mod8==Trainer.SHUFFLEMOD) {
			currentSubset.shuffle();
			return;
		}else if(!sort || sortInThread) {
			return;
		}else if(mod64==5) {
			currentSubset.sortSamples(1f, allowMTSort);
		}else if(sortAll || mod8==5) {
//			currentSubset.sortSamples(1f, allowMTSort);
			currentSubset.sortSamples2(1f, samplesThisEpoch, allowMTSort, flist);
		}else if(mod8==1) {
			currentSubset.sortSamples2(fractionPerEpoch0*6, samplesThisEpoch, allowMTSort, flist);
		}else if(mod8==3 || mod8==7){
//			currentSubset.sortSamples(fractionPerEpoch0*3, allowMTSort);
			currentSubset.sortSamples2(fractionPerEpoch0*3, samplesThisEpoch, allowMTSort, flist);
		}
	}
	
	/** Resets all statistical accumulators to zero for new epoch.
	 * Clears error sums and confusion matrix counters. */
	private void clearStats() {
		rawErrorSum=0;
		weightedErrorSum=0;
		tpSum=tnSum=fpSum=fnSum=0;
	}
	
	/**
	 * Accumulates statistics from completed job into epoch totals.
	 * Adds error sums and confusion matrix values from job results.
	 * @param job Completed job containing statistics to accumulate
	 */
	private void gatherStats(JobResults job) {
		assert(job.errorSum>=0) : job.errorSum;
		assert(job.weightedErrorSum>=0) : job.weightedErrorSum;
		rawErrorSum+=job.errorSum;
		weightedErrorSum+=job.weightedErrorSum;
		tpSum+=job.tpSum;
		tnSum+=job.tnSum;
		fpSum+=job.fpSum;
		fnSum+=job.fnSum;
	}
	
	/**
	 * Converts accumulated statistics into rates and percentages.
	 * Calculates error rates, false positive rates, false negative rates, etc.
	 * Updates network statistics with computed values.
	 * @param samples Total number of samples processed to normalize rates
	 */
	private void mergeStats(int samples) {
		final int outputs=(data!=null ? data.numOutputs() : validateSet.numOutputs());
		final double invSamples=1.0/samples;
		final double invOutputs=1.0/outputs;
//		final double e1=net0.errorSum*invSamples;
//		final double we1=net0.weightedErrorSum*invSamples;
		
//		assert(false) : "TP="+tpSum+", TN="+tnSum+", FP="+fpSum+", FN="+fnSum+"; sum="+(tpSum+tnSum+fpSum+fnSum);
		fpRate=fpSum*invSamples*invOutputs;
		fnRate=fnSum*invSamples*invOutputs;
		tpRate=tpSum*invSamples*invOutputs;
		tnRate=tnSum*invSamples*invOutputs;
		final double e3=rawErrorSum*invSamples;
		final double we3=weightedErrorSum*invSamples;
		assert(e3>=0) : e3+", "+we3+", "+samples;
		assert(we3>=0) : e3+", "+we3+", "+samples;
		
//		assert!Double.isNaN(e3) : invSamples;

		rawErrorRate=e3;//Tools.max(e1, e3);//, e2);
		weightedErrorRate=we3;//Tools.max(we1, we3);//, e2);
		setNetStats(net0);
	}
	
	/**
	 * Calculates network performance statistics using validation set.
	 * Determines optimal cutoff threshold based on target FPR, FNR, or crossover point.
	 * Sorts validation set and computes rates at calculated threshold.
	 * @param retainOldCutoff Whether to keep existing cutoff or recalculate
	 */
	void calcNetStats(boolean retainOldCutoff) {
		
		SampleSet set=validateSet;
		if(validateThisEpoch){
			if(retainOldCutoff) {
				set.sortByValue();
				lastCutoff=net0.cutoff;
				fpRate=set.calcFPRFromCutoff(lastCutoff);
				fnRate=set.calcFNRFromCutoff(lastCutoff);
			}else if(crossoverFpMult>0) {
				set.sortByValue();
				lastCutoff=set.calcCutoffFromCrossover(crossoverFpMult);
				fpRate=set.calcFPRFromCutoff(lastCutoff);
				fnRate=set.calcFNRFromCutoff(lastCutoff);
			}else if(targetFPR>=0) {
				set.sortByValue();
				fpRate=targetFPR;
				fnRate=set.calcFNRFromFPR(targetFPR);
				lastCutoff=set.calcCutoffFromFPR(fpRate);
				fpRate=set.calcFPRFromCutoff(lastCutoff);
			}else if(targetFNR>=0) {
				set.sortByValue();
				fnRate=targetFNR;
				fpRate=set.calcFPRFromFNR(targetFNR);
				lastCutoff=set.calcCutoffFromFNR(fnRate);//TODO: Test this function
				fnRate=set.calcFNRFromCutoff(lastCutoff);
			}else{
				lastCutoff=Trainer.cutoffForEvaluation;
//				fpRate=validateSet.calcFPRFromCutoff(lastCutoff);
//				fnRate=validateSet.calcFNRFromCutoff(lastCutoff);
			}
			tpRate=(set.numPositive/(double)set.samples.length)-fnRate;
			tnRate=(set.numNegative/(double)set.samples.length)-fpRate;
		}else {
			//Not sure what to do here if retainOldCutoff=true
			lastCutoff=Trainer.cutoffForEvaluation;
		}
		
//		assert(false) : crossoverFpMult+", "+lastCutoff+", "+fpRate+", "+fnRate+", "+targetFPR+", "+set.numPositive+", "+set.numNegative;
//		assert(fnRate<1) : fnRate+", "+targetFNR+", "+targetFPR;//This was added because one time I forgot to include positive samples
		setNetStats(net0);
	}
	
	/**
	 * Transfers calculated statistics into network object.
	 * Updates network's performance metrics and training parameters.
	 * @param net Network to update with current statistics
	 */
	private void setNetStats(CellNet net) {
		net.fpRate=(float) fpRate;
		net.fnRate=(float) fnRate;
		net.tpRate=(float) tpRate;
		net.tnRate=(float) tnRate;
		net.errorRate=(float) rawErrorRate;
		net.weightedErrorRate=(float) weightedErrorRate;
		net.alpha=(float) alpha;
		net.annealStrength=(float) annealStrength;
		net.epoch=currentEpoch;
		if(lastCutoff!=999) {net.setCutoff((float)lastCutoff);}
	}
	
	/*--------------------------------------------------------------*/
	
	/** Indicates whether training completed successfully without errors.
	 * @return true if no error state was encountered */
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parent trainer that created and manages this thread */
	private final Trainer parent;

	/** Primary network being trained by this thread */
	private final CellNet net0;//Basis network
	/** Copy of primary network used for worker thread initialization */
	private final CellNet net00;//A copy
	/** Array of network copies for parallel worker threads */
	private final CellNet[] subnets; //Copies for worker threads (if they don't make copies themselves)

	/*--------------------------------------------------------------*/
	
	/** Training sample set used by this thread */
	private SampleSet data;
	/** Validation sample set for performance evaluation */
	private SampleSet validateSet;
	
	/** Reusable float list for sorting operations */
	private final FloatList flist;
	/** Current training subset being processed */
	private Subset currentSubset;
	/** Array of samples from current subset */
	private Sample[] currentSamples;
	
	/** Read-write lock for coordinating access to shared sample sets */
	private final ReentrantReadWriteLock setLock;

	/*--------------------------------------------------------------*/
	
//	private final long annealSeed;
	/** Number of parallel jobs to launch per training epoch */
	private final int jobsPerEpoch;
//	private final int jobsPerBatch; //TODO: Change threads to this.
	
	/** Whether to process job results in submission order for determinism */
	private final boolean orderedJobs; //Without ordered, very very slight nondeterminism.
	/** Queue for collecting results from completed worker jobs */
	private final ArrayBlockingQueue<JobResults> jobResultsQueue;
	/** Queue for distributing jobs to worker threads */
	private final ArrayBlockingQueue<JobData> workerQueue;
	/** Queue for asynchronous job launching via LaunchThread */
	private final ArrayBlockingQueue<JobData> launchQueue;
	/** Profiler for timing different phases of training execution */
	final Profiler mprof=new Profiler("M", 13);
	
	/** Whether this thread performs training (true) or evaluation only (false) */
	private final boolean training;
	
	/*--------------------------------------------------------------*/
	
	/** Maximum number of training epochs to execute */
	final int maxEpochs;
	/** Target error rate for early stopping */
	final float targetError;
	/** Target false positive rate for threshold selection */
	final float targetFPR;
	/** Target false negative rate for threshold selection */
	final float targetFNR;
	/** Multiplier for finding FPR/FNR crossover point for threshold selection */
	final float crossoverFpMult;
	
	/*--------------------------------------------------------------*/

	/** Whether to sort all samples every epoch */
	final boolean sortAll;
	/** Whether to enable sample sorting */
	final boolean sort;
	/** Whether to perform sorting within worker threads */
	final boolean sortInThread;
	/** Whether to shuffle training subsets periodically */
	final boolean shuffleSubset; //Only if sortInThread is true
	/** Whether to launch jobs asynchronously via LaunchThread */
	final boolean launchInThread;
	/** Whether to allow multi-threaded sorting of samples */
	final boolean allowMTSort;
	
	/*--------------------------------------------------------------*/
	
	/** Initial learning rate at start of training */
	final double alphaZero;
	/** Learning rate multiplier (unused in current implementation) */
	final double alphaMult;
	/** Secondary learning rate multiplier (unused in current implementation) */
	final double alphaMult2;
	/** Epoch at which learning rate stops increasing and begins decaying */
	final int peakAlphaEpoch;
	
	/** Amount to increase learning rate each epoch until peak */
	final double alphaIncrease;
	/** Number of epochs over which to adjust learning rate (unused) */
	final int alphaEpochs;
	/** Multiplicative factor for learning rate decay after peak epoch */
	final double alphaDropoff;
	
	/** Initial strength of simulated annealing weight perturbations */
	final float annealStrength0;
	/** Probability of applying annealing in each eligible epoch */
	final float annealProb;
	/** Factor for calculating annealing strength decay over epochs */
	final float annealMult2;
	/** Initial dropoff rate for annealing strength decay */
	final double annealDropoff0;
	
	/** Earliest epoch at which annealing can begin */
	final int minAnnealEpoch;
	/** Latest epoch at which annealing can occur */
	final int maxAnnealEpoch;

	/** Initial fraction of training samples to use per epoch */
	private final float fractionPerEpoch0;
	/** Final fraction of training samples to use per epoch */
	private final float fractionPerEpoch2;
//	private float fractionPerEpoch;
	/** Multiplier applied to fraction per epoch to account for sample dumping */
	private float fpeMult=1.0f;
	/**
	 * Epoch at which to begin transitioning from initial to final fraction per epoch
	 */
	private final int fpeStart;
	
	/** Fraction of easy positive samples to remove during triage */
	private final float positiveTriage;
	/** Fraction of easy negative samples to remove during triage */
	private final float negativeTriage;
	/** Epoch at which to begin removing easy samples via triage */
	private final int startTriage;
	
	/** Whether to print training status messages */
	private final boolean printStatus;
	/** Number of epochs between status prints */
	private final int printInterval;
	
	/** Fraction of samples to discard during periodic dumping */
	private final float dumpRate;
	/** Epoch at which to perform sample dumping */
	private final int dumpEpoch;

	/** Epoch after which weight multiplier becomes 1.0 */
	private final int minWeightEpoch;
	/** Reciprocal of minimum weight epoch for efficient calculation */
	private final float minWeightEpochInverse;
	
	/*--------------------------------------------------------------*/
	
	/** Lowest error rate achieved during training */
	float bestErrorRate=999;
	/** Lowest false negative rate achieved during training */
	float bestFNR=999;

	/** Sum of raw errors accumulated across current epoch */
	double rawErrorSum=0;
	/** Sum of weighted errors accumulated across current epoch */
	double weightedErrorSum=0;
	long tpSum=0, tnSum=0, fpSum=0, fnSum=0;

	/** Raw error rate calculated from accumulated errors */
	double rawErrorRate=999f;
	/** Weighted error rate calculated from accumulated weighted errors */
	double weightedErrorRate=999f;
	
	double fpRate=0, fnRate=0, tpRate, tnRate;
	/** Decision threshold used for binary classification in last evaluation */
	double lastCutoff=999f;
	
	/** Current strength of annealing weight perturbations */
	double annealStrength;
	/** Current rate of annealing strength decay per epoch */
	double annealDropoff;
	/** Current learning rate for gradient descent updates */
	double alpha;
	
	/*--------------------------------------------------------------*/
	
	/** Next epoch number at which to print status information */
	private int nextPrintEpoch=1;
	
	/** Number of samples selected for processing in current epoch */
	private int samplesThisEpoch=-1;
	/** Whether validation should be performed in current epoch */
	private boolean validateThisEpoch=false;
	/** Current epoch number (0-based) */
	private int currentEpoch=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Random number generator for annealing weight perturbations */
	final Random randyAnneal;
	
//	private static final Sample[] poisonSamples=new Sample[0];
	
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
