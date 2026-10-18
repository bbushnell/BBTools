package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import shared.Random;
import shared.Shared;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import shared.Tools;
import structures.FloatList;

/** Trains the supplied net0 through bounded epochs, using copies as worker
 * scratch and publishing the trained network after its epoch loop.
 * @author Brian Bushnell
 * @date 2014 */
public class TrainerThread extends Thread{

	/** Creates a trainer thread from parent-owned queues and configuration.
	 * @param parent_ parent trainer supplying data, queues, and training options
	 * @param net0_ network to train */
	public TrainerThread(Trainer parent_, CellNet net0_){
		parent=parent_;
		net0=net0_;
		net00=net0.copy(false);
		randyAnneal=shared.Shared.random(net00.seed);

		jobsPerEpoch=parent.jobsPerEpoch;

		orderedJobs=parent.orderedJobs;
		simdFFMode=parent.simdFFMode;
		randySimd=shared.Shared.random(net00.seed^0x5EED5L); //Per-net deterministic RNG for the random-SIMD_FF mode
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

	/** Runs the epoch loop, then publishes net0 before queuing launcher poison and
	 * releasing the set write lock on the normal path. */
	@Override
	public void run(){
		if(setLock!=null){
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

		if(launchInThread){new LaunchThread().start();}

		runEpochs();

		parent.networkQueue.add(net0);///TODO: Ensure net's stats are correct
		if(launchInThread){launchQueue.add(JobData.POISON);}
		if(setLock!=null){
//			System.err.println("Final Unlock");
			setLock.writeLock().unlock();
		}
	}

	/** Advances epochs while below the cap and bestErrorRate exceeds targetError.
	 * This class never updates bestErrorRate, so that test is not a measured convergence
	 * check. Optionally dumps samples, trains, schedules validation on print epochs,
	 * and updates the parent best-network snapshot after validation.
	 * @return number of epochs completed */
	private int runEpochs(){
		while(currentEpoch<maxEpochs && (bestErrorRate>targetError)){
			mprof.reset();
			currentEpoch++;

			if(currentEpoch==dumpEpoch && dumpRate>0){
				dump(data);
			}

			if(training){
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

			if(validated){
				calcNetStats(!training || Trainer.setCutoffForEvaluation);
				parent.compareWithBest(net0.copy(false));
			}

			mprof.log();//11: 11078/10597
			//System.err.println("M finished epoch "+currentEpoch);
		}
		return currentEpoch;
	}

	/** Legacy dump path that independently retains positive and negative high-error
	 * samples, optionally keeps a partial tail, shuffles the survivors, and rebuilds
	 * subsets.
	 * @param data sample set to reduce
	 * @deprecated use {@link #dump(SampleSet)} */
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

		for(Sample s : data.samples){
			s.setPivot();//Necessary because the assertion failed once.  Usually works, though.
			assert(s.checkPivot()) : s.pivot+", "+s.calcPivot(); //TODO: Slow
			if(s.positive){plist.add(s);}
			else{nlist.add(s);}
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
		if(parent.partialDumpFraction<1){
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
		if(parent.partialDumpFraction<1){
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
		Shared.shuffle(list, shared.Shared.random(SampleSet.shuffleSeed+1));
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
			if(setsize>0){
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

	/** Retains high-error positive and negative samples, optionally keeps partial
	 * low-error tails, shuffles the survivors, rebuilds subsets, and adjusts the
	 * fraction-per-epoch multiplier.
	 * @param data sample set to reduce */
	private void dump(SampleSet data){
		final float retainFraction=(1-dumpRate);//Fraction completely retained
		final float retainFraction2=(parent.partialDumpFraction<1 ? 1-parent.partialDumpFraction*dumpRate : retainFraction);//Total fraction retained, including the partials
		runTestingInterval(data.samples);
		ArrayList<Sample> plist=new ArrayList<Sample>();
		ArrayList<Sample> nlist=new ArrayList<Sample>();

//		System.err.println("len="+data.samples.length+", ret="+retainFraction+", +ret2="+retainFraction2);

		for(Sample s : data.samples){
			s.setPivot();//Necessary because the assertion failed once.  Usually works, though.
			assert(s.checkPivot()) : s.pivot+", "+s.calcPivot(); //TODO: Slow
			if(s.positive){plist.add(s);}
			else{nlist.add(s);}
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
		Shared.shuffle(list, shared.Shared.random(SampleSet.shuffleSeed+1));
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
			if(setsize>0){
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

	/** Sorts one polarity by pivot, copies the retained prefix, and fractionally
	 * samples the low-error tail when configured.
	 * @param inList samples to sort and select from
	 * @param outList destination for retained samples
	 * @param retainCount prefix length to retain unconditionally */
	private void dumpList(ArrayList<Sample> inList, ArrayList<Sample> outList, int retainCount){
		Collections.sort(inList);
		for(int i=0; i<retainCount; i++){
			outList.add(inList.get(i));
			assert(i==0 || inList.get(i).pivot<=inList.get(i-1).pivot);
		}

		if(parent.partialDumpFraction<1){
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

	/** Selects the current training subset, launches worker jobs with backpropagation
	 * enabled, gathers worker deltas, applies accumulated changes, and updates annealing,
	 * alpha, and triage state. */
	private void runTrainingInterval(){
//		synchronized(LOCK){

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

		synchronized(net0){
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

	/** Returns the early-epoch square-root weight ramp, or 1 after minWeightEpoch.
	 * @return weight multiplier for the current epoch */
	private final float weightMult(){
		if(currentEpoch>=minWeightEpoch){return 1.0f;}
		return (float)Math.sqrt((currentEpoch+1)*minWeightEpochInverse);
	}

	/** Evaluates up to maxLinesV samples without backpropagation and merges their
	 * validation statistics.
	 * @param set validation sample array */
	private void runTestingInterval(Sample[] set){
		final int vlines=Tools.min(parent.maxLinesV, set.length);
//		synchronized(LOCK){
			clearStats();
			net0.clear();
			int jobs=launchJobs(net0, set, vlines, false, 1f, false);
			mprof.log();//8: 1738/1709
//		}

		gatherResults(null, jobResultsQueue, false, jobs);
		lock();
		//System.err.println("M done waiting for threads.");

//		synchronized(LOCK){
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

	/** Releases one read hold and acquires the write hold when setLock is enabled.
	 * Caller must hold the read lock; this is not an atomic mode conversion. */
	void lock(){
//		System.err.println("Lock");
		if(setLock!=null){
			setLock.readLock().unlock();
			setLock.writeLock().lock();
		}
	}

	/** Releases one write hold and acquires the read hold when setLock is enabled.
	 * Caller must own the write lock; this is not an atomic mode conversion. */
	void unlock(){
//		System.err.println("Unlock");
		if(setLock!=null){
			setLock.writeLock().unlock();
			setLock.readLock().lock();
		}
	}

	/** Resolves this batch's SIMD feed-forward choice from simdFFMode (0=off, 1=on, 2=random-per-batch). */
	private boolean computeSimdFF(){
		if(simdFFMode==1){return true;}
		if(simdFFMode==2){return randySimd.nextBoolean();}
		return false;
	}

	/*--------------------------------------------------------------*/

	/** Optionally anneals the trained network, then decays annealing strength and
	 * slows the decay when the initial dropoff has nearly exhausted the strength. */
	private void anneal(){
		if(currentEpoch>maxAnnealEpoch){annealStrength=0;}
		else if(currentEpoch>=minAnnealEpoch && annealStrength>0 &&
				annealProb>0 && currentEpoch+1<maxEpochs){
			if(annealProb>randyAnneal.nextFloat()){
				net0.anneal((float)annealStrength, randyAnneal);
				//annealDropoff=annealDropoff*0.999f;
			}
			annealStrength=annealStrength*annealDropoff;

			if(annealDropoff0==annealDropoff && annealStrength*40<annealStrength0){
				annealDropoff=(1-(1-annealDropoff)*0.25f);//Slow anneal dropoff
			}
		}
	}

	/** Raises alpha through peakAlphaEpoch, then applies exponential decay. */
	private void adjustAlpha(){
		if(currentEpoch<=peakAlphaEpoch){
			alpha+=alphaIncrease;
		}else{
			alpha*=alphaDropoff;
		}
	}

	/** Triage-removes easy positives and negatives only after startTriage and only
	 * after a full-subset epoch. */
	private void triage(){//Do this AFTER processing the epoch
		if(currentEpoch>=startTriage && samplesThisEpoch==currentSamples.length){
			currentSubset.triage(currentEpoch, startTriage, positiveTriage, negativeTriage);
		}
	}

	/** Converts aggregate launch requests into per-worker jobs. */
	private class LaunchThread extends Thread{

		//Called by start()
		@Override
		public void run(){
			for(JobData job=getJob(); job!=JobData.POISON; job=getJob()){
				launchJobsInner(job.immutableNet, job.set, job.maxSamples, job.epoch, job.alpha,
						job.backprop, job.weightMult, job.sort);
			}
		}

		/** Blocks until a launch request is available, reporting interruptions and
		 * retrying.
		 * @return next launch request, including poison jobs */
		JobData getJob(){
			JobData job=null;
			while(job==null){
				try{
					job=launchQueue.take();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			return job;
		}

	}

	/** Launches worker jobs directly or delegates launch to {@link LaunchThread}.
	 * @param net0 network to copy or share with workers
	 * @param set samples to distribute
	 * @param numSamples requested number of samples to process
	 * @param backprop whether workers run network backpropagation
	 * @param weightMult weight multiplier passed to workers
	 * @param sort whether workers should sort extra samples
	 * @return number of worker results expected */
	private int launchJobs(CellNet net0, Sample[] set, int numSamples,
			boolean backprop, float weightMult, boolean sort){
		if(launchInThread && setLock==null){
			JobData job=new JobData(net0, jobResultsQueue, currentEpoch, numSamples, alpha,
					backprop, weightMult, sort, true, null, set, setLock, 0, 0);
			launchQueue.add(job);
			return Tools.min(jobsPerEpoch, numSamples, set.length);
		}else{
			return launchJobsInner(net0, set, numSamples, currentEpoch, alpha, backprop, weightMult, sort);
		}
	}

	//Does not seem faster...
	/** Launches list-backed jobs when setLock is disabled and delegates to the
	 * set-backed launcher otherwise.
	 * @param net0 network to copy or refresh from for workers
	 * @param set samples to distribute
	 * @param numSamples_ requested number of samples to process
	 * @param epoch current training epoch
	 * @param alpha learning rate for this job
	 * @param backprop whether workers run network backpropagation
	 * @param weightMult weight multiplier passed to workers
	 * @param sort whether workers should receive extra sortable samples
	 * @return number of worker jobs launched */
	private int launchJobsInner(CellNet net0, Sample[] set, int numSamples_, int epoch, double alpha,
			boolean backprop, float weightMult, boolean sort){
		if(setLock!=null){return launchJobs_SetLock(net0, set, numSamples_, epoch, alpha, backprop, weightMult, sort);}
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
		for(int jid=0; jid<jobsPerEpoch && jid<numSamples; jid++){
			ArrayList<Sample> list=new ArrayList<Sample>(listLen);
			int idx=jid;
			while(idx<numSamples){
				list.add(set[idx]);
				idx+=jobsPerEpoch;
			}
			final int toProcess=list.size();
			if(sortFlag){
				while(idx<samplesToSend){
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
				if(Trainer.setNetInWorkerThread){
					job=new JobData(Trainer.setNetInWorkerThread ? immutableNet : null, jobResultsQueue, epoch, toProcess, alpha,
							backprop, weightMult, sort, false, list, null, setLock, jid, jobsPerEpoch);
				}else{
					synchronized(net0){
						synchronized(subnets[jid]){
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

//		if(sortFlag & shuffleSubset && ((epoch&7)==3)){
//			currentSubset.shuffle();
//		}
		return jobs;
	}

	/** Launches set-backed jobs that stride through the same sample array under
	 * the optional sample-set lock.
	 * @param net0 network to copy or refresh from
	 * @param set samples to distribute by stride
	 * @param numSamples_ requested number of samples to process
	 * @param epoch current training epoch
	 * @param alpha learning rate for this job
	 * @param backprop whether workers run network backpropagation
	 * @param weightMult weight multiplier passed to workers
	 * @param sort ignored by this set-backed path
	 * @return number of worker jobs launched */
	private int launchJobs_SetLock(CellNet net0, Sample[] set, int numSamples_, int epoch, double alpha,
			boolean backprop, float weightMult, boolean sort){

		final int numSamples=Tools.min(numSamples_, set.length);

//		assert(!sort);
		int sent=0;
		int jobs=0;

		//This does not seem to change anything...
//		final CellNet immutableNet=Trainer.copyNetInWorkerThread ||  Trainer.setNetInWorkerThread ? net0.copy(false) : null;
		final CellNet immutableNet=Trainer.copyNetInWorkerThread ||  Trainer.setNetInWorkerThread ? net00.setFrom(net0, false) : null;
		final boolean simdFFThis=computeSimdFF(); //Per-batch SIMD_FF choice for this launch
		unlock();
		for(int jid=0; jid<jobsPerEpoch; jid++){

			final int toProcess=(numSamples-jid+jobsPerEpoch-1)/jobsPerEpoch;//I think this is right
			final JobData job;
			if(Trainer.copyNetInWorkerThread){
				job=new JobData(immutableNet, jobResultsQueue, epoch, numSamples, alpha,
						backprop, weightMult, false, true, null, set, setLock, jid, jobsPerEpoch);
			}else{
				if(Trainer.setNetInWorkerThread){
					job=new JobData(Trainer.setNetInWorkerThread ? immutableNet : null, jobResultsQueue, epoch, numSamples, alpha,
							backprop, weightMult, false, false, null, set, setLock, jid, jobsPerEpoch);
				}else{
					synchronized(net0){
						synchronized(subnets[jid]){
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
			job.simdFF=simdFFThis; //Per-batch SIMD feed-forward flag (off/on/random per simdFFMode)
			workerQueue.add(job);
		}

		assert(sent==numSamples && jobs==jobsPerEpoch) :
			"sort="+sort+", sent="+sent+", samples="+numSamples+", samples_="+numSamples_+", "+
			", setlen="+set.length+", jobs="+jobsPerEpoch;
		return jobs;
	}

	/** Gathers the expected worker results in deterministic or arrival order.
	 * @param net0 network to accumulate deltas into when requested
	 * @param mq queue containing worker results
	 * @param accumulate whether to accumulate worker network deltas
	 * @param numJobs exact number of results to consume */
	private void gatherResults(final CellNet net0, final ArrayBlockingQueue<JobResults> mq,
			final boolean accumulate, final int numJobs){
		if(orderedJobs){
			gatherResultsOrdered(net0, mq, accumulate, numJobs);
		}else{
			gatherResultsDisordered(net0, mq, accumulate, numJobs);
		}
	}

	/** Consumes worker results as they arrive and optionally accumulates deltas,
	 * allowing slight nondeterminism from arrival order.
	 * @param net0 network to accumulate deltas into when requested
	 * @param mq queue containing worker results
	 * @param accumulate whether to accumulate worker network deltas
	 * @param numJobs exact number of results to consume */
	private void gatherResultsDisordered(final CellNet net0, final ArrayBlockingQueue<JobResults> mq,
			final boolean accumulate, final int numJobs){
		//System.err.println("M waiting for threads.");
		for(int i=0; i<numJobs; i++){
			JobResults job=null;
			while(job==null){
				try{job=mq.take();}
				catch(InterruptedException e){e.printStackTrace();}
			}
			assert(job.epoch==currentEpoch) : job.epoch+", "+currentEpoch+", "+job.tid;
			gatherStats(job);
			if(accumulate && job.net!=null){net0.accumulate(job.net);}
			else{assert(!accumulate || jobsPerEpoch>samplesThisEpoch);}
		}
	}

	/** Buffers worker results by job ID and processes the consecutive prefix so
	 * accumulation follows submission order.
	 * @param net0 network to accumulate deltas into when requested
	 * @param mq queue containing worker results
	 * @param accumulate whether to accumulate worker network deltas
	 * @param numJobs exact number of results to consume */
	private void gatherResultsOrdered(final CellNet net0, final ArrayBlockingQueue<JobResults> mq,
			final boolean accumulate, final int numJobs){
		JobResults[] results=new JobResults[numJobs];
		//System.err.println("M waiting for threads.");

		int next=0;
		for(int i=0; i<numJobs; i++){
			{//Get a job
				JobResults job=null;
				while(job==null){
					try{job=mq.take();}
					catch(InterruptedException e){e.printStackTrace();}
				}
				assert(job.epoch==currentEpoch) : job.epoch+", "+currentEpoch+", "+job.tid;
				results[job.jid]=job;
			}

			//Process as many consecutive jobs as are available
			//Can be moved outside of the loop to ensure read-write exclusion on net0, if needed
			while(next<numJobs && results[next]!=null){//This loop actually seems to take very little time
				final JobResults job=results[next];
				gatherStats(job);
				if(accumulate && job.net!=null){net0.accumulate(job.net);}
				else{assert(!accumulate || jobsPerEpoch>samplesThisEpoch);}
				next++;
			}
		}
	}

	/** Advances the legacy print schedule used here to select validation epochs.
	 * This method does not print.
	 * @return true for evaluation-only mode, a scheduled epoch, or the final epoch */
	private boolean handlePrintInterval(){
		boolean print=!training || currentEpoch==maxEpochs;
		if(/*printStatus && */currentEpoch>=nextPrintEpoch){
			print=true;
			if(currentEpoch<printInterval){
				nextPrintEpoch=nextPrintEpoch*4;
				if(nextPrintEpoch>printInterval){
					nextPrintEpoch=printInterval;
				}
			}
			if(currentEpoch>=nextPrintEpoch){
				nextPrintEpoch+=printInterval;
			}
			nextPrintEpoch=Tools.min(nextPrintEpoch, maxEpochs);
		}
		return print;
	}

	/** Computes the current fraction-per-epoch ramp and dump multiplier.
	 * @return fraction of subset samples to process this epoch */
	float calcFractionPerEpoch(){
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
	/** Computes samples for this epoch, forcing scheduled full passes and otherwise
	 * applying the same fraction ramp used by calcFractionPerEpoch.
	 * @param currentSubset subset being processed this epoch
	 * @return number of samples to process from the subset */
	int calcSamplesThisEpoch(Subset currentSubset){
		final int len=currentSubset.samples.length;
		if(currentEpoch>=currentSubset.nextFullPassEpoch){//This should really go outside the function.
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

	/** Selects the epoch subset, computes samplesThisEpoch, then shuffles on
	 * Trainer.SHUFFLEMOD, fully sorts every 64 epochs at offset 5, or runs the
	 * partial-sort cadence for other eligible odd epochs. */
	private void selectTrainingSubset(){
		currentSubset=data.currentSubset(currentEpoch);
		currentSamples=currentSubset.samples;
		samplesThisEpoch=calcSamplesThisEpoch(currentSubset);
		assert(setLock==null || setLock.writeLock().isHeldByCurrentThread());
		final int mod8=currentEpoch&7, mod64=currentEpoch&63;
		if(shuffleSubset && mod8==Trainer.SHUFFLEMOD){
			currentSubset.shuffle();
			return;
		}else if(!sort || sortInThread){
			return;
		}else if(mod64==5){
			currentSubset.sortSamples(1f, allowMTSort);
		}else if(sortAll || mod8==5){
//			currentSubset.sortSamples(1f, allowMTSort);
			currentSubset.sortSamples2(1f, samplesThisEpoch, allowMTSort, flist);
		}else if(mod8==1){
			currentSubset.sortSamples2(fractionPerEpoch0*6, samplesThisEpoch, allowMTSort, flist);
		}else if(mod8==3 || mod8==7){
//			currentSubset.sortSamples(fractionPerEpoch0*3, allowMTSort);
			currentSubset.sortSamples2(fractionPerEpoch0*3, samplesThisEpoch, allowMTSort, flist);
		}
	}

	/** Clears interval error sums and confusion counts. */
	private void clearStats(){
		rawErrorSum=0;
		weightedErrorSum=0;
		tpSum=tnSum=fpSum=fnSum=0;
	}

	/** Adds one worker result's errors and confusion counts to interval totals.
	 * @param job completed worker result */
	private void gatherStats(JobResults job){
		assert(job.errorSum>=0) : job.errorSum;
		assert(job.weightedErrorSum>=0) : job.weightedErrorSum;
		rawErrorSum+=job.errorSum;
		weightedErrorSum+=job.weightedErrorSum;
		tpSum+=job.tpSum;
		tnSum+=job.tnSum;
		fpSum+=job.fpSum;
		fnSum+=job.fnSum;
	}

	/** Converts accumulated interval totals into rates and writes them to net0.
	 * @param samples sample count used to normalize error sums */
	private void mergeStats(int samples){
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

	/** Selects or retains a cutoff from validation results and configured FPR,
	 * FNR, or crossover targets, then writes the resulting rates to net0.
	 * @param retainOldCutoff whether to score the existing network cutoff */
	void calcNetStats(boolean retainOldCutoff){

		SampleSet set=validateSet;
		if(validateThisEpoch){
			if(retainOldCutoff){
				set.sortByValue();
				lastCutoff=net0.cutoff;
				fpRate=set.calcFPRFromCutoff(lastCutoff);
				fnRate=set.calcFNRFromCutoff(lastCutoff);
			}else if(crossoverFpMult>0){
				set.sortByValue();
				lastCutoff=set.calcCutoffFromCrossover(crossoverFpMult);
				fpRate=set.calcFPRFromCutoff(lastCutoff);
				fnRate=set.calcFNRFromCutoff(lastCutoff);
			}else if(targetFPR>=0){
				set.sortByValue();
				fpRate=targetFPR;
				fnRate=set.calcFNRFromFPR(targetFPR);
				lastCutoff=set.calcCutoffFromFPR(fpRate);
				fpRate=set.calcFPRFromCutoff(lastCutoff);
			}else if(targetFNR>=0){
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
		}else{
			//Not sure what to do here if retainOldCutoff=true
			lastCutoff=Trainer.cutoffForEvaluation;
		}

//		assert(false) : crossoverFpMult+", "+lastCutoff+", "+fpRate+", "+fnRate+", "+targetFPR+", "+set.numPositive+", "+set.numNegative;
//		assert(fnRate<1) : fnRate+", "+targetFNR+", "+targetFPR;//This was added because one time I forgot to include positive samples
		setNetStats(net0);
	}

	/** Writes current rates, alpha, anneal strength, epoch, and cutoff onto a net.
	 * @param net network to update */
	private void setNetStats(CellNet net){
		net.fpRate=(float) fpRate;
		net.fnRate=(float) fnRate;
		net.tpRate=(float) tpRate;
		net.tnRate=(float) tnRate;
		net.errorRate=(float) rawErrorRate;
		net.weightedErrorRate=(float) weightedErrorRate;
		net.alpha=(float) alpha;
		net.annealStrength=(float) annealStrength;
		net.epoch=currentEpoch;
		if(lastCutoff!=999){net.setCutoff((float)lastCutoff);}
	}

	/*--------------------------------------------------------------*/

	/** Reports the public error flag, which this class does not set.
	 * @return true if {@link #errorState} is false */
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

	/** Optional lock guarding set-backed worker jobs. */
	private final ReentrantReadWriteLock setLock;

	/*--------------------------------------------------------------*/

//	private final long annealSeed;
	/** Number of parallel jobs to launch per training epoch */
	private final int jobsPerEpoch;
//	private final int jobsPerBatch; //TODO: Change threads to this.

	/** Whether to process job results in submission order for determinism */
	private final boolean orderedJobs; //Without ordered, very very slight nondeterminism.
	/** SIMD feed-forward mode: 0=off, 1=on, 2=random-per-batch (set from Trainer.simdFFMode) */
	private final int simdFFMode;
	/** Per-net RNG for simdFFMode==2 (random per-batch SIMD_FF choice) */
	private final shared.Random randySimd;
	/** Queue for collecting results from completed worker jobs */
	private final ArrayBlockingQueue<JobResults> jobResultsQueue;
	/** Queue for distributing jobs to worker threads */
	private final ArrayBlockingQueue<JobData> workerQueue;
	/** Queue for aggregate requests to the optional launch thread. */
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

	/** Whether remaining eligible epochs use sortSamples2 with a fraction of 1. */
	final boolean sortAll;
	/** Whether to enable sample sorting */
	final boolean sort;
	/** Whether to perform sorting within worker threads */
	final boolean sortInThread;
	/** Whether to shuffle training subsets at the configured cadence. */
	final boolean shuffleSubset;
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

	/** Copied print-status option; unused by this thread. */
	private final boolean printStatus;
	/** Interval used by the legacy print schedule to select validation epochs */
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

	/** Reserved error-rate field; TrainerThread does not update it. */
	float bestErrorRate=999;
	/** Reserved false-negative-rate field; TrainerThread does not update it. */
	float bestFNR=999;

	/** Sum of raw sample errors for the current interval. */
	double rawErrorSum=0;
	/** Sum of weighted sample errors for the current interval. */
	double weightedErrorSum=0;
	/** Confusion counts accumulated for the current interval. */
	long tpSum=0, tnSum=0, fpSum=0, fnSum=0;

	/** Raw error rate calculated from accumulated errors */
	double rawErrorRate=999f;
	/** Weighted error rate calculated from accumulated weighted errors */
	double weightedErrorRate=999f;

	/** Confusion interval rates. */
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

	/** Next scheduled validation epoch under the legacy print schedule */
	private int nextPrintEpoch=1;

	/** Number of samples selected for processing in current epoch */
	private int samplesThisEpoch=-1;
	/** Whether validation should be performed in current epoch */
	private boolean validateThisEpoch=false;
	/** Epoch counter: zero before start, one-based during training. */
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

	/** Unused status output stream placeholder. */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;

}
