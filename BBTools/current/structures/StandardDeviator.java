package structures;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.Scaffold;
import shared.Timer;
import shared.Tools;
import template.Accumulator;
import template.ThreadWaiter;

public class StandardDeviator implements Accumulator<StandardDeviator.ProcessThread> {

	public StandardDeviator(boolean stranded_, int strand_) {
		stranded=stranded_;
		strand=strand_;
	}
	
	//This one is to slow, maybe because short contigs need too much locking.
	//They would have to be put into lists of contigs.
//	public boolean calculateStuffQueued(int threads0, Collection<Scaffold> scaffolds, 
//			boolean calcStd, boolean calcMedian, double windowAvg, int window) {
//
//		final int threads=Tools.mid(1, threads0, scaffolds.size());
//		Timer t=new Timer();
//		
//		System.err.println("threads="+threads);
//		
//		{//First pass to sum depth
//			totalSum=0;
//			bins=0;
//			ArrayBlockingQueue<Scaffold> queue=new ArrayBlockingQueue<Scaffold>(threads*2);
//			//Fill a list with ProcessThreads
//			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
//			for(int i=0; i<threads; i++){
//				alpt.add(new ProcessThread(i, queue, false, false, -1, -1, -1));
//			}//TODO: For some reason this is running singlethreaded... 97% CPU.
//			//Start the threads and wait for them to finish
//			ThreadWaiter.startThreads(alpt);
//			for(Scaffold s : scaffolds) {addToQueue(queue, s);}
//			addToQueue(queue, POISON);
//			boolean success=ThreadWaiter.waitForThreads(alpt, this);
//			errorState&=!success;
//			t.stopAndStart("calcSumMT:");
//		}
//		
//		final double globalMean=totalSum/(double)bins;
//		{//Second pass for statistics
//			totalSum=0;
//			bins=0;
//			ArrayBlockingQueue<Scaffold> queue=new ArrayBlockingQueue<Scaffold>(threads*2);
//			//Fill a list with ProcessThreads
//			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
//			for(int i=0; i<threads; i++){
//				alpt.add(new ProcessThread(i, queue, calcStd, calcMedian, globalMean, windowAvg, window));
//			}//This part only uses 150% CPU
//			//Start the threads and wait for them to finish
//			ThreadWaiter.startThreads(alpt);
//			for(Scaffold s : scaffolds) {addToQueue(queue, s);}
//			addToQueue(queue, POISON);
//			boolean success=ThreadWaiter.waitForThreads(alpt, this);
//			errorState&=!success;
//			t.stopAndStart("calcStdMT:");
//		}
//		
//		return errorState;
//	}
	
	public boolean calculateStuff(int threads0, ArrayList<Scaffold> scaffolds, 
			boolean calcStd, boolean calcMedian, boolean calcHist, 
			int minDepth, int histMax, double windowAvg, int window) {

		final int threads=Tools.mid(1, threads0, scaffolds.size());
		Timer t=new Timer();
		
//		if(verbose) {System.err.println("threads="+threads);}
		
		{//First pass to sum depth for global stdev
			totalSum=0;
			bins=0;
			//Fill a list with ProcessThreads
			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
			for(int i=0; i<threads; i++){
				alpt.add(new ProcessThread(i, threads, scaffolds, false, false, false, minDepth, -1, -1, -1, -1));
			}
			//Start the threads and wait for them to finish
			ThreadWaiter.startThreads(alpt);
			boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
			errorState&=!success;
			if(verbose) {t.stopAndStart("calcSumMT:");}
		}
		
		final double globalMean=totalSum/(double)bins;
		{//Second pass for statistics
			totalSum=0;
			bins=0;
			//Fill a list with ProcessThreads
			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
			for(int i=0; i<threads; i++){
				alpt.add(new ProcessThread(i, threads, scaffolds,
						calcStd, calcMedian, calcHist, minDepth, histMax, globalMean, windowAvg, window));
			}
			//Start the threads and wait for them to finish
			ThreadWaiter.startThreads(alpt);
			boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
			errorState&=!success;
			if(verbose) {t.stopAndStart("calcStdMT:");}
		}
		
		return errorState;
	}

	public double[] standardDeviation(ArrayList<Scaffold> scaffolds, int minscaf0){
		totalSum=bins=0;
		mean=stdev=0;
		long scafsCounted=0;
		final int minscaf=Tools.max(minscaf0, 1);
		
		for(Scaffold scaf : scaffolds){
			if(scaf.length>=minscaf){
				final CoverageArray ca=(CoverageArray)(strand==1 ? scaf.obj1 : scaf.obj0);
				bins+=scaf.length;
				totalSum+=(ca==null ? 0 : ca.sum());
				scafsCounted++;
			}
		}
		
		if(bins<1){return new double[] {mean, stdev};}
		mean=totalSum/(double)bins;
//		System.err.println("sum="+totalSum+", bins="+bins+", mean="+mean+", scaffolds="+scafsCounted);
		double sumdev2=0;
		for(Scaffold scaf : scaffolds){
			if(scaf.length>=minscaf){
				final CoverageArray ca=(CoverageArray)(strand==1 ? scaf.obj1 : scaf.obj0);
				if(ca==null) {
					sumdev2+=(mean*mean*scaf.length);
				}else {
					sumdev2+=ca.devSum(mean);
				}
			}
		}

		stdev=Math.sqrt(sumdev2/bins);
//		System.err.println("sumdev2="+sumdev2+", stdev="+stdev);
		return new double[] {mean, stdev};
	}

	/** This works but is slow.  Just use unbinned std. */
	@Deprecated
	public double[] standardDeviationBinned(ArrayList<Scaffold> scaffolds, int binsize, int minscaf0){
		totalSumBinned=binsBinned=0;
		meanBinned=stdevBinned=0;
		final int minscaf=Tools.max(minscaf0, 1);
		
		for(Scaffold scaf : scaffolds){
			if(scaf.length>=minscaf){
				CoverageArray ca=(CoverageArray)(strand==1 ? scaf.obj1 : scaf.obj0);
				int bins=(int)((scaf.length+(long)binsize-1L)/binsize);
				assert(bins>0) : scaf.length+", "+binsize+", "+minscaf;
				binsBinned+=bins;
				if(ca==null) {
					totalSumBinned+=0;
				}else {
					totalSumBinned+=ca.sum()/bins;
				}
			}
		}
		
		if(binsBinned<1){return new double[] {meanBinned, stdevBinned};}
		meanBinned=totalSumBinned/(double)binsBinned;
		double sumdev2=0;
		
		for(Scaffold scaf : scaffolds){
			if(scaf.length>=minscaf){
				CoverageArray ca=(CoverageArray)(strand==1 ? scaf.obj1 : scaf.obj0);
				int lastPos=-1, nextPos=binsize-1;
				long tempSum=0;
				final int lim=scaf.length-1;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					tempSum+=x;
					if(i>=nextPos || i==lim){
						int bin=(i-lastPos);
						double depth=(tempSum/(double)bin);
						double dev=meanBinned-depth;
						sumdev2+=(dev*dev);
						nextPos+=binsize;
						lastPos=i;
						tempSum=0;
					}
				}
			}
		}
		
		stdevBinned=Math.sqrt(sumdev2/binsBinned);
		return new double[] {meanBinned, stdevBinned};
	}

	//Unsafe because it will fail if there are over 2 billion bins
	public double[] standardDeviationBinnedUnsafe(String fname, ArrayList<Scaffold> scaffolds, int binsize, int minscaf0){
		final int minscaf=Tools.max(minscaf0, 1);
		
		LongList list=new LongList();
		for(Scaffold scaf : scaffolds){
			CoverageArray ca=(CoverageArray)(strand==1 ? scaf.obj1 : scaf.obj0);
			int lastPos=-1, nextPos=binsize-1;
			long sum=0;
			final int lim=scaf.length-1;
			for(int i=0; i<scaf.length; i++){
				int x=(ca==null ? 0 : ca.get(i));
				sum+=x;
				if(i>=nextPos || i==lim){
					int bin=(i-lastPos);
					if(scaf.length>=minscaf){
						list.add((int)(10*(sum/(double)bin)));
					}
					nextPos+=binsize;
					lastPos=i;
					sum=0;
				}
			}
		}
		list.sort();
		
		mean=0.1*list.mean();
		double median=0.1*list.median();
		double mode=0.1*list.mode();
		stdev=0.1*list.stdev();
		return new double[] {mean, median, mode, stdev};
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Threading           ----------------*/
	/*--------------------------------------------------------------*/
	
//	static class ProcessThreadQueued extends Thread {
//		
//		ProcessThreadQueued(int tid_, ArrayBlockingQueue<Scaffold> queue_, boolean calcStd_, 
//				boolean calcMedian_, double globalMean_, double windowAvg_, int window_){
//			tid=tid_;
//			queue=queue_;
//			
//			calcStd=calcStd_;
//			calcMedian=calcMedian_;
//			globalMean=globalMean_;
//			windowAvg=windowAvg_;
//			window=window_;
//		}
//		
//		@Override
//		public void run() {
////			synchronized(this) {
//				for(Scaffold next=getNext(queue); next!=POISON; next=getNext(queue)){
////					synchronized(next) {
//						binsT+=next.length;
//						processCA((CoverageArray)next.obj1);
//						processCA((CoverageArray)next.obj2);
////					}
//				}
//				addToQueue(queue, POISON);//Send it to the next one
////			}
//		}
//		
//		private void processCA(CoverageArray ca) {
//			if(ca==null) {return;}
//			totalSumT+=ca.sum();
//			if(calcStd) {ca.standardDeviation();}
//			if(calcMedian) {ca.median();}
//			if(globalMean>-1) {ca.devSum(globalMean);}
//			if(windowAvg>0) {ca.basesUnderAverageCoverage(windowAvg, window);}
//		}
//		
//		long totalSumT=0;
//		long binsT=0;
//		
//		final int tid;
//		final ArrayBlockingQueue<Scaffold> queue;
//		
//		final boolean calcStd;
//		final boolean calcMedian;
//		final double globalMean;
//		final double windowAvg;
//		final int window;
//	}
	
	class ProcessThread extends Thread {
		
		ProcessThread(int tid_, int threads_, ArrayList<Scaffold> list_, 
				boolean calcStd_, boolean calcMedian_, boolean calcHist_, int minDepth_, int histMax_, 
				double globalMean_, double windowAvg_, int window_){
//			System.err.println(calcHist_);
			tid=tid_;
			threads=threads_;
			list=list_;
			
			calcStd=calcStd_;
			calcMedian=calcMedian_;
			calcHist=calcHist_;
			minDepth=minDepth_;
			histMax=histMax_;
			globalMean=globalMean_;
			windowAvg=windowAvg_;
			window=window_;
			if(calcHist) {
				assert(histMax>0) : histMax;
				hist0=new LongList(Tools.min(histMax+1, Character.MAX_VALUE+1));
				if(stranded) {hist1=new LongList(Tools.min(histMax+1, Character.MAX_VALUE+1));}
			}
		}
		
		@Override
		public synchronized void run() {
			for(int i=tid; i<list.size(); i+=threads){
				Scaffold scaf=list.get(i);
				synchronized(scaf) {
					binsT+=scaf.length;
					processCA((CoverageArray)scaf.obj0, 0, scaf.length, scaf.basehits);
					if(stranded) {processCA((CoverageArray)scaf.obj1, 1, scaf.length, scaf.basehits);}
				}
			}
		}
		
		/** In the case of preallocated ca's with 0 mapped reads, this goes faster */
		private void processCA(CoverageArray ca, int strand, int length, long mappedBases) {
			if(mappedBases>0) {
				if(calcHist) {addToHist(ca, strand, length);}
				if(ca==null) {return;}
				totalSumT+=ca.calcSumAndCovered(minDepth);
				if(calcStd) {ca.standardDeviation();}
				if(calcMedian) {ca.median();}
				if(globalMean>-1) {ca.devSum(globalMean);}
				if(windowAvg>0) {ca.basesUnderAverageCoverage(windowAvg, window);}
			}else {
				if(calcHist) {addToHist(null, strand, length);}
				if(ca==null) {return;}
				totalSumT+=0;
				ca.setSum(0);
				ca.setCovered(0);
				ca.setStdev(0);
				ca.setMedian(0);
				ca.setDevSum(mappedBases);
				if(globalMean>-1) {ca.setDevSum((globalMean*globalMean)*length);}
				if(windowAvg>0) {ca.setUnderWindowAverage(length);}
			}
		}
		
		private void addToHist(CoverageArray ca, int strand, int length) {
			LongList hist=(strand==0 ? hist0 : hist1);
			assert(hist!=null) : stranded+", "+strand+", "+length;
			if(ca==null) {
				hist.increment(0, length);
				return;
			}
			assert(ca.length()==length) : ca.length()+", "+length+", "+ca.maxIndex;
			for(int i=0, len=ca.length(); i<len; i++) {
				int d=Tools.min(ca.get(i), histMax);
				hist.increment(d);
			}
		}
		
		long totalSumT=0;
		long binsT=0;
		
		final int tid;
		final int threads;
		final ArrayList<Scaffold> list;
		LongList hist0;
		LongList hist1;
		
		final boolean calcStd;
		final boolean calcMedian;
		final boolean calcHist;
		final int minDepth;
		final int histMax;
		final double globalMean;
		final double windowAvg;
		final int window;
	}
	
	private static Scaffold getNext(ArrayBlockingQueue<Scaffold> queue) {
		Scaffold next=null;
		while(next==null) {
			try {
				next=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return next;
	}
	
	private static void addToQueue(ArrayBlockingQueue<Scaffold> queue, Scaffold s) {
		while(s!=null) {
			try {
				queue.put(s);
				s=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Override
	public void accumulate(ProcessThread t) {
		synchronized(t) {
			totalSum+=t.totalSumT;
			bins+=t.binsT;
			if(hist0==null) {hist0=t.hist0;}
			else if(t.hist0!=null) {hist0.incrementBy(t.hist0);}
			if(hist1==null) {hist1=t.hist1;}
			else if(t.hist1!=null) {hist1.incrementBy(t.hist1);}
		}
	}

	@Override
	public boolean success() {
		// TODO Auto-generated method stub
		return false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long totalSum;
	long bins;
	double mean;
	double stdev;
	
	long totalSumBinned;
	long binsBinned;
	double meanBinned;
	double stdevBinned;

	public LongList hist0;
	public LongList hist1;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	public final int strand;
	public final boolean stranded;
	public boolean errorState=false;
	
	private static final Scaffold POISON=new Scaffold("POISON", null, -1);
	public static boolean verbose=false;
	
}
