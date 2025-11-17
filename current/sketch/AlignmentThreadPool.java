package sketch;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import shared.Tools;

/**
 * Thread pool for processing Comparison objects that require SSU identity calculations.
 * Manages a pool of worker threads to perform alignment computations in parallel.
 * Uses a producer-consumer pattern with blocking queues for thread-safe job distribution.
 *
 * @author Brian Bushnell
 * @date July 2025
 */
public class AlignmentThreadPool {
	
	/** Creates a new thread pool with the specified maximum thread count.
	 * @param maxThreads_ Maximum number of worker threads to spawn */
	public AlignmentThreadPool(int maxThreads_) {
		maxThreads=maxThreads_;
		assert(maxThreads>0);
		tlist=new ArrayList<AlignmentThread>(maxThreads);
	}
	
	/**
	 * Processes a list of Comparison objects, submitting those needing alignment to worker threads.
	 * Only submits comparisons that return true from needsAlignment().
	 * Blocks until all submitted jobs are completed.
	 *
	 * @param list List of Comparison objects to process
	 * @param maxRecords Maximum number of records to process from the list
	 */
	public void addJobs(ArrayList<Comparison> list, int maxRecords){
		if(list==null || list.isEmpty() || maxRecords<1){return;}
		final int limit=Tools.min(list.size(), maxRecords);
		ArrayBlockingQueue<Comparison> dest=new ArrayBlockingQueue<Comparison>(limit);
		int added=0;
		for(int i=0; i<limit; i++){
			Comparison c=list.get(i);
			if(c.needsAlignment()){
				addJob(c, dest);
				added++;
			}
		}
		for(int i=0; i<added; i++){
			take(dest);
		}
	}
	
	/**
	 * Submits a single Comparison for alignment processing.
	 * Creates a new worker thread if needed and threads are available.
	 * @param c The Comparison object requiring alignment
	 * @param dest Queue to receive the processed Comparison
	 */
	public void addJob(Comparison c, ArrayBlockingQueue<Comparison> dest){
		if(tlist.size()<maxThreads){spawnThread();}
		assert(!poisoned);
		AlignmentJob job=new AlignmentJob(c, dest);
		put(job);
	}
	
	/** Creates a new worker thread if the current thread count is below maximum
	 * and all existing threads are busy. Thread-safe creation with synchronization. */
	private synchronized void spawnThread(){
		final int size=tlist.size();
		if(size<maxThreads && busy.get()>=size){
//			AlignmentThread alt=new AlignmentThread(source, busy);
			AlignmentThread alt=new AlignmentThread();
			tlist.add(alt);
			alt.start();
		}
	}
	
	/** Initiates shutdown of the thread pool by adding a poison pill to the job queue.
	 * Worker threads will terminate when they encounter the poison job. */
	synchronized void poison(){
		assert(!poisoned);
		if(poisoned){return;}
		put(poison);
		poisoned=true;
	}
	
	/**
	 * Adds a job to the shared work queue with busy counter increment.
	 * Blocks until the job can be successfully queued, handling interruptions.
	 * @param job The AlignmentJob to add to the queue
	 */
	private void put(AlignmentJob job){
		busy.incrementAndGet();
		boolean success=false;
		while(!success){
			try {
				source.put(job);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Retrieves an item from the specified blocking queue, handling interruptions.
	 * Blocks until an item becomes available.
	 *
	 * @param <X> Type of items in the queue
	 * @param queue The queue to retrieve from
	 * @return The retrieved item
	 */
	private final <X> X take(ArrayBlockingQueue<X> queue){
		X x=null;
		while(x==null){
			try {
				x=queue.take();
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		return x;
	}
	
	/** Worker thread that processes AlignmentJob objects from the shared queue.
	 * Continues processing until a poison pill is encountered, then propagates the poison. */
	private class AlignmentThread extends Thread {
		
		AlignmentThread(){}
		
		private final AlignmentJob next(){
			return take(source);
		}
		
		@Override
		public void run(){
			AlignmentJob job=null;
			for(job=next(); !job.isPoison(); job=next()){
				job.doWork();
				busy.decrementAndGet();
			}
			put(poison);
		}
		
//		private final ArrayBlockingQueue<AlignmentJob> source;
//		private final AtomicInteger busy;
//
//		private static final AlignmentJob poison=new AlignmentJob(null, null);
		
	}
	
	/** List of active worker threads in the pool */
	final ArrayList<AlignmentThread> tlist;
	/** Maximum number of worker threads allowed in the pool */
	final int maxThreads;
	/** Counter tracking the number of jobs currently being processed */
	final AtomicInteger busy=new AtomicInteger(0);
	/** Flag indicating whether the thread pool has been poisoned for shutdown */
	private boolean poisoned=false;

	/** Poison pill job used to signal worker threads to terminate */
	private static final AlignmentJob poison=new AlignmentJob(null, null);
	/** Shared queue containing AlignmentJob objects for worker thread processing */
	private static final ArrayBlockingQueue<AlignmentJob> source=new ArrayBlockingQueue<AlignmentJob>(4096);
	
}
