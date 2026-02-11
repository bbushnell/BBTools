package stream;

/**
 * Coordinates ordered processing with unordered worker threads,
 * accepting UNORDERED input and providing ORDERED output.
 * * This class is designed to follow a "scatter" process (like a
 * multithreaded Streamer) that produces results out-of-order.
 * It uses a JobQueue on the input to re-order the incoming
 * out-of-order jobs before distributing them to its own
 * internal worker threads, and a second JobQueue to re-order
 * their output.
 * * Producer (unordered) → JobQueueIn (orders) → Workers (unordered) → JobQueueOut (orders) → Consumer (ordered)
 * * @author Brian Bushnell
 * @contributor Gemini/Isla
 * @date November 16, 2025
 *
 * @param <I> Input job type (must implement HasID)
 * @param <O> Output job type (must implement HasID)
 */
public class OrderedQueueSystem2<I extends HasID, O extends HasID> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * @param numWorkers Number of internal worker threads this queue will feed
	 * @param orderedOutput If true, output is ordered; if false, output is unordered
	 * @param inputPrototype Prototype for creating input poison pills
	 * @param outputPrototype Prototype for creating output poison/last pills
	 */
	public OrderedQueueSystem2(int numWorkers, boolean orderedOutput, I inputPrototype_, O outputPrototype_){
		this(numWorkers+BUFFER_PADDING, (BUFFER_MULT*numWorkers)/2+BUFFER_PADDING, 
			numWorkers, orderedOutput, inputPrototype_, outputPrototype_);
	}

	/**
	 * @param capacityIn Capacity of the input re-ordering queue
	 * @param capacityOut Capacity of the output re-ordering queue
	 * @param numWorkers Number of internal worker threads this queue will feed
	 * @param orderedOutput If true, output is ordered; if false, output is unordered
	 * @param inputPrototype Prototype for creating input poison pills
	 * @param outputPrototype Prototype for creating output poison/last pills
	 */
	public OrderedQueueSystem2(int capacityIn, int capacityOut, int numWorkers_,
			boolean orderedOutput, I inputPrototype_, O outputPrototype_){
		// The input queue MUST be ordered to re-sort the incoming unordered data
		inq=new JobQueue<I>(capacityIn, true, true, 0); 
		outq=new JobQueue<O>(capacityOut, orderedOutput, true, 0);
		inputPrototype=inputPrototype_;
		outputPrototype=outputPrototype_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Producer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Add input job for processing.
	 * This is thread-safe and accepts out-of-order jobs.
	 */
	public boolean addInput(I job){
		if(job==null){return false;}
		assert(!job.last()) : "Use poison() to terminate";
		synchronized(this){
			final long id=job.id();
			assert(!lastSeen || job.poison());
			assert(id!=prevID || job.poison() || job.last()) : "Nonunique ID: "+prevID;
			prevID=id;
			maxSeenId=Math.max(id, maxSeenId);
		}
		// JobQueue.add() is blocking and handles its own wait/interrupt
		return inq.add(job);
	}

	/** Signal end of input - injects LAST to output and POISON to input. */
	@SuppressWarnings("unchecked")
	public synchronized void poison(){
		if(lastSeen){return;}
		if(verbose) {System.err.println("OQS2: poison()");}
		
		final long finalID = maxSeenId + 1;

		// Add ONE lastJob for the final consumer
		O lastJob=(O)outputPrototype.makeLast(finalID);
		outq.add(lastJob);

		// Add ONE poison pill for the worker threads.
		// The first worker to get it will exit its loop and re-inject it.
		I poisonJob=(I)inputPrototype.makePoison(finalID);
		inq.add(poisonJob);

		lastSeen=true;
	}

	/** Wait for processing to complete. */
	public synchronized void waitForFinish(){
		while(!finished){
			try{this.wait();}
			catch(InterruptedException e){e.printStackTrace();}
		}
	}

	/** Convenience: poison and wait. */
	public void poisonAndWait(){
		poison();
		waitForFinish();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker API           ----------------*/
	/*--------------------------------------------------------------*/

	/** Get next input job in-order (blocks). */
	public I getInput(){
		I job=inq.take(); // Pulls from the re-ordering input queue
		if(verbose) {System.err.println("OQS2: getInput I "+job.id()+": "+job.poison()+", "+job.last());}
		return job;
	}

	/** Add processed output job. */
	public void addOutput(O job){
		if(verbose) {System.err.println("OQS2: addOutput O "+job.id()+": "+job.poison()+", "+job.last());}
		outq.add(job);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Consumer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** Check if more output is coming. */
	public boolean hasMore(){
		return outq.hasMore();
	}

	/** Get next output job in order (blocks). */
	public O getOutput(){
		return outq.take();
	}

	/** Signal that processing is complete from the consumer side. */
	public synchronized void setFinished(boolean force){
		if(verbose) {System.err.println("OQS2: setFinished()");}
		finished=true;
		inq.poison((I)inputPrototype.makePoison(maxSeenId+1), force);  // Tell input queue to stop blocking any waiting producers
		outq.poison((O)outputPrototype.makePoison(maxSeenId+1), force); // Tell output queue to stop blocking any waiting workers
		this.notifyAll();
	}
	
	public synchronized boolean finished() {return finished;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final JobQueue<I> inq;
	private final JobQueue<O> outq;
	private final I inputPrototype;
	private final O outputPrototype;

	private long maxSeenId=-1;
	private long prevID=-1;
	private volatile boolean finished=false;
	private volatile boolean lastSeen=false;
	private static final boolean verbose=false;

	public static int BUFFER_PADDING=4;
	public static int BUFFER_MULT=3;

}