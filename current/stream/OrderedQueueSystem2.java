package stream;

/**
 * Supplies input and output queues for caller-owned workers; creates no threads.
 * Input arrival and worker completion may be out of order, but ordinary job IDs
 * must be unique and dense from zero, remain stable while queued, and have one
 * corresponding output per input ID. The input queue requests ordering; the output
 * request is configurable, but the current JobQueue implementation forces ordering
 * even when orderedOutput is false. That flag does not permit gaps in output IDs.
 *
 * Callers coordinate input registration before poison(), workers consume getInput()
 * and submit addOutput(), and a consumer/owner signals setFinished(). This wrapper
 * does not process jobs, join workers or infer completion from thread state.
 *
 * @author Brian Bushnell
 * @contributor Gemini/Isla
 * @date November 16, 2025
 *
 * @param <I> Input job type (must implement HasID)
 * @param <O> Output job type (must implement HasID)
 */
public class OrderedQueueSystem2<I extends HasID, O extends HasID>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Derives input capacity as numWorkers+BUFFER_PADDING and output capacity as
	 * (BUFFER_MULT*numWorkers)/2+BUFFER_PADDING, using integer arithmetic.
	 * No workers are created; the current static sizing values are captured by construction.
	 * @param numWorkers Worker-count hint used only in the capacity formulas
	 * @param orderedOutput Requested output policy; current JobQueue forces ordering regardless
	 * @param inputPrototype_ Factory for type-compatible poison markers with queue-compatible IDs
	 * @param outputPrototype_ Factory for type-compatible output LAST and poison markers
	 */
	public OrderedQueueSystem2(int numWorkers, boolean orderedOutput, I inputPrototype_, O outputPrototype_){
		this(numWorkers+BUFFER_PADDING, (BUFFER_MULT*numWorkers)/2+BUFFER_PADDING,
			numWorkers, orderedOutput, inputPrototype_, outputPrototype_);
	}

	/**
	 * Creates two bounded JobQueues starting at ID zero and retains marker factories.
	 * Capacities control JobQueue backpressure rather than a hard retained-job count.
	 * Factories are not null-checked; their markers must match the respective job type
	 * and ordered terminal IDs when used. This constructor does not start or count workers.
	 * @param capacityIn Input backpressure setting; JobQueue asserts it is greater than one
	 * @param capacityOut Output backpressure setting; JobQueue asserts it is greater than one
	 * @param numWorkers_ Unused compatibility parameter when capacities are explicit
	 * @param orderedOutput Requested output policy; current JobQueue still forces ordering
	 * @param inputPrototype_ Retained factory for input poison markers
	 * @param outputPrototype_ Retained factory for output LAST and poison markers
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
	 * Registers and enqueues a borrowed input job, blocking according to JobQueue policy.
	 * Ordinary IDs may arrive out of order but must be globally unique and dense from zero.
	 * This wrapper checks only the immediately preceding registration for duplicate IDs;
	 * it does not prove global uniqueness. Complete ordinary registrations before poison().
	 * @param job Input job; null is ignored, and LAST is rejected by an assertion
	 * @return false for null, otherwise the input queue's add result
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
		// JobQueue.add() is blocking and handles its own wait/interrupt.
		//Historical ordering rationale: metadata is registered under this monitor before enqueue.
		//With dense IDs, compatible marker IDs and no new ordinary registrations after poison(),
		//already-registered jobs sort before the terminal even if their blocking adds finish later.
		//This is a caller/queue protocol requirement, not a general concurrency or completion guarantee.
		return inq.add(job);
	}

	/** Seals the registered input range once, then requests an output LAST and input poison.
	 * Both factories receive maxSeenId+1; the output marker is enqueued first.
	 * Blocking enqueues occur outside this monitor. Subsequent callers return once
	 * lastSeen is set, even while the first caller is still publishing its markers. */
	@SuppressWarnings("unchecked")
	public void poison(){
		//[OQS2 deadlock FIXED 2026-06-20 (greenlit)] - the SAME latent bug found+fixed in OrderedQueueSystem
		//(via stream/bam/BgzfInputStreamMT3#003): poison() was `synchronized(this)` and did the BLOCKING
		//JobQueue enqueues (outq.add + inq.add) WHILE HOLDING the monitor. If outq/inq was full, poison()
		//blocked holding `this`, so setFinished() (also synchronized(this)) could never run to poison the queues
		//and wake the blocked enqueue -> circular deadlock. REACHABLE in the live writers (BamWriter, FastqWriter,
		//SamWriter, BgzfOutputStreamMT2): on close() main calls poison(), and if the WRITER thread errors on
		//out.write() (full disk / broken pipe) its catch calls setFinished(true) -> deadlock against the blocked
		//poison(). FIX = mirror addInput's own contract (L74): update lastSeen/read maxSeenId UNDER the lock, do
		//the blocking enqueues OUTSIDE it. (addInput already did this right; poison() was just inconsistent.)
		final long finalID;
		synchronized(this){
			if(lastSeen){return;}
			lastSeen=true;
			finalID=maxSeenId+1;
		}
		if(verbose){System.err.println("OQS2: poison()");}

		// Add ONE lastJob for the final consumer (blocking, but OUTSIDE the lock now)
		outq.add((O)outputPrototype.makeLast(finalID));

		// Add ONE poison pill for the worker threads (blocking, outside the lock).
		//JobQueue retains/reinserts poison and returns null to workers; this wrapper owns no worker loop.
		inq.add((I)inputPrototype.makePoison(finalID));
	}

	/** Waits for the flag set by an external setFinished() call, not for worker joins.
	 * Interrupted waits print a trace and continue; this wrapper does not restore the interrupt. */
	public synchronized void waitForFinish(){
		while(!finished){
			try{this.wait();}catch(InterruptedException e){e.printStackTrace();}
		}
	}

	/** Seals input with poison(), then waits for an external completion notification. */
	public void poisonAndWait(){
		poison();
		waitForFinish();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker API           ----------------*/
	/*--------------------------------------------------------------*/

	/** Takes the next input in queue order, blocking until the queue can supply it.
	 * @return Borrowed job, or null when JobQueue reports a terminal condition */
	public I getInput(){
		I job=inq.take();//Pulls from the re-ordering input queue.
		//TODO: Probable bug - STR266: enabling verbose makes this diagnostic dereference a null terminal.
		if(verbose){System.err.println("OQS2: getInput I "+job.id()+": "+job.poison()+", "+job.last());}
		return job;
	}

	/** Enqueues a borrowed processed job, blocking according to output queue policy.
	 * Ordinary output IDs must cover the input range densely; normal final signaling
	 * is supplied by poison(). This method does not copy or validate the payload.
	 * @param job Nonnull processed job with a stable corresponding input ID */
	public void addOutput(O job){
		if(verbose){System.err.println("OQS2: addOutput O "+job.id()+": "+job.poison()+", "+job.last());}
		outq.add(job);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Consumer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** Delegates the output queue's availability indicator; does not wait for completion.
	 * @return Current JobQueue.hasMore() result; not a guarantee of a nonnull next take */
	public boolean hasMore(){return outq.hasMore();}

	/** Takes the next output using the current queue ordering policy.
	 * @return Borrowed job, including the LAST marker itself, or null for a queue terminal
	 * condition. Consumers must recognize LAST or accept its marker payload. */
	public O getOutput(){return outq.take();}

	/** Sets the completion flag, requests termination on both queues and notifies waiters.
	 * Called by the external consumer/owner; it does not join workers or verify their results.
	 * @param force Forwarded unchanged to each JobQueue.poison call */
	//TODO: Possible bug - see JobQueue.poison(): force=false does not guarantee blocked threads wake
	//if a stream has an ID gap (dead worker); error paths should pass force=true.  Also see
	//JobQueue.hasMore(): after force=true, hasMore() stays true while take() returns null.
	@SuppressWarnings("unchecked")
	public synchronized void setFinished(boolean force){
		if(verbose){System.err.println("OQS2: setFinished()");}
		finished=true;
		inq.poison((I)inputPrototype.makePoison(maxSeenId+1), force);//Request input termination under the selected policy.
		outq.poison((O)outputPrototype.makePoison(maxSeenId+1), force);//Request output termination under the selected policy.
		this.notifyAll();
	}

	/** @return Stored external-completion flag, not an independently measured worker status */
	public synchronized boolean finished(){return finished;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input queue, constructed with ordering requested and initial ID zero. */
	private final JobQueue<I> inq;
	/** Output queue; the current JobQueue implementation forces ordering. */
	private final JobQueue<O> outq;
	/** Retained input marker factory, not a queued ordinary job. */
	private final I inputPrototype;
	/** Retained output marker factory, used for LAST and poison signals. */
	private final O outputPrototype;

	/** Largest registered input ID, used to choose terminal IDs under this monitor. */
	private long maxSeenId=-1;
	/** Immediately preceding registered ID; not a set of all previously submitted IDs. */
	private long prevID=-1;
	/** Completion notification set by the external consumer/owner. */
	private volatile boolean finished=false;
	/** Whether poison() has sealed input registration; not a queue-drained indicator. */
	private volatile boolean lastSeen=false;
	/** Compile-time diagnostic switch; disabled for normal builds. */
	private static final boolean verbose=false;

	/** Extra capacity added by the convenience constructor; does not resize existing queues. */
	public static int BUFFER_PADDING=4;
	/** Output sizing multiplier applied before integer division by two; configure before construction. */
	public static int BUFFER_MULT=3;

}
