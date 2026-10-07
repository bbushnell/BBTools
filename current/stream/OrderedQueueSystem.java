package stream;

import java.util.concurrent.ArrayBlockingQueue;

/**
 * Supplies FIFO input and ordered output queues for caller-owned workers.
 * Creates no threads and does not process jobs or join workers. Ordinary input
 * submission and final poison publication follow the single-producer contract;
 * the FIFO does not reorder or validate input IDs. Output IDs must be unique,
 * dense from zero and stable while queued, with a corresponding result per input ID.
 * The requested output policy is forwarded, but current JobQueue forces ordering
 * even when ordered is false; that flag does not permit missing output IDs.
 *
 * Workers receive the actual input poison object and own its handling/propagation.
 * Output LAST is also a returned marker payload; output poison maps to null.
 * An external consumer/owner calls setFinished() to notify completion. That call
 * requests output-queue termination only and does not stop or join input workers.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 25, 2025
 *
 * @param <I> Input job type (must implement HasID)
 * @param <O> Output job type (must implement HasID)
 */
public class OrderedQueueSystem<I extends HasID, O extends HasID>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Derives input capacity as threads+4 and output backpressure as (3*threads)/2+4.
	 * Uses integer arithmetic and delegates argument validation to the queue constructors.
	 * @param threads Capacity-sizing hint only; no threads are created
	 * @param ordered Requested output policy; current JobQueue still forces ordering
	 * @param inputPrototype_ Retained factory for type-compatible input poison markers
	 * @param outputPrototype_ Retained factory for output LAST and poison markers */
	public OrderedQueueSystem(int threads, boolean ordered, I inputPrototype_, O outputPrototype_){
		this(threads+4, (3*threads)/2+4, ordered, inputPrototype_, outputPrototype_);
	}

	/** Creates an input FIFO and an output JobQueue starting at ID zero.
	 * Prototypes are retained without null checks and must supply nonnull type-compatible
	 * markers when signaling; output markers must fit the queue's ordering contract.
	 * @param capacityIn Positive FIFO capacity required by ArrayBlockingQueue
	 * @param capacityOut Output backpressure setting; JobQueue asserts greater than one and applies a minimum of two
	 * @param ordered Requested output policy, currently forced true inside JobQueue
	 * @param inputPrototype_ Input marker factory, compatible with the caller's worker protocol
	 * @param outputPrototype_ Output marker factory, compatible with ordered terminal IDs */
	public OrderedQueueSystem(int capacityIn, int capacityOut,
			boolean ordered, I inputPrototype_, O outputPrototype_){
		inq=new ArrayBlockingQueue<I>(capacityIn);
		outq=new JobQueue<O>(capacityOut, ordered, true, 0);
		inputPrototype=inputPrototype_;
		outputPrototype=outputPrototype_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Producer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** Registers an input ID, then performs a blocking FIFO enqueue outside this monitor.
	 * Use a single ordinary producer and coordinate poison() after its submissions;
	 * the metadata lock does not make registration and FIFO insertion one atomic operation.
	 * This method does not copy jobs, reorder IDs or validate global uniqueness.
	 * @param job Borrowed input; null is ignored, LAST is rejected when assertions are enabled
	 * @return false for null, otherwise true after enqueue completes */
	public boolean addInput(I job){
		if(job==null){return false;}
		assert(!job.last()) : "Use poison() to terminate";
		//Historical threading contract: a single ordinary producer coordinates submission and poison().
		//Metadata registration is locked, but FIFO insertion follows outside the lock. Concurrent ordinary
		//producers could let poison overtake a registered job; metadata publication does not order their puts.
		synchronized(this){
			assert(!lastSeen || job.poison());
			maxSeenId=Math.max(job.id(), maxSeenId);
//			if(finished) {
//				if(job.poison()) {return inq.offer(job);}
//				return false;
//			}
		}
		return putJob(job);
	}

	/** Seals input registration once and requests output LAST followed by input POISON.
	 * Both factories receive maxSeenId+1. Ordinary producer submissions must be complete
	 * before this call. Blocking marker enqueues occur outside the metadata monitor;
	 * later callers can return while the first caller is still publishing its markers. */
	@SuppressWarnings("unchecked")
	public void poison(){
		//[stream/bam/BgzfInputStreamMT3#003] FIXED 2026-06-20 (greenlit): poison() was `synchronized(this)`
		//and did the BLOCKING enqueues (outq.add @ a full bounded JobQueue, and putJob->inq.put @ a full
		//ArrayBlockingQueue) WHILE HOLDING the monitor. If outq/inq was full, poison() blocked holding `this`,
		//so setFinished() (also synchronized(this)) could never run to set outq.poisoned=true + notifyAll and
		//wake the blocked enqueue -> circular deadlock (jstack-proven via MT3 on truncated input: producer in
		//poison()->outq.add holding the OQS monitor, main in close()->setFinished waiting for it). The 6 other
		//OQS users (BamStreamer, Fasta/Fastq/SamStreamer, ByteFile3/4) never hit it because a SEPARATE consumer
		//thread keeps draining outq until the LAST marker, so outq.add never blocks - but the hazard is real.
		//FIX = mirror addInput's already-blessed contract: update lastSeen/read maxSeenId UNDER the lock, do the
		//blocking enqueues OUTSIDE it. Now a full queue can't wedge the monitor that setFinished needs to unstick it.
		final long id;
		synchronized(this){
			if(lastSeen){return;}//idempotent: a no-op once already poisoned
			lastSeen=true;
			id=maxSeenId+1;
		}
		if(verbose){System.err.println("OQS: poison()");}

		//End-of-input protocol: two marker roles, both requested with id=maxSeenId+1.
		//LAST -> outq: its ordered ID follows real outputs 0..maxSeenId.
		//POISON -> inq: caller-owned workers handle it; a factory may use a fixed poison ID.
		outq.add((O)outputPrototype.makeLast(id));//Blocking, but outside the lock.
		putJob((I)inputPrototype.makePoison(id));//Blocking FIFO enqueue, also outside the lock.
	}

	/** Waits for the flag set by an external setFinished() call, not for worker joins.
	 * Interrupted waits print a trace and continue without restoring the interrupt. */
	public synchronized void waitForFinish(){
		while(!finished){
			try{this.wait();}catch(InterruptedException e){e.printStackTrace();}
		}
	}

	/** Publishes terminal markers, then waits for an external completion notification. */
	public void poisonAndWait(){
		poison();
		waitForFinish();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker API           ----------------*/
	/*--------------------------------------------------------------*/

	/** Takes the next FIFO entry, retrying interrupted waits after printing their traces.
	 * Input poison is returned unchanged; workers must recognize and propagate it
	 * according to their protocol. This wrapper does not recycle that marker itself.
	 * @return Nonnull borrowed job, including an input poison marker */
	//TODO: Possible bug - getInput/putJob swallow interrupts entirely (print + retry): workers are
	//unkillable except via the poison pill.  Consistent with the house poison-driven-shutdown
	//policy (see JobQueue.take()'s documented hazard), but a stray interrupt on a healthy worker
	//prints an alarming stack trace and is otherwise ignored.  If interrupt-based shutdown is ever
	//needed, make it poison-driven instead of adding a break here.
	public I getInput(){
		I job=null;
		while(job==null){
			try{
				job=inq.take();
			}catch(InterruptedException e){
//				synchronized(this) {
//					if(finished) {
//						//I'm not sure what to do... return null?
//					}
//				}
				e.printStackTrace();
			}
		}
		if(verbose){System.err.println("OQS: getInput I "+job.id()+": "+job.poison()+", "+job.last());}
		return job;
	}

	/** Enqueues a borrowed output job under JobQueue's blocking/order policy.
	 * Ordinary result IDs must densely cover the input range; poison() supplies LAST.
	 * @param job Nonnull output with a stable corresponding input ID */
	public void addOutput(O job){
		if(verbose){System.err.println("OQS: addOutput O "+job.id()+": "+job.poison()+", "+job.last());}
		outq.add(job);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Consumer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** Delegates output availability without waiting for a job or removing one (historical #001 doc fix).
	 * @return Current JobQueue.hasMore() result; not a guarantee of a nonnull next output */
	public boolean hasMore(){return outq.hasMore();}

	/** Takes the next output under the current ordered-queue contract.
	 * @return Borrowed job, including LAST itself, or null for a queue terminal condition.
	 * Consumers must recognize LAST or accept its marker payload. */
	public O getOutput(){return outq.take();}

	/** Sets the external-completion flag, requests output termination and notifies waiters.
	 * Does not enqueue input poison, stop workers or join them. Caller coordination
	 * determines when processing is complete; this method does not verify results.
	 * @param force Forwarded unchanged to the output JobQueue.poison call */
	//TODO: Possible bug - see JobQueue.poison(): force=false does not guarantee the consumer wakes
	//if the output stream has an ID gap (dead worker); error paths should pass force=true.  Also
	//see JobQueue.hasMore(): after force=true, hasMore() stays true while take() returns null.
	@SuppressWarnings("unchecked")
	public synchronized void setFinished(boolean force){
		if(verbose){System.err.println("OQS: setFinished()");}
		finished=true;
		outq.poison((O)outputPrototype.makePoison(maxSeenId+1), force);
		this.notifyAll();
	}

	/** @return Stored completion notification, not an independently measured worker status */
	public synchronized boolean finished(){return finished;}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Enqueues the supplied job, retrying after interrupted puts without restoring the interrupt.
	 * Does not consult finished or lastSeen; producer coordination belongs to its callers.
	 * @param job Nonnull job supplied by addInput or the input marker factory
	 * @return true after the enqueue completes */
	private boolean putJob(I job){
		if(verbose){System.err.println("OQS: putJob I "+job.id()+": "+job.poison()+", "+job.last());}
		while(job!=null){
			try{
				inq.put(job);
				job=null;
			}catch(InterruptedException e){
//				synchronized(this) {
//					if(finished) {
//						if(job.poison()) {
//							//return inq.offer(job);
//							continue;//Risk of it never getting inserted
//						}else {
//							//return false;
//						}
//					}
//				}
				e.printStackTrace();
			}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** FIFO input storage; no ID reordering or automatic marker handling. */
	private final ArrayBlockingQueue<I> inq;
	/** Output queue, currently ordered regardless of the requested policy. */
	private final JobQueue<O> outq;
	/** Retained input poison factory; callers supply compatible worker marker handling. */
	private final I inputPrototype;
	/** Retained output LAST/poison factory. */
	private final O outputPrototype;

	/** Largest registered input ID; metadata reads and writes use this monitor. */
	private long maxSeenId=-1;
	/** Completion flag written by the external consumer/owner. */
	private volatile boolean finished=false;
	/** Whether poison() sealed input registration, not whether its enqueues finished. */
	private volatile boolean lastSeen=false;
	/** Compile-time diagnostic switch, disabled in normal builds. */
	private static final boolean verbose=false;

}
