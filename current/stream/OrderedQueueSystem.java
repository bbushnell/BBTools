package stream;

import java.util.concurrent.ArrayBlockingQueue;

/**
 * Coordinates ordered processing with unordered worker threads.
 * 
 * Handles queue management, LAST/POISON propagation, and synchronization
 * for parallel processing pipelines that require ordered output.
 * 
 * Producer adds input jobs → Workers transform → Consumer gets ordered output
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 25, 2025
 *
 * @param <I> Input job type (must implement HasID)
 * @param <O> Output job type (must implement HasID)
 */
public class OrderedQueueSystem<I extends HasID, O extends HasID> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public OrderedQueueSystem(int threads, boolean ordered, I inputPrototype_, O outputPrototype_){
		this(threads+4, (3*threads)/2+4, ordered, inputPrototype_, outputPrototype_);
	}

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

	/** Add input job for processing. */
	public boolean addInput(I job){
		if(job==null){return false;}
		assert(!job.last()) : "Use poison() to terminate";
		//THREADING CONTRACT: SINGLE producer (javadoc: "Producer adds input jobs"). maxSeenId/lastSeen update under lock, but putJob (the inq.put) runs OUTSIDE the lock below - safe ONLY with one producer. The locks publish maxSeenId/lastSeen/finished to the WORKER and CONSUMER threads (e.g. setFinished() reads maxSeenId), NOT to order concurrent producers. With >1 producer, poison() could enqueue its input-POISON ahead of a counted-but-not-yet-put job -> a worker stops early -> that job is orphaned -> the ordered outq waits on its id forever.
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

	/** Signal end of input - injects LAST to output and POISON to input. */
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
		if(verbose) {System.err.println("OQS: poison()");}

		//End-of-input protocol: two DISTINCT markers, both stamped id=maxSeenId+1 (one past every real job, so the ORDERED outq sorts them strictly last):
		//  LAST  -> outq: marks the end of ordered output; the workers' real outputs (ids 0..maxSeenId) reorder ahead of it.
		//  POISON-> inq:  terminates the worker pool (a worker taking it stops pulling input).
		outq.add((O)outputPrototype.makeLast(id));        //blocking, but OUTSIDE the lock now
		putJob((I)inputPrototype.makePoison(id));          //blocking (inq.put), also outside the lock
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

	/** Get next input job (blocks). */
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
		if(verbose) {System.err.println("OQS: getInput I "+job.id()+": "+job.poison()+", "+job.last());}
		return job;
	}

	/** Add processed output job. */
	public void addOutput(O job){
		if(verbose) {System.err.println("OQS: addOutput O "+job.id()+": "+job.poison()+", "+job.last());}
		outq.add(job);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Consumer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** @return true if more ordered output remains to be consumed. Does NOT block or remove anything (despite the old copy-pasted javadoc - #001 doc fix). */
	public boolean hasMore(){
		return outq.hasMore();
	}

	/** Get next output job in order (blocks). */
	public O getOutput(){
		return outq.take();
	}

	/** Signal that processing is complete. */
	@SuppressWarnings("unchecked")
	public synchronized void setFinished(boolean force){
		if(verbose) {System.err.println("OQS: setFinished()");}
		finished=true;
		outq.poison((O)outputPrototype.makePoison(maxSeenId+1), force);
		this.notifyAll();
	}
	
	public synchronized boolean finished() {return finished;}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private boolean putJob(I job){
		if(verbose) {System.err.println("OQS: putJob I "+job.id()+": "+job.poison()+", "+job.last());}
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

	private final ArrayBlockingQueue<I> inq;
	private final JobQueue<O> outq;
	private final I inputPrototype;
	private final O outputPrototype;

	private long maxSeenId=-1;
	private volatile boolean finished=false;
	private volatile boolean lastSeen=false;
	private static final boolean verbose=false;

}