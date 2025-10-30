package stream;

import java.util.concurrent.ArrayBlockingQueue;

/**
 * OrderedQueueSystem2: interrupt-aware variant of OrderedQueueSystem.
 *
 * Differences vs. OrderedQueueSystem:
 * - Uses blocking put/take but honors interrupts to enable fast shutdown.
 * - Does not spin/offer-poll; callers can interrupt worker/producer threads
 *   to break out of waits immediately.
 * - Preserves LAST-to-output and POISON-to-input semantics for ordered pipelines.
 *
 * Producer adds input jobs → Workers transform → Consumer gets ordered output
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 30, 2025
 */
public class OrderedQueueSystem2<I extends HasID, O extends HasID> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public OrderedQueueSystem2(int capacity, boolean ordered,
		I inputPrototype_, O outputPrototype_){
		inq=new ArrayBlockingQueue<I>(Math.max(2, capacity));
		outq=new JobQueue<O>(Math.max(2, capacity), ordered, true, 0);
		inputPrototype=inputPrototype_;
		outputPrototype=outputPrototype_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Producer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** Add input job for processing. Honors interrupts; does not flip cancellation. */
	public void addInput(I job){
		if(job==null){return;}
		assert(!job.last()) : "Use poison() to terminate";
		synchronized(this){
			if(lastSeen || cancelled){return;}
			maxSeenId=Math.max(job.id(), maxSeenId);
		}
		try{
			inq.put(job);
		}catch(InterruptedException ie){
			// Preserve interrupt; caller decides whether to cancel/poison
			Thread.currentThread().interrupt();
		}
	}

	/** Signal end of input - injects LAST to output and POISON to input. */
	@SuppressWarnings("unchecked")
	public synchronized void poison(){
		if(lastSeen){return;}
		if(cancelled){lastSeen=true; notifyAll(); return;}

		// Use maxSeenId+1 as the terminal ID for both streams
		final long termId=maxSeenId+1;
		O lastJob=(O)outputPrototype.makeLast(termId);
		outq.add(lastJob);

		I poisonJob=(I)inputPrototype.makePoison(termId);
		try{
			inq.put(poisonJob);
		}catch(InterruptedException ie){
			// Preserve interrupt; consumers will still see LAST on outq
			Thread.currentThread().interrupt();
		}
		lastSeen=true;
		notifyAll();
	}

	/** Cancel processing (fast shutdown). Interrupt worker/producer threads externally. */
	public synchronized void cancel(){
		cancelled=true;
		notifyAll();
		// Optional: nudge any blocked takers by injecting a poison
		try{
			@SuppressWarnings("unchecked")
			I poison=(I)inputPrototype.makePoison(maxSeenId+1);
			inq.offer(poison);
		}catch(Throwable t){/* ignore */}
	}

	/** Wait for processing to complete. Interrupts will break early. */
	public synchronized void waitForFinish(){
		while(!finished && !cancelled){
			try{this.wait();}
			catch(InterruptedException ie){
				// If shutting down, exit; otherwise, continue waiting
				if(cancelled || lastSeen){break;}
				// Do not preserve interrupt here to avoid tight interrupt loops in callers
			}
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

	/** Get next input job (blocks). On interrupt during shutdown, returns a poison job. */
	@SuppressWarnings("unchecked")
	public I getInput(){
		while(true){
			try{
				I job=inq.take();
				if(verbose){System.err.println("OQS2: getInput I "+job.id()+": "+job.poison()+", "+job.last());}
				return job;
			}catch(InterruptedException ie){
				// If cancellation or last seen, hand back a poison to let worker exit
				if(cancelled || lastSeen){
					final long termId=Math.max(maxSeenId+1, 0);
					return (I)inputPrototype.makePoison(termId);
				}
				// Otherwise, ignore transient interrupt and continue waiting
				// (do not re-interrupt here to avoid immediate rethrow in take())
			}
		}
	}

	/** Add processed output job. */
	public void addOutput(O job){
		if(verbose){System.err.println("OQS2: addOutput O "+job.id()+": "+job.poison()+", "+job.last());}
		outq.add(job);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Consumer API          ----------------*/
	/*--------------------------------------------------------------*/

	/** Get next output job in order (blocks). */
	public O getOutput(){
		return outq.take();
	}

	/** Signal that processing is complete. */
	public synchronized void setFinished(){
		if(verbose){System.err.println("OQS2: setFinished()");}
		finished=true;
		this.notifyAll();
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
	private volatile boolean cancelled=false;
	private static final boolean verbose=false;
}
