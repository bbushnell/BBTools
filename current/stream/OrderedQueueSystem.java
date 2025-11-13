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
	public void addInput(I job){
		if(job==null){return;}
		assert(!job.last()) : "Use poison() to terminate";
		synchronized(this){
			assert(!lastSeen || job.poison());
			maxSeenId=Math.max(job.id(), maxSeenId);
		}
		putJob(job);
	}

	/** Signal end of input - injects LAST to output and POISON to input. */
	@SuppressWarnings("unchecked")
	public synchronized void poison(){
		if(lastSeen){return;}
		if(verbose) {System.err.println("OQS: poison()");}

		// Both use maxSeenId+1
		O lastJob=(O)outputPrototype.makeLast(maxSeenId+1);
		outq.add(lastJob);

		I poisonJob=(I)inputPrototype.makePoison(maxSeenId+1);
		putJob(poisonJob);

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

	/** Get next input job (blocks). */
	public I getInput(){
		I job=null;
		while(job==null){
			try{
				job=inq.take();
			}catch(InterruptedException e){
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

	/** Get next output job in order (blocks). */
	public boolean hasMore(){
		return outq.hasMore();
	}

	/** Get next output job in order (blocks). */
	public O getOutput(){
		return outq.take();
	}

	/** Signal that processing is complete. */
	public synchronized void setFinished(){
		if(verbose) {System.err.println("OQS: setFinished()");}
		finished=true;
		outq.poison();
		this.notifyAll();
	}
	
	public synchronized boolean finished() {return finished;}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private void putJob(I job){
		if(verbose) {System.err.println("OQS: putJob I "+job.id()+": "+job.poison()+", "+job.last());}
		while(job!=null){
			try{
				inq.put(job);
				job=null;
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
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