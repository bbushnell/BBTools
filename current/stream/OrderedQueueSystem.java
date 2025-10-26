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
	
	public OrderedQueueSystem(int capacity, boolean ordered, I inputPrototype_, O outputPrototype_){
		inq=new ArrayBlockingQueue<I>(capacity);
		outq=new JobQueue<O>(capacity, ordered, true, 0);
		inputPrototype=inputPrototype_;
		outputPrototype=outputPrototype_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Producer API          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Add input job for processing. */
	public void addInput(I job){
		if(job==null){return;}
		assert(!job.poison() && !job.last()) : "Use poison() to terminate";
		synchronized(this){
			assert(!lastSeen);
			maxSeenId=Math.max(job.id(), maxSeenId);
		}
		putJob(job);
	}
	
	/** Signal end of input - injects LAST to output and POISON to input. */
	@SuppressWarnings("unchecked")
	public synchronized void poison(){
	    if(lastSeen){return;}
	    
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
		return job;
	}
	
	/** Add processed output job. */
	public void addOutput(O job){
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
		finished=true;
		this.notifyAll();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	private void putJob(I job){
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
	
}