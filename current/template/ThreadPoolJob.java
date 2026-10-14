package template;

import java.util.concurrent.ArrayBlockingQueue;

import shared.KillSwitch;

/**
 * Base job for a package-local executor that returns inputs to a bounded queue.
 * A subclass supplies doWork; the default implementation terminates the JVM.
 * doJob stores its result and returns the input only after doWork returns
 * normally. It does not return the input after an exception. Callers handle
 * poison markers before execution: ArrayBlockingQueue cannot accept null inputs.
 *
 * @author Brian Bushnell
 * @date August 26, 2019
 */
public class ThreadPoolJob<X, Y>{

	/**
	 * Constructs a ThreadPoolJob with input data and a destination queue for returning completed jobs.
	 * @param x_ Input data to be processed by this job
	 * @param dest_ Queue where completed jobs are returned for coordination
	 */
	public ThreadPoolJob(X x_, ArrayBlockingQueue<X> dest_){
		x=x_;
		dest=dest_;
	}
	
	/** Executes the complete job workflow.
	 * Calls doWork(), stores its result, then returns the input to the queue.
	 * If doWork throws, cleanup is skipped; this is not a finally block. */
	final void doJob(){
		result=doWork();
		cleanup();
	}
	
	/**
	 * Customization method for job-specific processing.
	 * The concrete default calls KillSwitch.kill and does not return normally.
	 * @return Processed result of type Y
	 */
	public Y doWork(){
		KillSwitch.kill("Unimplemented Method");
		return null;
	}
	
	/**
	 * Returns the job input data to the destination queue after processing.
	 * Blocks until successful insertion, retrying on InterruptedException.
	 * Does not restore the interrupt flag or recognize poison markers here.
	 */
	final void cleanup(){
		boolean success=false;
		while(!success){
			try{
				dest.put(x);
				success=true;
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/** Checks whether this job is a poison pill for thread pool shutdown.
	 * @return true if the input payload is null, false otherwise */
	final boolean isPoison(){return x==null;}
	
	/** Input data to be processed by this job. */
	public final X x;
	/** Destination queue for returning completed jobs. */
	final ArrayBlockingQueue<X> dest;
	/** Result of job processing, set by doWork(). */
	public Y result; 
	
}
