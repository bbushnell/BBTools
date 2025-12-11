package template;

import java.util.concurrent.ArrayBlockingQueue;

import shared.KillSwitch;

/**
 * Template class for jobs executed in thread pools with producer-consumer pattern.
 * Manages job execution, result handling, and cleanup coordination between threads.
 * Uses generics where X is the input data type and Y is the result type.
 *
 * @author Brian Bushnell
 * @date August 26, 2019
 */
public class ThreadPoolJob<X, Y> {

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
	 * Calls doWork() to process the data, then cleanup() to return job to queue. */
	final void doJob(){
		result=doWork();
		cleanup();
	}
	
	/**
	 * Abstract method for subclasses to implement job-specific processing logic.
	 * Default implementation throws KillSwitch exception to force override.
	 * @return Processed result of type Y
	 */
	public Y doWork(){
		KillSwitch.kill("Unimplemented Method");
		return null;
	}
	
	/**
	 * Returns the job input data to the destination queue after processing.
	 * Blocks until successful insertion, retrying on InterruptedException.
	 * Ensures proper coordination between producer and consumer threads.
	 */
	final void cleanup(){
		boolean success=false;
		while(!success) {
			try {
				dest.put(x);
				success=true;
			} catch (InterruptedException e) {
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
