package template;

import java.util.concurrent.ArrayBlockingQueue;

import shared.KillSwitch;

/**
 * 
 * @author Brian Bushnell
 * @date August 26, 2019
 *
 */
public class ThreadPoolJob<X, Y> {

	/**
	 * Constructs a ThreadPoolJob with input data and destination queue.
	 * @param x_ Input data to be processed by this job
	 * @param dest_ Queue where completed jobs are returned for coordination
	 */
	public ThreadPoolJob(X x_, ArrayBlockingQueue<X> dest_){
		x=x_;
		dest=dest_;
	}
	
	/** Process a job */
	final void doJob(){
		result=doWork();
		cleanup();
	}
	
	/** Do whatever specific work needs to be done for this job */
	public Y doWork(){
		KillSwitch.kill("Unimplemented Method");
		return null;
	}
	
	/** Retire the job to the destination queue */
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
	
	/** Checks if this job is a poison pill for thread pool shutdown.
	 * @return true if x is null (poison pill), false otherwise */
	final boolean isPoison(){return x==null;}
	
	/** Input data to be processed by this job */
	public final X x;
	/** Destination queue for returning completed jobs */
	final ArrayBlockingQueue<X> dest;
	/** Result of job processing, set by doWork() method */
	public Y result; 
	
}
