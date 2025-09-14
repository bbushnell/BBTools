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
	
	final boolean isPoison(){return x==null;}
	
	public final X x;
	final ArrayBlockingQueue<X> dest;
	public Y result; 
	
}
