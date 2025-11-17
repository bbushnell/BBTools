package sketch;

import java.util.concurrent.ArrayBlockingQueue;

/**
 * Manages concurrent sequence alignment job processing with error resilience.
 * Executes sequence similarity calculations within a thread-safe, blocking queue-based
 * concurrent processing framework for the sketch comparison pipeline.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class AlignmentJob {
	
	/**
	 * Creates an alignment job with a comparison task and result destination.
	 * @param c_ The comparison task to process
	 * @param dest_ Queue to receive the processed comparison result
	 */
	AlignmentJob(Comparison c_, ArrayBlockingQueue<Comparison> dest_){
		c=c_;
		dest=dest_;
	}
	
	/**
	 * Executes the alignment job by processing the comparison and returning results.
	 * Performs SSU identity calculation on the comparison object with error handling,
	 * then places the result in the destination queue.
	 */
	void doWork(){
		assert(!isPoison());
		try {
			c.ssuIdentity();
		}catch (Throwable t){
			t.printStackTrace();
		}
		put();
	}
	
	/** Places the processed comparison into the destination queue with retry logic.
	 * Blocks until the queue accepts the result, handling interruptions gracefully. */
	private void put(){
		boolean success=false;
		while(!success){
			try {
				dest.put(c);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/** Returns true if this is a poison pill job used for thread termination */
	final boolean isPoison(){return c==null;}
	
	/** The comparison task to be processed by this alignment job */
	final Comparison c;
	/** Queue to receive the processed comparison result */
	final ArrayBlockingQueue<Comparison> dest;
	
}