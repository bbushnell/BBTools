package sketch;

import java.util.concurrent.ArrayBlockingQueue;

public class AlignmentJob {
	
	AlignmentJob(Comparison c_, ArrayBlockingQueue<Comparison> dest_){
		c=c_;
		dest=dest_;
	}
	
	void doWork(){
		assert(!isPoison());
		try {
			c.ssuIdentity();
		}catch (Throwable t){
			t.printStackTrace();
		}
		put();
	}
	
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
	
	final boolean isPoison(){return c==null;}
	
	final Comparison c;
	final ArrayBlockingQueue<Comparison> dest;
	
}