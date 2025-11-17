package fileIO;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Jan 2, 2013
 *
 */
public class LoadThread<X> extends Thread{
	
	/**
	 * Factory method to create and start a new LoadThread for the specified file.
	 * Automatically starts the thread and returns the handle for monitoring completion.
	 *
	 * @param fname Path to the file to load
	 * @param c Class type of the object to deserialize
	 * @return Started LoadThread instance
	 */
	public static <Y> LoadThread<Y> load(String fname, Class<Y> c){
		LoadThread<Y> lt=new LoadThread<Y>(fname, c);
		lt.start();
		return lt;
	}
	
	/**
	 * Constructs a new LoadThread for the specified file and class type.
	 * Registers the thread in the thread pool management system.
	 * @param fname_ Path to the file to load
	 * @param c_ Class type of the object to deserialize
	 */
	private LoadThread(String fname_, Class<X> c_){
		fname=fname_;
		c=c_;
		addThread(1);
	}
	
	@Override
	public void run(){
		addRunningThread(1);
		output=ReadWrite.read(c, fname, false);
		addRunningThread(-1);
		synchronized(this){this.notify();}
	}
	
	
	/**
	 * Adds threads to the pool management system with memory-aware limits.
	 * Maintains thread count invariants and handles pool capacity constraints.
	 * @param x Number of threads to add (positive) or remove (negative)
	 * @return Total number of active threads after the operation
	 */
	private static final int addThread(int x){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				activeThreads[0]+=x;
				activeThreads[1]+=x;
			}else{
				addRunningThread(x);
			}
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
					
			return activeThreads[0];
		}
	}
	
	/**
	 * Manages transitions between waiting and running thread states.
	 * Blocks when thread limit is reached to prevent resource exhaustion.
	 * Maintains thread pool state consistency and notifies waiting threads.
	 *
	 * @param x Number of threads changing state (positive to start, negative to stop)
	 * @return Current number of running threads
	 */
	private static final int addRunningThread(int x){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				assert(activeThreads[1]>=x);
				while(activeThreads[2]>=lim){
					try {
						activeThreads.wait();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				activeThreads[1]-=x; //Remove from waiting
			}else{
				activeThreads[0]+=x; //Remove from active
			}
			activeThreads[2]+=x; //Change number running
			
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
			
			if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notify();}
//			System.err.println(activeThreads[2]);
//			try {
//				activeThreads.wait(5000);
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
			return activeThreads[2];
		}
	}
	
	/**
	 * Returns the total number of active threads in the pool.
	 * Active threads include both waiting and currently running threads.
	 * @return Current count of active LoadThread instances
	 */
	public static final int countActiveThreads(){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
			return activeThreads[0];
		}
	}
	
	/**
	 * Blocks until all LoadThread instances have completed their work.
	 * Polls thread state with timeout to handle potential deadlocks.
	 * Notifies waiting threads when conditions change.
	 */
	public static final void waitForReadingToFinish(){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			while(activeThreads[0]>0){
				assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
						activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
				try {
					activeThreads.wait(8000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notify();}
			}
		}
	}
	
	/**
	 * Blocks until this specific LoadThread completes execution.
	 * Uses thread state monitoring and join operations for synchronization.
	 * Returns immediately if the output is already available.
	 */
	public final void waitForThisToFinish(){
		if(output==null){
			while(this.getState()!=State.TERMINATED){
				try {
					this.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/** {active, waiting, running} <br>
	 * Active means running or waiting.
	 */
	public static int[] activeThreads={0, 0, 0};
	
	/** Path to the file being loaded by this thread */
	private final String fname;
	/** Class type of the object to deserialize from the file */
	private final Class<X> c;
	/** The loaded object result, null until loading completes */
	public X output=null;
	
	/** Unused legacy field for tracking running threads */
	private static final int[] RUNNING=new int[1];
	/** Maximum number of concurrent LoadThread instances allowed */
	public static int LIMIT=Tools.min(12, Tools.max(Shared.threads(), 1));
	
}
