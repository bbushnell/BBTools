package fileIO;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;

/**
 * Generic thread for asynchronously loading objects from files.
 * Manages concurrent file I/O operations with thread pooling and synchronization.
 * Limits concurrent reads to prevent memory exhaustion and I/O contention.
 *
 * @author Brian Bushnell
 * @date Jan 2, 2013
 */
public class LoadThread<X> extends Thread{
	
	/** Creates, starts, and returns a LoadThread that asynchronously reads an object of type Y from a file. */
	public static <Y> LoadThread<Y> load(String fname, Class<Y> c){
		LoadThread<Y> lt=new LoadThread<Y>(fname, c);
		lt.start();
		return lt;
	}
	
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
			
			if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notifyAll();}
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
	
	/** @return The number of active load threads (waiting plus running). */
	public static final int countActiveThreads(){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
			return activeThreads[0];
		}
	}
	
	/** Blocks the calling thread until all active load threads have finished. */
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
				if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notifyAll();}
			}
		}
	}
	
	/** Blocks until this thread's load completes and {@link #output} is populated. */
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
	
	/** Shared load-thread counters indexed as [total, waiting, running]; invariant: total==waiting+running. Synchronized on to throttle concurrent reads to LIMIT. */
	//FIXED [fileIO/LoadThread#002]: single notify() -> notifyAll() (both sites). Worker threads (addRunningThread,
	//no timeout) and the main thread (waitForReadingToFinish, 8s timeout) both wait on this monitor; a finishing
	//thread's single notify() could wake the main thread (which can't progress while workers are queued) instead of a
	//ready worker, stalling that worker until the 8s timeout re-notified. notifyAll() wakes all; those that can't
	//proceed re-wait. Not a deadlock and output was always correct (the timeout self-heals), so this is a latency fix.
	//Same twin pattern fixed in ReadWrite.addRunningThread/waitForWritingToFinish. -Furina flagged 2026-06-16, G11 fixed 2026-07-04
	public static int[] activeThreads={0, 0, 0};

	private final String fname;
	private final Class<X> c;
	/** The loaded object; null until this thread finishes reading. */
	public X output=null;
	
	/** Maximum number of load threads permitted to run simultaneously. */
	public static int LIMIT=Tools.min(12, Tools.max(Shared.threads(), 1));
	
}
