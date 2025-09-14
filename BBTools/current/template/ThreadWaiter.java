package template;

import java.util.concurrent.locks.ReadWriteLock;

public class ThreadWaiter {
	
	/** Wait for all threads to start running */
	public static final <T extends Thread> boolean waitForThreadsToStart(Iterable<T> iter){

		//Wait for all threads to start running
		boolean success=true;
		for(T t : iter){
			//Wait until this thread has started
			while(t.getState()==Thread.State.NEW){
				Thread.yield();
			}
		}
		
		return success;
	}
	
	/** Wait for completion of all threads */
	public static final <T extends Thread> boolean waitForThreadsToFinish(Iterable<T> iter){

		//Wait for completion of all threads
		boolean success=true;
		for(T t : iter){

			//Wait until this thread has terminated
			while(t.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					t.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
		}
		
		return success;
	}
	
	public static final <T extends Thread> void startThreads(Iterable<T> iter){
		for(Thread t : iter){t.start();}
	}
	
	/**
	 * @param iter List of Threads.
	 * @return success
	 */
	public static final <T extends Thread> boolean startAndWait(Iterable<T> iter){
		startThreads(iter);
		boolean sr=waitForThreadsToStart(iter);
		boolean fr=waitForThreadsToFinish(iter);
		return fr && sr;
	}
	
	/**
	 * @param iter List of Threads.
	 * @return success
	 */
	public static final <T extends Thread> boolean startAndWait(Iterable<T> iter, 
			Accumulator<T> acc){
		final ReadWriteLock rwlock=acc.rwlock();
		if(rwlock!=null) {
//			rwlock.writeLock().lock();
			rwlock.readLock().lock();
		}
		startThreads(iter);
		boolean sr=waitForThreadsToStart(iter);
		if(rwlock!=null) {
//			rwlock.writeLock().unlock();
//			rwlock.readLock().lock();
		}
		boolean fr=waitForThreadsToFinish(iter);
		if(rwlock!=null) {
			rwlock.readLock().unlock();
			rwlock.writeLock().lock();
		}
		boolean ar=accumulate(iter, acc);
		if(rwlock!=null) {
			rwlock.writeLock().unlock();
		}
		return fr && sr && ar;
	}
	
	/** Wait for completion of all threads, and accumulate results */
	public static final <T extends Thread> boolean waitForThreadsToFinish(Iterable<T> iter, 
			Accumulator<T> acc){
		final ReadWriteLock rwlock=acc.rwlock();
		if(rwlock!=null) {
//			rwlock.writeLock().lock();
			rwlock.readLock().lock();
		}
//		startThreads(iter);
		boolean sr=waitForThreadsToStart(iter);
		if(rwlock!=null) {
//			rwlock.writeLock().unlock();
//			rwlock.readLock().lock();
		}
		boolean fr=waitForThreadsToFinish(iter);
		if(rwlock!=null) {
			rwlock.readLock().unlock();
			rwlock.writeLock().lock();
		}
		boolean ar=accumulate(iter, acc);
		if(rwlock!=null) {
			rwlock.writeLock().unlock();
		}
		return fr && sr && ar;
	}
	
	/** Accumulate results from all threads */
	private static final <T> boolean accumulate(Iterable<T> iter, Accumulator<T> acc){

		//Wait for completion of all threads
		for(T t : iter){
//			assert(t.getState()==Thread.State.TERMINATED);//Not strictly necessary; requires T to be a thread.

			//Accumulate per-thread statistics
			acc.accumulate(t);
		}
		
		return acc.success();
	}
	
}
