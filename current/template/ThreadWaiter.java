package template;

import java.util.concurrent.locks.ReadWriteLock;

/**
 * Utility helpers for starting, waiting on, and aggregating thread results (with optional locks).
 * Provides thread lifecycle helpers and accumulation via an Accumulator.
 * @author Brian Bushnell
 */
public class ThreadWaiter {
	
	/**
	 * Spins until all threads in the iterable leave the NEW state.
	 * @param iter Threads to monitor
	 * @return true (success)
	 */
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
	
	/**
	 * Waits for completion of all threads by monitoring their termination state.
	 * Uses join() operations with exception handling for interrupted threads.
	 * @param iter Collection of threads to wait for completion
	 * @return true (always returns success)
	 */
	public static final <T extends Thread> boolean waitForThreadsToFinish(Iterable<T> iter){

		//Wait for completion of all threads
		boolean success=true;
		final Thread self=Thread.currentThread();
		for(T t : iter){
			if(t==self) {continue;}
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
	
	public static final <T extends Thread> boolean waitForThreadsToFinish(T[] iter){

		//Wait for completion of all threads
		boolean success=true;
		final Thread self=Thread.currentThread();
		for(T t : iter){
			if(t==self) {continue;}
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
	
	/** Starts all threads in the provided iterable.
	 * @param iter Collection of threads to start */
	public static final <T extends Thread> void startThreads(Iterable<T> iter){
		for(Thread t : iter){t.start();}
	}
	
	/**
	 * Starts all threads and waits for them to complete.
	 * Combines thread startup and completion synchronization in a single operation.
	 * @param iter Collection of threads to start and wait for
	 * @return true if both starting and finishing operations succeed, false otherwise
	 */
	public static final <T extends Thread> boolean startAndWait(Iterable<T> iter){
		startThreads(iter);
		boolean sr=waitForThreadsToStart(iter);
		boolean fr=waitForThreadsToFinish(iter);
		return fr && sr;
	}
	
	/**
	 * Starts threads, waits for completion, and accumulates results using the provided Accumulator (with optional locks).
	 * @param iter Threads to start and wait for
	 * @param acc Accumulator for results
	 * @return true if start, wait, and accumulate all succeed
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
	
	/**
	 * Waits for thread completion and accumulates results with thread synchronization.
	 * Uses read-write locks to coordinate between thread execution and result accumulation.
	 * Similar to startAndWait but assumes threads are already started.
	 *
	 * @param iter Collection of threads to wait for completion
	 * @param acc Accumulator for collecting thread results
	 * @return true if waiting, finishing, and accumulation all succeed, false otherwise
	 */
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
	
	/**
	 * Accumulates results from all threads using the provided accumulator.
	 * Iterates through all threads and delegates result accumulation to the accumulator.
	 *
	 * @param iter Collection of threads to accumulate results from
	 * @param acc Accumulator interface for collecting thread results
	 * @return success status from the accumulator
	 */
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
