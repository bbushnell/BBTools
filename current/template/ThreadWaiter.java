package template;

import java.util.concurrent.locks.ReadWriteLock;

/**
 * Starts or joins caller-owned threads and optionally collects their results.
 * Joining establishes visibility of worker writes; it does not detect uncaught
 * worker failures. Implementations of Accumulator must inspect worker success flags.
 * Iterables must support repeated traversal and retain the same threads throughout
 * an operation. Callers coordinate startup; waiting on an unstarted thread can spin
 * indefinitely. Interruptions are printed and retried without restoring the flag.
 * @author Brian Bushnell
 */
public class ThreadWaiter{
	
	/**
	 * Spins until every thread leaves NEW. Another thread must start any NEW entries;
	 * this method does not start them or establish successful worker completion.
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
	 * Joins already-started threads, skipping the calling thread.
	 * Retries interrupted joins; does not inspect worker exceptions or success flags.
	 * @param iter Collection of threads to wait for completion
	 * @return true (always returns success)
	 */
	public static final <T extends Thread> boolean waitForThreadsToFinish(Iterable<T> iter){

		//Wait for completion of all threads
		boolean success=true;
		final Thread self=Thread.currentThread();
		for(T t : iter){
			if(t==self){continue;}
			//Wait until this thread has terminated
			while(t.getState()!=Thread.State.TERMINATED){
				try{
					//Attempt a join operation
					t.join();
				}catch(InterruptedException e){
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
		}
		
		return success;
	}
	
	/** Array form of the already-started-thread join; skips the calling thread.
	 * @param iter Stable array of nonnull, already-started threads
	 * @return true after all other threads terminate, regardless of worker failures */
	public static final <T extends Thread> boolean waitForThreadsToFinish(T[] iter){

		//Wait for completion of all threads
		boolean success=true;
		final Thread self=Thread.currentThread();
		for(T t : iter){
			if(t==self){continue;}
			//Wait until this thread has terminated
			while(t.getState()!=Thread.State.TERMINATED){
				try{
					//Attempt a join operation
					t.join();
				}catch(InterruptedException e){
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
		}
		
		return success;
	}
	
	/** Starts each thread in iteration order; an exception leaves earlier threads started.
	 * @param iter Collection of threads to start */
	public static final <T extends Thread> void startThreads(Iterable<T> iter){
		for(Thread t : iter){t.start();}
	}
	
	/**
	 * Starts all threads and waits for them to complete.
	 * Combines startup and joining; normal return does not certify worker success.
	 * @param iter Collection of threads to start and wait for
	 * @return true after normal startup and joining; worker failures are not inspected
	 */
	public static final <T extends Thread> boolean startAndWait(Iterable<T> iter){
		startThreads(iter);
		boolean sr=waitForThreadsToStart(iter);
		boolean fr=waitForThreadsToFinish(iter);
		return fr && sr;
	}
	
	/**
	 * Starts and joins workers, then collects their results on the calling thread.
	 * If supplied, the accumulator read lock covers startup/join and its write lock
	 * covers accumulation. Workers must not require that write lock to finish.
	 * Acquired locks are released in finally blocks before exceptions propagate.
	 * A startup exception can still leave earlier workers running; callers must not
	 * treat this as a general failure-recovery API.
	 * @param iter Threads to start and wait for
	 * @param acc Accumulator for results
	 * @return true if start, wait, and accumulate all succeed
	 */
	public static final <T extends Thread> boolean startAndWait(Iterable<T> iter, 
			Accumulator<T> acc){
		return waitAndAccumulate(iter, acc, true);
	}
	
	/**
	 * Joins already-started workers and collects results on the calling thread.
	 * Uses the same optional read-then-write lock protocol as startAndWait.
	 * Workers must not require the write lock to finish. Acquired locks are released
	 * in finally blocks before exceptions propagate. This method never starts NEW entries.
	 *
	 * @param iter Collection of threads to wait for completion
	 * @param acc Accumulator for collecting thread results
	 * @return true if waiting, finishing, and accumulation all succeed, false otherwise
	 */
	public static final <T extends Thread> boolean waitForThreadsToFinish(Iterable<T> iter, 
			Accumulator<T> acc){
		return waitAndAccumulate(iter, acc, false);
	}

	/** Runs the optional startup, join and accumulation phases under their locks.
	 * Exceptions propagate unchanged after releasing this invocation's lock holds.
	 * @param iter Stable collection of workers
	 * @param acc Accumulator whose optional lock protects each phase
	 * @param start Whether to start the workers before joining
	 * @return Combined wait and accumulator success */
	private static <T extends Thread> boolean waitAndAccumulate(final Iterable<T> iter,
			final Accumulator<T> acc, final boolean start){
		final ReadWriteLock rwlock=acc.rwlock();
		final boolean sr, fr;
		if(rwlock!=null){rwlock.readLock().lock();}
		try{
			if(start){startThreads(iter);}
			sr=waitForThreadsToStart(iter);
			fr=waitForThreadsToFinish(iter);
		}finally{
			if(rwlock!=null){rwlock.readLock().unlock();}
		}
		if(rwlock!=null){rwlock.writeLock().lock();}
		try{
			final boolean ar=accumulate(iter, acc);
			return fr && sr && ar;
		}finally{
			if(rwlock!=null){rwlock.writeLock().unlock();}
		}
	}

	/**
	 * Accumulates results from all threads using the provided accumulator.
	 * Delegates in iteration order after the caller has joined its workers.
	 * Success is determined by the accumulator, not by thread termination alone.
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
