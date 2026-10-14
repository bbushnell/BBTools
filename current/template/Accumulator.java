package template;

import java.util.concurrent.locks.ReadWriteLock;

/**
 * Collects caller-owned worker results after their completion.
 * ThreadWaiter invokes accumulate on its calling thread after joining workers;
 * the implementation is responsible for interpreting each worker success flag.
 * This interface itself supplies neither synchronization nor failure detection.
 *
 * @author Brian Bushnell
 * @date November 19, 2015
 * @param <T> Type of object containing statistics to accumulate
 */
public interface Accumulator<T>{
	
	/**
	 * Accumulates personal variables from finished threads.
	 * Merges thread-local statistics into the shared accumulator state.
	 * @param t Object containing thread-local statistics to merge
	 */
	public void accumulate(T t);
	
//	/** A shared lock preventing premature accumulation */
//	public ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/**
	 * Supplies the optional lock used by ThreadWaiter. Its read lock covers worker
	 * execution/waiting, and its write lock covers result accumulation. Workers
	 * must not require that write lock before they can terminate.
	 * @return Caller-owned coordination lock, or null to omit that protocol
	 */
	public ReadWriteLock rwlock();
	
	/** Reports success according to this accumulator's collected state.
	 * @return true if all operations completed without errors */
	public boolean success();
	
}
