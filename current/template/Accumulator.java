package template;

import java.util.concurrent.locks.ReadWriteLock;

/**
 * Interface for accumulating statistics captured by threads.
 * Provides a contract for collecting and consolidating thread-local data
 * with thread-safe access controls.
 *
 * @author Brian Bushnell
 * @date November 19, 2015
 * @param <T> Type of object containing statistics to accumulate
 */
public interface Accumulator<T> {
	
	/**
	 * Accumulates personal variables from finished threads.
	 * Merges thread-local statistics into the shared accumulator state.
	 * @param t Object containing thread-local statistics to merge
	 */
	public void accumulate(T t);
	
//	/** A shared lock preventing premature accumulation */
//	public ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/**
	 * Gets the shared lock preventing premature accumulation.
	 * Provides thread-safe coordination between statistics collection
	 * and accumulation phases.
	 * @return ReadWriteLock for coordinating thread access
	 */
	public ReadWriteLock rwlock();
	
	/** Reports if processing finished successfully.
	 * @return true if all operations completed without errors */
	public boolean success();
	
}
