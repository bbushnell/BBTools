package template;

import java.util.concurrent.locks.ReadWriteLock;

/**
 * Interface for accumulating statistics captured by threads.
 * 
 * @author Brian Bushnell
 * @date November 19, 2015
 *
 * @param <T>
 */
public interface Accumulator<T> {
	
	/** Accumulate personal variables from finished threads */
	public void accumulate(T t);
	
//	/** A shared lock preventing premature accumulation */
//	public ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/** A shared lock preventing premature accumulation */
	public ReadWriteLock rwlock();
	
	/** True if it finished successfully */
	public boolean success();
	
}
