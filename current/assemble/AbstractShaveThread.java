package assemble;

/**
 * Base worker for the removal phase of {@link Shaver#shave()}.
 *
 * <p>The thread owns one partition-removal loop and stores local counters that
 * the coordinating {@link Shaver} consumes after join. Uncaught failures are
 * retained for phase-level propagation instead of being lost on stderr.</p>
 *
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
public abstract class AbstractShaveThread extends Thread{

	/**
	 * Creates a removal worker and installs failure capture on this thread.
	 * @param id_ worker identifier assigned by the coordinating Shaver
	 */
	public AbstractShaveThread(int id_){
		id=id_;
		setUncaughtExceptionHandler((thread, throwable) -> failure=throwable);
	}

	/** Claims and removes table partitions until no work remains. */
	@Override
	public final void run(){
		while(processNextTable()){}
	}

	/**
	 * Claims and removes one table partition for the concrete kmer-table type.
	 * @return true if this worker should attempt another table
	 */
	protected abstract boolean processNextTable();

	/*--------------------------------------------------------------*/

	/** Per-worker removal count accumulated by Shaver after join. */
	long kmersRemovedT=0;

	/** Uncaught worker failure, captured so the coordinating Shaver can propagate it. */
	volatile Throwable failure=null;

	/** Immutable worker identifier assigned by the coordinating Shaver. */
	final int id;

}
