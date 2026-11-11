package assemble;

import java.util.ArrayList;

import shared.KillSwitch;
import stream.ConcurrentReadInputStream;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Shared state for Tadpole build workers. Concrete subclasses implement the
 * traversal; this base supplies worker-owned scratch, result collections and
 * counters, and captures uncaught failures for the coordinator to inspect.
 * Input streams and their array are borrowed. Callers and concrete workers
 * manage worker/stream startup, joining and resource closure; this base does
 * not perform those operations or provide a table-access synchronization policy.
 *
 * @author Brian Bushnell
 * @date Jul 18, 2015
 */
abstract class AbstractBuildThread extends Thread{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retains the worker's identity, operation mode and borrowed streams.
	 * Installs an uncaught-exception handler even when a subclass overrides run().
	 * The coordinator must join the worker and inspect failure before accepting
	 * its results; capturing a failure does not cancel other workers or close I/O.
	 *
	 * @param id_ Worker identifier used for k-mer ownership claims
	 * @param mode_ Tadpole operation mode, rather than a bitmask of independent flags
	 * @param crisa_ Borrowed input streams, or null for table-seeded assembly
	 */
	public AbstractBuildThread(final int id_, final int mode_, final ConcurrentReadInputStream[] crisa_){
		id=id_;
		crisa=crisa_;
		mode=mode_;
		setUncaughtExceptionHandler((thread, throwable) -> failure=throwable);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed stream array shared with other workers; this object does not close it. */
	final ConcurrentReadInputStream[] crisa;

	/** Selected operation interpreted by the concrete worker. */
	final int mode;
	/** Minimum seed depth for this worker's current table pass. */
	int minCountSeedCurrent;

	/** Reused worker-local counts for the four left extensions. */
	final int[] leftCounts=KillSwitch.allocInt1D(4);
	/** Reused worker-local counts for the four right extensions. */
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	/** Reused worker-local sequence construction buffer. */
	final ByteBuilder builderT=new ByteBuilder();

	/** Worker-local insert-size histogram; bin zero includes unresolved estimates. */
	final LongList insertSizes=new LongList();

	/** Contig objects collected by the worker and consumed by the coordinator after join. */
	ArrayList<Contig> contigs=new ArrayList<Contig>();

	/** Individual reads encountered by this worker, including both mates. */
	long readsInT=0;
	/** Bases in encountered reads, before assembly or filtering decisions. */
	long basesInT=0;
	/** Reads already marked discarded when the concrete read-processing path checks them. */
	long lowqReadsT=0;
	/** All bases in those discarded reads, not a count of individual low-quality bases. */
	long lowqBasesT=0;
	/** Uncaught worker failure, captured so the coordinating thread can propagate it. */
	volatile Throwable failure=null;
	/** Worker identifier supplied to the constructor. */
	final int id;

}
