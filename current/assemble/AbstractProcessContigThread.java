package assemble;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.KillSwitch;
import structures.ByteBuilder;

/**
 * Distributes a stable contig list across workers through a shared atomic index.
 * Each claimed contig is passed to the left hook, then the right hook. Subclasses
 * own traversal and graph updates; the atomic index alone does not synchronize
 * those updates. Scratch arrays, path state and statistics belong to this worker.
 * The coordinator starts and joins workers, checks captured failures and merges
 * their counters. This base neither closes resources nor cancels sibling workers.
 *
 * @author Brian Bushnell
 * @date July 12, 2018
 */
public abstract class AbstractProcessContigThread extends Thread{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows the stable contig list and shared next-index counter, normally zero
	 * before workers start. Installs failure capture inherited by concrete workers. */
	AbstractProcessContigThread(ArrayList<Contig> contigs_, AtomicInteger next_){
		contigs=contigs_;
		next=next_;
		setUncaughtExceptionHandler((thread, throwable) -> failure=throwable);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Processes this worker's borrowed list; uncaught failures are captured by Thread. */
	@Override
	public void run(){
		processContigs(contigs);
	}

	/**
	 * Claims successive list indexes and calls both end hooks with reusable scratch.
	 * Workers sharing next must use the same stable list and must not resize it or
	 * replace its entries during this phase. A worker instance may not be entered
	 * concurrently. Scratch is not cleared here; each hook establishes what it needs.
	 * Direct calls propagate exceptions to their caller; Thread failure capture
	 * applies only when the exception escapes a started worker's run method.
	 * @param contigs Borrowed list to distribute using this worker's shared counter
	 */
	public final void processContigs(ArrayList<Contig> contigs){
		//TODO: Probable bug ASM048 - for a near-Integer.MAX_VALUE list, extra worker claims can
		// wrap next negative and reach get() with an invalid index. A size bound or wrap-safe claim needs review.
		for(int cnum=next.getAndIncrement(); cnum<contigs.size(); cnum=next.getAndIncrement()){
			Contig c=contigs.get(cnum);
			processContigLeft(c, leftCounts, rightCounts, extraCounts, bb);
			processContigRight(c, leftCounts, rightCounts, extraCounts, bb);
		}
	}

	/**
	 * Processes the left end of a contig.
	 * Implementation defined by subclasses for specific assembly operations.
	 *
	 * @param c Contig to process
	 * @param leftCounts Array for counting left-end base occurrences
	 * @param rightCounts Array for counting right-end base occurrences
	 * @param extraCounts Array for additional counting operations
	 * @param bb ByteBuilder for sequence manipulation
	 */
	abstract void processContigLeft(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb);

	/**
	 * Processes the right end of a contig.
	 * Implementation defined by subclasses for specific assembly operations.
	 *
	 * @param c Contig to process
	 * @param leftCounts Array for counting left-end base occurrences
	 * @param rightCounts Array for counting right-end base occurrences
	 * @param extraCounts Array for additional counting operations
	 * @param bb ByteBuilder for sequence manipulation
	 */
	abstract void processContigRight(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb);

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reused four-base counts for leftward extensions. */
	final int[] leftCounts=KillSwitch.allocInt1D(4);
	/** Reused four-base counts for rightward extensions. */
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	/** Additional four-base scratch counts used by concrete traversals. */
	final int[] extraCounts=KillSwitch.allocInt1D(4);

	/** Borrowed contig list; its structure stays fixed while workers run. */
	final ArrayList<Contig> contigs;
	/** Shared next index; claiming work does not lock contigs or their graph edges. */
	final AtomicInteger next;

	/** Last traversal length in bases, including the first step from the source end. */
	int lastLength=-1;
	/** Last traversal's target contig id, or -1 when no target was accepted. */
	int lastTarget=-1;
	/** Last traversal termination code from ShaveObject. */
	int lastExitCondition=-1;
	/** Destination-side orientation bits for an accepted edge; not valid after failure. */
	int lastOrientation=-1;
	/** Reused path sequence builder; hooks reset it before constructing a path. */
	ByteBuilder bb=new ByteBuilder();
	/** Directed edges added by this worker. */
	long edgesMadeT=0;
	/** End classifications recomputed, and how many changed from their previous code. */
	long graphEndsRefreshedT=0, graphEndCodesChangedT=0;
	/** Examined endpoint seeds owned by this contig, absent, or owned by another contig. */
	long endpointSeedsUniqueT=0, endpointSeedsMissingT=0, endpointSeedsAmbiguousT=0;
	/** Graph traversals attempted when the concrete worker enables diagnostic accounting. */
	long traversalAttemptsT=0;
	/** Local branch alternatives ignored by the concrete traversal's wash policy. */
	long localHairsIgnoredT=0, localBubblesIgnoredT=0;
	/** Refreshed endpoint classifications, indexed by ShaveObject code. */
	final long[] graphEndCodeCountsT=new long[ShaveObject.MAX_CODE];
	/** Refreshed classifications for endpoints owned uniquely by their source contig. */
	final long[] uniqueEndCodeCountsT=new long[ShaveObject.MAX_CODE];
	/** Refreshed classifications for endpoint seeds absent from the table. */
	final long[] missingEndCodeCountsT=new long[ShaveObject.MAX_CODE];
	/** Refreshed classifications for endpoint seeds owned by another contig. */
	final long[] ambiguousEndCodeCountsT=new long[ShaveObject.MAX_CODE];
	/** Termination-code histogram for diagnostically counted traversals. */
	final long[] exitCountsT=new long[ShaveObject.MAX_CODE];
	/** Uncaught worker failure, captured so the coordinating thread can propagate it. */
	volatile Throwable failure=null;

}
