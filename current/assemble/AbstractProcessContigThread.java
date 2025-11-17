package assemble;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.KillSwitch;
import structures.ByteBuilder;

/**
 * Thread for exploring connectivity graph between contigs.
 * @author Brian Bushnell
 * @date July 12, 2018
 *
 */
public abstract class AbstractProcessContigThread extends Thread {

	/**
	 * Constructs a thread for processing contigs.
	 * @param contigs_ List of contigs to process
	 * @param next_ Atomic counter for thread-safe work distribution
	 */
	AbstractProcessContigThread(ArrayList<Contig> contigs_, AtomicInteger next_){
		contigs=contigs_;
		next=next_;
	}
	
	@Override
	public void run(){
		processContigs(contigs);
	}

	/**
	 * Processes all contigs assigned to this thread.
	 * Uses atomic indexing to claim contigs and processes both left and right
	 * ends of each contig. Thread-safe for concurrent execution.
	 * @param contigs List of contigs to process
	 */
	public final void processContigs(ArrayList<Contig> contigs){
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

	/** Array for counting base occurrences at left contig ends */
	final int[] leftCounts=KillSwitch.allocInt1D(4);
	/** Array for counting base occurrences at right contig ends */
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	/** Array for additional counting operations during contig processing */
	final int[] extraCounts=KillSwitch.allocInt1D(4);

	/** List of contigs assigned to this thread for processing */
	final ArrayList<Contig> contigs;
	/** Atomic counter for thread-safe distribution of work across threads */
	final AtomicInteger next;

	/** Length from the last processed contig or operation */
	int lastLength=-1;
	/** Target identifier from the last processed contig or operation */
	int lastTarget=-1;
	/** Exit condition code from the last processed contig or operation */
	int lastExitCondition=-1;
	/** Orientation value from the last processed contig or operation */
	int lastOrientation=-1;
	/** ByteBuilder for efficient sequence manipulation during contig processing */
	ByteBuilder bb=new ByteBuilder();
	/** Count of edges created by this thread during contig processing */
	long edgesMadeT=0;

}
