package assemble;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.KillSwitch;
import structures.ByteBuilder;

/**
 * Abstract base class for multithreaded contig processing during assembly.
 * Each thread processes contigs from a shared list using atomic indexing
 * for thread-safe work distribution. Subclasses implement specific
 * processing logic for left and right contig ends.
 *
 * @author Brian Bushnell
 * @date July 12, 2018
 */
public abstract class AbstractProcessContigThread extends Thread {

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

	final int[] leftCounts=KillSwitch.allocInt1D(4);
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	final int[] extraCounts=KillSwitch.allocInt1D(4);

	final ArrayList<Contig> contigs;
	final AtomicInteger next;

	int lastLength=-1;
	int lastTarget=-1;
	int lastExitCondition=-1;
	int lastOrientation=-1;
	ByteBuilder bb=new ByteBuilder();
	long edgesMadeT=0;

}
