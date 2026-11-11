package assemble;

import shared.KillSwitch;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Base worker for the exploration phase of {@link Shaver#shave()}.
 *
 * <p>Each instance owns one backing {@link Thread}, two thread-local kmer
 * workspaces, and per-worker statistics that the coordinating {@link Shaver}
 * reads only after the worker joins. Uncaught worker failures are captured for
 * phase-level propagation instead of being lost on stderr.</p>
 *
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
public abstract class AbstractExploreThread extends ShaveObject implements Runnable{

	/**
	 * Constructs thread-local kmer workspaces and the backing worker thread.
	 * @param id_ Thread identifier for tracking and coordination
	 * @param kbig_ Kmer length for graph structures being explored
	 */
	public AbstractExploreThread(int id_, int kbig_){
		id=id_;
		myKmer=new Kmer(kbig_);
		myKmer2=new Kmer(kbig_);
		thread=new Thread(this);
		thread.setUncaughtExceptionHandler((thread, throwable) -> failure=throwable);
	}

	/** Explores primary table cells and overflow-forest victims, then tallies bubbles. */
	@Override
	public final void run(){
		//TODO:

		//With processNextVictims enabled, the number of dead ends found drops from the first pass to the next, then stabilizes.
		//So, they are not being reset correctly.

		//Also, the number found - even with one thread - is nondeterministic if both are enabled.
		//Unstable whether or not processNextVictims is disabled.  But that's probably to be expected as the count is not exact.
		//What should be exact is the number of kmers removed for being dead ends.

		//The number is lower than expected.  65k for 600k reads with errors.  Most are bubbles, but 40% should be dead ends, or 240k.

		while(processNextTable(myKmer, myKmer2)){}
		while(processNextVictims(myKmer, myKmer2)){}

		for(int i=0; i<removeMatrixT.length; i++){
			for(int j=0; j<removeMatrixT.length; j++){
				if((i==F_BRANCH || i==B_BRANCH) && (j==F_BRANCH || j==B_BRANCH)){
					bubblesFoundT+=removeMatrixT[i][j];
				}
			}
		}
	}

	/** Uses this worker's private kmer pair to claim and process the next table. */
	boolean processNextTable(){return processNextTable(myKmer, myKmer2);}
	/**
	 * Claims and explores one table partition for the concrete kmer-table type.
	 * @param kmer primary thread-local kmer for graph traversal
	 * @param temp secondary scratch kmer kept distinct from live traversal state
	 * @return true if this worker should attempt another table
	 */
	protected abstract boolean processNextTable(final Kmer kmer, Kmer temp);

	/** Uses this worker's private kmer pair to traverse overflow-forest victims. */
	boolean processNextVictims(){return processNextVictims(myKmer, myKmer2);}//was (myKmer, myKmer): aliased scratch temp onto the live kmer; sibling processNextTable() above uses (myKmer, myKmer2). Dead convenience method (no callers) - fixed as landmine removal [assemble/AbstractExploreThread#001 FIXED]
	/**
	 * Traverses each table's overflow-forest victims without aliasing scratch
	 * state onto the primary traversal kmer.
	 */
	protected abstract boolean processNextVictims(final Kmer kmer, Kmer temp);

	/*--------------------------------------------------------------*/

	/** Starts the owned worker thread. */
	public final void start(){thread.start();}
	/** Returns the owned worker thread state for Shaver's join loop. */
	public final Thread.State getState(){return thread.getState();}
	/** Waits for the owned worker thread to finish. */
	public final void join() throws InterruptedException{thread.join();}

	/*--------------------------------------------------------------*/

	/** Number of primary or victim kmers tested by this worker. */
	long kmersTestedT=0;
	/** Aggregate count reported as dead ends; can include bubble paths selected by exploreAndMark(). */
	long deadEndsFoundT=0;
	/** Number of bubble paths inferred from this worker's remove matrix. */
	long bubblesFoundT=0;

	/** Immutable worker identifier assigned by the coordinating Shaver. */
	final int id;
	/** Primary and secondary thread-local kmers kept distinct during traversal. */
	final Kmer myKmer, myKmer2;

	/** Reused left-extension counts owned by this worker. */
	final int[] leftCounts=KillSwitch.allocInt1D(4);
	/** Reused right-extension counts owned by this worker. */
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	/** Reused path buffer owned by this worker. */
	final ByteBuilder builderT=new ByteBuilder();

	/** Per-worker transition matrix for all explored path termination pairs. */
	long[][] countMatrixT=new long[MAX_CODE+1][MAX_CODE+1];
	/** Per-worker transition matrix for paths selected for removal before ownership checks. */
	long[][] removeMatrixT=new long[MAX_CODE+1][MAX_CODE+1];

	/** Uncaught worker failure, captured so the coordinating Shaver can propagate it. */
	volatile Throwable failure=null;

	/** Backing thread whose uncaught exceptions populate {@link #failure}. */
	public final Thread thread;

}
