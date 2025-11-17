package assemble;

import shared.KillSwitch;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Searches for dead ends.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 *
 */
abstract class AbstractExploreThread extends ShaveObject implements Runnable {

	/**
	 * Constructor
	 */
	public AbstractExploreThread(int id_, int kbig_){
		id=id_;
		myKmer=new Kmer(kbig_);
		myKmer2=new Kmer(kbig_);
		thread=new Thread(this);
	}

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

	/**
	 * Processes the next table using the thread's default k-mer objects.
	 * Convenience method that delegates to the abstract processNextTable implementation.
	 * @return true if processing should continue, false if complete
	 */
	boolean processNextTable(){return processNextTable(myKmer, myKmer2);}
	abstract boolean processNextTable(final Kmer kmer, Kmer temp);

	/**
	 * Processes the next victim list using the thread's default k-mer objects.
	 * Convenience method that delegates to the abstract processNextVictims implementation.
	 * @return true if processing should continue, false if complete
	 */
	boolean processNextVictims(){return processNextVictims(myKmer, myKmer);}
	abstract boolean processNextVictims(final Kmer kmer, Kmer temp);

	/*--------------------------------------------------------------*/

	/** Starts the underlying thread for concurrent execution */
	public final void start(){thread.start();}
	/** Gets the current state of the underlying thread.
	 * @return The thread's current execution state */
	public final Thread.State getState(){return thread.getState();}
	/** Waits for the thread to complete execution.
	 * @throws InterruptedException If the current thread is interrupted while waiting */
	public final void join() throws InterruptedException{thread.join();}

	/*--------------------------------------------------------------*/
	
	/** Thread-local count of k-mers tested during exploration */
	long kmersTestedT=0;
	/** Thread-local count of dead ends discovered and removed */
	long deadEndsFoundT=0;
	/** Thread-local count of bubbles identified from branch intersections */
	long bubblesFoundT=0;
	
	/** Unique identifier for this thread within the exploration pool */
	final int id;
	final Kmer myKmer, myKmer2;

	/** Array for counting left-extending bases during k-mer exploration */
	final int[] leftCounts=KillSwitch.allocInt1D(4);
	/** Array for counting right-extending bases during k-mer exploration */
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	/** Thread-local byte buffer for constructing sequences during exploration */
	final ByteBuilder builderT=new ByteBuilder();

	/** Thread-local matrix tracking counts of each exploration code combination */
	long[][] countMatrixT=new long[MAX_CODE+1][MAX_CODE+1];
	/**
	 * Thread-local matrix tracking removal counts for each exploration code combination
	 */
	long[][] removeMatrixT=new long[MAX_CODE+1][MAX_CODE+1];
	
	/** The underlying Java thread object for concurrent execution */
	public final Thread thread;
	
}
