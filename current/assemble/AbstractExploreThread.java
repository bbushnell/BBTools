package assemble;

import shared.KillSwitch;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Abstract thread class for exploring k-mer graph structures to identify dead ends and bubbles.
 * Part of the assembly shaving pipeline that removes spurious branches from k-mer graphs.
 * Extends ShaveObject to inherit exploration codes and branch constants.
 *
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
abstract class AbstractExploreThread extends ShaveObject implements Runnable {

	/**
	 * Constructs an exploration thread with specified ID and k-mer size.
	 * Initializes k-mer working objects and creates the underlying thread.
	 * @param id_ Thread identifier for tracking and coordination
	 * @param kbig_ K-mer size for the graph structures being explored
	 */
	public AbstractExploreThread(int id_, int kbig_){
		id=id_;
		myKmer=new Kmer(kbig_);
		myKmer2=new Kmer(kbig_);
		thread=new Thread(this);
	}

	/** Main thread execution that iteratively processes tables and victim lists,
	 * then accumulates bubble statistics from removal matrices after exploration. */
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
	/**
	 * Processes the next k-mer table section for dead end detection.
	 * Abstract method that must be implemented by concrete subclasses.
	 *
	 * @param kmer Primary k-mer object for graph traversal
	 * @param temp Temporary k-mer object for computations
	 * @return true if more table processing is needed, false if complete
	 */
	abstract boolean processNextTable(final Kmer kmer, Kmer temp);

	boolean processNextVictims(){return processNextVictims(myKmer, myKmer);}
	abstract boolean processNextVictims(final Kmer kmer, Kmer temp);

	/*--------------------------------------------------------------*/

	public final void start(){thread.start();}
	public final Thread.State getState(){return thread.getState();}
	public final void join() throws InterruptedException{thread.join();}

	/*--------------------------------------------------------------*/
	
	long kmersTestedT=0;
	long deadEndsFoundT=0;
	long bubblesFoundT=0;
	
	final int id;
	/** Secondary thread-local k-mer used for temporary computations. */
	final Kmer myKmer, myKmer2;

	final int[] leftCounts=KillSwitch.allocInt1D(4);
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	final ByteBuilder builderT=new ByteBuilder();

	long[][] countMatrixT=new long[MAX_CODE+1][MAX_CODE+1];
	long[][] removeMatrixT=new long[MAX_CODE+1][MAX_CODE+1];
	
	public final Thread thread;
	
}
