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

	boolean processNextTable(){return processNextTable(myKmer, myKmer2);}
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
	final Kmer myKmer, myKmer2;

	final int[] leftCounts=KillSwitch.allocInt1D(4);
	final int[] rightCounts=KillSwitch.allocInt1D(4);
	final ByteBuilder builderT=new ByteBuilder();

	long[][] countMatrixT=new long[MAX_CODE+1][MAX_CODE+1];
	long[][] removeMatrixT=new long[MAX_CODE+1][MAX_CODE+1];
	
	public final Thread thread;
	
}
