package assemble;

/**
 * @author Brian Bushnell
 * @date Jul 20, 2015
 *
 */
abstract class AbstractShaveThread extends Thread{

	/** Creates a new shave thread with the specified identifier.
	 * @param id_ Thread identifier for tracking and debugging */
	public AbstractShaveThread(int id_){
		id=id_;
	}
	
	/** Main thread execution method that processes k-mer tables until exhausted.
	 * Continuously calls processNextTable() until no more tables remain. */
	@Override
	public final void run(){
		while(processNextTable()){}
	}
	
	/**
	 * Processes the next available k-mer table, removing dead-end k-mers.
	 * Implementation varies by subclass based on specific k-mer table structure.
	 * @return true if a table was processed, false if no more tables available
	 */
	abstract boolean processNextTable();
	
	/*--------------------------------------------------------------*/
	
	long kmersRemovedT=0;
	
	final int id;
	
}