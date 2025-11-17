package assemble;

/**
 * @author Brian Bushnell
 * @date Jul 20, 2015
 *
 */
/**
 * Removes dead-end kmers.
 */
abstract class AbstractShaveThread extends Thread{

	/**
	 * Constructor
	 */
	public AbstractShaveThread(int id_){
		id=id_;
	}
	
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
	
	/** Count of k-mers removed by this thread during processing */
	long kmersRemovedT=0;
	
	/** Unique identifier for this thread instance */
	final int id;
	
}