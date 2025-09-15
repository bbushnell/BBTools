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
	
	abstract boolean processNextTable();
	
	/*--------------------------------------------------------------*/
	
	long kmersRemovedT=0;
	
	final int id;
	
}