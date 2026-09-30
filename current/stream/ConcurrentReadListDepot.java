package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

/**
 * Preallocates reusable ArrayList buffers and exposes empty/ready blocking queues.
 * Queue operations support concurrent access; mutable list contents and multi-step
 * ownership transfers remain the caller's responsibility. This class does not clear
 * lists, move them between queues, validate membership or implement a terminal protocol.
 * Empty/full are caller-maintained roles, not checks on list contents. Callers may
 * insert other lists; the private array retains references only to the original pool.
 *
 * @author Brian Bushnell
 * @param <K> Element type stored in the buffer lists
 */
public class ConcurrentReadListDepot<K>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates numBufs distinct empty lists and places them in empty; full starts empty.
	 * Each queue has capacity numBufs+1. List capacity is an allocation hint, not a size
	 * limit, and this class performs no additional argument validation.
	 * @param bufSize Initial capacity of each ArrayList
	 * @param numBufs Number of lists initially allocated, not a membership limit
	 */
	public ConcurrentReadListDepot(int bufSize, int numBufs){
		bufferSize=bufSize;
		bufferCount=numBufs;

		lists=new ArrayList[numBufs];
		empty=new ArrayBlockingQueue<ArrayList<K>>(numBufs+1);
		full=new ArrayBlockingQueue<ArrayList<K>>(numBufs+1);

		for(int i=0; i<lists.length; i++){
			lists[i]=new ArrayList<K>(bufSize);
			empty.add(lists[i]);
		}

	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Available-list queue, initially containing all newly allocated empty lists. */
	public final ArrayBlockingQueue<ArrayList<K>> empty;
	/** Ready-list queue, initially empty; enqueued lists need not be filled to capacity. */
	public final ArrayBlockingQueue<ArrayList<K>> full;

	/** Initial list capacity, not a maximum list size. */
	public final int bufferSize;
	/** Initial allocation count; not current queue occupancy or enforced membership. */
	public final int bufferCount;

	/** Retained references to the initially allocated lists, regardless of their location. */
	private final ArrayList<K>[] lists;

}
