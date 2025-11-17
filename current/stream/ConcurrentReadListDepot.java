package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

/**
 * Thread-safe depot for managing a pool of reusable ArrayList buffers.
 * Maintains empty and full buffer queues using ArrayBlockingQueue for concurrent access.
 * Useful for producer-consumer patterns where multiple threads need access to
 * pre-allocated ArrayList instances.
 *
 * @author Brian Bushnell
 */
public class ConcurrentReadListDepot<K> {
	
	
	
	/**
	 * Constructs a ConcurrentReadListDepot with specified buffer configuration.
	 * Pre-allocates the specified number of ArrayList instances and places them
	 * in the empty queue ready for use.
	 *
	 * @param bufSize Initial capacity for each ArrayList buffer
	 * @param numBufs Number of ArrayList buffers to create and manage
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
	
	
	/** Queue of empty ArrayList buffers available for use */
	public final ArrayBlockingQueue<ArrayList<K>> empty;
	/** Queue of full ArrayList buffers ready for processing */
	public final ArrayBlockingQueue<ArrayList<K>> full;
	
	/** Initial capacity for each ArrayList buffer */
	public final int bufferSize;
	/** Total number of ArrayList buffers managed by this depot */
	public final int bufferCount;
	
	
	/** Array holding references to all managed ArrayList instances */
	private final ArrayList<K>[] lists;
	
}
