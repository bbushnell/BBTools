package stream;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

/**
 * Thread-safe depot for managing concurrent access to pre-allocated buffer pools.
 * Uses blocking queues to coordinate empty and full buffer exchanges between
 * producer and consumer threads. Designed for high-throughput scenarios where
 * buffer allocation overhead needs to be minimized.
 *
 * @author Brian Bushnell
 * @param <K> Type of objects stored in the buffer lists
 */
public class ConcurrentDepot<K> {
	
	/**
	 * Creates a depot with pre-allocated buffer pools for concurrent access.
	 * Initializes empty queue with all buffers and full queue as empty.
	 * @param bufSize Initial capacity for each buffer list
	 * @param numBufs Number of buffer lists to pre-allocate
	 */
	@SuppressWarnings("unchecked")
	public ConcurrentDepot(int bufSize, int numBufs){
		bufferSize=bufSize;
		bufferCount=numBufs;
		
		lists=new ArrayList[numBufs];
		empty=new ArrayBlockingQueue<ArrayList<K>>(numBufs+1, fair);
		full=new ArrayBlockingQueue<ArrayList<K>>(numBufs+1, fair);
		
		for(int i=0; i<lists.length; i++){
			lists[i]=new ArrayList<K>(bufSize);
			empty.add(lists[i]);
		}
		
	}
	
	
	/** Queue of empty buffer lists available for producers to fill */
	public final ArrayBlockingQueue<ArrayList<K>> empty;
	/** Queue of full buffer lists ready for consumers to process */
	public final ArrayBlockingQueue<ArrayList<K>> full;
	
	/** Initial capacity allocated for each buffer list */
	public final int bufferSize;
	/** Total number of buffer lists in the depot */
	public final int bufferCount;
	
	/** Whether blocking queues use fair ordering for waiting threads */
	public static boolean fair=false;
	
	/** Internal array holding references to all pre-allocated buffer lists */
	private final ArrayList<K>[] lists;
	
}
