package stream;

import java.util.Comparator;
import java.util.PriorityQueue;

import shared.Tools;

/**
 * Thread-safe job queue with optional ordering and capacity bounds.
 * Simplifies multithreaded producer-consumer patterns by handling synchronization,
 * ordering, and backpressure automatically.
 * 
 * Supports two primary modes:
 * - Ordered: Jobs are retrieved in sequential ID order, using the heap as a reordering buffer
 * - Unordered: Jobs are retrieved as available, prioritized by ID but not strictly ordered
 * 
 * Bounded mode prevents memory issues by blocking producers when the queue reaches capacity,
 * while still allowing jobs matching nextID to be added immediately to prevent deadlocks.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 23, 2025
 * 
 * @param <K> Job type implementing HasID for identification and ordering
 */
//TODO: Make high-speed version with 4 heaps using id()&3 to select heap and reduce lock contention
public class JobQueue<K extends HasID>{
	
	public JobQueue(int capacity_){this(capacity_, true, true, 0);}
	public JobQueue(int capacity_, boolean ordered_){this(capacity_, ordered_, true, 0);}
	
	/**
	 * Creates a new JobQueue with specified behavior.  Suggested capacity is 1+(2.5*threads).
	 * 
	 * @param capacity_ Maximum number of jobs allowed in queue before blocking producers (must be >1)
	 * @param ordered_ If true, jobs are released in strict ID order; if false, jobs released as available
	 * @param bounded_ If true, producers block when queue reaches capacity; if false, unbounded growth
	 * @param firstID Expected ID of the first job (typically 0)
	 */
	public JobQueue(int capacity_, boolean ordered_, boolean bounded_, long firstID){
		assert(capacity_>1) : "Capacity is too small: "+capacity_;
		capacity=Math.max(capacity_, 2);
		half=(capacity+1)/2; // Used for lazy notification optimization
		ordered=ordered_;
		bounded=bounded_;
		nextID=firstID;
		maxSeen=firstID-1;
		heap=new PriorityQueue<K>(Tools.mid(1, capacity, 96), new HasIDComparator<K>());
	}
	
	/**
	 * Adds a job to the queue, blocking if necessary to respect capacity bounds.
	 * In bounded mode, blocks when heap is at capacity UNLESS this job's ID matches nextID,
	 * which prevents deadlocks by allowing the exact job the consumer needs to be added.
	 * 
	 * @param job Job to add to the queue
	 * @throws InterruptedException if interrupted while waiting for capacity
	 * @return True when the add is successful.
	 */
	public boolean add(K job) {
		final long id=job.id();
		synchronized(heap){
			// Block if bounded, at capacity, and this isn't the job the consumer is waiting for
			// The id>heap.peek().id() check prevents deadlock by letting nextID through
			while(bounded && heap.size()>=capacity && id>=nextID+half && id>heap.peek().id()){
				try {
					heap.wait();
				} catch (InterruptedException e){
					Thread.currentThread().interrupt(); // Preserve interrupt status for caller
				}
			}
			heap.add(job);
			maxSeen=Math.max(maxSeen, job.id());
			if(verbose){
				System.err.println("Worker: added job " + id +
					" to heap (heap size now " + heap.size() + ")");
			}
			// Lazy notify: only wake consumer if this is the job they need or heap was empty
			if(id==nextID || (!ordered && heap.size()==1)){heap.notifyAll();}
		}
		return true;
	}
	
	/**
	 * Retrieves the next job from the queue, waiting if necessary.
	 * In ordered mode, waits for jobs in strict sequential ID order.
	 * In unordered mode, returns jobs as they become available.
	 * Returns null after receiving a job marked as last().
	 * 
	 * @return Next job to process, or null if processing is complete
	 */
	public K take(){
		K job=null;
		if(verbose){System.err.println("Consumer waiting for "+nextID);}
		synchronized(heap){
			while(job==null && !lastSeen){
				// Wait if heap is empty or (in ordered mode) next job isn't ready yet
				while(heap.isEmpty() || (ordered && heap.peek().id()>nextID)){
					if(verbose){System.err.println("Consumer waiting; heap.size()="+heap.size());}
					try {
						heap.wait();
					} catch (InterruptedException e){
						Thread.currentThread().interrupt(); // Preserve interrupt status
						// Don't return null here - wait for explicit last signal
					}
				}
				job=heap.poll();
				if(verbose){System.err.println("Consumer fetched "+job.id());}
				assert(job.id()<=nextID || !ordered); // Defensive check for ordering
				nextID++; // Advance to next expected ID
				lastSeen=lastSeen || job.last(); // Check for shutdown signal
				if(job.last()) {heap.add((K)job.makePoison(job.id()+1));}
				else if(job.poison()) {heap.add(job);}
				final int size=heap.size();
				// Lazy notify: wake producers only when necessary to reduce overhead
				// Skip notification when heap is mostly full (more jobs coming soon anyway)
				if(size==half || size==0 || (ordered && heap.peek().id()!=nextID) || job.poison() || job.last()){
					heap.notifyAll();
				}
			}
		}
		return job;
	}
	
	public boolean hasMore(){
		synchronized(heap){return !lastSeen;}
	}
	
	public long nextID(){
		synchronized(heap){return nextID;}
	}
	
	public long maxSeen(){
		synchronized(heap){return maxSeen;}
	}

	/** Comparator for sorting jobs */
	private static class HasIDComparator<K extends HasID> implements Comparator<K> {

		@Override
		public int compare(K a, K b){return Long.compare(a.id(), b.id());}

	}

	/** Next expected job ID in ordered mode */
	private long nextID;
	/** Highest ID seen */
	private long maxSeen;
	/** True once a job marked last() has been seen */
	private boolean lastSeen=false;
	/** Priority queue storing jobs, ordered by ID */
	private final PriorityQueue<K> heap;
	/** If true, jobs are released in strict ID order */
	private final boolean ordered;
	/** If true, producers block when queue reaches capacity */
	private final boolean bounded;
	/** Maximum jobs allowed in queue before blocking */
	private final int capacity;
	/** Half of capacity, used for lazy notification optimization */
	private final int half;//TODO: Change so high numbers can only add under half
	/** Enable debug output */
	private static final boolean verbose=false;
	
}