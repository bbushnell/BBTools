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
 * CONTRACT: Only supports ordered streams with dense, complete IDs (firstID, firstID+1, ...,
 * no skipped IDs).  A missing ID stalls the consumer until a poison/last pill arrives.
 * Producers of sparse streams must insert empty placeholder jobs for filtered-out IDs.
 * The ordered_ constructor parameter is currently ADVISORY-ONLY (see constructor): every
 * queue runs in ordered mode, because unordered mode violated assumptions elsewhere.
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
	 * Creates a new JobQueue with specified behavior.  Suggested capacity is 3+(1.5*threads).
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
		quarter=(half+1)/2;//Anyone can add under quarter full
		ordered=ordered_ || true;//TODO: Review all cases where this can legitimately be set to false.
		//TODO: Possible bug - ordered_ is silently ignored (forced true).  Callers passing false
		//(OQS 'ordered', OQS2 'orderedOutput') get strict ordering plus its hidden dense-ID
		//requirement: any skipped ID hangs the consumer until poison.  Intent was "ordered=false
		//makes things cheaper", but correctness was hard to ensure; disordered mode violated
		//assumptions elsewhere and caused hangs (Brian).  Safe reading: ordered=true means "I
		//need ordering"; ordered=false does NOT mean "gap-tolerant".  Class javadoc updated.
		bounded=bounded_;
		nextID=firstID;
		maxSeen=firstID-1;
		heap=new PriorityQueue<K>(Tools.mid(1, capacity, 96), new HasIDComparator<K>());
	}
	
	/**
	 * Adds a job to the queue, blocking if necessary to respect capacity bounds.
	 * In bounded mode, blocks until id is within capacity of nextID.
	 * The capacity bound is SOFT backpressure, not a hard limit: if the calling thread is
	 * INTERRUPTED while waiting for capacity, it stops waiting and adds the job anyway (relaxing
	 * the bound by at most one job per interrupted thread). Interrupts are treated as a shutdown
	 * signal across all callers, so this never grows unbounded. (Does NOT throw; the interrupt
	 * status is preserved for the caller.)
	 *
	 * @param job Job to add to the queue
	 * @return True when the add is successful (always, once it returns).
	 */
	public boolean add(K job) {
		final long id=job.id();
		final long ticket=id-capacity;
		boolean warn=verbose2;
		if(verbose2){System.err.println(name+" Worker: got ticket "+ticket+" for job "+id);}
		synchronized(heap){
			//Old version:
			// Block if bounded, at capacity, and this isn't the job the consumer is waiting for
			// The id>heap.peek().id() check prevents deadlock by letting nextID through
			// while(bounded && heap.size()>=capacity && id>=nextID+half && id>heap.peek().id() && !poisoned){
			
			//New version:  Take a ticket.
			if(!bounded) {//skip wait
			}else if(ordered) {
//				while(bounded && heap.size()>=capacity && id>=nextID+half && id>heap.peek().id() && !poisoned){
				while(ticket>nextID && heap.size()>quarter && !poisoned){
					if(verbose2 && warn) {
						warn=false;
						System.err.println(name+" Worker can't add "+id+": ticket "+ticket+">"+nextID);
					}
					try {
						heap.wait();
					} catch (InterruptedException e){
						Thread.currentThread().interrupt(); // Preserve interrupt status for caller
						break; //#001 fix: do NOT loop back into wait() with the interrupt flag re-armed - wait() then throws InterruptedException immediately every iteration -> 100% CPU busy-spin holding the heap monitor (the BgzfInputStreamMT2 clean-close race: worker parked here, interrupted by close(), poisoned still false). Stop waiting; fall through to heap.add(job). Capacity is relaxed only during this interrupted shutdown window (no data loss; interrupt status preserved for the caller).
					}
				}
			}else {
				while(heap.size()>=capacity && !poisoned){
					if(verbose2) {System.err.println(name+" Worker can't add "+id+": size "+heap.size()+">="+capacity);}
					try {
						heap.wait();
					} catch (InterruptedException e){
						Thread.currentThread().interrupt(); // Preserve interrupt status for caller
						break; //#001 fix: do NOT loop back into wait() with the interrupt flag re-armed - wait() then throws InterruptedException immediately every iteration -> 100% CPU busy-spin holding the heap monitor (the BgzfInputStreamMT2 clean-close race: worker parked here, interrupted by close(), poisoned still false). Stop waiting; fall through to heap.add(job). Capacity is relaxed only during this interrupted shutdown window (no data loss; interrupt status preserved for the caller).
					}
				}
			}
			heap.add(job);
			maxSeen=Math.max(maxSeen, job.id());
			if(verbose2){
				System.err.println(name+" Worker: added job " + toString(job) +
					" to heap (heap size now " + heap.size() + ")");
			}
			// Lazy notify: only wake consumer if this is the job they need or heap was empty
			if(id==nextID || (!ordered && heap.size()==1)){
				if(verbose2) {System.err.println(name+" Worker notify.");}
				heap.notifyAll();
			}
		}
		return true;
	}
	
	private final String toString(K k) {
		if(k==null) {return "null";}
		String s="id="+k.id();
		if(k.poison()) {s+=" poison";}
		if(k.last()) {s+=" last";}
		return s;
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
		if(verbose2){System.err.println(name+" Consumer waiting for "+nextID);}
		synchronized(heap){
			while(job==null && !lastSeen && !poisoned){
				// Wait if heap is empty or (in ordered mode) next job isn't ready yet
				while(!heapReady() && !lastSeen && !poisoned){
					if(verbose2){
						System.err.println(name+" Consumer waiting for ("+nextID+"); heap.size()="+heap.size()+
							(heap.isEmpty() ? "" : ": "+toString(heap.peek())));
					}
					try {
						heap.wait();
					} catch (InterruptedException e){
						Thread.currentThread().interrupt(); // Preserve interrupt status
						// Don't return null here - wait for explicit last signal
						// CONTRACT/HAZARD (Furina 2026-06-25): this DELIBERATELY ignores interrupts and keeps waiting
						// for a real terminal (heapReady / lastSeen / poisoned). Note this is the SAME re-arm-then-
						// loop-back-to-wait() pattern that busy-spun in add() (#001 fix): if a thread is ever
						// interrupted while parked here, wait() will throw immediately every iteration -> 100% CPU
						// RUNNABLE spin holding the heap monitor. It is left AS-IS ON PURPOSE because the spin is
						// currently UNREACHABLE: no JobQueue.take() caller is ever interrupted while in take()
						// (verified across all users - interrupts target producers/workers, which call add(), and the
						// waitForFinish-join thread, never the take()-consumer). To unblock a parked consumer, the
						// correct mechanism is poison()/last (which this loop's !poisoned/!lastSeen conditions honor),
						// NEVER an interrupt. If a future caller needs to interrupt a consumer, do NOT just add a
						// `break` here (it would fall through to heap.poll() on a not-ready heap -> NPE or ordering
						// violation); make it poison-driven instead, or apply the swallow-don't-re-arm variant.
					}
				}
				if(lastSeen || poisoned) {return null;}
				job=heap.poll();
				if(verbose2){System.err.println(name+" Consumer fetched "+toString(job));}
				assert(job.id()<=nextID || !ordered); // Defensive check for ordering
				nextID++; // Advance to next expected ID
				lastSeen=lastSeen || job.last(); // Check for shutdown signal
				if(job.last()) {
					if(verbose) {System.err.println(name+" Consumer fetched last and added poison.");}
					heap.add((K)job.makePoison(job.id()+1));
				}else if(job.poison()) {
					if(verbose) {System.err.println(name+" Consumer fetched and reinserted poison.");}
					heap.add(job);
				}
				final int size=heap.size();
				// Lazy notify: wake producers only when necessary to reduce overhead
				// Skip notification when heap is mostly full (more jobs coming soon anyway)
				if(size==half || size==0 || (ordered && heap.peek().id()!=nextID) || job.poison() || job.last()){
					heap.notifyAll();
					if(verbose2) {System.err.println(name+" Consumer notify.");}
				}
			}
		}
		//TODO: Possible bug - a last() marker is RETURNED to the consumer as a normal job here
		//(only poison maps to null; the null arrives on the NEXT call).  Every consumer must
		//either tolerate the makeLast() payload as processable (empty) cargo or check job.last()
		//itself - nothing in getOutput()/take() docs says so.  If any job type's makeLast()
		//carries non-empty or half-initialized payload, it gets processed as data.
		return job==null || job.poison() ? null : job;
	}
	
	/**
	 * Retrieves the next job if it is ready (in order). 
	 * Returns null immediately if the queue is empty or the next ordered ID is missing.
	 * Non-blocking version of take().
	 */
	public K poll(){
		K job=null;
		synchronized(heap){
			if(!heapReady()){return null;} // Return immediately if not ready
			
			job=heap.poll();
			if(verbose2){System.err.println(name+" Consumer polled "+toString(job));}
			assert(job.id()<=nextID || !ordered); 
			nextID++; 
			lastSeen=lastSeen || job.last();
			if(job.last()) {
				if(verbose) {System.err.println(name+" Consumer polled last and added poison.");}
				heap.add((K)job.makePoison(job.id()+1));
			}else if(job.poison()) {
				if(verbose) {System.err.println(name+" Consumer polled and reinserted poison.");}
				heap.add(job);
			}
			final int size=heap.size();

			if(size==half || size==0 || (ordered && heap.peek().id()!=nextID) || job.poison() || job.last()){
				heap.notifyAll();
				if(verbose2) {System.err.println(name+" Consumer notify.");}
			}
		}
		return job==null || job.poison() ? null : job;
	}
	
	private boolean heapReady() {
		synchronized(heap) {
			if(heap.isEmpty()) {return lastSeen;}
			K k=heap.peek();
			if(verbose2) {System.err.println("heapReady found "+k.id()+"; nextID="+nextID);}
			if(k.id()<=nextID) {return true;}//Poison may be lower than expected
			return !ordered && !k.last() && !k.poison();//TODO: Add a normal() function.
		}
	}
	
	public boolean hasMore(){
		//TODO: Possible bug - ignores 'poisoned'.  After a FORCED shutdown (poison(pill, force=true),
		//e.g. writer error paths calling setFinished(true)), take() returns null forever while
		//hasMore() stays true - a `while(hasMore()){process(getOutput());}` consumer spins or NPEs
		//on the exact path where a thread just died.  Candidate fix is `!lastSeen && !poisoned`,
		//but that changes visible behavior for every user under nondeterministic timing, so it
		//needs real multi-workflow validation, not a drive-by edit.  NOTE: the GRACEFUL path is
		//unaffected - poisoned is only ever set by force=true, so with the intended never-forced
		//protocol hasMore()==!lastSeen is correct.
		synchronized(heap){return !lastSeen;}
	}
	
	public long nextID(){
		synchronized(heap){return nextID;}
	}
	
	public long maxSeen(){
		synchronized(heap){return maxSeen;}
	}
	
	//TODO: Possible bug - with force=false the pill's id (maxSeenId+1) sorts strictly LAST, so if
	//the stream has a genuine gap (a worker died holding job k), heapReady stays false at nextID=k
	//and the consumer never wakes despite the notifyAll: only force=true guarantees liveness after
	//job loss.  Fine under the graceful never-lose-a-job protocol; callers on ERROR paths must use
	//force=true.  ALSO NOTE (deliberate, verified): OQS shutdown can put LAST@m+1 and POISON@m+1
	//into the same heap - equal ids, arbitrary tie-break.  Both orders terminate correctly today
	//(traced 2026-07-14); do not "fix" the tie without re-tracing both interleavings.
	public void poison(K poison, boolean force) {
		assert(poison!=null && poison.poison()) : poison;
		synchronized(heap){
			if(verbose2) {System.err.println(name+" poison().");}
			poisoned=poisoned||force;
			heap.add(poison);
			heap.notifyAll();
		}
	}
	
	public void notifyHeap() {
		synchronized(heap) {heap.notifyAll();}
	}

	/** Comparator for sorting jobs */
	private static class HasIDComparator<K extends HasID> implements Comparator<K> {

		@Override
		public int compare(K a, K b){return Long.compare(a.id(), b.id());}

	}
	
	public String name="";

	/** Next expected job ID in ordered mode */
	private long nextID;
	/** Highest ID seen */
	private long maxSeen;
	/** True once a job marked last() has been seen */
	private boolean lastSeen=false;
	/** Threads interacting with this should shut down */
	private boolean poisoned=false;
	/** Priority queue storing jobs, ordered by ID */
	private final PriorityQueue<K> heap;
	/** If true, jobs are released in strict ID order */
	private final boolean ordered;
	/** If true, producers block when queue reaches capacity */
	private final boolean bounded;
	/** Maximum jobs allowed in queue before blocking */
	private final int capacity;
	/** Half of capacity, used for lazy notification optimization */
	private final int half;
	private final int quarter;
	/** Enable debug output */
	private static final boolean verbose=false;//Should be for important events like thread death
	private static final boolean verbose2=false;
	
}