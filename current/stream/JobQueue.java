package stream;

import java.util.Comparator;
import java.util.PriorityQueue;

import shared.Tools;

/**
 * Priority queue using the heap monitor for queue and counter operations.
 * Every instance currently forces ordered mode, including when ordered_ is false.
 * Intended input has dense, complete IDs starting at firstID; readiness accepts a head
 * ID less than or equal to nextID, and each removal increments nextID once.
 * Job references are retained; callers must keep queued IDs and marker flags stable.
 *
 * Admission uses soft ticket/occupancy gating rather than a hard heap-size limit.
 * take() returns a LAST-marked job unless it is also poison, then stops on later calls.
 * Poison jobs map to null. Enqueuing a terminal marker does not bypass ID readiness;
 * the explicit forced-stop flag is handled separately. poll() and hasMore() have
 * their own documented state checks and are not interchangeable with take().
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 23, 2025
 * @param <K> Job type supplying stable identification and marker properties
 */
//TODO: Make high-speed version with 4 heaps using id()&3 to select heap and reduce lock contention
public class JobQueue<K extends HasID>{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a bounded, ordered queue whose first expected ID is zero. */
	public JobQueue(int capacity_){this(capacity_, true, true, 0);}
	/** Creates a bounded queue starting at zero; ordered_ is currently advisory only. */
	public JobQueue(int capacity_, boolean ordered_){this(capacity_, ordered_, true, 0);}
	
	/** Creates the queue and initializes its expected-ID and admission counters.
	 * Suggested capacity is 3+(1.5*threads); capacity is a soft window, not a heap limit.
	 * @param capacity_ Requested window; asserted greater than one, then clamped to at least two
	 * @param ordered_ Advisory flag; current implementation always selects ordered mode
	 * @param bounded_ Enables admission waits when true; false skips those waits
	 * @param firstID Initial expected ID; maxSeen starts at firstID-1 */
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
	
	/*--------------------------------------------------------------*/
	/*-------------------        Methods        --------------------*/
	/*--------------------------------------------------------------*/

	/** Retains and inserts a job, waiting according to the current admission policy.
	 * In bounded ordered mode, waits while id-capacity exceeds nextID, heap size exceeds
	 * quarter, and no forced stop is recorded. This permits the needed IDs through and
	 * does not enforce a hard maximum heap size. An interrupted wait restores interrupt
	 * status and leaves the wait loop before insertion. Unbounded mode skips admission waits.
	 * Updates maxSeen from job.id(); notifications follow the existing readiness policy.
	 * @param job Nonnull job whose ID and marker properties remain stable while queued
	 * @return true after insertion */
	public boolean add(K job){
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
			if(!bounded){//skip wait
			}else if(ordered){
//				while(bounded && heap.size()>=capacity && id>=nextID+half && id>heap.peek().id() && !poisoned){
				while(ticket>nextID && heap.size()>quarter && !poisoned){
					if(verbose2 && warn){
						warn=false;
						System.err.println(name+" Worker can't add "+id+": ticket "+ticket+">"+nextID);
					}
					try{
						heap.wait();
					}catch(InterruptedException e){
						Thread.currentThread().interrupt(); // Preserve interrupt status for caller
						break; //#001 fix: do NOT loop back into wait() with the interrupt flag re-armed - wait() then throws InterruptedException immediately every iteration -> 100% CPU busy-spin holding the heap monitor (the BgzfInputStreamMT2 clean-close race: worker parked here, interrupted by close(), poisoned still false). Stop waiting; fall through to heap.add(job). Capacity is relaxed only during this interrupted shutdown window (no data loss; interrupt status preserved for the caller).
					}
				}
			}else{
				while(heap.size()>=capacity && !poisoned){
					if(verbose2){System.err.println(name+" Worker can't add "+id+": size "+heap.size()+">="+capacity);}
					try{
						heap.wait();
					}catch(InterruptedException e){
						Thread.currentThread().interrupt(); // Preserve interrupt status for caller
						break; //#001 fix: do NOT loop back into wait() with the interrupt flag re-armed - wait() then throws InterruptedException immediately every iteration -> 100% CPU busy-spin holding the heap monitor (the BgzfInputStreamMT2 clean-close race: worker parked here, interrupted by close(), poisoned still false). Stop waiting; fall through to heap.add(job). Capacity is relaxed only during this interrupted shutdown window (no data loss; interrupt status preserved for the caller).
					}
				}
			}
			heap.add(job);
			maxSeen=Math.max(maxSeen, job.id());
			if(verbose2){
				System.err.println(name+" Worker: added job "+toString(job)+
					" to heap (heap size now "+heap.size()+")");
			}
			// Lazy notify: only wake consumer if this is the job they need or heap was empty
			if(id==nextID || (!ordered && heap.size()==1)){
				if(verbose2){System.err.println(name+" Worker notify.");}
				heap.notifyAll();
			}
		}
		return true;
	}
	
	/** Formats an ID and marker flags for diagnostics; null becomes the literal null string. */
	private final String toString(K k){
		if(k==null){return "null";}
		String s="id="+k.id();
		if(k.poison()){s+=" poison";}
		if(k.last()){s+=" last";}
		return s;
	}
	
	/** Waits for readiness, an observed LAST job, or a forced stop under the heap monitor.
	 * Returns a removed LAST job as data unless it is also poison; subsequent takes return
	 * null. A removed poison job is reinserted and maps to null. Every removal advances
	 * nextID once. A LAST job also supplies a poison marker through makePoison(id+1).
	 * Completion uses markers/forced-stop state, not consumer interruption; retain the
	 * documented historical interruption constraint in the implementation below.
	 * @return Ready job, including an ordinary LAST job, or null at a terminal condition */
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
					try{
						heap.wait();
					}catch(InterruptedException e){
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
				if(lastSeen || poisoned){return null;}
				job=heap.poll();
				if(verbose2){System.err.println(name+" Consumer fetched "+toString(job));}
				assert(job.id()<=nextID || !ordered); // Defensive check for ordering
				nextID++; // Advance to next expected ID
				lastSeen=lastSeen || job.last(); // Check for shutdown signal
				if(job.last()){
					if(verbose){System.err.println(name+" Consumer fetched last and added poison.");}
					heap.add((K)job.makePoison(job.id()+1));
				}else if(job.poison()){
					if(verbose){System.err.println(name+" Consumer fetched and reinserted poison.");}
					heap.add(job);
				}
				final int size=heap.size();
				// Lazy notify: wake producers only when necessary to reduce overhead
				// Skip notification when heap is mostly full (more jobs coming soon anyway)
				if(size==half || size==0 || (ordered && heap.peek().id()!=nextID) || job.poison() || job.last()){
					heap.notifyAll();
					if(verbose2){System.err.println(name+" Consumer notify.");}
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
	
	/** Removes a ready job without waiting and applies the same removal/marker accounting.
	 * Readiness is checked through heapReady(); this method does not first apply take's
	 * lastSeen/forced-stop exit checks. A removed LAST job is returned unless also poison.
	 * @return Ready job, or null when not ready or when a poison job was removed */
	public K poll(){
		K job=null;
		synchronized(heap){
			if(!heapReady()){return null;} // Return immediately if not ready
			
			job=heap.poll();
			if(verbose2){System.err.println(name+" Consumer polled "+toString(job));}
			assert(job.id()<=nextID || !ordered); 
			nextID++; 
			lastSeen=lastSeen || job.last();
			if(job.last()){
				if(verbose){System.err.println(name+" Consumer polled last and added poison.");}
				heap.add((K)job.makePoison(job.id()+1));
			}else if(job.poison()){
				if(verbose){System.err.println(name+" Consumer polled and reinserted poison.");}
				heap.add(job);
			}
			final int size=heap.size();

			if(size==half || size==0 || (ordered && heap.peek().id()!=nextID) || job.poison() || job.last()){
				heap.notifyAll();
				if(verbose2){System.err.println(name+" Consumer notify.");}
			}
		}
		return job==null || job.poison() ? null : job;
	}
	
	/** Tests head readiness under the monitor. Empty heaps report lastSeen; ordered heads
	 * are ready when {@code id<=nextID}. The inactive unordered branch excludes terminal markers. */
	private boolean heapReady(){
		synchronized(heap){
			if(heap.isEmpty()){return lastSeen;}
			K k=heap.peek();
			if(verbose2){System.err.println("heapReady found "+k.id()+"; nextID="+nextID);}
			if(k.id()<=nextID){return true;}//Poison may be lower than expected
			return !ordered && !k.last() && !k.poison();//TODO: Add a normal() function.
		}
	}
	
	/** Returns a synchronized snapshot of !lastSeen; does not inspect heap contents or forced-stop state. */
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
	
	/** Returns the synchronized expected-ID counter, advanced once for each removed job. */
	public long nextID(){synchronized(heap){return nextID;}}
	
	/** Returns the synchronized maximum recorded by add(), initially firstID-1. */
	public long maxSeen(){synchronized(heap){return maxSeen;}}
	
	//TODO: Possible bug - with force=false the pill's id (maxSeenId+1) sorts strictly LAST, so if
	//the stream has a genuine gap (a worker died holding job k), heapReady stays false at nextID=k
	//and the consumer never wakes despite the notifyAll: only force=true guarantees liveness after
	//job loss.  Fine under the graceful never-lose-a-job protocol; callers on ERROR paths must use
	//force=true.  ALSO NOTE (deliberate, verified): OQS shutdown can put LAST@m+1 and POISON@m+1
	//into the same heap - equal ids, arbitrary tie-break.  Both orders terminate correctly today
	//(traced 2026-07-14); do not "fix" the tie without re-tracing both interleavings.
	/** Inserts the supplied marker and notifies waiters without updating maxSeen.
	 * force=true latches the forced-stop state; false does not clear an existing latch.
	 * Ordinary marker insertion remains subject to the consumer's ID-readiness policy.
	 * @param poison Nonnull marker asserted to report poison()==true
	 * @param force Latch forced-stop state when true */
	public void poison(K poison, boolean force){
		assert(poison!=null && poison.poison()) : poison;
		synchronized(heap){
			if(verbose2){System.err.println(name+" poison().");}
			poisoned=poisoned || force;
			heap.add(poison);
			heap.notifyAll();
		}
	}
	
	/** Notifies all heap waiters without changing queue contents or readiness state. */
	public void notifyHeap(){synchronized(heap){heap.notifyAll();}}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes        -----------------*/
	/*--------------------------------------------------------------*/

	/** Compares job IDs only; supplies no additional tie-breaker. */
	private static class HasIDComparator<K extends HasID> implements Comparator<K>{

		/** Compares the two supplied job IDs. */
		@Override
		public int compare(K a, K b){return Long.compare(a.id(), b.id());}

	}
	
	/*--------------------------------------------------------------*/
	/*--------------------        Fields        --------------------*/
	/*--------------------------------------------------------------*/

	/** Caller-supplied prefix used by diagnostic messages. */
	public String name="";

	/** Next expected job ID in ordered mode */
	private long nextID;
	/** Highest ID recorded through add(), including its marker jobs; direct poison() insertion is excluded. */
	private long maxSeen;
	/** True after removing a job whose last() flag is set. */
	private boolean lastSeen=false;
	/** Forced-stop latch set by poison(..., true); consulted by add() and take(). */
	private boolean poisoned=false;
	/** Priority queue storing jobs, ordered by ID */
	private final PriorityQueue<K> heap;
	/** Ordered readiness selection; currently forced true by every constructor. */
	private final boolean ordered;
	/** Enables soft admission waiting when true. */
	private final boolean bounded;
	/** Window subtracted from an incoming ID to form its admission ticket. */
	private final int capacity;
	/** (capacity+1)/2, used as a lazy-notification threshold. */
	private final int half;
	/** (half+1)/2; ordered admission waits require heap size greater than this threshold. */
	private final int quarter;

	/*--------------------------------------------------------------*/
	/*------------------        Constants        -------------------*/
	/*--------------------------------------------------------------*/

	/** Enables diagnostic output for selected terminal events. */
	private static final boolean verbose=false;//Should be for important events like thread death
	/** Enables detailed queue-operation diagnostics. */
	private static final boolean verbose2=false;
	
}
