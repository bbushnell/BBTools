package fileIO;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;

/**
 * Creates one thread per serialized-object load and limits how many reads are admitted.
 * Waiting requests still have their own threads; this is not a reusable thread pool
 * or a bound on object size or aggregate memory. LOW_MEMORY selects one admitted
 * reader; otherwise LIMIT applies. Configure these settings before starting loads.
 * Successful reads assign output and release their counters. Failures are not
 * propagated to callers here and can leave counters registered after a worker exits.
 * Counter-based waiting and per-thread waiting have different completion guarantees.
 *
 * @author Brian Bushnell
 * @date Jan 2, 2013
 */
public class LoadThread<X> extends Thread{
	
	/** Registers a request, starts its thread and returns the handle.
	 * Admission to the read happens inside run(), after the new thread starts.
	 * @param <Y> Requested static result type; no runtime type check is performed here
	 * @param fname Input name interpreted by ReadWrite.read
	 * @param c Type token forwarded to ReadWrite.read, which currently ignores it
	 * @return Started load handle; output may still be null after termination
	 */
	public static <Y> LoadThread<Y> load(String fname, Class<Y> c){
		LoadThread<Y> lt=new LoadThread<Y>(fname, c);
		//TODO: Probable bug - FIO-070: a failed start leaves the constructor's waiting registration behind.
		lt.start();
		return lt;
	}
	
	/** Stores the request and counts it as waiting before the factory starts its thread. */
	private LoadThread(String fname_, Class<X> c_){
		fname=fname_;
		c=c_;
		addThread(1);
	}
	
	/** Waits for admission, reads with allowSubprocess=false and assigns the result.
	 * Successful null results are allowed. Assignment precedes counter release and
	 * this thread's final notification; no exception/result envelope is provided.
	 */
	@Override
	public void run(){
		addRunningThread(1);
		//TODO: Probable bug - FIO-070: escaping read failures skip counter release and notification; see align2.Block.read's dependency TODO.
		output=ReadWrite.read(c, fname, false);
		addRunningThread(-1);
		synchronized(this){this.notify();}
	}
	
	
	/** Registers a waiting request for the internal +1 transition; -1 releases a running request.
	 * @return Total registered requests after the change
	 */
	private static final int addThread(int x){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				activeThreads[0]+=x;
				activeThreads[1]+=x;
			}else{
				addRunningThread(x);
			}
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
					
			return activeThreads[0];
		}
	}
	
	/** Admits one request for the internal +1 transition, or releases one for -1.
	 * Counter updates and wakeups share the activeThreads monitor. Interrupted waits
	 * print the interruption and retry; they do not cancel a request.
	 * @return Running count after the change
	 */
	private static final int addRunningThread(int x){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(x!=0);
			if(x>0){
				assert(activeThreads[1]>=x);
				while(activeThreads[2]>=lim){
					try{
						activeThreads.wait();
					}catch(InterruptedException e){
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				activeThreads[1]-=x; //Remove from waiting
			}else{
				activeThreads[0]+=x; //Remove from active
			}
			activeThreads[2]+=x; //Change number running
			
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
			
			if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notifyAll();}
//			System.err.println(activeThreads[2]);
//			try {
//				activeThreads.wait(5000);
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
			return activeThreads[2];
		}
	}
	
	/** Reads the registered waiting-plus-running count under its monitor.
	 * @return Registered requests, which may include stale entries after failures
	 */
	public static final int countActiveThreads(){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
					activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
			return activeThreads[0];
		}
	}
	
	/** Waits until the registered count reaches zero; does not join each thread.
	 * Successful workers decrement the count before their final notification/exit.
	 * Stale failure registrations can prevent this wait from returning. Interrupted
	 * waits are printed and retried; the eight-second interval is a polling wait,
	 * not an overall timeout or failure detector.
	 */
	public static final void waitForReadingToFinish(){
		final int lim=(Shared.LOW_MEMORY ? 1 : LIMIT);
		synchronized(activeThreads){
			while(activeThreads[0]>0){
				assert(activeThreads[0]==(activeThreads[1]+activeThreads[2]) && activeThreads[0]>=0 && activeThreads[1]>=0 &&
						activeThreads[2]>=0 && activeThreads[2]<=lim) : Arrays.toString(activeThreads);
				try{
					activeThreads.wait(8000);
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(activeThreads[2]==0 || (activeThreads[2]<lim && activeThreads[1]>0)){activeThreads.notifyAll();}
			}
		}
	}
	
	/** Joins until termination only when the initial output check observes null.
	 * A nonnull output can cause return before counter release or thread termination;
	 * that plain-field check supplies no synchronization of its own. Null can remain
	 * after a successful null read or worker failure. Exceptions are not rethrown to
	 * the caller, and interrupted joins are printed and retried.
	 */
	public final void waitForThisToFinish(){
		if(output==null){
			while(this.getState()!=State.TERMINATED){
				try{
					this.join();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/** Shared [total, waiting, running] registrations, guarded by this array's monitor.
	 * Total equals waiting plus running. External code must not replace or mutate
	 * this bookkeeping during use; the limit covers reads, not all created threads.
	 */
	//FIXED [fileIO/LoadThread#002]: single notify() -> notifyAll() (both sites). Worker threads (addRunningThread,
	//no timeout) and the main thread (waitForReadingToFinish, 8s timeout) both wait on this monitor; a finishing
	//thread's single notify() could wake the main thread (which can't progress while workers are queued) instead of a
	//ready worker, stalling that worker until the 8s timeout re-notified. notifyAll() wakes all; those that can't
	//proceed re-wait. Not a deadlock and output was always correct (the timeout self-heals), so this is a latency fix.
	//Same twin pattern fixed in ReadWrite.addRunningThread/waitForWritingToFinish. -Furina flagged 2026-06-16, G11 fixed 2026-07-04
	public static int[] activeThreads={0, 0, 0};

	/** Input name passed to ReadWrite.read. */
	private final String fname;
	/** Forwarded type token; the current read helper does not validate against it. */
	private final Class<X> c;
	/** Plain result field assigned before counter release; null also represents a null result or failure. */
	public X output=null;
	
	/** Read-admission limit when LOW_MEMORY is false; configure before starting loads. */
	public static int LIMIT=Tools.min(12, Tools.max(Shared.threads(), 1));
	
}
