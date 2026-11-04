package prok;

import java.util.ArrayList;

import stream.Streamer;
import template.ThreadWaiter;

/** Starts an owned producer cohort and joins every started worker before return.
 * Callers inspect worker success and accumulate results after joining. Declare
 * this scope before borrowed reader/writer owners, so its close restores observed
 * interruption only after native I/O cleanup. It owns no streams or global settings.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaWorkers implements AutoCloseable{

	/** Captures interruption before native I/O operations can consume it. */
	public NcrnaWorkers(){restoreInterrupt=Thread.interrupted();}

	/** Starts fresh workers, preserving the original startup failure and interrupts.
	 * If none starts, drains the optional already-started FASTA input without doing
	 * feature/alignment work. The selected ZT/ST readers expose real terminal
	 * consumption through hasMore; an interrupted nonterminal null must be retried.
	 * With a partial cohort, started consumers finish before outputs can close.
	 * Native input startup failure is outside this consumer-lifetime contract.
	 * @param threads Fresh, nonnull workers owned by this invocation
	 * @param input Borrowed started FASTA stream, or null for nonstreaming workers
	 */
	public void startAndJoin(final ArrayList<? extends Thread> threads, final Streamer input){
		assert(threads!=null && !threads.isEmpty()) : "The validated worker budget creates a nonempty owned cohort";
		assert(Thread.currentThread()==owner) : "The worker scope retains interruption for its creating thread: "+owner.getName();
		restoreInterrupt|=Thread.interrupted();
		try{
			ThreadWaiter.startThreads(threads);
		}catch(RuntimeException | Error failure){
			boolean started=false;
			for(final Thread worker:threads){started|=worker.getState()!=Thread.State.NEW;}
			if(!started && input!=null){
				try{
					while(true){
						if(input.nextList()==null){
							if(!input.hasMore()){break;}
							restoreInterrupt=true;
						}
						restoreInterrupt|=Thread.interrupted();
					}
				}catch(RuntimeException | Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
			}
			throw failure;
		}finally{
			restoreInterrupt|=Thread.interrupted();
			for(final Thread worker:threads){
				if(worker.getState()==Thread.State.NEW){continue;}
				while(worker.isAlive()){
					try{worker.join();}catch(InterruptedException interrupted){
						restoreInterrupt=true;
						interrupted.printStackTrace();
					}
				}
			}
		}
	}

	/** Called after I/O resource owners close; workers have already joined in startAndJoin. */
	@Override
	public void close(){
		assert(Thread.currentThread()==owner) : "Only the creating thread can restore its pending interruption: "+owner.getName();
		if(restoreInterrupt){owner.interrupt();}
	}

	private final Thread owner=Thread.currentThread();
	private boolean restoreInterrupt;
}
