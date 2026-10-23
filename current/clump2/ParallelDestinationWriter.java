package clump2;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLongArray;

import stream.Read;
import stream.Writer;

/**
 * Reuses a bounded worker set to write one split KmerSort3 round across
 * zero-threaded destination writers. Each round gives every destination at most
 * one add opportunity and waits for every worker to observe the round before
 * the submitting coordinator returns.
 */
final class ParallelDestinationWriter{
	
	/**
	 * Starts a bounded pool for the supplied destination array.
	 * @param destinations_ Output writers indexed like submitted rounds
	 * @param ids_ Dense per-destination output IDs
	 * @param maxWorkers Maximum worker threads to keep for the helper lifetime
	 */
	ParallelDestinationWriter(final Writer[] destinations_, final AtomicLongArray ids_, final int maxWorkers){
		if(destinations_==null || destinations_.length<1){
			throw new IllegalArgumentException("Parallel writer requires at least one destination");
		}
		if(ids_==null || ids_.length()<destinations_.length){
			throw new IllegalArgumentException("Need one ID counter per destination");
		}
		final int workerCount=Math.max(1, Math.min(destinations_.length, maxWorkers));
		destinations=destinations_;
		ids=ids_;
		workers=new Worker[workerCount];
		for(int i=0; i<workers.length; i++){
			workers[i]=new Worker(i);
		}
		
		int started=0;
		try{
			for(int i=0; i<workers.length; i++){
				workers[i].start();
				started++;
			}
		}catch(RuntimeException|Error e){
			stopStartedWorkers(started);
			throw e;
		}
	}
	
	/**
	 * Writes one split round, returning only after every worker has observed this
	 * round or after a worker failure has been observed. A destination writer may
	 * queue a list behind a missing earlier ID, so the caller must preserve the
	 * list and Read payloads until the destination writers finish.
	 * @param round One list per destination; null and empty lists are skipped
	 */
	void submitRound(final ArrayList<Read>[] round){
		if(round==null || round.length!=destinations.length){
			throw new IllegalArgumentException("Round length must match destinations: "
					+(round==null ? -1 : round.length)+", "+destinations.length);
		}
		synchronized(lock){
			rethrow(failure);
			if(closed){throw new IllegalStateException("Parallel destination writer is closed");}
			if(activeRound!=null){throw new IllegalStateException("Concurrent submitRound calls are not supported");}
			activeRound=round;
			remaining=workers.length;
			nextDestination.set(0);
			roundGeneration++;
			lock.notifyAll();
			
			while(activeRound!=null){
				try{
					lock.wait();
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
					recordFailure(e);
					throw new RuntimeException(e);
				}
			}
			rethrow(failure);
		}
	}
	
	/** Stops the worker pool after any current round has completed. */
	void close(){
		synchronized(lock){
			while(activeRound!=null){
				try{
					lock.wait();
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
					recordFailure(e);
					break;
				}
			}
			closed=true;
			lock.notifyAll();
		}
		joinWorkers();
	}
	
	/** Retains the first worker failure and prevents further per-destination claims. */
	private void recordFailure(final Throwable cause){
		if(failure==null){failure=(cause==null ? new RuntimeException("Aborted") : cause);}
		nextDestination.set(destinations.length);
	}
	
	/** Marks one worker's round complete and wakes the coordinator after the last one. */
	private void finishWorkerRound(){
		synchronized(lock){
			remaining--;
			if(remaining==0){
				activeRound=null;
				lock.notifyAll();
			}
		}
	}
	
	/** Joins all workers, preserving interruption for the caller. */
	private void joinWorkers(){
		joinWorkers(workers.length);
	}
	
	/** Joins the requested live prefix, preserving interruption for the caller. */
	private void joinWorkers(final int count){
		boolean interrupted=false;
		for(int i=0; i<count; i++){
			final Worker worker=workers[i];
			while(worker.getState()!=Thread.State.TERMINATED){
				try{
					worker.join();
				}catch(InterruptedException e){
					interrupted=true;
				}
			}
		}
		if(interrupted){Thread.currentThread().interrupt();}
	}
	
	/** Stops already-started workers after a constructor start failure. */
	private void stopStartedWorkers(final int started){
		synchronized(lock){
			closed=true;
			lock.notifyAll();
		}
		joinWorkers(started);
	}
	
	/** Rethrows the original unchecked worker failure on the coordinator. */
	private static void rethrow(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof Error){throw (Error)failure;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		throw new RuntimeException(failure);
	}
	
	/** One reusable destination writer. */
	private final class Worker extends Thread{
		
		Worker(final int id){super("KmerSort3Writer-"+id);}
		
		@Override
		public void run(){
			int observedGeneration=0;
			while(true){
				final ArrayList<Read>[] round;
				synchronized(lock){
					while(activeRound==null || observedGeneration==roundGeneration){
						if(activeRound==null && (closed || failure!=null)){return;}
						try{
							lock.wait();
						}catch(InterruptedException e){
							Thread.currentThread().interrupt();
							if(activeRound==null || observedGeneration==roundGeneration){
								recordFailure(e);
								return;
							}
						}
					}
					if(activeRound==null){return;}
					round=activeRound;
					observedGeneration=roundGeneration;
				}
				
				try{
					for(int destination=nextDestination.getAndIncrement(); destination<destinations.length;
							destination=nextDestination.getAndIncrement()){
						final ArrayList<Read> list=round[destination];
						if(list!=null && !list.isEmpty()){
							destinations[destination].add(list, ids.getAndIncrement(destination));
						}
					}
				}catch(Throwable t){
					synchronized(lock){recordFailure(t);}
				}finally{
					finishWorkerRound();
				}
			}
		}
	}
	
	private final Writer[] destinations;
	private final AtomicLongArray ids;
	private final Worker[] workers;
	private final AtomicInteger nextDestination=new AtomicInteger(0);
	private final Object lock=new Object();
	
	private ArrayList<Read>[] activeRound;
	private int remaining=0;
	private int roundGeneration=0;
	private boolean closed=false;
	private Throwable failure=null;
}
