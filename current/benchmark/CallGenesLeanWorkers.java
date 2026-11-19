package benchmark;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.IntFunction;

/** Ordered, bounded genome dispatch on a pool distinct from classifier workers.
 * Every active worker joins before this method returns or throws; cancellation
 * must never let the owner close shared biological resources under a caller.
 * Results accumulate for the entire cohort; only submitted pending tasks are bounded.
 * @author Sayu, Fischl
 */
public final class CallGenesLeanWorkers{

	/** Dispatch each ordinal in [0, count) and return nonnull results in input order.
	 * Submit at most twice the worker count ahead of the next input-order result.
	 * Each task owns its readers and mutable callers; callers retain shared resources
	 * until this method exits. Tasks check the first-failure flag before their body;
	 * already-started tasks finish and close before any return or exception.
	 * Interruption is restored after the pool joins; failed or interrupted runs throw
	 * instead of returning a partial result list.
	 * @param count Positive number of tasks
	 * @param workers Pool size, from one through count
	 * @param task Concurrent task function receiving a zero-based input ordinal
	 * @param <T> Result type
	 * @return New list of all successful results in ordinal order
	 */
	public static <T> ArrayList<T> run(int count, int workers, IntFunction<T> task){
		if(count<1 || workers<1 || workers>count || task==null){
			throw new IllegalArgumentException("Ordered dispatch requires a nonempty cohort and bounded positive workers");
		}
		final ExecutorService pool=Executors.newFixedThreadPool(workers);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ArrayDeque<Future<T>> pending=new ArrayDeque<Future<T>>();
		final ArrayList<T> result=new ArrayList<T>(count);
		final int window=(int)Math.min(count, 2L*workers);
		int submitted=0;
		try{
			while(result.size()<count){
				while(submitted<count && pending.size()<window && failure.get()==null){
					final int ordinal=submitted++;
					pending.add(pool.submit(()->{
						if(failure.get()!=null){return null;}
						try{
							final T value=task.apply(ordinal);
							if(value==null){
								throw new IllegalStateException("Missing genome result at ordinal "+ordinal);
							}
							return value;
						}catch(RuntimeException | Error problem){
							failure.compareAndSet(null, problem);
							throw problem;
						}
					}));
				}
				if(failure.get()!=null){throw new IllegalStateException("Genome worker failed", failure.get());}
				final T value=pending.remove().get();
				if(value==null){throw new IllegalStateException("Genome cohort is incomplete", failure.get());}
				result.add(value);
			}
		}catch(InterruptedException interrupted){
			failure.compareAndSet(null, interrupted);
			Thread.currentThread().interrupt();
			throw new IllegalStateException("Interrupted while joining genome workers", interrupted);
		}catch(ExecutionException failed){
			throw new IllegalStateException("Genome worker failed", failed.getCause());
		}finally{
			pool.shutdown();
			boolean interrupted=Thread.interrupted();
			while(!pool.isTerminated()){
				try{
					pool.awaitTermination(1, TimeUnit.DAYS);
				}catch(InterruptedException again){
					interrupted=true;
					failure.compareAndSet(null, again);
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
		}
		if(Thread.currentThread().isInterrupted() || failure.get()!=null){
			throw new IllegalStateException("Genome dispatch failed or was interrupted during final join", failure.get());
		}
		return result;
	}
	/** Static dispatcher; no instance state. */
	private CallGenesLeanWorkers(){}
}
