package prot;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

/** Reference-owned bounded worker pool for independent stop-tail batches.
 * Callers retain immutable peptide bytes until assign returns. Null query slots
 * denote preclassified unsupported stops and remain null. Each result preserves
 * its input slot; worker completion order never controls gene-candidate order.
 * @author Keqing
 */
public final class ProteinTailPool implements AutoCloseable{
	public ProteinTailPool(final ProteinTailAssigner reference, final int threads_){
		if(reference==null || threads_<1){throw new IllegalArgumentException("Tail pool requires a reference and positive worker count");}
		threads=threads_;
		workers=ThreadLocal.withInitial(reference::newWorker);
		final AtomicInteger names=new AtomicInteger();
		pool=threads==1 ? null : Executors.newFixedThreadPool(threads, task->{
			final Thread worker=new Thread(task, "CallGenes-HBM-"+names.incrementAndGet());
			worker.setDaemon(true); return worker;
		});
	}

	/** At most threads coarse tasks per call, never one Future or lock per protein.
	 * Concurrent callers share the same pool. Every task completes before return,
	 * including before any failed task is rethrown to the owning caller.
	 */
	public ProteinTailAssigner.Hit[] assign(final byte[][] queries){
		if(queries==null){throw new IllegalArgumentException("Null stop-tail batch");}
		final ProteinTailAssigner.Hit[] results=new ProteinTailAssigner.Hit[queries.length];
		if(queries.length==0){return results;}
		final int chunks=1+(queries.length-1)/BATCH;
		final AtomicInteger next=new AtomicInteger();
		final Callable<Void> task=()->{
			final ProteinTailAssigner.Worker worker=workers.get();
			for(int chunk=next.getAndIncrement(); chunk<chunks; chunk=next.getAndIncrement()){
				final int first=chunk*BATCH, end=first+Math.min(BATCH, queries.length-first);
				for(int i=first; i<end; i++){
					if(queries[i]!=null){results[i]=worker.assign("stop", queries[i]);}
				}
			}
			return null;
		};
		try{
			if(pool==null){task.call();}
			else{
				final ArrayList<Callable<Void>> tasks=new ArrayList<Callable<Void>>();
				for(int t=0; t<Math.min(threads, chunks); t++){tasks.add(task);}
				for(Future<Void> result : pool.invokeAll(tasks)){result.get();}
			}
		}catch(InterruptedException failure){
			Thread.currentThread().interrupt();
			throw new IllegalStateException("Interrupted while waiting for stop-tail workers", failure);
		}catch(ExecutionException failure){
			throw new IllegalStateException("A stop-tail worker failed", failure.getCause());
		}catch(RuntimeException | Error failure){throw failure;}
		catch(Exception failure){throw new IllegalStateException("Stop-tail batch failed", failure);}
		for(int i=0; i<queries.length; i++){
			assert((queries[i]==null)==(results[i]==null)) : "Every searchable tail must publish one result at its original index: "+i;
		}
		return results;
	}

	/** Normal callers have already joined every batch; daemon workers cannot hold a failed CLI open. */
	@Override
	public void close(){if(pool!=null){pool.shutdown();}}
	private final int threads;
	private final ExecutorService pool;
	private final ThreadLocal<ProteinTailAssigner.Worker> workers;
	private static final int BATCH=32;
}
