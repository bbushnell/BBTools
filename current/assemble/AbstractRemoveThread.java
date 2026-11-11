package assemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import kmer.AbstractKmerTableSet;
import kmer.HashArray1D;
import kmer.KmerNode;
import kmer.KmerTableSet;
import shared.Timer;
import ukmer.HashArrayU1D;
import ukmer.KmerNodeU;
import ukmer.KmerTableSetU;

/**
 * Removes k-mers with counts outside a retained inclusive range.
 * Shared worker scaffold for standard and unlimited k-mer tables.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
public abstract class AbstractRemoveThread extends Thread{

	/**
	 * Creates a worker sharing the caller's table counter.
	 *
	 * @param id_ worker identifier
	 * @param min_ minimum count to retain
	 * @param max_ maximum count to retain
	 * @param nextTable_ shared next-table counter, initially zero
	 */
	public AbstractRemoveThread(int id_, int min_, int max_, AtomicInteger nextTable_){
		id=id_;
		min=min_;
		max=max_;
		nextTable=nextTable_;
		assert(nextTable.get()==0);
	}
	
	/**
	 * Processes tables until the shared counter reaches the table count.
	 * Retains any worker failure for the coordinating process() call.
	 */
	@Override
	public final void run(){
		try{
			while(processNextTable()){}
		}catch(Throwable t){
			failure=t;
		}
	}
	
	/**
	 * Claims and prunes one table when any remain.
	 * Implementations exclusively mutate count storage in their claimed shard,
	 * clear ownership state, and regenerate table metadata after zeroing all
	 * values outside the retained count range.
	 *
	 * @return true if a table was processed; false when no tables remain
	 */
	abstract boolean processNextTable();
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Runs k-mer pruning workers across a standard or unlimited table set.
	 * Uses {@link KmerTableSet} for that exact class and casts all other inputs to
	 * {@link KmerTableSetU}; callers must
	 * provide a positive worker count and must not read or mutate the table set
	 * concurrently while workers clear ownership and regenerate processed shards.
	 * Interrupted joins are logged and retried until each worker terminates.
	 *
	 * @param threads positive number of worker threads
	 * @param min minimum count to retain
	 * @param max maximum count to retain
	 * @param tables table set to process
	 * @param print whether to emit timing and removal statistics
	 * @return total k-mers removed across all tables
	 * @throws RuntimeException first worker runtime failure, after all workers join
	 * @throws Error first worker error, after all workers join
	 */
	public static long process(final int threads, final int min, final int max, AbstractKmerTableSet tables, boolean print){
		Timer t=new Timer();
		
		final AtomicInteger nextTable=new AtomicInteger(0);
		long kmersRemoved=0;
		Throwable failure=null;
		
		/* Create Removethreads */
		ArrayList<AbstractRemoveThread> alpt=new ArrayList<AbstractRemoveThread>(threads);
		for(int i=0; i<threads; i++){
			final AbstractRemoveThread art;
			if(tables.getClass()==KmerTableSet.class){
				art=new RemoveThread1(i, min, max, nextTable, (KmerTableSet)tables);
			}else{
				art=new RemoveThread2(i, min, max, nextTable, (KmerTableSetU)tables);
			}
			alpt.add(art);
		}
		for(AbstractRemoveThread pt : alpt){pt.start();}
		
		for(AbstractRemoveThread pt : alpt){
			while(true){
				try{
					pt.join();
					break;
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			
			kmersRemoved+=pt.kmersRemovedT;
			if(failure==null && pt.failure!=null){failure=pt.failure;}
		}
		rethrow(failure);

		t.stop();
		if(print){
			outstream.println("Removed "+kmersRemoved+" kmers.");
			outstream.println("Remove time: "+t);
		}
		
		return kmersRemoved;
	}

	/** Throws a worker failure after all workers have joined. */
	private static void rethrow(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new RuntimeException(failure);
	}
	
	/*--------------------------------------------------------------*/
	
	private static class RemoveThread1 extends AbstractRemoveThread{

		/** Creates a standard-table removal worker. */
		public RemoveThread1(int id_, int min_, int max_, AtomicInteger nextTable_, KmerTableSet tables_){
			super(id_, min_, max_, nextTable_);
			tables=tables_;
		}
		
		@Override
		boolean processNextTable(){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArray1D table=tables.getTable(tnum);
			final int[] values=table.values();
			final int lim=table.arrayLength();
			for(int cell=0; cell<lim; cell++){
				final int value=values[cell];
				if(value<min || value>max){values[cell]=0;}
			}
			for(KmerNode kn : table.victims().array()){
				if(kn!=null){traverseKmerNode(kn);}
			}
			
			table.clearOwnership();
			kmersRemovedT+=table.regenerate(0);
			return true;
		}
		
		/** Clears out-of-range counts from a binary-tree victim chain. */
		private void traverseKmerNode(KmerNode kn){
			if(kn==null){return;}
			final int value=kn.count();
			if(value<min || value>max){kn.set(0);}
			traverseKmerNode(kn.left());
			traverseKmerNode(kn.right());
		}
		
		/** Concrete standard table set mutated by this worker. */
		private final KmerTableSet tables;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private static class RemoveThread2 extends AbstractRemoveThread{

		/** Creates an unlimited-table removal worker. */
		public RemoveThread2(int id_, int min_, int max_, AtomicInteger nextTable_, KmerTableSetU tables_){
			super(id_, min_, max_, nextTable_);
			tables=tables_;
		}
		
		@Override
		boolean processNextTable(){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final int[] values=table.values();
			final int lim=table.arrayLength();
			for(int cell=0; cell<lim; cell++){
				final int value=values[cell];
				if(value<min || value>max){values[cell]=0;}
			}
			for(KmerNodeU kn : table.victims().array()){
				if(kn!=null){traverseKmerNode(kn);}
			}
			
			table.clearOwnership();
			kmersRemovedT+=table.regenerate(0);
			return true;
		}
		
		/** Clears out-of-range counts from a binary-tree victim chain. */
		private void traverseKmerNode(KmerNodeU kn){
			if(kn==null){return;}
			final int value=kn.count();
			if(value<min || value>max){kn.set(0);}
			traverseKmerNode(kn.left());
			traverseKmerNode(kn.right());
		}
		
		/** Concrete unlimited table set mutated by this worker. */
		private final KmerTableSetU tables;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Number of k-mers removed by this worker. */
	long kmersRemovedT=0;

	/** Throwable raised by this worker, visible after join. */
	Throwable failure=null;
	
	/** Thread identifier retained by concrete workers. */
	final int id;
	/** Minimum k-mer count threshold for retention. */
	final int min;
	/** Maximum k-mer count threshold for retention. */
	final int max;
	
	/** Shared atomic counter indicating the next table index to process. */
	final AtomicInteger nextTable;
	
	/** Print stream used for status messages and timing output. */
	static PrintStream outstream=System.err;
	
}
