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
 * Removes k-mers whose counts fall outside a specified range.
 * Base class for threaded k-mer pruning across standard and unlimited table sets.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
public abstract class AbstractRemoveThread extends Thread{

	/**
	 * Creates a removal thread with bounds and shared table counter.
	 * Initializes ID, min/max thresholds, and the atomic table index tracker.
	 *
	 * @param id_ Thread identifier
	 * @param min_ Minimum k-mer count to retain
	 * @param max_ Maximum k-mer count to retain
	 * @param nextTable_ Shared atomic counter for table iteration
	 */
	public AbstractRemoveThread(int id_, int min_, int max_, AtomicInteger nextTable_){
		id=id_;
		min=min_;
		max=max_;
		nextTable=nextTable_;
		assert(nextTable.get()==0);
	}
	
	/** Executes removal across tables until none remain.
	 * Loops calling processNextTable() until exhaustion. */
	@Override
	public final void run(){
		while(processNextTable()){}
	}
	
	/**
	 * Processes the next table shard for k-mer removal.
	 * Implemented by subclasses for specific table types.
	 * @return true if a table was processed; false when no tables remain
	 */
	abstract boolean processNextTable();
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Runs multi-threaded k-mer pruning across a table set.
	 * Spawns worker threads (standard or unlimited), joins them, and aggregates removed counts.
	 *
	 * @param threads Number of worker threads
	 * @param min Minimum count to retain
	 * @param max Maximum count to retain
	 * @param tables Table set to process
	 * @param print Whether to emit timing/removal stats
	 * @return Total k-mers removed across all tables
	 */
	public static long process(final int threads, final int min, final int max, AbstractKmerTableSet tables, boolean print){
		Timer t=new Timer();
		
		final AtomicInteger nextTable=new AtomicInteger(0);
		long kmersRemoved=0;
		
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
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			kmersRemoved+=pt.kmersRemovedT;
		}

		t.stop();
		if(print){
			outstream.println("Removed "+kmersRemoved+" kmers.");
			outstream.println("Remove time: "+t);
		}
		
		return kmersRemoved;
	}
	
	/*--------------------------------------------------------------*/
	
	private static class RemoveThread1 extends AbstractRemoveThread{

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
		
		private void traverseKmerNode(KmerNode kn){
			if(kn==null){return;}
			final int value=kn.count();
			if(value<min || value>max){kn.set(0);}
			traverseKmerNode(kn.left());
			traverseKmerNode(kn.right());
		}
		
		private final KmerTableSet tables;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private static class RemoveThread2 extends AbstractRemoveThread{

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
		
		private void traverseKmerNode(KmerNodeU kn){
			if(kn==null){return;}
			final int value=kn.count();
			if(value<min || value>max){kn.set(0);}
			traverseKmerNode(kn.left());
			traverseKmerNode(kn.right());
		}
		
		private final KmerTableSetU tables;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Number of k-mers removed by this thread. */
	long kmersRemovedT=0;
	
	/** Thread identifier for coordination and logging. */
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