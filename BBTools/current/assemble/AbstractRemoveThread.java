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
 * Removes kmers with counts outside a certain range.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 */
public abstract class AbstractRemoveThread extends Thread{

	/**
	 * Constructor
	 */
	public AbstractRemoveThread(int id_, int min_, int max_, AtomicInteger nextTable_){
		id=id_;
		min=min_;
		max=max_;
		nextTable=nextTable_;
		assert(nextTable.get()==0);
	}
	
	@Override
	public final void run(){
		while(processNextTable()){}
	}
	
	abstract boolean processNextTable();
	
	/*--------------------------------------------------------------*/
	
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

		/**
		 * Constructor
		 */
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

		/**
		 * Constructor
		 */
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
	
	long kmersRemovedT=0;
	
	final int id;
	final int min;
	final int max;
	
	final AtomicInteger nextTable;
	
	/** Print messages to this stream */
	static PrintStream outstream=System.err;
	
}