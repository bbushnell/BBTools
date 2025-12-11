package ukmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import kmer.DumpThread;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Thread class for parallel dumping of k-mers from AbstractKmerTableU arrays.
 * Each thread processes multiple tables sequentially while coordinating with other threads
 * to ensure all tables are processed exactly once.
 *
 * @author Brian Bushnell
 * @date Nov 16, 2015
 */
public class DumpThreadU extends Thread{
	
	/**
	 * Launches multiple DumpThreadU workers to dump k-mers with shared writer and atomic counters.
	 * @param k K-mer length
	 * @param mincount Minimum count threshold
	 * @param maxcount Maximum count threshold
	 * @param tables Tables to process
	 * @param bsw Output writer
	 * @param remaining Atomic remaining counter
	 * @return true if all threads succeeded
	 */
	public static boolean dump(final int k, final int mincount, final int maxcount, final AbstractKmerTableU[] tables, final ByteStreamWriter bsw, AtomicLong remaining){
		final int threads=DumpThread.NUM_THREADS>0 ? DumpThread.NUM_THREADS : Tools.min(tables.length, (Tools.mid(1, Shared.threads()-1, 6)));
		final AtomicInteger lock=new AtomicInteger(0);
		final ArrayList<DumpThreadU> list=new ArrayList<DumpThreadU>(threads);
		for(int i=0; i<threads; i++){
			list.add(new DumpThreadU(k, mincount, maxcount, lock, tables, bsw, remaining));
		}
		for(DumpThreadU t : list){t.start();}
		boolean success=true;
		for(DumpThreadU t : list){
			while(t.getState()!=Thread.State.TERMINATED){
				try {
					t.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			success&=t.success;
		}
		return success;
	}
	
	/**
	 * Worker constructor wiring k-mer thresholds, shared table index, table array, writer, and remaining counter.
	 * @param k_ K-mer length
	 * @param mincount_ Minimum count
	 * @param maxcount_ Maximum count
	 * @param nextTable_ Atomic table index for coordination
	 * @param tables_ Tables to dump
	 * @param bsw_ Output writer
	 * @param remaining_ Atomic remaining counter
	 */
	public DumpThreadU(final int k_, final int mincount_, final int maxcount_, final AtomicInteger nextTable_, final AbstractKmerTableU[] tables_, final ByteStreamWriter bsw_, AtomicLong remaining_){
		k=k_;
		mincount=mincount_;
		maxcount=maxcount_;
		nextTable=nextTable_;
		tables=tables_;
		bsw=bsw_;
		remaining=remaining_;
	}
	
	/**
	 * Thread execution method that dumps k-mers from assigned tables.
	 * Uses atomic counter to claim tables sequentially and processes each table
	 * by calling dumpKmersAsBytes_MT. Buffers output in ByteBuilder for efficiency.
	 */
	@Override
	public void run(){
		final ByteBuilder bb=new ByteBuilder(16300);
		for(int i=nextTable.getAndIncrement(); i<tables.length; i=nextTable.getAndIncrement()){
			AbstractKmerTableU t=tables[i];
			t.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
		}
		if(bb.length()>0){
			synchronized(bsw){bsw.addJob(bb);}
		}
		success=true;
	}
	
	final int k;
	final int mincount;
	final int maxcount;
	final AtomicInteger nextTable;
	final AtomicLong remaining;
	final AbstractKmerTableU[] tables;
	final ByteStreamWriter bsw;
	boolean success=false;
	
}
