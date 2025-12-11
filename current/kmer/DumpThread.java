package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Multithreaded worker for dumping k-mers from k-mer tables to output streams.
 * Each thread processes a subset of k-mer tables, extracting k-mers within specified count ranges.
 * Coordinates with other DumpThread instances via atomic counters for work distribution.
 *
 * @author Brian Bushnell
 * @date Nov 16, 2015
 */
public class DumpThread extends Thread{
	
	/**
	 * Launches multiple DumpThread instances to process k-mer tables using a shared writer and atomic counters.
	 * Determines thread count, distributes table work, waits for completion, and returns success.
	 * @param k K-mer length for output
	 * @param mincount Minimum count threshold
	 * @param maxcount Maximum count threshold
	 * @param tables K-mer tables to process
	 * @param bsw Output writer
	 * @param remaining Atomic counter tracking remaining k-mers
	 * @return true if all threads completed successfully
	 */
	public static boolean dump(final int k, final int mincount, final int maxcount, final AbstractKmerTable[] tables, final ByteStreamWriter bsw, AtomicLong remaining){
		final int threads=NUM_THREADS>0 ? NUM_THREADS : Tools.min(tables.length, (Tools.mid(1, Shared.threads()-1, 6)));
		final AtomicInteger lock=new AtomicInteger(0);
		final ArrayList<DumpThread> list=new ArrayList<DumpThread>(threads);
		for(int i=0; i<threads; i++){
			list.add(new DumpThread(k, mincount, maxcount, lock, tables, bsw, remaining));
		}
		for(DumpThread t : list){t.start();}
		boolean success=true;
		for(DumpThread t : list){
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
	 * Constructs a worker with thresholds, shared table index, table array, output writer, and remaining counter.
	 * @param k_ K-mer length
	 * @param mincount_ Minimum count
	 * @param maxcount_ Maximum count
	 * @param nextTable_ Atomic table index for work sharing
	 * @param tables_ Tables to dump
	 * @param bsw_ Output writer
	 * @param toDump_ Atomic remaining counter
	 */
	public DumpThread(final int k_, final int mincount_, final int maxcount_, final AtomicInteger nextTable_, final AbstractKmerTable[] tables_, final ByteStreamWriter bsw_, final AtomicLong toDump_){
		k=k_;
		mincount=mincount_;
		maxcount=maxcount_;
		nextTable=nextTable_;
		tables=tables_;
		bsw=bsw_;
		remaining=toDump_;
	}
	
	/**
	 * Processes k-mer tables by dumping k-mers within count range to output stream.
	 * Uses atomic work distribution to claim tables, then delegates to table-specific dump methods.
	 * Synchronizes final buffer writes to prevent output corruption.
	 */
	@Override
	public void run(){
		final ByteBuilder bb=new ByteBuilder(16300);
		for(int i=nextTable.getAndIncrement(); i<tables.length; i=nextTable.getAndIncrement()){
			AbstractKmerTable t=tables[i];
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
	final AtomicLong remaining;
	final AtomicInteger nextTable;
	final AbstractKmerTable[] tables;
	final ByteStreamWriter bsw;
	boolean success=false;
	
	public static int NUM_THREADS=-1;
	
}
