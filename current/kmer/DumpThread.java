package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date Nov 16, 2015
 *
 */
public class DumpThread extends Thread{
	
	/**
	 * Launches multiple threads to dump k-mers from k-mer tables to a stream writer.
	 * Creates optimal number of threads based on table count and system resources.
	 * Each thread processes different tables using atomic work distribution.
	 *
	 * @param k K-mer length for output formatting
	 * @param mincount Minimum k-mer count to include in output
	 * @param maxcount Maximum k-mer count to include in output
	 * @param tables Array of k-mer tables to process
	 * @param bsw ByteStreamWriter for output serialization
	 * @param remaining AtomicLong counter tracking remaining k-mers to dump
	 * @return true if all threads completed successfully, false otherwise
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
	 * Creates a new DumpThread with parameters for k-mer extraction.
	 *
	 * @param k_ K-mer length for output formatting
	 * @param mincount_ Minimum k-mer count threshold
	 * @param maxcount_ Maximum k-mer count threshold
	 * @param nextTable_ Atomic counter for work distribution among threads
	 * @param tables_ Array of k-mer tables to process
	 * @param bsw_ ByteStreamWriter for serialized output
	 * @param toDump_ AtomicLong tracking remaining k-mers to process
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
	
	/** K-mer length for output formatting */
	final int k;
	/** Minimum k-mer count threshold for inclusion in output */
	final int mincount;
	/** Maximum k-mer count threshold for inclusion in output */
	final int maxcount;
	/**
	 * Atomic counter tracking remaining k-mers to be processed across all threads
	 */
	final AtomicLong remaining;
	/** Atomic index for work distribution - indicates next table to process */
	final AtomicInteger nextTable;
	/** Array of k-mer tables to process for k-mer extraction */
	final AbstractKmerTable[] tables;
	/** ByteStreamWriter for serializing k-mer output to stream */
	final ByteStreamWriter bsw;
	/** Flag indicating whether this thread completed processing successfully */
	boolean success=false;
	
	/** Override for thread count; if negative, uses automatic thread calculation */
	public static int NUM_THREADS=-1;
	
}
