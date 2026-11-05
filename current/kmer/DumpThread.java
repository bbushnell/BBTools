package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Partitions stable k-mer tables among workers that submit unordered text blocks.
 * Each worker owns its buffer; atomic counters distribute tables and an optional
 * record allowance. The caller owns the started writer and must finish it and
 * check its error state separately. Tables must not change during traversal.
 *
 * @author Brian Bushnell
 * @date Nov 16, 2015
 */
public class DumpThread extends Thread{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Starts and joins table workers; does not close or join the output writer.
	 * Output order depends on worker scheduling. Interrupted joins are logged and retried.
	 * @param k K-mer length for output
	 * @param mincount Inclusive minimum count forwarded to each table implementation
	 * @param maxcount Inclusive maximum count forwarded to each table implementation
	 * @param tables K-mer tables to process
	 * @param bsw Already-started writer, borrowed without transferring lifecycle ownership
	 * @param remaining Shared record allowance, or null for no limit; may become negative
	 * @return Whether every worker returned normally, not whether asynchronous output succeeded;
	 * delegated table return values are not accumulated
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
				try{
					t.join();
				}catch(InterruptedException e){
					//Preserve the existing policy of waiting for all workers despite interruption.
					e.printStackTrace();
				}
			}
			success&=t.success;
		}
		return success;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows tables, writer and shared counters; construction does not start the worker.
	 * @param k_ K-mer length
	 * @param mincount_ Minimum count
	 * @param maxcount_ Maximum count
	 * @param nextTable_ Atomic table index for work sharing
	 * @param tables_ Tables to dump
	 * @param bsw_ Output writer
	 * @param toDump_ Shared record allowance, or null for unlimited output
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

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Claims tables and delegates record formatting, then transfers its remaining buffer
	 * to the writer. An uncaught failure leaves success false; the coordinator joins
	 * this worker before reading that flag.
	 */
	@Override
	public void run(){
		//Table dump methods copy full blocks before reusing bb; the final block is transferred without reuse.
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

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Decoding width forwarded to table serializers. */
	final int k;
	/** Inclusive lower count bound forwarded to table implementations. */
	final int mincount;
	/** Inclusive upper count bound forwarded to table implementations. */
	final int maxcount;
	/** Shared allowance, or null for unlimited output. */
	final AtomicLong remaining;
	/** Next unclaimed table index; shared by this dump's workers. */
	final AtomicInteger nextTable;
	/** Borrowed tables, stable for the complete dump. */
	final AbstractKmerTable[] tables;
	/** Caller-owned started writer. */
	final ByteStreamWriter bsw;
	/** Published by normal completion and read after thread termination. */
	boolean success=false;

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Positive worker override for short and long dumps; otherwise use at most six and at most the table count. */
	public static int NUM_THREADS=-1;

}
