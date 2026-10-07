package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * FASTQ, FASTA, HEADER and SCARF writer driven by the caller or one worker thread.
 * Construction opens the output; start creates the optional worker and is also
 * called lazily on the first submission. A ByteBuilder holds one formatted batch.
 * Ordered or threaded configurations use JobQueue, which currently always requires
 * dense batch IDs starting at zero, including empty batches for skipped work.
 * Lists, wrappers and Read payloads remain borrowed until formatting completes;
 * queued caller-driven batches may still be retained when their submission returns.
 * Direct unordered, unthreaded use requires one submitting caller: the stream lock
 * protects byte writes, not conversion, counters or all lifecycle state.
 * Finish producers before normal waitForFinish or poisonAndWait completion and
 * inspect the returned error state. Counters are observations, not completion barriers.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 18, 2025
 */
public class FastqWriterST2 implements Writer{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Opens ordered caller-driven output with queue capacity three and no append.
	 * @param out_ Output name; its extension selects format, defaulting to FASTQ
	 * @param writeR1_ Select first mates
	 * @param writeR2_ Select second mates
	 * @param overwrite Permit replacing an existing file */
	public FastqWriterST2(String out_, boolean writeR1_, boolean writeR2_, boolean overwrite){
		this(FileFormat.testOutput(out_, FileFormat.FASTQ, null, true, overwrite, false, true), 
			writeR1_, writeR2_, false, 3);
	}
	
	/** Opens caller-driven output with ordering from the format and capacity three.
	 * @param ffout_ Output format and opening policy
	 * @param writeR1_ Select pair-number-zero roots
	 * @param writeR2_ Select second mates */
	public FastqWriterST2(FileFormat ffout_, boolean writeR1_, boolean writeR2_){
		this(ffout_, writeR1_, writeR2_, false, 3);
	}

	/**
	 * Opens output with explicit worker use; at least one mate selection is required.
	 * @param ffout_ Output file format/path
	 * @param writeR1_ Write read 1
	 * @param writeR2_ Write read 2
	 * @param threaded_ If true, a separate thread will handle writing (Producer-Consumer model)
	 * @param queueCapacity_ Soft queue capacity; must exceed one when queued
	 */
	public FastqWriterST2(FileFormat ffout_, boolean writeR1_, boolean writeR2_, 
			boolean threaded_, int queueCapacity_){
		this(ffout_, writeR1_, writeR2_, threaded_, queueCapacity_, false);
	}

	/** Opens output with optional compression on the leaf's formatting thread only.
	 * @param ffout_ FASTQ/FASTA/HEADER/SCARF output; unknown format defaults to FASTQ
	 * @param writeR1_ Select first mates
	 * @param writeR2_ Select second mates
	 * @param threaded_ Use one formatting/writing worker when true
	 * @param queueCapacity_ Soft queue capacity, greater than one when queued
	 * @param lightweight_ Use native lightweight backends without subprocesses */
	FastqWriterST2(FileFormat ffout_, boolean writeR1_, boolean writeR2_,
			boolean threaded_, int queueCapacity_, boolean lightweight_){
		lightweight=lightweight_;
		ffout=ffout_;
		fname=ffout_.name();
		writeR1=writeR1_;
		writeR2=writeR2_;
		format=(ffout.format()==UNKNOWN ? FASTQ : ffout.format());
		assert(format==FASTQ || format==FASTA || format==HEADER || format==FileFormat.SCARF) : ffout;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Config
		ordered=ffout.ordered();
		threaded=threaded_;
		
		// Create a queue when the format requests ordering or a worker was requested.
		if(ordered || threaded){
			// JobQueue handles ordering, capacity bounds, and backpressure.
			queue=new JobQueue<ListNum<Read>>(queueCapacity_, ordered, true, 0);
			queue.name="*FQWriterST2";
		}else{
			// Pure lightweight ST/unordered mode (needs external synchronization if MT)
			queue=null;
		}
		
		// Open output stream
		//ffout.append() must be honored: app=t previously truncated (hardcoded false; replicated via stream.sh 2026-09-05)
		outstream=(lightweight ? LightweightOutputStream.open(ffout)
			: ReadWrite.getOutputStream(fname, ffout.append(), true, false));
		if(verbose){outstream2.println("Made FastqWriterST (Ordered: "+ordered+", Threaded: "+threaded+")");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Starts the optional writer worker once; direct mode only records startup. */
	@Override
	public void start(){
		synchronized(this){
			if(started){return;}
			started=true;
			if(threaded && queue!=null){
				writerThread=new Thread(new WriterRunnable());
				writerThread.start();
			}
			if(verbose){outstream2.println("Started "+getClass().getName());}
		}
	}
	
	/** Ends submissions after producers finish and drains queued caller-driven work.
	 * Does not join a worker or close the output stream. */
	@Override
	public synchronized void poison(){
		if(verbose){outstream2.println("Called poison "+getClass().getName());}
		// Ensure only one thread executes the shutdown sequence
		if(poisoned){return;}
		poisoned=true;
		
		if(queue!=null){
			// 1. Calculate the ID for the poison pill. It must be > maxSeen to guarantee order.
			long poisonID=queue.maxSeen()+1;
			
			// 2. Create the poison pill (ListNum acts as the HasID object)
			ListNum<Read> poison=new ListNum<Read>(null, poisonID, ListNum.POISON);
			
			// 3. Inject the poison pill into the queue using the explicit JobQueue API.
			queue.poison(poison, false); 
			
			// 4. CRITICAL: If running in Host-Driven (unthreaded) mode, the producer thread 
			// calling poison() MUST perform the final drain.
			if(!threaded){
				// Use take() to block and ensure everything is drained until the pill is found.
				while(true){
					// take() will return null when it retrieves the poison pill.
					ListNum<Read> job=queue.take(); 
					
					if(job==null){break;}
					
					writeReads(job.list);
				}
			}
			// If threaded, the worker thread handles the take() loop and exits gracefully.
		}
		if(verbose){outstream2.println("Finished poison "+getClass().getName());}
	}
	
	/** Requests normal termination and finalizes output after waiting.
	 * @return Sticky reported error state */
	@Override
	public synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/** Marks abandonment without draining, joining the worker or closing the stream.
	 * Subsequent waitForFinish calls return the cached error without finalizing output. */
	@Override
	public synchronized void finishError(){
		errorState=true;
		if(queue!=null){
			poisoned=true;
			queue.poison(new ListNum<Read>(null, queue.maxSeen()+1, ListNum.POISON), true);
		}
		finished=true;
	}
	
	/** Wraps a borrowed nonnull Read list and submits it without copying.
	 * @param list Batch roots; null entries are skipped during formatting
	 * @param id Dense zero-based ID for queued configurations */
	@Override
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}
	
	/** Starts lazily and writes directly or retains the original wrapper in the queue.
	 * Caller-driven ordered operation drains ready batches after adding this one.
	 * Keep queued list contents and payloads stable until they have been formatted.
	 * @param reads Null wrapper is ignored; normal wrappers require a nonnull list */
	@Override
	public void addReads(ListNum<Read> reads){
		if(verbose){outstream2.println("Called addReads "+(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));}
		if(reads==null){return;}
		if(!started){start();}
		
		if(queue==null){
			// LIGHTWEIGHT MODE (Unordered/Unthreaded): Direct synchronized write.
			writeReads(reads.list);
		}else{
			if(verbose){outstream2.println("Adding to queue "+(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));}
			// 1. Add to queue (blocks if full/backpressure enabled)
			queue.add(reads);
			
			if(!threaded){
				// HOST-DRIVEN MODE (Ordered/Unthreaded): 
				// The producer acts as the draining consumer immediately after adding.
				// Uses the non-blocking poll() to drain any contiguous jobs that are ready.

				assert(ordered);
				synchronized(queue){//Essential to maintain ordering
					for(ListNum<Read> job=queue.poll(); job!=null && !job.poison(); job=queue.poll()){
						if(verbose){
							outstream2.println("Draining queue "+
								(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));
						}
						writeReads(job.list);
					}
				}
			}
			// If threaded, the background thread handles the draining.
		}
		if(verbose){outstream2.println("Finished addReads "+(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));}
	}
	
	/** Converts SAM sequence, quality and name fields into an unlinked Read batch.
	 * Payload arrays are shared and pair numbers retained for mate selection. Numeric
	 * IDs become -1; other alignment metadata and ListNum poison/last flags are not
	 * forwarded to the new wrapper. The batch ID is retained.
	 * @param lines Null wrapper is ignored; otherwise its list and records must be nonnull */
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines){
			//Fixed STR242: preserve SAM mate identity before writeReads selects R1 or R2.
			Read r=new Read(sl.seq, sl.qual, sl.qname, -1, false);
			r.setPairnum(sl.pairnum());
			reads.add(r);
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}
	
	/** @return Observed count of formatted selected reads, including HEADER records */
	@Override
	public long readsWritten(){return readsWritten;}
	
	/** @return Formatted sequence bases; HEADER output does not increment this count */
	@Override
	public long basesWritten(){return basesWritten;}
	
	/** Poisons, joins the optional worker and finalizes output after producers finish.
	 * An interrupted join restores interruption and still attempts finalization, so
	 * that path does not guarantee worker termination. Repeated calls return cached state.
	 * @return Sticky reported error state */
	@Override
	public synchronized boolean waitForFinish(){

		if(verbose){
			outstream2.println("FastqWriterST close() 1");
			new Exception().printStackTrace();
		}
		if(finished){return errorState;}
		poison(); 
		// poison() handles draining the queue if unthreaded.
		if(verbose){outstream2.println("FastqWriterST close() 2");}

		// If threaded, wait for the worker to process the poison pill and exit.
		if(threaded && writerThread!=null){
			try{
				writerThread.join();
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
			}
		}
		if(verbose){outstream2.println("FastqWriterST close() 3");}

		boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess() && !lightweight);
		if(verbose){outstream2.println("FastqWriterST close() 4");}
		finished=true;
		return errorState|=b;
	}
	
	/** @return Current error flag; this query is not a completion barrier */
	@Override
	public boolean errorState(){return errorState;}
	
	/** @return Whether finished is set and the cached error flag is clear */
	@Override
	public boolean finishedSuccessfully(){return !errorState && finished;}
	
	/** @return Output name captured at construction */
	@Override
	public final String fname(){return fname;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Logic          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Formats a whole batch in a fresh builder, then writes it under the stream lock.
	 * @param reads Borrowed nonnull roots; each formatter skips null entries */
	private void writeReads(ArrayList<Read> reads){
		ByteBuilder bb=new ByteBuilder();
		
		// Format reads
		if(format==FASTQ){
			writeFastq(reads, bb);
		}else if(format==FASTA){
			writeFasta(reads, bb);
		}else if(format==HEADER){
			writeHeader(reads, bb);
		}else if(format==FileFormat.SCARF){
			writeScarf(reads, bb);
		}else{
			throw new RuntimeException("Bad format: "+format);
		}

		// The actual write to stream must be synchronized
		write(bb);
		bb=null;
	}
	
	/** Writes nonempty bytes and clears the builder after a successful stream write.
	 * @param bb Buffer owned by the current formatting invocation */
	private void write(ByteBuilder bb){
		if(bb.length()==0){return;}//was <0 (never fired); skip empty buffers
		byte[] array=bb.toBytes();
		if(verbose){outstream2.println("FQWST write("+array.length+")");}
		try{
			synchronized(outstream){outstream.write(array);}//Works when synchronized on stream, hangs when synchronized on this
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
		if(verbose){outstream2.println("FQWST write("+array.length+") finished");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Writer Thread         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Optional single consumer that formats queued batches and writes their bytes. */
	private class WriterRunnable implements Runnable{
		/** Drains normal batches until the queue supplies its terminal signal. */
		@Override
		public void run(){
			if(verbose){outstream2.println("WR Thread started");}
			try{
				for(ListNum<Read> job=queue.take(); job!=null && !job.poison(); job=queue.take()){
					if(verbose){outstream2.println("WR Thread got "+job.id+", "+job.last());}
					// Blocking wait for next job (or poison pill)
					writeReads(job.list);
				}
			}catch(Throwable t){
				//[stream/FastqWriterST2#001] A write failure (write() rethrows IOException) must not strand the
				//producer. Set errorState (so waitForFinish reports loud), then FORCE-poison the queue
				//(force=true sets JobQueue.poisoned so a producer blocked in add()'s capacity-wait unblocks
				//instead of hanging), then surface loud. The blessed SWEEP-A writer-death fix.
				errorState=true;
				if(queue!=null){
					try{queue.poison(new ListNum<Read>(null, queue.maxSeen()+1, ListNum.POISON), true);}catch(Throwable t2){}
				}
				throw new RuntimeException("FastqWriterST2 writer thread failed; output may be incomplete.", t);
			}
			if(verbose){outstream2.println("WR Thread finished");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Appends selected mates as FASTQ and advances formatting counters.
	 * @param reads Roots in output order; null entries are skipped
	 * @param bb Destination byte builder */
	private void writeFastq(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				r1.toFastq(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r1.length();
			}
			if(writeR2 && r2!=null){
				r2.toFastq(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r2.length();
			}
		}
	}
	
	/** Appends selected mates as FASTA and advances formatting counters.
	 * @param reads Roots in output order; null entries are skipped
	 * @param bb Destination byte builder */
	private void writeFasta(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				r1.toFasta(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r1.length();
			}
			if(writeR2 && r2!=null){
				r2.toFasta(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r2.length();
			}
		}
	}
	
	/** Serializes selected mates with Read's existing SCARF formatter.
	 * @param reads Roots in output order; null entries are skipped
	 * @param bb Destination byte builder */
	private void writeScarf(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){r1.toScarf(bb).nl(); readsWritten++; basesWritten+=r1.length();}
			if(writeR2 && r2!=null){r2.toScarf(bb).nl(); readsWritten++; basesWritten+=r2.length();}
		}
	}

	/** Appends selected mate IDs and increments readsWritten, but not basesWritten.
	 * @param reads Roots in output order; null entries are skipped
	 * @param bb Destination byte builder */
	private void writeHeader(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				bb.appendln(r1.id);
				readsWritten++;
			}
			if(writeR2 && r2!=null){
				bb.appendln(r2.id);
				readsWritten++;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file path */
	public final String fname;
	/** Output file format */
	final FileFormat ffout;
	/** Output file format as an int */
	public final int format;
	/** Output opened by the constructor and finalized through ReadWrite. */
	OutputStream outstream;
	/** Write R1 reads (pairnum==0) */
	final boolean writeR1;
	/** Write R2 reads (pairnum==1 or mate) */
	final boolean writeR2;
	
	/** Format ordering flag; current JobQueue also orders threaded unordered requests. */
	final boolean ordered;
	/** Configuration: Use internal thread */
	final boolean threaded;
	/** Restricts output compression to the leaf's formatting thread. */
	private final boolean lightweight;
	
	/** Queue for ordering/threading. Null if unordered+unthreaded. */
	final JobQueue<ListNum<Read>> queue;
	/** Internal writer thread. Null if unthreaded. */
	Thread writerThread;
	
	/** Selected records counted during formatting, before output completion. */
	protected long readsWritten=0;
	/** Formatted sequence bases; stays zero for HEADER output. */
	protected long basesWritten=0;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** True after start() called */
	private boolean started=false;
	/** Normal finalization or explicit abandonment has been recorded. */
	private boolean finished=false;
	/** Normal or forced termination has been requested. */
	private boolean poisoned=false;

	/*--------------------------------------------------------------*/
	/*----------------       Output Settings        ----------------*/
	/*--------------------------------------------------------------*/

	/** Supported serialization identifiers and the fallback sentinel. */
	private static final int FASTQ=FileFormat.FASTQ;
	private static final int FASTA=FileFormat.FASTA;
	private static final int HEADER=FileFormat.HEADER;
	private static final int UNKNOWN=FileFormat.UNKNOWN;
	
	/** Compile-time diagnostic output switch. */
	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	protected PrintStream outstream2=System.err;
	
}
