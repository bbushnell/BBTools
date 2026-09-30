package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Queued FASTA and numeric QUAL writer with one formatting/output thread.
 * Batches use the ordered JobQueue contract, with dense IDs beginning at zero.
 * Submitted lists, reads and their payloads are borrowed until consumed.
 * Each batch writes FASTA before QUAL; the two writes are not atomic.
 * Output selection follows Read.pairnum and mate links; entries are not deduplicated.
 * Lifecycle flags and counters alone do not certify arbitrary concurrent use.
 *
 * @author Collei, Brian Bushnell
 * @date November 22, 2025
 */
public class FastaQualWriterST implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens both outputs and creates the queue; processing starts separately.
	 * QUAL inherits the FASTA append choice; its descriptor enables subprocesses
	 * and requests overwrite separately. At least one mate selection is required.
	 * @param ffFa Nonnull FASTA output descriptor
	 * @param qf Nonnull QUAL output filename
	 * @param writeR1_ Emit entries whose pair number is zero
	 * @param writeR2_ Emit pair-number-one entries or the mate of other entries
	 */
	public FastaQualWriterST(FileFormat ffFa, String qf,
			boolean writeR1_, boolean writeR2_){
		ffoutFa=ffFa;
		//Qual file must mirror the fasta's append mode or the two files desynchronize on app=t
		ffoutQual=FileFormat.testOutput(qf, FileFormat.QUAL, null, true, true, ffFa.append(), false);
		fnameFa=ffFa.name();
		fnameQual=qf;

		writeR1=writeR1_;
		writeR2=writeR2_;

		assert(writeR1 || writeR2) : "Must write at least one mate";

		// Always configure ordered batches for this writer.
		ordered=true;
		queueCapacity=4;

		// JobQueue handles ordering
		queue=new JobQueue<ListNum<Read>>(queueCapacity, ordered, true, 0);
		queue.name="*FastaQualWriterST";

		// Open output streams
		//ff.append() must be honored: app=t previously truncated (hardcoded false; replicated via stream.sh 2026-09-05)
		outstreamFa=ReadWrite.getOutputStream(fnameFa, ffFa.append(), true, ffFa.allowSubprocess());
		outstreamQual=ReadWrite.getOutputStream(fnameQual, ffoutQual.append(), true, ffoutQual.allowSubprocess());

		if(verbose){outstream.println("Made FastaQualWriterST for "+fnameFa);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the writer thread once; repeated calls do not create a new thread. */
	@Override
	public void start(){
		synchronized(this){
			if(started){return;}
			started=true;
			writerThread=new Thread(new WriterRunnable());
			writerThread.start();
			if(verbose){outstream.println("Started "+getClass().getName());}
		}
	}

	/** Wraps a borrowed nonnull list and submits it with its ordered batch ID.
	 * @param list Reads retained until consumed, with null entries skipped
	 * @param id Dense batch ID beginning at zero */
	@Override
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}

	/**
	 * Submits a borrowed batch, starting the writer if needed; null batches are ignored.
	 * Queue backpressure can block submission. Do not modify submitted data before
	 * consumption, or submit ordinary data after poisoning.
	 * @param reads Batch with an ordered ID and nonnull list
	 */
	@Override
	public void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		if(!started){start();}

		// Add to queue (blocks if full/backpressure enabled)
		queue.add(reads);
	}

	/**
	 * Converts SAM entries to Reads and submits them under the original batch ID.
	 * Names, bases and qualities are retained without copying their payload arrays;
	 * numeric IDs become -1. SAM pair numbers are retained; mate links are not constructed.
	 * @param lines Batch of nonnull SAM entries; a null batch is ignored
	 */
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines){
			// Fixed #002: retain the SAM mate number so writeReads applies R1/R2 selection.
			//The constructor alone leaves both mates as pairnum0 with no mate link.
			Read r=new Read(sl.seq, sl.qual, sl.qname, -1, false);
			r.setPairnum(sl.pairnum());
			reads.add(r);
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}

	/** Queues one end marker after the largest observed ID; does not start or join the writer. */
	@Override
	public synchronized void poison(){
		if(poisoned){return;}
		poisoned=true;

		// Calculate poison ID
		long poisonID=queue.maxSeen()+1;
		ListNum<Read> poison=new ListNum<Read>(null, poisonID, ListNum.POISON);
		queue.poison(poison, false);

		if(verbose){outstream.println("Poisoned "+getClass().getName());}
	}

	/**
	 * Sends poison if needed, attempts to join an existing thread, then finishes both outputs.
	 * Does not start a thread. An interrupted join is ignored without a retry here.
	 * Marks closed after finishing streams; an already-closed call returns the stored flag.
	 * @return Accumulated writer error and output-finishing status
	 */
	@Override
	public synchronized boolean waitForFinish(){
		if(closed){return errorState;}

		// Ensure poison is sent if not already
		if(!poisoned){poison();}

		// Wait for worker thread to finish draining
		if(writerThread!=null){
			try{
				writerThread.join();
			}catch(InterruptedException e){
				// Ignore
			}
		}

		boolean b=ReadWrite.finishWriting(null, outstreamFa, fnameFa, ffoutFa.allowSubprocess());
		boolean b2=ReadWrite.finishWriting(null, outstreamQual, fnameQual, ffoutQual.allowSubprocess());

		outstreamFa=null;
		outstreamQual=null;
		closed=true;

		return errorState|=(b || b2);
	}

	/** Sends the end marker and delegates to waitForFinish.
	 * @return Accumulated error status */
	@Override
	public synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/** Marks an error, force-poisons the queue and marks closed without joining or closing streams. */
	@Override
	public synchronized void finishError(){
		errorState=true;
		poisoned=true;
		queue.poison(new ListNum<Read>(null, queue.maxSeen()+1, ListNum.POISON), true);
		closed=true;
	}

	/** ORs a status into the error flag, printing a stack trace on the first true transition. */
	private boolean setError(boolean b){
		if(b && !errorState){
			new RuntimeException("Triggered error state. this="+toString()).printStackTrace(outstream);
		}
		errorState|=b;
		return errorState;
	}

	/** Returns the observed error flag; this call does not wait for completion. */
	@Override
	public boolean errorState(){return errorState;}

	/** Returns the observed closed-without-error flags, without joining the writer. */
	@Override
	public boolean finishedSuccessfully(){return !errorState && closed;}

	/** Returns both output names in a parenthesized, comma-separated string. */
	@Override
	public final String fname(){return "("+fnameFa+","+fnameQual+")";}

	/** Returns output names and the observed closed/poisoned flags. */
	@Override
	public final String toString(){
		return "FastaQualWriterST"+fname()+" closed="+closed+", poisoned="+poisoned;
	}

	/** Returns the observed count of formatted individual records, incremented before batch writes. */
	@Override
	public long readsWritten(){return readsWritten;}

	/** Returns the observed formatted-base count; it is not an output-flush receipt. */
	@Override
	public long basesWritten(){return basesWritten;}

	/*--------------------------------------------------------------*/
	/*----------------        Writer Thread         ----------------*/
	/*--------------------------------------------------------------*/

	/** Consumes ordered batches and formats them on the single writer thread. */
	private class WriterRunnable implements Runnable{
		/** Drains until an end marker or null; retains the historical writer-error handling below. */
		@Override
		public void run(){
			try{
				// JobQueue handles ordering. We just take() the next valid job.
				for(ListNum<Read> job=queue.take(); job!=null && !job.poison(); job=queue.take()){
					writeReads(job.list);
				}
			}catch(Throwable t){
				//[stream/FastaQualWriterST#001] A write failure (write() rethrows IOException) must not strand
				//the producer. Set errorState (so waitForFinish reports loud), then FORCE-poison the queue
				//(force=true sets JobQueue.poisoned so a producer blocked in add()'s capacity-wait unblocks
				//instead of hanging), then surface loud. The blessed SWEEP-A writer-death fix.
				errorState=true;
				if(queue!=null){
					try{queue.poison(new ListNum<Read>(null, queue.maxSeen()+1, ListNum.POISON), true);}catch(Throwable t2){}
				}
				throw new RuntimeException("FastaQualWriterST writer thread failed; output may be incomplete.", t);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Logic          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Formats selected records into two fresh batch buffers, then writes FASTA followed by QUAL.
	 * Null entries are skipped. R1 selects pairnum0; R2 selects pairnum1 directly or
	 * a different entry's mate. Repeated records or links can therefore be emitted repeatedly.
	 */
	private void writeReads(ArrayList<Read> reads){
		ByteBuilder bbFa=new ByteBuilder();
		ByteBuilder bbQual=new ByteBuilder();

		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);

			if(writeR1 && r1!=null){
				writeRead(r1, bbFa, bbQual);
			}
			if(writeR2 && r2!=null){
				writeRead(r2, bbFa, bbQual);
			}
		}

		// Write both buffers to their respective streams
		// Single thread, so no race condition between Fa and Qual writes here
		write(bbFa, outstreamFa);
		write(bbQual, outstreamQual);
	}

	/**
	 * Appends FASTA using current Shared.FASTA_WRAP and one unwrapped numeric QUAL line.
	 * Existing signed quality bytes are printed as integers; null qualities use current
	 * Shared.FAKE_QUAL once per base. Both headers use r.id, falling back to numericID
	 * for a null name. Increments counters before either buffer is written.
	 */
	private void writeRead(Read r, ByteBuilder bbFa, ByteBuilder bbQual){
		// Write Fasta
		r.toFasta(bbFa);
		bbFa.nl();

		// Fixed #003: use Read.toFasta's numeric-ID fallback so paired output headers agree.
		//Appending a null String directly instead emitted the literal null in QUAL.
		// Write Qual header
		bbQual.append('>');
		if(r.id==null){bbQual.append(r.numericID);}else{bbQual.append(r.id);}
		bbQual.nl();

		// Write Qual scores
		byte[] quals=r.quality;
		if(quals!=null){
			for(byte b : quals){
				int q=b;
				bbQual.append(q).append(' ');
			}
			if(quals.length>0){bbQual.length--;} //Trim trailing space
		}else{
			int fake=Shared.FAKE_QUAL;
			int len=r.length();
			for(int i=0; i<len; i++){
				bbQual.append(fake).append(' ');
			}
			if(len>0){bbQual.length--;}
		}
		bbQual.nl();

		readsWritten++;
		basesWritten+=r.length();
	}

	/** Writes a nonempty buffer, clearing it after success; wraps IOException in RuntimeException. */
	private void write(ByteBuilder bb, OutputStream os){
		if(bb.length()==0){return;}//was <0 (never fired); skip empty buffers
		byte[] array=bb.toBytes();
		try{
			os.write(array); //No synchronized(this) needed; only writerThread calls this
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** FASTA output name retained from the descriptor. */
	private final String fnameFa;
	/** QUAL output name supplied by the caller. */
	private final String fnameQual;
	/** Borrowed FASTA output descriptor. */
	private final FileFormat ffoutFa;
	/** QUAL descriptor with matching append mode and separately configured options. */
	private final FileFormat ffoutQual;

	/** FASTA stream opened at construction, nulled by normal finishing. */
	private OutputStream outstreamFa;
	/** QUAL stream opened at construction, nulled by normal finishing. */
	private OutputStream outstreamQual;

	/** Whether to format pair-number-zero entries. */
	private final boolean writeR1;
	/** Whether to format second mates selected from entries or mate links. */
	private final boolean writeR2;

	/** Always true for this writer’s queue configuration. */
	private final boolean ordered;
	/** Configured queue capacity, fixed at four. */
	private final int queueCapacity;

	/** Ordered, backpressured queue configured to begin at batch ID zero. */
	private final JobQueue<ListNum<Read>> queue;
	/** Single formatting/output thread created by start. */
	private Thread writerThread;

	/** Formatted individual records; updated before batch output. */
	private long readsWritten=0;
	/** Bases in formatted records; updated before batch output. */
	private long basesWritten=0;

	/** Observed error flag, also publicly mutable for legacy callers. */
	public boolean errorState=false;
	/** Whether start has created the writer thread. */
	private boolean started=false;
	/** Whether an end marker has been requested. */
	private boolean poisoned=false;
	/** Set by normal finishing or finishError, not itself proof of stream closure. */
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Compile-time diagnostic logging switch. */
	public static final boolean verbose=false;

	/** Print status messages to this output stream */
	private final PrintStream outstream=System.err;

}
