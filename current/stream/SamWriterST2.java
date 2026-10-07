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
 * SAM/BAM writer that formats on the calling thread or one optional writer thread.
 * Construction opens the output; start writes its header and starts the worker if requested.
 * Ordered or threaded configurations use JobQueue, which currently always requires
 * dense batch IDs starting at zero, including empty batches for skipped work.
 * Direct unordered, unthreaded use requires one submitting caller: the stream lock
 * protects individual byte writes, not all conversion, counters or batch ordering.
 * Coordinate producers before normal poisonAndWait completion and inspect its result.
 * Submitted SAM lists and records remain borrowed until their writes complete.
 * Read conversion creates a new list but may reuse and normalize attached SamLines.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 18, 2025
 */
public class SamWriterST2 implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens caller-driven output with both mates selected and queue capacity three.
	 * Ordering comes from the output format.
	 * @param ffout_ Output format and opening policy
	 * @param header_ Optional borrowed SAM header lines
	 * @param useSharedHeader_ Prefer the shared input header over the supplied header */
	public SamWriterST2(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_){
		this(ffout_, header_, useSharedHeader_, false, 3, true, true);
	}
	
	/**
	 * Opens output with both mates selected and configurable worker use.
	 * @param ffout_ Output file format
	 * @param header_ Optional header lines
	 * @param useSharedHeader_ If true, pulls header from SamReadInputStream
	 * @param threaded_ If true, a separate thread handles writing (Producer-Consumer)
	 * @param queueCapacity_ Soft queue capacity; must exceed one when a queue is needed
	 */
	public SamWriterST2(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_,
			boolean threaded_, int queueCapacity_){
		this(ffout_, header_, useSharedHeader_, threaded_, queueCapacity_, true, true);
	}

	/** Opens output with explicit mate selection.
	 * @param ffout_ Output format, append policy and ordering flag
	 * @param header_ Optional borrowed SAM header lines
	 * @param useSharedHeader_ Prefer the shared input header
	 * @param threaded_ Format and write queued batches on one worker when true
	 * @param queueCapacity_ Soft queue capacity, greater than one when queued
	 * @param writeR1_ Select pair-number-zero roots from Read batches
	 * @param writeR2_ Select second mates from Read batches */
	public SamWriterST2(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_,
			boolean threaded_, int queueCapacity_, boolean writeR1_, boolean writeR2_){
		this(ffout_, header_, useSharedHeader_, threaded_, queueCapacity_, writeR1_, writeR2_, false);
	}

	/** Internal lightweight path always supplies BAM dictionary lines, overriding SAM
	 * header suppression. The BAM backend suppresses binary header re-emission on append.
	 * @param ffout_ Output format, append policy and ordering flag
	 * @param header_ Optional borrowed SAM header lines
	 * @param useSharedHeader_ Prefer the shared input header
	 * @param threaded_ Use one writer thread when true
	 * @param queueCapacity_ Soft queue capacity, greater than one when queued
	 * @param writeR1_ Select first mates from Read batches
	 * @param writeR2_ Select second mates from Read batches
	 * @param lightweight_ Use native lightweight output backends without subprocesses */
	SamWriterST2(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_,
			boolean threaded_, int queueCapacity_, boolean writeR1_, boolean writeR2_, boolean lightweight_){
		lightweight=lightweight_;
		
		ffout=ffout_;
		fname=ffout.name();
		header=header_;
		useSharedHeader=useSharedHeader_;
		writeR1=writeR1_;
		writeR2=writeR2_;
		final boolean nativeDictionary=(lightweight && ffout.bam());
		// Native BAM consumes dictionary lines even when append suppresses binary header emission.
		supressHeader=(!nativeDictionary && (ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists())));
		supressHeaderSequences=(!nativeDictionary && (ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader));
		
		// Config
		ordered=ffout.ordered(); // Derive ordering requirement from file format
		threaded=threaded_;
		
		// Only create queue if we need ordering or threading
		if(ordered || threaded){
			// JobQueue handles ordering, capacity bounds, and backpressure.
			// The job type is ListNum<SamLine> since SamWriter deals in lines.
			queue=new JobQueue<ListNum<SamLine>>(queueCapacity_, ordered, true, 0);
			queue.name="*SamWriterST2";
		}else{
			// Pure lightweight ST/unordered mode (needs external synchronization if MT)
			queue=null;
		}

		if(lightweight){
			outstream=LightweightOutputStream.open(ffout);
		}else if(ffout.bam()){
			outstream=ReadWrite.getBamOutputStream(fname, ffout.append());
		}else {
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}
		if(verbose){outstream2.println("Made SamWriterST2 (Ordered: "+ordered+", Threaded: "+threaded+")");}
		
//		System.err.println(ffout.ordered()+", "+ordered+", "+threaded+", "+(queue!=null));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes the header once and starts the optional worker; repeated calls return. */
	@Override
	public void start(){
		synchronized(this){
			if(started){return;}
			writeHeader();
			started=true;
			if(threaded && queue!=null){
				writerThread=new Thread(new WriterRunnable());
				writerThread.start();
			}
		}
	}

	/** Ends submissions after producers finish; queued caller-driven output drains here.
	 * Does not close the output stream or join an optional worker. */
	@Override
	public synchronized void poison(){
		// Ensure only one thread executes the shutdown sequence
		if(poisoned){return;}
		poisoned=true;
		
		if(queue!=null){
			// 1. Calculate the ID for the poison pill. It must be > maxSeen to guarantee order.
			long poisonID=queue.maxSeen()+1;
			
			// 2. Create the poison pill (ListNum acts as the HasID object)
			// SamWriter uses ListNum<SamLine>
			ListNum<SamLine> poison=new ListNum<SamLine>(null, poisonID, ListNum.POISON);
			
			// 3. Inject the poison pill into the queue using the explicit JobQueue API.
			queue.poison(poison, false); 
			
			// 4. CRITICAL: If running in Host-Driven (unthreaded) mode, the producer thread 
			// calling poison() MUST perform the final drain using the blocking take().
			if(!threaded){
				// Drain loop terminates when queue.take() returns null (the poison pill)
				while(true){
					ListNum<SamLine> job=queue.take(); 
					
					if(job==null){break;}
					
					writeLines(job);
				}
			}
			// If threaded, the worker thread handles the take() loop and exits gracefully.
		}
	}

	/** Poisons, waits for normal completion and finalizes the output.
	 * @return Sticky reported error state */
	@Override
	public synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/** Marks abandonment without draining, joining the worker or closing the stream.
	 * Setting closed also makes subsequent waitForFinish return the cached error. */
	@Override
	public synchronized void finishError(){
		errorState=true;
		if(queue!=null){
			poisoned=true;
			queue.poison(new ListNum<SamLine>(null, queue.maxSeen()+1, ListNum.POISON), true);
		}
		closed=true;
	}
	
	/** Joins an existing worker and finalizes output; call poison first for queued use.
	 * An interrupted join restores interruption and still attempts finalization, so
	 * that path is not a guarantee that the worker has terminated.
	 * @return Sticky reported error state, or its cached value when already closed */
	@Override
	public final synchronized boolean waitForFinish(){
		if(closed){return errorState;}
		// If threaded, wait for the worker to process the poison pill and exit.
		if(threaded && writerThread!=null){
			try{
				writerThread.join();
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
			}
		}
		
		boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess() && !lightweight);
		closed=true;
		return errorState|=b;
	}
	
	/** Converts and submits a nonnull Read list; null entries are skipped.
	 * @param list Borrowed Read roots, possibly paired
	 * @param id Dense zero-based batch ID in queued configurations */
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}

	/** Converts selected reads on the submitting thread and submits a fresh list.
	 * Attached SAM records may remain borrowed; ListNum poison/last flags are not forwarded.
	 * @param reads Null wrapper is ignored; otherwise its list must be nonnull */
	@Override
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		ArrayList<SamLine> lines=toSamLines(reads.list, writeR1, writeR2);
		addLines(new ListNum<SamLine>(lines, reads.id));
	}

	/** Submits SAM records directly, without applying Read mate-selection flags.
	 * Starts lazily; queued lists and records must stay stable until consumed.
	 * Caller-driven ordered use drains currently contiguous jobs after submission.
	 * @param lines Null wrapper is ignored; otherwise a nonnull borrowed list */
	@Override
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		if(!started){start();}
		
		if(queue==null){
			// LIGHTWEIGHT MODE (Unordered/Unthreaded): Direct synchronized write.
			writeLines(lines);
		}else{
			// 1. Add to queue (blocks if full/backpressure enabled)
			queue.add(lines);
			
			if(!threaded){
				// HOST-DRIVEN MODE (Ordered/Unthreaded): 
				// Producer acts as the draining consumer immediately after adding.
				// Uses the non-blocking poll() to drain any contiguous jobs that are ready.
				assert(ordered);
				synchronized(queue){//Essential to maintain ordering
					for(ListNum<SamLine> job=queue.poll(); job!=null; job=queue.poll()){
						if(verbose){System.err.println("Writing job "+job.id+", expected="+expected);}
						writeLines(job);
					}
				}
			}
			// If threaded, the background thread handles the draining.
		}
	}
	
	/** Serializes a batch and counts nonnull SAM records, including secondary alignments.
	 * Counts advance during formatting and do not certify durable output.
	 * @param lines Batch selected by the active direct or queued path */
	private void writeLines(ListNum<SamLine> lines){
		assert(!ordered || lines.id()==expected++) : lines.id+", "+expected;
		ByteBuilder bb=new ByteBuilder();
		for(SamLine sl : lines){
			if(sl==null){continue;}
			sl.toBytes(bb);
			bb.nl();
			readsWritten++;
			basesWritten+=sl.length();

			if(bb.length()>=BUFFER_SIZE){
				write(bb);
			}
		}
		if(bb.length()>0){
			write(bb);
		}
	}
	
	/** Writes buffered bytes under the stream lock, then clears the builder.
	 * @param bb Buffer owned by this formatting invocation */
	private void write(ByteBuilder bb){
		if(bb.length()==0){return;}//was <0 (never fired); skip empty buffers
		byte[] array=bb.toBytes();
		try{
			// CRITICAL: Synchronize on the shared stream object to prevent file corruption
			synchronized(outstream){outstream.write(array);}
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}
	
	/** @return Observed count of formatted SAM records, including secondary alignments */
	@Override
	public long readsWritten(){return readsWritten;}

	/** @return Sum of formatted SAM record lengths, including secondary alignments */
	@Override
	public long basesWritten(){return basesWritten;}

	/** @return Current error flag; this accessor is not a completion barrier */
	@Override
	public boolean errorState(){return errorState;}
	
	/** @return Whether closed is set and the cached error flag is clear */
	@Override
	public boolean finishedSuccessfully(){return !errorState && closed;}
	
	/** @return Output filename supplied by the format */
	@Override
	public final String fname(){return fname;}

	/*--------------------------------------------------------------*/
	/*----------------        Writer Thread         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Consumes queued SAM batches on the optional single formatting/writing thread. */
	private class WriterRunnable implements Runnable{
		/** Drains normal batches until the queue supplies its terminal signal. */
		@Override
		public void run(){
			try{
				// Loop uses the blocking take() and terminates when it retrieves the poison pill (null)
				for(ListNum<SamLine> job=queue.take(); job!=null; job=queue.take()){
					writeLines(job);
				}
			}catch(Throwable t){
				//[stream/SamWriterST2#001] A write failure (write()/writeLines rethrow IOException) must not
				//strand the producer. Set errorState (so waitForFinish reports loud), then FORCE-poison the
				//queue (force=true sets JobQueue.poisoned so a producer blocked in add()'s capacity-wait
				//unblocks instead of hanging), then surface loud. The blessed SWEEP-A writer-death fix.
				errorState=true;
				if(queue!=null){
					try{queue.poison(new ListNum<SamLine>(null, queue.maxSeen()+1, ListNum.POISON), true);}catch(Throwable t2){}
				}
				throw new RuntimeException("SamWriterST2 writer thread failed; output may be incomplete.", t);
			}
		}
	}


	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts all selected records, retaining attached SAM metadata when enabled.
	 * @param reads Batch roots; null entries are skipped
	 * @return SAM records for both mate selections */
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads){
		return toSamLines(reads, true, true);
	}

	/** Converts roots using the same pair-number selection as FASTQ writers.
	 * R1 is a pairnum-zero root; R2 is a pairnum-one root or the mate of an R1 root.
	 * Attached SAM records may be reused; optional secondary alignments are added.
	 * @param reads Batch roots; null entries are skipped
	 * @param writeR1 Select first mates
	 * @param writeR2 Select second mates
	 * @return Selected SAM records in root order, with R1 before R2 when both occur */
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads, boolean writeR1, boolean writeR2){
		ArrayList<SamLine> samLines=new ArrayList<SamLine>();
		for(final Read r : reads){
			if(r==null){continue;}
			// Fixed STR240: a standalone R2 occupies a batch slot but is not an R1.
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			SamLine sl1=(r1==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE && r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE && r2.samline!=null ? r2.samline : new SamLine(r2, 1)));
			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
				sl2.qname=sl1.qname;
			}
			assert(sl1!=null || sl2!=null) : "A nonnull Read must supply one SAM record before mate selection: "+r;
			if(writeR1){addSamLine(r1, sl1, samLines);}
			if(writeR2){addSamLine(r2, sl2, samLines);}
		}
		return samLines;
	}

	/** Appends a primary and enabled secondary alignments without copying the primary.
	 * @param r Source read; null is ignored
	 * @param primary Primary SAM record; null is ignored
	 * @param samLines Destination list in output order */
	private static void addSamLine(Read r, SamLine primary, ArrayList<SamLine> samLines){
		if(r==null || primary==null){return;}
		assert(!ReadStreamWriter.ASSERT_CIGAR || !r.mapped() || primary.cigar!=null) : r;
		samLines.add(primary);
		ArrayList<SiteScore> list=r.sites;
		if(ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
			final Read clone=r.clone();
			for(int i=1; i<list.size(); i++){
				SiteScore ss=list.get(i);
				clone.match=null;
				clone.setFromSite(ss);
				clone.setSecondary(true);
				SamLine secondary=new SamLine(clone, r.pairnum());
				assert(!secondary.nonSecondary());
				assert(!ReadStreamWriter.USE_ATTACHED_SAMLINE || secondary.cigar!=null) : r;
				samLines.add(secondary);
			}
		}
	}
	
	/** Selects shared, supplied or generated headers, in that priority order.
	 * A requested but absent shared header becomes an empty list with a warning.
	 * @return Borrowed header list, or a newly generated/empty list */
	ArrayList<byte[]> getHeader(){
		ArrayList<byte[]> headerLines;
		if(useSharedHeader){
			headerLines=SamReadInputStream.getSharedHeader(true);
		}else if(header!=null){
			headerLines=header;
		}else {
			headerLines=SamHeader.makeHeaderList(supressHeaderSequences, 
				ReadStreamWriter.MINCHROM, ReadStreamWriter.MAXCHROM);
		}
		if(headerLines==null){
			outstream2.println("Warning: Header was null, creating empty header");
			headerLines=new ArrayList<byte[]>();
		}
		return headerLines;
	}

	/** Writes unsuppressed header lines once; native BAM receives dictionary text too. */
	protected void writeHeader(){
		if(headerWritten || supressHeader){return;}
		ArrayList<byte[]> headerLines=getHeader();
		
		ByteBuilder bb=new ByteBuilder();
		try{
			for(byte[] line : headerLines){
				bb.append(line).nl();
				if(bb.length()>=16384){
					outstream.write(bb.toBytes());
					bb.clear();
				}
			}
			if(bb.length()>=1){
				outstream.write(bb.toBytes());
				bb.clear();
			}
		}catch(IOException e){
			throw new RuntimeException(e);
		}
		headerWritten=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output path and opening policy captured at construction. */
	final String fname;
	final FileFormat ffout;
	/** Prefer shared input metadata when choosing the header. */
	final boolean useSharedHeader;
	/** Suppress text header output except required native BAM dictionary input. */
	final boolean supressHeader;
	/** Omit sequence entries when generating a header. */
	final boolean supressHeaderSequences;
	/** Whether this instance has completed writing its unsuppressed header. */
	boolean headerWritten=false;
	/** Optional caller-owned explicit header; shared metadata takes precedence. */
	final ArrayList<byte[]> header;
	/** Open output backend, finalized through ReadWrite. */
	final OutputStream outstream;
	
	/** Output-format ordering flag; JobQueue also orders threaded unordered requests. */
	final boolean ordered;
	/** Whether one worker consumes queued SAM batches. */
	final boolean threaded;
	/** Restricts output compression to the leaf's formatting thread. */
	private final boolean lightweight;
	/** Present for ordered or threaded operation; null for direct caller-driven output. */
	final JobQueue<ListNum<SamLine>> queue;
	/** Optional worker, created by start. */
	Thread writerThread;
	
	/** Formatted records and their lengths; uncoordinated reads are only observations. */
	private long readsWritten=0;
	private long basesWritten=0;
	/** Sticky error indication from worker, abandonment or output finalization. */
	private boolean errorState=false;
	/** Write R1 records (pairnum 0). */
	private final boolean writeR1;
	/** Write R2 records (pairnum 1). */
	private final boolean writeR2;
	/** Header/start transition has completed. */
	private boolean started=false;
	/** Completion or explicit abandonment has been recorded. */
	private boolean closed=false;
	/** Normal or forced termination has been requested. */
	private boolean poisoned=false;

	/** Expected next batch number, advanced only by enabled ordered assertions. */
	private long expected=0;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Byte threshold at which a formatted SAM batch flushes its builder. */
	private static final int BUFFER_SIZE=65536;
	/** Compile-time diagnostic output switch. */
	public static final boolean verbose=false;
	/** Destination for constructor and missing-header diagnostics. */
	protected PrintStream outstream2=System.err;

}
