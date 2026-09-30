package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Formats SAM text with parallel workers and a separate output thread.
 *
 * Construction opens the selected SAM/BAM output through ReadWrite. Workers retain
 * submitted Read/SamLine references until conversion and serialization complete;
 * keep those objects, header references and shared SAM settings stable while in use.
 * OrderedQueueSystem2 receives the descriptor's ordering preference, but current
 * JobQueue forces strict ordering for both queues. Supply dense batch IDs starting
 * at zero even when the descriptor is unordered; empty batches can preserve IDs.
 *
 * Start once, submit batches, then coordinate normal poisonAndWait completion.
 * There is no duplicate-start guard. Public counters accumulate worker formatting
 * totals during normal output finalization, rather than reporting live write progress.
 *
 * @author Brian Bushnell, Isla
 * @date October 25, 2025
 */
public class SamWriter implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens output using default queue capacities and selects list entries plus their mates.
	 * @param ffout_ Nonnull output descriptor
	 * @param threads_ Requested workers; values below one use DEFAULT_THREADS before clamping
	 * @param header_ Optional borrowed header list
	 * @param useSharedHeader_ Prefer the shared input header when available
	 */
	public SamWriter(FileFormat ffout_, int threads_,
		ArrayList<byte[]> header_, boolean useSharedHeader_){
		this(ffout_, threads_, header_, useSharedHeader_, true, true);
	}

	/** Opens output with explicit conversion selection and default queue capacities.
	 * @param ffout_ Nonnull output descriptor
	 * @param threads_ Requested workers; values below one use DEFAULT_THREADS before clamping
	 * @param header_ Optional borrowed header list
	 * @param useSharedHeader_ Prefer shared header over provided/generated headers
	 * @param writeR1_ Emit list-entry records during Read conversion
	 * @param writeR2_ Emit attached-mate records during Read conversion
	 */
	public SamWriter(FileFormat ffout_, int threads_,
		ArrayList<byte[]> header_, boolean useSharedHeader_, boolean writeR1_, boolean writeR2_){
		this(ffout_, threads_, header_, useSharedHeader_, writeR1_, writeR2_, 0, 0);
	}

	/** Captures configuration, constructs queues and opens output without starting threads.
	 * Worker count is clamped between one and Shared.threads. Header suppression captures
	 * current globals and append-to-existing-file status; header contents are selected later.
	 * Output opening delegates to ReadWrite's BAM path for BAM, otherwise its buffered path.
	 * @param ffout_ Nonnull output descriptor; its append and ordering settings are retained
	 * @param threads_ Requested workers; values below one select DEFAULT_THREADS
	 * @param header_ Optional borrowed header list and byte arrays
	 * @param useSharedHeader_ Request shared-header lookup before provided/generated headers
	 * @param writeR1_ Select list entries, irrespective of their stored pair-number bit
	 * @param writeR2_ Select attached mates during Read conversion
	 * @param inputCapacity_ Positive explicit capacity; nonpositive chooses the worker-based default
	 * @param outputCapacity_ Positive explicit capacity; nonpositive chooses the worker-based default
	 * @throws IllegalArgumentException If either resulting queue capacity is less than two
	 */
	public SamWriter(FileFormat ffout_, int threads_,
		ArrayList<byte[]> header_, boolean useSharedHeader_, boolean writeR1_, boolean writeR2_,
		int inputCapacity_, int outputCapacity_){
		ffout=ffout_;
		fname=ffout.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		header=header_;
		useSharedHeader=useSharedHeader_;
		writeR1=writeR1_;
		writeR2=writeR2_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);

		//Create prototype jobs for OrderedQueueSystem2
		SamWriterInputJob inputProto=new SamWriterInputJob(null, null, ListNum.PROTO, -1);
		SamWriterOutputJob outputProto=new SamWriterOutputJob(-1, null, ListNum.PROTO);

		final int inputCapacity=(inputCapacity_>0 ? inputCapacity_ :
				threads+OrderedQueueSystem2.BUFFER_PADDING);
		final int outputCapacity=(outputCapacity_>0 ? outputCapacity_ :
				(OrderedQueueSystem2.BUFFER_MULT*threads)/2+OrderedQueueSystem2.BUFFER_PADDING);
		if(inputCapacity<2 || outputCapacity<2){
			throw new IllegalArgumentException("Writer queue capacities must be at least 2: "+
					inputCapacity+", "+outputCapacity);
		}
		oqs=new OrderedQueueSystem2<SamWriterInputJob, SamWriterOutputJob>(
			inputCapacity, outputCapacity, threads, ffout.ordered(), inputProto, outputProto);

		if(ffout.bam()){
			outstream=ReadWrite.getBamOutputStream(fname, ffout.append());
		}else{
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}
		if(verbose){System.err.println("outstream="+outstream.getClass());}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts a fresh set of workers and one output thread; call once per writer instance.
	 * Repeated calls are not rejected and do not reset queues, headers or counters.
	 */
	@Override
	public void start(){spawnThreads();}

	/** Wraps and submits borrowed reads with the queue's dense ID sequence starting at zero.
	 * ListNum construction may assign Read.rand under its global deterministic-number option.
	 * @param list Nonnull input list kept stable until worker consumption completes
	 * @param id Batch ID; submit an empty list rather than omitting an expected ID
	 */
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}

	/** Enqueues retained reads for later worker conversion and serialization.
	 * Uses a NORMAL input job with the wrapper ID, not the wrapper's marker type.
	 * Does not start threads; null wrappers are ignored. Queue submission may block.
	 * @param reads Wrapper with a nonnull list, or null to do nothing
	 */
	@Override
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		SamWriterInputJob job=new SamWriterInputJob(reads, null, ListNum.NORMAL, reads.id);
		oqs.addInput(job);
	}

	/** Enqueues borrowed SAM records for worker serialization without Read conversion.
	 * Bypasses writeR1/writeR2 selection and secondary expansion; supplied records are
	 * serialized as listed, skipping null entries. Wrapper marker type is not forwarded.
	 * @param lines Wrapper with a stable nonnull list, or null to do nothing
	 */
	@Override
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		SamWriterInputJob job=new SamWriterInputJob(null, lines, ListNum.NORMAL, lines.id);
		oqs.addInput(job);
	}

	/** Delegates normal end-of-input signaling to OQS2 after submissions are coordinated.
	 * Queues its final output/input markers; does not itself wait for output completion.
	 */
	@Override
	public final void poison(){
		oqs.poison();
	}

	/** Waits for OQS2's finished flag, then reads this writer's error flag.
	 * Does not initiate poison. Normal completion follows worker aggregation and cleanup;
	 * forced completion can publish the queue flag without joining every thread here.
	 * @return Current writer error flag
	 */
	public final boolean waitForFinish(){
		oqs.waitForFinish();
		return errorState();
	}

	/** Signals normal termination and waits according to waitForFinish's contract.
	 * @return Current writer error flag
	 */
	public final boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/** Marks an external pipeline failure and requests OQS2 forced completion.
	 * Does not join worker/output threads or perform stream cleanup itself.
	 */
	@Override
	public synchronized void finishError(){
		setErrorState(true);
		oqs.setFinished(true);
	}

	/** Returns the observed worker-formatting count accumulated during normal finalization. */
	@Override
	public long readsWritten(){return readsWritten;}

	/** Returns the observed SamLine-length total accumulated during normal finalization. */
	@Override
	public long basesWritten(){return basesWritten;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces the thread list, creates output ID zero and worker IDs one through threads,
	 * then starts each. No restart/reset or duplicate-start protection is performed.
	 */
	void spawnThreads(){
		final int totalThreads=threads+1; //Workers plus writer

		alpt=new ArrayList<ProcessThread>(totalThreads);
		for(int i=0; i<totalThreads; i++){alpt.add(new ProcessThread(i));}

		for(ProcessThread pt : alpt){pt.start();}
	}

	/** Converts all list entries and attached mates, including configured secondary alignments.
	 * @param reads Nonnull list; null entries are skipped
	 * @return New list of generated or reused SAM records
	 */
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads){
		return toSamLines(reads, true, true);
	}

	/** Builds SAM records and applies positional entry/mate selection.
	 * May reuse attached SamLines under USE_ATTACHED_SAMLINE. Before selection, KEEP_NAMES=false
	 * can replace the mate record's qname with the entry's qname, including on reused objects.
	 * Selected records include configured secondary alignments through addSamLine.
	 * @param reads Nonnull source list with optional attached mates; null entries are skipped
	 * @param writeR1 Include list-entry records, not a filter on Read.pairnum()
	 * @param writeR2 Include attached-mate records
	 * @return Newly allocated list containing generated and/or borrowed SamLines
	 */
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads, boolean writeR1, boolean writeR2){
		ArrayList<SamLine> samLines=new ArrayList<SamLine>();

		for(final Read r1 : reads){
			if(r1==null){continue;}
			Read r2=(r1==null ? null : r1.mate);

			SamLine sl1=(r1==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE
				&& r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE
				&& r2.samline!=null ? r2.samline : new SamLine(r2, 1)));

			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) ||
				!sl2.qname.equals(sl1.qname))){
				sl2.qname=sl1.qname;
			}
			assert(sl1!=null) : r1;
			if(writeR1){addSamLine(r1, sl1, samLines);}
			if(writeR2){addSamLine(r2, sl2, samLines);}
		}
		return samLines;
	}

	/** Appends a primary record and optional secondary records built from later mapping sites.
	 * Secondary conversion uses a cloned Read; a null Read or primary is ignored.
	 * @param r Source read supplying mapping sites
	 * @param primary Generated or attached primary SAM record
	 * @param samLines Destination list
	 */
	private static void addSamLine(Read r, SamLine primary, ArrayList<SamLine> samLines){
		if(r==null || primary==null){return;}

		assert(!ReadStreamWriter.ASSERT_CIGAR || !r.mapped() || primary.cigar!=null) : r;
		samLines.add(primary);

		// Handle secondary alignments
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

	/** Selects shared, then provided, then generated header lines; an empty list is accepted.
	 * Shared lookup may wait for a publisher. Sequence suppression and current chromosome
	 * bounds apply only to generation. Retained/shared lists are returned without copying.
	 * @return Nonnull header list; an unexpected null generation result becomes an empty list
	 */
	ArrayList<byte[]> getHeader(){
		if(verbose){System.err.println("Fetching header: "+useSharedHeader+","+(header!=null));}
		ArrayList<byte[]> headerLines=null;
		//getSharedHeader may now return null when no SAM input exists (non-sam -> sam); fall through to
		//an explicit header, then to a generated one, rather than blocking or emitting an empty header.
		if(useSharedHeader){headerLines=SamReadInputStream.getSharedHeader(true);}
		if(headerLines==null && header!=null){headerLines=header;}
		if(headerLines==null){
			headerLines=SamHeader.makeHeaderList(supressHeaderSequences,
				ReadStreamWriter.MINCHROM, ReadStreamWriter.MAXCHROM);
		}
		if(headerLines==null){
			System.err.println("Warning: Header was null, creating empty header");
			headerLines=new ArrayList<byte[]>();
		}
		if(verbose){System.err.println("Fetched header: "+(headerLines==null ? "null" : headerLines.size()));}
		return headerLines;
	}

	/** Writes selected header lines once unless all-header suppression was captured.
	 * Appends newlines and writes accumulated bytes at a 16384-byte threshold. Marks written
	 * only after output succeeds; suppressed calls return without changing headerWritten.
	 */
	protected synchronized void writeHeader(){
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
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input job retaining one Read or SamLine wrapper for NORMAL jobs, or a marker. */
	static class SamWriterInputJob implements HasID{

		/** Retains wrappers; enabled assertions check normal-job exclusivity, marker type and IDs. */
		public SamWriterInputJob(ListNum<Read> reads_, ListNum<SamLine> lines_, int type_, long id_){
			reads=reads_;
			lines=lines_;
			type=type_;
			id=id_;
			assert(type!=ListNum.NORMAL || ((reads==null)!=(lines==null)));
			assert(type==ListNum.NORMAL || type==ListNum.LAST || type==ListNum.POISON || type==ListNum.PROTO);
			assert(reads==null || id==reads.id);
			assert(lines==null || id==lines.id);
		}

		/** Returns the queue ID. */
		@Override
		public long id(){return id;}

		/** Reports the input poison marker type. */
		@Override
		public boolean poison(){return type==ListNum.POISON;}

		/** Reports the last-job marker type. */
		@Override
		public boolean last(){return type==ListNum.LAST;}

		/** Creates a poison marker with the requested ID and no input wrapper. */
		@Override
		public SamWriterInputJob makePoison(long id_){
			return new SamWriterInputJob(null, null, ListNum.POISON, id_);
		}

		/** Creates a last marker with the requested ID and no input wrapper. */
		@Override
		public SamWriterInputJob makeLast(long id_){
			return new SamWriterInputJob(null, null, ListNum.LAST, id_);
		}

		/** Optional borrowed Read wrapper; factory-created markers omit it. */
		public final ListNum<Read> reads;
		/** Optional borrowed SamLine wrapper; factory-created markers omit it. */
		public final ListNum<SamLine> lines;
		/** NORMAL, LAST, POISON or PROTO. */
		public final int type;
		/** ID copied from a normal wrapper or supplied for a marker. */
		public final long id;
	}

	/** Formatted output bytes with an ordering ID, or a payload-free marker. */
	static class SamWriterOutputJob implements HasID{

		/** Retains bytes/ID/type; enabled assertions check marker type and payload consistency. */
		public SamWriterOutputJob(long id_, byte[] bytes_, int type_){
			id=id_;
			bytes=bytes_;
			type=type_;
			assert((type==ListNum.NORMAL)==(bytes!=null));
			assert(type==ListNum.NORMAL || type==ListNum.LAST || type==ListNum.PROTO || type==ListNum.POISON);
		}

		/** Returns the output ordering ID. */
		@Override
		public long id(){return id;}

		/** Reports the poison marker type. */
		@Override
		public boolean poison(){return type==ListNum.POISON;}

		/** Reports normal end-of-output marker type. */
		@Override
		public boolean last(){return type==ListNum.LAST;}

		/** Creates a payload-free poison marker with the requested ID. */
		@Override
		public SamWriterOutputJob makePoison(long id_){
			return new SamWriterOutputJob(id_, null, ListNum.POISON);
		}

		/** Creates a payload-free last marker with the requested ID. */
		@Override
		public SamWriterOutputJob makeLast(long id_){
			return new SamWriterOutputJob(id_, null, ListNum.LAST);
		}

		/** ID preserved from an input job or supplied for a marker. */
		public final long id;
		/** Owned formatted byte snapshot; null for markers. */
		public final byte[] bytes;
		/** NORMAL, LAST, POISON or PROTO. */
		public final int type;
	}

	/** Processing thread - converts reads/lines to SAM text or writes output. */
	private class ProcessThread extends Thread{

		/** Names a thread for its role; ID zero writes output and positive IDs format jobs.
		 * @param tid_ Role ID assigned by spawnThreads
		 */
		ProcessThread(final int tid_){
			tid=tid_;
			setName("SamLineWriter-"+(tid==0 ? "Output" : "Worker-"+tid));
		}

		/** Executes the assigned role while holding this thread's monitor.
		 * A returned role method sets success. Throwables retain the historical error handling
		 * below: record failure, request forced queue completion and rethrow with the thread ID.
		 */
		@Override
		public void run(){
			try{
				synchronized(this){
					if(tid==0){
						writeOutput(); //Writer thread
					}else{
						processJobs(); //Worker thread
					}
					success=true;
				}
			}catch(Throwable t){
				//[stream/SamWriter#001 FIXED 2026-06-21]: a thread death without poisoning/finishing the OQS2 HANGS -- a writer death
				//(write error in writeOutput/writeHeader) leaves outq undrained so workers block in addOutput + `finished` unset so main
				//hangs in waitForFinish; a worker death (a throw in processJobs/toSamLines/toBytes) leaves its ordered job undelivered so
				//the writer's getOutput blocks forever on the gap. Escalate to a LOUD finish: setErrorState + oqs.setFinished(true)
				//[force=true poisons inq+outq, waking writer+workers+main], then rethrow so the thread dies loud (stack trace) instead of
				//hanging. errorState surfaces via poisonAndWait's return so the caller decides the exit. Mirrors BamWriter#001/FastqWriter#001.
				setErrorState(true);
				oqs.setFinished(true);
				throw new RuntimeException("SamWriter thread "+tid+" failed", t);
			}
		}

		/** Writes the header and queued byte blocks, then joins other threads and aggregates counts.
		 * ThreadWaiter skips this output thread. Normal flow combines worker success and cleanup
		 * status before publishing non-forced queue completion. Counts come from worker formatting.
		 */
		void writeOutput(){
			//Write header first
			writeHeader();

			//Write ordered data blocks
			SamWriterOutputJob job=oqs.getOutput();
			while(job!=null && !job.last()){
				try{
					outstream.write(job.bytes);
				}catch(Exception e){
					throw new RuntimeException("Error writing output", e);
				}

				job=oqs.getOutput();
			}

			//Wait for other threads and accumulate statistics
			ThreadWaiter.waitForThreadsToFinish(alpt);
			synchronized(SamWriter.this){
				for(ProcessThread pt : alpt){
					if(pt!=this){
						synchronized(pt){
							readsWritten+=pt.readsWrittenT;
							basesWritten+=pt.basesWrittenT;
							setErrorState(!pt.success);
						}
					}
				}
			}

			if(verbose){System.err.println("Consumer finished accumulating.");}
			boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			errorState|=b;
			if(verbose){System.err.println("Consumer finished writing.");}
			oqs.setFinished(false);
			if(verbose){System.err.println("Consumer set oqs finished.");}
		}

		/** Converts queued Read wrappers or serializes supplied SamLines into copied byte blocks.
		 * Direct SamLine jobs bypass mate selection/secondary expansion. Counts each nonnull SAM
		 * record and SamLine.length before output submission, clears the reused builder afterward,
		 * and returns a received poison job to the input queue for remaining workers.
		 */
		void processJobs(){
			final ByteBuilder bb=new ByteBuilder();

			SamWriterInputJob job=oqs.getInput();
			while(job!=null && !job.poison()){
				//Convert to SamLines if needed
				ArrayList<SamLine> lines;
				if(job.lines!=null){
					lines=job.lines.list;
				}else{
					lines=toSamLines(job.reads.list, writeR1, writeR2);
				}

				//Format SamLines to bytes and count
				for(SamLine sl : lines){
					if(sl==null){continue;}
					sl.toBytes(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=sl.length();
				}

				//Create output job
				SamWriterOutputJob outJob=new SamWriterOutputJob(job.id(), bb.toBytes(), ListNum.NORMAL);
				oqs.addOutput(outJob);
				bb.clear();

				job=oqs.getInput();
			}

			//Re-inject poison for other workers
			if(job!=null){oqs.addInput(job);}
		}

		/** Number of SAM records formatted by this worker, including secondary records. */
		protected long readsWrittenT=0;
		/** Sum of formatted SamLine lengths. */
		protected long basesWrittenT=0;
		/** Set true when the selected role method returns without throwing. */
		boolean success=false;
		/** Role ID: zero for output, positive for conversion workers. */
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Getters and Setters      ----------------*/
	/*--------------------------------------------------------------*/

	/** Accumulates an error observation without clearing an earlier error.
	 * @param b Error observation to combine
	 */
	synchronized void setErrorState(boolean b){
		errorState|=b;
	}

	/** Returns the current local error flag while holding this writer's monitor. */
	@Override
	public synchronized boolean errorState(){return errorState;}

	/** Reads local error and queue completion flags; does not join or finalize threads. */
	@Override
	public boolean finishedSuccessfully(){return !errorState && oqs.finished();}

	/** Returns the output name captured during construction. */
	@Override
	public final String fname(){return fname;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output file name. */
	final String fname;
	/** Retained output descriptor, also consulted during cleanup. */
	final FileFormat ffout;
	/** Selected conversion-worker count, excluding the output thread. */
	final int threads;
	/** Whether shared-header lookup precedes provided/generated headers. */
	final boolean useSharedHeader;
	/** Captured all-header suppression, including append to existing output. */
	final boolean supressHeader;
	/** Captured suppression of sequence lines in generated headers only. */
	final boolean supressHeaderSequences;
	/** Set after successful header output; suppression alone does not set it. */
	boolean headerWritten=false;
	/** Optional borrowed header list and byte arrays. */
	final ArrayList<byte[]> header;
	/** Queue coordinator; current JobQueue requires dense IDs starting at zero. */
	final OrderedQueueSystem2<SamWriterInputJob, SamWriterOutputJob> oqs;
	/** Output stream opened during construction. */
	final OutputStream outstream;

	/** Output and worker threads created by the latest start call. */
	private ArrayList<ProcessThread> alpt;

	/** Public accumulator of worker-formatted SAM records, merged during normal finalization. */
	public long readsWritten=0;
	/** Public accumulator of worker SamLine lengths, merged during normal finalization. */
	public long basesWritten=0;
	/** Local accumulated error flag. */
	private boolean errorState=false;
	/** Include source-list entries during Read conversion; not a pair-number filter. */
	private final boolean writeR1;
	/** Include attached mates during Read conversion; direct SamLine jobs bypass selection. */
	private final boolean writeR2;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Default requested worker count before clamping to Shared.threads. */
	public static int DEFAULT_THREADS=6;

	/** Shared diagnostic switch. */
	public static final boolean verbose=false;

	/** Retained diagnostic stream; current implementation prints directly to System.err. */
	protected PrintStream outstream2=System.err;

}
