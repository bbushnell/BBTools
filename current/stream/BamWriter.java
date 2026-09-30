package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import stream.bam.SamToBamConverter;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Writes BAM binary files with parallel conversion and ordered output.
 * 
 * Workers convert borrowed Read/SamLine objects to binary records. The output thread
 * selects a header, creates a converter and writes queued byte blocks. Construction
 * opens output through ReadWrite's configured BGZF path, which may open backend resources.
 * Keep input/header references and shared conversion settings stable while in use.
 *
 * Current JobQueue forces ordered processing even when the descriptor requests unordered
 * output. Supply dense IDs starting at zero; empty batches can preserve expected IDs.
 * Start once, coordinate submissions, then use normal poisonAndWait completion. Repeated
 * start is not guarded. Counters merge worker conversion totals during normal finalization.
 * 
 * @author Brian Bushnell, Isla
 * @date October 25, 2025
 */
public class BamWriter implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens output and selects list entries plus attached mates during Read conversion.
	 * @param ffout_ Nonnull output descriptor
	 * @param threads_ Requested workers; below one uses DEFAULT_THREADS before clamping
	 * @param header_ Optional borrowed header used only outside shared-header mode
	 * @param useSharedHeader_ Select shared-header mode; unavailable shared data becomes empty
	 */
	public BamWriter(FileFormat ffout_, int threads_,
		ArrayList<byte[]> header_, boolean useSharedHeader_){
		this(ffout_, threads_, header_, useSharedHeader_, true, true);
	}

	/** Captures configuration, creates queues and opens output with append forwarded.
	 * Worker count is clamped between one and Shared.threads. Header suppression captures
	 * current globals and append-to-existing-file status. Appending relies on the selected
	 * header matching the existing reference dictionary. Does not start this writer's threads.
	 * @param ffout_ Nonnull descriptor supplying output name, append and ordering preference
	 * @param threads_ Requested workers; below one uses DEFAULT_THREADS
	 * @param header_ Borrowed header list and arrays, used when shared-header mode is false
	 * @param useSharedHeader_ Select shared-header lookup rather than provided/generated headers
	 * @param writeR1_ Select list-entry records during Read conversion, not a pair-bit filter
	 * @param writeR2_ Select attached-mate records during Read conversion
	 */
	public BamWriter(FileFormat ffout_, int threads_,
		ArrayList<byte[]> header_, boolean useSharedHeader_, boolean writeR1_, boolean writeR2_){
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

		oqs=new OrderedQueueSystem2<SamWriterInputJob, SamWriterOutputJob>(
			threads, ffout.ordered(), inputProto, outputProto);
		
		//ffout.append() must be honored: app=t previously TRUNCATED while supressHeader (line above) simultaneously
		//suppressed the header, yielding a headerless corrupt BAM (replicated via stream.sh 2026-09-05)
		outstream=ReadWrite.getBgzipStream(fname, ffout.append());
		if(verbose){System.err.println("outstream="+outstream.getClass());}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts a fresh output thread and worker set; call once per instance.
	 * Repeated calls are not rejected and do not reset counters, queues or header state.
	 */
	@Override
	public void start(){spawnThreads();}
	
	/** Wraps and submits retained reads with dense IDs starting at zero.
	 * ListNum construction may assign Read.rand under its deterministic-number option.
	 * @param list Nonnull list kept stable until worker consumption completes
	 * @param id Expected batch ID; use empty lists instead of omitting IDs
	 */
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}

	/** Enqueues borrowed reads for worker conversion and BAM serialization.
	 * Uses a NORMAL job with the wrapper ID, not its marker type. Does not start threads.
	 * @param reads Wrapper with a stable nonnull list, or null to do nothing
	 */
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		SamWriterInputJob job=new SamWriterInputJob(reads, null, ListNum.NORMAL, reads.id);
		oqs.addInput(job);
	}

	/** Enqueues borrowed SAM records for binary serialization, skipping Read conversion.
	 * Bypasses mate selection and secondary expansion; wrapper marker type is not forwarded.
	 * @param lines Wrapper with a stable nonnull list, or null to do nothing
	 */
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		SamWriterInputJob job=new SamWriterInputJob(null, lines, ListNum.NORMAL, lines.id);
		oqs.addInput(job);
	}

	/** Delegates normal end-of-input markers to OQS2 after submissions are coordinated.
	 * Does not itself wait for output completion.
	 */
	public final void poison(){
		oqs.poison();
	}

	/** Waits for OQS2's finished flag, then reports the local error flag.
	 * Does not initiate poison. Forced queue completion is not a join of every thread here.
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

	/** Records external failure, notifies converter waiters and requests forced queue completion.
	 * Does not join threads or perform stream cleanup itself.
	 */
	@Override
	public synchronized void finishError(){
		//Workers may be parked in getConverter() on this monitor rather than in
		//OQS2; wake them using the same header-failure condition as failHeader().
		headerFailed=true;
		setErrorState(true);
		this.notifyAll();
		oqs.setFinished(true);
	}

	/** Returns observed worker-record totals accumulated during normal finalization. */
	@Override
	public long readsWritten(){return readsWritten;}

	/** Returns observed worker SamLine-length totals accumulated during normal finalization. */
	@Override
	public long basesWritten(){return basesWritten;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces the thread list, creates output role zero plus conversion workers, and starts all.
	 * No duplicate-start protection or lifecycle-state reset is performed.
	 */
	void spawnThreads(){
		final int totalThreads=threads+1; //Workers plus writer

		alpt=new ArrayList<ProcessThread>(totalThreads);
		for(int i=0; i<totalThreads; i++){alpt.add(new ProcessThread(i));}

		for(ProcessThread pt : alpt){pt.start();}
	}

	/** Builds converter metadata, delegates binary-header output and publishes the converter.
	 * Extracts SN keys from selected SQ lines and collects matching full-name aliases. Converter
	 * setup still occurs when header bytes are suppressed. On success marks headerWritten and
	 * notifies waiters; selected setup/output exceptions use failHeader before propagation.
	 * Dictionary agreement concerns are recorded separately below (#003).
	 */
	protected synchronized void writeHeader(){
		if(headerWritten){return;}

		ArrayList<byte[]> headerLines=getHeader();

		// Extract reference names for converter
		ArrayList<String> refNames=new ArrayList<String>();
		for(byte[] line : headerLines){
			if(line.length>3 && line[0]=='@' && line[1]=='S' && line[2]=='Q'){
				String lineStr=new String(line, StandardCharsets.US_ASCII);
				String[] fields=lineStr.split("\\t");
				for(int i=1; i<fields.length; i++){
					if(fields[i].startsWith("SN:")){
						refNames.add(fields[i].substring(3));
						break;
					}
				}
			}
		}
		//TODO: Probable bug #003 - these converter IDs include every extracted SN, but BamWriterHelper
		//emits dictionary entries only for SN with positive LN. Sequence suppression emits zero entries
		//when a header is written; all-header suppression emits no header. Dictionary agreement needs
		//a separate review before changing behavior.
		//Validate reference-name mapping and aliases before emitting any header bytes.
		final SamToBamConverter converter;
		try{
			converter=new SamToBamConverter(refNames.toArray(new String[0]),
					collectReferenceAliases());
		}catch(RuntimeException|Error e){
			failHeader();
			throw e;
		}

		// Write BAM header using helper
		stream.bam.BamWriterHelper writer=new stream.bam.BamWriterHelper(outstream);
		try{
			writer.writeHeaderFromLines(headerLines, supressHeader, supressHeaderSequences);
		}catch(IOException e){
			//[stream/BamWriter#001 FIXED 2026-06-20 (greenlit + verified)]: header-write failure here dies BEFORE sharedConverter
			//is set + notifyAll -> workers would block forever in getConverter(). Set headerFailed + notifyAll so those workers
			//wake and abort loudly instead of hanging, THEN crash loud. (setFinished(true) alone cannot wake them: getConverter
			//waits on the BamWriter monitor, not the OQS2 monitor that setFinished notifies.)
			failHeader();
			throw new RuntimeException(e);
		}

		// Create converter for workers
		sharedConverter=converter;
		headerWritten=true;
		this.notifyAll();
	}

	/** Records header failure, notifies converter waiters and delegates output cleanup.
	 * Called with this writer's monitor held by writeHeader; cleanup may block or throw.
	 */
	private void failHeader(){
		headerFailed=true;
		setErrorState(true);
		this.notifyAll();
		ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
	}

	/** Collects available full scaffold names containing whitespace when TRIM_RNAME is enabled.
	 * The converter associates an alias only when its shortened name exists in the SN map.
	 * @return Collected aliases, or null when disabled or no qualifying names exist
	 */
	private static String[] collectReferenceAliases(){
		if(!Shared.TRIM_RNAME || Data.scaffoldNames==null){return null;}
		final ArrayList<String> aliases=new ArrayList<String>();
		for(byte[][] names : Data.scaffoldNames){
			if(names==null){continue;}
			for(byte[] name : names){
				if(name!=null && Tools.indexOfWhitespace(name)>=0){
					aliases.add(new String(name, StandardCharsets.US_ASCII));
				}
			}
		}
		return aliases.isEmpty() ? null : aliases.toArray(new String[0]);
	}
	
	/** Waits for a converter or the header-failure flag; printed interruptions are retried.
	 * Throws if no converter is available afterward, otherwise returns a shallow clone.
	 * Clones share the reference map built during initialization, not independent dictionaries.
	 * @return Clone of the published converter
	 */
	private SamToBamConverter getConverter(){
		synchronized(this){
			//[stream/BamWriter#001 FIXED] workers park here until writeHeader() builds sharedConverter+notifyAll; the
			//`&& !headerFailed` exit + the throw below mean a failed header wakes them to a LOUD abort instead of a permanent hang.
			while(sharedConverter==null && !headerFailed){
				try{this.wait();}
				catch(InterruptedException e){e.printStackTrace();}
			}
			if(sharedConverter==null){//initialization failed or an external finishError() woke this worker before conversion setup
				throw new RuntimeException("BAM writer initialization did not complete; aborting worker thread.");
			}
			return (SamToBamConverter)sharedConverter.clone();
		}
	}
	
	/** Converts list entries and mates, including configured secondary alignments.
	 * @param reads Nonnull list; null entries are skipped
	 * @return New list of generated or reused SAM records
	 */
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads){
		return toSamLines(reads, true, true);
	}

	/** Builds SAM records and applies positional entry/mate selection.
	 * May reuse attached SamLines. Before selection, KEEP_NAMES=false may normalize a
	 * mate qname to the list entry's qname, including on reused objects.
	 * @param reads Nonnull list with optional attached mates; null entries are skipped
	 * @param writeR1 Include list-entry records, not a filter on stored pair-number bits
	 * @param writeR2 Include attached-mate records
	 * @return New list containing generated and/or borrowed records and configured secondaries
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

	/** Appends a primary record and configured secondary records generated from later sites.
	 * Uses a Read clone for secondary conversion; ignores a null Read or primary.
	 * @param r Read supplying mapping sites
	 * @param primary Generated or attached primary record
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
	
	/** Selects shared-header mode, otherwise a provided header, otherwise generated lines.
	 * Shared lookup may wait; an unavailable selected header becomes an empty list.
	 * Unlike SamWriter, this method does not fall back from a null shared result (#002).
	 * @return Nonnull selected header list without defensive copying
	 */
	ArrayList<byte[]> getHeader(){
		if(verbose){System.err.println("Fetching header: "+useSharedHeader+","+(header!=null));}
		ArrayList<byte[]> headerLines;
		//TODO: Probable bug #002 - a null shared header skips the provided/generated alternatives
		//and becomes empty below, even if a supplied reference header exists. Expected fallback
		//policy and runtime consequences need separate review.
		if(useSharedHeader){
			headerLines=SamReadInputStream.getSharedHeader(true);
		}else if(header!=null){
			headerLines=header;
		}else{
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

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input job retaining one Read or SamLine wrapper for NORMAL jobs, or a marker. */
	static class SamWriterInputJob implements HasID{

		/** Retains wrappers; enabled assertions check normal exclusivity, marker type and IDs. */
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

		/** Reports the last marker type. */
		@Override
		public boolean last(){return type==ListNum.LAST;}

		/** Creates a payload-free poison marker with the requested ID. */
		@Override
		public SamWriterInputJob makePoison(long id_){
			return new SamWriterInputJob(null, null, ListNum.POISON, id_);
		}

		/** Creates a payload-free last marker with the requested ID. */
		@Override
		public SamWriterInputJob makeLast(long id_){
			return new SamWriterInputJob(null, null, ListNum.LAST, id_);
		}

		/** Optional borrowed Read wrapper. */
		public final ListNum<Read> reads;
		/** Optional borrowed SamLine wrapper. */
		public final ListNum<SamLine> lines;
		/** NORMAL, LAST, POISON or PROTO. */
		public final int type;
		/** Wrapper ID or supplied marker ID. */
		public final long id;
	}

	/** Binary-record bytes with an ordering ID, or a payload-free marker. */
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

		/** Reports the normal end-of-output marker type. */
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

		/** Input job ID or supplied terminal ID. */
		public final long id;
		/** Owned binary-record snapshot, null for markers. */
		public final byte[] bytes;
		/** NORMAL, LAST, POISON or PROTO. */
		public final int type;
	}

	/** Processing thread - converts reads/lines to BAM binary or writes output. */
	private class ProcessThread extends Thread{

		/** Names this thread for output role zero or a positive worker role.
		 * @param tid_ Role ID assigned by spawnThreads
		 */
		ProcessThread(final int tid_){
			tid=tid_;
			setName(tid==0 ? "BamWriter-Output" : "BamWriter-Worker-"+tid);
			if(verbose){System.err.println("tid "+tid+" created.");}
		}

		/** Executes the selected role, retaining historical error handling on a throwable.
		 * If the role returns, restores the -1 counter offsets and marks success.
		 */
		@Override
		public void run(){
			try{
				if(tid==0){
					writeOutput(); //Writer thread
					if(verbose){System.err.println("Consumer "+tid+" finished.");}
				}else{
					processJobs(); //Worker thread
					if(verbose){System.err.println("Worker "+tid+" finished.");}
				}
			}catch(Throwable t){
				//[stream/BamWriter#001 FIXED 2026-06-20 (greenlit + verified)]: ANY thread failure here used to die WITHOUT
				//poisoning/finishing the OQS2 -> a writer death left outq undrained (workers block in addOutput) and `finished`
				//unset (main HANGS in waitForFinish); a worker death left its ordered job id undelivered (the writer's getOutput
				//blocks forever on the gap). Now escalate to a LOUD finish: setErrorState + oqs.setFinished(true) [force=true
				//poisons inq+outq, waking the blocked writer/workers AND main], then rethrow so it crashes loud instead of hanging.
				setErrorState(true);
				oqs.setFinished(true);
				throw new RuntimeException("BamWriter thread "+tid+" failed", t);
			}

			synchronized(this){
				readsWrittenT++;
				basesWrittenT++;
				success=true;
			}
		}

		/** Writes the header and queued bytes, then joins other threads and merges their counters.
		 * ThreadWaiter skips the current thread. Accumulation includes unsuccessful workers and
		 * records their error state. Normal cleanup precedes non-forced queue completion.
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
					//[stream/BamWriter#001 FIXED] write error (disk-full/broken-pipe) -> rethrow; run()'s catch escalates to
					//oqs.setFinished(true), so workers+main wake to a loud finish instead of hanging on the undrained outq.
					throw new RuntimeException("Error writing output", e);
				}

				job=oqs.getOutput();
			}

			//Wait for other threads and accumulate statistics
			ThreadWaiter.waitForThreadsToFinish(alpt);
			synchronized(BamWriter.this){
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
			//Close output stream and signal completion
			boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			errorState|=b;
			if(verbose){System.err.println("Consumer finished writing.");}
			oqs.setFinished(false);
			if(verbose){System.err.println("Consumer set oqs finished.");}
		}

		/** Obtains a converter clone and serializes input batches into copied binary-record bytes.
		 * Direct SamLine jobs bypass mate selection/secondary expansion. Counts records and
		 * SamLine lengths before submitting output, then clears the reused builder. Returns a
		 * received poison job to the input queue for remaining workers.
		 */
		void processJobs(){
			final ByteBuilder bb=new ByteBuilder(65536);
			final SamToBamConverter converter=getConverter();

			SamWriterInputJob job=oqs.getInput();
			while(job!=null && !job.poison()){
				assert(bb.length()==0 && bb.array.length>=65536);
				//Convert to SamLines if needed
				ArrayList<SamLine> lines;
				if(job.lines!=null){
					lines=job.lines.list;
				}else{
					lines=toSamLines(job.reads.list, writeR1, writeR2);
				}

				//Format SamLines to BAM bytes and count
				for(SamLine sl : lines){
					if(sl==null){continue;}
					converter.appendAlignment(sl, bb);
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

		/** SAM-record counter with a -1 offset restored when run's role returns. */
		//Normal accounting: -1 plus the final run() increment yields the formatted record count;
		//the base counter uses the same offset around SamLine.length sums. Aggregation includes
		//every non-output thread and separately records !success, so failed workers are not
		//automatically excluded. This offset rationale does not establish complete failure-path totals.
		protected long readsWrittenT=-1;
		/** SamLine-length sum with the same normal-completion offset. */
		protected long basesWrittenT=-1;
		/** Set after the selected role returns and counter offsets are restored. */
		boolean success=false;
		/** Role ID: zero for output, positive for conversion workers. */
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Getters and Setters      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Accumulates an error observation without clearing earlier errors.
	 * @param b Error observation to combine
	 */
	synchronized void setErrorState(boolean b){
		errorState|=b;
	}
	
	/** Returns the current local error flag under this writer's monitor. */
	@Override
	public synchronized boolean errorState(){return errorState;}
	
	/** Reads local error and queue-completion flags without joining or finalizing threads. */
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
	/** Conversion-worker count, excluding the output thread and compression backend resources. */
	final int threads;
	/** Selects shared-header mode instead of provided/generated headers. */
	final boolean useSharedHeader;
	/** Captured all-header suppression, including append to an existing output. */
	final boolean supressHeader;
	/** Suppresses generated sequence lines and the helper's emitted binary reference table. */
	final boolean supressHeaderSequences;
	/** Set after converter setup and delegated header output, even when bytes were suppressed. */
	boolean headerWritten=false;
	/** Optional borrowed header list and byte arrays. */
	final ArrayList<byte[]> header;
	/** Queue coordinator; current JobQueue requires dense IDs starting at zero. */
	final OrderedQueueSystem2<SamWriterInputJob, SamWriterOutputJob> oqs;
	/** Stream opened through the configured BGZF output path during construction. */
	final OutputStream outstream;
	
	/** Output and worker threads created by the latest start call. */
	private ArrayList<ProcessThread> alpt;
	/** Published converter; shallow worker clones share its reference-name mapping. */
	private SamToBamConverter sharedConverter;
	
	/** Public accumulator of worker-record totals merged during normal finalization. */
	public long readsWritten=0;
	/** Public accumulator of worker SamLine-length totals merged during normal finalization. */
	public long basesWritten=0;
	/** Local accumulated error flag. */
	private boolean errorState=false;
	/** Select source-list entries during Read conversion; not a pair-bit filter. */
	private final boolean writeR1;
	/** Select attached mates during Read conversion; direct SamLine jobs bypass selection. */
	private final boolean writeR2;
	/** Header-failure/external-stop condition observed by converter waiters; see historical #001. */
	private volatile boolean headerFailed=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Default requested conversion-worker count before clamping to Shared.threads. */
	public static int DEFAULT_THREADS=7; // BAM benefits from more threads

	/** Shared diagnostic switch. */
	public static final boolean verbose=false;
	
	/** Retained diagnostic stream; current implementation prints directly to System.err. */
	protected PrintStream outstream2=System.err;

}
