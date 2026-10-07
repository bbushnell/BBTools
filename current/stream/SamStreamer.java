package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * SAM reader with one input thread and one or more parsing workers coordinated by
 * OrderedQueueSystem. ByteFile may add background threads. Construction does not
 * start reading; configure sampling before one start and consume nextLines or
 * nextReads. Shared header publication affects other readers in the JVM.
 * Global ByteFile force flags take precedence over this reader's local BF2 preference.
 * Input limits count nonheader records before sampling; discarded records are not
 * parsed. Output batches can be empty after sampling. Ordering follows the queue
 * implementation, not a guarantee that an ordered=false request disables ordering.
 * Terminal consumption does not itself wait for final counter aggregation.
 *
 * @author Brian Bushnell
 * @contributor Isla, Shinobu
 * @date November 4, 2016
 */
public class SamStreamer implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves a SAM-fallback descriptor and configures an unstarted reader.
	 *
	 * @param fname_ Input path; subprocess input is permitted
	 * @param threads_ Parsing workers, clamped to 1..Shared.threads(); values below one use DEFAULT_THREADS
	 * @param saveHeader_ Request shared-header publication
	 * @param ordered_ Ordering request passed to the queue system
	 * @param maxReads_ Maximum nonheader input records before sampling; negative means unlimited
	 * @param makeReads_ True to convert SamLine entries to Read objects
	 */
	public SamStreamer(final String fname_, final int threads_, final boolean saveHeader_, final boolean ordered_,
			final long maxReads_, final boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_,
			saveHeader_, ordered_, maxReads_, makeReads_);
	}

	/**
	 * Captures input configuration and creates the queues without starting readers.
	 * Header retention disabled here leaves any earlier shared publication untouched.
	 *
	 * @param ffin_ Nonnull input descriptor passed to ByteFile
	 * @param threads_ Parsing workers, clamped to 1..Shared.threads(); values below one use DEFAULT_THREADS
	 * @param saveHeader_ Request shared-header publication
	 * @param ordered_ Ordering request passed to the queue system
	 * @param maxReads_ Maximum nonheader input records before sampling; negative means unlimited
	 * @param makeReads_ True to convert SamLine entries to Read objects
	 */
	public SamStreamer(final FileFormat ffin_, final int threads_, final boolean saveHeader_, final boolean ordered_,
			final long maxReads_, final boolean makeReads_){
		fname=ffin_.name();
		ffin=ffin_;
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		saveHeader=saveHeader_;
		header=(saveHeader ? new ArrayList<byte[]>() : null);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		makeReads=makeReads_;

		// Create OQS with prototypes for LAST/POISON generation
		ListNum<byte[]> inputPrototype=new ListNum<byte[]>(null, 0, ListNum.PROTO);
		ListNum<SamLine> outputPrototype=new ListNum<SamLine>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<ListNum<byte[]>, ListNum<SamLine>>(
			threads, ordered_, inputPrototype, outputPrototype);

	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resets reported counters and starts the input/parsing threads. Call once;
	 * queues, header, error and sampling state are not reset for reuse.
	 */
	@Override
	public void start(){
		if(verbose){outstream.println("SamStreamer.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Process the reads in separate threads
		spawnThreads();

		if(verbose){outstream.println("Started.");}
	}

	/** Marks the output queue finished using its forced terminal mechanism.
	 * Does not join workers or close the input thread's local ByteFile. This is not
	 * a guarantee that all background work has stopped when the method returns.
	 */
	@Override
	public synchronized void close(){
		//[stream/SamStreamer#002 partial fix 2026-09-05] was a no-op; see BamStreamer.close for the rationale.
		//The input thread free-runs the remaining file after this and exits via its own poison protocol.
		oqs.setFinished(true);
	}

	/** @return Input filename being streamed. */
	@Override
	public String fname(){return fname;}

	/** @return Delegated queue terminal state; not an occupancy or worker-completion check */
	@Override
	public boolean hasMore(){return oqs.hasMore();}

	/** @return false; this reader does not link records into paired Read batches */
	@Override
	public boolean paired(){return false;}

	/** @return 0 because SamStreamer does not track pair numbers. */
	@Override
	public int pairnum(){return 0;}

	/** @return Currently observed retained-record total, aggregated by the input thread after waiting for workers */
	@Override
	public long readsProcessed(){return readsProcessed;}

	/** @return Currently observed retained sequence-base total; terminal consumption is not an aggregation barrier */
	@Override
	public long basesProcessed(){return basesProcessed;}

	/**
	 * Configures sampling before start. Workers share the outer generator, so a seed
	 * alone does not guarantee the same retained records across worker schedules.
	 * Sampling precedes parsing; discarded records are not validated.
	 * @param rate Keep threshold against nextFloat; values at least one disable sampling
	 * @param seed Seed for the shared generator created by Shared.threadLocalRandom
	 */
	@Override
	public void setSampleRate(final float rate, final long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	/** @return Next list of reads, delegating to {@link #nextReads()}. */
	@Override
	public ListNum<Read> nextList(){return nextReads();}

	/** Converts a nextLines batch to new containers reusing its attached Read objects.
	 * Requires construction with makeReads=true. Sampling may leave an empty batch.
	 * Read numeric IDs are the input chunk's firstRecordNum plus retained ordinal,
	 * not necessarily each original record's position after sampling compaction.
	 * @return New Read batch with the parsed batch's ID, or null when nextLines returns null
	 */
	public ListNum<Read> nextReads(){
		assert(makeReads);
		ListNum<SamLine> lines=nextLines();
		if(lines==null){return null;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		if(!lines.isEmpty()){
			for(SamLine line : lines){
				assert(line.obj!=null);
				reads.add((Read)line.obj);
			}
		}
		ListNum<Read> ln=new ListNum<Read>(reads, lines.id);
		return ln;
	}

	/** Takes parsed output through the queue, waiting as required by that implementation.
	 * A LAST batch marks output finished. Null/LAST checks the currently observed error
	 * flag and aborts the JVM if it is set; this call does not join or await counter
	 * aggregation. Nonterminal empty batches remain distinct from end of input.
	 * @return Parsed batch, possibly empty, or null on a terminal queue result
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		ListNum<SamLine> list=oqs.getOutput();
		if(verbose){
			if(list==null){outstream.println("Consumer got null.");}else{outstream.println("Consumer got list "+list.id()+" type "+list.type);}
		}
		if(list==null || list.last()){
			if(list!=null && list.last()){
				oqs.setFinished(true);
			}
			//[stream/SamStreamer#001 FIXED] crash LOUD on a producer/worker error instead of silently truncating: a thread
			//death set errorState + force-poisoned the OQS (so we got null/last here). A bare `return null` would look like
			//clean EOF -> wrong/partial results. KillSwitch.kill exits loudly (BBTools contract: crash, never silently wrong).
			if(errorState){KillSwitch.kill("Error reading SAM file (corrupt or truncated): "+fname);}
			return null;
		}
		return list;
	}

	/** @return Currently observed input/worker error flag, without waiting for completion */
	@Override
	public boolean errorState(){return errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Starts one input thread followed by the configured number of parsing workers.
	 * The input thread retains the thread list for later waiting and aggregation.
	 */
	void spawnThreads(){
		//Determine how many threads may be used
		final int threads=this.threads+1;

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt));
		}
		if(verbose){outstream.println("Spawned threads.");}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Worker that either pulls byte chunks from disk (tid 0) or converts them to SamLine/Read objects (other tids).
	 */
	private class ProcessThread extends Thread{

		/** Captures a thread role without starting it.
		 * @param tid_ Zero for input, positive for a parsing worker
		 * @param alpt_ Thread list retained only by the input thread
		 */
		ProcessThread(final int tid_, final ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("SamStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}

		/** Runs the selected role, marking success only when that role returns normally. */
		@Override
		public void run(){
			//Process the reads
			if(tid==0){
				processInputThread();
			}else{
				makeReads();
			}

			//Indicate successful exit status
			success=true;
			if(verbose){outstream.println("tid "+tid+" terminated.");}
		}

		/**
		 * Reads input blocks, records thrown input failures and calls queue poison in finally.
		 * Then waits for the thread list and aggregates each other thread's counts/success.
		 * Output terminal consumption alone does not wait for this aggregation.
		 */
		void processInputThread(){
			//[stream/SamStreamer#001 FIXED 2026-06-20]: the try/finally GUARANTEES oqs.poison() runs even when processBytes() throws
			//(empty-line AIOOBE at L281/L344, or any read error), so the workers in getInput() + the consumer in getOutput() wake
			//instead of hanging forever; the catch records errorState so nextLines crashes LOUD via KillSwitch instead of silently
			//truncating. Input-thread death leaves NO ordered gap (chunks 0..k all delivered, LAST reachable) so plain poison() suffices
			//here (the worker side needs setFinished(true), see makeReads). Mirrors the greenlit BamStreamer#001.
			try{
				processBytes();
			}catch(Throwable t){
				errorState=true;
				outstream.println("SamStreamer: error reading "+fname+": "+t);
			}finally{
				oqs.poison();//ALWAYS poison so workers (getInput) + consumer (getOutput) wake
			}
			if(verbose){outstream.println("tid "+tid+" done with processBytes + poisoning.");}

			//Wait for completion of all threads
			boolean allSuccess=true;
			ThreadWaiter.waitForThreadsToFinish(alpt);
			for(ProcessThread pt : alpt){
				//Wait until this thread has terminated
				if(pt!=this){
					//Accumulate per-thread statistics
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					allSuccess&=pt.success;
				}
			}
			if(verbose){outstream.println("tid "+tid+" noted all process threads finished.");}

			//Track whether any threads failed
			if(!allSuccess){errorState=true;}
			if(verbose){outstream.println("tid "+tid+" finished! Error="+errorState);}
		}

		/** Reads SAM lines into byte-array batches for the queue.
		 * Prefers BF2 locally, honoring global force flags and backend allowances
		 * through the ByteFile factory. Nonheader records count toward
		 * maxReads before sampling. Batch thresholds count records and full line bytes
		 * without newlines, and are checked after adding each record.
		 * Requested headers publish before the first record or on normal loop exit.
		 * Normal completion closes the local ByteFile and folds its reported error;
		 * exceptional exits are handled by processInputThread.
		 */
		void processBytes(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

			ByteFile bf=ByteFile.makeByteFileWithPreference(ffin, 2);

			long listNumber=0;
			long reads=0;
			int bytes=0;
			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<byte[]> ln=new ListNum<byte[]>(new ArrayList<byte[]>(slimit), listNumber++);
			ln.firstRecordNum=reads;

			for(byte[] line=bf.nextLine(); line!=null && reads<maxReads; line=bf.nextLine()){
				if(line[0]=='@'){
					if(header!=null){
						if(Shared.TRIM_READ_DESCRIPTION){line=SamReadInputStream.trimHeaderSQ(line);}
						header.add(line);
					}
				}else{
					//comprehension: the FIRST non-@ line flushes the accumulated header to the shared store exactly once (header=null
					//after), so workers see the full header via SamReadInputStream before any record is parsed. The duplicate flush
					//after normal input-loop exit covers header-only files (no records => the in-loop flush never fires).
					if(header!=null){
						SamReadInputStream.setSharedHeader(header);
						header=null;
					}
					reads++;
					bytes+=line.length;
					ln.add(line);
					if(ln.size()>=slimit || bytes>=blimit){
						oqs.addInput(ln);
						ln=new ListNum<byte[]>(new ArrayList<byte[]>(slimit), listNumber++);
						ln.firstRecordNum=reads;
						bytes=0;
					}
				}
			}

			if(header!=null){
				SamReadInputStream.setSharedHeader(header);
				header=null;
			}
			if(verbose){outstream.println("tid "+tid+" ran out of input.");}
			if(ln.size()>0){
				oqs.addInput(ln);
			}
			ln=null;
			if(verbose){outstream.println("tid "+tid+" done reading bytes.");}
			errorState|=bf.close();//Fold the reader's error state (truncated/corrupt input) so it isn't silently dropped at the streamer boundary
			if(verbose){outstream.println("tid "+tid+" closed stream.");}
		}

		/**
		 * Samples and compacts each input batch, then parses retained records.
		 * Preserves batch IDs, including for batches made empty by sampling. Converted
		 * Read IDs start at the original chunk offset and increase for retained records.
		 * On failure, sets the error flag, marks output finished and rethrows; normal
		 * poison consumption reinserts that marker for the other workers.
		 */
		void makeReads(){
			//[stream/SamStreamer#001 FIXED 2026-06-20]: a parse throw (SamLine ctor / line[0] AIOOBE / sl.toRead) used to kill this
			//worker with its ordered job UNDELIVERED -> the ordered consumer blocked forever on the gap (a plain LAST marker sorts
			//AFTER the gap, so it can't release it; verified in JobQueue.take/heapReady). The catch force-poisons outq via
			//oqs.setFinished(true) (sets JobQueue.poisoned -> take() returns null past the gap) + records errorState, so the consumer
			//wakes and crashes LOUD in nextLines. NOT oqs.poison() here: poison() sets lastSeen, which would trip addInput's assert in
			//the still-reading input thread. Rethrow so run() skips success=true (this dead worker reports failure to the input thread).
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}

			final LineParser1 lp=new LineParser1('\t');
			try{
				ListNum<byte[]> list=oqs.getInput();
				while(list!=null && !list.poison()){
					if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}

					// Apply subsampling if needed
					if(samplerate<1f && randy!=null){
						int nulled=0;
						for(int i=0; i<list.size(); i++){
							if(randy.nextFloat()>=samplerate){
								list.list.set(i, null);
								nulled++;
							}
						}
						if(nulled>0){Tools.condenseStrict(list.list);}
					}

					ListNum<SamLine> reads=new ListNum<SamLine>(
						new ArrayList<SamLine>(list.size()), list.id);
					long readID=list.firstRecordNum;
					for(byte[] line : list){
						if(line[0]=='@'){
							//Ignore header lines
						}else{
							SamLine sl=new SamLine(lp.set(line));
							reads.add(sl);
							if(makeReads){
								Read r=sl.toRead(FASTQ.PARSE_CUSTOM);
								sl.obj=r;
								r.samline=sl;
								r.numericID=readID++;
								if(!r.validated()){r.validate(true);}
							}
							readsProcessedT++;
							basesProcessedT+=(sl.seq==null ? 0 : sl.length());
						}
					}
					oqs.addOutput(reads);
					list=oqs.getInput();
				}
				if(verbose){outstream.println("tid "+tid+" done making reads.");}
				//Re-inject poison for other workers
				if(list!=null){oqs.addInput(list);}
			}catch(Throwable t){
				errorState=true;
				oqs.setFinished(true);//force-poison outq -> release the consumer past THIS worker's undelivered-job gap
				throw new RuntimeException("SamStreamer worker "+tid+" failed: "+fname, t);
			}
		}

		/** Retained records parsed by this worker. */
		protected long readsProcessedT=0;
		/** Retained sequence bases parsed by this worker. */
		protected long basesProcessedT=0;
		/** Whether the selected role returned normally; input failures also use the outer error flag. */
		boolean success=false;
		/** Zero for the input role, positive for a parsing worker. */
		final int tid;

		/** Thread list used by the input role for waiting; null for parsing workers. */
		ArrayList<ProcessThread> alpt;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path. */
	public final String fname;

	/** FileFormat descriptor for the input file. */
	final FileFormat ffin;

	/** OrderedQueueSystem coordinating byte input and SamLine output queues. */
	final OrderedQueueSystem<ListNum<byte[]>, ListNum<SamLine>> oqs;

	/** Parsing worker count; spawnThreads adds one input thread. */
	final int threads;
	/** Captured request for shared-header publication. */
	final boolean saveHeader;
	/** Whether to create attached Read objects in addition to parsed SamLine objects. */
	final boolean makeReads;

	/** Header accumulator; null when retention is disabled or after publication. */
	ArrayList<byte[]> header;

	/** Retained record total aggregated by the input thread after worker waiting. */
	protected long readsProcessed=0;
	/** Retained sequence-base total aggregated by the input thread after worker waiting. */
	protected long basesProcessed=0;

	/** Maximum nonheader input records before sampling; negative constructor limits become Long.MAX_VALUE. */
	final long maxReads;

	/** Output stream for status and verbose messages. */
	protected PrintStream outstream=System.err;
	/** Reported input/worker errors; reading the flag does not wait for completion. */
	public boolean errorState=false;
	/** Keep threshold configured before start; values at least one disable sampling. */
	float samplerate=1f;
	/** Shared generator used by parsing workers; no per-worker random stream is allocated here. */
	shared.Random randy=null;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Record threshold captured by the input thread; initialize before start. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Full SAM line-byte threshold excluding newlines; initialize before start. */
	public static int TARGET_LIST_BYTES=shared.Shared.bufferSize();
	/** Default number of worker threads when none is specified. */
	public static int DEFAULT_THREADS=3;
	/** Enables verbose logging when true. */
	public static final boolean verbose=false;

}
