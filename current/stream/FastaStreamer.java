package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * FASTA reader with one input thread and configurable parallel conversion workers.
 * An ordered queue returns batches in input-job order; a ByteFile backend may add
 * its own input threads. Configure once before start, then drain through terminal
 * and close. Restart and arbitrary concurrent close are not supported protocols.
 * Limits apply to input before positional sampling; interleaved limits normally
 * count pairs. Returned batch IDs follow input jobs, including empty sampled batches.
 * Counter getters report asynchronously aggregated worker totals and can lag EOF.
 * Names use US-ASCII decoding. Read validation may apply other configured changes.
 * 
 * @author Brian Bushnell
 * @contributor Shinobu (documentation and lifecycle maintenance)
 * @date November 5, 2025
 */
public class FastaStreamer implements Streamer{

	/** Runs the factory-selected reader and prints counts from returned data to stderr.
	 * The factory can choose another reader class. Reported errors reject normal
	 * statistics; ordinary Java unwinding attempts close, but VM halt bypasses it.
	 * @param args FASTA path followed by an optional default conversion-thread count
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		String fname=args[0];
		if(args.length>1){DEFAULT_THREADS=Integer.parseInt(args[1]);}

		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		Streamer st=StreamerFactory.makeStreamer(ff, 0, true, -1, true, true);
		long reads=0, bases=0;
		//[stream/FastaStreamer#002] Factory-selected ST/2ST readers can report errors at EOF.
		try{
			st.start();
			for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
				for(Read r : ln){
					reads+=r.pairCount();
					bases+=r.pairLength();
				}
			}
		}finally{st.close();}
		if(st.errorState()){
			throw new RuntimeException("FastaStreamer failed while reading input; see preceding diagnostic.");
		}

		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Captures input settings and creates queues without opening a backend.
	 * @param fname_ Input path; descriptor permits subprocess input
	 * @param threads_ Conversion workers; below 1 selects DEFAULT_THREADS
	 * @param pairnum_ 0 for unpaired/R1/interleaved, 1 for separate R2 input
	 * @param maxReads_ Input limit before sampling; negative means unlimited
	 */
	public FastaStreamer(String fname_, int threads_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), threads_, pairnum_, maxReads_);
	}

	/** Captures descriptor, flags and input limit, and creates ordered queues.
	 * Conversion workers are clamped between 1 and Shared.threads(); start adds
	 * one input thread. Direct-constructor values below 1 differ from factory hint 0.
	 * @param ffin_ Nonnull descriptor, including interleaving/amino settings
	 * @param threads_ Conversion workers; below 1 selects DEFAULT_THREADS
	 * @param pairnum_ 0 for interleaved input, otherwise 0 or 1; checked by assertion
	 * @param maxReads_ Input limit before sampling; normally pairs when interleaved,
	 * otherwise individual headers; negative means unlimited
	 */
	public FastaStreamer(FileFormat ffin_, int threads_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		pairnum=pairnum_;
		flag=(ffin.amino() || Shared.AMINO_IN ? Read.AAMASK : 0);
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		// Create OQS with prototypes for LAST/POISON generation
		ListNum<byte[]> inputPrototype=new ListNum<byte[]>(null, 0, ListNum.PROTO);
		ListNum<Read> outputPrototype=new ListNum<Read>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<ListNum<byte[]>, ListNum<Read>>(
			threads, true, inputPrototype, outputPrototype);

		if(verbose){outstream.println("Made FastaStreamer-"+threads);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resets exposed counters and launches input/conversion threads.
	 * Call once on a fresh instance: this does not reset queues or error state for
	 * restart. Keep public flags and global parsing configuration stable while running.
	 */
	@Override
	public void start(){
		if(verbose){outstream.println("FastaStreamer.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Process the reads in separate threads
		spawnThreads();

		if(verbose){outstream.println("FastaStreamer started.");}
	}

	/** Attempts backend cleanup, then forced output finish even if cleanup throws.
	 * This does not join workers or guarantee cancellation/terminal delivery.
	 * Producer completion uses backend-only cleanup and ordinary poison instead.
	 */
	@Override
	public void close(){
		//[stream/FastaStreamer#003] Preserve reported errors and the finish attempt on every close path.
		try{closeInput();}finally{oqs.setFinished(true);}
	}

	/** Returns the input name captured from the descriptor. */
	@Override
	public String fname(){return fname;}

	/** Returns the output queue's availability hint without consuming a batch.
	 * This is not a completion barrier; forced finish can leave the hint true
	 * while nextList returns null. Drive iteration with nextList itself.
	 */
	@Override
	public boolean hasMore(){return oqs.hasMore();}

	/** Returns observed error status; this plain read does not wait for completion. */
	@Override
	public boolean errorState(){return errorState;}

	/** Returns the interleaving setting captured at construction. */
	@Override
	public boolean paired(){return interleaved;}

	/** Returns the input-side pair marker used for unpaired conversion. */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns the currently aggregated worker-read total, which may lag returned data or EOF.
	 * Unpaired workers count retained constructed Reads; paired workers count all
	 * constructed individual Reads before sampling. Failed/unpublished work can count.
	 * The input thread adds worker totals after joins; this getter does not wait for it.
	 */
	@Override
	public synchronized long readsProcessed(){return readsProcessed;}

	/** Returns the corresponding aggregated base total; see readsProcessed(). */
	@Override
	public synchronized long basesProcessed(){return basesProcessed;}

	/** Configures positional sampling before start; keep settings stable while running.
	 * Unpaired sampling uses input-record positions; interleaved sampling uses pair
	 * positions, retaining both mates together. Dropped records still advance IDs.
	 * @param rate Retained fraction, normally in [0,1]; not validated here
	 * @param seed Nonnegative for reproducibility; negative resolves once to a random seed
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		sampleSeed=Streamer.resolveSampleSeed(seed);
	}

	/** Takes the next ordered output batch, which may be empty after sampling.
	 * Batch IDs follow input jobs; firstRecordNum remains unspecified (-1).
	 * Returned lists/Read data are not recycled here. A terminal/null queue result
	 * returns null unless an error was reported, in which case KillSwitch halts the VM.
	 * Terminal observation does not wait for worker joins or counter aggregation.
	 */
	@Override
	public ListNum<Read> nextList(){
		ListNum<Read> list=oqs.getOutput();
		if(verbose){
			if(list==null){outstream.println("Consumer got null.");}else{outstream.println("Consumer got list "+list.id()+" type "+list.type);}
		}
		if(list==null || list.last()){
			if(list!=null && list.last()){
				oqs.setFinished(true);
			}
			//[stream/FastaStreamer#001 FIXED] crash LOUD on a producer/worker error instead of silently truncating: a thread death
			//set errorState + force-poisoned the OQS (so we reach here with null/last). A bare `return null` would look like clean EOF ->
			//wrong/partial results downstream. KillSwitch.kill exits loudly (BBTools contract: crash, never silently wrong).
			if(errorState){KillSwitch.kill("Error reading FASTA file (corrupt or truncated): "+fname);}
			return null;
		}
		return list;
	}

	/** Always throws UnsupportedOperationException; FASTA has no SamLine representation. */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Closes an established backend without abandoning queued output.
	 * Clears the reference only after normal return. Only writes true to the shared
	 * error flag so cleanup cannot overwrite a worker's concurrent failure report.
	 * Throws propagate; this adds no concurrent-close or physical-closure guarantee.
	 */
	private void closeInput(){
		boolean completed=false;
		try{
			if(bf!=null){
				final boolean error=bf.close();
				if(error){errorState=true;}
				bf=null;
			}
			completed=true;
		}finally{if(!completed){errorState=true;}}
	}

	/** Creates and starts all conversion workers plus input thread zero.
	 * Only the input thread retains the group for subsequent joins and aggregation.
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
		for(ProcessThread pt : alpt){pt.start();}
		if(verbose){outstream.println("Started threads.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input coordinator (tid 0) or a converter; each owns its local parse buffers/counts. */
	private class ProcessThread extends Thread{

		/** Creates a named role without starting it.
		 * @param tid_ Zero for input/coordinator, positive for conversion
		 * @param alpt_ Complete group; retained only by the input coordinator
		 */
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("FastaStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}

		/** Dispatches this role and marks a normal return, independently of errorState.
		 * Conversion failure propagates after flagging/force-finish. The input path
		 * can catch a processing failure and still return normally with errorState set.
		 */
		@Override
		public void run(){
			//Process the reads
			synchronized(this){
				if(tid==0){
					processBytes();
				}else{
					if(interleaved){
						makeReadsInterleaved();
					}else{
						makeReadsSingle();
					}
				}
			}

			//Indicate successful exit status
			success=true;
			if(verbose){outstream.println("tid "+tid+" terminated.");}
		}

		/** Reads jobs, attempts backend cleanup and ordinary poison, then joins converters.
		 * Aggregation occurs after poison, so output EOF can precede final counters.
		 * Cleanup/poison failure can skip joins and aggregation; exceptions may replace
		 * earlier failures, and blocking can prevent completion.
		 */
		void processBytes(){
			//[stream/FastaStreamer#001 FIXED 2026-06-21] Input failure must record status and attempt ordinary poison.
			//Input jobs retain their ordered IDs; worker failures need force-finish for their undelivered gaps.
			//The catch records errorState so nextList fails loudly. Cleanup and poison may block or throw;
			//nested finally preserves the poison attempt, not a universal wakeup guarantee.
			try{
				processBytes0();
			}catch(Throwable t){
				errorState=true;
				outstream.println("FastaStreamer: error reading "+fname+": "+t);
			}finally{
				//[stream/FastaStreamer#004] Close even on producer failure, before ordinary queue completion.
				//Do not use public close here: its forced output finish can abandon valid queued data.
				try{closeInput();}finally{oqs.poison();}
			}
			if(verbose){outstream.println("tid "+tid+" done with processBytes0 + poisoning.");}

			//Wait for completion of all threads
			boolean allSuccess=true;
			ThreadWaiter.waitForThreadsToFinish(alpt);
			for(ProcessThread pt : alpt){
				//Wait until this thread has terminated
				if(pt!=this){
					synchronized(pt){
						synchronized(FastaStreamer.this){
							//Accumulate per-thread statistics
							readsProcessed+=pt.readsProcessedT;
							basesProcessed+=pt.basesProcessedT;
							allSuccess&=pt.success;
						}
					}
				}
			}
			if(verbose){outstream.println("tid "+tid+" noted all process threads finished.");}

			//Track whether any threads failed
			if(!allSuccess){errorState=true;}
			if(verbose){outstream.println("tid "+tid+" finished! Error="+errorState);}
		}

		/** Opens the selected ByteFile and groups nonempty physical lines into input jobs.
		 * Snapshots Shared.bufferLen()/bufferData() targets. Sequence-line bytes count
		 * toward the byte target; splitting is checked at headers and paired boundaries.
		 * Targets are not hard caps. Header-led input produces whole-record jobs.
		 * Limits count input headers before sampling. Interleaved finite limits below
		 * Long.MAX_VALUE/2 are doubled; larger limits use the stored value directly.
		 * Loop transitions can read beyond the final included record. The caller owns cleanup.
		 */
		private void processBytes0(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

			bf=ByteFile.makeByteFile(ffin);

			long listNumber=0;
			long totalReads=0;

			int headersInList=0;
			int bytesInList=0;

			final int slimit=Shared.bufferLen();
			final int blimit=Shared.bufferData();
			ListNum<byte[]> ln=new ListNum<byte[]>(new ArrayList<byte[]>(), listNumber++);
			ln.firstRecordNum=totalReads;
			final long limit=maxReads*(interleaved && maxReads<Long.MAX_VALUE/2 ? 2 : 1);
			
			for(byte[] line=bf.nextLine(); line!=null && totalReads<=limit; line=bf.nextLine()){
				if(line.length>0){
					if(line[0]!='>'){
						ln.add(line);
						bytesInList+=line.length;
					}else{
						//Found a header.
						if((headersInList>=slimit || bytesInList>=blimit) && 
							(!interleaved || ((headersInList&1)==0))){
							oqs.addInput(ln);
							ln=new ListNum<byte[]>(new ArrayList<byte[]>(), listNumber++);
							ln.firstRecordNum=totalReads;
							headersInList=0;
							bytesInList=0;
						}
						if(totalReads<limit){ln.add(line);}
						headersInList++;
						totalReads++;
					}
				}
			}
			if(verbose){outstream.println("tid "+tid+" ran out of input.");}
			if(ln.size()>0){
				oqs.addInput(ln);
			}
			ln=null;
			if(verbose){outstream.println("tid "+tid+" done reading bytes.");}
		}

		/** Converts unpaired jobs, sampling by input position before Read construction.
		 * Selected Reads are explicitly validated regardless of VALIDATE_IN_CONSTRUCTOR.
		 * IDs advance over discarded records; counters include retained constructed Reads.
		 * A headerless job throws. Empty sampled output still retains its ordered job ID.
		 * On conversion failure, records error status, force-finishes output and rethrows.
		 */
		void makeReadsSingle(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}

			//[stream/FastaStreamer#001 FIXED 2026-06-21]: worker death (the live "No header for record" throw below, or a throw in
			//Read.validate) used to leave its ordered job UNDELIVERED -> the ordered consumer blocked forever on the gap (a plain LAST
			//marker sorts AFTER the gap, so it can't release it; verified in JobQueue.take/heapReady). The catch force-poisons outq via
			//oqs.setFinished(true) (sets JobQueue.poisoned -> take() returns null past the gap) + records errorState, so the consumer wakes
			//and crashes LOUD in nextList. NOT oqs.poison() here: poison() sets lastSeen, tripping addInput's assert in the still-reading
			//input thread. Rethrow so run() skips success=true. Mirrors the greenlit SamStreamer/FastqStreamer workers.
			final ByteBuilder bb=new ByteBuilder(4096);
			try{
				ListNum<byte[]> list=oqs.getInput();
				while(list!=null && !list.poison()){
					if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}

					ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>(50), list.id());
					long readID=list.firstRecordNum;

					// Parse lines into reads using ByteBuilder
					byte[] header=null;

					for(byte[] line : list){
						if(line.length>0 && line[0]=='>'){
							// Save previous record if exists
							if(header!=null){
								//Positional sampling (Streamer.sampleKeep): thread-safe + reproducible; the shared
								//PRNG raced across workers. readID now advances for DROPPED records too, so kept
								//reads carry their true file position (twin-file mates then share numericID).
								if(samplerate>=1f || Streamer.sampleKeep(readID, sampleSeed, samplerate)){
									Read r=new Read(bb.toBytes(), null,
										new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, StandardCharsets.US_ASCII), readID, flag, true);
									r.setPairnum(pairnum);
									if(!r.validated()){r.validate(true);}
									reads.add(r);
									readsProcessedT++;
									basesProcessedT+=r.length();
								}
								readID++;
							}
							header=line;
							bb.clear();
						}else{
							bb.append(line);
						}
					}

					// Save final record
					if(header!=null){
						if(samplerate>=1f || Streamer.sampleKeep(readID, sampleSeed, samplerate)){
							Read r=new Read(bb.toBytes(), null,
								new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, StandardCharsets.US_ASCII), readID, flag, true);
							r.setPairnum(pairnum);
							if(!r.validated()){r.validate(true);}
							reads.add(r);
							readsProcessedT++;
							basesProcessedT+=r.length();
						}
					}else{
						throw new RuntimeException("No header for record "+readID+
							" length "+bb.length()+" in "+fname);
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
				throw new RuntimeException("FastaStreamer worker "+tid+" failed: "+fname, t);
			}
		}

		/** Constructs and validates all job records, then pairs and samples by pair position.
		 * Counters include constructed individual Reads before sampling, unlike the
		 * unpaired path. Mates share positional numeric IDs; output contains first mates.
		 * Odd jobs report an error and only pair the even prefix. Empty sampled output
		 * still occupies its job ID. Conversion failure flags, force-finishes and rethrows.
		 */
		void makeReadsInterleaved(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}

			//[stream/FastaStreamer#001 FIXED 2026-06-21]: same worker-death gap fix as makeReadsSingle -- a throw in Read.validate (or the
			//odd-pair path) leaves an undelivered ordered job -> consumer hang; the catch force-poisons outq via setFinished(true) + records
			//errorState so the consumer wakes and crashes LOUD in nextList. The odd-pair case is softened to errorState (below) so it no longer
			//relies on an -ea-only worker assert. Mirrors the greenlit SamStreamer/FastqStreamer workers.
			final ByteBuilder bb=new ByteBuilder(4096);
			try{
				ListNum<byte[]> list=oqs.getInput();
				while(list!=null && !list.poison()){
					if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}

					ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>(), list.id());
					long readID=list.firstRecordNum/2;

					// Parse lines into reads using ByteBuilder
					//STR-072/IO021: Honor trd consistently with legacy FASTA (Brian, 2026-09-29).
					//CutGff2 also trims locally after its 2026-09-09 lookup-miss report; Read.FIX_HEADER is separate.
					ArrayList<Read> allReads=new ArrayList<Read>();
					byte[] header=null;

					for(byte[] line : list){
						if(line.length>0 && line[0]=='>'){
							// Save previous record if exists
							if(header!=null){
								Read r=new Read(bb.toBytes(), null, new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1,
									StandardCharsets.US_ASCII), 0, flag, true);
								if(!r.validated()){r.validate(true);}
								allReads.add(r);
								readsProcessedT++;
								basesProcessedT+=r.length();
							}
							header=line;
							bb.clear();
						}else{
							bb.append(line);
						}
					}

					// Save final record
					if(header!=null){
						Read r=new Read(bb.toBytes(), null, new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1,
							StandardCharsets.US_ASCII), 0, flag, true);
						if(!r.validated()){r.validate(true);}
						allReads.add(r);
						readsProcessedT++;
						basesProcessedT+=r.length();
					}

					// Pair them up -- an odd count means a truncated/malformed interleaved FASTA (legit-but-bad input); record
					// errorState (-> nextList crashes LOUD) and bound the loop to even, instead of an -ea-only worker assert. Matches FastqStreamer.
					if((allReads.size()&1)!=0){
						System.err.println("Incomplete pair near pairnum "+readID+" in "+fname);
						errorState=true;
					}
					final int lim=(allReads.size()|1)-1;
					//Positional sampling by pair index; the pair is kept or dropped as a unit, and readID
					//advances for dropped pairs too so ids reflect true file position.
					for(int i=0; i<lim; i+=2){
						if(samplerate>=1f || Streamer.sampleKeep(readID, sampleSeed, samplerate)){
							Read r1=allReads.get(i);
							Read r2=allReads.get(i+1);
							r1.setPairnum(0);
							r2.setPairnum(1);
							r1.mate=r2;
							r2.mate=r1;
							reads.add(r1);
							r1.numericID=readID;
							r2.numericID=readID;
						}
						readID++;
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
				throw new RuntimeException("FastaStreamer worker "+tid+" failed: "+fname, t);
			}
		}

		/** Constructed Reads: retained only in unpaired mode, before sampling in paired mode. */
		protected long readsProcessedT=0;
		/** Total lengths of this worker's counted Reads, including unpublished work on failure. */
		protected long basesProcessedT=0;
		/** True when run's processing call returns normally; reported errorState may still be true. */
		boolean success=false;
		/** Zero for input/coordinator, positive for a converter. */
		final int tid;

		/** Complete group for coordinator joins/aggregation; null for converters. */
		ArrayList<ProcessThread> alpt;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;

	/** Primary input file */
	final FileFormat ffin;
	
	/** Backend opened by the input thread and shared with lifecycle cleanup. */
	public ByteFile bf;//TODO: Consider restricting public visibility; cleanup currently needs the shared reference.

	/** Ordered conversion jobs/results and ordinary versus forced termination protocol. */
	final OrderedQueueSystem<ListNum<byte[]>, ListNum<Read>> oqs;

	/** Conversion-worker count; excludes the input thread and any backend workers. */
	final int threads;
	/** Input-side marker for unpaired conversion; interleaved input requires zero. */
	final int pairnum;
	/** Captured pairing mode for adjacent records. */
	final boolean interleaved;

	/** Aggregate of converter read totals; populated asynchronously after worker joins. */
	protected long readsProcessed=0;
	/** Aggregate of converter base totals, potentially incomplete at output EOF. */
	protected long basesProcessed=0;

	/** Input limit before sampling; negative constructor arguments normalize to Long.MAX_VALUE. */
	final long maxReads;
	/** Read-constructor flags; configure before start and keep stable while converting. */
	public int flag;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

//	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
//	public static int TARGET_LIST_BYTES=262144;
	/** Default conversion-worker request used by constructors and FASTA factory selection. */
	public static int DEFAULT_THREADS=3;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** Reported input, conversion or cleanup failure; getters do not wait for publication. */
	public boolean errorState=false;
	/** Positional retention threshold, configured before startup. */
	private float samplerate=1f;
	/** Seed for positional sampling (Streamer.sampleKeep); resolved from setSampleRate's seed */
	private long sampleSeed=17;

}
