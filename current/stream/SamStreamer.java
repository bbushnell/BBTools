package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;
import shared.Tools;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Multithreaded SAM reader built on OrderedQueueSystem.
 * Creates worker threads that read SAM records in order, optionally materializing
 * Read objects, and tracks processed read/base counts.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 4, 2016
 */
public class SamStreamer implements Streamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Builds a streamer for an input filename, creating a FileFormat internally.
	 * Configures thread count, ordering, optional header retention, maximum reads to
	 * consume, and whether to construct Read objects alongside SamLine records.
	 *
	 * @param fname_ Input SAM filename
	 * @param threads_ Desired worker threads (clamped to available cores)
	 * @param saveHeader_ True to retain SAM header lines
	 * @param ordered_ True to preserve read order through the pipeline
	 * @param maxReads_ Maximum reads to process (-1 for unlimited)
	 * @param makeReads_ True to convert SamLine entries to Read objects
	 */
	public SamStreamer(String fname_, int threads_, boolean saveHeader_, boolean ordered_, 
			long maxReads_, boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_, 
			saveHeader_, ordered_, maxReads_, makeReads_);
	}
	
	/**
	 * Builds a streamer from a preconstructed FileFormat.
	 * Initializes header retention, maximum reads, and the ordered queue system used
	 * to move byte blocks and parsed lines between threads.
	 *
	 * @param ffin_ Input file description
	 * @param threads_ Desired worker threads (clamped to available cores)
	 * @param saveHeader_ True to retain SAM header lines
	 * @param ordered_ True to preserve read order through the pipeline
	 * @param maxReads_ Maximum reads to process (-1 for unlimited)
	 * @param makeReads_ True to convert SamLine entries to Read objects
	 */
	public SamStreamer(FileFormat ffin_, int threads_, boolean saveHeader_, boolean ordered_, 
			long maxReads_, boolean makeReads_){
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
		
//		if(verbose || true){outstream.println("Made SamStreamer-"+threads);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resets counters and spawns processing threads to begin streaming. */
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

	/** Closes the streamer. Currently a placeholder for resource cleanup. */
	@Override
	public synchronized void close(){
		//TODO: Unimplemented
	}
	
	/** @return Input filename being streamed. */
	@Override
	public String fname() {return fname;}
	
	/** @return True while the ordered queue still holds output. */
	@Override
	public boolean hasMore() {return oqs.hasMore();}
	
	/** SamStreamer produces unpaired reads; always returns false. */
	@Override
	public boolean paired(){return false;}

	/** @return 0 because SamStreamer does not track pair numbers. */
	@Override
	public int pairnum(){return 0;}
	
	/** @return Number of reads processed across all threads. */
	@Override
	public long readsProcessed() {return readsProcessed;}
	
	/** @return Number of bases processed across all threads. */
	@Override
	public long basesProcessed() {return basesProcessed;}
	
	/**
	 * Configures subsampling of reads during streaming.
	 * @param rate Fraction of reads to keep (1 for all)
	 * @param seed Seed for the per-thread random sampler
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	/** @return Next list of reads, delegating to {@link #nextReads()}. */
	@Override
	public ListNum<Read> nextList(){return nextReads();}
	
	/** Retrieves the next batch of parsed reads from the output queue.
	 * Blocks until data are available or streaming is complete. */
	public ListNum<Read> nextReads(){
		assert(makeReads);
		ListNum<SamLine> lines=nextLines();
		if(lines==null){return null;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		if(!lines.isEmpty()) {
			for(SamLine line : lines){
				assert(line.obj!=null);
				reads.add((Read)line.obj);
			}
		}
		ListNum<Read> ln=new ListNum<Read>(reads, lines.id);
		return ln;
	}
	
	/** Retrieves the next batch of parsed SamLine objects from the output queue.
	 * Handles termination markers and marks the queue finished when the last list is seen. */
	@Override
	public ListNum<SamLine> nextLines(){
		ListNum<SamLine> list=oqs.getOutput();
		if(verbose){
			if(list==null) {outstream.println("Consumer got null.");}
			else {outstream.println("Consumer got list "+list.id()+" type "+list.type);}
		}
		if(list==null || list.last()){
			if(list!=null && list.last()){
				oqs.setFinished(true);
			}
			return null;
		}
		return list;
	}
	
	/** @return True if any worker thread reported an error. */
	@Override
	public boolean errorState() {return errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates and starts the producer/consumer threads that perform streaming work.
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
	private class ProcessThread extends Thread {
		
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("SamStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}
		
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
		 * Drives the input side: reads byte slices, queues them, and aggregates thread statistics.
		 */
		void processInputThread(){
			processBytes();
			if(verbose){outstream.println("tid "+tid+" done with processBytes.");}
			
			// Signal completion via OQS
			oqs.poison();
			if(verbose){outstream.println("tid "+tid+" done poisoning.");}
			
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
		
		/** Reads SAM lines into ListNum<byte[]> chunks, preserving headers when requested,
		 * and feeds them into the ordered queue system. */
		void processBytes(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

			ByteFile.FORCE_MODE_BF2=true;
			ByteFile bf=ByteFile.makeByteFile(ffin);
			
			long listNumber=0;
			long reads=0;
			int bytes=0;
			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<byte[]> ln=new ListNum<byte[]>(new ArrayList<byte[]>(slimit), listNumber++);
			ln.firstRecordNum=reads;
			
			for(byte[] line=bf.nextLine(); line!=null && reads<maxReads; line=bf.nextLine()){
				if(line[0]=='@'){
					if(header!=null) { 
						if(Shared.TRIM_RNAME){line=SamReadInputStream.trimHeaderSQ(line);}
						header.add(line);
					}
				}else{
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
			bf.close();
			if(verbose){outstream.println("tid "+tid+" closed stream.");}
		}
		
		/**
		 * Parses queued byte slices into SamLine objects and optionally materializes Read objects.
		 */
		void makeReads(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}
			
			final LineParser1 lp=new LineParser1('\t');
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
					if(nulled>0) {Tools.condenseStrict(list.list);}
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
			if(list!=null) {oqs.addInput(list);}
		}

		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		boolean success=false;
		final int tid;
		
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
	
	/** Number of worker threads (plus one input thread). */
	final int threads;
	/** True to collect SAM headers before record processing. */
	final boolean saveHeader;
	/** True to convert SamLines to Read objects instead of leaving raw lines. */
	final boolean makeReads;
	
	/** Stored SAM header lines; null when headers are not being retained. */
	ArrayList<byte[]> header;
	
	/** Number of reads processed across all threads. */
	protected long readsProcessed=0;
	/** Number of bases processed across all threads. */
	protected long basesProcessed=0;
	
	/** Stop after processing this many reads (Long.MAX_VALUE for unlimited). */
	final long maxReads;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Target number of records per chunk before pushing to the queue. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Target number of bytes per chunk before pushing to the queue. */
	public static int TARGET_LIST_BYTES=shared.Shared.bufferSize();
	/** Default number of worker threads when none is specified. */
	public static int DEFAULT_THREADS=3;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status and verbose messages. */
	protected PrintStream outstream=System.err;
	/** Enables verbose logging when true. */
	public static final boolean verbose=false;
	/** Flag set when any processing thread encounters an error. */
	public boolean errorState=false;
	/** Fraction of reads to retain when subsampling is enabled. */
	float samplerate=1f;
	/** Random generator used for subsampling decisions. */
	shared.Random randy=null;
	
}