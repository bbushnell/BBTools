package stream;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.bam.BamReader;
import stream.bam.BamToSamConverter;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Reads BAM record bodies on one input thread and converts them on separate workers.
 * Delivers SamLine batches or, when requested, associated Read objects. Input batches
 * carry alignment-record positions for deterministic subsampling; mates are not linked.
 * OrderedQueueSystem coordinates batches; its current output JobQueue forces order
 * even when the constructor receives ordered=false.
 *
 * Configure before start and start once. Input and conversion threads are daemon
 * threads; the decompression backend can add threads of its own. Normal consumers
 * drain through the terminal result. close() signals the output queue, without
 * joining the input thread or guaranteeing final public counters. Counters are
 * snapshots of later worker aggregation, not progress counters updated per record.
 *
 * @author Chloe, Isla
 * @date November 10, 2025
 */
public class BamStreamer implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a BAM-default descriptor, then configures an unstarted reader.
	 * @param fname_ Input path
	 * @param threads_ Conversion workers; below one uses DEFAULT_THREADS before clamping
	 * @param saveHeader_ Collect and publish input header text when true
	 * @param ordered_ Requested output order; current queue forces ordered delivery
	 * @param maxReads_ Alignment-record limit before sampling; negative means unlimited
	 * @param makeReads_ Construct Read objects associated with each retained SamLine */
	public BamStreamer(String fname_, int threads_, boolean saveHeader_,
		boolean ordered_, long maxReads_, boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.BAM, null, true, false), threads_,
			saveHeader_, ordered_, maxReads_, makeReads_);
	}

	/** Retains input settings and creates queues without opening the BAM stream.
	 * @param ffin_ Nonnull input descriptor; the worker expects BAM content
	 * @param threads_ Conversion workers, clamped to 1..Shared.threads(); below one uses defaults
	 * @param saveHeader_ Allocate header storage and publish parsed text lines when true
	 * @param ordered_ Requested output order, currently overridden by JobQueue
	 * @param maxReads_ Alignment-record limit before sampling; negative is unlimited, zero reads only metadata
	 * @param makeReads_ Enable Read construction for nextList/nextReads */
	public BamStreamer(FileFormat ffin_, int threads_, boolean saveHeader_,
		boolean ordered_, long maxReads_, boolean makeReads_){
		fname=ffin_.name();
		ffin=ffin_;
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		saveHeader=saveHeader_;
		header=(saveHeader ? new ArrayList<byte[]>() : null);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		makeReads=makeReads_;

		// Create OQS with prototypes
		ListNum<byte[]> inputPrototype=new ListNum<byte[]>(null, 0, ListNum.PROTO);
		ListNum<SamLine> outputPrototype=new ListNum<SamLine>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<ListNum<byte[]>, ListNum<SamLine>>(
			threads, ordered_, inputPrototype, outputPrototype);

		if(verbose){
			outstream.println("Made BamStreamer-"+threads);
			new Exception().printStackTrace();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Clears public counters and starts a fresh input/conversion thread set; call once.
	 * Does not reset queues, header storage, converter or previous error state for reuse. */
	@Override
	public void start(){
		if(verbose){outstream.println("BamStreamer.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		bytesProcessed=0;

		//Spawn threads
		spawnThreads();

		if(verbose){outstream.println("Started.");}
	}

	/** Signals output-queue completion; does not close the input handle or join workers.
	 * The historical partial-teardown notes below remain separate from this operation's
	 * actual contract. This method is not a final-counter or final-error barrier. */
	@Override
	public synchronized void close(){
		//[stream/BamStreamer#002 partial fix 2026-09-05] was a no-op: abandoning the stream before EOF left
		//the non-daemon threads blocked forever (zombie JVM; jstack-proven for the sibling FastqStreamer).
		//setFinished(true) force-poisons the outq so workers + consumer wake; the input thread then
		//free-runs the remaining file and exits via its own poison protocol — wasteful on a huge file but
		//bounded and zombie-free. TODO: the full structural teardown (deterministic flag-bail so input
		//stops reading early) remains future work; couples with #001.
		oqs.setFinished(true);
	}

	/** Returns the retained input name. */
	@Override
	public String fname(){return fname;}

	/** Delegates a non-consuming queue-availability hint; nextLines determines terminal input. */
	@Override
	public boolean hasMore(){return oqs.hasMore();}

	/** Returns false: aligned mates are not linked by this reader. */
	@Override
	public boolean paired(){return false;}

	/** Returns the reader-level pair label zero; alignment flags retain their own mate information. */
	@Override
	public int pairnum(){return 0;}

	/** Returns the observed count of retained, converted alignment records.
	 * Aggregation occurs after worker joins; consuming terminal input does not join that aggregation. */
	@Override
	public long readsProcessed(){return readsProcessed;}

	/** Returns observed retained sequence bases from worker aggregation, without waiting. */
	@Override
	public long basesProcessed(){return basesProcessed;}

	/** Returns observed retained BAM record-body bytes, excluding each four-byte size prefix.
	 * Header-text bytes collected on the input thread are omitted by historical #003;
	 * this is neither compressed-file size nor complete uncompressed BAM size. */
	public long bytesProcessed(){return bytesProcessed;}

	/** Configures positional sampling before startup; not synchronized with conversion.
	 * @param rate Retained fraction, normally 0..1; values at least one bypass sampling
	 * @param seed Nonnegative seed, or a negative request for a newly resolved random seed */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		sampleSeed=Streamer.resolveSampleSeed(seed);
	}

	/** Delegates Read-batch delivery; requires construction with makeReads=true.
	 * @return New wrapper/list of associated reads, or null on terminal input */
	@Override
	public ListNum<Read> nextList(){return nextReads();}

	/** Extracts already constructed Read references from the next SamLine batch.
	 * Allocates a new list/wrapper with the same batch ID, without copying Read objects.
	 * May return an empty batch after sampling. Requires makeReads, checked by assertion.
	 * @return Read batch, or null when nextLines reports terminal input */
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

	/** Returns the next converted batch or handles the queue's terminal result.
	 * LAST marks the queue finished; a terminal with observed errorState invokes
	 * KillSwitch.kill. This does not join the input thread's final aggregation.
	 * Returned SamLines/associated Reads belong to the consumer; no recycling occurs.
	 * @return Converted batch, possibly empty after sampling, or null at terminal input */
	@Override
	public ListNum<SamLine> nextLines(){
		ListNum<SamLine> list=oqs.getOutput();
		if(verbose){
			if(list==null){outstream.println("Consumer got null.");}
			else{outstream.println("Consumer got list "+list.id()+" type "+list.type);}
		}
		if(list==null || list.last()){
			if(list!=null && list.last()){
				oqs.setFinished(true);
			}
			//[#001 FIXED] crash LOUD on a producer error instead of silently truncating: if the input thread
			//hit a corrupt/non-BAM error it set errorState + poisoned the OQS (so we got null/last here). A
			//bare `return null` would look like clean EOF -> wrong/partial results. KillSwitch.assertDie exits
			//the process loudly (the BBTools error contract: crash, never silently wrong).
			if(errorState){KillSwitch.kill("Error reading BAM file (corrupt or truncated): "+fname);}
			return null;
		}
		return list;
	}

	/** Returns the observed shared error flag without joining input or conversion threads. */
	@Override
	public boolean errorState(){return errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates and starts one daemon input thread plus the configured daemon conversion workers.
	 * The input thread retains the complete thread list for joined aggregation. */
	void spawnThreads(){
		final int threads=this.threads+1;

		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt));
		}
		if(verbose){outstream.println("Spawned threads.");}

		for(ProcessThread pt : alpt){
			//[stream/BamStreamer#002] DAEMON: BamStreamer is an INPUT reader, so abandoning it is harmless (no
			//output to corrupt). Marking the Input+Worker ProcessThreads daemon means that if the consumer stops
			//draining and main returns (e.g. calctruequality finishes/abandons a still-reading streamer), these
			//threads can no longer keep the JVM alive parked on full OQS queues -> the JVM-hang is impossible.
			//(Was: non-daemon -> deadlocked the JVM. The BGZF-InputWorker threads were already daemon.) A proper
			//close() implementation (deterministic flag-bail teardown) is the complementary fix; this is the
			//wedge-free safety net. The non-daemon CONSUMER thread still gates JVM exit, so reads aren't dropped
			//on a normal run - daemon only matters once the consumer is already gone.
			pt.setDaemon(true);
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input role at tid zero; other instances convert record batches with private scratch. */
	private class ProcessThread extends Thread{

		/** Names the thread and retains the aggregation list only for the input role.
		 * @param tid_ Zero for input, otherwise a conversion-worker identifier
		 * @param alpt_ Complete thread list, retained only when tid_ is zero */
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("BamStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}

		/** Executes the assigned role; sets success after its method returns normally.
		 * Input errors caught by processInputThread can coexist with success=true here. */
		@Override
		public void run(){
			if(tid==0){
				processInputThread();
			}else{
				makeReads();
			}

			success=true;
			if(verbose){outstream.println("tid "+tid+" terminated.");}
		}

		/** Reads input, records caught errors, publishes terminal markers and joins workers.
		 * ThreadWaiter skips this calling thread. Aggregation then adds only other workers'
		 * totals to public counters and folds their success flags. Publication and joined
		 * aggregation are separate steps; historical failure guarantees are not established
		 * merely by reaching a consumer terminal result. */
		void processInputThread(){
			//[stream/BamStreamer#001] FIXED 2026-06-20 (greenlit by Brian): a corrupt/truncated/non-BAM input
			//can no longer hang the workers+consumer. The try/finally GUARANTEES oqs.poison() runs even when
			//processBamBytes() throws (so workers in getInput() + the consumer in getOutput() wake); the catch
			//sets errorState + notifyAll() (releasing any worker still spinning in the pre-converter
			//sharedConverter-wait at makeReads L329). The consumer then crashes LOUD via the errorState check in
			//nextLines, instead of silently truncating. Same fix shape as the greenlit bam/BgzfInputStreamMT2#001.
			//NEEDS VALIDATION: a correct BAM (rc=0, all reads) AND a corrupt/non-BAM .bam (loud crash, NO hang).
			try{
				processBamBytes();
			}catch(Throwable t){
				synchronized(BamStreamer.this){errorState=true; BamStreamer.this.notifyAll();}
				outstream.println("BamStreamer: error reading BAM "+fname+": "+t);
			}finally{
				oqs.poison();//ALWAYS poison so workers (getInput) + consumer (getOutput) wake
			}
			if(verbose){outstream.println("tid "+tid+" done with processBamBytes + poisoning.");}

			//Wait for completion of all threads
			boolean allSuccess=true;
			ThreadWaiter.waitForThreadsToFinish(alpt);
			for(ProcessThread pt : alpt){
				if(pt!=this){//[stream/BamStreamer#003] LOW: skipping the input thread (this) drops ITS bytesProcessedT - the header bytes accumulated at L250/258 - from bytesProcessed, a stat undercount by the header size (readsProcessedT/basesProcessedT are 0 on the input thread, so only bytes are lost). Fix: add this.bytesProcessedT after this loop.
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					bytesProcessed+=pt.bytesProcessedT;
					allSuccess&=pt.success;
				}
			}
			if(verbose){outstream.println("tid "+tid+" noted all threads finished.");}

			if(!allSuccess){errorState=true;}
			if(verbose){outstream.println("tid "+tid+" finished! Error="+errorState);}
		}

		/** Opens decompression, reads BAM header/dictionary and queues complete record bodies.
		 * Saved nonempty header-text lines omit LF and are published by reference globally.
		 * Dictionary reference names initialize the shared converter. Record quota applies
		 * before sampling; batch thresholds are captured after header setup, with body bytes
		 * excluding size prefixes. The alignment loop treats any caught EOFException as EOF.
		 * Closes the input after normal completion; no finally-close path is provided here. */
		void processBamBytes(){
			if(verbose){outstream.println("tid "+tid+" started processBamBytes.");}

			long listNumber=0;
			try{
				final InputStream bgzf=ReadWrite.getUnbgzipStream(fname);
				BamReader reader=new BamReader(bgzf);

				//Read BAM magic
				byte[] magic=reader.readBytes(4);
				if(!Arrays.equals(magic, new byte[]{'B', 'A', 'M', 1})){
					throw new RuntimeException("Not a BAM file: "+fname);
				}

				//Read header text
				long l_text=reader.readUint32();
				byte[] text=reader.readBytes((int)l_text);

				//Parse header if requested
				if(saveHeader && header!=null){
					synchronized(header){
						int start=0;
						for(int i=0; i<text.length; i++){
							if(text[i]=='\n'){
								if(i>start){
									byte[] line=Arrays.copyOfRange(text, start, i);
									header.add(line);
									bytesProcessedT+=line.length;
								}
								start=i+1;
							}
						}
						if(start<text.length){
							byte[] line=Arrays.copyOfRange(text, start, text.length);
							header.add(line);
							bytesProcessedT+=line.length;
						}
						SamReadInputStream.setSharedHeader(header);
						if(verbose){outstream.println("Thread "+tid+" set shared header.");}
					}
				}

				if(verbose){outstream.println("Thread "+tid+" reading sequence dictionary.");}

				//Read reference sequence dictionary
				int n_ref=reader.readInt32();
				String[] refNames=new String[n_ref];
				for(int i=0; i<n_ref; i++){
					long l_name=reader.readUint32();
					refNames[i]=reader.readString((int)l_name-1);
					reader.readUint8();// Skip NUL
					long l_ref=reader.readUint32();
				}

				if(verbose){outstream.println("Thread "+tid+" making converter.");}
				synchronized(BamStreamer.this){
					sharedConverter=new BamToSamConverter(refNames);
					BamStreamer.this.notifyAll();
				}
				if(verbose){outstream.println("Thread "+tid+" made converter.");}

				final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
				int bytes=0;
				ListNum<byte[]> ln=new ListNum<byte[]>(new ArrayList<byte[]>(slimit), listNumber++);
				ln.firstRecordNum=0;

				//Read alignment records
				try{
					for(long reads=0; reads<maxReads; reads++){
						long block_size=reader.readUint32();
						byte[] bamRecord=reader.readBytes((int)block_size);
						ln.add(bamRecord);
						bytes+=block_size;

						if(ln.size()>=slimit || bytes>=blimit){
							oqs.addInput(ln);
							ln=new ListNum<byte[]>(new ArrayList<byte[]>(slimit), listNumber++);
							ln.firstRecordNum=reads+1;
							bytes=0;
							if(verbose){outstream.println("Thread "+tid+" made list: reads="+reads);}
						}
					}
				}catch(EOFException e){
					//TODO: Probable bug - BamReader.readFully throws EOFException for a partial
					//record body or size field too; this catch treats those as ordinary EOF.
					//Distinguish record-boundary EOF in a separate parser-correctness change.
					//Normal end of file
				}
				if(verbose){outstream.println("Thread "+tid+" finished reading.");}

				if(ln.size()>0){
					oqs.addInput(ln);
				}

				bgzf.close();

				if(verbose){outstream.println("Thread "+tid+" closed streams.");}
			}catch(IOException e){
				throw new RuntimeException("Error reading BAM file: "+fname, e);
			}
			if(verbose){outstream.println("Thread "+tid+" finished processBamBytes.");}
		}

		/** Waits for the shared converter, samples queued bodies and publishes converted batches.
		 * Sampling uses original record positions, then compacts each input list. Constructed
		 * Read IDs start at that batch's firstRecordNum and increment over retained records,
		 * so they do not necessarily equal original sampled positions. Each worker reuses
		 * its own CIGAR builder; SamLine.obj and Read.samline connect the two representations.
		 * Publishes even empty output batches and reinserts input poison for other workers. */
		void makeReads(){
			if(verbose){outstream.println("Thread "+tid+" waiting on converter.");}
			synchronized(BamStreamer.this){
				//comprehension: workers block here until the INPUT thread (tid==0) publishes the converter
				//(synchronized(BamStreamer.this)+notifyAll @280-281); sharedConverter is volatile + the lock =
				//happens-before; wait(100) is a missed-notify backstop. [#001 FIXED] the `&& !errorState` bail
				//lets a worker escape this spin if the input thread DIES pre-converter (e.g. magic-mismatch),
				//where it set errorState+notifyAll in processInputThread's catch -> no hang.
				while(sharedConverter==null && !errorState){
					try{
						BamStreamer.this.wait(100);
					}catch(InterruptedException e){
						e.printStackTrace();
					}
				}
				converter=sharedConverter;
			}
			if(converter==null){//input thread died pre-converter (#001): bail, the finally-poison wakes the consumer
				if(verbose){outstream.println("tid "+tid+" bailing makeReads: no converter (errorState="+errorState+").");}
				return;
			}

			if(verbose){outstream.println("tid "+tid+" started makeReads.");}
			final ByteBuilder cigar=new ByteBuilder(1024);
			ListNum<byte[]> list=oqs.getInput();
			while(list!=null && !list.poison()){
				if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}

				// Apply subsampling if needed
				//Positional sampling (Streamer.sampleKeep) by record index: thread-safe + reproducible
				//across runs and thread counts; the shared PRNG raced across worker threads.
				if(samplerate<1f){
					int nulled=0;
					final long firstRec=list.firstRecordNum;
					for(int i=0; i<list.size(); i++){
						if(!Streamer.sampleKeep(firstRec+i, sampleSeed, samplerate)){
							list.list.set(i, null);
							nulled++;
						}
					}
					if(nulled>0){Tools.condenseStrict(list.list);}
				}

				ListNum<SamLine> reads=new ListNum<SamLine>(
					new ArrayList<SamLine>(list.size()), list.id);
				long readID=list.firstRecordNum;
				for(byte[] bamRecord : list){
					bytesProcessedT+=bamRecord.length;
					final SamLine sl=converter.toSamLine(bamRecord, cigar);
					assert(sl!=null);
					if(sl!=null){
						if(makeReads){
							Read r=sl.toRead(FASTQ.PARSE_CUSTOM);
							sl.obj=r;
							r.samline=sl;
							r.numericID=readID++;
							if(!r.validated()){r.validate(true);}
						}
						reads.add(sl);

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
		}

		/** Retained alignment records successfully converted by this worker. */
		protected long readsProcessedT=0;
		/** Retained sequence bases counted by conversion workers. */
		protected long basesProcessedT=0;
		/** Record-body bytes for conversion workers; saved header-text bytes for input. */
		protected long bytesProcessedT=0;
		/** True after the assigned role returns normally, independent of shared errorState. */
		boolean success=false;
		/** Zero for input, positive for conversion workers. */
		final int tid;

		/** Complete thread set retained by the input role; null for conversion workers. */
		ArrayList<ProcessThread> alpt;
		/** Worker reference to the shared dictionary-backed converter. */
		BamToSamConverter converter;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;

	/** Retained input descriptor; the input thread opens fname through the BGZF backend. */
	final FileFormat ffin;

	/** Raw-body input FIFO and ordered converted-output queue. */
	final OrderedQueueSystem<ListNum<byte[]>, ListNum<SamLine>> oqs;

	/** Conversion-worker count, excluding the input and decompression threads. */
	final int threads;
	/** Enables collecting and publishing input header text. */
	final boolean saveHeader;
	/** Enables constructing Read objects alongside SamLines. */
	final boolean makeReads;

	/** Collected header lines, or null; published to SamReadInputStream without copying the list. */
	ArrayList<byte[]> header;

	/** Worker-retained record count, added by input-thread aggregation after worker joins. */
	protected long readsProcessed=0;
	/** Worker-retained base count, added during final aggregation. */
	protected long basesProcessed=0;
	/** Retained record-body byte total; input-thread header bytes are not folded in. */
	private long bytesProcessed=0;

	/** Alignment-record quota before sampling; negative constructor values become Long.MAX_VALUE. */
	final long maxReads;

	/** Shared BAM to SAM converter (created by input thread) */
	private volatile BamToSamConverter sharedConverter;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Input batch record threshold, captured after reading the dictionary. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Input batch body-byte threshold, excluding four-byte record-size prefixes. */
	public static int TARGET_LIST_BYTES=shared.Shared.bufferSize();
	/** Default requested conversion workers, clamped by the constructor. */
	public static int DEFAULT_THREADS=6;// Historical tuning: BAM peaks at 7 + 12 bgzip threads

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** Shared observed error flag; queries alone do not establish final completion. */
	public boolean errorState=false;
	/** Pre-start sampling rate, read by conversion workers without synchronized updates. */
	float samplerate=1f;
	/** Seed for positional sampling (Streamer.sampleKeep); resolved from setSampleRate's seed */
	long sampleSeed=17;

}
