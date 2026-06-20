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
 * Multithreaded BAM file reader using OrderedQueueSystem.
 * Input thread reads BAM binary and converts to intermediate format.
 * Worker threads convert to SamLine objects.
 *
 * @author Chloe, Isla
 * @date November 10, 2025
 */
public class BamStreamer implements Streamer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public BamStreamer(String fname_, int threads_, boolean saveHeader_, 
		boolean ordered_, long maxReads_, boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.BAM, null, true, false), threads_, 
			saveHeader_, ordered_, maxReads_, makeReads_);
	}

	/** Constructor. */
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

	@Override
	public synchronized void close(){
		//TODO: Unimplemented
		//TODO: Possible bug [stream/BamStreamer#002] - LOW/MED: close() is a no-op, so abandoning the stream before EOF leaves the non-daemon ProcessThreads (1 input + N workers) blocked forever (workers in oqs.getInput(), input in ThreadWaiter) with no termination signal. Couples with #001. Structural: close() should poison the OQS + interrupt/join the threads. Author-acknowledged TODO.
	}
	
	@Override
	public String fname() {return fname;}
	
	@Override
	public boolean hasMore() {return oqs.hasMore();}
	
	@Override
	public boolean paired(){return false;}

	@Override
	public int pairnum(){return 0;}
	
	@Override
	public long readsProcessed(){return readsProcessed;}
	
	@Override
	public long basesProcessed(){return basesProcessed;}
	
	public long bytesProcessed(){return bytesProcessed;}
	
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	@Override
	public ListNum<Read> nextList(){return nextReads();}
	
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
			//[#001 FIXED] crash LOUD on a producer error instead of silently truncating: if the input thread
			//hit a corrupt/non-BAM error it set errorState + poisoned the OQS (so we got null/last here). A
			//bare `return null` would look like clean EOF -> wrong/partial results. KillSwitch.assertDie exits
			//the process loudly (the BBTools error contract: crash, never silently wrong).
			if(errorState){KillSwitch.kill("Error reading BAM file (corrupt or truncated): "+fname);}
			return null;
		}
		return list;
	}
	
	@Override
	public boolean errorState() {return errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawn process threads */
	void spawnThreads(){
		final int threads=this.threads+1;

		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt));
		}
		if(verbose){outstream.println("Spawned threads.");}

		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private class ProcessThread extends Thread {

		/** Constructor */
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("BamStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}

		/** Called by start() */
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
					reader.readUint8(); //Skip NUL
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

		/** Worker threads convert BAM records to SamLines */
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
			if(list!=null) {oqs.addInput(list);}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** Number of bytes processed by this thread */
		protected long bytesProcessedT=0;
		/** True only if this thread has completed successfully */
		boolean success=false;
		/** Thread ID */
		final int tid;

		ArrayList<ProcessThread> alpt;
		BamToSamConverter converter;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;
	
	/** Primary input file */
	final FileFormat ffin;
	
	final OrderedQueueSystem<ListNum<byte[]>, ListNum<SamLine>> oqs;
	
	final int threads;
	final boolean saveHeader;
	final boolean makeReads;
	
	ArrayList<byte[]> header;
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	private long bytesProcessed=0;
	
	/** Quit after processing this many input reads */
	final long maxReads;
	
	/** Shared BAM to SAM converter (created by input thread) */
	private volatile BamToSamConverter sharedConverter;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	public static int TARGET_LIST_BYTES=shared.Shared.bufferSize();
	public static int DEFAULT_THREADS=6; // BAM benefits from more threads; peaks at 7 + 12 bgzip threads
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	float samplerate=1f;
	shared.Random randy=null;
	
}