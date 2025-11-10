package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.LineParser1;
import shared.Shared;
import shared.Tools;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Multithreaded SAM file reader using OrderedQueueSystem.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 4, 2016
 */
public class SamStreamer implements Streamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public SamStreamer(String fname_, int threads_, boolean saveHeader_, boolean ordered_, 
			long maxReads_, boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), threads_, 
			saveHeader_, ordered_, maxReads_, makeReads_);
	}
	
	/** Constructor. */
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

	@Override
	public synchronized void close(){
		//TODO: Unimplemented
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
	public long readsProcessed() {return readsProcessed;}
	
	@Override
	public long basesProcessed() {return basesProcessed;}
	
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
				oqs.setFinished();
			}
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
	
	private class ProcessThread extends Thread {
		
		/** Constructor */
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			setName("SamStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}
		
		/** Called by start() */
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
		
		/** 
		 * Input thread reads lines from file and produces byte[] lists.
		 */
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
		
		/** Worker threads parse byte[] into SamLines */
		void makeReads(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}
			
			final LineParser1 lp=new LineParser1('\t');
			ListNum<byte[]> list=oqs.getInput();
			while(!list.poison()){
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
			oqs.addInput(list);
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** True only if this thread has completed successfully */
		boolean success=false;
		/** Thread ID */
		final int tid;
		
		ArrayList<ProcessThread> alpt;
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
	
	/** Quit after processing this many input reads */
	final long maxReads;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int TARGET_LIST_SIZE=200;
	public static int TARGET_LIST_BYTES=250000;
	public static int DEFAULT_THREADS=3;
	
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
	java.util.Random randy=null;
	
}