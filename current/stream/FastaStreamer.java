package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Loads FASTA files rapidly with multiple threads.
 * 
 * @author Brian Bushnell
 * @date November 5, 2025
 */
public class FastaStreamer implements Streamer {

	public static void main(String[] args) {
		Timer t=new Timer();
		String fname=args[0];
		if(args.length>1) {DEFAULT_THREADS=Integer.parseInt(args[1]);}

		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		Streamer st=StreamerFactory.makeStreamer(ff, 0, true, -1, true, true);
		st.start();
		long reads=0, bases=0;
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
			for(Read r : ln) {
				reads+=r.pairCount();
				bases+=r.pairLength();
			}
		}
		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public FastaStreamer(String fname_, int threads_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), threads_, pairnum_, maxReads_);
	}

	/** Constructor. */
	public FastaStreamer(FileFormat ffin_, int threads_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		pairnum=pairnum_;
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

	@Override
	public void close(){
		if(bf!=null) {bf.close(); bf=null;}
	}

	@Override
	public String fname() {return fname;}

	@Override
	public boolean hasMore(){
		return oqs.hasMore();
	}

	@Override
	public boolean errorState() {return errorState;}

	@Override
	public boolean paired(){return interleaved;}

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public synchronized long readsProcessed() {return readsProcessed;}

	@Override
	public synchronized long basesProcessed() {return basesProcessed;}

	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	@Override
	public ListNum<Read> nextList(){
		ListNum<Read> list=oqs.getOutput();
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
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

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
			setName("FastaStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}

		/** Called by start() */
		@Override
		public void run(){
			//Process the reads
			synchronized(this) {
				if(tid==0){
					processBytes();
				}else{
					if(interleaved) {
						makeReadsInterleaved();
					}else {
						makeReadsSingle();
					}
				}
			}

			//Indicate successful exit status
			success=true;
			if(verbose){outstream.println("tid "+tid+" terminated.");}
		}

		void processBytes(){
			processBytes0();
			if(verbose){outstream.println("tid "+tid+" done with processBytes0.");}

			// Signal completion via OQS
			oqs.poison();
			if(verbose){outstream.println("tid "+tid+" done poisoning.");}

			//Wait for completion of all threads
			boolean allSuccess=true;
			ThreadWaiter.waitForThreadsToFinish(alpt);
			for(ProcessThread pt : alpt){
				//Wait until this thread has terminated
				if(pt!=this){
					synchronized(pt) {
						synchronized(FastaStreamer.this) {
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

		/** 
		 * Thread 0 reads the actual file and produces lists of byte[] (raw lines).
		 * Each list starts with a '>' line and ends just before the next '>'.
		 * Lists are sent when they reach 200 headers or 200kb, whichever comes first.
		 */
		private void processBytes0(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}

			bf=ByteFile.makeByteFile(ffin);

			long listNumber=0;
			long totalReads=0;

			int headersInList=0;
			int bytesInList=0;

			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<byte[]> ln=new ListNum<byte[]>(new ArrayList<byte[]>(), listNumber++);
			ln.firstRecordNum=totalReads;
			final long limit=maxReads*(interleaved && maxReads<Long.MAX_VALUE/2 ? 2 : 1);
			
			for(byte[] line=bf.nextLine(); line!=null && totalReads<=limit; line=bf.nextLine()){
				if(line.length>0) {
					if(line[0]!='>'){
						ln.add(line);
						bytesInList+=line.length;
					}else {
						//Found a header.
						if((headersInList>=slimit || bytesInList>=blimit) && 
							(!interleaved || ((headersInList&1)==0))){
							oqs.addInput(ln);
							ln=new ListNum<byte[]>(new ArrayList<byte[]>(), listNumber++);
							ln.firstRecordNum=totalReads;
							headersInList=0;
							bytesInList=0;
						}
						if(totalReads<limit) {ln.add(line);}
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
			bf.close();
			if(verbose){outstream.println("tid "+tid+" closed stream.");}
		}

		/** Iterate through the reads */
		void makeReadsSingle(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}

			final ByteBuilder bb=new ByteBuilder(4096);
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
							if(samplerate>=1f || randy.nextFloat()<samplerate){
								Read r=new Read(bb.toBytes(), null, 
									new String(header, 1, header.length-1, StandardCharsets.US_ASCII), readID++, true);
								r.setPairnum(pairnum);
								if(!r.validated()){r.validate(true);}
								reads.add(r);
								readsProcessedT++;
								basesProcessedT+=r.length();
							}
						}
						header=line;
						bb.clear();
					}else{
						bb.append(line);
					}
				}

				// Save final record
				if(header!=null){
					if(samplerate>=1f || randy.nextFloat()<samplerate){
						Read r=new Read(bb.toBytes(), null, 
							new String(header, 1, header.length-1, StandardCharsets.US_ASCII), readID++, true);
						r.setPairnum(pairnum);
						if(!r.validated()){r.validate(true);}
						reads.add(r);
						readsProcessedT++;
						basesProcessedT+=r.length();
					}
				}else {
					throw new RuntimeException("No header for record "+readID+
						" length "+bb.length()+" in "+fname);
				}

				oqs.addOutput(reads);
				list=oqs.getInput();
			}
			if(verbose){outstream.println("tid "+tid+" done making reads.");}
			//Re-inject poison for other workers
			if(list!=null) {oqs.addInput(list);}
		}

		/** Iterate through the reads */
		void makeReadsInterleaved(){
			if(verbose){outstream.println("tid "+tid+" started makeReads.");}

			ListNum<byte[]> list=oqs.getInput();
			final ByteBuilder bb=new ByteBuilder(4096);
			while(list!=null && !list.poison()){
				if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}

				ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>(), list.id());
				long readID=list.firstRecordNum/2;

				// Parse lines into reads using ByteBuilder
				ArrayList<Read> allReads=new ArrayList<Read>();
				byte[] header=null;

				for(byte[] line : list){
					if(line.length>0 && line[0]=='>'){
						// Save previous record if exists
						if(header!=null){
							Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1, StandardCharsets.US_ASCII), 0, true);
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
					Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1, StandardCharsets.US_ASCII), 0, true);
					if(!r.validated()){r.validate(true);}
					allReads.add(r);
					readsProcessedT++;
					basesProcessedT+=r.length();
				}

				// Pair them up
				assert((allReads.size()&1)==0) : "Odd number of reads for interleaved list: "+allReads.size();
				for(int i=0, lim=allReads.size(); i<lim; i+=2){
					if(samplerate>=1f || randy.nextFloat()<samplerate){
						Read r1=allReads.get(i);
						Read r2=allReads.get(i+1);
						r1.setPairnum(0);
						r2.setPairnum(1);
						r1.mate=r2;
						r2.mate=r1;
						reads.add(r1);
						r1.numericID=readID;
						r2.numericID=readID++;
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
	
	public ByteFile bf;//TODO: Should not be a field, just internal.

	final OrderedQueueSystem<ListNum<byte[]>, ListNum<Read>> oqs;

	final int threads;
	final int pairnum;
	final boolean interleaved;

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
	public static int TARGET_LIST_BYTES=262144;
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
	private float samplerate=1f;
	private java.util.Random randy=null;

}