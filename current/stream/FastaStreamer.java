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
			//[stream/FastaStreamer#001 FIXED 2026-06-21]: try/finally GUARANTEES oqs.poison() even when processBytes0() throws (corrupt/
			//truncated FASTA from bf.nextLine, OOM, etc.) so the workers in getInput() + the consumer in getOutput() wake instead of hanging;
			//the catch records errorState so nextList crashes LOUD via KillSwitch. Input-thread death leaves NO ordered gap (chunks delivered
			//in order, LAST reachable) so plain poison() suffices here (workers need setFinished(true), see makeReads*). Same OQS-thread-death
			//fix as the greenlit SamStreamer/FastqStreamer/BamWriter#001.
			try{
				processBytes0();
			}catch(Throwable t){
				errorState=true;
				outstream.println("FastaStreamer: error reading "+fname+": "+t);
			}finally{
				oqs.poison();// Signal completion via OQS -- ALWAYS, so workers (getInput) + consumer (getOutput) wake
			}
			if(verbose){outstream.println("tid "+tid+" done with processBytes0 + poisoning.");}

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

			final int slimit=Shared.bufferLen();
			final int blimit=Shared.bufferData();
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
								if(samplerate>=1f || randy.nextFloat()<samplerate){
									Read r=new Read(bb.toBytes(), null,
										new String(header, 1, header.length-1, StandardCharsets.US_ASCII), readID++, flag, true);
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
								new String(header, 1, header.length-1, StandardCharsets.US_ASCII), readID++, flag, true);
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
			}catch(Throwable t){
				errorState=true;
				oqs.setFinished(true);//force-poison outq -> release the consumer past THIS worker's undelivered-job gap
				throw new RuntimeException("FastaStreamer worker "+tid+" failed: "+fname, t);
			}
		}

		/** Iterate through the reads */
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
					ArrayList<Read> allReads=new ArrayList<Read>();
					byte[] header=null;

					for(byte[] line : list){
						if(line.length>0 && line[0]=='>'){
							// Save previous record if exists
							if(header!=null){
								Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1,
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
						Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1,
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
					for(int i=0; i<lim; i+=2){
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
			}catch(Throwable t){
				errorState=true;
				oqs.setFinished(true);//force-poison outq -> release the consumer past THIS worker's undelivered-job gap
				throw new RuntimeException("FastaStreamer worker "+tid+" failed: "+fname, t);
			}
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
	public int flag;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

//	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
//	public static int TARGET_LIST_BYTES=262144;
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
	private shared.Random randy=null;

}