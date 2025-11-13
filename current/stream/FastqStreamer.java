package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Loads FASTQ files rapidly with multiple threads.
 * 
 * @author Isla
 * @date October 30, 2025
 */
public class FastqStreamer implements Streamer {
	
	public static void main(String[] args) {
		Timer t=new Timer();
		String fname=args[0];
		if(args.length>1) {DEFAULT_THREADS=Integer.parseInt(args[1]);}
		if(args.length>2) {Shared.SIMD=true;}
		if(args.length>3) {Read.VALIDATE_VECTOR=Parse.parseBoolean(args[3]);}
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
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
	public FastqStreamer(String fname_, int threads_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTQ, null, true, false), threads_, pairnum_, maxReads_);
	}
	
	/** Constructor. */
	public FastqStreamer(FileFormat ffin_, int threads_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		
		// Create OQS with prototypes for LAST/POISON generation
		ListNum<byte[][]> inputPrototype=new ListNum<byte[][]>(null, 0, ListNum.PROTO);
		ListNum<Read> outputPrototype=new ListNum<Read>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<ListNum<byte[][]>, ListNum<Read>>(
			threads, true, inputPrototype, outputPrototype);
		
		if(verbose){outstream.println("Made FastqStreamer-"+threads);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		if(verbose){outstream.println("FastqStreamer.start() called.");}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads();
		
		if(verbose){outstream.println("FastqStreamer started.");}
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
	public long readsProcessed() {return readsProcessed;}
	
	@Override
	public long basesProcessed() {return basesProcessed;}
	
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
		throw new UnsupportedOperationException("FASTQ does not support SamLine");
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
			setName("FastqStreamer-"+(tid==0 ? "Input" : "Worker-"+tid));
			alpt=(tid==0 ? alpt_ : null);
		}
		
		/** Called by start() */
		@Override
		public void run(){
			//Process the reads
			if(tid==0){
				processBytes();
			}else{
				if(interleaved) {
					makeReadsInterleaved();
				}else {
					makeReadsSingle();
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
		 * Thread 0 reads the actual file and produces lists of byte[][] (4 lines per read).
		 */
		private void processBytes0(){
			if(verbose){outstream.println("tid "+tid+" started processBytes.");}
			
			bf=ByteFile.makeByteFile(ffin);
			
			long listNumber=0;
			long reads=0;
			int bytes=0;
			
			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<byte[][]> ln=new ListNum<byte[][]>(new ArrayList<byte[][]>(slimit), listNumber++);
			ln.firstRecordNum=reads;
			
			while(reads<maxReads){
				// Read 4 lines per FASTQ record
				byte[] header=bf.nextLine();
				if(header==null){break;}
				byte[] bases=bf.nextLine();
				byte[] plus=bf.nextLine();
				byte[] quals=bf.nextLine();
				bytes+=2*bases.length;//Ignore header, usually short
				
				if(bases==null || plus==null || quals==null){
					// Incomplete record at end of file
					break;
				}
				
				reads++;
				byte[][] record=new byte[][]{header, bases, plus, quals};
				ln.add(record);
				
				if(ln.size()>=slimit || bytes>=blimit){
					oqs.addInput(ln);
					ln=new ListNum<byte[][]>(new ArrayList<byte[][]>(slimit), listNumber++);
					ln.firstRecordNum=reads;
					bytes=0;
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
			
			ListNum<byte[][]> list=oqs.getInput();
			while(list!=null && !list.poison()){
				if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}
				
				ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>(list.size()), list.id());
				long readID=list.firstRecordNum;
				
				if(samplerate>=1f){
					for(byte[][] quad : list){
						Read r=quadToRead(quad, pairnum, readID++);
						reads.add(r);
					}
				}else{
					for(byte[][] quad : list){
						if(randy.nextFloat()<samplerate){
							Read r=quadToRead(quad, pairnum, readID++);
							reads.add(r);
						}
					}
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
			
			ListNum<byte[][]> list=oqs.getInput();
			while(list!=null && !list.poison()){
				if(verbose){outstream.println("tid "+tid+" grabbed blist "+list.id());}
				
				ListNum<Read> reads=new ListNum<Read>(new ArrayList<Read>((list.size()+1)/2), list.id());
				long readID=list.firstRecordNum/2;
				ArrayList<byte[][]> quads=list.list;
				assert((quads.size()&1)==0) : "Odd number of quads for interleaved list: "+quads.size();

				if(samplerate>=1f){
					for(int i=0, lim=quads.size(); i<lim; i+=2){
						byte[][] quad1=quads.get(i);
						byte[][] quad2=quads.get(i+1);
						Read r1=quadToRead(quad1, 0, readID);
						Read r2=quadToRead(quad2, 1, readID++);
						r1.mate=r2;
						r2.mate=r1;
						reads.add(r1);
					}
				}else{
					for(int i=0, lim=quads.size(); i<lim; i+=2){
						if(randy.nextFloat()<samplerate){
							byte[][] quad1=quads.get(i);
							byte[][] quad2=quads.get(i+1);
							Read r1=quadToRead(quad1, 0, readID);
							Read r2=quadToRead(quad2, 1, readID++);
							r1.mate=r2;
							r2.mate=r1;
							reads.add(r1);
						}
					}
				}
				
				oqs.addOutput(reads);
				list=oqs.getInput();
			}
			if(verbose){outstream.println("tid "+tid+" done making reads.");}
			//Re-inject poison for other workers
			if(list!=null) {oqs.addInput(list);}
		}
		
		private Read quadToRead(byte[][] quad, int pairnum, long id) {
//			Read r=FASTQ.quadToRead_slow(quad, false, null, readID, 0);
			
			Read r=FASTQ.quadToReadVec(quad, id, 0, fname);
			r.setPairnum(pairnum);
			
			if(!r.validated()){r.validate(true);}

			readsProcessedT++;
			basesProcessedT+=r.length();
			return r;
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
	
	/** Input source */
	private ByteFile bf;
	
	final OrderedQueueSystem<ListNum<byte[][]>, ListNum<Read>> oqs;
	
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
	public static int DEFAULT_THREADS=2;
	
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