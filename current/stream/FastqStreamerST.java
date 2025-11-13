package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Single-threaded FASTQ loader with simple buffering.
 * Simpler alternative to FastqStreamer for cases where threading overhead isn't worth it.
 * 
 * @author Brian Bushnell, Isla
 * @date November 10, 2025
 */
public class FastqStreamerST implements Streamer {
	
	public static void main(String[] args) {
		Timer t=new Timer();
		String fname=args[0];
		if(args.length>1) {Shared.SIMD=Parse.parseBoolean(args[1]);}
		if(args.length>2) {Read.VALIDATE_VECTOR=Parse.parseBoolean(args[2]);}
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		FastqStreamerST st=new FastqStreamerST(ff, 0, -1);
		st.start();
		long reads=0, bases=0;
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
			for(Read r : ln) {
				reads+=r.pairCount();
				bases+=r.pairLength();
			}
		}
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public FastqStreamerST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTQ, null, true, false), pairnum_, maxReads_);
	}
	
	/** Constructor. */
	public FastqStreamerST(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		
		// Simple output queue
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);
		
		if(verbose){outstream.println("Made FastqStreamerST");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		if(verbose){outstream.println("FastqStreamerST.start() called.");}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Start processing thread
		thread=new ProcessThread();
		thread.start();
		
		if(verbose){outstream.println("FastqStreamerST started.");}
	}
	
	@Override
	public void close(){
		if(bf!=null) {bf.close(); bf=null;}
	}
	
	@Override
	public String fname() {return fname;}
	
	@Override
	public boolean hasMore(){
		return !finished;
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
		try{
			ListNum<Read> list=outputQueue.take();
			if(verbose){
				if(list==null || list.last()) {outstream.println("Consumer got terminal list.");}
				else {outstream.println("Consumer got list "+list.id());}
			}
			if(list!=null && list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				errorState=!thread.success;
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}
	
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTQ does not support SamLine");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		/** Constructor */
		ProcessThread(){
			setName("FastqStreamerST-Worker");
		}
		
		/** Called by start() */
		@Override
		public void run(){
			try{
				if(interleaved) {
					processInterleaved();
				}else {
					processSingle();
				}
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				// Send terminal list
				try{
					ListNum<Read> terminal=new ListNum<Read>(null, -1, ListNum.LAST);
					outputQueue.put(terminal);
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("ProcessThread terminated.");}
		}
		
		void processSingle() throws InterruptedException{
			if(verbose){outstream.println("Started processSingle.");}

			bf=ByteFile.makeByteFile(ffin);
			
			long listNumber=0;
			long readID=0;
			int bytes=0;
			
			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<Read> ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);
			
			while(readID<maxReads){
				// Read 4 lines per FASTQ record
				byte[] header=bf.nextLine();
				if(header==null){break;}
				byte[] bases=bf.nextLine();
				byte[] plus=bf.nextLine();
				byte[] quals=bf.nextLine();
				
				if(bases==null || plus==null || quals==null){break;}
				
				bytes+=2*bases.length;
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){
					byte[][] quad=new byte[][]{header, bases, plus, quals};
					Read r=quadToRead(quad, pairnum, readID);
					ln.add(r);
				}
				readID++;
				
				if(ln.size()>=slimit || bytes>=blimit){
					outputQueue.put(ln);
					ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);
					bytes=0;
				}
			}
			
			if(ln.size()>0){
				outputQueue.put(ln);
			}
			bf.close();
			if(verbose){outstream.println("Finished processSingle.");}
		}
		
		void processInterleaved() throws InterruptedException{
			if(verbose){outstream.println("Started processInterleaved.");}

			ByteFile bf=ByteFile.makeByteFile(ffin);
			
			long listNumber=0;
			long readID=0;
			int bytes=0;
			
			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<Read> ln=new ListNum<Read>(new ArrayList<Read>(slimit/2), listNumber++);
			
			while(readID<maxReads){
				// Read 8 lines per interleaved pair
				byte[] header1=bf.nextLine();
				if(header1==null){break;}
				byte[] bases1=bf.nextLine();
				byte[] plus1=bf.nextLine();
				byte[] quals1=bf.nextLine();
				
				byte[] header2=bf.nextLine();
				byte[] bases2=bf.nextLine();
				byte[] plus2=bf.nextLine();
				byte[] quals2=bf.nextLine();
				
				if(bases1==null || quals1==null || bases2==null || quals2==null){break;}
				
				bytes+=2*(bases1.length+bases2.length);
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){
					byte[][] quad1=new byte[][]{header1, bases1, plus1, quals1};
					byte[][] quad2=new byte[][]{header2, bases2, plus2, quals2};
					Read r1=quadToRead(quad1, 0, readID);
					Read r2=quadToRead(quad2, 1, readID);
					r1.mate=r2;
					r2.mate=r1;
					ln.add(r1);
				}
				readID++;
				
				if(ln.size()>=slimit || bytes>=blimit){
					outputQueue.put(ln);
					ln=new ListNum<Read>(new ArrayList<Read>(slimit/2), listNumber++);
					bytes=0;
				}
			}
			
			if(ln.size()>0){
				outputQueue.put(ln);
			}
			bf.close();
			if(verbose){outstream.println("Finished processInterleaved.");}
		}
		
		private Read quadToRead(byte[][] quad, int pairnum, long id) {
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
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	public final String fname;
	
	/** Primary input file */
	final FileFormat ffin;
	
	/** Output queue */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;
	
	/** Processing thread */
	private ProcessThread thread;
	
	/** Input source */
	private ByteFile bf;
	
	final int pairnum;
	final boolean interleaved;
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	
	/** Quit after processing this many input reads */
	final long maxReads;
	
	/** Set when terminal list is received */
	private boolean finished=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int TARGET_LIST_SIZE=200;
	public static int TARGET_LIST_BYTES=262144;
	private static final int QUEUE_SIZE=4;
	
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