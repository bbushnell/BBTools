package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Reads GFA formatted sequence data and converts it to Read objects.
 * @author Brian Bushnell
 * @date November 21, 2025
 */
public class GfaStreamerST implements Streamer {
	
	public static void main(String[] args) {
		Timer t=new Timer();
		String fname=args[0];
		if(args.length>1) {Shared.SIMD=Parse.parseBoolean(args[1]);}
		if(args.length>2) {Read.VALIDATE_VECTOR=Parse.parseBoolean(args[2]);}
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.GFA, null, true, true);
		GfaStreamerST st=new GfaStreamerST(ff, 0, -1);
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
	public GfaStreamerST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.GFA, null, true, false), pairnum_, maxReads_);
	}
	
	/** Constructor. */
	public GfaStreamerST(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		
		// Simple output queue
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);
		
		if(verbose){outstream.println("Made GfaStreamerST");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		if(verbose){outstream.println("GfaStreamerST.start() called.");}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Start processing thread
		thread=new ProcessThread();
		thread.start();
		
		if(verbose){outstream.println("GfaStreamerST started.");}
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
	public boolean paired(){return false;}

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
			assert(list!=null) : "Pulled null list.";//Should never happen
			if(verbose){
				if(list==null || list.last()) {outstream.println("Consumer got terminal list.");}
				else {outstream.println("Consumer got list "+list.id());}
			}
			if(list==null || list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				errorState=!thread.success;
				outputQueue.add(list);//Re-inject
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
		throw new UnsupportedOperationException("GFA does not support SamLine");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		/** Constructor */
		ProcessThread(){
			setName("GfaStreamerST-Worker");
		}
		
		/** Called by start() */
		@Override
		public void run(){
			try{
				processSingle();
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
				byte[] line=bf.nextLine();
				if(line==null){break;}
				if(line.length<1 || line[0]!='S') {continue;}
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){
					Read r=toRead(line, pairnum, readID);
					ln.add(r);
					bytes+=r.length();
				}
				
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
		
		private Read toRead(byte[] line, int pairnum, long id) {
			lp.set(line);
			String h=lp.parseString(1);
			byte[] bases=lp.parseByteArray(2);
			Read r=new Read(bases, null, h, id);
			r.setPairnum(pairnum);

			readsProcessedT++;
			basesProcessedT+=r.length();
			return r;
		}

		private final LineParser1 lp=new LineParser1('\t');
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

	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
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
	private shared.Random randy=null;
	
}