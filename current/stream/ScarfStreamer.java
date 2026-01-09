package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.LineParser1;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import structures.ListNum;

/**
 * Single-threaded Scarf loader.
 * Handles the ancient "Header:Sequence:Qualities" format.
 * Assumes ASCII-64 quality encoding (Phred+64).
 * @author Collei, Brian Bushnell
 * @date November 21, 2025
 */
public class ScarfStreamer implements Streamer {
	
	public static void main(String[] args) {
		Timer t=new Timer();
		String fname=args[0];
		if(args.length>1) {Shared.SIMD=Parse.parseBoolean(args[1]);}
		if(args.length>2) {Read.VALIDATE_VECTOR=Parse.parseBoolean(args[2]);}
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.SCARF, null, true, true);
		ScarfStreamer st=new ScarfStreamer(ff, 0, -1);
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
	public ScarfStreamer(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.SCARF, null, true, false), pairnum_, maxReads_);
	}
	
	/** Constructor. */
	public ScarfStreamer(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		
		// Simple output queue
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);
		
		if(verbose){outstream.println("Made ScarfStreamerST");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		if(verbose){outstream.println("ScarfStreamerST.start() called.");}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Start processing thread
		thread=new ProcessThread();
		thread.start();
		
		if(verbose){outstream.println("ScarfStreamerST started.");}
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
	public boolean paired(){return false;} //Scarf is typically single-ended or in separate files

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
		throw new UnsupportedOperationException("SCARF does not support SamLine");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		/** Constructor */
		ProcessThread(){
			setName("ScarfStreamerST-Worker");
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
				// Read 1 line per Scarf record
				byte[] line=bf.nextLine();
				if(line==null){break;}
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){
					Read r=scarfToRead(line, readID);
					if(r!=null){
						ln.add(r);
						bytes+=r.length(); //Estimate size
					}
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

		private Read scarfToRead(byte[] line, long id) {

			//Scarf format: Header:Sequence:Qualities
			//Parse from right to left to allow colons in header

			int a=-1, b=-1;
			final byte colon=':';
			for(int i=line.length-1; i>=0; i--){
				if(line[i]==colon){
					if(b<0){b=i;}
					else{
						assert(a<0);
						a=i;
						break;
					}
				}
			}

			if(a<0 || b<0){
				//Malformed line; ignore or crash? 
				//FastqStreamer crashes on bad format, so we assert validity.
				if(verbose){System.err.println("Skipping malformed scarf line: "+new String(line));}
				return null; 
			}

			//Copy arrays to create new Read object
			String header=new String(line, 0, a);
			byte[] bases=Arrays.copyOfRange(line, a+1, b);
			byte[] quals=Arrays.copyOfRange(line, b+1, line.length);

			//Convert Qualities from Phred+64 (ASCII) to Phred+0 (Numeric)
			//Internal BBTools standard is numeric 0-41.
			Vector.applyQualOffset(quals, bases, -ASCII_OFFSET);

			Read r=new Read(bases, quals, header, id, true);//True forces validation
			r.setPairnum(pairnum);

			readsProcessedT++;
			basesProcessedT+=r.length();
			return r;
		}
		
		/** Parses colons */
		private final LineParser1 lp=new LineParser1(':');
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
	private static final int ASCII_OFFSET=64; //Scarf is Phred+64
	
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