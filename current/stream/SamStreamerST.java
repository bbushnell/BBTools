package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.LineParser1;
import shared.Shared;
import structures.ListNum;

/**
 * Single-threaded SAM line loader with simple buffering.
 * Simpler alternative to SamLineStreamer for cases where threading overhead isn't worth it.
 * 
 * @author Brian Bushnell, Isla
 * @date November 10, 2025
 */
public class SamStreamerST implements Streamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	SamStreamerST(String fname_, boolean saveHeader_, long maxReads_, boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false), 
			saveHeader_, maxReads_, makeReads_);
	}
	
	/** Constructor. */
	SamStreamerST(FileFormat ffin_, boolean saveHeader_, long maxReads_, boolean makeReads_){
		fname=ffin_.name();
		ffin=ffin_;
		saveHeader=saveHeader_;
		header=(saveHeader ? new ArrayList<byte[]>() : null);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		makeReads=makeReads_;
		outq=new ArrayBlockingQueue<ListNum<SamLine>>(QUEUE_SIZE);
		if(verbose){outstream.println("Made SamLineStreamerST");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		if(verbose){outstream.println("SamLineStreamerST.start() called.");}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Start processing thread
		thread=new ProcessThread();
		thread.start();
		
		if(verbose){outstream.println("SamLineStreamerST started.");}
	}

	@Override
	public synchronized void close(){
		//TODO: Unimplemented
	}
	
	@Override
	public String fname() {return fname;}
	
	@Override
	public boolean hasMore() {return !finished;}
	
	@Override
	public boolean errorState() {return errorState;}
	
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
	
	@Override
	public ListNum<SamLine> nextLines(){
		try{
			ListNum<SamLine> list=outq.take();
			if(verbose && list!=null){
				if(list.last()) {outstream.println("Consumer got terminal list.");}
				else {outstream.println("Consumer got list "+list.id()+" size "+list.size());}
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
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		/** Constructor */
		ProcessThread(){
			setName("SamLineStreamerST-Worker");
		}
		
		/** Called by start() */
		@Override
		public void run(){
			try{
				processFileDirectly();
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				// Send terminal list
				try{
					ListNum<SamLine> terminal=new ListNum<SamLine>(null, -1, false, true);
					outq.put(terminal);
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("ProcessThread terminated.");}
		}
		
		/** Read file directly and process lines in same thread */
		void processFileDirectly() throws InterruptedException{
			if(verbose){outstream.println("Started processFileDirectly.");}
			
			ByteFile.FORCE_MODE_BF2=true;
			ByteFile bf=ByteFile.makeByteFile(ffin);
			
			final LineParser1 lp=new LineParser1('\t');
			long listNumber=0;
			long readID=0;
			int bytes=0;
			
			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<SamLine> ln=new ListNum<SamLine>(new ArrayList<SamLine>(slimit), listNumber++);
			
			for(byte[] line=bf.nextLine(); line!=null && readID<maxReads; line=bf.nextLine()){
				if(line[0]=='@'){
					// Handle header
					if(header!=null) { 
						if(Shared.TRIM_RNAME){line=SamReadInputStream.trimHeaderSQ(line);}
						header.add(line);
					}
				}else{
					// First non-header line: save header
					if(header!=null){
						SamReadInputStream.setSharedHeader(header);
						header=null;
					}
					
					// Apply subsampling if needed
					if(samplerate>=1f || randy==null || randy.nextFloat()<samplerate){
						SamLine sl=new SamLine(lp.set(line));
						ln.add(sl);
						bytes+=(sl.seq==null ? 0 : 2*sl.length());
						
						if(makeReads){
							Read r=sl.toRead(FASTQ.PARSE_CUSTOM);
							sl.obj=r;
							r.samline=sl;
							r.numericID=readID;
							if(!r.validated()){r.validate(true);}
						}
						
						readsProcessedT++;
						basesProcessedT+=(sl.seq==null ? 0 : sl.length());
					}
					readID++;
					
					if(ln.size()>=slimit || bytes>=blimit){
						outq.put(ln);
						ln=new ListNum<SamLine>(new ArrayList<SamLine>(slimit), listNumber++);
						bytes=0;
					}
				}
			}
			
			// Handle leftover header if file had no reads
			if(header!=null){
				SamReadInputStream.setSharedHeader(header);
				header=null;
			}
			
			if(ln.size()>0){
				outq.put(ln);
			}
			
			bf.close();
			if(verbose){outstream.println("Finished processFileDirectly.");}
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
	
	final ArrayBlockingQueue<ListNum<SamLine>> outq;
	private ProcessThread thread;
	private boolean finished=false;
	
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
	private static final int QUEUE_SIZE=8;
	
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