package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
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
		//[SamStreamerST#001] implemented (was an empty TODO): close+fold the reader so a worker-exception path — which skips the worker's own bf.close() — doesn't leak the handle. On the normal path the worker nulls bf, so this is a clean no-op then. Same field-based close as FastaStreamerST.
		if(bf!=null){errorState|=bf.close(); bf=null;}
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
			assert(list!=null) : "Pulled null list.";//Should never happen
			if(verbose){
				if(list==null || list.last()) {outstream.println("Consumer got terminal list.");}
				else {outstream.println("Consumer got list "+list.id());}
			}
			if(list==null || list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				errorState|=!thread.success;//errorState-fold-clobber [family sweep 2026-06-22, same as stream/FastaStreamerST#004]: |= not =, else a plain '=' OVERWRITES the worker's reader-fold (errorState|=bf.close() on truncated input parsed to EOF→success=true) → truncation silently dropped. Crash-loud-on-truncation restored.
				outq.add(list);//Re-inject
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
				//COMPREHENSION [OQS-death-safe, same pattern as FastaStreamerST]: terminal pushed in the FINALLY → delivered on EVERY worker exit (incl. an Error that catch(Exception) misses) → the consumer's outq.take() in nextLines never strands; errorState|=!success there reads the worker's outcome.
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
			
			ByteFile.FORCE_MODE_BF2=true;//[SamStreamerST#003 QUESTION — FLOAT to Brian] permanent JVM-wide global mutation: forces BF2 mode for EVERY ByteFile created afterward and never resets. Presumably the intended perf default for SAM streaming, but it's a side-effect on unrelated readers. Intended?
			bf=ByteFile.makeByteFile(ffin);//[SamStreamerST#001] use the bf FIELD (not a local) so the implemented close() can reach it on the worker-exception path; matches FastaStreamerST#002
			
			final LineParser1 lp=new LineParser1('\t');
			long listNumber=0;
			long readID=0;
			int bytes=0;
			
			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<SamLine> ln=new ListNum<SamLine>(new ArrayList<SamLine>(slimit), listNumber++);
			
			for(byte[] line=bf.nextLine(); line!=null && readID<maxReads; line=bf.nextLine()){
				if(line.length==0){continue;}//[SamStreamerST#002] skip blank lines: SAM has none, but ByteFile can return a length-0 line (concat / hand-edit artifact) → line[0] below would AIOOBE on the worker thread. The canonical FastaStreamerST guards length too.
				if(line[0]=='@'){
					// Handle header
					if(header!=null) { 
						if(Shared.TRIM_READ_DESCRIPTION){line=SamReadInputStream.trimHeaderSQ(line);}
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
			
			errorState|=bf.close(); bf=null;//Fold the reader's error state, then null the field so the external close() can't double-close (matches FastaStreamerST)
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
	/** Input source — a FIELD so close() can reach it on the worker-exception path (#001) */
	private ByteFile bf;
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
	
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	public static int TARGET_LIST_BYTES=shared.Shared.bufferSize();
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
	private shared.Random randy=null;
	
}