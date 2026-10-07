package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.Parse;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import structures.ListNum;

/**
 * Loads SCARF records on one parsing worker and publishes newly allocated batches.
 * Each line contains Header:Sequence:Qualities; the two rightmost colons delimit
 * the fields, allowing colons in the header. Quality text is decoded as Phred+64
 * through Vector.applyQualOffset, followed by Read validation.
 * Configure sampling before starting, start once, and consume through nextList()
 * until terminal input. This reader assigns a pair number but does not link mates
 * or detect interleaving. The input-line limit applies before sampling; public
 * counters describe produced reads and are copied when a terminal batch is consumed.
 * The caller owns delivered batches; this class does not recycle them.
 * Concurrent configuration, repeated startup and close during parsing are not
 * coordinated by this implementation. ByteFile may add its own input threads.
 * @author Collei, Brian Bushnell
 * @date November 21, 2025
 */
public class ScarfStreamer implements Streamer{
	
	/** Counts delivered records from a SCARF path and prints elapsed-time statistics.
	 * @param args Input path, optional SIMD flag and optional vector-validation flag */
	public static void main(String[] args){
		Timer t=new Timer();
		String fname=args[0];
		//Explicit requests still require the optional module and a supported vector width.
		if(args.length>1){Shared.SIMD=Parse.parseBoolean(args[1]) && Vector.simd256;}
		if(args.length>2){Read.VALIDATE_VECTOR=Parse.parseBoolean(args[2]);}
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.SCARF, null, true, true);
		ScarfStreamer st=new ScarfStreamer(ff, 0, -1);
		st.start();
		long reads=0, bases=0;
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			for(Read r : ln){
				reads+=r.pairCount();
				bases+=r.pairLength();
			}
		}
		//TODO: Probable bug - the diagnostic does not check st.errorState() after
		//terminal consumption, so a recorded read error can leave a successful exit.
		t.stop();//Tools formats stored elapsed time; an unstopped timer yields infinite rates.
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resolves a SCARF-default input descriptor; opening is deferred to the worker.
	 * @param fname_ Input path, resolved with subprocess support enabled
	 * @param pairnum_ Pair label, zero or one; no mate linking is performed
	 * @param maxReads_ Maximum input lines before sampling; negative means unlimited */
	public ScarfStreamer(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.SCARF, null, true, false), pairnum_, maxReads_);
	}
	
	/** Retains input settings and allocates a four-batch output queue without starting.
	 * @param ffin_ Nonnull input descriptor, interpreted as SCARF by this reader
	 * @param pairnum_ Pair label, zero or one, checked when assertions are enabled
	 * @param maxReads_ Maximum input lines before sampling; negative is unlimited, zero reads none */
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
	
	/** Clears public counters and starts a new parsing worker; call once per reader.
	 * Does not reset the queue, finished flag or prior error state for reuse. */
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
	
	/** Closes and clears the currently assigned input, if present.
	 * Does not join/cancel the worker, publish a terminal or fold this close result.
	 * Use after normal consumption; concurrent worker access is not coordinated. */
	@Override
	public void close(){
		if(bf!=null){bf.close(); bf=null;}
	}
	
	/** Returns the retained input path. */
	@Override
	public String fname(){return fname;}
	
	/** Returns true while nextList has not consumed the terminal marker; does not inspect input.
	 * True before startup and after an interrupted take that has not seen a terminal. */
	@Override
	public boolean hasMore(){
		return !finished;
	}
	
	/** Returns the currently recorded error flag without waiting for completion. */
	@Override
	public boolean errorState(){return errorState;}
	
	/** Returns false: records receive a pair label but no linked mate. */
	@Override
	public boolean paired(){return false;}// Scarf is typically single-ended or in separate files

	/** Returns the configured zero/one pair label. */
	@Override
	public int pairnum(){return pairnum;}
	
	/** Returns produced-read count copied on terminal consumption, otherwise initially zero.
	 * Skipped input lines do not contribute to this count. */
	@Override
	public long readsProcessed(){return readsProcessed;}
	
	/** Returns produced bases copied on terminal consumption, otherwise initially zero. */
	@Override
	public long basesProcessed(){return basesProcessed;}
	
	/** Configures line sampling; call before start, not concurrently with parsing.
	 * Sampling precedes parsing, so skipped lines are not validated.
	 * @param rate Retention probability; values at least one bypass the random draw
	 * @param seed Seed passed to Shared.threadLocalRandom when sampling is enabled */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}
	
	/** Waits for a produced batch or terminal marker after startup.
	 * A terminal copies worker totals, folds worker success, marks finished and is
	 * requeued for later calls. Interrupted waits set errorState and return null
	 * without restoring interrupt status or marking finished. Intended for a single
	 * consuming caller; no worker join or batch recycling occurs here.
	 * @return Owned batch, or null for a terminal marker or interrupted wait */
	@Override
	public ListNum<Read> nextList(){
		try{
			ListNum<Read> list=outputQueue.take();
			assert(list!=null) : "Pulled null list.";//Should never happen
			if(verbose){
				if(list==null || list.last()){outstream.println("Consumer got terminal list.");}
				else{outstream.println("Consumer got list "+list.id());}
			}
			if(list==null || list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				errorState|=!thread.success;//errorState-fold-clobber [family sweep 2026-06-22, same as stream/FastaStreamerST#004]: |= not =, else a plain '=' OVERWRITES the worker's reader-fold (errorState|=bf.close() on truncated input parsed to EOF→success=true) → truncation silently dropped. Crash-loud-on-truncation restored.
				outputQueue.add(list);//Re-inject
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}
	
	/** SCARF has no direct SAM-line delivery path.
	 * @throws UnsupportedOperationException Always */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("SCARF does not support SamLine");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Owns ordinary input parsing and publishes batches through the bounded queue. */
	private class ProcessThread extends Thread{
		
		/** Names the single parsing worker; does not open input. */
		ProcessThread(){
			setName("ScarfStreamerST-Worker");
		}
		
		/** Processes input and attempts to publish LAST in finally.
		 * Caught Exceptions set the shared error flag; normal return sets success even
		 * if input close reported an error. Errors are not caught by this handler. */
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
		
		/** Opens input, samples/parses up to the line quota and publishes batches.
		 * Read IDs count all input lines; output batch IDs count emitted batches.
		 * Size/byte limits are captured on entry; the byte estimate counts bases only.
		 * Normal completion closes input and folds its close status. This method has
		 * no finally-close path; terminal publication belongs to run().
		 * @throws InterruptedException If a blocking batch publication is interrupted */
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
						bytes+=r.length();// Estimate size
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
			errorState|=bf.close();//Fold the reader's error state (truncated/corrupt input) so it isn't silently dropped at the streamer boundary
			if(verbose){outstream.println("Finished processSingle.");}
		}

		/** Copies a line's fields, converts qualities and validates a newly allocated Read.
		 * Parses from the right so header colons are preserved; header decoding uses the
		 * platform charset. The missing-separator path retains historical #001 behavior.
		 * Counts only successfully constructed reads and labels each with pairnum.
		 * @param line Nonnull SCARF line without its line terminator
		 * @param id Zero-based source-line ID, including skipped lines
		 * @return New unlinked read, or null for missing separators with assertions off */
		private Read scarfToRead(byte[] line, long id){
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
				//[stream/ScarfStreamer#001] Malformed scarf line (needs >=2 colons for Header:Seq:Qual).
				//Crash-loud under -ea (KillSwitch.assertDie halts the VM with the offending line) so bad
				//input is never silently dropped; under -da the assert vanishes and we fall through to the
				//best-effort skip below. Resolves the author's inline "ignore or crash?" per crash-loud policy.
				assert(false) : KillSwitch.assertDie("Malformed scarf line (need >=2 colons): "+new String(line));
				if(verbose){System.err.println("Skipping malformed scarf line: "+new String(line));}
				return null;
			}

			//Copy arrays to create new Read object
			String header=new String(line, 0, a);
			byte[] bases=Arrays.copyOfRange(line, a+1, b);
			byte[] quals=Arrays.copyOfRange(line, b+1, line.length);

			//Convert Phred+64 text to numeric qualities using the shared offset policy;
			//subsequent Read validation applies the configured quality/base rules.
			Vector.applyQualOffset(quals, bases, -ASCII_OFFSET);

			Read r=new Read(bases, quals, header, id, true);//True forces validation
			r.setPairnum(pairnum);

			readsProcessedT++;
			basesProcessedT+=r.length();
			return r;
		}
		
		/** Number of reads constructed after sampling, not total source lines. */
		protected long readsProcessedT=0;
		/** Bases in constructed reads after validation. */
		protected long basesProcessedT=0;
		/** True when processSingle returned normally; separate close errors can still exist. */
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	public final String fname;
	
	/** Retained input descriptor; opening occurs on the parsing worker. */
	final FileFormat ffin;
	
	/** Bounded FIFO of produced batches and the requeued LAST marker. */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;
	
	/** Most recently started worker; startup is intended once per instance. */
	private ProcessThread thread;
	
	/** Worker-assigned input; public close also accesses it without synchronization. */
	private ByteFile bf;
	
	/** Zero/one label applied to each produced read; does not establish pairing. */
	final int pairnum;
	
	/** Produced reads copied from the worker when a consumer receives LAST. */
	protected long readsProcessed=0;
	/** Produced bases copied from the worker when a consumer receives LAST. */
	protected long basesProcessed=0;
	
	/** Input-line quota before sampling; negative constructor values become Long.MAX_VALUE. */
	final long maxReads;
	
	/** Set when terminal list is received */
	private boolean finished=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Batch entry threshold captured when parsing begins; configure a positive value. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Batch base-count threshold captured when parsing begins; checked after each line. */
	public static int TARGET_LIST_BYTES=262144;
	/** Maximum queued data/terminal batches. */
	private static final int QUEUE_SIZE=4;
	/** Fixed SCARF text quality offset. */
	private static final int ASCII_OFFSET=64;// Scarf is Phred+64
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** Error accumulator shared by worker and consumer; not a completion barrier. */
	public boolean errorState=false;
	/** Pre-start sampling probability; values at least one bypass sampling. */
	private float samplerate=1f;
	/** Random source retained when sampling is configured below one. */
	private shared.Random randy=null;
	
}
