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
 * Streams GFA segment sequences through one parser worker and a bounded output queue.
 * ByteFile may use additional backend threads. Start once and configure sampling and
 * batch targets beforehand; restart and arbitrary concurrent close are unsupported.
 * StreamerFactory currently selects this implementation for GFA regardless of its thread hint.
 *
 * Lines beginning with S consume positional IDs and the input limit before sampling.
 * Other lines are ignored; discarded S lines are not parsed. Retained lines supply
 * tab fields 1 and 2 as an ASCII name and copied bases, with null qualities and
 * shared Read-constructor validation. Extra tags and graph relationships are ignored.
 * This prefix-based sequence reader is not a complete GFA validator.
 *
 * Consuming the terminal copies successful-conversion totals and failure status.
 * Counts can include records in an unpublished batch discarded after a later failure.
 * An interrupted wait also returns null, without that normal terminal update.
 * Callers must inspect errorState rather than treating every null as successful EOF.
 * Cleanup and terminal publication are attempted in finally, but blocking or
 * interruption can prevent delivery; this is not a cancellation protocol.
 *
 * @author Brian Bushnell
 * @contributor Shinobu (correctness, lifecycle and documentation)
 * @date November 21, 2025
 */
public class GfaStreamerST implements Streamer{
	
	/** Drains input, closes in finally and reports delivered reads and bases after checking status.
	 * Explicit SIMD requests honor the normal capability gate. The third argument controls
	 * Read.VALIDATE_VECTOR, which does not disable all parser/backend SIMD operations.
	 * @param args Input path, optionally followed by SIMD and vector-validation booleans
	 * @throws RuntimeException If the reader reports an error after draining and closing
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		String fname=args[0];
		//[stream/GfaStreamerST#004] Explicit true still requires the normal module and CPU capability gate.
		if(args.length>1){Shared.SIMD=Parse.parseBoolean(args[1]) && simd.Vector.simd256;}
		if(args.length>2){Read.VALIDATE_VECTOR=Parse.parseBoolean(args[2]);}
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.GFA, null, true, true);
		GfaStreamerST st=new GfaStreamerST(ff, 0, -1);
		long reads=0, bases=0;
		try{
			st.start();
			for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
				for(Read r : ln){
					reads+=r.pairCount();
					bases+=r.pairLength();
				}
			}
		}finally{st.close();}
		//[stream/GfaStreamerST#002] A consumed terminal does not imply successful input.
		if(st.errorState()){
			throw new RuntimeException("GfaStreamerST encountered an error while reading "+fname);
		}
		t.stop();//[stream/GfaStreamerST#003] Publish elapsed time before calculating throughput.
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resolves an input descriptor without opening the reader's backend.
	 * @param fname_ Input path; descriptor permits subprocess input
	 * @param pairnum_ Side marker, 0 or 1; this reader does not link mates
	 * @param maxReads_ S-prefix line limit before sampling; negative means unlimited
	 */
	public GfaStreamerST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.GFA, null, true, false), pairnum_, maxReads_);
	}
	
	/** Captures the descriptor, side marker and input limit; the worker opens input later.
	 * @param ffin_ Input descriptor
	 * @param pairnum_ Side marker, asserted to be 0 or 1
	 * @param maxReads_ S-prefix line limit before sampling; negative means unlimited
	 */
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
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resets exposed counters and launches a parser worker; call once on a fresh reader.
	 * Repeated calls do not reset the queue, finished flag or error state.
	 */
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
	
	/** Closes the backend and latches returned errors or exceptional completion.
	 * Clears the backend reference only after normal return; repeated successful close
	 * is then a no-op. Does not join the worker, mark finished or publish terminal.
	 * Coordinate with processing; this unsynchronized method is not cancellation.
	 */
	@Override
	public void close(){
		if(bf==null){return;}
		boolean completed=false;
		try{
			//[stream/GfaStreamerST#005] Latch backend errors without overwriting an existing worker error.
			if(bf.close()){errorState=true;}
			completed=true;
			bf=null;
		}finally{
			if(!completed){errorState=true;}
		}
	}
	
	/** Returns the captured input name. */
	@Override
	public String fname(){return fname;}
	
	/** Returns !finished without inspecting the queue or worker; not a nonblocking readiness test. */
	@Override
	public boolean hasMore(){return !finished;}
	
	/** Returns the currently observed error flag, without joining or synchronizing. */
	@Override
	public boolean errorState(){return errorState;}
	
	/** Returns false; records are not linked into pairs by this reader. */
	@Override
	public boolean paired(){return false;}

	/** Returns the configured side marker, independently of paired(). */
	@Override
	public int pairnum(){return pairnum;}
	
	/** Returns successful retained conversions copied when a normal terminal is consumed.
	 * May include records never delivered because a later parse failure discarded their batch.
	 */
	@Override
	public long readsProcessed(){return readsProcessed;}
	
	/** Returns converted bases copied at normal terminal consumption, not necessarily delivered bases. */
	@Override
	public long basesProcessed(){return basesProcessed;}
	
	/** Configures sampling before start; rates at least one bypass RNG and retain every S line.
	 * Smaller rates obtain a fresh generator from Shared.threadLocalRandom(seed).
	 * Discarded S lines still consume IDs and limits but bypass parsing and Read validation.
	 * @param rate Retention probability, normally in [0, 1]
	 * @param seed Seed passed to the generator factory
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}
	
	/** Takes the next data batch, blocking until data, terminal or interruption.
	 * Data batch IDs start at zero; firstRecordNum remains -1. With positive batch targets,
	 * no empty data batches are emitted. Terminal consumption marks finished, copies
	 * conversion totals, folds worker status and reinserts the terminal for later calls.
	 * An interrupted take returns null and sets errorState without those terminal updates.
	 * @return A data batch, or null for terminal/interruption; inspect errorState
	 */
	@Override
	public ListNum<Read> nextList(){
		try{
			ListNum<Read> list=outputQueue.take();
			assert(list!=null) : "Pulled null list.";//Should never happen
			if(verbose){
				if(list==null || list.last()){outstream.println("Consumer got terminal list.");}else{outstream.println("Consumer got list "+list.id());}
			}
			if(list==null || list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				//errorState-fold-clobber [family sweep 2026-06-22, same as stream/FastaStreamerST#004]:
				//Fold, do not overwrite: successful processing may still have a backend close error.
				//Preserve that status for API callers and main's error rejection.
				errorState|=!thread.success;
				outputQueue.add(list);//Re-inject
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}
	
	/** GFA input does not expose SAM-line batches.
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException("GFA does not support SamLine");}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Owns parsing, sampling draws and conversion counters for one start invocation. */
	private class ProcessThread extends Thread{
		
		/** Sets a diagnostic worker name without opening input. */
		ProcessThread(){setName("GfaStreamerST-Worker");}
		
		/** Parses input and attempts cleanup before terminal publication.
		 * Processing Exceptions are logged and flagged; processing Errors escape after finally.
		 * Cleanup exceptions can replace an active failure and are outside the processing catch.
		 * A blocked close or blocked/interrupted terminal put can prevent delivery. A put
		 * interruption is printed without setting an additional error flag.
		 */
		@Override
		public void run(){
			try{
				processSingle();
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				//[stream/GfaStreamerST#006] Attempt cleanup on normal and exceptional processing exits.
				try{close();}finally{
					//Retain the terminal attempt even if backend cleanup throws.
					try{
						ListNum<Read> terminal=new ListNum<Read>(null, -1, ListNum.LAST);
						outputQueue.put(terminal);
					}catch(InterruptedException e){
						e.printStackTrace();
					}
				}
			}
			if(verbose){outstream.println("ProcessThread terminated.");}
		}
		
		/** Opens the backend and parses retained S-prefix lines up to the pre-sampling limit.
		 * Captures batch targets for retained entries/bases; flushes after adding
		 * a whole record and emits a final nonempty batch. Cleanup belongs to run's finally.
		 * @throws InterruptedException If publishing a data batch is interrupted
		 */
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
				if(line.length<1 || line[0]!='S'){continue;}
				
				if(samplerate>=1f || randy.nextFloat()<samplerate){
					Read r=toRead(line, pairnum, readID);
					ln.add(r);
					bytes+=r.length();
				}
				readID++;//[stream/GfaStreamerST#001] count every S-line so maxReads terminates the loop
				//and each read gets a unique numericID; formerly readID stayed 0 -> maxReads ignored + all ids 0.

				if(ln.size()>=slimit || bytes>=blimit){
					outputQueue.put(ln);
					ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);
					bytes=0;
				}
			}
			
			if(ln.size()>0){
				outputQueue.put(ln);
			}
			if(verbose){outstream.println("Finished processSingle.");}
		}
		
		/** Converts tab fields 1/2 and increments totals after successful Read construction.
		 * Uses ASCII names, copied bases and null qualities; Read.VALIDATE_IN_CONSTRUCTOR controls validation.
		 * Does not validate complete GFA syntax, tags or graph relationships.
		 * @param line Retained S-prefix line
		 * @param pairnum Configured side marker
		 * @param id Position among S-prefix input lines, including discarded lines
		 * @return The constructed unlinked record
		 */
		private Read toRead(byte[] line, int pairnum, long id){
			lp.set(line);
			String h=lp.parseString(1);
			byte[] bases=lp.parseByteArray(2);
			Read r=new Read(bases, null, h, id);
			r.setPairnum(pairnum);

			readsProcessedT++;
			basesProcessedT+=r.length();
			return r;
		}

		/** Reused field parser owned by this worker. */
		private final LineParser1 lp=new LineParser1('\t');
		/** Successful retained conversions, including any unpublished records. */
		protected long readsProcessedT=0;
		/** Bases in successful retained conversions, including any unpublished records. */
		protected long basesProcessedT=0;
		/** True after normal processing return; cleanup may still fail afterward. */
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Captured input name. */
	public final String fname;
	
	/** Descriptor used when the worker opens its backend. */
	final FileFormat ffin;
	
	/** Bounded handoff queue for data batches and a reinserted terminal marker. */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;
	
	/** Worker from the latest start call; repeated start is not supported. */
	private ProcessThread thread;
	
	/** Backend opened by the worker and cleared only after normal close return. */
	private ByteFile bf;
	
	/** Side marker placed on each unlinked record. */
	final int pairnum;
	
	/** Converted read count copied from the worker on normal terminal consumption. */
	protected long readsProcessed=0;
	/** Converted base count copied from the worker on normal terminal consumption. */
	protected long basesProcessed=0;
	
	/** S-prefix input-line limit before sampling; negative constructor values become Long.MAX_VALUE. */
	final long maxReads;
	
	/** Set when nextList consumes a terminal; an interrupted take does not set it. */
	private boolean finished=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained-entry batch target; default captured at class initialization, value captured by processing.
	 * Configure a positive value before start; positivity is not checked here.
	 */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Retained-base batch target, captured by processing; a whole record may exceed it.
	 * Configure a positive value before start.
	 */
	public static int TARGET_LIST_BYTES=262144;
	/** Maximum queued batches, including terminal. */
	private static final int QUEUE_SIZE=4;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Destination for optional verbose diagnostics. */
	protected PrintStream outstream=System.err;
	/** Compile-time switch for verbose diagnostics. */
	public static final boolean verbose=false;
	/** Latched processing, close or consumer-interruption status; not a synchronized live view. */
	public boolean errorState=false;
	/** Retention probability for S-prefix lines; configure before start. */
	private float samplerate=1f;
	/** Fresh sampling generator for rates below one; null when sampling is bypassed. */
	private shared.Random randy=null;
	
}
