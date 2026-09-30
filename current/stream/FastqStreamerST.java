package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Four-line FASTQ reader with one parser worker and a bounded output queue.
 * Construction captures input configuration; start opens input on the worker.
 * The selected ByteFile backend may add its own threads. Configure shared parser
 * settings and sampling before starting, then drain and close from coordinated
 * caller code. Start once; restart and arbitrary concurrent close are not supported.
 *
 * Limits count physical records, or pairs when interleaved, before sampling.
 * Retained reads keep positional numeric IDs; paired mates share an ID. Sampling
 * can produce empty nonterminal batches and bypasses decoder checks for discarded
 * records. Retained records use the nucleotide decoder and shared quality settings;
 * this class is not a complete FASTQ validator.
 *
 * Consuming the terminal batch copies successful-conversion totals and failure
 * status. Those totals can include an unpublished first mate if its second mate
 * fails conversion. An interrupted wait also returns null, but does not perform
 * that terminal update. Callers must inspect errorState rather than treating every
 * null result as successful exhaustion. Current StreamerFactory does not select
 * this implementation; it remains available directly and through its diagnostic main.
 * 
 * @author Brian Bushnell, Isla
 * @contributor Shinobu (correctness, lifecycle and documentation)
 * @date November 10, 2025
 */
public class FastqStreamerST implements Streamer{
	
	/**
	 * Drains the input, closes in finally and reports delivered reads, bases and batches.
	 * The optional SIMD request honors the normal capability gate; the third argument
	 * controls Read.VALIDATE_VECTOR, which does not disable all SIMD operations.
	 * @param args Input path, optionally followed by SIMD and vector-validation booleans
	 * @throws RuntimeException If the reader reports an error after draining and closing
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		String fname=args[0];
		//[stream/FastqStreamerST#004] Explicit true still requires the normal module and CPU capability gate.
		if(args.length>1){Shared.SIMD=Parse.parseBoolean(args[1]) && simd.Vector.simd256;}
		if(args.length>2){Read.VALIDATE_VECTOR=Parse.parseBoolean(args[2]);}
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		FastqStreamerST st=new FastqStreamerST(ff, 0, -1);
		long reads=0, bases=0, lists=0;
		try{
			st.start();
			for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
				lists++;//[stream/FastqStreamerST#003] Count every returned batch once.
				for(Read r : ln){
					reads+=r.pairCount();
					bases+=r.pairLength();
				}
			}
		}finally{st.close();}
		//[stream/FastqStreamerST#002] A consumed terminal does not imply successful input.
		if(st.errorState()){
			throw new RuntimeException("FastqStreamerST encountered an error while reading "+fname);
		}
		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
		System.err.println("lists="+lists+", reads="+reads+", bases="+bases);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Resolves a FASTQ descriptor without reading file contents, then captures configuration.
	 * @param fname_ Input path; descriptor permits subprocess input
	 * @param pairnum_ Side marker, 0 or 1; interleaved input requires 0 by assertion
	 * @param maxReads_ Physical record/pair limit before sampling; negative means unlimited
	 */
	public FastqStreamerST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTQ, null, true, false), pairnum_, maxReads_);
	}
	
	/**
	 * Captures the input descriptor, pairing and limit without opening input or starting work.
	 * @param ffin_ Input descriptor supplying name, pairing and backend options
	 * @param pairnum_ Side marker for unpaired input, 0 or 1; interleaved input requires 0
	 * @param maxReads_ Physical record/pair limit before sampling; negative means unlimited
	 */
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
	
	/** Resets published counters and starts one parser worker.
	 * Call once: queue, completion and error state are not reset for a restart.
	 * Backend opening and batch-target capture occur on the worker.
	 */
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
	
	/** Closes the current backend and latches returned errors or exceptional completion.
	 * Clears the backend reference only after normal return; repeated successful close
	 * is a no-op. Does not join the worker, publish a terminal or mark finished. Coordinate
	 * with processing; this unsynchronized method is not a cancellation protocol.
	 */
	@Override
	public void close(){
		if(bf==null){return;}
		boolean completed=false;
		try{
			//[stream/FastqStreamerST#005] Preserve backend errors, including those reported after truncated input.
			//Only latch true: a successful close must not overwrite an error already reported by the worker.
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
	
	/** Returns !finished; it does not inspect the queue, input or worker state. */
	@Override
	public boolean hasMore(){return !finished;}
	
	/** Returns the currently observed error flag, without a join or synchronization. */
	@Override
	public boolean errorState(){return errorState;}
	
	/** Reports whether adjacent input records are linked as interleaved mates. */
	@Override
	public boolean paired(){return interleaved;}

	/** Returns the captured side marker; interleaved input uses zero here. */
	@Override
	public int pairnum(){return pairnum;}
	
	/** Returns successful retained individual-read conversions copied at terminal consumption.
	 * Until that update this value remains at its start-reset value. It can include
	 * a converted first mate whose failed second mate prevented delivery of the pair.
	 */
	@Override
	public long readsProcessed(){return readsProcessed;}
	
	/** Returns bases in successful conversions, copied with readsProcessed at terminal consumption. */
	@Override
	public long basesProcessed(){return basesProcessed;}
	
	/** Configures positional sampling before start; no range validation is performed here.
	 * The input-record or pair index determines the decision, and discarded positions
	 * still advance the numeric ID. Discarded records bypass decoder validation.
	 * @param rate Requested retained fraction, normally in [0,1]; zero keeps none, one all
	 * @param seed Nonnegative deterministic seed; negative requests a random seed resolved once
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		sampleSeed=Streamer.resolveSampleSeed(seed);
	}
	
	/** Takes the next batch, including possibly empty nonterminal batches.
	 * Data batch IDs start at zero; firstRecordNum remains unspecified (-1).
	 * A terminal batch copies worker counts and failure status, marks finished and is
	 * reinserted for later calls. Interrupted waiting instead sets errorState and returns
	 * null without updating finished or counters; it does not restore interrupt status.
	 * Terminal delivery depends on the worker completing its publication attempt.
	 * @return Next batch, or null after terminal consumption or interrupted waiting
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
				//Family sweep 2026-06-22: |= preserves worker/backend errors even when processing
				//returned normally (for example, a truncated record). Plain assignment erased them.
				//That review described factory dispatch as commented out; current factory does not select this class.
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
	
	/** Rejects SAM-line access because this reader produces FASTQ Read objects.
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTQ does not support SamLine");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Owns parsing and conversion counts; cleanup precedes its terminal-publication attempt. */
	private class ProcessThread extends Thread{
		
		/** Sets a diagnostic worker name without opening input. */
		ProcessThread(){
			setName("FastqStreamerST-Worker");
		}
		
		/** Selects the captured pairing mode and attempts cleanup before terminal publication.
		 * Processing exceptions are logged and flagged; processing Errors escape after finally.
		 * Cleanup exceptions can replace an active failure. A blocked close or interrupted queue insertion can
		 * prevent terminal delivery, so this is not a general cancellation guarantee.
		 */
		@Override
		public void run(){
			try{
				if(interleaved){
					processInterleaved();
				}else{
					processSingle();
				}
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				//[stream/FastqStreamerST#006] Attempt cleanup on normal and exceptional processing exits.
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
		
		/** Parses four physical lines per input record and samples before conversion.
		 * Limits and IDs count completed physical records, including discarded records.
		 * The size threshold counts retained entries; the byte budget counts twice base
		 * lengths before sampling and can flush an empty batch. Incomplete trailing input
		 * flags an error and retains prior completed data. Cleanup belongs to run's finally.
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
				// Read 4 lines per FASTQ record
				byte[] header=bf.nextLine();
				if(header==null){break;}
				byte[] bases=bf.nextLine();
				byte[] plus=bf.nextLine();
				byte[] quals=bf.nextLine();
				
				if(bases==null || quals==null){
					System.err.println("Incomplete record at end of file: "+new String(header));
					errorState=true;
					break;
				}
				//[stream/FastqStreamerST#007] Normalize +description only; preserve invalid prefixes for the decoder assertion.
				if(plus!=null && plus.length>1 && plus[0]=='+'){plus=PLUS;}

				bytes+=2*bases.length;

				//Positional sampling (Streamer.sampleKeep): same subset as the MT streamer for the same seed
				if(samplerate>=1f || Streamer.sampleKeep(readID, sampleSeed, samplerate)){
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
			if(verbose){outstream.println("Finished processSingle.");}
		}
		
		/** Parses eight physical lines per pair and samples each pair as a unit.
		 * Limits and IDs count physical pairs before sampling. Both retained mates are
		 * converted before linking/delivery; failed R2 conversion can leave R1 counted but
		 * unpublished. Size threshold counts pair entries, despite the half-sized initial
		 * capacity; the byte budget counts twice both sequence lengths before sampling.
		 * Incomplete trailing records/pairs flag an error. Final flush sends only nonempty data.
		 * @throws InterruptedException If publishing a data batch is interrupted
		 */
		void processInterleaved() throws InterruptedException{
			if(verbose){outstream.println("Started processInterleaved.");}

			//Assign the FIELD bf (not a local) so close() can reach it - matches processSingle.
			//[stream/FastqStreamerST#001] interleaved formerly shadowed bf with a local, leaving the
			//field null, so close() (the consumer's BF4-hang-prevention on limited reads) was a no-op.
			bf=ByteFile.makeByteFile(ffin);
			
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
				if(bases1==null || quals1==null || header2==null || bases2==null || quals2==null){
					if(bases1==null || quals1==null || bases2==null || quals2==null){
						System.err.println("Incomplete record at end of file: "+new String(header1));
					}
					if(header2==null){
						System.err.println("Incomplete pair at end of file, pairnum "+readID);
					}
					errorState=true;
					break;
				}
				//[stream/FastqStreamerST#007] Preserve invalid separator prefixes on either mate for the decoder assertion.
				if(plus1!=null && plus1.length>1 && plus1[0]=='+'){plus1=PLUS;}
				if(plus2!=null && plus2.length>1 && plus2[0]=='+'){plus2=PLUS;}

				bytes+=2*(bases1.length+bases2.length);

				//Positional sampling by pair index (Streamer.sampleKeep); the pair is a unit
				if(samplerate>=1f || Streamer.sampleKeep(readID, sampleSeed, samplerate)){
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
			if(verbose){outstream.println("Finished processInterleaved.");}
		}
		
		/** Converts and validates a retained nucleotide record, then increments worker totals.
		 * Shared quality detection and normalization apply; encoding conversion still occurs
		 * when CHANGE_QUALITY is false. Counters advance only after successful conversion
		 * and validation, independently of subsequent batch publication or mate conversion.
		 * @param quad Header, bases, separator and quality lines
		 * @param pairnum Side marker to assign to the converted read
		 * @param id Physical input-record or pair position
		 * @return Converted and validated read
		 */
		private Read quadToRead(byte[][] quad, int pairnum, long id){
			Read r=FASTQ.quadToReadVec(quad, id, 0, fname);
			r.setPairnum(pairnum);
			
			if(!r.validated()){r.validate(true);}

			readsProcessedT++;
			basesProcessedT+=r.length();
			return r;
		}

		/** Successful retained individual-read conversions, including any unpublished mate. */
		protected long readsProcessedT=0;
		/** Bases in successful retained conversions, independently of delivery. */
		protected long basesProcessedT=0;
		/** True when processing returned normally; subsequent cleanup may still report an error. */
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Captured input path/name. */
	public final String fname;
	
	/** Captured descriptor used when the worker opens its backend. */
	final FileFormat ffin;
	
	/** Bounded data/terminal queue; empty data batches are distinct from terminal batches. */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;
	
	/** Parser worker created by start; its totals are copied on terminal consumption. */
	private ProcessThread thread;
	
	/** Worker-opened backend, also reached by caller close; clear only after normal closure. */
	private ByteFile bf;
	
	/** Side marker for unpaired input; interleaved construction requires zero. */
	final int pairnum;
	/** Captured descriptor pairing mode. */
	final boolean interleaved;
	
	/** Successful-conversion total published by terminal consumption, not a live worker counter. */
	protected long readsProcessed=0;
	/** Base total published alongside readsProcessed at terminal consumption. */
	protected long basesProcessed=0;
	
	/** Maximum physical input records or pairs before sampling; negative arguments become Long.MAX_VALUE. */
	final long maxReads;
	
	/** Set by terminal consumption, not by close, worker exit or an interrupted take. */
	private boolean finished=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained-entry batch target, initialized from Shared.bufferLen and captured by processing. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Batch budget of twice input sequence lengths, including discarded records; excludes headers. */
	public static int TARGET_LIST_BYTES=262144;
	/** Maximum queued data/terminal batches. */
	private static final int QUEUE_SIZE=4;
	/** Canonical separator replacing a valid plus-prefixed description before decoding. */
	private static final byte[] PLUS=new byte[]{(byte)'+'};
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Destination for optional verbose diagnostics; main/errors also write directly to System.err. */
	protected PrintStream outstream=System.err;
	/** Compile-time verbose diagnostic switch. */
	public static final boolean verbose=false;
	/** Recorded processing, cleanup or waiting failure; not an independent completion barrier. */
	public boolean errorState=false;
	/** Retained fraction configured before start; no range validation in this class. */
	private float samplerate=1f;
	/** Positional sampling seed; setter resolves a negative request once before processing. */
	private long sampleSeed=17;
	
}
