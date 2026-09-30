package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.ByteFile1Fc;
import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Shared;
import structures.IntList;
import structures.ListNum;

/**
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * FASTA reader with one parser worker and a bounded output queue.
 * ByteFile1Fc joins sequence lines into record blocks; this reader copies bases
 * and decodes US-ASCII names into unpaired Reads with null qualities.
 * Interleaved input is unsupported. Limits, numeric IDs and counts use retained
 * records after sampling; discarded blocks can leave gaps in batch IDs.
 * Configure before start, then consume through the terminal result and close.
 * Start once; close does not cancel or join the worker. Status getters are plain
 * observations, not cross-thread completion barriers. The normal factory path
 * selects this backend only with its SIMD/version-two conditions satisfied.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 12, 2025
 */
public class FastaStreamer2ST implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Configures an input path and queue without opening the backend.
	 * @param fname_ FASTA path; descriptor permits subprocess input
	 * @param pairnum_ Pair marker (0 or 1), without linking mates
	 * @param maxReads_ Retained-record limit; negative means unlimited
	 */
	public FastaStreamer2ST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Configures a noninterleaved descriptor and queue without opening input.
	 * Initial Read flags reflect descriptor/global amino settings.
	 * @param ffin_ Nonnull input descriptor
	 * @param pairnum_ Pair marker (0 or 1), checked by assertion
	 * @param maxReads_ Retained-record limit; negative means unlimited
	 */
	public FastaStreamer2ST(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		flag=(ffin.amino() || Shared.AMINO_IN ? Read.AAMASK : 0);
		pairnum=pairnum_;
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(!interleaved) : "FastaStreamer2ST does not support interleaved files";
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		// Simple output queue
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);
		if(verbose){outstream.println("Made FastaStreamer2ST");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the parser worker and resets exposed counters.
	 * Call once: this precondition is not enforced, and terminal/error/queue state
	 * is not reset for a restart. Configure sampling and flags before calling.
	 */
	@Override
	public void start(){
		if(verbose){outstream.println("FastaStreamer2ST.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Start processing thread
		thread=new ProcessThread();
		thread.start();

		if(verbose){outstream.println("FastaStreamer2ST started.");}
	}

	/** Closes an established backend, retaining reported errors and exceptional cleanup.
	 * The reference clears only after normal return; callers may retry after a throw.
	 * Does not mark consumption finished or join the worker. Call after consuming
	 * the terminal result; this is not a cancellation or concurrent-close protocol.
	 */
	@Override
	public void close(){
		boolean completed=false;
		try{
			//[stream/FastaStreamer2ST#001] Preserve backend close/read errors on every close path.
			if(bf!=null){errorState|=bf.close(); bf=null;}
			completed=true;
		}finally{errorState|=!completed;}
	}

	/** Returns the configured input path. */
	@Override
	public String fname(){return fname;}

	/** Reports whether terminal consumption is pending, not whether data is queued. */
	@Override
	public boolean hasMore(){return !finished;}

	/** Returns cached backend, worker, consumer-interruption or cleanup failure status.
	 * This plain read is not a completion barrier or a general validation result.
	 */
	@Override
	public boolean errorState(){return errorState;}

	/** Returns false; this reader does not link mates. */
	@Override
	public boolean paired(){return false;}

	/** Returns the pair marker applied to retained reads. */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns retained-read totals copied when this consumer takes the terminal marker.
	 * Remains zero during ordinary data consumption. A failed batch or queue handoff
	 * can leave counted reads that were never delivered to the consumer.
	 */
	@Override
	public long readsProcessed(){return readsProcessed;}

	/** Returns the corresponding retained-base total; see readsProcessed(). */
	@Override
	public long basesProcessed(){return basesProcessed;}

	/** Sets the retention threshold and creates a fresh generator for sampled input.
	 * Configure before start and keep stable while the worker runs.
	 * @param rate Retention probability; at least 1 disables sampling, 0 discards all
	 * @param seed Nonnegative for reproducibility; negative selects a time-based seed
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	/** Takes the next nonempty data batch, or consumes/reinserts the terminal marker.
	 * Data arrays and lists remain caller-owned. Backend blocks determine batch
	 * sizes; IDs count blocks, including discarded ones, and firstRecordNum counts
	 * retained reads before that block. Terminal consumption copies worker totals
	 * and folds failure status before returning null.
	 * An interrupted take also returns null and sets errorState, without setting
	 * finished, copying counters or restoring the interrupt flag. Queue operations
	 * may block; terminal delivery is not a cancellation guarantee.
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

	/** Always throws UnsupportedOperationException; FASTA has no SamLine representation. */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns parsing state and publishes fresh lists, then attempts cleanup and termination. */
	private class ProcessThread extends Thread{

		/** Creates the named worker; input opens during processing. */
		ProcessThread(){setName("FastaStreamer2ST-Worker");}

		/** Runs parsing, catches Exceptions and attempts cleanup before terminal publication.
		 * Errors still pass through finally. Cleanup may replace a parser exception;
		 * blocking or interruption can prevent terminal delivery. An input constructor
		 * failure before assigning bf cannot be cleaned through that reference.
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
				//[stream/FastaStreamer2ST#002] Attempt cleanup before terminal publication, including failure.
				//Cleanup can replace a parser exception; blocking/interruption can prevent terminal delivery.
				try{close();}finally{
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

		/** Opens the record-block backend and publishes retained records until EOF/limit.
		 * Uses meaningful newline boundaries, not backing-array capacity. Header/base
		 * slices are copied before sampling; only retained records become Reads.
		 * Read constructor settings determine validation; no extra validate call is made.
		 * Counts advance before batch handoff; run handles cleanup on exit.
		 * @throws InterruptedException If publication of a data batch is interrupted
		 */
		void processSingle() throws InterruptedException{
			if(verbose){outstream.println("Started processSingle.");}

			bf=new ByteFile1Fc(ffin);
			IntList newlines=new IntList(256);
			//Correctness rests on the ByteFile1Fc contract: each block starts with '>', holds
			//complete records, and carries EXACTLY 2 '\n' per record (header-\n, sequence-\n) with
			//\r and internal sequence newlines stripped (it unwraps multi-line FASTA via SIMD).
			//That is why the 2-newline-per-record walk below and block[nl0+1]=='>' are correct;
			//wrapped FASTA is NOT a bug here - it is unwrapped upstream. A non-'>' start = corrupt
			//block -> the assert crashes loud (-ea) / best-effort (-da).

			long listNumber=0;
			while(readsProcessedT<maxReads){
				// Get block of records with newline positions
				byte[] block=bf.nextLine(newlines);
				if(block==null || block.length==0){break;}

				ArrayList<Read> readList=new ArrayList<Read>();
				ListNum<Read> reads=new ListNum<Read>(readList, listNumber++);
				reads.firstRecordNum=readsProcessedT;

				for(int i=0, nl0=-1; i<newlines.size() && readsProcessedT<maxReads; i++){
					int nl1=newlines.get(i);
					int nl2=(newlines.size()>i+1 ? newlines.get(i+1) : nl1);
					assert(block[nl0+1]=='>') : nl0+", "+(char)block[nl0+1];
					final byte[] header=KillSwitch.copyOfRange(block, nl0+2, nl1);
					final byte[] bases=(nl2>nl1 ? KillSwitch.copyOfRange(block, nl1+1, nl2) : null);
					
					if(samplerate>=1f || randy.nextFloat()<samplerate){
						Read r=new Read(bases, null, new String(header, 0, ReadHeader.end(header, 0, Shared.TRIM_READ_DESCRIPTION), StandardCharsets.US_ASCII), readsProcessedT, flag);
						r.setPairnum(pairnum);
						readList.add(r);
						readsProcessedT++;
						basesProcessedT+=r.length();
					}
					
					if(bases!=null){
						i++;
						nl0=nl2;
					}else{nl0=nl1;}
				}

				if(readList.size()>0){outputQueue.put(reads);}
			}

			if(verbose){outstream.println("Finished processSingle.");}
		}

		/** Retained Reads constructed and added to local batches, including unpublished ones. */
		protected long readsProcessedT=0;
		/** Summed lengths of the worker's counted retained Reads. */
		protected long basesProcessedT=0;
		/** True after normal parsing/handoff return; cleanup status is recorded separately. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;

	/** Primary input file */
	final FileFormat ffin;

	/** Data batches and terminal marker; capacity counts lists, not records or bytes. */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;

	/** Worker created by start; not joined by close. */
	private ProcessThread thread;
	
	/** Concrete backend, opened by the worker and cleared after normal close returns. */
	private ByteFile1Fc bf;

	/** Pair marker attached to retained reads. */
	final int pairnum;
	/** Descriptor setting, asserted false during construction. */
	final boolean interleaved;

	/** Worker retained-read count copied at terminal consumption. */
	protected long readsProcessed=0;
	/** Worker retained-base count copied at terminal consumption. */
	protected long basesProcessed=0;

	/** Retained-record limit, normalized to Long.MAX_VALUE for a negative argument. */
	final long maxReads;
	/** Read-constructor flags; configure before start and keep stable during parsing. */
	public int flag;

	/** Set when terminal list is received */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum queued lists, including the terminal marker; not a memory or read bound. */
	private static final int QUEUE_SIZE=2;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** Accumulated reported backend, worker, consumer-interruption and cleanup failures. */
	public boolean errorState=false;
	/** Retention threshold, configured before starting the worker. */
	private float samplerate=1f;
	/** Fresh seeded generator for sampling; null when sampling is disabled. */
	private shared.Random randy=null;

}
