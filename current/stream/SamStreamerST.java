package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;
import structures.ListNum;

/** SAM parser with one background worker and a bounded queue of parsed line batches.
 * Construction does not start reading; call start once, then consume nextLines or
 * nextReads from one consumer. The ByteFile backend may use additional threads.
 * Sampling precedes parsing, so discarded records are not validated. Input limits
 * count nonheader records before sampling; reported counts include retained records.
 * Header publication and ByteFile selection involve JVM-wide settings. Configure
 * sampling before start; concurrent reconfiguration/lifecycle use is not supported.
 * @author Brian Bushnell, Isla, Shinobu
 * @date November 10, 2025
 */
public class SamStreamerST implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a SAM-fallback descriptor without content probing; does not start reading.
	 * @param fname_ Input path or standard-input name; subprocess input is permitted
	 * @param saveHeader_ Request shared-header publication
	 * @param maxReads_ Maximum input nonheader records before sampling; negative is unlimited
	 * @param makeReads_ Attach Read objects to parsed SamLine objects
	 */
	SamStreamerST(final String fname_, final boolean saveHeader_, final long maxReads_, final boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.SAM, null, true, false),
			saveHeader_, maxReads_, makeReads_);
	}

	/** Captures configuration and allocates the queue without opening the input.
	 * Header retention disabled here leaves any previous shared publication untouched.
	 * @param ffin_ Nonnull input descriptor passed to ByteFile when the worker starts
	 * @param saveHeader_ Collect and publish header lines through SamReadInputStream
	 * @param maxReads_ Maximum input nonheader records before sampling; negative is unlimited
	 * @param makeReads_ Build Read objects for use by nextReads/nextList
	 */
	SamStreamerST(final FileFormat ffin_, final boolean saveHeader_, final long maxReads_, final boolean makeReads_){
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

	/** Resets reported counters and launches the SAM worker. Call once per instance;
	 * this method does not reset the queue, header, error or terminal-consumed state.
	 */
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

	/** Closes the currently observed ByteFile and retains its reported error flag.
	 * Does not join the SAM worker or mark terminal consumption; backend close may
	 * wait for its own threads. This is not a concurrent teardown guarantee.
	 */
	@Override
	public synchronized void close(){
		//[SamStreamerST#001] Close/fold the field after an exception skips the worker's close.
		//Normal completion nulls bf, making this a no-op then; same pattern as FastaStreamerST.
		if(bf!=null){errorState|=bf.close(); bf=null;}
	}

	/** @return Input name captured from the descriptor */
	@Override
	public String fname(){return fname;}

	/** @return Whether the consumer has not yet consumed a terminal batch; not a queue-size check */
	@Override
	public boolean hasMore(){return !finished;}

	/** @return Currently observed error flag, without waiting for the worker */
	@Override
	public boolean errorState(){return errorState;}

	/** @return false; this streamer does not link records into paired Read batches */
	@Override
	public boolean paired(){return false;}

	/** @return Zero; there is no separate mate-input index for this streamer */
	@Override
	public int pairnum(){return 0;}

	/** @return Retained SAM record count copied from the worker on terminal consumption */
	@Override
	public long readsProcessed(){return readsProcessed;}

	/** @return Retained sequence-base count copied from the worker on terminal consumption */
	@Override
	public long basesProcessed(){return basesProcessed;}

	/** Configures sequential random sampling before parsing; call before start.
	 * Rates at least one disable sampling. Otherwise uses a new seeded generator
	 * from Shared.threadLocalRandom, not Streamer's positional sampling helper.
	 * @param rate Keep threshold compared with nextFloat; values are not validated here
	 * @param seed Generator seed passed to Shared.threadLocalRandom
	 */
	@Override
	public void setSampleRate(final float rate, final long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	/** @return Next Read batch from nextReads, or null on terminal/interruption */
	@Override
	public ListNum<Read> nextList(){return nextReads();}

	/** Takes the next parsed batch, waiting for the worker as needed.
	 * Terminal consumption marks finished, copies retained counters and folds worker
	 * failure into errorState, then reinserts the terminal for subsequent calls.
	 * Interruption sets errorState and returns null without those terminal actions
	 * or restoring the interrupt flag. Returned data batches are not copied.
	 * @return Parsed line batch, or null on terminal/interruption
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		try{
			ListNum<SamLine> list=outq.take();
			assert(list!=null) : "Pulled null list.";//Should never happen
			if(verbose){
				if(list==null || list.last()){outstream.println("Consumer got terminal list.");}
				else{outstream.println("Consumer got list "+list.id());}
			}
			if(list==null || list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				//errorState-fold-clobber [family sweep 2026-06-22; FastaStreamerST#004]: |= preserves
				//bf.close errors even when normal return sets success=true. Callers must check the flag.
				errorState|=!thread.success;
				outq.add(list);//Re-inject
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}

	/** Consumes a parsed batch and collects its attached Read objects without copying them.
	 * Allocates a new list and ListNum wrapper. Requires construction with makeReads=true;
	 * numeric IDs follow input positions and may have gaps after sampling.
	 * @return Read batch, or null when nextLines returns null
	 */
	public ListNum<Read> nextReads(){
		assert(makeReads);
		ListNum<SamLine> lines=nextLines();
		if(lines==null){return null;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		if(!lines.isEmpty()){
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

	/** Parses input and publishes batches through the enclosing instance's queue. */
	private class ProcessThread extends Thread{

		/** Sets the existing diagnostic thread name; does not start the worker. */
		ProcessThread(){
			setName("SamLineStreamerST-Worker");
		}

		/** Runs the parser, reports caught Exceptions, then attempts terminal publication.
		 * success records normal parser return, independently of reported close errors.
		 * Finally also runs for uncaught Errors; an interrupted terminal put is printed
		 * and not retried. Terminal consumption elsewhere folds !success into errorState.
		 */
		@Override
		public void run(){
			try{
				processFileDirectly();
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				//COMPREHENSION [OQS-death-safe; FastaStreamerST]: finally attempts a terminal even
				//after an uncaught Error. Delivery still depends on put completing without interruption.
				try{
					ListNum<SamLine> terminal=new ListNum<SamLine>(null, -1, false, true);
					outq.put(terminal);
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("ProcessThread terminated.");}
		}

		/** Opens ByteFile, parses retained records and publishes batches/header on normal progress.
		 * Uses a local BF2 preference, honoring global force flags without changing them.
		 * Blank lines are skipped;
		 * lines beginning with {@code @} are headers. Requested publication occurs before the first
		 * nonheader record or on normal loop exit. The input limit counts nonheader records
		 * before sampling; loop initialization/update may fetch beyond that logical limit.
		 * Batch byte limits estimate twice retained sequence length, not input SAM bytes.
		 * Normal exit publishes any remaining batch and folds/closes the ByteFile.
		 * @throws InterruptedException If a data-batch queue put is interrupted
		 */
		void processFileDirectly() throws InterruptedException{
			if(verbose){outstream.println("Started processFileDirectly.");}

			//[SamStreamerST#003; STR-023] Prefer BF2 for this reader without changing later readers' defaults.
			//Explicit force flags retain their priority; an unavailable preference uses automatic selection.
			bf=ByteFile.makeByteFileWithPreference(ffin, 2);//[SamStreamerST#001] use the bf FIELD (not a local) so the implemented close() can reach it on the worker-exception path; matches FastaStreamerST#002

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
					if(header!=null){
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

		/** Retained records parsed by this worker; excludes sampled-out records. */
		protected long readsProcessedT=0;
		/** Sequence bases in retained records; absent sequences contribute zero. */
		protected long basesProcessedT=0;
		/** True after normal processFileDirectly return; reported close errors may still exist. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input name captured from the descriptor. */
	public final String fname;

	/** Input descriptor used when the worker opens ByteFile. */
	final FileFormat ffin;

	/** Bounded data/terminal queue; terminal consumers reinsert the marker. */
	final ArrayBlockingQueue<ListNum<SamLine>> outq;
	/** Worker created by start; counters are read after terminal consumption. */
	private ProcessThread thread;
	/** Input source — a FIELD so close() can reach it on the worker-exception path (#001) */
	private ByteFile bf;
	/** Consumer-side terminal-consumed flag; close/interruption does not set it. */
	private boolean finished=false;

	/** Captured header request, used to allocate header storage at construction. */
	final boolean saveHeader;
	/** Whether each retained SamLine receives an attached validated Read. */
	final boolean makeReads;

	/** Header under construction, then null after publication; null when retention is disabled. */
	ArrayList<byte[]> header;

	/** Retained record count, updated on terminal consumption rather than continuously. */
	protected long readsProcessed=0;
	/** Retained base count, updated on terminal consumption rather than continuously. */
	protected long basesProcessed=0;

	/** Limit on input nonheader records before sampling; negative requests become Long.MAX_VALUE. */
	final long maxReads;

	/** Stream for optional verbose messages; exception traces use their default stream. */
	protected PrintStream outstream=System.err;
	/** Reported input, worker or consumer-interruption error; callers must inspect it. */
	public boolean errorState=false;
	/** Sampling threshold configured before start. */
	private float samplerate=1f;
	/** Sequential generator for rates below one; null disables sampling. */
	private shared.Random randy=null;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Batch record threshold, copied from Shared at class initialization. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Batch threshold for twice retained sequence length, copied from Shared at class initialization. */
	public static int TARGET_LIST_BYTES=shared.Shared.bufferSize();
	/** Maximum number of queued parsed batches/terminal markers. */
	private static final int QUEUE_SIZE=8;
	/** Compile-time verbose diagnostics switch. */
	public static final boolean verbose=false;

}
