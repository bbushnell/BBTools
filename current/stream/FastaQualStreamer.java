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
import structures.ByteBuilder;
import structures.ListNum;

/** Queued FASTA plus numeric QUAL reader with one parser worker.
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * ByteFile selection may add backend input workers. Configure before starting,
 * start once, consume batches to a terminal result, and explicitly close. This
 * class supplies neither a restart nor a concurrent cancellation protocol.
 *
 * Records follow file order; header names are not compared and extra QUAL records
 * are not checked. Sequence lines are concatenated without parser case conversion.
 * Scores use literal-space integer fields narrowed to bytes; this is not general
 * whitespace or numeric-range validation. Length equality is checked before
 * sampling, while only retained reads undergo configurable Read validation.
 * Retained header names use the JVM's default charset.
 *
 * Returned lists and copied base/quality arrays are not recycled by this class.
 * Terminal-marker consumption publishes processed counters and worker status.
 * An interrupted consumer wait also returns null, but does not finalize counters
 * or availability state. Check errorState as well as the terminal result.
 * @author Collei, Brian Bushnell, Shinobu
 * @date November 21, 2025
 */
public class FastaQualStreamer implements Streamer{

	/** Drains the reader, rejects reported failure and prints individual-read/base statistics.
	 * @param args FASTA path, QUAL path, and optional SIMD boolean subject to capability
	 * @throws RuntimeException If the reader reports failure
	 */
	public static void main(final String[] args){
		final Timer t=new Timer();
		final String faName=args[0];
		final String qualName=args[1];
		//[stream/FastaQualStreamer#005] Explicit true still requires the normal capability gate; false disables SIMD.
		if(args.length>2){Shared.SIMD=Parse.parseBoolean(args[2]) && simd.Vector.simd256;}

		final FileFormat ffFa=FileFormat.testInput(faName, FileFormat.FASTA, null, true, true);
		//Retain generic input detection for the diagnostic main's quality descriptor.
		final FileFormat ffQual=FileFormat.testInput(qualName, FileFormat.UNKNOWN, null, true, true);

		final FastaQualStreamer st=new FastaQualStreamer(ffFa, ffQual, 0, -1);
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
		//[stream/FastaQualStreamer#002] A terminal batch can also report worker failure.
		if(st.errorState()){throw new RuntimeException("FastaQualStreamer failed while reading input; see preceding diagnostic.");}
		t.stop();//[stream/FastaQualStreamer#006] Timer.elapsed stays zero until stop; rates otherwise become nonfinite.
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates descriptors and queue state without opening input or starting workers.
	 * The QUAL descriptor inherits the FASTA descriptor's subprocess permission.
	 * @param ff FASTA input descriptor
	 * @param qf Numeric QUAL input path
	 * @param pairnum_ Side marker, 0 or 1; does not create linked mates
	 * @param maxReads_ Maximum input records before sampling; negative means unlimited
	 */
	public FastaQualStreamer(final FileFormat ff, final String qf, final int pairnum_, final long maxReads_){
		this(ff, FileFormat.testInput(qf, FileFormat.QUAL, null, ff.allowSubprocess(), false),
			pairnum_, maxReads_);
	}

	/** Retains input descriptors and configuration; input opens in the parser worker.
	 * @param ffinFa_ FASTA descriptor
	 * @param ffinQual_ Numeric QUAL descriptor whose records correspond by order
	 * @param pairnum_ Side marker, 0 or 1; asserted when assigned to a retained Read
	 * @param maxReads_ Maximum input records before sampling; negative means unlimited
	 */
	public FastaQualStreamer(final FileFormat ffinFa_, final FileFormat ffinQual_, final int pairnum_, final long maxReads_){
		ffinFa=ffinFa_;
		ffinQual=ffinQual_;
		fname=ffinFa_.name();//Primary name is the FASTA file.
		pairnum=pairnum_;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);
		if(verbose){outstream.println("Made FastaQualStreamer for "+fname+" + "+ffinQual_.name());}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the single parser worker and resets exposed counters.
	 * Call once per instance. This does not reset the queue, finished flag or error
	 * state and therefore does not implement restart. Input/configuration globals
	 * and sampling settings must remain stable while processing.
	 */
	@Override
	public void start(){
		if(verbose){outstream.println("FastaQualStreamer.start() called.");}
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		//Start processing thread
		thread=new ProcessThread();
		thread.start();
		if(verbose){outstream.println("FastaQualStreamer started.");}
	}

	/** Attempts both input closes, retaining reported errors and marking exceptional cleanup.
	 * References clear only after normal close returns. If both closes throw, the QUAL
	 * exception replaces the FASTA exception; an uncleared reference may be retried later.
	 * Backend close may wait or throw. This is not a cancellation or concurrent-close
	 * protocol; ordinary callers close after consuming the worker's terminal result.
	 */
	@Override
	public void close(){
		boolean completed=false;
		try{
			try{
				//[stream/FastaQualStreamer#001,#004] Retain BOTH backend flags without clearing a prior error.
				if(bfFa!=null){errorState|=bfFa.close(); bfFa=null;}
			}finally{
				if(bfQual!=null){errorState|=bfQual.close(); bfQual=null;}
			}
			completed=true;
		}finally{errorState|=!completed;}
	}

	/** Returns the primary FASTA input name.
	 * @return Name captured from the FASTA descriptor
	 */
	@Override
	public String fname(){return fname;}

	/** Reports whether a terminal marker has not yet been consumed.
	 * True is possible before start, while the queue is empty, or after an interrupted
	 * wait. This is an availability hint, not a data or completion guarantee.
	 * @return Inverse of the terminal-consumption flag
	 */
	@Override
	public boolean hasMore(){return !finished;}

	/** Reports recorded parsing, consumer-wait or cleanup errors.
	 * Consult after terminal consumption for the ordinary completed-reader status;
	 * this getter does not join the worker or supply a visibility barrier itself.
	 * @return Currently recorded error flag
	 */
	@Override
	public boolean errorState(){return errorState;}

	/** Reports that this reader never links interleaved mates.
	 * @return false
	 */
	@Override
	public boolean paired(){return false;}

	/** Returns the side marker assigned to retained reads.
	 * @return Configured pair marker, normally 0 or 1
	 */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns the processed-record total copied at terminal-marker consumption.
	 * Fresh readers expose zero until then. Counts include sampled-out records and
	 * may include records in an unpublished partial batch. A record that fails parsing
	 * or retained-Read validation is not counted; a failed handoff can follow counting.
	 * @return Last published processed-record count
	 */
	@Override
	public long readsProcessed(){return readsProcessed;}

	/** Returns input bases for the records included in readsProcessed().
	 * Uses concatenated sequence lengths, before any retained-Read validation changes.
	 * @return Last published processed-base count
	 */
	@Override
	public long basesProcessed(){return basesProcessed;}

	/** Sets PRNG-based retention before start; rates at least one retain every record.
	 * Rates below one allocate a new generator. No rate validation or positional
	 * sampling is supplied; dropped records still undergo parsing and length checks,
	 * but bypass Read construction/validation. Do not reconfigure while processing.
	 * @param rate Requested retention fraction, normally in [0,1]
	 * @param seed Nonnegative deterministic seed, or negative for time-based seeding
	 */
	@Override
	public void setSampleRate(final float rate, final long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	/** Waits for the next retained-read batch after start.
	 * Data batch IDs increase from zero; firstRecordNum remains unspecified (-1).
	 * Taking terminal publishes worker counts/status, sets finished and reinserts
	 * terminal so later calls can also return null. An interrupted take instead
	 * records error and returns null without copying counts or setting finished;
	 * the interrupt flag is not restored here. Inspect errorState after null.
	 * @return Next retained batch, or null for terminal or an interrupted wait
	 */
	@Override
	public ListNum<Read> nextList(){
		try{
			final ListNum<Read> list=outputQueue.take();
			if(list==null || list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				errorState|=!thread.success;//[stream/FastaQualStreamer#001] |= not =, else this clobbers the
				//worker's reader-fold (errorState|=bfFa/bfQual.close()) when truncated input parses to EOF→success=true.
				if(list!=null){outputQueue.add(list);}//Re-inject terminal
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}

	/** Rejects SAM-line representation for this sequence reader.
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FastaQualStreamer does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns parser scratch state, input lookahead and generated counters. */
	private class ProcessThread extends Thread{

		/** Creates the named worker without starting input processing. */
		ProcessThread(){setName("FastaQualStreamer-Worker");}

		/** Parses input, catches Exception, and attempts cleanup then terminal publication.
		 * Errors are not caught, but still pass through finally. Cleanup exceptions
		 * can replace an active parser exception; blocking/interruption can prevent
		 * terminal delivery. A normal process return alone sets success.
		 */
		@Override
		public void run(){
			try{
				process();
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				//[stream/FastaQualStreamer#003] Cleanup precedes terminal publication, including parse/open failures.
				try{close();}finally{
					//Attempt publication even if close throws; blocking/interruption behavior is unchanged.
					try{
						ListNum<Read> terminal=new ListNum<Read>(null, -1, ListNum.LAST);
						outputQueue.put(terminal);
					}catch(InterruptedException e){
						e.printStackTrace();
					}
				}
			}
		}

		/** Opens both inputs and parses corresponding records before sampling.
		 * Snapshots batch targets, expected to be positive, without validating them.
		 * Leading nonheader lines are skipped independently; headers are never compared.
		 * FASTA controls termination, so extra
		 * QUAL records are unchecked. Numeric fields are split on literal spaces
		 * and narrowed to bytes without a general range check. Length mismatch fails
		 * before sampling. Only retained records are constructed and validated.
		 * Copied arrays and fresh data lists are handed off without recycling here.
		 * @throws InterruptedException If a data-batch queue insertion is interrupted
		 */
		void process() throws InterruptedException{
			bfFa=ByteFile.makeByteFile(ffinFa);
			bfQual=ByteFile.makeByteFile(ffinQual);

			long listNumber=0;
			long readID=0;
			int bytes=0;

			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<Read> ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);

			// Using LineParser to split space-delimited qualities
			final LineParser1 lp=new LineParser1(' ');
			final ByteBuilder bbBases=new ByteBuilder();
			final ByteBuilder bbQuals=new ByteBuilder();

			// Lookahead buffers
			byte[] nextFaLine=bfFa.nextLine();
			byte[] nextQualLine=bfQual.nextLine();

			while(readID<maxReads && nextFaLine!=null){

				// 1. Process Header
				final byte[] faHeader=nextFaLine;
				if(faHeader==null || faHeader.length==0 || faHeader[0]!='>'){
					// Skip garbage/blank lines between records if any, or sync issues
					if(nextFaLine==null){break;}//EOF
					nextFaLine=bfFa.nextLine();
					continue;
				}

				// 2. Process Bases
				bbBases.clear();
				while(true){
					nextFaLine=bfFa.nextLine();
					if(nextFaLine==null){break;}
					if(nextFaLine.length>0 && nextFaLine[0]=='>'){
						break;//Hold the next header in nextFaLine.
					}
					bbBases.append(nextFaLine);
				}

				// 3. Process Qual Header
				// Skip garbage until we find a header
				while(nextQualLine!=null && (nextQualLine.length==0 || nextQualLine[0]!='>')){
					nextQualLine=bfQual.nextLine();
				}

				if(nextQualLine==null){
					throw new RuntimeException("FASTA file has more records than QUAL file (EOF reached in QUAL).");
				}

				//For speed, records are assumed to correspond by order; header names are not checked.

				// 4. Process Qual Scores
				bbQuals.clear();
				while(true){
					nextQualLine=bfQual.nextLine();
					if(nextQualLine==null){break;}
					if(nextQualLine.length>0 && nextQualLine[0]=='>'){
						break;//Found next header.
					}

					// Parse line of space-delimited integers
					lp.set(nextQualLine);
					final int terms=lp.terms();
					for(int i=0; i<terms; i++){
						final int q=lp.parseInt(i);
						// Store as raw byte value (numeric Phred score)
						bbQuals.append((byte)q);
					}
				}

				// 5. Integrity Check
				if(bbBases.length!=bbQuals.length){
					throw new RuntimeException("Read "+readID+": Sequence length ("+bbBases.length+
						") != Quality length ("+bbQuals.length+").\nHeader: "+new String(faHeader));
				}

				if(samplerate>=1f || randy.nextFloat()<samplerate){
					// Fasta headers start with '>', strip it for the Read object
					final String id=new String(faHeader, 1, ReadHeader.end(faHeader, 1, Shared.TRIM_READ_DESCRIPTION)-1);
					final Read r=new Read(bbBases.toBytes(), bbQuals.toBytes(), id, readID);
					r.setPairnum(pairnum);
					if(!r.validated()){r.validate(true);}//Validate even when constructor validation is disabled.

					ln.add(r);
					bytes+=r.length();
				}
				readID++;

				readsProcessedT++;
				basesProcessedT+=bbBases.length;

				if(ln.size()>=slimit || bytes>=blimit){
					outputQueue.put(ln);
					ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);
					bytes=0;
				}
			}

			if(ln.size()>0){outputQueue.put(ln);}
		}

		/** Successfully processed input records, including sampled-out records. */
		protected long readsProcessedT=0;
		/** Concatenated input bases in the successfully processed records. */
		protected long basesProcessedT=0;
		/** True when parsing/handoff returns normally; cleanup failures also set outer errorState. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary FASTA filename captured at construction. */
	public final String fname;
	/** FASTA descriptor; opened by the worker. */
	final FileFormat ffinFa;
	/** Corresponding numeric QUAL descriptor; opened by the worker. */
	final FileFormat ffinQual;

	/** Bounded handoff of retained-read lists and the reinserted terminal marker. */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;
	/** Parser worker installed by start; not joined by close. */
	private ProcessThread thread;

	/** Worker-opened FASTA backend, cleared when its close returns normally. */
	private ByteFile bfFa;
	/** Worker-opened QUAL backend, cleared when its close returns normally. */
	private ByteFile bfQual;

	/** Side marker assigned to each retained unpaired read. */
	final int pairnum;

	/** Record count copied from the worker when terminal is consumed; initially zero. */
	protected long readsProcessed=0;
	/** Input-base count copied from the worker when terminal is consumed; initially zero. */
	protected long basesProcessed=0;

	/** Input-record limit before sampling; negative constructor limits normalize to MAX_VALUE. */
	final long maxReads;

	/** Set only by terminal-marker consumption, not by close or an interrupted take. */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained-record batch target, initially captured from Shared at class initialization.
	 * Configure a positive value before start; each worker snapshots it when processing begins.
	 */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Retained-base batch target; a record can exceed it. Configure positive before start. */
	public static int TARGET_LIST_BYTES=262144;
	/** Maximum queued lists, including a terminal marker while present. */
	private static final int QUEUE_SIZE=4;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination for optional verbose lifecycle messages. */
	protected PrintStream outstream=System.err;
	/** Enables optional lifecycle messages when recompiled as true. */
	public static final boolean verbose=false;
	/** Accumulated worker/consumer/cleanup status; no independent join or visibility barrier. */
	public boolean errorState=false;
	/** Retention fraction read during parsing; configure before start. */
	private float samplerate=1f;
	/** Per-setter PRNG for rates below one; unused/null when all records are retained. */
	private shared.Random randy=null;

}
