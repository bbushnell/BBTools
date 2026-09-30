package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import structures.ByteBuilder;
import structures.ListNum;

/** FASTA reader with one parsing worker and a bounded queue of Read batches.
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * The selected ByteFile backend may add input threads. StreamerFactory uses this
 * worker fallback where its zero-worker and SIMD alternatives are not selected.
 * Configure before start, start once on a fresh instance, consume to terminal and close.
 * Closure does not interrupt or join the worker, and restart is not implemented.
 *
 * Headers are decoded as ASCII; sequence lines are concatenated into Reads with
 * null qualities and configured Read validation. Limits count fragments: individual
 * unpaired reads or completed interleaved pairs. Worker counters count individual reads.
 * Sampling rejects headers before construction. Counter getters update when a
 * terminal marker is consumed, and failed processing may count unpublished records.
 * The worker attempts cleanup and then terminal publication in finally, including
 * after an assertion failure. Cleanup or queue insertion can block; publication
 * can be interrupted. Cleanup exceptions may replace an active parsing exception.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @contributor Shinobu (documentation, formatting and interleaved limits)
 * @date November 12, 2025
 */
public class FastaStreamerST implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an unstarted reader using inferred FASTA format settings.
	 * @param fname_ Input filename
	 * @param pairnum_ 0 for unpaired/R1, or 1 for separate R2 input
	 * @param maxReads_ Fragment limit (unpaired reads or interleaved pairs); negative for unlimited
	 */
	public FastaStreamerST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Captures format, pairing, amino flags and the fragment limit.
	 * Input is opened later by the worker. Interleaved limits stop after complete pairs.
	 * With assertions enabled, odd input before the limit fails the incomplete-pair assertion.
	 * @param ffin_ Input descriptor, including interleaving/amino settings
	 * @param pairnum_ 0 for unpaired/interleaved/R1, or 1 for separate R2 input
	 * @param maxReads_ Fragment limit (unpaired reads or interleaved pairs); negative for unlimited
	 */
	public FastaStreamerST(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		flag=(ffin.amino() || Shared.AMINO_IN ? Read.AAMASK : 0);
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(pairnum==0 || !interleaved);
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		// Simple output queue
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);

		if(verbose){outstream.println("Made FastaStreamerST");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resets outer counters and launches a parsing worker; call once on a fresh instance.
	 * Repeated calls launch additional workers without resetting the queue or finished flag.
	 */
	@Override
	public void start(){
		if(verbose){outstream.println("FastaStreamerST.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Start processing thread
		thread=new ProcessThread();
		thread.start();

		if(verbose){outstream.println("FastaStreamerST started.");}
	}

	/** Closes the currently available ByteFile, folding reported errors and latching thrown failures.
	 * The reference clears only after normal return, allowing a later retry after a throw.
	 * Does not interrupt/join the worker, prevent later opening, publish a terminal
	 * marker or mark consumption finished. Prefer calling after terminal consumption.
	 */
	@Override
	public void close(){
		boolean completed=false;
		try{
			//[stream/FastaStreamerST#005] Keep returned and thrown cleanup failures before terminal publication.
			if(bf!=null){errorState|=bf.close(); bf=null;}
			completed=true;
		}finally{errorState|=!completed;}
	}

	/** Returns the captured input name. */
	@Override
	public String fname(){return fname;}

	/** Reports whether normal terminal consumption is still pending; not a queue inspection. */
	@Override
	public boolean hasMore(){
		return !finished;
	}

	/** Returns the accumulated worker, reader-close or consumer-interruption status. */
	@Override
	public boolean errorState(){return errorState;}

	/** Returns the interleaving setting captured at construction. */
	@Override
	public boolean paired(){return interleaved;}

	/** Returns the configured input-side pair number. */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns constructed individual reads copied from the worker on terminal consumption.
	 * Counts both mates and may include records not published when processing failed.
	 */
	@Override
	public long readsProcessed(){return readsProcessed;}

	/** Returns constructed bases copied from the worker on terminal consumption. */
	@Override
	public long basesProcessed(){return basesProcessed;}

	/** Configures PRNG header selection before the worker starts.
	 * Every interleaved rate below one, including zero, is rejected with assertions
	 * enabled. With assertions disabled, that mode does not guarantee pair alignment.
	 * This is not positional Streamer.sampleKeep sampling.
	 * @param rate Requested retention threshold, normally in [0,1]
	 * @param seed Seed passed to Shared.threadLocalRandom when sampling is active
	 * @throws AssertionError If assertions are enabled and interleaved rate is below one
	 */
	@Override
	public void setSampleRate(float rate, long seed){
		//[stream/FastaStreamerST#001] Crash-loud guard (Brian 2026-06-22): fractional sampling of INTERLEAVED FASTA desyncs read pairs (the start-read1 roll has no file-parity guard). Unsupported weird corner — best-effort only. Under -ea (default) crash loud with the workaround; under -da the determined user gets silent best-effort. No hot-loop parity fix for a path nobody uses; crash-don't-corrupt.
		assert(!(interleaved && rate<1f)) : "Fractional sampling of interleaved FASTA is unsupported (read pairs would desync). Workaround: convert to FASTQ, subsample, then convert back to FASTA. ["+fname+"]";
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	/** Takes a batch, or consumes/reinserts the terminal and copies worker counters/status.
	 * Normal terminal reinsertion supports repeated EOF calls. Consumer interruption
	 * returns null with an error flag, without marking normal terminal consumption.
	 * @return Next data batch, or null at terminal input or consumer interruption
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
				errorState|=!thread.success;//[stream/FastaStreamerST#004] |= not =: a plain '=' OVERWRITES the worker's reader-error fold (errorState|=bf.close() fires on truncated/corrupt input that still parsed to completion → success=true), silently dropping the truncation signal and defeating the 2b reader-fold. OR-ing keeps BOTH the reader fold and a worker-thread failure.
				outputQueue.add(list);//Re-inject
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}

	/** Rejects the unsupported SAM-line view.
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses selected headers/sequences and publishes batches to the bounded output queue. */
	private class ProcessThread extends Thread{

		/** Creates a named worker; input opening occurs during run. */
		ProcessThread(){
			setName("FastaStreamerST-Worker");
		}

		/** Processes the selected mode, then attempts cleanup and terminal publication.
		 * Errors also pass through finally. Cleanup may replace a parsing exception;
		 * blocked cleanup or an interrupted/blocked put can prevent terminal delivery.
		 * success records normal parsing/handoff return; cleanup status is separate.
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
				//[stream/FastaStreamerST#005] Both parsing modes need cleanup even after failure.
				//Historical worker-error rationale: finally also runs after an AssertionError, which catch(Exception) misses.
				//If the terminal is delivered, nextList folds !success without erasing prior reader errors.
				//Keep a terminal attempt even when cleanup throws; neither operation guarantees completion.
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

		/** Concatenates sequence lines for selected headers and publishes unpaired Reads.
		 * Rejected headers do not count toward the construction limit or counters. IDs count
		 * constructed reads; loop transitions may fetch beyond the final retained record.
		 * run() attempts backend cleanup after this method returns or throws.
		 * @throws InterruptedException If batch publication is interrupted
		 */
		void processSingle() throws InterruptedException{
			if(verbose){outstream.println("Started processSingle.");}

			bf=ByteFile.makeByteFile(ffin);

			long listNumber=0;
			int readsInList=0;
			int bytesInList=0;

			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<Read> ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);
			ln.firstRecordNum=readsProcessedT;

			final ByteBuilder bb=new ByteBuilder(4096);
			byte[] header=null;
			byte[] line=null;

			for(line=bf.nextLine(); line!=null && readsProcessedT<maxReads; line=bf.nextLine()){

				if(line.length>0 && line[0]=='>'){
					if(header!=null){
						Read r=new Read(bb.toBytes(), null, new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, 
							StandardCharsets.US_ASCII), readsProcessedT, flag);
						r.setPairnum(pairnum);
						ln.add(r);
						readsProcessedT++;
						basesProcessedT+=r.length();
						readsInList++;
						bytesInList+=r.length();
					}
					header=null;
					bb.clear();

					if(samplerate>=1f || randy.nextFloat()<samplerate){header=line;}
					
					if(readsInList>=slimit || bytesInList>=blimit){
						outputQueue.put(ln);
						ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);
						ln.firstRecordNum=readsProcessedT;
						readsInList=0;
						bytesInList=0;
					}
				}else if(header!=null){
					bb.append(line);
				}
			}

			//Do not flush a lookahead header after the requested read limit has been reached.
			if(line==null && header!=null && readsProcessedT<maxReads){
				Read r=new Read(bb.toBytes(), null, 
					new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, StandardCharsets.US_ASCII), readsProcessedT, flag);
				r.setPairnum(pairnum);
				ln.add(r);
				readsProcessedT++;
				basesProcessedT+=r.length();
			}

			if(ln.size()>0){
				outputQueue.put(ln);
			}
			if(verbose){outstream.println("Finished processSingle.");}
		}

		/** Constructs adjacent records, links pairs and publishes first-mate entries.
		 * Names are not checked for pair agreement. Limits count pairs; counters count
		 * individual reads. With assertions enabled, odd input before the limit fails the final assertion.
		 * Supported interleaved use requires samplerate at least one under assertions.
		 * run() attempts backend cleanup after this method returns or throws.
		 * @throws InterruptedException If batch publication is interrupted
		 */
		void processInterleaved() throws InterruptedException{
			if(verbose){outstream.println("Started processInterleaved.");}

			bf=ByteFile.makeByteFile(ffin);//Use the field so run's finally/close can reach it after the odd-pair assertion; matches processSingle.

			long listNumber=0;
			long readID=0;
			int readsInList=0;
			int bytesInList=0;

			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<Read> ln=new ListNum<Read>(new ArrayList<Read>(slimit/2), listNumber++);
			ln.firstRecordNum=readID;

			final ByteBuilder bb=new ByteBuilder(4096);
			byte[] header=null;
			Read pending=null;
			byte[] line=null;

			for(line=bf.nextLine(); line!=null && readID<maxReads; line=bf.nextLine()){

				if(line.length>0 && line[0]=='>'){
					if(header!=null){
						// Finish current read
						Read r=new Read(bb.toBytes(), null, 
							new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, StandardCharsets.US_ASCII), 0, flag);
						readsProcessedT++;
						basesProcessedT+=r.length();
						bb.clear();

						if(pending==null){
							// This is read1
							pending=r;
							pending.setPairnum(0);
							header=line; // Start building read2
						}else{
							// This is read2
							r.setPairnum(1);
							pending.mate=r;
							r.mate=pending;
							pending.numericID=readID;
							r.numericID=readID++;

							ln.add(pending);
							readsInList+=2;
							bytesInList+=pending.length()+r.length();
							pending=null;

							// Decide whether to start building next pair
							header=(samplerate>=1f || randy.nextFloat()<samplerate) ? line : null;

							// Check if we should ship current list
							if(readsInList>=slimit || bytesInList>=blimit){
								outputQueue.put(ln);
								ln=new ListNum<Read>(new ArrayList<Read>(slimit/2), listNumber++);
								ln.firstRecordNum=readID;
								readsInList=0;
								bytesInList=0;
							}
						}
					}else{
						//TODO: Possible bug under -da [stream/FastaStreamerST#001] - sampling can desynchronize interleaved pairs.
						//The setSampleRate assertion normally rejects this unsupported mode; retain that guard and its rationale.
						//COMPREHENSION [the deviation]: this start-read1 roll has NO file-parity guard. Under sampling (samplerate<1), if a pair's read1 is rejected (header stays null) and THIS branch then accepts the very next '>' — that pair's read2 — read2 is labeled read1 (pending, pairnum 0) and mate-paired with the NEXT pair's read1; mate-links+pairnums desync for the rest of the file. Under samplerate>=1f (common case) header is always set, strict read1/read2 alternation holds → no desync. `pending` tracks mid-pair but NOT file parity.
						// Not currently building - decide whether to start
						header=(samplerate>=1f || randy.nextFloat()<samplerate) ? line : null;
						bb.clear();
					}
				}else if(header!=null){
					// Accumulate bases
					bb.append(line);
				}
			}

			//Flush a final mate only before the pair limit; a lookahead header may belong to the next pair.
			if(line==null && header!=null && readID<maxReads){
				Read r=new Read(bb.toBytes(), null, 
					new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1, StandardCharsets.US_ASCII), 0, flag);
				readsProcessedT++;
				basesProcessedT+=r.length();

				if(pending==null){
					// File ends on read1 - incomplete pair
					pending=r;
					pending.setPairnum(0);
					ln.add(pending);
				}else{
					// This is read2 - complete the pair
					r.setPairnum(1);
					pending.mate=r;
					r.mate=pending;
					pending.numericID=readID;
					r.numericID=readID++;
					ln.add(pending);
					pending=null;
				}
			}

			// Handle incomplete pair at end
			//Historical intentional dual-mode behavior: under -ea a pending mate fails loudly; under -da EOF may
			//publish a lone read best-effort. Here AssertionError escapes catch(Exception), but run's finally attempts
			//terminal publication. On successful delivery, nextList records !success; interrupted publication is not guaranteed.
			assert(pending==null) : "Odd number of reads in interleaved FASTA file: "+fname;

			if(ln.size()>0){
				outputQueue.put(ln);
			}
			if(verbose){outstream.println("Finished processInterleaved.");}
		}

		/** Individual Reads constructed by this worker, including both mates and possibly unpublished data. */
		protected long readsProcessedT=0;
		/** Bases in constructed Reads, possibly including unpublished data on failure. */
		protected long basesProcessedT=0;
		/** True after normal parsing/handoff return; cleanup status is recorded separately. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input name captured at construction. */
	public final String fname;

	/** Input descriptor used when the worker opens a ByteFile. */
	final FileFormat ffin;

	/** Bounded handoff of data batches and a terminal marker. */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;

	/** Worker assigned by start; a fresh instance is required for each lifecycle. */
	private ProcessThread thread;
	
	/** Worker-opened reader, cleared after normal closure; also reached by outer close. */
	private ByteFile bf;

	/** Configured input-side pair marker; interleaved input requires zero under assertions. */
	final int pairnum;
	/** Whether adjacent constructed records are paired. */
	final boolean interleaved;

	/** Constructed individual-read total copied when a terminal marker is consumed. */
	protected long readsProcessed=0;
	/** Constructed base total copied when a terminal marker is consumed. */
	protected long basesProcessed=0;

	/** Fragment limit: constructed unpaired reads or completed interleaved pairs. */
	final long maxReads;
	/** Flags passed to Read construction, initially including the amino flag where requested. */
	public int flag;

	/** Set by normal terminal consumption, not by close or consumer interruption. */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructed individual-read batching threshold, initialized from Shared. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Constructed-base batching threshold, excluding headers and object overhead. */
	public static int TARGET_LIST_BYTES=262144;
	/** Maximum number of queued batches/terminal markers. */
	private static final int QUEUE_SIZE=2;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination for optional verbose messages; caught exceptions use their default stream. */
	protected PrintStream outstream=System.err;
	/** Compile-time optional diagnostics switch. */
	public static final boolean verbose=false;
	/** Accumulated worker, reader-close or consumer-interruption failure. */
	public boolean errorState=false;
	/** Header-retention threshold configured before start. */
	private float samplerate=1f;
	/** PRNG for rates below one; assigned by setSampleRate before start. */
	private shared.Random randy=null;

}
