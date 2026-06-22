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

/**
 * Single-thread (ST) FASTA file loader: spawns exactly ONE worker ProcessThread that reads+parses the whole file
 * and hands ListNum&lt;Read&gt; chunks to the consumer over a bounded ArrayBlockingQueue (QUEUE_SIZE=2) — cf.
 * FastaStreamerZT, which is host-driven and spawns no worker. Selected by StreamerFactory for threads==1 / cores&lt;8
 * on fasta input (not the SIMD path). The worker's run() pushes its terminal list in a FINALLY, so a worker-thread
 * death — including an AssertionError from the odd-interleaved assert, which catch(Exception) does not catch — still
 * delivers a terminal to the consumer: structurally immune to the OQS-thread-death-hang class despite having a worker.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 12, 2025
 */
public class FastaStreamerST implements Streamer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public FastaStreamerST(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/** Constructor. */
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

	@Override
	public void close(){
		if(bf!=null) {errorState|=bf.close(); bf=null;}//Fold the reader's error state (matches ZT + the reader-fold family); only bites if the worker threw before its own close, leaving the field bf open
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
	public boolean paired(){return interleaved;}

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public long readsProcessed() {return readsProcessed;}

	@Override
	public long basesProcessed() {return basesProcessed;}

	@Override
	public void setSampleRate(float rate, long seed){
		//[stream/FastaStreamerST#001] Crash-loud guard (Brian 2026-06-22): fractional sampling of INTERLEAVED FASTA desyncs read pairs (the start-read1 roll has no file-parity guard). Unsupported weird corner — best-effort only. Under -ea (default) crash loud with the workaround; under -da the determined user gets silent best-effort. No hot-loop parity fix for a path nobody uses; crash-don't-corrupt.
		assert(!(interleaved && rate<1f)) : "Fractional sampling of interleaved FASTA is unsupported (read pairs would desync). Workaround: convert to FASTQ, subsample, then convert back to FASTA. ["+fname+"]";
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

	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private class ProcessThread extends Thread {

		/** Constructor */
		ProcessThread(){
			setName("FastaStreamerST-Worker");
		}

		/** Called by start() */
		@Override
		public void run(){
			try{
				if(interleaved) {
					processInterleaved();
				}else {
					processSingle();
				}
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				//COMPREHENSION [OQS-death-safe — the headline]: this terminal push is in a FINALLY, not the catch — it runs on EVERY exit, including an AssertionError from the odd-interleaved assert in processInterleaved (an Error, which catch(Exception) above does NOT catch). So whatever kills the worker, the consumer's outputQueue.take() in nextList() always gets this terminal and reads errorState=!success → a worker death CANNOT strand the consumer. The structurally-correct OQS-death-hang fix, via a bounded queue.
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
					if(header!=null) {
						Read r=new Read(bb.toBytes(), null, new String(header, 1, header.length-1, 
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
					
					if(readsInList>=slimit || bytesInList>=blimit) {
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

			// Handle EOF
			if(line==null && header!=null) {
				Read r=new Read(bb.toBytes(), null, 
					new String(header, 1, header.length-1, StandardCharsets.US_ASCII), readsProcessedT, flag);
				r.setPairnum(pairnum);
				ln.add(r);
				readsProcessedT++;
				basesProcessedT+=r.length();
			}

			if(ln.size()>0){
				outputQueue.put(ln);
			}
			errorState|=bf.close(); bf=null;//Fold the reader's error state, then null the field so a later close() can't double-close (close() then no-ops cleanly)
			if(verbose){outstream.println("Finished processSingle.");}
		}

		void processInterleaved() throws InterruptedException{
			if(verbose){outstream.println("Started processInterleaved.");}

			bf=ByteFile.makeByteFile(ffin);//Use the bf FIELD (not a local) so close() can reach it on the -ea odd-interleaved crash path (the assert below throws before this method's own close); matches processSingle

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

			for(line=bf.nextLine(); line!=null && readsProcessedT<maxReads; line=bf.nextLine()){

				if(line.length>0 && line[0]=='>'){
					if(header!=null){
						// Finish current read
						Read r=new Read(bb.toBytes(), null, 
							new String(header, 1, header.length-1, StandardCharsets.US_ASCII), 0, flag);
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
						//TODO: Possible bug [stream/FastaStreamerST#001] - interleaved+sampling pair desync (LOW; shared with FastaStreamerZT; FLOAT: is interleaved-FASTA fractional sampling an intended workflow?)
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

			// Handle EOF
			if(line==null && header!=null) {
				Read r=new Read(bb.toBytes(), null, 
					new String(header, 1, header.length-1, StandardCharsets.US_ASCII), 0, flag);
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
			//COMPREHENSION [INTENTIONAL dual-mode, per Brian — NOT a bug; same pattern as FastaStreamerZT]: the EOF block above populates the lone read1 of an odd interleaved file AND this assert guards it, deliberately. Under -ea (default) it CRASHES LOUD; under -da it is skipped and the lone read ships best-effort. ST-vs-ZT wrinkle: here it runs on the WORKER thread, so the AssertionError (an Error) escapes catch(Exception) in run() — but run()'s FINALLY still delivers the terminal, so the consumer ends clean (errorState=true), not stranded.
			assert(pending==null) : "Odd number of reads in interleaved FASTA file: "+fname;

			if(ln.size()>0){
				outputQueue.put(ln);
			}
			errorState|=bf.close(); bf=null;//Fold the reader's error state, then null the field so a later close() can't double-close (close() then no-ops cleanly)
			if(verbose){outstream.println("Finished processInterleaved.");}
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

	/** Output queue */
	final ArrayBlockingQueue<ListNum<Read>> outputQueue;

	/** Processing thread */
	private ProcessThread thread;
	
	/** Input source */
	private ByteFile bf;

	final int pairnum;
	final boolean interleaved;

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads */
	final long maxReads;
	public int flag;

	/** Set when terminal list is received */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	public static int TARGET_LIST_BYTES=262144;
	private static final int QUEUE_SIZE=2;

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