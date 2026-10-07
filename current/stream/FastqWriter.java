package stream;

import java.io.OutputStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Writes FASTQ, FASTA, SCARF, headers, or attachments with parallel formatting.
 * 
 * Workers format selected Read mates; a separate output thread writes blocks
 * ordered by batch ID through OrderedQueueSystem2. IDs must be unique and dense,
 * starting at zero, even when producers submit batches out of order.
 * <p>The constructor opens the output stream. Call start once, submit all batches,
 * then poisonAndWait. Submitted lists and their reads must remain unchanged until
 * writing completes. Null list entries are skipped; empty batches still preserve
 * their position in the ID sequence. Normal completion finalizes output and
 * publishes counters; error completion may abandon queued work. Bare standard
 * output/error streams are flushed and left open by ReadWrite.finishWriting.
 * <p>R1 is an entry with pairnum 0. R2 is either an entry with pairnum 1 or the
 * linked mate of an R1 entry. Each enabled mate is emitted independently. Header
 * and attachment output count emitted reads but do not add to the base counter.
 * 
 * @author Isla
 * @date October 30, 2025
 */
public class FastqWriter implements Writer{
	
	/**
	 * Copies sequence input through the streamer/writer factories and reports counts.
	 * @param args Input, optional output, formatting threads, presence-based SIMD request and BGZF flag
	 */
	public static void main(final String[] args){
		Timer t=new Timer();
		String in=args[0];
		String out=(args.length<2 || args[1].equalsIgnoreCase("null") ? null : args[1]);
		int threads=DEFAULT_THREADS;
		if(args.length>2){threads=Integer.parseInt(args[2]);}
		//STR-006: Preserve the argument-presence trigger while honoring JVM/hardware capability.
		if(args.length>3){Shared.SIMD=simd.Vector.simd256;}
		if(args.length>4){
			ReadWrite.ALLOW_NATIVE_BGZF=ReadWrite.PREFER_NATIVE_BGZF_IN=
				ReadWrite.PREFER_NATIVE_BGZF_OUT=Parse.parseBoolean(args[4]);
		}

//		ByteFile.FORCE_MODE_BF4=false;
		FileFormat ffin=FileFormat.testInput(in, FileFormat.FASTQ, null, true, true);
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.FASTQ, null, true, true, false, true);
		
		SamLine.SET_FROM_OK=ffin.samOrBam();
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=(ffout!=null && ffout.samOrBam() && ffin.samOrBam());

		final boolean inputReads=(!ffin.samOrBam());
		final boolean outputReads=(ffout!=null && !ffout.samOrBam());
		
		Streamer st=StreamerFactory.makeStreamer(ffin, 0, true, -1, true, outputReads);
		Writer fw=WriterFactory.makeWriter(ffout, true, true, threads, null, ffin.samOrBam());
//		assert(false) : "\n"+ffin+"\n"+ffout+"\n"+st.getClass()+"\n"+fw.getClass();
		process(st, fw, t, inputReads || outputReads);
	}
	
	/** Drains reads or SAM lines, finishes optional output, and prints processing counts. */
	private static void process(final Streamer st, final Writer fw, final Timer t, final boolean readMode){
		st.start();
		if(fw!=null){fw.start();}
		long reads=0, bases=0;
		if(readMode){
			for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
				for(Read r : ln){
					reads+=r.pairCount();
					bases+=r.pairLength();
				}
				if(fw!=null){fw.addReads(ln);}
			}
		}else{
			for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()){
				for(SamLine sl : ln){
					reads++;
					bases+=sl.length();
				}
				if(fw!=null){fw.addLines(ln);}
			}
		}
		if(fw!=null){
			fw.poisonAndWait();
			assert(reads==fw.readsWritten());
			assert(bases==fw.basesWritten());
			reads=fw.readsWritten();
			bases=fw.basesWritten();
		}
		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Opens a named output, inferring its format with FASTQ as the default.
	 * @param out_ Output filename
	 * @param threads_ Formatting workers; values below one use DEFAULT_THREADS
	 * @param writeR1_ Emit entries marked as R1, including unpaired reads
	 * @param writeR2_ Emit entries marked as R2 or linked R2 mates
	 * @param overwrite Allow replacement of an existing file
	 */
	public FastqWriter(final String out_, final int threads_, final boolean writeR1_,
		final boolean writeR2_, final boolean overwrite){
		this(FileFormat.testOutput(out_, FileFormat.FASTQ, null, true, overwrite, false, true), 
			threads_, writeR1_, writeR2_);
	}
	
	/**
	 * Opens output and creates queues; threads are started separately by start().
	 * The format must be FASTQ, FASTA, SCARF, HEADER or ATTACHMENT; UNKNOWN uses FASTQ.
	 * @param ffout_ Nonnull output descriptor; its append flag is honored
	 * @param threads_ Requested formatting workers, clamped to 1..Shared.threads();
	 * values below one first select DEFAULT_THREADS
	 * @param writeR1_ Emit R1 entries; at least one mate selection must be enabled
	 * @param writeR2_ Emit R2 entries or linked R2 mates
	 */
	public FastqWriter(final FileFormat ffout_, final int threads_,
		final boolean writeR1_, final boolean writeR2_){
		ffout=ffout_;
		fname=ffout_.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		writeR1=writeR1_;
		writeR2=writeR2_;
		format=(ffout.format()==UNKNOWN ? FASTQ : ffout.format());
		assert(format==FASTQ || format==FASTA || format==HEADER || format==SCARF || format==ATTACHMENT) : ffout;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Create OQS
		FastqWriterInputJob inputProto=new FastqWriterInputJob(null, 0, ListNum.PROTO);
		FastqWriterOutputJob outputProto=new FastqWriterOutputJob(0, null, ListNum.PROTO);
		oqs=new OrderedQueueSystem2<FastqWriterInputJob, FastqWriterOutputJob>(
			threads, true, inputProto, outputProto);
		
		// Open output stream
		//ffout.append() must be honored: app=t previously truncated (hardcoded false; replicated via stream.sh 2026-09-05)
		outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, false);
//		System.err.println("os class: "+outstream.getClass());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Starts formatting workers and the output thread; call once before submission. */
	@Override
	public void start(){spawnThreads();}
	
	/** Returns emitted-read count accumulated during normal output completion. */
	@Override
	public long readsWritten(){return readsWritten;}
	
	/** Returns emitted bases after normal completion; HEADER and ATTACHMENT count zero. */
	@Override
	public long basesWritten(){return basesWritten;}
	
	/**
	 * Wraps and submits a read list without copying it; may block for queue space.
	 * @param list Nonnull list, possibly empty or containing null entries
	 * @param id Unique batch ID in a dense sequence starting at zero
	 */
	@Override
	public final void add(final ArrayList<Read> list, final long id){addReads(new ListNum<Read>(list, id));}
	
	/**
	 * Submits a batch without copying its contents; may block for queue space.
	 * Producers must finish submission before poison(); retain empty batches so
	 * the dense ID sequence has no gaps. Do not mutate submitted reads or lists.
	 * @param reads Nonnull batch with a nonnull list and unique ID starting at zero
	 */
	@Override
	public void addReads(final ListNum<Read> reads){
		FastqWriterInputJob job=new FastqWriterInputJob(reads, reads.id(), ListNum.NORMAL);
		oqs.addInput(job);
	}
	
	/**
	 * Converts each SAM line into an unlinked Read and submits the original batch ID.
	 * Sequence and quality arrays are shared; SAM mate numbers are retained for selection.
	 * Alignment fields, mate links and attachment objects are not retained.
	 * @param lines Nonnull batch of nonnull SAM lines; IDs follow addReads' contract
	 */
	@Override
	public void addLines(final ListNum<SamLine> lines){
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines){
			//STR-367: Preserve SAM mate identity before formatters select R1 or R2.
			final Read r=new Read(sl.seq, sl.qual, sl.qname, -1, false);
			r.setPairnum(sl.pairnum());
			reads.add(r);
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}
	
	/** Signals the end of submission; may block while queuing terminal markers. */
	@Override
	public void poison(){oqs.poison();}
	
	/**
	 * Waits for the queue's completion signal after poison or an error.
	 * Normal completion includes output finalization and counter accumulation. Forced
	 * completion releases the wait without guaranteeing that every thread has exited.
	 * @return True if this writer recorded an error
	 */
	@Override
	public boolean waitForFinish(){
		oqs.waitForFinish();
		return errorState;
	}
	
	/** Signals end of input and waits for completion; returns true on recorded error. */
	@Override
	public boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/** Force-finishes after external failure; queued output may be abandoned. */
	@Override
	public synchronized void finishError(){
		errorState=true;
		oqs.setFinished(true);
	}
	
	/** Returns the recorded error flag; use waitForFinish for normal completion. */
	@Override
	public boolean errorState(){return errorState;}
	
	/** Returns whether queues are finished and no writer error is currently recorded. */
	@Override
	public boolean finishedSuccessfully(){return !errorState && oqs.finished();}
	
	/** Returns the output filename from the constructor's descriptor. */
	@Override
	public final String fname(){return fname;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn worker and writer threads. */
	void spawnThreads(){
		final int totalThreads=threads+1;// Workers plus writer
		
		alpt=new ArrayList<ProcessThread>(totalThreads);
		for(int i=0; i<totalThreads; i++){
			alpt.add(new ProcessThread(i));
		}
		
		for(ProcessThread pt : alpt){
			pt.start();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Read batch or terminal marker for the ordered input queue. */
	private static class FastqWriterInputJob implements HasID{//TODO: This class should just be a ln
		/** Creates a normal job with reads, or a terminal/prototype job with null reads. */
		FastqWriterInputJob(final ListNum<Read> reads_, final long id_, final int type_){
			reads=reads_;
			id=id_;
			type=type_;
		}
		
		@Override public long id(){return id;}
		@Override public boolean poison(){return type==ListNum.POISON;}
		@Override public boolean last(){return type==ListNum.LAST;}
		@Override public FastqWriterInputJob makePoison(final long id){
			return new FastqWriterInputJob(null, id, ListNum.POISON);
		}
		@Override public FastqWriterInputJob makeLast(final long id){
			return new FastqWriterInputJob(null, id, ListNum.LAST);
		}
		
		final ListNum<Read> reads;
		final long id;
		final int type;
	}
	
	/** Formatted byte block or terminal marker for the ordered output queue. */
	private static class FastqWriterOutputJob implements HasID{
		/** Creates a data job with bytes, or a terminal/prototype job with null bytes. */
		FastqWriterOutputJob(final long id_, final byte[] bytes_, final int type_){
			id=id_;
			bytes=bytes_;
			type=type_;
		}
		
		@Override public long id(){return id;}
		@Override public boolean poison(){return type==ListNum.POISON;}
		@Override public boolean last(){return type==ListNum.LAST;}
		@Override public FastqWriterOutputJob makePoison(final long id){
			return new FastqWriterOutputJob(id, null, ListNum.POISON);
		}
		@Override public FastqWriterOutputJob makeLast(final long id){
			return new FastqWriterOutputJob(id, null, ListNum.LAST);
		}
		
		final long id;
		final byte[] bytes;
		final int type;
	}
	
	/** Formats reads in the selected format, or writes the ordered byte blocks. */
	private class ProcessThread extends Thread{
		
		/** Assigns worker identity; ID zero is the sole output thread. */
		ProcessThread(final int tid_){
			tid=tid_;
			setName("FastqWriter-"+(tid==0 ? "Output" : "Worker-"+tid));
		}
		
		/** Runs one role; a failure records an error and releases queues before rethrowing. */
		@Override
		public void run(){
			try{
				synchronized(this){
					if(tid==0){
						writeOutput();// Writer thread
					}else{
						processJobs();// Worker thread
					}
					success=true;
				}
			}catch(Throwable t){
				//[stream/FastqWriter#001 FIXED 2026-06-21]: a thread death without poisoning/finishing the OQS2 HANGS -- a writer death
				//(write error in writeOutput) leaves outq undrained so workers block in addOutput + `finished` unset so main hangs in
				//waitForFinish; a worker death (a throw in processJobs/format) leaves its ordered job undelivered so the writer's getOutput
				//blocks forever on the gap. Escalate to a LOUD finish: errorState + oqs.setFinished(true) [force=true poisons inq+outq,
				//waking writer+workers+main], then rethrow so the thread dies loud (stack trace) instead of hanging. errorState surfaces via
				//poisonAndWait's return so the caller decides the exit. Mirrors the greenlit BamWriter#001.
				errorState=true;
				oqs.setFinished(true);
				throw new RuntimeException("FastqWriter thread "+tid+" failed", t);
			}
		}
		
		/** Writer thread - outputs ordered blocks to disk. */
		void writeOutput(){
			// Write ordered data blocks
			FastqWriterOutputJob job=oqs.getOutput();
			while(job!=null && !job.last()){
				try{
					outstream.write(job.bytes);
				}catch(Exception e){
					throw new RuntimeException("Error writing output", e);
				}
				
				job=oqs.getOutput();
			}
			
			// Wait for other threads and accumulate statistics
			ThreadWaiter.waitForThreadsToFinish(alpt);
			synchronized(FastqWriter.this){
				for(ProcessThread pt : alpt){
					if(pt!=this){// This thread not successful yet!
						synchronized(pt){
							readsWritten+=pt.readsWrittenT;
							basesWritten+=pt.basesWrittenT;
							errorState|=!pt.success;
						}
					}
				}
			}
			
			// Close output stream and signal completion
			boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			errorState|=b;
			oqs.setFinished(true);
		}
		
		/** Worker thread - converts reads to formatted bytes. */
		void processJobs(){
			final ByteBuilder bb=new ByteBuilder();
			
			FastqWriterInputJob job=oqs.getInput();
			while(job!=null && !job.poison()){
				assert(job.reads!=null) : job.last()+", "+job.poison();
				ArrayList<Read> reads=job.reads.list;
				
				// Format reads
				if(format==FASTQ){
					writeFastq(reads, bb);
				}else if(format==FASTA){
					writeFasta(reads, bb);
				}else if(format==HEADER){
					writeHeader(reads, bb);
				}else if(format==SCARF){
					writeScarf(reads, bb);
				}else if(format==ATTACHMENT){
					writeAttachment(reads, bb);
				}else{
					throw new RuntimeException("Bad format: "+format);
				}
				
				// Create output job
				FastqWriterOutputJob outJob=new FastqWriterOutputJob(job.id(), bb.toBytes(), ListNum.NORMAL);
				oqs.addOutput(outJob);
				bb.clear();
				
				job=oqs.getInput();
			}
			
			// Re-inject poison for other workers
			if(job!=null){oqs.addInput(job);}
		}
		
		/** Appends enabled mates as FASTQ, skipping null entries and counting reads/bases. */
		private void writeFastq(final ArrayList<Read> reads, final ByteBuilder bb){
			for(Read r : reads){
				if(r==null){continue;}
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					r1.toFastq(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r1.length();
				}
				if(writeR2 && r2!=null){
					r2.toFastq(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r2.length();
				}
			}
		}
		
		/** Appends enabled mates as FASTA, skipping null entries and counting reads/bases. */
		private void writeFasta(final ArrayList<Read> reads, final ByteBuilder bb){
			for(Read r : reads){
				if(r==null){continue;}
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					r1.toFasta(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r1.length();
				}
				if(writeR2 && r2!=null){
					r2.toFasta(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r2.length();
				}
			}
		}
		
		/** Appends enabled mates' IDs and newlines; increments reads but not bases. */
		private void writeHeader(final ArrayList<Read> reads, final ByteBuilder bb){
			for(Read r : reads){
				if(r==null){continue;}
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					bb.appendln(r1.id);
					readsWrittenT++;
				}
				if(writeR2 && r2!=null){
					bb.appendln(r2.id);
					readsWrittenT++;
				}
			}
		}

		/**
		 * Appends enabled mates' attachment text, with one newline per emitted read.
		 * Null entries are skipped; selected reads must have nonnull attachment objects.
		 * Increments the read counter but not the base counter.
		 * @param reads Batch supplied by processJobs
		 * @param bb Reusable output buffer
		 */
		private void writeAttachment(final ArrayList<Read> reads, final ByteBuilder bb){
			for(Read r : reads){
				if(r==null){continue;}
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					bb.append(r1.obj.toString()).nl();
					readsWrittenT++;
				}
				if(writeR2 && r2!=null){
					bb.append(r2.obj.toString()).nl();
					readsWrittenT++;
				}
			}
		}
		
		/**
		 * Appends enabled mates as SCARF, skipping null entries and updating emitted counters.
		 * R2 is either the entry itself (pairnum 1) or the linked mate of an R1 entry.
		 * @param reads Batch supplied by processJobs
		 * @param bb Reusable output buffer, with one newline appended per emitted read
		 */
		private void writeScarf(final ArrayList<Read> reads, final ByteBuilder bb){
			assert(reads!=null && bb!=null) : "processJobs must supply a read list and its reusable formatting buffer";
			for(Read r : reads){
				if(r==null){continue;}
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					r1.toScarf(bb).nl();
					readsWrittenT++;
					basesWrittenT+=r1.length();
				}
				if(writeR2 && r2!=null){
					r2.toScarf(bb).nl();
					readsWrittenT++;
					basesWrittenT+=r2.length();
				}
			}
		}
		
		/** Number of selected reads formatted by this worker. */
		protected long readsWrittenT=0;
		/** Bases formatted as FASTQ, FASTA or SCARF; zero for headers and attachments. */
		protected long basesWrittenT=0;
		/** True only if this thread completed successfully. */
		boolean success=false;
		/** Thread ID. */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file path */
	public final String fname;
	/** Output file format */
	final FileFormat ffout;
	/** Output file format as an int */
	public final int format;
	/** Output stream */
	OutputStream outstream;
	/** OQS for coordinating workers and writer */
	final OrderedQueueSystem2<FastqWriterInputJob, FastqWriterOutputJob> oqs;
	/** Number of worker threads */
	final int threads;
	/** Write R1 reads (pairnum==0) */
	final boolean writeR1;
	/** Write R2 reads (pairnum==1 or mate) */
	final boolean writeR2;
	/** Thread list for accumulation */
	private ArrayList<ProcessThread> alpt;
	/** Emitted reads, accumulated by the output thread during normal completion. */
	protected long readsWritten=0;
	/** Emitted sequence bases; HEADER and ATTACHMENT do not contribute. */
	protected long basesWritten=0;
	/** True if an error was encountered */
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final int FASTQ=FileFormat.FASTQ;
	private static final int FASTA=FileFormat.FASTA;
	private static final int HEADER=FileFormat.HEADER;
	private static final int SCARF=FileFormat.SCARF;
	private static final int ATTACHMENT=FileFormat.ATTACHMENT;
	private static final int UNKNOWN=FileFormat.UNKNOWN;
	
	/** Default number of formatting workers, before clamping to Shared.threads(). */
	public static int DEFAULT_THREADS=3;
	/** Compile-time diagnostic flag. */
	public static final boolean verbose=false;
	
}
