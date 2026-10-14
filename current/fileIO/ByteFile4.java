package fileIO;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.HasID;
import stream.OrderedQueueSystem;
import structures.IntList;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Byte-line reader with one input producer and parallel line-extraction workers.
 * Construction opens the input and starts the pipeline. The producer copies chunks
 * ending at the last LF; workers use Vector helpers to locate delimiters and copy
 * individual lines. OrderedQueueSystem returns chunks in their original ID order.
 * LF and its preceding CR are omitted; an unterminated final line is retained,
 * with a trailing CR also omitted. Blank lines share an empty byte array.
 *
 * Batch consumption is the fast path: nextList and nextListAsync return the same
 * queue payloads without ByteFile's record-origin numbering or outer monitor.
 * Scalar consumption has a separate, unsynchronized cache and requires one caller.
 * Switching from scalar reads or pushBack to batches is unsupported: queued batches
 * bypass any unread scalar cache (FIO-016). Scalar reads also consume internal chunk
 * IDs, so later batch IDs need not begin at zero. Pushback is a legacy scalar-only
 * facility, not a supported way to prepare a batch-reader header boundary.
 *
 * Use one consumption mode for a reader's lifetime. Returned nonempty line arrays
 * and batch lists are not recycled by this reader. EOF does not close the input;
 * close reports the latched read/close error flag and requests queue termination,
 * but does not join the threads. Reset and reopening a closed pipeline are unsupported.
 *
 * @author Brian Bushnell
 * @contributor Isla, Shinobu
 * @date November 2, 2025
 */
public final class ByteFile4 extends ByteFile{

	/** Prints a selected scalar line range or measures batch consumption.
	 * Print mode normalizes line endings to LF; speed mode consumes all batches.
	 * Reports read/write/close errors after the diagnostic returns.
	 * @param args Input, range, worker, buffer and common I/O options */
	public static void main(String[] args){
		String fname=null;
		long first=0, last=-1;
		boolean speedtest=false;
		int threads=DEFAULT_THREADS;
		int bufferSize=bufferlen;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("input")){
				fname=b;
			}else if(a.equals("speedtest")){
				speedtest=Parse.parseBoolean(b);
				if(speedtest){
					first=0;
					last=Long.MAX_VALUE;
				}
			}else if(a.equals("first") || a.equals("start")){
				first=Long.parseLong(b);
			}else if(a.equals("last") || a.equals("stop")){
				last=Long.parseLong(b);
			}else if(a.equals("lines")){
				first=0;
				last=Long.parseLong(b);
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("buffersize") || a.equals("buffer")){
				bufferSize=(int)Parse.parseKMG(b);
			}else if(a.equals("verbose")){
//				verbose=Parse.parseBoolean(b);
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//Do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//Do nothing
			}else if(arg.indexOf('=')<0 && fname==null){
				// First non-flag argument is the filename
				fname=arg;
			}else{
				System.err.println("Unknown parameter "+arg);
				assert(false) : "Unknown parameter "+arg;
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}

		if(fname==null){fname="stdin";}

		bufferlen=bufferSize;//FIO-014: construction starts the producer, so configure its buffer first.
		ByteFile4 bf=new ByteFile4(fname, true, threads);
		boolean error=speedtest(bf, first, last, !speedtest);
		error|=bf.close();
		if(error){throw new RuntimeException("ByteFile4: error reading, writing or closing "+fname);}
	}

	/** Runs the diagnostic and reports output errors; the caller checks reader closure separately.
	 * @param bf Already-started reader consumed by this diagnostic
	 * @param first Zero-based first line to print; ignored in batch speed mode
	 * @param last Exclusive print limit; a negative value means no limit
	 * @param reprint Whether to write scalar lines to borrowed stdout
	 * @return Whether printing detected an output error */
	private static boolean speedtest(ByteFile4 bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		boolean error=false;
		if(last<0){last=Long.MAX_VALUE;}
		bf.start();

		if(reprint){
			BufferedOutputStream bos=new BufferedOutputStream(System.out, 65536);
			for(long i=0; i<first; i++){bf.nextLine();}
			try{
				for(long i=first; i<last; i++){
					byte[] s=bf.nextLine();
					if(s==null){break;}

					lines++;
					bytes+=s.length+1;
					//				System.out.println(new String(s));
					bos.write(s);
					bos.write('\n');
				}
				bos.flush();//FIO-012: emit the final buffer without closing borrowed stdout.
				error|=System.out.checkError();//PrintStream records failures instead of throwing IOException.
			}catch(IOException e){
				error=true;
				e.printStackTrace();
			}
		}else{
			for(ListNum<byte[]> ln=bf.nextListAsync(); ln!=null; ln=bf.nextListAsync()){
				for(byte[] line : ln.list){
					lines++;
					bytes+=line.length+1;
				}
			}
		}
		t.stop();

		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
			System.err.println("Bytes: "+bytes);
		}
		return error;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens a named input with the default worker request.
	 * @param fname Input name or supported stream alias
	 * @param allowSubprocess_ Whether the descriptor may use subprocess decoding */
	public ByteFile4(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}

	/** Opens the supplied descriptor with the default worker request.
	 * @param ff Input descriptor retained by ByteFile */
	public ByteFile4(FileFormat ff){
		this(ff, DEFAULT_THREADS);
	}

	/** Resolves an input descriptor, opens it and starts the pipeline.
	 * @param fname Input name or supported stream alias
	 * @param allowSubprocess_ Whether subprocess decoding is allowed
	 * @param threads_ Requested extraction workers; values below one use DEFAULT_THREADS */
	public ByteFile4(String fname, boolean allowSubprocess_, int threads_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false), threads_);
	}

	/** Creates ordered queues and starts one producer plus the effective workers.
	 * @param ff Input descriptor retained by ByteFile
	 * @param threads_ Requested workers, defaulted below one and clamped to 1..Shared.threads() */
	public ByteFile4(FileFormat ff, int threads_){
		super(ff);
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());

		// Create OQS with prototypes
		BufferJob inputPrototype=new BufferJob(null, 0, BufferJob.PROTO);
		ListNum<byte[]> outputPrototype=new ListNum<byte[]>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<BufferJob, ListNum<byte[]>>(
			threads, true, inputPrototype, outputPrototype);

		if(verbose){System.err.println("Made ByteFile4 with "+threads+" threads");}
		start();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens input and starts threads once; already-open calls do nothing.
	 * Calling this after close cannot re-arm the final queue system and is unsupported. */
	public synchronized void start(){
		if(open){return;}
		open();
		spawnThreads();
	}

	/** Rejects rewind by reporting an error, closing and resetting only base counters. */
	@Override
	public final void reset(){
		//ByteFile4's producer/worker pipeline and its final, non-re-armable OrderedQueueSystem (setFinished is
		//permanent) cannot be restarted, so reset() cannot restore a readable state. Rather than silently return
		//an empty stream (contract violation), fail loudly: warn and latch errorState so close() reports a nonzero
		//exit. Callers needing reset() should use ByteFile1/ByteFile2 instead.
		System.err.println("Error: ByteFile4 does not support reset() for '"+name()+
			"'; the multithreaded read pipeline cannot be restarted. Use ByteFile1 or ByteFile2 if reset is required.");
		errorState=true;
		close();
		superReset();
	}

	/** Closes the underlying input and requests forced output-queue termination.
	 * Resets the delivered-line counter but retains the scalar cache; callers stop
	 * consuming before closure. Does not join producer/workers or make the reader reusable.
	 * @return Sticky read/close error flag, including a rejected reset request */
	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open);}
		if(!open){return errorState;}
		open=false;

		if(is!=null){
			errorState|=ReadWrite.finishReading(is, name(), ff.subprocess);
			is=null;
		}
		
		oqs.setFinished(true);

		lineNum.set(-1);
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name());}
		return errorState;
	}

	/** Returns one cached line, fetching the next ordered chunk when needed.
	 * Single-caller scalar mode only; do not mix with batch consumption.
	 * @return Line without its terminator, or null at the output terminal marker */
	@Override
	public final byte[] nextLine(){
		// Get next line from current list, or fetch new list
		if(currentList==null || listPos>=currentList.size()){
			currentList=oqs.getOutput();
			listPos=0;

			if(currentList==null || currentList.last()){
				if(currentList!=null && currentList.last()){
					oqs.setFinished(true);
				}
				return null;
			}
		}

		lineNum.incrementAndGet();
		return currentList.get(listPos++);
	}

	/*
	 * Performance note (2026-08-21): Deleting this override would let ByteFile.nextList()
	 * assign ListNum.firstRecordNum, but would also add its nlSync monitor to this hot path.
	 * ByteFile4's normal speedtest calls nextListAsync(), so two temporary A/B builds changed
	 * only that loop to call nextList().  On Dori ln004, using a 37,336,632,684-byte raw
	 * FASTQ in /tmp, t=2, and five ABBA-balanced runs per build, the direct path had an
	 * 8.795-second median (8.683-8.964) and the inherited synchronized path had an
	 * 8.986-second median (8.826-10.345).  Node-wide load was not recorded, so the exact
	 * slowdown is not conclusive.  The test did not establish that synchronization was free,
	 * therefore this override remains.  Callers needing exact record origins should assign
	 * them at an existing ordered boundary, as TaxTreeText.BatchSource does.
	 */
	/** Returns nextListAsync directly, without assigning firstRecordNum.
	 * @return Ordered batch, or null at the output terminal marker */
	@Override
	public final ListNum<byte[]> nextList(){
		return nextListAsync();//Avoid synchronization in ByteFile's default implementation
	}

	/** Takes an ordered output chunk directly from the queue.
	 * IDs correspond to producer chunks and are dense from zero in batch-only use.
	 * No scalar-cache reconciliation or record-origin numbering is performed.
	 * @return Owned batch/list payload, or null for LAST or queue termination */
	@Override
	public final ListNum<byte[]> nextListAsync(){
		//TODO: Probable bug (FIO-016) - nextLine can leave unread entries in currentList.
		//Fetching directly from oqs bypasses that scalar-read tail, including a pushed-back header boundary.
		//Deferred under IO-037: deprecated pushback/mixed-mode support must not slow the batch hot path.
		ListNum<byte[]> list=oqs.getOutput();
		if(list==null || list.last()){
			if(list!=null && list.last()){
				oqs.setFinished(true);
			}
			return null;
		}
		lineNum.addAndGet(list.size());
		return list;
	}

	/** Legacy scalar-cache rewind/insertion; not supported as a scalar-to-batch bridge.
	 * With a consumed cache entry, rewinds that entry and ignores the supplied array.
	 * Otherwise inserts the supplied array; does not decrement the delivered-line counter.
	 * @param line Array inserted only when there is no consumed entry to rewind */
	@Override
	public void pushBack(byte[] line){
		if(currentList==null){
			// Shouldn't happen, but be defensive
			currentList=new ListNum<byte[]>(new ArrayList<byte[]>(), -1);
			currentList.add(line);
			listPos=0;
		}else if(listPos>0){
			listPos--;
		}else{
			// listPos==0, insert at beginning
			assert(listPos==0 && currentList!=null);
			currentList.list.add(0, line);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts one producer and the configured number of extraction workers. */
	private void spawnThreads(){
		final int totalThreads=threads+1; // Workers plus producer

		alpt=new ArrayList<ProcessThread>(totalThreads);
		for(int i=0; i<totalThreads; i++){
			alpt.add(new ProcessThread(i));
		}

		for(ProcessThread pt : alpt){
			pt.start();
		}
	}

	/** Marks the reader open and acquires its input stream.
	 * @return Stream owned by this reader until close delegates to finishReading */
	private synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened ByteFile4 "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		return is;
	}

	/** Computes the next buffer capacity without allocating storage.
	 * FIO-015: long arithmetic prevents overflow; the cap follows ByteFile1.
	 * @param current Current array length
	 * @return Twice the current length, capped at Integer.MAX_VALUE-8
	 * @throws IOException If the full buffer is already at the cap, using the producer's read-error path */
	static int nextBufferSize(final int current) throws IOException{
		final int maxCap=Integer.MAX_VALUE-8;
		if(current>=maxCap){
			throw new IOException("ByteFile4: line exceeds reader capacity ("+maxCap+" bytes); current="+current);
		}
		return (int)Math.min(2L*current, maxCap);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/** @return Stored open flag, not a worker-liveness or EOF test */
	@Override
	public boolean isOpen(){return open;}

	/** @return Borrowed input stream; null after close */
	@Override
	public final InputStream is(){return is;}

	/** @return Delivered lines minus one; pushback replay counts again, close resets to -1 */
	@Override
	public final long lineNum(){return lineNum.longValue();}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Producer at tid zero, otherwise an extraction worker sharing the ordered queues. */
	private class ProcessThread extends Thread{

		/** @param tid_ Zero for the producer; positive values identify extraction workers */
		ProcessThread(final int tid_){
			tid=tid_;
			setName("ByteFile4-"+(tid==0 ? "Producer" : "Worker-"+tid));
		}

		/** Executes the assigned role; success records only a normal method return. */
		@Override
		public void run(){
			if(tid==0){
				produceBufferJobs();
				if(verbose){System.err.println("Producer finished.");}
			}else{
				processBufferJobs();
				if(verbose){System.err.println("Worker "+tid+" finished.");}
			}
			success=true;
		}

		/** Reads and copies complete-line chunks, retaining a partial suffix for the next read.
		 * EOF submits a final unterminated suffix. IOException latches an error; bytes
		 * not yet submitted are not delivered. Publishes poison, then joins other threads. */
		void produceBufferJobs(){
			long jobId=0;
			byte[] buffer=new byte[bufferlen];
			int bstop=0;

			try{
				while(true){
					if(verbose){System.err.println("Producer attempting read.");}
					// Read data into buffer
					final int r;
					synchronized(ByteFile4.this){
						if(is==null){break;}
						r=is.read(buffer, bstop, buffer.length-bstop);
					}
					if(r>0){
						bstop+=r;
					}

					// Handle EOF
					if(r<0){
						if(bstop>0){
							// Package remaining data
							byte[] copy=Arrays.copyOf(buffer, bstop);
							BufferJob job=new BufferJob(copy, jobId++);
							if(verbose){System.err.println("Producer adding EOF job "+job.id+".");}
							oqs.addInput(job);
						}
						break;
					}

					// Find last newline
					int lastNL=Vector.findLastSymbol(buffer, bstop, slashn);

					// If we have a newline, package job
					if(lastNL>=0){
						byte[] copy=Arrays.copyOf(buffer, lastNL+1);
						BufferJob job=new BufferJob(copy, jobId++);
						if(verbose){System.err.println("Producer adding normal job "+job.id+".");}
						oqs.addInput(job);

						// Shift remaining bytes to start
						int remaining=bstop-lastNL-1;
						if(remaining>0){
							System.arraycopy(buffer, lastNL+1, buffer, 0, remaining);
							bstop=remaining;
						}else{
							bstop=0;
						}
					}else{
						// No newline found - expand buffer if full
						if(bstop==buffer.length){
							buffer=Arrays.copyOf(buffer, nextBufferSize(buffer.length));
						}
					}
				}
			}catch(IOException e){
				if(!errorState){//Warn loudly once; errorState doubles as the spam-guard
					System.err.println("Error: truncated or corrupt input for '"+name()+"'; data may be incomplete.");
					if(verbose){e.printStackTrace();}
				}
				errorState=true;//Already submitted chunks may be delivered; close reports the incomplete input.
			}

			// Signal completion
			if(verbose){System.err.println("Producer adding poison.");}
			oqs.poison();

			// Wait for workers
			if(verbose){System.err.println("Producer waiting for finish.");}
			ThreadWaiter.waitForThreadsToFinish(alpt);
		}

		/** Copies lines from each input chunk and publishes one list with its original ID.
		 * Consumes LF/CRLF and a final unterminated suffix; reinserts input poison for peers. */
		void processBufferJobs(){
			if(verbose){System.err.println("Worker "+tid+" getInput()");}
			BufferJob job=oqs.getInput();
			if(verbose){System.err.println("Worker "+tid+" got job "+(job==null ? -1 : job.id));}
			IntList positions=new IntList();

			while(job!=null && !job.poison()){
				positions.clear();

				// Find all newlines in this buffer chunk
				Vector.findSymbols(job.buffer, 0, job.buffer.length, slashn, positions);

				ListNum<byte[]> output=new ListNum<byte[]>(
					new ArrayList<byte[]>(positions.size()), job.id());

				int start=0;
				for(int i=0; i<positions.size(); i++){
					int nlpos=positions.get(i);
					int limit=(nlpos>0 && job.buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;

					if(start==limit){
						output.add(blankLine);
					}else{
						byte[] line=Arrays.copyOfRange(job.buffer, start, limit);
						output.add(line);
					}
					start=nlpos+1;
				}

				// BufferJob should end exactly at a newline, no residual data
//				assert(start==job.buffer.length) : "Incomplete line in BufferJob! start="+start+
//					", length="+job.buffer.length;
				
				//Fix for files not terminated by newlines
				if(start<job.buffer.length){
					int limit=job.buffer.length;
					if(limit>start && job.buffer[limit-1]==slashr){limit--;}

					if(start==limit){output.add(blankLine);}else{
						byte[] line=Arrays.copyOfRange(job.buffer, start, limit);
						output.add(line);
					}
				}

				if(verbose){System.err.println("Worker "+tid+" addOutput("+output.id+")");}
				oqs.addOutput(output);
				if(verbose){System.err.println("Worker "+tid+" getInput()");}
				job=oqs.getInput();
				if(verbose){System.err.println("Worker "+tid+" got job "+(job==null ? -1 : job.id));}
			}
			if(verbose){System.err.println("Worker "+tid+" exited on "+
				(job==null ? "null" : "poison="+job.poison()));}

			// Re-inject poison for other workers
			if(job!=null){oqs.addInput(job);}
			if(verbose){System.err.println("Worker "+tid+" reinjected poison.");}
		}

		/** Normal-return marker; the reader does not aggregate this field. */
		boolean success=false;
		/** Role selector: zero is the producer. */
		final int tid;
	}

	/** Immutable chunk descriptor; normal jobs own a copied byte buffer, markers use null. */
	private static final class BufferJob implements HasID{

		/** @param buffer_ Copied chunk handed to one extraction worker
		 * @param id_ Zero-based producer chunk ID */
		BufferJob(byte[] buffer_, long id_){
			this(buffer_, id_, NORMAL);
		}

		/** @param buffer_ Chunk payload, or null for a prototype/terminal marker
		 * @param id_ Queue ordering ID
		 * @param type_ NORMAL, PROTO, LAST or POISON */
		BufferJob(byte[] buffer_, long id_, int type_){
			buffer=buffer_;
			id=id_;
			type=type_;
		}

		/** @return Ordering ID retained at construction */
		@Override
		public long id(){return id;}

		/** @return Whether this is an input poison marker */
		@Override
		public boolean poison(){return type==POISON;}

		/** @return Whether this is a LAST marker */
		@Override
		public boolean last(){return type==LAST;}

		/** @param id_ Terminal ordering ID
		 * @return New payload-free poison marker */
		@Override
		public BufferJob makePoison(long id_){
			return new BufferJob(null, id_, POISON);
		}

		/** @param id_ Terminal ordering ID
		 * @return New payload-free last marker */
		@Override
		public BufferJob makeLast(long id_){
			return new BufferJob(null, id_, LAST);
		}

		/** Copied chunk payload; null on markers. */
		final byte[] buffer;
		/** Stable producer/terminal ordering ID. */
		final long id;
		/** Marker or normal-payload type. */
		final int type;

		/** Factory prototype type. */
		static final int PROTO=-1;
		/** Ordinary chunk type. */
		static final int NORMAL=0;
		/** Ordered last-result marker type. */
		static final int LAST=3;
		/** Input termination marker type. */
		static final int POISON=4;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Single-producer input FIFO and ordered worker output; never re-armed. */
	private final OrderedQueueSystem<BufferJob, ListNum<byte[]>> oqs;
	/** Effective extraction-worker count, excluding the producer. */
	private final int threads;
	/** Producer followed by extraction workers; used by the producer to join peers. */
	private ArrayList<ProcessThread> alpt;

	/** Input lifecycle flag; not a worker-completion flag. */
	private boolean open=false;
	/** Input stream accessed under this monitor by producer/open/close. */
	private InputStream is;
	/** Successful scalar returns plus emitted batch sizes, minus one. */
	private AtomicLong lineNum=new AtomicLong(-1);

	/** Scalar-only cache, separate from direct batch consumption. */
	private ListNum<byte[]> currentList=null;
	/** Next scalar-cache position. */
	private int listPos=0;

	/** Shared zero-length representation of an empty line. */
	private static final byte[] blankLine=new byte[0];

	/** Compile-time diagnostic logging switch. */
	public static final boolean verbose=false;
	/** Input buffering request passed to ReadWrite when opening. */
	public static boolean BUFFERED=false;
	/** Initial producer buffer size, read when the producer starts. */
	public static int bufferlen=65536;
	/** Default extraction-worker request before clamping to Shared.threads(). */
	public static int DEFAULT_THREADS=2;

	/** Sticky producer read, input close or rejected-reset error. */
	private boolean errorState=false;
}
