package fileIO;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import parse.Parse;
import parse.Parser;
import shared.KillSwitch;
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
 * Extracts byte lines using one input producer and configured extraction workers.
 * Construction opens input and starts all threads. The producer copies chunks and
 * newline offsets into jobs; workers copy individual lines into ordered output lists.
 * Returned line arrays remain available to consumers, while blank lines share one
 * empty array. CR before LF and a final unterminated CR are both stripped.
 *
 * Scalar consumption caches a chunk list. Batch consumption returns any unread
 * cached tail before taking another output chunk, retaining its producer chunk ID. Batch
 * sizes follow producer chunks, not the base TARGET_LIST_SIZE. Inherited nextList
 * serializes its own callers and numbers record origins; direct scalar/async calls
 * are not coordinated with that monitor. Async lists retain the default origin -1.
 *
 * Consuming terminal output signals queue completion but does not close input.
 * close closes the stream and clears the open flag, without stopping/joining the
 * pipeline or clearing scalar cache (FIO-022). Reset is unsupported and terminates
 * the JVM with a diagnostic even when assertions are disabled. Configure shared
 * settings before construction and keep these lifecycle
 * and mixed-mode limitations separate from ordinary full-drain reading.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @contributor Shinobu
 * @date November 2, 2025
 */
public final class ByteFile3 extends ByteFile{

	/** Prints a selected line range or reports throughput through the existing pipeline.
	 * Closes input after consumption and reports its read/close error flags.
	 * Print-mode stdout failures are reported after read/close errors. The buffer option
	 * is assigned before construction starts the producer.
	 * @param args Input path or in=, first/last/lines, speedtest, threads, buffer size
	 * and shared parser options; absent input defaults to stdin */
	public static void main(String[] args){
		String fname=null;
		long first=0, last=100;
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
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(Parser.parseCommonStatic(arg, a, b)){
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
		
		//FIO-021: publish the parsed setting before construction starts the producer.
		ByteFile3.bufferlen=bufferSize;
		ByteFile3 bf=new ByteFile3(fname, true, threads);
		bf.start();
		speedtest(bf, first, last, !speedtest);
		final boolean outputError=!speedtest && System.out.checkError();
		
		//FIO-020/023: report close errors without invoking unsupported reset.
		final boolean error=bf.close();
		if(error){throw new RuntimeException("ByteFile3: error reading or closing "+bf.name());}
		if(outputError){throw new RuntimeException("ByteFile3: error writing output");}
	}

	/** Skips scalar lines, then prints a range or drains batches through nextList.
	 * The speed path ignores last and counts one delimiter byte per returned line.
	 * @param bf Reader consumed by the diagnostic
	 * @param first Number of scalar lines to skip
	 * @param last Exclusive print limit
	 * @param reprint Whether to print decoded line Strings to borrowed stdout */
	private static void speedtest(ByteFile3 bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){bf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length+1;
				//FIO-024: main checks the borrowed PrintStream after consumption, before closing input.
				System.out.println(new String(s));
			}

			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
//			for(long i=first; i<last; i++){
//				byte[] s=bf.nextLine();
//				if(s==null){break;}
//				lines++;
//				bytes+=s.length+1;
//			}
			for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()){
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
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a TEXT-fallback descriptor and starts with the default worker request.
	 * @param fname Input filename or supported stream alias
	 * @param allowSubprocess_ Whether the descriptor allows subprocess decoding */
	public ByteFile3(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}

	/** Resolves a descriptor and starts the producer and requested extraction workers.
	 * @param fname Input filename or supported stream alias
	 * @param allowSubprocess_ Whether the descriptor allows subprocess decoding
	 * @param threads_ Worker request, bounded by the descriptor constructor */
	public ByteFile3(String fname, boolean allowSubprocess_, int threads_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false), threads_);
	}

	/** Starts reading the supplied descriptor with DEFAULT_THREADS requested workers.
	 * @param ff Borrowed input descriptor retained by the base class */
	public ByteFile3(FileFormat ff){
		this(ff, DEFAULT_THREADS);
	}
	
	/** Creates ordered job queues, opens input and starts one producer plus workers.
	 * Uses DEFAULT_THREADS for a nonpositive request, then bounds the worker count
	 * between one and Shared.threads(). Queue capacity is derived from that count.
	 * @param ff Borrowed input descriptor
	 * @param threads_ Requested extraction-worker count, excluding the producer */
	public ByteFile3(FileFormat ff, int threads_){
		super(ff);
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());

		// Create OQS with prototypes
		BufferJob inputPrototype=new BufferJob(null, null, 0, BufferJob.PROTO);
		ListNum<byte[]> outputPrototype=new ListNum<byte[]>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<BufferJob, ListNum<byte[]>>(
			threads, true, inputPrototype, outputPrototype);

		if(verbose){System.err.println("Made ByteFile3 with "+threads+" threads");}
		start();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens and starts threads when the open flag is false; otherwise does nothing.
	 * Reuses the existing queue object, so this is not a pipeline-reset operation. */
	public synchronized void start(){
		if(open){return;}
		open();
		spawnThreads();
	}

	/** Rejects unsupported rewind by terminating the JVM, regardless of assertion mode.
	 * The fixed queue pipeline cannot be restarted; use ByteFile1 or ByteFile2 for rereading. */
	@Override
	public final void reset(){
		KillSwitch.kill("ByteFile3 does not support reset. Use ByteFile1 or ByteFile2 for rereading: "+name());
	}

	/** For an open reader, closes the stream, clears open and resets the line counter.
	 * An already-closed call returns the current error flag without further cleanup.
	 * Does not coordinate producer access, clear cached output, or stop/join workers.
	 * @return Current sticky read/close flag, not a final joined-pipeline status */
	@Override
	public synchronized final boolean close(){
		//TODO: Probable bug (FIO-022) - closing input does not signal the queues or join workers;
		//an early consumer stop can leave output producers blocked and input still in use.
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open);}
		if(!open){return errorState;}
		open=false;

		if(is!=null){
			errorState|=ReadWrite.finishReading(is, name(), allowSubprocess());
			is=null;
		}

		lineNum=-1;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name());}
		return errorState;
	}

	/** Returns cached lines before taking another ordered chunk; counts delivered lines.
	 * Terminal output signals queue completion. Requires serialized scalar consumption.
	 * @return Copied line bytes, a shared empty array for blanks, or null at terminal output */
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

		lineNum++;
		return currentList.get(listPos++);
	}
	
	/** Returns a cached scalar-list remainder before taking another producer chunk.
	 * A remainder gets a new list retaining its chunk ID and the original line arrays.
	 * Counts delivered lines and signals queue completion on null/LAST output.
	 * @return Chunk-ID-bearing line list with default origin -1, or null at terminal output */
	@Override
	public final ListNum<byte[]> nextListAsync(){
		//FIO-026: preserve the unread scalar tail before advancing the output queue.
		if(currentList!=null && listPos<currentList.size()){
			final ArrayList<byte[]> tail=new ArrayList<byte[]>(currentList.list.subList(listPos, currentList.size()));
			final ListNum<byte[]> ln=new ListNum<byte[]>(tail, currentList.id);
			currentList=null;
			listPos=0;
			lineNum+=ln.size();
			return ln;
		}
		ListNum<byte[]> ln=oqs.getOutput();
		if(ln==null || ln.last()){
			oqs.setFinished(true);
			return null;
		}

		lineNum+=ln.size();
		return ln;
	}

	/** Rejects pushback unconditionally.
	 * @param line Unused line argument
	 * @throws UnsupportedOperationException always */
	@Override
	public void pushBack(byte[] line){
		throw new UnsupportedOperationException("pushBack not supported in ByteFile3");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces the thread list and starts producer index zero followed by workers. */
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

	/** Sets open and acquires a stream through ReadWrite; leaves queue/cache state intact.
	 * @return Newly acquired input stream
	 * @throws RuntimeException if the open flag is already set */
	private synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened ByteFile3 "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		return is;
	}

	/** Computes the next buffer capacity without allocating storage.
	 * FIO-025: long arithmetic prevents overflow; the cap follows ByteFile1 and ByteFile4.
	 * @param current Current array length
	 * @return Twice the current length, capped at Integer.MAX_VALUE-8
	 * @throws IOException If the full buffer is already at the cap, using the producer's read-error path */
	static int nextBufferSize(final int current) throws IOException{
		final int maxCap=Integer.MAX_VALUE-8;
		if(current>=maxCap){
			throw new IOException("ByteFile3: line exceeds reader capacity ("+maxCap+" bytes); current="+current);
		}
		return (int)Math.min(2L*current, maxCap);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/** @return Stored open flag; EOF consumption does not clear it */
	@Override
	public boolean isOpen(){return open;}

	/** @return Borrowed producer stream, or null after close; direct reads are uncoordinated */
	@Override
	public final InputStream is(){return is;}

	/** @return Stored delivered-line count minus one; closing open input resets it to -1 */
	@Override
	public final long lineNum(){return lineNum;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Index zero produces chunk jobs; positive indices extract copied line arrays. */
	private class ProcessThread extends Thread{

		/** Assigns this thread's role and diagnostic name.
		 * @param tid_ Zero for the producer, positive for an extraction worker */
		ProcessThread(final int tid_){
			tid=tid_;
			setName("ByteFile3-"+(tid==0 ? "Producer" : "Worker-"+tid));
		}

		/** Executes the assigned role, then marks success if the role returns normally. */
		@Override
		public void run(){
			if(tid==0){
				produceBufferJobs();
			}else{
				processBufferJobs();
			}
			success=true;
		}

		/** Reads input, scans LF positions and submits copied chunks with dense IDs.
		 * Preserves a newline-free suffix between fills, growing the buffer when full.
		 * EOF submits a remaining suffix; caught IOException drops unsubmitted residual
		 * bytes and latches errorState. Then publishes terminal markers and waits for
		 * workers; ThreadWaiter skips this calling producer rather than joining itself. */
		void produceBufferJobs(){
			long jobId=0;
			IntList positions=new IntList();
			byte[] buffer=new byte[bufferlen];
			int bstop=0;

			try{
				while(true){
					// Fill buffer
					int r=is.read(buffer, bstop, buffer.length-bstop);
					if(r>0){
						bstop+=r;
					}

					// Find newlines
					positions.clear();
					Vector.findSymbols(buffer, 0, bstop, slashn, positions);

					// If no newlines and buffer is full, expand
					if(positions.size()==0 && bstop==buffer.length){
						//FIO-025: a full capped buffer follows the existing IOException error path.
						buffer=Arrays.copyOf(buffer, nextBufferSize(buffer.length));
						continue;
					}

					// If we hit EOF with no newlines, package what we have
					if(r<0){
						if(bstop>0){
							byte[] copy=Arrays.copyOf(buffer, bstop);
							BufferJob job=new BufferJob(copy, positions.copy(), jobId++);
							oqs.addInput(job);
						}
						break;
					}

					// Package job with data up to last newline
					if(positions.size()>0){
						int lastNL=positions.get(positions.size()-1);
						byte[] copy=Arrays.copyOf(buffer, lastNL+1);
						BufferJob job=new BufferJob(copy, positions.copy(), jobId++);
						oqs.addInput(job);

						// Shift remaining bytes to start
						int remaining=bstop-lastNL-1;
						System.arraycopy(buffer, lastNL+1, buffer, 0, remaining);
						bstop=remaining;
					}
				}
			}catch(IOException e){
				if(!errorState){//Warn loudly once; errorState doubles as the spam-guard
					System.err.println("Error: truncated or corrupt input for '"+name()+"'; data may be incomplete.");
					if(verbose){e.printStackTrace();}
				}
				errorState=true;//Already submitted jobs can finish; unsubmitted residual bytes are not enqueued here.
			}

			// Signal completion
			oqs.poison();

			// Wait for workers
			ThreadWaiter.waitForThreadsToFinish(alpt);
		}

		/** Converts chunk jobs to copied lines, preserving chunk IDs and within-chunk order.
		 * Strips CR before LF and a terminal CR in an unterminated suffix. Reuses blankLine
		 * for empty records, then reinserts input poison for the remaining workers. */
		void processBufferJobs(){
			BufferJob job=oqs.getInput();
			while(job!=null && !job.poison()){
				ListNum<byte[]> output=new ListNum<byte[]>(
					new ArrayList<byte[]>(job.positions.size()+1), job.id());

				int start=0;
				for(int i=0; i<job.positions.size(); i++){
					int nlpos=job.positions.get(i);
					int limit=(nlpos>0 && job.buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;

					if(start==limit){
						output.add(blankLine);
					}else{
						byte[] line=Arrays.copyOfRange(job.buffer, start, limit);
						output.add(line);
					}
					start=nlpos+1;
				}

//				// Handle any remaining bytes after last newline
//				if(start<job.buffer.length){
//					byte[] line=Arrays.copyOfRange(job.buffer, start, job.buffer.length);
//					if(line.length>0){
//						output.add(line);
//					}
//				}

				//Handle any remaining bytes after last newline
				//Also handles files not terminated with newlines
				if(start<job.buffer.length){
					int limit=job.buffer.length;
					if(limit>start && job.buffer[limit-1]==slashr){limit--;}
					if(start==limit){output.add(blankLine);}
					else{
						byte[] line=Arrays.copyOfRange(job.buffer, start, limit);
						output.add(line);
					}
				}

				oqs.addOutput(output);
				job=oqs.getInput();
			}

			// Re-inject poison for other workers
			if(job!=null){oqs.addInput(job);}
		}

		/** Set on normal run completion; this reader does not aggregate the flag. */
		boolean success=false;
		/** Zero selects input production; positive indices select extraction. */
		final int tid;
	}

	/** Owns copied chunk bytes/newline offsets, or null payloads for control markers. */
	private static final class BufferJob implements HasID{

		/** Creates an ordinary chunk job.
		 * @param buffer_ Copied chunk bytes
		 * @param positions_ Copied LF offsets within the chunk
		 * @param id_ Dense producer-assigned job ID */
		BufferJob(byte[] buffer_, IntList positions_, long id_){
			this(buffer_, positions_, id_, NORMAL);
		}

		/** Creates a data job or payload-free marker/prototype.
		 * @param buffer_ Chunk bytes, or null for a marker/prototype
		 * @param positions_ LF offsets, or null for a marker/prototype
		 * @param id_ Job or marker ID
		 * @param type_ NORMAL, LAST, POISON or PROTO */
		BufferJob(byte[] buffer_, IntList positions_, long id_, int type_){
			buffer=buffer_;
			positions=positions_;
			id=id_;
			type=type_;
		}

		/** @return Stored job ID */
		@Override
		public long id(){return id;}

		/** @return Whether this is an input poison marker */
		@Override
		public boolean poison(){return type==POISON;}

		/** @return Whether this is a LAST marker */
		@Override
		public boolean last(){return type==LAST;}

		/** @param id_ Marker ID to retain
		 * @return New payload-free poison marker */
		@Override
		public BufferJob makePoison(long id_){
			return new BufferJob(null, null, id_, POISON);
		}

		/** @param id_ Marker ID to retain
		 * @return New payload-free LAST marker */
		@Override
		public BufferJob makeLast(long id_){
			return new BufferJob(null, null, id_, LAST);
		}

		/** Copied chunk bytes; null on markers/prototypes. */
		final byte[] buffer;
		/** Copied LF offsets relative to buffer; null on markers/prototypes. */
		final IntList positions;
		/** Stable producer/marker ID. */
		final long id;
		/** One of the job-type constants below. */
		final int type;

		/** Prototype used only to manufacture queue markers. */
		static final int PROTO=-1;
		/** Ordinary data-bearing chunk. */
		static final int NORMAL=0;
		/** End-of-output marker type, available through the HasID factory. */
		static final int LAST=3;
		/** Input termination marker forwarded between workers. */
		static final int POISON=4;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Final queue pipeline; constructed once and not restartable. */
	private final OrderedQueueSystem<BufferJob, ListNum<byte[]>> oqs;
	/** Extraction-worker count, excluding the producer. */
	private final int threads;
	/** Producer at index zero followed by workers; replaced by spawnThreads. */
	private ArrayList<ProcessThread> alpt;

	/** Stream-open flag, not an EOF or worker-completion flag. */
	private boolean open=false;
	/** Stream shared with the producer; close does not join before clearing it. */
	private InputStream is;
	/** Delivered lines minus one; scalar adds one, batch adds list size, close resets. */
	private long lineNum=-1;

	/** Scalar consumer's current output chunk; batching first returns any unread tail. */
	private ListNum<byte[]> currentList=null;
	/** Next scalar index in currentList; retained across close. */
	private int listPos=0;

	/** Shared array returned for empty lines. */
	private static final byte[] blankLine=new byte[0];

	/** Diagnostic logging switch. */
	public static boolean verbose=false;
	/** Whether ReadWrite should buffer the input stream in addition to producer chunks. */
	public static boolean BUFFERED=false;
	/** Initial producer buffer capacity, read when its thread begins. */
	public static int bufferlen=131072;
	/** Default requested extraction-worker count. */
	public static int DEFAULT_THREADS=1;

	/** Sticky producer read/stream-close error flag; no joined completion guarantee. */
	private boolean errorState=false;
}
