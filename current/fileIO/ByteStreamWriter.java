package fileIO;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;

import assemble.Contig;
import dna.AminoAcid;
import dna.Data;
import kmer.AbstractKmerTable;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import ukmer.AbstractKmerTableU;



/**
 * Writes byte-builder jobs through a background thread to an output opened by
 * the constructor. Read formatting is optional; already formatted bytes pass
 * through unchanged. Compression and BAM conversion are delegated to ReadWrite.
 *
 * <p>Call {@link #start()} before submitting jobs or buffered output. In unordered
 * mode, one producer may use the print/println methods and the mutable buffer.
 * Those methods are not synchronized. Alternatively, producers may submit owned
 * builders with {@link #add(ByteBuilder, long)}; do not mix job submission with
 * buffered printing. Ordered mode has no print buffer and requires exactly one
 * job for every ID from zero through the final ID, including empty jobs.
 *
 * <p>Job submission transfers the builder and its backing array without copying;
 * do not mutate or reuse them afterward. Stop and join all producers before
 * {@link #poisonAndWait()}, then inspect its error result. Direct forcePrint
 * calls bypass the queue and need exclusive access to the output, normally for
 * headers before startup. Caught direct-write errors are retained for shutdown.
 * Missing ordered IDs fail loudly after shutdown is requested; poisonAndWait
 * waits for writer termination before propagating that failure.
 *
 * <p>Buffered string/character print and println overloads retain only the low
 * eight bits of each UTF-16 character; they do not encode Unicode text. ASCII
 * is preserved, but characters above U+00FF lose information. For a chosen text
 * encoding, encode explicitly and submit bytes. In contrast, forcePrint(String)
 * uses the platform default charset, so these paths are not interchangeable.
 *
 * @author Brian Bushnell, Shinobu
 * @date Oct 21, 2014
 *
 */
public class ByteStreamWriter extends Thread{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Performance test program for ByteStreamWriter.
	 * Creates a test byte array and writes it repeatedly to measure throughput.
	 * @param args Command-line arguments: [filename] [iterations]
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		final int alen=1000;
		byte[] array=new byte[alen];
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.numberToBase[i&3];
		}
		array[array.length-1]='\n';
		long iters=Long.parseLong(args[1]);
		String fname=args[0];
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, true);
		bsw.start();
		for(long i=0; i<iters; i++){
			bsw.print(array);
		}
		bsw.poisonAndWait();
		t.stop();
		System.err.println("MB/s: \t"+Tools.format("%.2f", ((alen*iters)/(t.elapsed/1000.0))));
		System.err.println("Time: \t"+t);
	}
	
	/** Opens unordered output with filename-based format detection and a TEXT fallback.
	 * Does not start the worker; see {@link #ByteStreamWriter(FileFormat)}.
	 * @param fname_ Output path
	 * @param overwrite_ Permission to replace an existing output file
	 * @param append_ Whether to append to existing raw or gzip output
	 * @param allowSubprocess_ Requested permission to use output subprocesses
	 */
	public ByteStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_){
		this(fname_, overwrite_, append_, allowSubprocess_, 0);
	}
	
	/** Opens unordered output with an optional format override (zero means detect).
	 * Does not start the worker; see {@link #ByteStreamWriter(FileFormat)}.
	 * @param fname_ Output path
	 * @param overwrite_ Permission to replace an existing output file
	 * @param append_ Whether to append to existing raw or gzip output
	 * @param allowSubprocess_ Requested permission to use output subprocesses
	 * @param format Optional FileFormat type override, or zero to autodetect
	 */
	public ByteStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_, int format){
		this(FileFormat.testOutput(fname_, FileFormat.TEXT, format, 0, allowSubprocess_, overwrite_, append_, false));
	}
	
	/**
	 * Opens output immediately and creates a five-job queue, but does not start
	 * the worker. Ordered mode uses a map of pending IDs; unordered mode allocates
	 * a print buffer. Append is accepted only for raw or gzip FileFormats here.
	 * @param ff Non-null output description, including ordering and overwrite policy
	 */
	public ByteStreamWriter(FileFormat ff){
		FASTQ=ff.fastq() || ff.text();
		FASTA=ff.fasta();
		BREAD=ff.bread();
		SAM=ff.samOrBam();
		BAM=ff.bam();
		SITES=ff.sites();
		INFO=ff.attachment();
		OTHER=(!FASTQ && !FASTA && !BREAD && !SAM && !BAM && !SITES && !INFO);
		
		
		fname=ff.name();
		overwrite=ff.overwrite();
		append=ff.append();
		allowSubprocess=ff.allowSubprocess();
		ordered=ff.ordered();
		assert(!(overwrite&append));
		assert(ff.canWrite()) : "File "+fname+" exists "+(new File(ff.name()).canWrite() ? 
				("and overwrite="+overwrite+".\nPlease add the flag ow to overwrite the file.\n") : 
					"and is read-only.");
		if(append && !(ff.raw() || ff.gzip())){throw new RuntimeException("Can't append to compressed files.");}
		
		if(!BAM || !Data.BAM_SUPPORT_OUT()){
			outstream=ReadWrite.getOutputStream(fname, append, true, allowSubprocess);
		}else{
			outstream=ReadWrite.getBamOutputStream(fname, append);
		}
		
		queue=new ArrayBlockingQueue<ByteBuilder>(5);
		if(ordered){
			buffer=null;
			map=new HashMap<Long, ByteBuilder>(MAX_CAPACITY);
		}else{
			buffer=new ByteBuilder(initialLen);
			map=null;
		}
	}
	
	/**
	 * Creates and starts a ByteStreamWriter unless fname is null.
	 * @param fname Output path
	 * @param ow Permission to replace an existing output file
	 * @param append Whether to append to existing raw or gzip output
	 * @param allowSub Requested permission to use output subprocesses
	 * @return Started writer, or null for a null path
	 */
	public static ByteStreamWriter makeBSW(String fname, boolean ow, boolean append, boolean allowSub){
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, 0, 0, allowSub, ow, append, false);
		return makeBSW(ff);
	}

	/**
	 * Creates and starts a ByteStreamWriter unless ff is null.
	 * @param ff Output description, or null
	 * @return Started writer, or null if ff is null
	 */
	public static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Primary Method        ----------------*/
	/*--------------------------------------------------------------*/

	
	/** Announces startup, consumes jobs through the terminal, then finishes output.
	 * Invoked by the thread machinery; callers use {@link #start()}. */
	@Override
	public void run(){
		if(verbose){System.err.println("running");}
		assert(open) : fname;
		
		synchronized(this){
			started=true;
			this.notify();
		}

		if(verbose){System.err.println("waiting for jobs");}
		
		processJobs();
		
		if(verbose){System.err.println("null/poison job");}
//		assert(false);
		open=false;
		errorState=ReadWrite.finishWriting(null, outstream, fname, allowSubprocess)|errorState;
		if(verbose){System.err.println("finish writing");}
		synchronized(this){notifyAll();}
		if(verbose){System.err.println("done");}
	}
	
	/**
	 * Core job processing loop that dequeues ByteBuilder jobs and writes them
	 * to the output stream. Continues until POISON2 sentinel is received.
	 * On a queued-write IOException, invokes KillSwitch to terminate the program;
	 * there is no write retry. Interrupted queue waits print the exception and retry.
	 * Only the writer thread may consume this queue; this method does not close output.
	 */
	public void processJobs(){
		
		ByteBuilder job=null;
		while(job==null){
			try{
				job=queue.take();
//				job.list=queue.take();
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(verbose){System.err.println("processing jobs");}
		while(job!=null && job!=POISON2){
			if(job.length()>0){
				try{
					outstream.write(job.array, 0, job.length());
				}catch(IOException e){
					//Safest option is to exit here to avoid a corrupt file; there is no retry.
					System.err.println("Program is exiting due to a failed write.");
					KillSwitch.exceptionKill(e);
				}
			}
			
			job=null;
			while(job==null){
				try{
					job=queue.take();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Control and Helpers     ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Starts the worker and waits for its startup handshake.
	 * A writer closed by closeBeforeStart is terminal and cannot be started.
	 * Buffered writes must not precede this call. Interrupted waits are logged
	 * and retried. */
	@Override
	public synchronized void start(){
		if(closedBeforeStart){
			throw new IllegalStateException("ByteStreamWriter cannot start after closeBeforeStart(): "+fname);
		}
		//ADDED [ByteStreamWriter start-guard]: catch the write-before-start bug class. Buffered print/println before start()
		//force an addJob that asserts(started) (-> -ea crash once output>maxLen), or are lost if the writer is never started.
		//forcePrint() writes straight to outstream and never touches buffer, so pre-start header writes stay exempt.
		//Motivated by tax/ImgRecord#001 (never started/closed) + tax/AnalyzeAccession#001 (started after the writes).
		assert(buffer==null || buffer.isEmpty()) : "ByteStreamWriter.start() called after "+buffer.length()+" bytes were already buffered via print/println; start() must precede all buffered writes. Use forcePrint() for pre-start direct writes (e.g. headers).";
		super.start();
		if(verbose){System.err.println(this.getState());}
		synchronized(this){
			while(!started){
				try{
					this.wait(20);
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Makes one terminal attempt to close constructor-opened output before the
	 * worker starts. Repeated calls perform no cleanup and return the current error state.
	 * This destructive close discards buffered print bytes, queued jobs and ordered
	 * map entries. Callers must stop all producers and all direct forcePrint users
	 * before closing; forcePrint bytes already handed to the backend are preserved
	 * on successful backend cleanup. Unchecked cleanup failures are rethrown after
	 * recording both the error and terminal state. Calling this after start() is invalid.
	 * @return True if direct writes or backend cleanup reported an error
	 * @throws IllegalStateException If the worker has already been started
	 */
	public synchronized boolean closeBeforeStart(){
		if(closedBeforeStart){return errorState;}
		if(started || getState()!=Thread.State.NEW){
			throw new IllegalStateException("closeBeforeStart() requires a NEW, unstarted writer: state="+
					getState()+", started="+started+", fname="+fname);
		}
		buffer=null;
		queue.clear();
		if(map!=null){map.clear();}
		open=false;
		try{
			errorState=ReadWrite.finishWriting(null, outstream, fname, allowSubprocess)|errorState;
		}catch(RuntimeException|Error e){
			errorState=true;
			throw e;
		}finally{
			closedBeforeStart=true;
			notifyAll();
		}
		return errorState;
	}

	
	/**
	 * Queues the remaining buffered bytes and a terminal after all producers finish.
	 * Requires a separately started worker unless closeBeforeStart already closed
	 * output; this method waits for startup and returns immediately for a preclosed
	 * writer. Normal return requests shutdown but does not wait for closure.
	 * Ordered IDs must be complete: a detected gap marks errorState, discards the
	 * unwritable suffix and queues shutdown before throwing. The accepted contiguous
	 * prefix is still written; the incomplete output must not be treated as successful.
	 * @throws RuntimeException If ordered jobs remain behind a missing ID, after shutdown is queued
	 */
	public synchronized void poison(){
		//Don't allow thread to shut down before it has started
		while(!closedBeforeStart && (!started || this.getState()==Thread.State.NEW)){
			try{
				this.wait(20);
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(closedBeforeStart || !open){return;}
		
		OrderedOutputException failure=null;
		if(ordered){
			//FIO-004: diagnose the missing ID before inserting any terminal into the ordered map.
			//All producers have stopped, so only the already-queued contiguous prefix is writable.
			if(!map.isEmpty()){
				failure=new OrderedOutputException(nextJobID, map.size());
				errorState=true;
				map.clear();
			}
		}else{
			if(buffer!=null){addJob(buffer);}
		}
		buffer=null;
//		System.err.println("Poisoned!");
//		assert(false);
		
//		assert(false) : open+", "+this.getState()+", "+started;
		open=false;
		addJob(POISON2);
		if(failure!=null){throw failure;}
	}
	
	/** 
	 * Waits for thread termination after a successful shutdown request.
	 * Returns immediately for output closed by closeBeforeStart. Does not start or
	 * poison the worker. Interrupted joins are logged and retried; the interrupted
	 * status is not restored here.
	 */
	public void waitForFinish(){
		while(!closedBeforeStart && this.getState()!=Thread.State.TERMINATED){
			try{
				this.join(1000);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Requests shutdown and joins the worker after producers have finished.
	 * A missing ordered ID is rethrown after joining the worker, which finishes output.
	 * Other poison failures propagate immediately because shutdown may not be queued.
	 * Includes recorded direct-write
	 * failures when direct calls finish before shutdown; see {@link #forcePrint(byte[])}.
	 * @return The errorState accumulated by direct writes, ordered gaps, output closure or a caller
	 * @throws RuntimeException If ordered jobs remain behind a missing ID
	 */
	public boolean poisonAndWait(){
		OrderedOutputException failure=null;
		try{poison();}
		catch(OrderedOutputException e){failure=e;}
		waitForFinish();
		if(failure!=null){throw failure;}
		return errorState;
	}

	/** Identifies the gap failure thrown only after a terminal has been queued.
	 * This permits poisonAndWait to join without catching unrelated pre-terminal failures. */
	private static final class OrderedOutputException extends RuntimeException{
		/** Captures the missing ID and buffered data-job count before the suffix is discarded.
		 * @param missingID Next ordered job ID required for a contiguous prefix
		 * @param bufferedJobs Number of out-of-order buffered jobs that will be discarded
		 */
		OrderedOutputException(final long missingID, final int bufferedJobs){
			super("ByteStreamWriter: ordered output incomplete - jobID "+missingID+
				" was never supplied (a producer failed or emitted nothing); "+bufferedJobs+
				" buffered job(s) could not be written. Every jobID in [0,maxJobID] must be add()ed exactly once; "+
				"a producer that yields no bytes must still add an empty ByteBuilder for its id.");
		}
		/** Serialization identity for this diagnostic exception. */
		private static final long serialVersionUID=1L;
	}
	
	//Synchronized to share this monitor with add()/addOrdered()/poison() and to gate on the started handshake; queue.put() is itself thread-safe, so the lock is belt-and-suspenders.
	/**
	 * Adds a ByteBuilder job to the processing queue.
	 * Blocks until queue space is available and ensures job is successfully queued.
	 * Requires that start() has completed before jobs can be submitted.
	 * Bypasses ID ordering and the print buffer; ordered producers should use add.
	 * No bytes are copied. Do not mutate the job or its array after submission.
	 * Interrupted queue waits are logged and retried.
	 * @param bb Non-null builder whose ownership is transferred to the writer
	 */
	public synchronized void addJob(ByteBuilder bb){
//		System.err.println("Got job "+(j.list==null ? "null" : j.list.size()));
		
		assert(started) : "Wait for start() to return before using the writer.";
//		while(!started || this.getState()==Thread.State.NEW){
//			try {
//				this.wait(20);
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		
		boolean success=false;
		while(!success){
			try{
				queue.put(bb);
				success=true;
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
				assert(!queue.contains(bb)); //Hopefully it was not added.
			}
		}
	}
	
	/** Queues a nonempty unordered print buffer. Does not flush the underlying
	 * OutputStream or wait for the queued bytes to be written. */
	public final void forceFlushBuffer(){
		flushBuffer(true);
	}
	
	/** Queues the unordered print buffer when its byte threshold is met, or when
	 * forced and nonempty. Replaces submitted buffers; callers must stop using old
	 * references. The threshold is checked after an append and may be exceeded.
	 * @param force Whether to submit a nonempty buffer below the threshold */
	public final void flushBuffer(boolean force){
		final int x=buffer.length();
		if(x>=maxLen || (force && x>0)){
			addJob(buffer);
			buffer=new ByteBuilder(initialLen);
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------           Ordering           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a job with ordering support for maintaining output sequence.
	 * Uses ordered or unordered processing depending on configuration.
	 * Ordered flow control limits queued job counts, not bytes: a later ID arriving
	 * with at least 128 pending entries waits until fewer than 32 remain or its ID
	 * becomes next. The next required ID is allowed through to make progress.
	 *
	 * <p><b>Ordered-mode contract:</b> jobIDs must be supplied contiguously over 0..maxJobID,
	 * each exactly once. A producer that yields no bytes for its id must STILL add an empty
	 * ByteBuilder for that id; a skipped jobID stalls the stream, and {@link #poison()} then
	 * reports a gap by throwing after queuing shutdown and discarding the unwritable suffix.
	 * Stop all producers before requesting shutdown. Submitted builders and arrays
	 * must not be mutated or reused; their bytes are not copied.
	 *
	 * @param job Non-null builder transferred to this writer, even if empty
	 * @param jobID Dense nonnegative sequence number; ignored in unordered mode
	 */
	public synchronized void add(ByteBuilder job, long jobID){
		
		if(ordered){
			int size=map.size();
//			System.err.print(size+", ");
//			System.err.println("A.Adding job "+jobID+"; next="+nextJobID+"; max="+maxJobID+", map="+map.keySet());
			final boolean flag=(size>=HALF_LIMIT);
			if(jobID>nextJobID && size>=ADD_LIMIT){
//				if(printBufferNotification){
//					System.err.println("Output buffer became full; key "+jobID+" waiting on "+nextJobID+".");
//					printBufferNotification=false;
//				}
				while(jobID>nextJobID && size>=HALF_LIMIT){
					try{
						this.wait(2000);
					}catch(InterruptedException e){
						e.printStackTrace();
					}
					size=map.size();
				}
//				if(printBufferNotification){
//					System.err.println("Output buffer became clear for key "+jobID+"; next="+nextJobID+", size="+size);
//				}
			}
//			System.err.println("B.Adding ordered job "+jobID+"; next="+nextJobID+"; max="+maxJobID);
			addOrdered(job, jobID);
			assert(jobID!=nextJobID);
			if(flag && jobID<nextJobID){this.notifyAll();}
		}else{
			addDisordered(job);
		}
	}
	
	/**
	 * Stores a job by ID and forwards any contiguous prefix to the worker queue.
	 * Does not write the bytes itself or copy the builder. An empty map wakes
	 * waiting producers; the caller supplies the outer flow-control policy.
	 *
	 * @param job ByteBuilder containing data to write
	 * @param jobID Unique sequence number no smaller than the next expected ID
	 */
	private synchronized void addOrdered(ByteBuilder job, long jobID){
//		System.err.println("addOrdered "+jobID+"; nextJobID="+nextJobID);
//		assert(false);
		assert(ordered);
		assert(job!=null) : jobID;
		assert(jobID>=nextJobID) : jobID+", "+nextJobID;
		maxJobID=Tools.max(maxJobID, jobID);
		ByteBuilder old=map.put(jobID, job);
		assert(old==null);
//		System.err.println("C.Adding ordered job "+jobID+"; next="+nextJobID+"; max="+maxJobID+", map="+map.keySet());
		
		if(jobID==nextJobID){
			do{
				ByteBuilder value=map.remove(nextJobID);
				//			System.err.println("Removing and queueing "+nextJobID+": "+value.toString());
				addJob(value);
				nextJobID++;
				//			System.err.println("D.nextJobID="+nextJobID);
			}while(map.containsKey(nextJobID));
			
			if(map.isEmpty()){notifyAll();}
		}else{
			assert(!map.containsKey(nextJobID));
		}
	}
	
	/**
	 * Adds a job for unordered processing.
	 * Queues the job without sequence tracking. Requires an empty print buffer;
	 * buffered print calls must not be mixed with job submission.
	 * @param job ByteBuilder containing data to write
	 */
	private synchronized void addDisordered(ByteBuilder job){
		assert(!ordered);
		assert(buffer==null || buffer.isEmpty());
		addJob(job);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Print             ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Writes directly using the platform default charset; see the byte[] overload
	 * for exclusive-access and error reporting. Does not append a newline.
	 * @param s Non-null string to encode and write
	 */
	public void forcePrint(String s){
		forcePrint(s.getBytes());
	}
	
	/** 
	 * Writes directly, bypassing the queue and buffer, without flushing or adding a
	 * newline. Intended for headers before startup or while the worker is idle; do not
	 * overlap with queued writes, other direct calls or shutdown. Caught IOExceptions
	 * are printed and latched in errorState for poisonAndWait, rather than rethrown.
	 * Buffered output can defer a failure until flush/close, so inspect shutdown status.
	 * @param b Non-null bytes to write synchronously to the underlying stream
	 */
	public void forcePrint(byte[] b){
		try{
			outstream.write(b, 0, b.length);
		}catch(IOException e){
			//FIO-005: retain direct-write failure even if later flush/close succeeds.
			errorState=true;
			e.printStackTrace();
		}
	}
	

	/**
	 * Borrows the current unordered print buffer for the single producer.
	 * After direct appends, call flushBuffer to apply the threshold. Do not retain
	 * or mutate the reference after a flush submits it or shutdown begins.
	 * @return Current mutable buffer; unavailable in ordered mode
	 */
	public ByteBuilder getBuffer(){
		assert(open);
		assert(buffer!=null);
		return buffer;
	}
	
	/** Appends the low eight bits of each character, without charset encoding.
	 * Prefer explicitly encoded bytes for text that is not ASCII.
	 * @param x Characters to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(CharSequence x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends the low eight bits of each character and conditionally flushes.
	 * Does not encode Unicode text; prefer explicitly encoded bytes for non-ASCII text.
	 * @param x StringBuilder to append
	 * @return This ByteStreamWriter for method chaining
	 */
	@Deprecated
	public ByteStreamWriter print(StringBuilder x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/** Appends the low eight bits of each character, without charset encoding.
	 * Unlike forcePrint(String), this does not use the platform default charset.
	 * @param x String to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(String x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/** Appends the low eight bits of each character followed by a tab.
	 * Does not encode Unicode text; see {@link #print(String)}.
	 * @param x String to append before the tab
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printt(String x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		buffer.append('\t');
		flushBuffer(false);
		return this;
	}

	/**
	 * Appends tab character to buffer.
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter tab(){return print('\t');}
	/**
	 * Appends newline character to buffer.
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter nl(){return print('\n');}
	
	/**
	 * Appends boolean value to buffer and conditionally flushes.
	 * @param x Boolean value to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(boolean x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends integer value to buffer and conditionally flushes.
	 * @param x Integer value to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(int x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends integer followed by tab character to buffer.
	 * @param x Integer value to append before tab
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printt(int x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		buffer.append((byte)'\t');
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends long value to buffer and conditionally flushes.
	 * @param x Long value to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(long x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends long followed by tab character to buffer.
	 * @param x Long value to append before tab
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printt(long x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x);
		buffer.append((byte)'\t');
		flushBuffer(false);
		return this;
	}
	
//	public ByteStreamWriter print(float x){
//		if(verbose){System.err.println("Added line '"+(x)+"'");}
//		assert(open) : x;
//		buffer.appendSlow(x);
//		flushBuffer(false);
//		return this;
//	}
//	
//	public ByteStreamWriter print(double x){
//		if(verbose){System.err.println("Added line '"+(x)+"'");}
//		assert(open) : x;
//		buffer.appendSlow(x);
//		flushBuffer(false);
//		return this;
//	}
	
	/**
	 * Appends float with specified decimal places to buffer.
	 * @param x Float value to append
	 * @param decimals Number of decimal places
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(float x, int decimals){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x, decimals);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends float with decimal places followed by tab to buffer.
	 * @param x Float value to append
	 * @param decimals Number of decimal places
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printt(float x, int decimals){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x, decimals);
		buffer.append('\t');
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends double with specified decimal places to buffer.
	 * @param x Double value to append
	 * @param decimals Number of decimal places
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(double x, int decimals){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x, decimals);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends double with decimal places followed by tab to buffer.
	 * @param x Double value to append
	 * @param decimals Number of decimal places
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printt(double x, int decimals){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : x;
		buffer.append(x, decimals);
		buffer.append('\t');
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends one literal byte, not its decimal representation, and conditionally flushes.
	 * @param x Byte value to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(byte x){
		if(verbose){System.err.println("Added line '"+((char)x)+"'");}
		assert(open) : ((char)x);
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends the low eight bits of a character and conditionally flushes.
	 * @param x Character to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(char x){
		if(verbose){System.err.println("Added line '"+(x)+"'");}
		assert(open) : (x);
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends byte array to buffer and conditionally flushes.
	 * @param x Byte array to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(byte[] x){
		if(verbose){System.err.println("Added line '"+new String(x)+"'");}
		assert(open) : new String(x);
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends byte array followed by newline to buffer.
	 * @param x Byte array to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(byte[] x){
		if(verbose){System.err.println("Added line '"+new String(x)+"'");}
		assert(open) : new String(x);
		buffer.append(x).nl();
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends partial byte array to buffer and conditionally flushes.
	 * @param x Byte array to append
	 * @param len Maximum number of bytes from the start of x to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(byte[] x, int len){
		if(verbose){System.err.println("Added line '"+new String(x)+"'");}
		assert(open) : new String(x);
		buffer.append(x, len);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends characters cast to bytes and conditionally flushes; this is not a charset encoder.
	 * @param x Character array to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(char[] x){
		if(verbose){System.err.println("Added line '"+new String(x)+"'");}
		assert(open) : new String(x);
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Copies ByteBuilder contents to the current buffer and conditionally flushes.
	 * The supplied builder may be reused after return.
	 * @param x ByteBuilder to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(ByteBuilder x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends ByteBuilder to buffer with optional direct job submission.
	 * When destroy is true and current buffer is empty, submits the ByteBuilder
	 * directly as a job rather than copying to buffer.
	 * In that case, ownership of x and its backing array passes to the writer;
	 * neither may be mutated or reused afterward. This overload requires unordered mode.
	 *
	 * @param x ByteBuilder to append
	 * @param destroy Whether to submit directly when buffer is empty
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(ByteBuilder x, boolean destroy){
		if(!destroy || buffer.length()>0){print(x);}
		else{
			if(verbose){System.err.println("Added line '"+x+"'");}
			assert(open) : x;
			addJob(x);
		}
		return this;
	}
	
	/**
	 * Formats and appends a Read object according to the configured output format.
	 * TEXT uses FASTQ formatting; BAM uses SAM text for downstream conversion.
	 * BREAD uses native Read text; SITES and INFO delegate to their Read formatters.
	 * Adds no extra newline and does not print the mate. Requires unordered mode and exclusive
	 * access to r while formatting, including any temporary native-text match mutation.
	 * @param r Read object to format and append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(Read r){
		assert(!OTHER);
		ByteBuilder x=(FASTQ ? r.toFastq(buffer) : FASTA ? r.toFasta(FASTA_WRAP, buffer) : SAM ? r.toSam(buffer) :
			SITES ? r.toSites(buffer) : INFO ? r.toInfo(buffer) : r.toText(true, buffer));
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Formats and appends a Contig object as FASTA format.
	 * @param c Contig object to format and append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter print(Contig c){
		assert(!OTHER);
		c.toFasta(FASTA_WRAP, buffer);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Delegates k-mer/count formatting to AbstractKmerTable, then checks the buffer threshold.
	 *
	 * @param kmer K-mer encoded as long
	 * @param count K-mer occurrence count
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printKmer(long kmer, long count, int k){
		AbstractKmerTable.toBytes(kmer, count, k, buffer);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Formats and appends a k-mer with integer array values.
	 *
	 * @param kmer K-mer encoded as long
	 * @param values Associated integer values array
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printKmer(long kmer, int[] values, int k){
		AbstractKmerTable.toBytes(kmer, values, k, buffer);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Delegates multiword k-mer/count formatting to AbstractKmerTableU.
	 *
	 * @param array K-mer encoded as long array
	 * @param count K-mer occurrence count
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printKmer(long[] array, long count, int k){
		AbstractKmerTableU.toBytes(array, count, k, buffer);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Formats and appends a long k-mer with integer values array.
	 *
	 * @param array K-mer encoded as long array
	 * @param values Associated integer values array
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printKmer(long[] array, int[] values, int k){
		AbstractKmerTableU.toBytes(array, values, k, buffer);
		flushBuffer(false);
		return this;
	}
	
//	public ByteStreamWriter printKmer(long kmer, long[] values, int k){
//		kmer64.AbstractKmerTable64.toBytes(kmer, values, k, buffer);
//		flushBuffer(false);
//		return this;
//	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------           Println            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Appends newline character to buffer.
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(){return print('\n');}
	/**
	 * Appends CharSequence followed by newline to buffer.
	 * Characters are narrowed to bytes as in {@link #print(CharSequence)}.
	 * @param x CharSequence to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(CharSequence x){print(x); return print('\n');}
	/**
	 * Appends String followed by newline to buffer.
	 * Characters are narrowed to bytes as in {@link #print(String)}.
	 * @param x String to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(String x){print(x); return print('\n');}
	/**
	 * Appends StringBuilder followed by newline to buffer.
	 * Characters are narrowed to bytes as in {@link #print(StringBuilder)}.
	 * @param x StringBuilder to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(StringBuilder x){print(x); return print('\n');}
	/**
	 * Appends integer followed by newline to buffer.
	 * @param x Integer to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(int x){print(x); return print('\n');}
	/**
	 * Appends long followed by newline to buffer.
	 * @param x Long to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(long x){print(x); return print('\n');}
//	public void println(float x){print(x); print('\n');}
//	public void println(double x){print(x); print('\n');}
	/**
	 * Appends float with decimal places followed by newline to buffer.
	 * @param x Float value to append
	 * @param d Number of decimal places
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(float x, int d){print(x, d); return print('\n');}
	/**
	 * Appends double with decimal places followed by newline to buffer.
	 * @param x Double value to append
	 * @param d Number of decimal places
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(double x, int d){print(x, d); return print('\n');}
	/**
	 * Appends byte followed by newline to buffer.
	 * @param x Byte to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(byte x){print(x); return print('\n');}
	/**
	 * Appends character followed by newline to buffer.
	 * Retains only the low eight bits of the character, as in {@link #print(char)}.
	 * @param x Character to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(char x){print(x); return print('\n');}
//	public ByteStreamWriter println(byte[] x){print(x); return print('\n');}
	/**
	 * Appends character array followed by newline to buffer.
	 * Characters are narrowed to bytes as in {@link #print(char[])}.
	 * @param x Character array to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(char[] x){print(x); return print('\n');}
	/**
	 * Appends ByteBuilder followed by newline to buffer.
	 * @param x ByteBuilder to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(ByteBuilder x){print(x); return print('\n');}
	/**
	 * Copies a ByteBuilder followed by newline. When destroy is true, appends the
	 * newline to x itself before copying it; this overload never submits x as a job.
	 * @param x ByteBuilder to append
	 * @param destroy Whether x may be mutated by appending the newline
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(ByteBuilder x, boolean destroy){
		if(destroy){return print(x.append('\n'));}else{print(x); return print('\n');}
	}
	/**
	 * Formats and appends k-mer with count followed by newline.
	 *
	 * @param kmer K-mer encoded as long
	 * @param count K-mer occurrence count
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printlnKmer(long kmer, int count, int k){printKmer(kmer, count, k); return print('\n');}
	/**
	 * Formats and appends k-mer with values followed by newline.
	 *
	 * @param kmer K-mer encoded as long
	 * @param values Associated integer values
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printlnKmer(long kmer, int[] values, int k){printKmer(kmer, values, k); return print('\n');}
	/**
	 * Formats and appends long k-mer with count followed by newline.
	 *
	 * @param array K-mer encoded as long array
	 * @param count K-mer occurrence count
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printlnKmer(long[] array, int count, int k){printKmer(array, count, k); return print('\n');}
	/**
	 * Formats and appends long k-mer with values followed by newline.
	 *
	 * @param array K-mer encoded as long array
	 * @param values Associated integer values
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printlnKmer(long[] array, int[] values, int k){printKmer(array, values, k); return print('\n');}
	/**
	 * Formats and appends Read followed by newline.
	 * @param r Read object to format and append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(Read r){print(r); return print('\n');}
	/**
	 * Formats and appends Contig followed by newline.
	 * @param c Contig object to format and append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(Contig c){print(c); return print('\n');}

	/**
	 * Formats and appends k-mer with long count followed by newline.
	 *
	 * @param kmer K-mer encoded as long
	 * @param count K-mer occurrence count as long
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printlnKmer(long kmer, long count, int k){printKmer(kmer, count, k); return print('\n');}
//	public ByteStreamWriter printlnKmer(long kmer, long[] values, int k){printKmer(kmer, values, k); return print('\n');}
	/**
	 * Formats and appends long k-mer with long count followed by newline.
	 *
	 * @param array K-mer encoded as long array
	 * @param count K-mer occurrence count as long
	 * @param k K-mer length
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter printlnKmer(long[] array, long count, int k){printKmer(array, count, k); return print('\n');}
//	public ByteStreamWriter printlnKmer(long[] array, long[] values, int k){printKmer(array, values, k); return print('\n');}
	

	
	/**
	 * Formats and appends Read with optional mate pair.
	 * @param r Read object to format and append
	 * @param paired Whether to also print the mate read if present
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(Read r, boolean paired){
		println(r);
		if(paired && r.mate!=null){println(r.mate);}
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Inherited          ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return "BSW for "+fname;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Single producer's mutable print buffer; null in ordered mode and after shutdown is queued. */
	private ByteBuilder buffer;
	
	/** Capacity used for newly allocated print buffers, not a maximum size. */
	public int initialLen=36000;
	/** Soft byte threshold for queue submission, checked after each buffered append. */
	public int maxLen=32768;
	/** Whether to overwrite existing output files */
	public final boolean overwrite;
	/** Whether to append to existing output files */
	public final boolean append;
	/** Whether to allow subprocess for file compression */
	public final boolean allowSubprocess;
	/** Output file name */
	public final String fname;
	/** Whether output must maintain sequential ordering */
	public final boolean ordered;
	/** Output stream for writing data */
	private final OutputStream outstream;
	/** Queue for passing write jobs to the worker thread */
	private final ArrayBlockingQueue<ByteBuilder> queue;
	
	/** Ordered jobs waiting for missing earlier IDs; null in unordered mode. */
	private final HashMap<Long, ByteBuilder> map;
	/** Next expected job ID for ordered processing */
	private long nextJobID=0;
	/** Highest job ID seen for ordered processing */
	private long maxJobID=-1;
	
	/** Set false when shutdown is queued or run completes; not a thread-safe public status API. */
	private boolean open=true;
	/** Whether the worker thread has started */
	private volatile boolean started=false;
	/** Whether closeBeforeStart made its terminal cleanup attempt. */
	private volatile boolean closedBeforeStart=false;
	
	/** Accumulates direct-write, ordered-gap, finishWriting and caller-set errors. Read after joining.
	 * Direct calls must finish before shutdown; queued write failures invoke KillSwitch. */
	public boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	/** Whether output format is BAM */
	private final boolean BAM;
	/** Whether Read formatting uses SAM text, including BAM conversion output. */
	private final boolean SAM;
	/** Whether Read formatting uses FASTQ, including a generic TEXT FileFormat. */
	private final boolean FASTQ;
	/** Whether output format is FASTA */
	private final boolean FASTA;
	/** Whether output format is BREAD (native Read text). */
	private final boolean BREAD;
	/** Whether output format is SITES */
	private final boolean SITES;
	/** Whether output format is INFO */
	private final boolean INFO;
	/** Whether output format is not a recognized bioinformatics format */
	private final boolean OTHER;
	
	/** FASTA line wrap length captured from Shared when this writer is constructed. */
	private final int FASTA_WRAP=Shared.FASTA_WRAP;
	
	/*--------------------------------------------------------------*/

//	private static final ByteBuilder POISON=new ByteBuilder("POISON_ByteStreamWriter");
	/** Sentinel object to signal thread termination */
	private static final ByteBuilder POISON2=new ByteBuilder(1);
	
	/** Whether to print verbose debugging information */
	public static boolean verbose=false;
	/** Initial HashMap capacity, not a hard limit on jobs or buffered bytes. */
	private final int MAX_CAPACITY=256;
	/** Job count threshold for triggering flow control in ordered mode */
	private final int ADD_LIMIT=MAX_CAPACITY/2;
	/** Pending-count low-water threshold checked by waiting ordered producers. */
	private final int HALF_LIMIT=ADD_LIMIT/4;
	
}
