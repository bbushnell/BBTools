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
 * @author Brian Bushnell
 * @date Oct 21, 2014
 *
 */
public class ByteStreamWriter extends Thread {
	
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
	
	/** @See primary constructor */
	public ByteStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_){
		this(fname_, overwrite_, append_, allowSubprocess_, 0);
	}
	
	/** @See primary constructor */
	public ByteStreamWriter(String fname_, boolean overwrite_, boolean append_, boolean allowSubprocess_, int format){
		this(FileFormat.testOutput(fname_, FileFormat.TEXT, format, 0, allowSubprocess_, overwrite_, append_, false));
	}
	
	/**
	 * Create a ByteStreamWriter for this FileFormat.
	 * @param ff Contains information about the file name, output format, etc.
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
	
	/** Creates and starts a ByteStreamWriter unless fname is null */
	public static ByteStreamWriter makeBSW(String fname, boolean ow, boolean append, boolean allowSub){
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, 0, 0, allowSub, ow, append, false);
		return makeBSW(ff);
	}

	/** Creates and starts a ByteStreamWriter unless ff is null */
	public static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Primary Method        ----------------*/
	/*--------------------------------------------------------------*/

	
	@Override
	public void run() {
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
		ReadWrite.finishWriting(null, outstream, fname, allowSubprocess);
		if(verbose){System.err.println("finish writing");}
		synchronized(this){notifyAll();}
		if(verbose){System.err.println("done");}
	}
	
	/**
	 * Core job processing loop that dequeues ByteBuilder jobs and writes them
	 * to the output stream. Continues until POISON2 sentinel is received.
	 * Includes retry logic for IOException handling with program termination
	 * on persistent write failures to prevent corrupt output files.
	 */
	public void processJobs() {
		
		ByteBuilder job=null;
		while(job==null){
			try {
				job=queue.take();
//				job.list=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(verbose){System.err.println("processing jobs");}
		while(job!=null && job!=POISON2){
			if(job.length()>0){
				for(int i=0; job!=null && i<100; i++) {
					try {
						outstream.write(job.array, 0, job.length());
						job=null;
					} catch (IOException e) {
						//Safest option is to exit here to avoid corrupt files.
						if(true) {
							System.err.println("Program is exiting due to a failed write.");
							KillSwitch.exceptionKill(e);
						}
						e.printStackTrace();
						
						System.err.println("Retrying...");
					}
				}
			}
			
			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Control and Helpers     ----------------*/
	/*--------------------------------------------------------------*/
	
	
	@Override
	public synchronized void start(){
		super.start();
		if(verbose){System.err.println(this.getState());}
		synchronized(this){
			while(!started){
				try {
					this.wait(20);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}

	
	/**
	 * Signals the writer thread to shut down gracefully.
	 * Waits for thread startup if necessary, flushes remaining buffer content,
	 * and sends POISON2 sentinel to terminate the processing loop.
	 */
	public synchronized void poison(){
		//Don't allow thread to shut down before it has started
		while(!started || this.getState()==Thread.State.NEW){
			try {
				this.wait(20);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(!open){return;}
		
		if(ordered){
			addOrdered(POISON2, maxJobID+1);
		}else{
			if(buffer!=null){addJob(buffer);}
		}
		buffer=null;
//		System.err.println("Poisoned!");
//		assert(false);
		
//		assert(false) : open+", "+this.getState()+", "+started;
		open=false;
		addJob(POISON2);
	}
	
	/** 
	 * Wait for this object's thread to terminate.
	 * Should be poisoned first.
	 */
	public void waitForFinish(){
		while(this.getState()!=Thread.State.TERMINATED){
			try {
				this.join(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Poison the thread, and wait for it to terminate.
	 * @return true if there was an error, false otherwise
	 */
	public boolean poisonAndWait(){
		poison();
		waitForFinish();
		return errorState;
	}
	
	//TODO Why is this synchronized?
	/**
	 * Adds a ByteBuilder job to the processing queue.
	 * Blocks until queue space is available and ensures job is successfully queued.
	 * Requires that start() has completed before jobs can be submitted.
	 * @param bb ByteBuilder containing data to write
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
			try {
				queue.put(bb);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				assert(!queue.contains(bb)); //Hopefully it was not added.
			}
		}
	}
	
	/** Forces immediate buffer flush regardless of current buffer size */
	public final void forceFlushBuffer(){
		flushBuffer(true);
	}
	
	/** Called after every write to the buffer */
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
	 * Implements flow control to prevent excessive memory usage in ordered mode.
	 *
	 * @param job ByteBuilder containing data to write
	 * @param jobID Sequence number for ordered processing
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
					try {
						this.wait(2000);
					} catch (InterruptedException e) {
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
	 * Processes jobs in sequential order using a HashMap to buffer out-of-order jobs.
	 * Immediately processes jobs when they arrive in the expected sequence,
	 * otherwise stores them until their turn arrives.
	 *
	 * @param job ByteBuilder containing data to write
	 * @param jobID Expected sequence number for this job
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
	 * Immediately queues the job without sequence tracking.
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
	 * Skip the  buffers and print directly.
	 * Mainly for headers with ordered streams.
	 * @param s String to print.
	 */
	public void forcePrint(String s){
		forcePrint(s.getBytes());
	}
	
	/** 
	 * Skip the  buffers and print directly.
	 * Mainly for headers with ordered streams.
	 * @param b Data to print.
	 */
	public void forcePrint(byte[] b){
		try {
			outstream.write(b, 0, b.length);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

	/**
	 * Returns the current buffer for direct manipulation.
	 * Requires that the writer is open and buffer exists.
	 * @return Current ByteBuilder buffer
	 */
	public ByteBuilder getBuffer() {
		assert(open);
		assert(buffer!=null);
		return buffer;
	}
	
	/** Avoid using this if possible. */
	public ByteStreamWriter print(CharSequence x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/**
	 * Appends StringBuilder to buffer and conditionally flushes.
	 * @param x StringBuilder to append
	 * @return This ByteStreamWriter for method chaining
	 */
	@Deprecated
	/** Avoid using this if possible. */
	public ByteStreamWriter print(StringBuilder x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/** Avoid using this if possible. */
	public ByteStreamWriter print(String x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		flushBuffer(false);
		return this;
	}
	
	/** Avoid using this if possible. */
	public ByteStreamWriter printt(String x){
		if(verbose){System.err.println("Added line '"+x+"'");}
		assert(open) : x;
		buffer.append(x);
		buffer.append('\t');
		flushBuffer(false);
		return this;
	}

	/** Appends tab character to buffer */
	public ByteStreamWriter tab(){return print('\t');}
	/** Appends newline character to buffer */
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
	 * Appends byte value to buffer and conditionally flushes.
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
	 * Appends character to buffer and conditionally flushes.
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
	 * @param len Number of bytes to append
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
	 * Appends character array to buffer and conditionally flushes.
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
	 * Appends ByteBuilder contents to current buffer and conditionally flushes.
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
	 * Supports FASTQ, FASTA, SAM, SITES, and INFO formats.
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
	 * Formats and appends a k-mer with count in standardized format.
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
	 * Formats and appends a long k-mer with count for k-mers longer than 32 bases.
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
	
	
	/** Appends newline character to buffer */
	public ByteStreamWriter println(){return print('\n');}
	/**
	 * Appends CharSequence followed by newline to buffer.
	 * @param x CharSequence to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(CharSequence x){print(x); return print('\n');}
	/**
	 * Appends String followed by newline to buffer.
	 * @param x String to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(String x){print(x); return print('\n');}
	/**
	 * Appends StringBuilder followed by newline to buffer.
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
	 * @param x Character to append
	 * @return This ByteStreamWriter for method chaining
	 */
	public ByteStreamWriter println(char x){print(x); return print('\n');}
//	public ByteStreamWriter println(byte[] x){print(x); return print('\n');}
	/**
	 * Appends character array followed by newline to buffer.
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
	 * Appends ByteBuilder followed by newline with optional direct submission.
	 * @param x ByteBuilder to append
	 * @param destroy Whether to append newline directly to x and submit as job
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
	
	/** Current buffer for accumulating output data */
	private ByteBuilder buffer;
	
	/** Initial buffer capacity when creating new buffers */
	public int initialLen=36000;
	/** Maximum buffer size before automatic flushing */
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
	
	/** For ordered output */
	private final HashMap<Long, ByteBuilder> map;
	/** Next expected job ID for ordered processing */
	private long nextJobID=0;
	/** Highest job ID seen for ordered processing */
	private long maxJobID=-1;
	
	/** Whether the writer is open and accepting data */
	private boolean open=true;
	/** Whether the worker thread has started */
	private volatile boolean started=false;
	
	/** TODO */
	public boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	/** Whether output format is BAM */
	private final boolean BAM;
	/** Whether output format is SAM */
	private final boolean SAM;
	/** Whether output format is FASTQ */
	private final boolean FASTQ;
	/** Whether output format is FASTA */
	private final boolean FASTA;
	/** Whether output format is BREAD (BBTools binary read format) */
	private final boolean BREAD;
	/** Whether output format is SITES */
	private final boolean SITES;
	/** Whether output format is INFO */
	private final boolean INFO;
	/** Whether output format is not a recognized bioinformatics format */
	private final boolean OTHER;
	
	/** Line wrap length for FASTA format output */
	private final int FASTA_WRAP=Shared.FASTA_WRAP;
	
	/*--------------------------------------------------------------*/

//	private static final ByteBuilder POISON=new ByteBuilder("POISON_ByteStreamWriter");
	/** Sentinel object to signal thread termination */
	private static final ByteBuilder POISON2=new ByteBuilder(1);
	
	/** Whether to print verbose debugging information */
	public static boolean verbose=false;
	/** Number of lists held before the stream blocks */
	private final int MAX_CAPACITY=256;
	/** Job count threshold for triggering flow control in ordered mode */
	private final int ADD_LIMIT=MAX_CAPACITY/2;
	/** Job count threshold for waking waiting threads in ordered mode */
	private final int HALF_LIMIT=ADD_LIMIT/4;
	
}
