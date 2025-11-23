package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Single-threaded FASTQ writer with simple buffering.
 * Now supports optional internal threading and ordered/unordered modes via JobQueue,
 * allowing MT producers while maintaining low memory overhead.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 18, 2025
 */
public class FastqWriterST2 implements Writer{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	// Default constructor maintains lightweight, unthreaded, but can optionally order
	public FastqWriterST2(String out_, boolean writeR1_, boolean writeR2_, boolean overwrite){
		this(FileFormat.testOutput(out_, FileFormat.FASTQ, null, true, overwrite, false, true), 
			writeR1_, writeR2_, false, 3);
	}
	
	public FastqWriterST2(FileFormat ffout_, boolean writeR1_, boolean writeR2_){
		this(ffout_, writeR1_, writeR2_, false, 3);
	}

	/**
	 * Full constructor with configuration flags.
	 * * @param ffout_ Output file format/path
	 * @param writeR1_ Write read 1
	 * @param writeR2_ Write read 2
	 * @param ordered_ If true, output will be ordered by ID (requires JobQueue overhead)
	 * @param threaded_ If true, a separate thread will handle writing (Producer-Consumer model)
	 * @param queueCapacity_ Size of the buffer (JobQueue)
	 */
	public FastqWriterST2(FileFormat ffout_, boolean writeR1_, boolean writeR2_, 
			boolean threaded_, int queueCapacity_){
		ffout=ffout_;
		fname=ffout_.name();
		writeR1=writeR1_;
		writeR2=writeR2_;
		format=(ffout.format()==UNKNOWN ? FASTQ : ffout.format());
		assert(format==FASTQ || format==FASTA || format==HEADER) : ffout;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Config
		ordered=ffout.ordered();
		threaded=threaded_;
		
		// Only create queue if we need ordering, threading, or need to handle MT producers.
		if(ordered || threaded){
			// JobQueue handles ordering, capacity bounds, and backpressure.
			queue=new JobQueue<ListNum<Read>>(queueCapacity_, ordered, true, 0);
			queue.name="*FQWriterST2";
		} else{
			// Pure lightweight ST/unordered mode (needs external synchronization if MT)
			queue=null;
		}
		
		// Open output stream
		outstream=ReadWrite.getOutputStream(fname, false, true, false);
		if(verbose){outstream2.println("Made FastqWriterST (Ordered: "+ordered+", Threaded: "+threaded+")");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		synchronized(this) {
			if(started) {return;}
			started=true;
			if(threaded && queue != null){
				writerThread=new Thread(new WriterRunnable());
				writerThread.start();
			}
			if(verbose){outstream2.println("Started "+getClass().getName());}
		}
	}
	
	@Override
	public synchronized void poison(){
		if(verbose){outstream2.println("Called poison "+getClass().getName());}
		// Ensure only one thread executes the shutdown sequence
		if(poisoned){return;}
		poisoned=true;
		
		if(queue != null){
			// 1. Calculate the ID for the poison pill. It must be > maxSeen to guarantee order.
			long poisonID=queue.maxSeen() + 1;
			
			// 2. Create the poison pill (ListNum acts as the HasID object)
			ListNum<Read> poison=new ListNum<Read>(null, poisonID, ListNum.POISON);
			
			// 3. Inject the poison pill into the queue using the explicit JobQueue API.
			queue.poison(poison, false); 
			
			// 4. CRITICAL: If running in Host-Driven (unthreaded) mode, the producer thread 
			// calling poison() MUST perform the final drain.
			if(!threaded){
				// Use take() to block and ensure everything is drained until the pill is found.
				while(true){
					// take() will return null when it retrieves the poison pill.
					ListNum<Read> job=queue.take(); 
					
					if(job==null) break;
					
					writeReads(job.list);
				}
			}
			// If threaded, the worker thread handles the take() loop and exits gracefully.
		}
		if(verbose){outstream2.println("Finished poison "+getClass().getName());}
	}
	
	@Override
	public synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}
	
	@Override
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}
	
	@Override
	public void addReads(ListNum<Read> reads){
		if(verbose){outstream2.println("Called addReads "+(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));}
		if(reads==null){return;}
		if(!started){start();}
		
		if(queue==null){
			// LIGHTWEIGHT MODE (Unordered/Unthreaded): Direct synchronized write.
			writeReads(reads.list);
		} else{
			if(verbose){outstream2.println("Adding to queue "+(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));}
			// 1. Add to queue (blocks if full/backpressure enabled)
			queue.add(reads);
			
			if(!threaded){
				// HOST-DRIVEN MODE (Ordered/Unthreaded): 
				// The producer acts as the draining consumer immediately after adding.
				// Uses the non-blocking poll() to drain any contiguous jobs that are ready.

				assert(ordered);
				synchronized(queue) {//Essential to maintain ordering
					for(ListNum<Read> job=queue.poll(); job!=null && !job.poison(); job=queue.poll()){
						if(verbose){
							outstream2.println("Draining queue "+
								(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));
						}
						writeReads(job.list);
					}
				}
			}
			// If threaded, the background thread handles the draining.
		}
		if(verbose){outstream2.println("Finished addReads "+(reads==null ? "null" : reads.id+", "+reads.poison()+", "+reads.last()));}
	}
	
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines){
			reads.add(new Read(sl.seq, sl.qual, sl.qname, -1, false));
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}
	
	@Override
	public long readsWritten(){return readsWritten;}
	
	@Override
	public long basesWritten(){return basesWritten;}
	
	@Override
	public synchronized boolean waitForFinish(){

		if(verbose){
			outstream2.println("FastqWriterST close() 1");
			new Exception().printStackTrace();
		}
		if(finished) {return errorState;}
		poison(); 
		// poison() handles draining the queue if unthreaded.
		if(verbose){outstream2.println("FastqWriterST close() 2");}

		// If threaded, wait for the worker to process the poison pill and exit.
		if(threaded && writerThread != null){
			try{
				writerThread.join();
			} catch (InterruptedException e){
				Thread.currentThread().interrupt();
			}
		}
		if(verbose){outstream2.println("FastqWriterST close() 3");}

		boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
		if(verbose){outstream2.println("FastqWriterST close() 4");}
		finished=true;
		return errorState|=b;
	}
	
	@Override
	public boolean errorState(){return errorState;}
	
	@Override
	public boolean finishedSuccessfully(){return !errorState && finished;}
	
	@Override
	public final String fname() {return fname;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Logic          ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeReads(ArrayList<Read> reads){
		ByteBuilder bb=new ByteBuilder();
		
		// Format reads
		if(format==FASTQ){
			writeFastq(reads, bb);
		}else if(format==FASTA){
			writeFasta(reads, bb);
		}else if(format==HEADER){
			writeHeader(reads, bb);
		}else{
			throw new RuntimeException("Bad format: "+format);
		}

		// The actual write to stream must be synchronized
		write(bb);
		bb=null;
	}
	
	private void write(ByteBuilder bb){
		if(bb.length()<0){return;}
		byte[] array=bb.toBytes();
		if(verbose){outstream2.println("FQWST write("+array.length+")");}
		try{
			synchronized(outstream){outstream.write(array);}//Works when synchronized on stream, hangs when synchronized on this
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
		if(verbose){outstream2.println("FQWST write("+array.length+") finished");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Writer Thread         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class WriterRunnable implements Runnable{
		@Override
		public void run(){
			if(verbose){outstream2.println("WR Thread started");}
			for(ListNum<Read> job=queue.take(); job!=null && !job.poison(); job=queue.take()){
				if(verbose){outstream2.println("WR Thread got "+job.id+", "+job.last());}
				// Blocking wait for next job (or poison pill)
				writeReads(job.list);
			}
			if(verbose){outstream2.println("WR Thread finished");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeFastq(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				r1.toFastq(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r1.length();
			}
			if(writeR2 && r2!=null){
				r2.toFastq(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r2.length();
			}
		}
	}
	
	private void writeFasta(ArrayList<Read> reads, ByteBuilder bb){
		// ... (Fasta writing logic) ...
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				r1.toFasta(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r1.length();
			}
			if(writeR2 && r2!=null){
				r2.toFasta(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r2.length();
			}
		}
	}
	
	private void writeHeader(ArrayList<Read> reads, ByteBuilder bb){
		// ... (Header writing logic) ...
		for(Read r : reads){
			if(r==null){continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				bb.appendln(r1.id);
				readsWritten++;
			}
			if(writeR2 && r2!=null){
				bb.appendln(r2.id);
				readsWritten++;
			}
		}
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
	OutputStream outstream;
	/** Write R1 reads (pairnum==0) */
	final boolean writeR1;
	/** Write R2 reads (pairnum==1 or mate) */
	final boolean writeR2;
	
	/** Configuration: Enforce ordering */
	final boolean ordered;
	/** Configuration: Use internal thread */
	final boolean threaded;
	
	/** Queue for ordering/threading. Null if unordered+unthreaded. */
	final JobQueue<ListNum<Read>> queue;
	/** Internal writer thread. Null if unthreaded. */
	Thread writerThread;
	
	protected long readsWritten=0;
	protected long basesWritten=0;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** True after start() called */
	private boolean started=false;
	private boolean finished=false;
	private boolean poisoned=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final int FASTQ=FileFormat.FASTQ;
	private static final int FASTA=FileFormat.FASTA;
	private static final int HEADER=FileFormat.HEADER;
	private static final int UNKNOWN=FileFormat.UNKNOWN;
	
	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	protected PrintStream outstream2=System.err;
	
}