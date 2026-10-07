package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

import stream.JobQueue;

/**
 * BGZF output stream with producer-side buffering, compression workers and one writer.
 * Construction starts the helper threads immediately. One caller owns write/flush/close
 * operations and the accumulation buffer; this class does not serialize concurrent callers.
 * Each submitted buffer becomes a job with an ascending ID. Workers create raw deflate
 * payloads and footer metadata, and the writer obtains jobs in order through JobQueue.
 * The default uncompressed block size is 65280 bytes, leaving room for BGZF framing.
 *
 * flush() submits buffered data and flushes the sink without awaiting queued compression
 * or writes. The writer owns normal final EOF emission and closing of the supplied sink.
 * close() requests a last job, joins the writer and interrupts remaining workers without
 * joining them. See the method contracts for the existing submission and interruption
 * behavior; historical shutdown comments below are not a general completion guarantee.
 *
 * @author Chloe
 * @contributor Isla
 * @date October 18, 2025
 */
public class BgzfOutputStreamMT extends OutputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts one compressor and one writer, using level 6 and 65280-byte input blocks. */
	public BgzfOutputStreamMT(OutputStream out){this(out, 1, 6, DEFAULT_BLOCK_SIZE);}

	/** Starts the requested compressors and one writer with 65280-byte input blocks. */
	public BgzfOutputStreamMT(OutputStream out, int threads, int compressionLevel){this(out, threads, compressionLevel, DEFAULT_BLOCK_SIZE);}

	/** Allocates queues and the accumulation buffer, then starts compression and writer threads.
	 * Queue capacity is 3+(3*threads)/2; the output JobQueue uses ordered, bounded mode.
	 * @param out Nonnull sink retained for writing, flushing and normal final closure
	 * @param threads Compressor count, asserted within 1 through 32; excludes the writer
	 * @param compressionLevel Deflate level, asserted within 0 through 9
	 * @param blockSize Input-buffer size, asserted within 1 through DEFAULT_BLOCK_SIZE */
	public BgzfOutputStreamMT(OutputStream out, int threads, int compressionLevel, int blockSize){
		assert(out!=null) : "Null output stream";
		assert(threads>0 && threads<=32) : "Invalid thread count: "+threads;
		assert(compressionLevel>=0 && compressionLevel<=9) : 
			"Invalid compression level: "+compressionLevel;
		assert(blockSize>0 && blockSize<=DEFAULT_BLOCK_SIZE) : 
			"Invalid BGZF block size: "+blockSize;

		this.out=out;
		this.workerThreads=threads;
		this.compressionLevel=compressionLevel;
		this.maxBlockSize=blockSize;

		// Queue sizes allow some buffering
		final int queueSize=3+(3*workerThreads)/2;
		this.inputQueue=new ArrayBlockingQueue<>(queueSize);
		this.jobQueue=new JobQueue<BgzfJob>(queueSize, true, true, 0);

		this.buffer=new byte[maxBlockSize];

		startThreads();

		assert(repOK()) : "Constructor postcondition failed";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the configured compressors and one writer; threads inherit the creator's daemon status. */
	private void startThreads(){
		assert(workers==null) : "Workers already started";
		assert(writer==null) : "Writer already started";

		// Start compression threads; daemon status is inherited from the creating thread.
		workers=new Thread[workerThreads];
		for(int i=0; i<workerThreads; i++){
			final int threadNum=i;
			workers[i]=new Thread(new Runnable(){
				public void run(){workerLoop();}
			}, "BGZF-Compressor-"+threadNum);
			workers[i].start();
		}

		// Start the writer with the same inherited daemon-status policy.
		writer=new Thread(new Runnable(){
			public void run(){writerLoop();}
		}, "BGZF-Writer");
		writer.start();
	}

	/** Buffers the low eight bits, submitting a block when full.
	 * Open state is asserted; stored worker errors are checked only when a block is submitted.
	 * @param b Value whose low byte is appended
	 * @throws IOException If block submission reports an error */
	@Override
	public void write(int b) throws IOException{
		assert(!closed) : "Stream closed";

		buffer[bufferPos++]=(byte)b;
		if(bufferPos>=maxBlockSize){
			flushBlock(false); // Not the last block
		}

		assert(bufferPos<maxBlockSize) : "Buffer overflow: "+bufferPos;
	}

	/** Copies the requested range into owned buffers and submits each full block.
	 * Null/range/open-state conditions are asserted. Checks stored worker errors even
	 * for an empty range; completion of this call does not await output completion.
	 * @param b Source bytes, not retained
	 * @param off First source position
	 * @param len Number of bytes to append
	 * @throws IOException If a recorded worker error or submission error is observed */
	@Override
	public void write(byte[] b, int off, int len) throws IOException{
		assert(b!=null) : "Null buffer";
		assert(off>=0 && len>=0 && len<=b.length-off) : 
			"Invalid offset/length: off="+off+", len="+len+", buf.length="+b.length;
		assert(!closed) : "Stream closed";

		if(verbose && len>0){
			System.err.println("write(): writing "+len+" bytes (bufferPos="+bufferPos+")");
		}

		// Check for worker errors
		if(workerError!=null){throw workerError;}

		while(len>0){
			int available=maxBlockSize-bufferPos;
			int toWrite=Math.min(available, len);

			assert(toWrite>0) : "toWrite should be positive: "+toWrite;
			assert(bufferPos+toWrite<=maxBlockSize) : 
				"Write would overflow: pos="+bufferPos+", toWrite="+toWrite;

			System.arraycopy(b, off, buffer, bufferPos, toWrite);
			bufferPos+=toWrite;
			off+=toWrite;
			len-=toWrite;

			if(bufferPos>=maxBlockSize){
				if(verbose){
					System.err.println("write(): buffer full, calling flushBlock(false)");
				}
				flushBlock(false); // Not the last block
			}
		}

		assert(bufferPos<=maxBlockSize) : "Buffer overflow: "+bufferPos;
	}

	/** Transfers a nonempty accumulation buffer to a job, then allocates its replacement.
	 * Empty buffers return immediately. Ordinary jobs use blocking put; final jobs try
	 * brief offers, then timed offers that recheck workerError. Replacement occurs only
	 * after submission. The job retains the old array without another payload copy.
	 * @param isLast Whether the submitted data job also marks final output
	 * @throws IOException If a stored error or interrupted submission is observed */
	private void flushBlock(boolean isLast) throws IOException{
		if(bufferPos==0){
			if(verbose){
				System.err.println("flushBlock(isLast="+isLast+"): bufferPos=0, nothing to flush");
			}
			return;
		}

		assert(bufferPos>0 && bufferPos<=maxBlockSize) : 
			"Invalid buffer position: "+bufferPos;

		if(verbose){
			System.err.println("flushBlock(isLast="+isLast+"): flushing "+bufferPos+
				" bytes as job "+nextJobId);
		}

		// Check for errors before submitting
		if(workerError!=null){throw workerError;}

		// Create job
		BgzfJob job=new BgzfJob(nextJobId++, buffer, null, isLast);
		synchronized(job){job.decompressedSize=bufferPos;}

		assert(!verbose || job.repOK()) : "flushBlock created invalid job";

		// Submit to input queue
		if(isLast){
			// Try brief offers first; the fallback below keeps retrying the real final block.
			boolean enqueued=false;
			for(int attempts=0; attempts<16 && !enqueued; attempts++){
				if(inputQueue.offer(job)){enqueued=true; break;}
				try{Thread.sleep(1);}catch(InterruptedException ie){Thread.currentThread().interrupt(); break;}
			}
			//[stream/bam/BgzfOutputStreamMT#001] FIXED 2026-06-20 (greenlit by Brian; adversarial-Sonnet
			//CONFIRMED HIGH before fix). OLD fallback SILENTLY DROPPED the real final block and injected an
			//EMPTY last-marker stamped jobQueue.nextID() (the writer's CURRENT gap) -> if the writer was
			//behind under output-sink backpressure, that marker dequeued AHEAD of in-flight blocks -> premature
			//writeEOF()+close() -> silently TRUNCATED .bam with a valid-looking EOF marker. NEW: never drop the
			//real block. After the quick non-blocking offers fail (transient backpressure: slow sink -> jobQueue
			//full -> workers blocked -> inputQueue full), BLOCK on the REAL job via a timed offer-loop that
			//rechecks workerError each tick -> waits out the (transient) backpressure as the writer drains the
			//sink and workers free inputQueue, but crashes LOUD if a worker/writer actually died (never hangs,
			//never drops). The real job carries isLast=true, so the writer still writes EOF AFTER it (correct
			//order). NEEDS VALIDATION: a slow-sink close (the offers-exhausted path) + byte-for-byte round-trip.
			if(!enqueued){
				boolean done=false;
				while(!done){
					if(workerError!=null){throw workerError;}//crash loud on real worker/writer death
					try{
						done=inputQueue.offer(job, 100, java.util.concurrent.TimeUnit.MILLISECONDS);
					}catch(InterruptedException ie){
						Thread.currentThread().interrupt();
						throw new IOException("Interrupted while submitting final BGZF block", ie);
					}
				}
				if(verbose){System.err.println("flushBlock: blocking-put final block "+job.id);}
			}
		}else{
			try{
				inputQueue.put(job);
				if(verbose){
					System.err.println("flushBlock: submitted job "+job.id+" to inputQueue");
				}
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
				throw new IOException("Interrupted while submitting job");
			}
		}

		// Allocate new buffer
		buffer=new byte[maxBlockSize];
		bufferPos=0;

		assert(bufferPos==0) : "Buffer position not reset";
	}

	/** Compresses jobs with thread-local raw Deflater and CRC32 instances.
	 * Publishes the deflate payload in compressed and the eight-byte CRC/ISIZE footer
	 * in decompressed. Empty jobs are passed through without compression. Exceptions
	 * caught by the loop are recorded in workerError; the deflater is ended on exit. */
	private void workerLoop(){
		Deflater deflater=new Deflater(compressionLevel, true); // true=nowrap mode
		if(FILTERED_BGZF){deflater.setStrategy(Deflater.FILTERED);}
		CRC32 crc=new CRC32();

		try{
			while(true){
				BgzfJob job;
				try{
					job=inputQueue.take();
				}catch(InterruptedException ie){
					// Graceful exit on interrupt during shutdown
					Thread.currentThread().interrupt();
					break;
				}

				if(job.isPoisonPill()){
					if(verbose){
						System.err.println("Worker: received POISON, exiting");
					}
					break;
				}

				if(verbose){
					System.err.println("Worker: processing job "+job.id+
						(job.lastJob ? " (LAST JOB)" : "")+" ("+job.decompressedSize+
						" bytes)");
				}

				// If this is the last job, inject poison for other workers (nonblocking)
				//Poison accounting (why no worker hangs at shutdown): the worker holding the lastJob
				//injects (workerThreads-1) poisons here, and close() puts exactly 1 more (~L483) =
				//workerThreads poisons total. The empty-marker lastJob path breaks WITHOUT consuming a
				//poison (~L262), leaving 1 spare (harmless); the non-empty lastJob path loops back and
				//consumes 1, balancing exactly. Either way every worker's inputQueue.take() eventually
				//returns a poison -> no worker blocks forever. (Math.max(1,..) over-injects 1 spare when
				//workerThreads==1; harmless.)
				if(job.lastJob){
					if(verbose){
						System.err.println("Worker: saw LAST JOB, injecting POISON");
					}
					for(int s=0, toSignal=Math.max(1, workerThreads-1); s<toSignal; s++){
						if(!inputQueue.offer(BgzfJob.POISON_PILL)){Thread.yield();}
					}
				}

				assert(job.decompressed!=null) : 
					"Worker received job with null decompressed data";

				synchronized(job){
					// Handle empty lastJob marker
					if(job.decompressedSize==0){
						if(verbose){
							System.err.println("Worker: received empty lastJob marker, "+
								"adding to queue without compressing");
						}
						job.compressed=new byte[0];
						job.compressedSize=0;

						if(!jobQueue.add(job)){break;}
						if(verbose){
							System.err.println("Worker: added empty lastJob to queue");
						}
						break; // Exit worker loop
					}

					assert(job.decompressedSize>0) : "Worker received job with zero size";
					assert(!verbose || job.repOK()) : "Worker received invalid job";

					// Compress into a payload buffer (no header/footer), writer will assemble
					deflater.reset();
					deflater.setInput(job.decompressed, 0, job.decompressedSize);
					deflater.finish();

					int cap=(maxBlockSize+1024);
					job.compressed=new byte[cap];
					int compressedSize=0;
					int guard=0;
					while(!deflater.finished()){
						int n=deflater.deflate(job.compressed, compressedSize, cap-compressedSize);
						if(n==0){
							// If no progress and not finished, avoid spinning forever
							if(deflater.needsInput()){break;}
							guard++;
							if(guard>4){break;}
							continue;
						}
						compressedSize+=n;
						// Grow buffer if extremely compressed output exceeds our conservative cap
						if(compressedSize==cap && !deflater.finished()){
							int newCap=cap*2;
							byte[] bigger=new byte[newCap];
							System.arraycopy(job.compressed, 0, bigger, 0, compressedSize);
							job.compressed=bigger;
							cap=newCap;
						}
					}

					assert(compressedSize>0) : "Deflate produced zero bytes";
					assert(deflater.finished()) : "Deflate not finished";

					// Calculate CRC32 for footer
					crc.reset();
					crc.update(job.decompressed, 0, job.decompressedSize);
					long crcValue=crc.getValue();

					// Store footer metadata (CRC32 + ISIZE) for writer
					byte[] meta=new byte[8];
					int mpos=0;
					mpos=writeInt32(meta, mpos, (int)crcValue);
					mpos=writeInt32(meta, mpos, job.decompressedSize);
					job.decompressed=meta;
					job.decompressedSize=mpos; // 8

					// Set compressed payload size
					job.compressedSize=compressedSize;

					assert(!verbose || job.repOK()) : "Worker produced invalid job";

					// Add to job queue
					if(!jobQueue.add(job)){break;}
				}
			}
		}catch(Exception e){
			workerError=new IOException("Worker thread failed", e);
		}finally{
			deflater.end();
		}
	}

	/** Writes each ordered job as header, deflate payload and CRC/ISIZE footer.
	 * A last job triggers EOF emission and sink closure after its payload. A null or
	 * poison result exits without that finalization. Caught exceptions are stored in
	 * workerError; header/footer scratch arrays are retained by this writer. */
	private void writerLoop(){
		try{
			while(true){
				BgzfJob job=jobQueue.take();
				if(job==null){
					if(verbose){
						System.err.println("Writer: JobQueue returned null, exiting");
					}
					break;
				}

				if(job.isPoisonPill()){//Should never happen
					if(verbose){
						System.err.println("Writer: received poison pill, exiting");
					}
					break;
				}

				if(verbose){
					System.err.println("Writer: received job "+job.id+
						(job.lastJob ? " (LAST JOB)" : "")+" (size="+job.compressedSize+
						" bytes)");
				}

				synchronized(job){
					// Write to output stream (skip if empty lastJob marker)
					if(job.compressedSize>0){
						assert(job.compressed!=null) : 
							"Writer received job with null compressed data";
						// Compute BSIZE and write header
						final int bsize=18+job.compressedSize+8-1;
						assert(bsize>=27 && bsize<=65535);
						if(writerHeader==null){writerHeader=new byte[18];}
						System.arraycopy(GZIP_HEADER_TEMPLATE, 0, writerHeader, 0, 18);
						writerHeader[16]=(byte)(bsize&0xFF);
						writerHeader[17]=(byte)((bsize>>8)&0xFF);
						out.write(writerHeader, 0, 18);
						// Write compressed payload
						out.write(job.compressed, 0, job.compressedSize);
						// Write footer from meta (CRC32 + ISIZE)
						if(writerFooter==null){writerFooter=new byte[8];}
						assert(job.decompressed!=null && job.decompressedSize==8);
						System.arraycopy(job.decompressed, 0, writerFooter, 0, 8);
						out.write(writerFooter, 0, 8);
					}
				}

				// If this was the last job, write EOF marker, close stream, and exit
				//The writer terminates on the FIRST dequeued job with lastJob==true and writes the BGZF
				//EOF marker exactly once here. Because jobQueue is ordered (heapReady only releases the
				//min id <= nextID, JobQueue.heapReady/take), in the NORMAL path the real lastJob is the
				//highest id and is dequeued strictly last -> EOF is correct. This same property is what
				//makes #001's fallback dangerous: a last-marker stamped with an EARLIER id (nextID at
				//fallback time) is dequeued early and triggers this EOF+close before later blocks ship.
				if(job.lastJob){
					if(verbose){
						System.err.println("Writer: lastJob written, writing EOF marker");
					}
					writeEOFMarker();

					// Close underlying stream - writer owns this
					try{
						out.close();
						if(verbose){
							System.err.println("Writer: closed underlying stream, exiting");
						}
					}catch(IOException e){
						workerError=e;
					}
					break; // Exit writer loop
				}
			}
		}catch(Exception e){
			workerError=new IOException("Writer thread failed", e);
		}
	}

	/** Submits buffered data and flushes the sink without waiting for queued jobs to finish.
	 * Returns immediately when already closed.
	 * @throws IOException If submission or flushing the underlying sink fails */
	@Override
	public void flush() throws IOException{flush(false);}

	/** Submits buffered data or attempts to enqueue an empty last marker, then flushes the sink.
	 * The empty-last path makes at most 16 nonblocking offers; it has no blocking fallback
	 * or explicit rejection when those offers fail. Neither path awaits writer completion.
	 * @param isLast Whether to request final output, including an empty marker when needed
	 * @throws IOException If data submission or the underlying flush reports an error */
	private void flush(boolean isLast) throws IOException{
		if(closed && !isLast){return;}

		if(verbose){
			System.err.println("flush(isLast="+isLast+"): bufferPos="+bufferPos);
		}

		// Flush current buffer
		if(bufferPos>0){
			flushBlock(isLast);
			if(verbose){
				System.err.println("flush: flushed buffer as job "+(nextJobId-1)+
					(isLast ? " (LAST)" : ""));
			}
		}else if(isLast){
			// Buffer empty but this is close() - create empty lastJob marker
			BgzfJob emptyJob=new BgzfJob(nextJobId++, new byte[0], null, true);

			if(verbose){
				System.err.println("flush: created empty lastJob (id="+emptyJob.id+")");
			}

			// Best-effort, non-blocking enqueue of last marker
			boolean enqueued=false;
			for(int attempts=0; attempts<16 && !enqueued; attempts++){
				if(inputQueue.offer(emptyJob)){enqueued=true; break;}
				try{Thread.sleep(1);}catch(InterruptedException ie){Thread.currentThread().interrupt(); break;}
			}
			//TODO: Probable bug [stream/bam/BgzfOutputStreamMT#002] - if these offers all fail,
			//flush(true) proceeds without a last marker, while close() subsequently joins the
			//writer that expects it. Source-only concern; separate validation and repair are
			//deferred, not part of this documentation pass.
		}

		// Flush underlying stream
		out.flush();
	}

//	@Override
//	public void close() throws IOException{
//		if(closed){return;} // Already closed - idempotent
//
//		if(verbose){System.err.println("close(): starting shutdown");}
//
//		// Submit last marker before marking closed, so in-flight jobs are still processed
//		flush(true);
//
//		// Now mark closed and nudge workers
//		closed=true;
//		inputQueue.offer(BgzfJob.POISON_PILL);
//		for(Thread t : workers){ if(t!=null) t.interrupt(); }
//
//		// Wait briefly for writer thread to finish
//		if(writer!=null){
//			try{writer.join(10);}catch(InterruptedException ie){Thread.currentThread().interrupt();}
//		}
//
//		// If still alive, keep nudging but do not block
//		if(writer!=null && writer.isAlive()){
//			inputQueue.offer(BgzfJob.POISON_PILL);
//		}
//
//		// If any error captured, surface it
//		if(workerError!=null){throw workerError;}
//	}
	
	/** Requests final output, marks closed, submits one worker poison and joins the writer.
	 * Already-closed calls return immediately. Interrupted poison submission or joining
	 * restores interrupt status and continues cleanup. Remaining workers are interrupted
	 * without joining; a recorded workerError is thrown afterward. Normal sink closure
	 * occurs in the writer, not directly here. See flush(boolean) for empty-marker limits.
	 * @throws IOException If final submission or flushing fails, or a stored worker error is present */
	@Override
	public void close() throws IOException{
		if(closed){return;} // Already closed - idempotent

		if(verbose){System.err.println("close(): starting shutdown");}

		// Submit last marker
		flush(true);

		closed=true;
		
		// Submit one additional poison; the worker loop also offers markers on a last job.
//		for(int i=0; i<workerThreads; i++){
		try{
			inputQueue.put(BgzfJob.POISON_PILL);
		}catch(InterruptedException e){
			Thread.currentThread().interrupt(); // Restore flag
			// If main thread is interrupted, we can't do much, but we shouldn't kill workers yet
		}
//		}
		
		// Wait for writer to finish writing everything
		if(writer!=null){
			try{
				writer.join(); // Wait for writer to finish (it exits when it sees lastJob)
			}catch(InterruptedException ie){
				Thread.currentThread().interrupt();
			}
		}
		
		// Cleanup: Interrupt workers only if they are still alive (stuck)
		for(Thread t : workers){
			if(t!=null && t.isAlive()){
				t.interrupt(); 
			}
		}

		// If any error captured, surface it
		if(workerError!=null){throw workerError;}
	}

	/** Writes the fixed 28-byte empty BGZF block, then flushes the underlying stream. */
	private void writeEOFMarker() throws IOException{
		// Standard 28-byte EOF marker
		byte[] eof=new byte[]{
			0x1f, (byte)0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
			0x00, (byte)0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
			0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00
		};

		out.write(eof);
		out.flush();
	}

	/** Writes the low 16 bits in little-endian order; returns the position after two bytes. */
	private int writeInt16(byte[] buf, int pos, int val){
		buf[pos++]=(byte)(val&0xFF);
		buf[pos++]=(byte)((val>>8)&0xFF);
		return pos;
	}

	/** Writes all 32 bits in little-endian order; returns the position after four bytes. */
	private int writeInt32(byte[] buf, int pos, int val){
		buf[pos++]=(byte)(val&0xFF);
		buf[pos++]=(byte)((val>>8)&0xFF);
		buf[pos++]=(byte)((val>>16)&0xFF);
		buf[pos++]=(byte)((val>>24)&0xFF);
		return pos;
	}

	/** Precomputed 18-byte gzip header template with BC subfield; BSIZE patched per block. */
	private static final byte[] GZIP_HEADER_TEMPLATE=new byte[]{
		31, (byte)139, 8, 4, // ID1, ID2, CM, FLG(FEXTRA)
		0, 0, 0, 0,          // MTIME (0)
		0, (byte)255,        // XFL=0, OS=255
		6, 0,                // XLEN=6
		66, 67,              // SI1='B', SI2='C'
		2, 0,                // SLEN=2
		0, 0                 // BSIZE (patched)
	};

	/** Copies the template, patches the low 16 BSIZE bits and writes 18 bytes to dest. */
	private static void writeHeaderWithBsize(byte[] dest, int bsize){
		// NUMA-friendly: allocate a local header copy in the calling thread
		final byte[] header=Arrays.copyOf(GZIP_HEADER_TEMPLATE, GZIP_HEADER_TEMPLATE.length);
		header[16]=(byte)(bsize&0xFF);
		header[17]=(byte)((bsize>>8)&0xFF);
		System.arraycopy(header, 0, dest, 0, 18);
	}

	/** Checks configured ranges, required references and accumulation-buffer bounds.
	 * Does not inspect job payloads, thread completion, queued work or sink state. */
	private boolean repOK(){
		if(out==null){return false;}
		if(workerThreads<=0 || workerThreads>32){return false;}
		if(compressionLevel<0 || compressionLevel>9){return false;}
		if(inputQueue==null || jobQueue==null){return false;}
		if(buffer==null || buffer.length!=maxBlockSize){return false;}
		if(bufferPos<0 || bufferPos>maxBlockSize){return false;}
		if(maxBlockSize<=0 || maxBlockSize>DEFAULT_BLOCK_SIZE){return false;}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of compressor threads; excludes the separate writer. */
	private final int workerThreads;
	/** Compression level (0-9, default 6) */
	private final int compressionLevel;
	/** Input queue for jobs to be compressed */
	private final ArrayBlockingQueue<BgzfJob> inputQueue;
	/** Job queue maintaining sequential output order */
	private final JobQueue<BgzfJob> jobQueue;
	/** Configured input-block capacity, at most DEFAULT_BLOCK_SIZE. */
	private final int maxBlockSize;
	/** Worker threads compressing blocks */
	private Thread[] workers;
	/** Writer thread outputting compressed blocks */
	private Thread writer;
	/** Next job ID to assign */
	private long nextJobId=0;
	/** Producer-owned accumulation buffer, transferred to a job on submission. */
	private byte[] buffer;
	/** Position in accumulation buffer */
	private int bufferPos=0;
	/** Retained sink; writer emits data/closes it, while caller-side flush also flushes it. */
	private final OutputStream out;
	/** Most recently assigned IOException from a compressor or writer; null initially. */
	private volatile IOException workerError=null;
	/** Per-writer small reusable header buffer (18B) */
	private byte[] writerHeader;
	/** Per-writer small reusable footer buffer (8B) */
	private byte[] writerFooter;
	/** Caller-side close latch, set after final flush and before waiting for helpers. */
	private volatile boolean closed=false;
	/** Compile-time disabled diagnostics; the system-property expression is commented out. */
	private static final boolean verbose=false;//Boolean.getBoolean("bgzf.debug");
	/** Default input-block capacity, 65280 bytes. */
	//[block-size lead FIXED 2026-06-20 (greenlit)] 0xff00 (65280), NOT 65536: the BGZF spec/samtools cap the
	//UNCOMPRESSED block here so the compressed block + 26B overhead always fits the 16-bit BSIZE field. At
	//65536 an incompressible full block deflates to ~65556 (empirically) -> bsize>65535 -> wraps -> corrupt
	//framing. Read buffers stay 65536 (other writers may emit up to the spec max), only WRITE size is capped.
	public static final int DEFAULT_BLOCK_SIZE=65280;
	/** Mutable strategy preference sampled by each compressor at startup; set before construction. */
	public static boolean FILTERED_BGZF=false;
}
