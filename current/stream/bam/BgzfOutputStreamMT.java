package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

import stream.JobQueue;

/**
 * Multithreaded BGZF (Blocked GZIP Format) output stream.
 *
 * Architecture:
 * - Producer (main thread): Accumulates data into 64KB blocks, creates jobs with ascending IDs
 * - Worker thread(s): Compresses blocks in parallel
 * - Writer thread: Uses JobQueue to retrieve jobs in sequential order and write to output
 *
 * Jobs are automatically ordered by JobQueue to maintain sequential output
 * even when workers complete out of order.
 *
 * @author Chloe
 * @contributor Isla
 * @date October 18, 2025
 */
public class BgzfOutputStreamMT extends OutputStream {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create a BGZF output stream with default thread count (1) and 64KB block size. */
	public BgzfOutputStreamMT(OutputStream out){
		this(out, 1, 6, DEFAULT_BLOCK_SIZE); // 1 worker, compression level 6
	}

	/** Create a BGZF output stream with custom thread count and compression level. */
	public BgzfOutputStreamMT(OutputStream out, int threads, int compressionLevel){
		this(out, threads, compressionLevel, DEFAULT_BLOCK_SIZE);
	}

	/** Create a BGZF output stream with fully custom threading and block size. */
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
		int queueSize=Math.max(workerThreads*2, 2);
		this.inputQueue=new ArrayBlockingQueue<>(queueSize);
		this.jobQueue=new JobQueue<BgzfJob>(queueSize, true, true, 0);

		this.buffer=new byte[maxBlockSize];

		startThreads();

		assert(repOK()) : "Constructor postcondition failed";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Start worker and writer threads. */
	private void startThreads(){
		assert(workers==null) : "Workers already started";
		assert(writer==null) : "Writer already started";

		// Start worker threads (NON-daemon - JVM waits for them)
		workers=new Thread[workerThreads];
		for(int i=0; i<workerThreads; i++){
			final int threadNum=i;
			workers[i]=new Thread(new Runnable(){
				public void run(){workerLoop();}
			}, "BGZF-Compressor-"+threadNum);
			workers[i].start();
		}

		// Start writer thread (NON-daemon - JVM waits for it)
		writer=new Thread(new Runnable(){
			public void run(){writerLoop();}
		}, "BGZF-Writer");
		writer.start();
	}

	@Override
	public void write(int b) throws IOException{
		assert(!closed) : "Stream closed";

		buffer[bufferPos++]=(byte)b;
		if(bufferPos>=maxBlockSize){
			flushBlock(false); // Not the last block
		}

		assert(bufferPos<maxBlockSize) : "Buffer overflow: "+bufferPos;
	}

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

	/** Submit current buffer as a job for compression. */
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

		// Submit to input queue (blocks if queue full)
		try{
			inputQueue.put(job);
			if(verbose){
				System.err.println("flushBlock: submitted job "+job.id+" to inputQueue");
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new IOException("Interrupted while submitting job");
		}

		// Allocate new buffer
		buffer=new byte[maxBlockSize];
		bufferPos=0;

		assert(bufferPos==0) : "Buffer position not reset";
	}

	/** Worker thread: Compress BGZF blocks. */
	private void workerLoop(){
		Deflater deflater=new Deflater(compressionLevel, true); // true=nowrap mode
		CRC32 crc=new CRC32();

		try{
			while(true){
				BgzfJob job=inputQueue.take();

				if(job.isPoisonPill()){
					if(verbose){
						System.err.println("Worker: received POISON, re-injecting and exiting");
					}
					inputQueue.put(BgzfJob.POISON_PILL); // Re-inject for other workers
					break;
				}

				if(verbose){
					System.err.println("Worker: processing job "+job.id+
						(job.lastJob ? " (LAST JOB)" : "")+" ("+job.decompressedSize+
						" bytes)");
				}

				// If this is the last job, inject poison for other workers
				if(job.lastJob){
					if(verbose){
						System.err.println("Worker: saw LAST JOB, injecting POISON");
					}
					inputQueue.put(BgzfJob.POISON_PILL);
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

					// Compress the data
					deflater.reset();
					deflater.setInput(job.decompressed, 0, job.decompressedSize);
					deflater.finish();

					byte[] compressed=new byte[maxBlockSize+1024];
					int compressedSize=deflater.deflate(compressed);

					assert(compressedSize>0) : "Deflate produced zero bytes";
					assert(deflater.finished()) : "Deflate not finished";

					// Calculate CRC32
					crc.reset();
					crc.update(job.decompressed, 0, job.decompressedSize);
					long crcValue=crc.getValue();

					// Calculate BSIZE (total block size minus 1)
					// Total=header(18)+compressed+footer(8)
					int bsize=18+compressedSize+8-1;

					assert(bsize>=27) : "BSIZE too small: "+bsize;
					assert(bsize<=65535) : "BSIZE exceeds uint16 max: "+bsize;

					// Build complete BGZF block
					job.compressed=new byte[bsize+1];
					int pos=0;

					// Write gzip header with BC subfield
					job.compressed[pos++]=31;          // ID1
					job.compressed[pos++]=(byte)139;   // ID2
					job.compressed[pos++]=8;           // CM (DEFLATE)
					job.compressed[pos++]=4;           // FLG (FEXTRA)
					pos=writeInt32(job.compressed, pos, 0);      // MTIME
					job.compressed[pos++]=0;           // XFL
					job.compressed[pos++]=(byte)255;   // OS (unknown)
					pos=writeInt16(job.compressed, pos, 6);      // XLEN
					job.compressed[pos++]=66;          // SI1 'B'
					job.compressed[pos++]=67;          // SI2 'C'
					pos=writeInt16(job.compressed, pos, 2);      // SLEN
					pos=writeInt16(job.compressed, pos, bsize);  // BSIZE

					assert(pos==18) : "Header size wrong: "+pos;

					// Write compressed data
					System.arraycopy(compressed, 0, job.compressed, pos, compressedSize);
					pos+=compressedSize;

					// Write footer (CRC32+ISIZE)
					pos=writeInt32(job.compressed, pos, (int)crcValue);
					pos=writeInt32(job.compressed, pos, job.decompressedSize);

					assert(pos==bsize+1) : "Total size wrong: expected "+(bsize+1)+
						", got "+pos;

					job.compressedSize=pos;

					// Clear decompressed data to free memory
					job.decompressed=null;
					job.decompressedSize=0;

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

	/** Writer thread: Write compressed blocks in sequential order. */
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
						out.write(job.compressed, 0, job.compressedSize);
					}
				}

				// If this was the last job, write EOF marker, close stream, and exit
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

	@Override
	public void flush() throws IOException{
		flush(false); // Not the last job
	}

	/** Internal flush implementation. */
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

			try{
				inputQueue.put(emptyJob);
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
				throw new IOException("Interrupted submitting empty lastJob");
			}
		}

		// Flush underlying stream
		out.flush();
	}

	@Override
	public void close() throws IOException{
		if(closed){return;} // Already closed - idempotent

		if(verbose){System.err.println("close(): starting shutdown");}

		// Flush with lastJob=true - this submits the lastJob marker
		flush(true);

		// Mark as closed
		closed=true;

		if(verbose){System.err.println("close(): waiting for writer thread to finish");}

		// Wait for writer thread to complete (ensures EOF is written and stream closed)
		if(writer!=null){
			while(writer.getState()!=Thread.State.TERMINATED){
				try{
					writer.join(1000);
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}

		// Check for errors
		if(workerError!=null){throw workerError;}

		if(verbose){System.err.println("close(): writer finished, close complete");}
	}

	/** Write BGZF EOF marker (28-byte empty block). */
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

	/** Write 16-bit integer in little-endian format. */
	private int writeInt16(byte[] buf, int pos, int val){
		buf[pos++]=(byte)(val & 0xFF);
		buf[pos++]=(byte)((val>>8) & 0xFF);
		return pos;
	}

	/** Write 32-bit integer in little-endian format. */
	private int writeInt32(byte[] buf, int pos, int val){
		buf[pos++]=(byte)(val & 0xFF);
		buf[pos++]=(byte)((val>>8) & 0xFF);
		buf[pos++]=(byte)((val>>16) & 0xFF);
		buf[pos++]=(byte)((val>>24) & 0xFF);
		return pos;
	}

	/** Validate internal state for debugging. */
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

	/** Number of worker threads (start with 1 for correctness, then scale up) */
	private final int workerThreads;
	/** Compression level (0-9, default 6) */
	private final int compressionLevel;
	/** Input queue for jobs to be compressed */
	private final ArrayBlockingQueue<BgzfJob> inputQueue;
	/** Job queue maintaining sequential output order */
	private final JobQueue<BgzfJob> jobQueue;
	/** Maximum uncompressed block size (<=64KB) */
	private final int maxBlockSize;
	/** Worker threads compressing blocks */
	private Thread[] workers;
	/** Writer thread outputting compressed blocks */
	private Thread writer;
	/** Next job ID to assign */
	private long nextJobId=0;
	/** Current accumulation buffer */
	private byte[] buffer;
	/** Position in accumulation buffer */
	private int bufferPos=0;
	/** Underlying output stream */
	private final OutputStream out;
	/** Error state from worker threads */
	private volatile IOException workerError=null;
	/** Stream closed flag (only accessed by main thread) */
	private boolean closed=false;
	/** Debug flag (enable with -Dbgzf.debug=true) */
	private static final boolean verbose=false;//Boolean.getBoolean("bgzf.debug");
	/** Default maximum uncompressed block size (64KB) */
	public static final int DEFAULT_BLOCK_SIZE=65536;
}