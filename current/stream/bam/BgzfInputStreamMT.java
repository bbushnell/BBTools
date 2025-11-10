package stream.bam;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

import stream.JobQueue;

/**
 * Multithreaded BGZF (Blocked GZIP Format) input stream.
 *
 * Architecture:
 * - Producer thread: Reads BGZF blocks from file, creates jobs with ascending IDs
 * - Worker thread(s): Decompresses blocks in parallel
 * - Consumer (main thread): Uses JobQueue to retrieve jobs in sequential order
 *
 * Jobs are automatically ordered by JobQueue to maintain sequential output
 * even when workers complete out of order.
 *
 * @author Chloe
 * @contributor Isla
 * @date October 23, 2025
 */
public class BgzfInputStreamMT extends InputStream {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public BgzfInputStreamMT(InputStream in){
		this(in, 4);
	}

	public BgzfInputStreamMT(InputStream in, int threads){
		assert(in!=null) : "Null input stream";
		assert(threads>0 && threads<=32) : "Invalid thread count: "+threads;

		this.in=in;
		this.workerThreads=threads;

		//Queue size: 3+workers*2 allows some buffering without excessive memory
		final int queueSize=3+(3*workerThreads)/2;
		this.inputQueue=new ArrayBlockingQueue<>(queueSize);
		this.jobQueue=new JobQueue<BgzfJob>(queueSize, true, true, 0);

		startThreads();

		assert(repOK()) : "Constructor postcondition failed";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Start producer and worker threads. */
	private void startThreads(){
		assert(producer==null) : "Threads already started";
		assert(workers==null) : "Workers already started";

		//Start producer thread
		producer=new Thread(new Runnable(){
			public void run(){producerLoop();}
		}, "BGZF-InputProducer");
		producer.setDaemon(true);
		producer.start();

		//Start worker threads
		workers=new Thread[workerThreads];
		for(int i=0; i<workerThreads; i++){
			final int threadNum=i;
			workers[i]=new Thread(new Runnable(){
				public void run(){workerLoop();}
			}, "BGZF-InputWorker-"+threadNum);
			workers[i].setDaemon(true);
			workers[i].start();
		}
	}

	/** Producer thread: Read BGZF blocks and create jobs. */
	private void producerLoop(){
		try{
			while(!closed){
				BgzfJob job=readNextBlock();
				if(job==null){
					if(!lastJobQueued){
						BgzfJob eofJob=new BgzfJob(nextJobId++, new byte[0], null, true);
						eofJob.decompressedSize=0;
						while(eofJob!=null) {
							try{
								inputQueue.put(eofJob);
								eofJob=null;
							}catch(InterruptedException e){
								// If closing, stop trying to enqueue and exit
								synchronized(BgzfInputStreamMT.this) {
									if(closed) {
										inputQueue.offer(eofJob);
										producerFinished=true;
										return;
									}
								}
							}
						}
						lastJobQueued=true;
					}
					break;
				}

				assert(job.compressed!=null || job.decompressed!=null) : 
					"Producer created job without data";
				if(job.compressed!=null){
					assert(job.compressedSize>0) : 
						"Producer created job with zero compressed size";
				}else{
					assert(job.decompressedSize>0) : 
						"Producer created job with zero decompressed size";
				}
				assert(!DEBUG || job.repOK()) : "Producer created invalid job";

				// Submit to input queue (be robust to interrupts during close)
				while(job!=null) {
					try{
						inputQueue.put(job);
						job=null;
					}catch(InterruptedException e){
						// If we are closing, stop trying to enqueue and exit loop
						synchronized(BgzfInputStreamMT.this) {
							if(closed) {
								job=null;
								producerFinished=true;
								return;
							}
						}
					}
				}
			}
		}catch(IOException e){
			workerError=e;
		}finally{
			producerFinished=true;
			if(!lastJobQueued){
				// Best-effort poison without blocking during shutdown
				inputQueue.offer(BgzfJob.POISON_PILL);
			}
		}
	}

	/** Read next BGZF block from input stream and create job. Returns null on EOF. */
	private BgzfJob readNextBlock() throws IOException{
		while(true){
			if(readingPlainGzip){
				BgzfJob job=readPlainGzipChunk();
				if(job!=null){return job;}
				closePlainStream(); //Plain gzip member exhausted; switch back to BGZF
				continue;
			}

			//Read gzip header (minimum 10 bytes)
			byte[] header=new byte[10];
			int bytesRead=readFully(header, 0, header.length);
			if(bytesRead==0){return null;} //EOF
			if(bytesRead<header.length){
				throw new EOFException("Truncated BGZF block header");
			}

			//Verify gzip signature
			assert((header[0] & 0xFF)==31 && (header[1] & 0xFF)==139) : 
				"Not a gzip file: "+(header[0] & 0xFF)+", "+(header[1] & 0xFF);
			if((header[0] & 0xFF)!=31 || (header[1] & 0xFF)!=139){
				throw new IOException("Not a gzip file");
			}

			//Check compression method (should be 8=DEFLATE)
			if(header[2]!=8){
				throw new IOException("Unsupported compression method: "+header[2]);
			}

			//Check flags - FEXTRA must be set for BGZF
			int flags=header[3] & 0xFF;
			boolean fextra=(flags & 0x04)!=0;
			if(!fextra){
				startPlainGzip(header, null, null);
				continue;
			}

			//Read XLEN (2 bytes, little-endian)
			byte[] xlenBytes=new byte[2];
			if(readFully(xlenBytes, 0, 2)<2){
				throw new EOFException("Truncated XLEN");
			}
			int xlen=((xlenBytes[1] & 0xFF)<<8) | (xlenBytes[0] & 0xFF);

			//Read extra field and find BC subfield
			byte[] extra=new byte[xlen];
			if(readFully(extra, 0, xlen)<xlen){
				throw new EOFException("Truncated extra field");
			}

			int bsize=findBsizeInExtra(extra, xlen);
			if(bsize<0){
				startPlainGzip(header, xlenBytes, extra);
				continue;
			}

			//Calculate compressed data length
			int alreadyRead=10+2+xlen;
			int remaining=(bsize+1)-alreadyRead;
			if(remaining<8){throw new IOException("Invalid BSIZE: "+bsize);}

			final int compressedSize=remaining-8; //Subtract CRC32 and ISIZE

			//Read complete block into job
			final byte[] compressed=new byte[compressedSize];
			if(readFully(compressed, 0, compressedSize)<compressedSize){
				throw new EOFException("Truncated compressed data");
			}

			//Read trailer (CRC32+ISIZE)
			byte[] trailer=new byte[8];
			if(readFully(trailer, 0, 8)<8){
				throw new EOFException("Truncated block trailer");
			}

			//Store trailer in job for worker validation
			ByteBuffer bb=ByteBuffer.wrap(trailer).order(ByteOrder.LITTLE_ENDIAN);
			long expectedCrc=bb.getInt() & 0xFFFFFFFFL;
			int expectedSize=bb.getInt();

			//Store expected values for worker validation
			final byte[] decompressed=new byte[8+65536]; //8 bytes metadata+max block
			ByteBuffer meta=ByteBuffer.wrap(decompressed).order(ByteOrder.LITTLE_ENDIAN);
			meta.putInt((int)expectedCrc);
			meta.putInt(expectedSize);

			if(compressedSize==0 && expectedSize==0){lastJobQueued=true;}

			BgzfJob job=new BgzfJob(nextJobId++, decompressed, compressed, lastJobQueued);
			synchronized(job){job.compressedSize=compressedSize;}
			if(DEBUG && job.lastJob){
				System.err.println("Producer: job "+job.id+" marked LAST (compressed="+
					compressedSize+", expected="+expectedSize+")");
			}

			assert(!DEBUG || job.repOK()) : "readNextBlock created invalid job";
			return job;
		}
	}

	/** Find BC subfield in gzip extra field and extract BSIZE. */
	private int findBsizeInExtra(byte[] extra, int xlen){
		int pos=0;
		while(pos+4<=xlen){
			int si1=extra[pos] & 0xFF;
			int si2=extra[pos+1] & 0xFF;
			int slen=((extra[pos+3] & 0xFF)<<8) | (extra[pos+2] & 0xFF);

			if(si1==66 && si2==67){ //'B' 'C'
				if(slen==2 && pos+6<=xlen){
					return ((extra[pos+5] & 0xFF)<<8) | (extra[pos+4] & 0xFF);
				}
			}

			pos+=4+slen;
		}
		return -1;
	}

	private void startPlainGzip(byte[] header, byte[] xlenBytes, byte[] extra) 
		throws IOException{
		if(readingPlainGzip){return;}
		byte[] prefix=buildPrefix(header, xlenBytes, extra);
		plainGzipStream=new GZIPInputStream(new PrefixedInputStream(prefix, in), 65536);
		readingPlainGzip=true;
	}

	private BgzfJob readPlainGzipChunk() throws IOException{
		if(plainGzipStream==null){return null;}

		byte[] buffer=new byte[65536];
		int total=0;
		while(total<buffer.length){
			int n=plainGzipStream.read(buffer, total, buffer.length-total);
			if(n<0){break;}
			if(n==0){continue;}
			total+=n;
			break;
		}

		if(total<=0){return null;}

		BgzfJob job=new BgzfJob(nextJobId++, buffer, null, false);
		job.decompressedSize=total;
		return job;
	}

	private void closePlainStream() throws IOException{
		if(plainGzipStream!=null){
			plainGzipStream.close();
			plainGzipStream=null;
		}
		readingPlainGzip=false;
	}

	private byte[] buildPrefix(byte[] header, byte[] xlenBytes, byte[] extra){
		int prefixLen=header.length;
		if(xlenBytes!=null){prefixLen+=xlenBytes.length;}
		if(extra!=null){prefixLen+=extra.length;}

		byte[] prefix=new byte[prefixLen];
		int pos=0;
		System.arraycopy(header, 0, prefix, pos, header.length);
		pos+=header.length;
		if(xlenBytes!=null){
			System.arraycopy(xlenBytes, 0, prefix, pos, xlenBytes.length);
			pos+=xlenBytes.length;
		}
		if(extra!=null && extra.length>0){
			System.arraycopy(extra, 0, prefix, pos, extra.length);
		}
		return prefix;
	}

	private static final class PrefixedInputStream extends InputStream{
		private final byte[] prefix;
		private int position=0;
		private final InputStream tail;

		PrefixedInputStream(byte[] prefix, InputStream tail){
			this.prefix=prefix;
			this.tail=tail;
		}

		@Override
		public int read() throws IOException{
			if(position<prefix.length){return prefix[position++] & 0xFF;}
			return tail.read();
		}

		@Override
		public int read(byte[] b, int off, int len) throws IOException{
			if(position<prefix.length){
				int toCopy=Math.min(len, prefix.length-position);
				System.arraycopy(prefix, position, b, off, toCopy);
				position+=toCopy;
				return toCopy;
			}
			return tail.read(b, off, len);
		}

		@Override
		public void close(){
			//Do not close tail stream; caller manages lifecycle
		}
	}

	/** Worker thread: Decompress BGZF blocks. */
	private void workerLoop(){
		Inflater inflater=new Inflater(true); //true=nowrap mode for raw deflate

		try{
			while(!closed){
				BgzfJob job=inputQueue.take();

				if(DEBUG){
					System.err.println("Worker: dequeued job "+job.id+
						(job.isPoisonPill() ? " (POISON)" : ""));
				}

				if(job.isPoisonPill()){
					if(DEBUG){
						System.err.println("Worker: got POISON, terminating");
					}
					// Do not re-inject poison to avoid potential blocking; close() interrupts others
					break;
				}

				synchronized(job){
					if(job.compressed==null){
						assert(job.decompressed!=null) : 
							"Job "+job.id+" missing decompressed data";
						if(DEBUG){
							System.err.println("Worker: job "+job.id+
								" already decompressed ("+job.decompressedSize+" bytes)");
						}
					}else{
						assert(!DEBUG || job.repOK()) : "Worker received invalid job";

						//Extract expected CRC and size from metadata
						ByteBuffer meta=ByteBuffer.wrap(job.decompressed).order(
							ByteOrder.LITTLE_ENDIAN);
						long expectedCrc=meta.getInt() & 0xFFFFFFFFL;
						int expectedSize=meta.getInt();

						//Decompress into same array (after metadata)
						inflater.reset();
						inflater.setInput(job.compressed, 0, job.compressedSize);

						if(DEBUG){System.err.println("Worker: inflating job "+job.id);}

						try{
							job.decompressedSize=inflater.inflate(job.decompressed, 8, 65536);
							if(DEBUG){
								System.err.println("Worker: inflated job "+job.id+" to "+
									job.decompressedSize+" bytes");
							}
						}catch(DataFormatException e){
							job.error=new IOException("Decompression failed for job "+
								job.id, e);
							if(!jobQueue.add(job)){break;}
							continue;
						}

						//Validate decompressed size
						assert(job.decompressedSize==expectedSize) : 
							"Size mismatch for job "+job.id+": expected "+expectedSize+
							", got "+job.decompressedSize;
						if(job.decompressedSize!=expectedSize){
							job.error=new IOException("Uncompressed size mismatch for job "+
								job.id+": expected "+expectedSize+", got "+
								job.decompressedSize);
							if(!jobQueue.add(job)){break;}
							continue;
						}

						//Verify CRC32
						CRC32 crc=new CRC32();
						crc.update(job.decompressed, 8, job.decompressedSize);
						long actualCrc=crc.getValue();

						assert(actualCrc==expectedCrc) : 
							"CRC32 mismatch for job "+job.id+": expected "+expectedCrc+
							", got "+actualCrc;
						if(actualCrc!=expectedCrc){
							job.error=new IOException("CRC32 mismatch for job "+job.id);
							if(!jobQueue.add(job)){break;}
							continue;
						}

						//Move decompressed data to start of array (remove metadata)
						System.arraycopy(job.decompressed, 8, job.decompressed, 0, 
							job.decompressedSize);

						//Clear compressed data to free memory
						job.compressed=null;
						job.compressedSize=0;

						assert(!DEBUG || job.repOK()) : "Worker produced invalid job";
					}
				}

				//Add to job queue
				if(!jobQueue.add(job)){break;}

				if(job.lastJob){
					if(DEBUG){
						System.err.println("Worker: job "+job.id+
							" marked LAST, injecting POISON to wake others");
					}
					// Try to wake remaining workers without blocking
					final int toSignal=Math.max(1, workerThreads-1);
					int signaled=0;
					while(signaled<toSignal){
						if(inputQueue.offer(BgzfJob.POISON_PILL)) {signaled++;}
						else {Thread.yield();}
					}
				}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
		}finally{
			inflater.end();
		}
	}

	@Override
	public int read() throws IOException{
		byte[] b=new byte[1];
		int n=read(b, 0, 1);
		return n<0 ? -1 : (b[0] & 0xFF);
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException{
		assert(b!=null) : "Null buffer";
		assert(off>=0 && len>=0 && len<=b.length-off) : 
			"Invalid offset/length: off="+off+", len="+len+", buf.length="+b.length;

		if(closed){throw new IOException("Stream closed");}
		if(len==0){return 0;}
		if(eofReached){return -1;}

		//Check for worker errors
		if(workerError!=null){throw workerError;}

		int totalRead=0;
		while(totalRead<len){
			//Need new block?
			if(currentBlockPos>=currentBlockSize){
				BgzfJob nextJob=jobQueue.take();
				if(nextJob==null){
					eofReached=true;
					return totalRead==0 ? -1 : totalRead;
				}

				synchronized(nextJob){
					if(DEBUG){
						System.err.println("Consumer: received job "+nextJob.id+
							", decompressedSize="+nextJob.decompressedSize+
							", decompressed="+(nextJob.decompressed!=null ? 
								"not null" : "NULL"));
					}

					//Check for job error
					if(nextJob.error!=null){
						throw new IOException("Decompression failed", nextJob.error);
					}

					assert(nextJob.decompressed!=null) : 
						"Job "+nextJob.id+" has null decompressed data";

					if(nextJob.lastJob){
						if(DEBUG){
							System.err.println("Consumer: job "+nextJob.id+
								" marked LAST, returning EOF");
						}
						eofReached=true;
						return totalRead==0 ? -1 : totalRead;
					}

					//Zero-length blocks (non-final) carry no data; skip and keep reading
					if(nextJob.decompressedSize==0){
						if(DEBUG){
							System.err.println("Consumer: job "+nextJob.id+
								" has 0 decompressed bytes, skipping");
						}
						continue;
					}

					assert(nextJob.decompressedSize>0) : 
						"Job "+nextJob.id+" has zero decompressed size";

					currentBlock=nextJob.decompressed;
					currentBlockSize=nextJob.decompressedSize;
				}
				currentBlockPos=0;
			}

			//Copy from current block
			int available=currentBlockSize-currentBlockPos;
			int toCopy=Math.min(available, len-totalRead);

			assert(toCopy>0) : "toCopy should be positive: "+toCopy;
			assert(currentBlockPos+toCopy<=currentBlockSize) : 
				"Copy would exceed block: pos="+currentBlockPos+", toCopy="+toCopy+
				", size="+currentBlockSize;

			System.arraycopy(currentBlock, currentBlockPos, b, off+totalRead, toCopy);
			currentBlockPos+=toCopy;
			totalRead+=toCopy;
		}

		assert(totalRead==len) : "Read wrong amount: expected "+len+", got "+totalRead;
		return totalRead;
	}

	/** Read exactly n bytes from input stream. Returns 0 on immediate EOF, n on success. */
	private int readFully(byte[] b, int off, int len) throws IOException{
		int total=0;
		while(total<len){
			int n=in.read(b, off+total, len-total);
			if(n<0){return total;}
			total+=n;
		}
		return total;
	}

	@Override
	public void close() throws IOException{
		if(closed){return;}

		closed=true;

		// Close streams early to unblock any I/O
		try {closePlainStream();} catch(IOException ignore) {}
		try {in.close();} catch(IOException ignore) {}

		// Best-effort: nudge workers via poison without blocking
		inputQueue.offer(BgzfJob.POISON_PILL);

		// Interrupt threads to break out of take()/put()
		if(producer!=null){producer.interrupt();}
		if(workers!=null){
			for(Thread worker : workers){worker.interrupt();}
		}

		// Wait briefly for threads to notice and exit
		try{
			if(producer!=null){producer.join(10);} // target ~10ms shutdown
			if(workers!=null){
				for(Thread worker : workers){worker.join(10);} // short joins; threads are daemon
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
		}
	}

	/** Validate internal state for debugging. */
	private boolean repOK(){
		if(in==null){return false;}
		if(workerThreads<=0 || workerThreads>32){return false;}
		if(inputQueue==null || jobQueue==null){return false;}
		if(currentBlockPos<0 || currentBlockSize<0){return false;}
		if(currentBlockPos>currentBlockSize){return false;}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of worker threads (start with 1 for correctness, then scale up) */
	private final int workerThreads;
	/** Input queue for jobs to be processed */
	private final ArrayBlockingQueue<BgzfJob> inputQueue;
	/** Job queue maintaining sequential output order */
	private final JobQueue<BgzfJob> jobQueue;
	/** Producer thread reading BGZF blocks */
	private Thread producer;
	/** Worker threads decompressing blocks */
	private Thread[] workers;
	/** Next job ID to assign by producer */
	private long nextJobId=0;
	/** Current decompressed block being read from */
	private byte[] currentBlock;
	/** Position in current block */
	private int currentBlockPos=0;
	/** Size of current block */
	private int currentBlockSize=0;
	/** Underlying input stream */
	private final InputStream in;
	/** Error state from worker threads */
	private volatile IOException workerError=null;
	/** Producer finished flag */
	private volatile boolean producerFinished=false;
	/** Stream closed flag */
	private volatile boolean closed=false;
	/** Whether a last-job marker has already been queued */
	private boolean lastJobQueued=false;
	/** Whether EOF has already been delivered to the caller */
	private boolean eofReached=false;
	/** Debug flag (enable with -Dbgzf.debug=true) */
	private static final boolean DEBUG=false;//Boolean.getBoolean("bgzf.debug");
	/** Plain gzip stream currently being drained by the producer (null when reading BGZF) */
	private GZIPInputStream plainGzipStream=null;
	/** True while the producer is emitting pre-decompressed plain gzip chunks */
	private boolean readingPlainGzip=false;
}
