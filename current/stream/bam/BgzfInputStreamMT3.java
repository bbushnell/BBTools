package stream.bam;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

import stream.OrderedQueueSystem;
import structures.BinaryByteWrapperLE;
import template.ThreadWaiter;

/**
 * BGZF reader using OrderedQueueSystem's input FIFO and ordered output queue.
 * Construction starts a daemon producer and daemon decompression workers. One caller
 * owns read state and closure; synchronized close does not serialize concurrent reads.
 * The producer assigns ascending job IDs, workers submit decoded results, and the caller
 * consumes them in order. Ordinary gzip fallback uses a replay prefix over borrowed input.
 *
 * Bulk reads try to fill the requested length across jobs until EOF. Stored producer IOExceptions,
 * normalized worker failures and per-job errors are checked at documented points. close() closes
 * inputs best-effort, signals OQS completion and waits for workers through ThreadWaiter;
 * it does not join the producer. Historical repair/validation notes are retained.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 15, 2025
 */
public class BgzfInputStreamMT3 extends InputStream{
	
	/** Reads a file, optionally copies decoded bytes to stdout, and reports bulk-call statistics.
	 * Arguments are filename, optional workers, any third argument to enable writing, and
	 * optional maximum counted calls. The fourth value limits counted/output bulk reads,
	 * not physical BGZF blocks; the loop reads before checking that limit, so an additional
	 * call can consume bytes without counting or writing them. Missing filename exits 1.
	 * @throws IOException If opening, reading or closing through the declared APIs fails */
	public static void main(String[] args) throws IOException{
		if(args.length<1){
			System.err.println("Usage: BgzfInputStreamMT3 <file.gz>");
			System.exit(1);
		}
		int threads=(args.length>1 ? Integer.parseInt(args[1]) : BgzfSettings.READ_THREADS);
		boolean write=args.length>2;
		long maxBlocks=(args.length>3 ? Long.parseLong(args[3]) : Long.MAX_VALUE);
		
		
		String filename=args[0];
		byte[] buffer=new byte[131072];
		
		long totalReads=0;
		long totalBytes=0;
		long startTime=System.nanoTime();
		
		try(InputStream fis=new java.io.FileInputStream(filename);
			InputStream bis=new BufferedInputStream(fis, 131072);
			InputStream bgzf=new BgzfInputStreamMT3(bis, threads)){
			
			int bytesRead;
			while((bytesRead=bgzf.read(buffer))>=0 && totalReads<maxBlocks){
				totalReads++;
				totalBytes+=bytesRead;
				if(write){System.out.write(buffer, 0, bytesRead);}
			}
		}
		long endTime=System.nanoTime();
		float seconds=(endTime-startTime)/1e9f;
		
		System.err.println("Read operations: "+totalReads);
		System.err.println("Total bytes:     "+totalBytes);
		System.err.println("Time:            "+String.format("%.3f", seconds)+" seconds");
		System.err.println("Throughput:      "+String.format("%.2f", totalBytes/seconds/1e6)+" MB/s");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts helpers using the current BgzfSettings.READ_THREADS value. */
	public BgzfInputStreamMT3(InputStream in){this(in, BgzfSettings.READ_THREADS);}

	/** Retains input, creates marker prototypes and OQS, then starts one producer plus workers.
	 * OQS uses input capacity threads+4 and output backpressure (3*threads)/2+4.
	 * @param in Nonnull input retained for producer reads and close()
	 * @param threads Decompression worker count, asserted within 1 through 32 */
	public BgzfInputStreamMT3(InputStream in, int threads){
		assert(in!=null) : "Null input stream";
		assert(threads>0 && threads<=32) : "Invalid thread count: "+threads;

		this.in=in;
		this.workerThreads=threads;

		// Create OQS with prototypes
		BgzfInputJob inputPrototype=new BgzfInputJob(0, null, 0, 0, false);
		BgzfInputJob outputPrototype=new BgzfInputJob(0, null, 0, 0, false);
		this.oqs=new OrderedQueueSystem<BgzfInputJob, BgzfInputJob>(
			threads, true, inputPrototype, outputPrototype);

		startThreads();

		assert(repOK()) : "Constructor postcondition failed";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the producer and all workers, explicitly marking every helper daemon. */
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

	/** Reads jobs using reusable framing scratch and submits them to the OQS input FIFO.
	 * Caught IOExceptions are stored in workerError. Finally marks producerFinished and
	 * calls OQS poison() to publish output LAST and input poison markers. */
	private void producerLoop(){
		// Declare temp arrays outside loop for reuse
		final byte[] fields=new byte[10];
		final byte[] xlenBytes=new byte[2];
		final byte[] extra=new byte[1024];
		
		try{
			while(!closed){
				BgzfInputJob job=readNextBlock(fields, xlenBytes, extra);
				if(job==null){
					break;
				}

				assert(job.compressed!=null || job.decompressed!=null) : 
					"Producer created job without data";
				if(job.compressed!=null){
					assert(job.compressed.length>0) : 
						"Producer created job with zero compressed size";
				}else{
					assert(job.decompressedSize>0) : 
						"Producer created job with zero decompressed size";
				}
				assert(!DEBUG || job.repOK()) : "Producer created invalid job";

				oqs.addInput(job);
			}
		}catch(IOException e){
			workerError=e;
		}finally{
			producerFinished=true;
			oqs.poison();
		}
	}

	/** Reads framing into a compressed BGZF job or returns a predecoded ordinary-gzip job.
	 * Absent FEXTRA or no usable BC field selects JDK gzip fallback. CRC/ISIZE metadata
	 * is stored in the job, with inflation left to workers. Both compressed size and
	 * expected size being zero mark the returned job LAST. Returns null at exhaustion;
	 * selected framing and truncation checks throw, without promising full validation.
	 * @param fields Ten-byte header scratch, also used for the eight-byte trailer
	 * @param xlenBytes Two-byte XLEN scratch
	 * @param extra Extra-field scratch; any larger replacement remains local */
	private BgzfInputJob readNextBlock(byte[] fields, byte[] xlenBytes, byte[] extra)
			throws IOException{
		while(true){
			if(readingPlainGzip){
				BgzfInputJob job=readPlainGzipChunk();
				if(job!=null){return job;}
				closePlainStream();
				continue;
			}

			//Read gzip header (minimum 10 bytes)
			int bytesRead=readFully(fields, 0, fields.length);
			if(bytesRead==0){return null;} //EOF
			if(bytesRead<fields.length){
				throw new EOFException("Truncated BGZF block header");
			}

			//Verify gzip signature
			assert((fields[0]&0xFF)==31 && (fields[1]&0xFF)==139) :
				"Not a gzip file: "+(fields[0]&0xFF)+", "+(fields[1]&0xFF);
			if((fields[0]&0xFF)!=31 || (fields[1]&0xFF)!=139){
				throw new IOException("Not a gzip file");
			}

			//Check compression method (should be 8=DEFLATE)
			if(fields[2]!=8){
				throw new IOException("Unsupported compression method: "+fields[2]);
			}

			//Check flags - FEXTRA must be set for BGZF
			int flags=fields[3]&0xFF;
			boolean fextra=(flags&0x04)!=0;
			if(!fextra){
				startPlainGzip(fields, null, null, 0);
				continue;
			}

			//Read XLEN (2 bytes, little-endian)
			if(readFully(xlenBytes, 0, 2)<2){
				throw new EOFException("Truncated XLEN");
			}
			int xlen=((xlenBytes[1]&0xFF)<<8)|(xlenBytes[0]&0xFF);

			//Read extra field and find BC subfield
			if(xlen>extra.length){extra=new byte[xlen];}
			if(readFully(extra, 0, xlen)<xlen){
				throw new EOFException("Truncated extra field");
			}

			int bsize=findBsizeInExtra(extra, xlen);
			if(bsize<0){
				startPlainGzip(fields, xlenBytes, extra, xlen);
				continue;
			}

			//Calculate compressed data length
			int alreadyRead=10+2+xlen;
			int remaining=(bsize+1)-alreadyRead;
			if(remaining<8){throw new IOException("Invalid BSIZE: "+bsize);}

			final int compressedSize=remaining-8;

			//Read compressed data
			final byte[] compressed=new byte[compressedSize];
			if(readFully(compressed, 0, compressedSize)<compressedSize){
				throw new EOFException("Truncated compressed data");
			}

			//Read trailer (CRC32+ISIZE)
			if(readFully(fields, 0, 8)<8){
				throw new EOFException("Truncated block trailer");
			}

			//Extract expected CRC and size using wrapper
			BinaryByteWrapperLE wrapper=new BinaryByteWrapperLE(fields);
			long expectedCrc=wrapper.getInt()&0xFFFFFFFFL;
			int expectedSize=wrapper.getInt();

			boolean isLast=(compressedSize==0 && expectedSize==0);

			BgzfInputJob job=new BgzfInputJob(nextJobId++, compressed, 
				expectedCrc, expectedSize, isLast);
			
			if(verbose && job.lastJob){
				System.err.println("Producer: job "+job.id+" marked LAST (compressed="+
					compressedSize+", expected="+expectedSize+")");
			}

			assert(!DEBUG || job.repOK()) : "readNextBlock created invalid job";
			return job;
		}
	}
	
	/** Scans subfields with six bytes available, returning the first length-two BC value or -1. */
	private final int findBsizeInExtra(final byte[] extra, final int xlen){
		for(int pos=0, lim=xlen-6; pos<=lim;){
			final int slen=((extra[pos+3]&0xFF)<<8)|(extra[pos+2]&0xFF);
			if(extra[pos]=='B' && extra[pos+1]=='C' && slen==2){
				return ((extra[pos+5]&0xFF)<<8)|(extra[pos+4]&0xFF);
			}
			pos+=4+slen;
		}
		return -1;
	}

	/** Starts JDK gzip decoding by replaying consumed framing ahead of the borrowed input.
	 * Returns unchanged if fallback is already active; the replay wrapper leaves in open. */
	private void startPlainGzip(byte[] header, byte[] xlenBytes, byte[] extra, int xlen)
		throws IOException{
		if(readingPlainGzip){return;}
		byte[] prefix=buildPrefix(header, xlenBytes, extra, xlen);
		plainGzipStream=new GZIPInputStream(new PrefixedInputStream(prefix, in), 65536);
		readingPlainGzip=true;
	}

	/** Copies header, optional XLEN bytes and exactly xlen extra bytes into one replay prefix. */
	private byte[] buildPrefix(byte[] header, byte[] xlenBytes, byte[] extra, final int xlen){
		int prefixLen=header.length;
		if(xlenBytes!=null){prefixLen+=xlenBytes.length;}
		prefixLen+=xlen;

		byte[] prefix=new byte[prefixLen];
		int pos=0;
		System.arraycopy(header, 0, prefix, pos, header.length);
		pos+=header.length;
		if(xlenBytes!=null){
			System.arraycopy(xlenBytes, 0, prefix, pos, xlenBytes.length);
			pos+=xlenBytes.length;
		}
		if(xlen>0){
			System.arraycopy(extra, 0, prefix, pos, xlen);
		}
		return prefix;
	}

	/** Returns the first positive gzip read in a fresh 65536-byte buffer as a predecoded job.
	 * Returns null without a decoder, at exhaustion, or when a caught exception coincides
	 * with closed=true. Other read exceptions propagate. */
	private BgzfInputJob readPlainGzipChunk() throws IOException{
		if(plainGzipStream==null){return null;}

		final byte[] buffer=new byte[65536];
		int total=0;
		while(total<buffer.length){
			int n=0;
			try{//Protects from a closing race condition
				n=plainGzipStream.read(buffer, total, buffer.length-total);
			}catch(Exception e){
				if(closed){return null;} // Expected shutdown error
				throw e; // Real error
			}
			if(n<0){break;}
			if(n==0){continue;}
			total+=n;
			break;
		}

		if(total<=0){return null;}

		BgzfInputJob job=new BgzfInputJob(nextJobId++, null, 0, 0, false);
		synchronized(job){
			job.decompressed=buffer;
			job.decompressedSize=total;
		}
		return job;
	}

	/** Closes the fallback decoder and clears its fields on success; the wrapper leaves in open. */
	private void closePlainStream() throws IOException{
		if(plainGzipStream!=null){
			plainGzipStream.close();
			plainGzipStream=null;
		}
		readingPlainGzip=false;
	}

	/** Replays borrowed prefix bytes before delegating to borrowed input, without owning its closure. */
	private static final class PrefixedInputStream extends InputStream{
		/** Borrowed already-consumed framing bytes. */
		private final byte[] prefix;
		/** Next unread prefix position. */
		private int position=0;
		/** Borrowed underlying input, closed by the outer reader. */
		private final InputStream tail;

		/** Retains prefix and tail without copying or closing them. */
		PrefixedInputStream(byte[] prefix, InputStream tail){
			this.prefix=prefix;
			this.tail=tail;
		}

		/** Reports unread prefix and tail bytes, saturating the sum to avoid overflow.
		 * Older GZIPInputStream implementations use available() to find concatenated members;
		 * inherited zero can discard a buffered next-member header and cause false corruption.
		 */
		@Override
		public int available() throws IOException{
			assert(position>=0 && position<=prefix.length) : "Gzip header replay position outside prefix: "+position+" of "+prefix.length;
			return (int)Math.min(Integer.MAX_VALUE, (long)(prefix.length-position)+tail.available());
		}

		/** Returns an unsigned prefix byte, or delegates when the prefix is exhausted. */
		@Override
		public int read() throws IOException{
			if(position<prefix.length){return prefix[position++]&0xFF;}
			return tail.read();
		}

		/** Copies only remaining prefix bytes when present, otherwise delegates to tail. */
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

		/** Leaves the borrowed tail open for the outer reader to close. */
		@Override
		public void close(){}//Do not close tail stream; caller manages lifecycle
	}

	/** Takes FIFO jobs, handles input poison by re-enqueuing it, and submits decoded results to OQS.
	 * Compressed jobs inflate into fresh 65536-byte arrays; predecoded jobs bypass inflation.
	 * Size/CRC assertions precede explicit mismatch branches. Those branches and caught
	 * DataFormatException attach per-job errors when reached. The outer catch normalizes
	 * other caught failures into workerError and calls setFinished(true); finally ends the inflater. */
	private void workerLoop(){
		final Inflater inflater=new Inflater(true);

		try{
			while(!closed){
				BgzfInputJob job=oqs.getInput();

				if(verbose){
					System.err.println("Worker: dequeued job "+job.id+
						(job.isPoisonPill() ? " (POISON)" : ""));
				}

				if(job.isPoisonPill()){
					if(verbose){
						System.err.println("Worker: got POISON, terminating");
					}
					//Re-inject for other workers
					oqs.addInput(job);
					break;
				}

				if(job.compressed==null){
					assert(job.decompressed!=null || job.last()) : 
						"Job "+job.id+" missing decompressed data";
					if(verbose){
						System.err.println("Worker: job "+job.id+
							" already decompressed ("+job.decompressedSize+" bytes)");
					}
				}else{
					assert(!DEBUG || job.repOK()) : "Worker received invalid job";

					//Decompress
					inflater.reset();
					inflater.setInput(job.compressed, 0, job.compressed.length);

					if(verbose){System.err.println("Worker: inflating job "+job.id);}
					final byte[] decompressed=new byte[65536];
					
					try{
						synchronized(job){
							job.decompressedSize=inflater.inflate(decompressed, 0, 65536);
							job.decompressed=decompressed;
						}
						if(verbose){
							System.err.println("Worker: inflated job "+job.id+" to "+
								job.decompressedSize+" bytes");
						}
					}catch(DataFormatException e){
						job.error=new IOException("Decompression failed for job "+
							job.id, e);
						oqs.addOutput(job);
						continue;
					}

					//Validate decompressed size
					assert(job.decompressedSize==job.expectedSize) : 
						"Size mismatch for job "+job.id+": expected "+job.expectedSize+
						", got "+job.decompressedSize;
					if(job.decompressedSize!=job.expectedSize){
						job.error=new IOException("Uncompressed size mismatch for job "+
							job.id+": expected "+job.expectedSize+", got "+
							job.decompressedSize);
						oqs.addOutput(job);
						continue;
					}

					//Verify CRC32
					CRC32 crc=new CRC32();
					crc.update(job.decompressed, 0, job.decompressedSize);
					long actualCrc=crc.getValue();

					assert(actualCrc==job.expectedCrc) : 
						"CRC32 mismatch for job "+job.id+": expected "+job.expectedCrc+
						", got "+actualCrc;
					if(actualCrc!=job.expectedCrc){
						job.error=new IOException("CRC32 mismatch for job "+job.id);
						oqs.addOutput(job);
						continue;
					}

					assert(!DEBUG || job.repOK()) : "Worker produced invalid job";
				}

				//Add to output queue
				oqs.addOutput(job);
			}
		}catch(Throwable e){
			//[stream/bam/BgzfInputStreamMT3#002] FIXED 2026-06-20 (greenlit, MT3 kept+fixed). Was
			//`catch(Exception){printStackTrace}` which SWALLOWED an unexpected worker death (no workerError, no
			//terminal) -> the ordered consumer could park forever on the dead worker's job id. Now: record the
			//error + force OQS shutdown so the consumer WAKES and throws LOUD (via #001's getOutput-null
			//workerError check). Throwable (not Exception) also covers AssertionError/Error, matching the MT2
			//producer fix. The inline DataFormatException path (L407) is unchanged - it addOutputs job.error in
			//order, the precise/preferred route; this catch is the backstop for everything else.
			workerError=(e instanceof IOException ? (IOException)e : new IOException("BGZF worker failed", e));
			oqs.setFinished(true);
		}finally{
			inflater.end();
		}
	}

	/** Reads through a fresh one-byte array, returning an unsigned byte or -1 at EOF. */
	@Override
	public int read() throws IOException{
		byte[] b=new byte[1];
		int n=read(b, 0, 1);
		return n<0 ? -1 : (b[0]&0xFF);
	}

	/** Fills the requested range across ordered outputs, or returns the bytes preceding EOF.
	 * Asserts buffer/range validity; closed is checked before zero length, then latched EOF
	 * before workerError. Null outputs recheck workerError; the LAST branch does not repeat
	 * that check. Both terminal-return branches latch EOF and call OQS setFinished(true).
	 * Per-job errors precede payload adoption; empty nonfinal jobs are skipped and LAST
	 * payload is not adopted.
	 * @param b Destination buffer
	 * @param off First destination position
	 * @param len Requested byte count
	 * @return len on a full read, a partial count at EOF, or -1 if EOF precedes any bytes
	 * @throws IOException If closed or an observed shared/per-job error is reported */
	@Override
	public int read(byte[] b, int off, int len) throws IOException{
		assert(b!=null) : "Null buffer";
		assert(off>=0 && len>=0 && len<=b.length-off) : 
			"Invalid offset/length: off="+off+", len="+len+", buf.length="+b.length;

		if(closed){throw new IOException("Stream closed");}
		if(len==0){return 0;}
		if(eofReached){return -1;}

		if(workerError!=null){throw workerError;}

		int totalRead=0;
		while(totalRead<len){
			if(currentBlockPos>=currentBlockSize){
				BgzfInputJob nextJob=oqs.getOutput();
				if(nextJob==null){
					//[stream/bam/BgzfInputStreamMT3#001] FIXED + VALIDATED 2026-06-20 (greenlit - Brian wants MT3
					//KEPT+fixed to transition to OQS). Re-check workerError on the getOutput-null path and throw
					//LOUD instead of returning silent EOF (mirrors the MT2 fix's consumer half). Now OBSERVABLE
					//and confirmed: once #003 (the shared-OQS poison() deadlock) was fixed, MT3 on a truncated
					//.gz crashes loud here (EOFException, no hang); a correct .gz reads clean. See #003.
					if(workerError!=null){throw workerError;}
					eofReached=true;
					oqs.setFinished(true);
					return totalRead==0 ? -1 : totalRead;
				}

				if(verbose){
					System.err.println("Consumer: received job "+nextJob.id+
						", decompressedSize="+nextJob.decompressedSize+
						", decompressed="+(nextJob.decompressed!=null ? 
							"not null" : "NULL"));
				}

				if(nextJob.error!=null){
					throw new IOException("Decompression failed", nextJob.error);
				}

				assert(nextJob.decompressed!=null || nextJob.lastJob) : 
					"Job "+nextJob.id+" has null decompressed data";

				if(nextJob.lastJob){
					if(verbose){
						System.err.println("Consumer: job "+nextJob.id+
							" marked LAST, returning EOF");
					}
					eofReached=true;
					oqs.setFinished(true);
					return totalRead==0 ? -1 : totalRead;
				}

				if(nextJob.decompressedSize==0){
					if(verbose){
						System.err.println("Consumer: job "+nextJob.id+
							" has 0 decompressed bytes, skipping");
					}
					continue;
				}

				assert(nextJob.decompressedSize>0) : 
					"Job "+nextJob.id+" has zero decompressed size";

				currentBlock=nextJob.decompressed;
				currentBlockSize=nextJob.decompressedSize;
				currentBlockPos=0;
			}

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

	/** Accumulates up to len bytes from in, returning zero or a partial count if EOF arrives first. */
	private int readFully(byte[] b, int off, int len) throws IOException{
		int total=0;
		while(total<len){
			int n=in.read(b, off+total, len-total);
			if(n<0){return total;}
			total+=n;
		}
		return total;
	}

	/** Under this monitor, latches closed, attempts stream closure and signals OQS completion.
	 * Stream-close IOExceptions are ignored. ThreadWaiter then waits for the worker array
	 * without a timeout, retrying interrupted joins; it does not join the producer. This
	 * method does not explicitly interrupt helpers or report workerError. Repeated calls
	 * return once they acquire the monitor and observe closed. */
	@Override
	public synchronized void close() throws IOException{
		if(closed){return;}
		closed=true;//[#003] part 1 FIXED: was `if(verbose)\n closed=true;` - a brace-less `if(verbose)` DANGLED
		//over this assignment, so `closed` was NEVER set (verbose=false) -> close() non-idempotent + the workers'
		//`while(!closed)` never terminated via the flag. The close()-time DEADLOCK itself (#003 part 2) is fixed
		//in shared OrderedQueueSystem.poison() (it no longer holds the OQS monitor across a blocking enqueue, so
		//this setFinished() below can always acquire it). Both validated: truncated .gz now crashes loud, not hangs.
		if(verbose){System.err.println("Called close.");}
		try{closePlainStream();}catch(IOException ignore){}
		try{in.close();}catch(IOException ignore){}
		if(verbose){System.err.println("Calling setFinished.");}
		oqs.setFinished(true);
		if(verbose){System.err.println("Wiating for workers.");}
		ThreadWaiter.waitForThreadsToFinish(workers);//Not strictly needed for daemons
		if(verbose){System.err.println("Close finished.");}
	}

	/** Checks required references, worker-count range and current block counters only. */
	private boolean repOK(){
		if(in==null){return false;}
		if(workerThreads<=0 || workerThreads>32){return false;}
		if(oqs==null){return false;}
		if(currentBlockPos<0 || currentBlockSize<0){return false;}
		if(currentBlockPos>currentBlockSize){return false;}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Decompression worker count, excluding the producer. */
	private final int workerThreads;
	/** Input FIFO, ordered output queue and external completion signaling. */
	private final OrderedQueueSystem<BgzfInputJob, BgzfInputJob> oqs;
	/** Daemon producer reading framing or ordinary gzip chunks. */
	private Thread producer;
	/** Daemon workers included in close()'s wait. */
	private Thread[] workers;
	/** Next job ID assigned sequentially by the producer. */
	private long nextJobId=0;
	/** Decoded buffer currently consumed by the caller. */
	private byte[] currentBlock;
	/** Next unread position in currentBlock. */
	private int currentBlockPos=0;
	/** Used bytes in currentBlock. */
	private int currentBlockSize=0;
	/** Retained input, closed by the outer reader rather than the replay wrapper. */
	private final InputStream in;
	/** Caught producer IOException or normalized worker failure; per-job errors are separate. */
	private volatile IOException workerError=null;
	/** Producer-exit observation, separate from consumer EOF state. */
	private volatile boolean producerFinished=false;
	/** Close latch also consulted by producer and worker loops. */
	private volatile boolean closed=false;
	/** Caller-side latch set by the null/LAST return paths. */
	private boolean eofReached=false;
	/** Compile-time disabled diagnostic output. */
	private static final boolean verbose=false;
	/** Compile-time disabled payload representation assertions. */
	private static final boolean DEBUG=false;
	/** Fallback decoder while ordinary gzip is being drained, otherwise null. */
	private GZIPInputStream plainGzipStream=null;
	/** Whether the producer is using the ordinary-gzip fallback. */
	private boolean readingPlainGzip=false;
}
