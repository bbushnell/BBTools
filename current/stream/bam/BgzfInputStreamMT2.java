package stream.bam;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

import dna.Data;
import fileIO.ReadWrite;
import stream.JobQueue;
import structures.BinaryByteWrapperLE;

/**
 * BGZF reader with daemon producer/workers and ordered, caller-side consumption.
 * Construction starts helpers immediately and asserts the existing native-input policy.
 * One caller owns read counters, the current block and closure; concurrent consumers
 * are not serialized. The producer assigns ascending IDs and retains payload arrays
 * in BgzfInputJob objects; workers supply decoded buffers through an ordered JobQueue.
 *
 * Ordinary gzip fallback replays consumed framing bytes through GZIPInputStream.
 * Bulk reads fill across jobs until the requested length or EOF, counting bytes copied
 * into caller buffers. Producer failures have a separate shared error slot checked
 * on entry and on null queue results; explicit worker errors travel with their jobs.
 * close() attempts stream closure, signals/interrupts helpers and uses short joins;
 * it is not a guarantee that every helper has exited. Historical fix notes are retained.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 14, 2025
 */
public class BgzfInputStreamMT2 extends InputStream{
	
	/** Drains a file, optionally writing decoded bytes to stdout, and reports I/O call statistics.
	 * Arguments are filename, optional worker count, and any optional third argument to
	 * enable writing (its value is not parsed). Missing filename exits 1. Read operations
	 * count bulk calls, not sequence records; timing includes optional output writes.
	 * @throws IOException If opening, reading or closing through the declared APIs fails */
	public static void main(String[] args) throws IOException{
		if(args.length<1){
			System.err.println("Usage: BgzfInputStreamMT2 <file.gz>");
			System.exit(1);
		}
		int threads=(args.length>1 ? Integer.parseInt(args[1]) : BgzfSettings.READ_THREADS);
		boolean write=args.length>2;
		
		String filename=args[0];
		byte[] buffer=new byte[131072];
		
		long totalReads=0;
		long totalBytes=0;
		long startTime=System.nanoTime();
		
		try(InputStream fis=new java.io.FileInputStream(filename);
			InputStream bis=new BufferedInputStream(fis, 131072);//reduces sys time (10%+)
			InputStream bgzf=new BgzfInputStreamMT2(bis, threads);
//			OutputStream bos=(write ? new BufferedOutputStream(System.out, 131072) : null)//slower
				){
			
			int bytesRead;
			while((bytesRead=bgzf.read(buffer))>=0){
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
	public BgzfInputStreamMT2(InputStream in){this(in, BgzfSettings.READ_THREADS);}

	/** Retains input, allocates queues and starts one producer plus the requested workers.
	 * Capacity is 3+(3*threads)/2. Existing policy assertions require ALLOW_NATIVE_BGZF
	 * and {@code PREFER_NATIVE_BGZF_IN || !Data.BGZIP()}.
	 * @param in Nonnull input retained for producer reads and close()
	 * @param threads Worker count, asserted within 1 through 32 */
	public BgzfInputStreamMT2(InputStream in, int threads){
		assert(in!=null) : "Null input stream";
		assert(threads>0 && threads<=32) : "Invalid thread count: "+threads;
		assert(ReadWrite.ALLOW_NATIVE_BGZF);
		assert(ReadWrite.PREFER_NATIVE_BGZF_IN || !Data.BGZIP()) : Data.BGZIP()+", "+Data.BGZIP_THREADED();
		this.in=in;
		this.workerThreads=threads;

		//Queue capacity is 3+(3*workers)/2, with integer arithmetic.
		final int queueSize=3+(3*workerThreads)/2;
		this.inputQueue=new ArrayBlockingQueue<>(queueSize);
		this.jobQueue=new JobQueue<BgzfInputJob>(queueSize, true, true, 0);
		
		startThreads();

		assert(repOK()) : "Constructor postcondition failed";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the producer and workers, explicitly marking every helper daemon. */
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

	/** Reads jobs with reusable framing scratch and publishes them to the worker queue.
	 * Supplies an empty LAST job at exhaustion when needed. Caught producer failures are
	 * normalized to IOException; finally records producerFinished and performs existing terminal
	 * signaling. Historical validation notes below are retained as their original evidence. */
	private void producerLoop(){
		// Declare temp arrays outside loop for reuse
		final byte[] header=new byte[10];
		final byte[] xlenBytes=new byte[2];
		final byte[] extra=new byte[1024];
		
		try{
			while(!closed){
				BgzfInputJob job=readNextBlock(header, xlenBytes, extra);
				if(job==null){
					if(!lastJobQueued){
						BgzfInputJob eofJob=new BgzfInputJob(nextJobId++, null, 0, 0, true);
						while(eofJob!=null){
							try{
								inputQueue.put(eofJob);
								eofJob=null;
							}catch(InterruptedException e){
								synchronized(BgzfInputStreamMT2.this){
									if(closed){
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
					assert(job.compressed.length>0) : 
						"Producer created job with zero compressed size";
				}else{
					assert(job.decompressedSize>0) : 
						"Producer created job with zero decompressed size";
				}
				assert(!DEBUG || job.repOK()) : "Producer created invalid job";

				// Submit to input queue
				while(job!=null){
					try{
						inputQueue.put(job);
						job=null;
					}catch(InterruptedException e){
						synchronized(BgzfInputStreamMT2.this){
							if(closed){
								job=null;
								producerFinished=true;
								return;
							}
						}
					}
				}
			}
		}catch(Throwable e){
			//#001-fix EXTENDED 2026-06-20 (greenlit): was `catch(IOException e)`, which MISSED the magic-check
			//AssertionError at readNextBlock L218 (under -ea, the assert fires BEFORE the IOException at L221).
			//An AssertionError is an Error, not an IOException -> it escaped the catch -> workerError stayed null
			//-> the finally's poison was skipped -> consumer HUNG on a non-gzip input (e.g. raw bytes as .bam).
			//Catching Throwable closes that: ANY producer death (IOException, AssertionError, OOM, RuntimeException)
			//now sets workerError and routes through the poison-the-consumer path below. Validated: random-bytes
			//.bam now crashes LOUD (was: 25s hang). Wrap non-IOExceptions so workerError stays an IOException.
			workerError=(e instanceof IOException ? (IOException)e : new IOException("BGZF producer failed", e));
		}finally{
			producerFinished=true;
			if(!lastJobQueued){
				inputQueue.offer(BgzfInputJob.POISON_PILL);
			}
			//#001-fix [stream/bam/BgzfInputStreamMT2#001]: on a producer error (e.g. truncated/corrupt gzip, or a
			//not-gzip magic-assertion failure) the only terminal signal previously reached the WORKER (inputQueue
			//POISON), never the CONSUMER -> read() parked forever in jobQueue.take(). Poison the jobQueue too so
			//the blocked consumer wakes (take() returns null) and read() throws workerError LOUD instead of hanging.
			//Cold path only (runs once, at producer exit, and only when errored) -> zero hot-path cost.
			if(workerError!=null){
				jobQueue.poison(BgzfInputJob.POISON_PILL, true);
			}
		}
	}

	/** Reads framing and creates a compressed BGZF job or a predecoded ordinary-gzip job.
	 * Absent FEXTRA or no usable BC field selects JDK gzip fallback. BGZF trailer values
	 * are stored separately in the job; payload inflation happens in workers. The current
	 * zero-compressed-size/zero-expected-size condition latches LAST; other empty decoded
	 * blocks are handled by the consumer. Returns null at exhaustion; selected framing
	 * and truncation checks throw. This is not a full format-validation contract.
	 * @param fields Ten-byte header scratch, reused for the eight-byte trailer
	 * @param xlenBytes Two-byte XLEN scratch
	 * @param extra Reused extra-field buffer; any larger replacement remains local */
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
			//To here

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

			if(compressedSize==0 && expectedSize==0){lastJobQueued=true;}//comprehension: degenerate terminator only - a valid deflate stream is >=2 bytes, so valid blocks have compressedSize>=2. The REAL BGZF EOF marker is compressedSize=2/expectedSize=0 -> inflates to 0 bytes -> skipped by the consumer (decompressedSize==0 -> continue), so concatenated BGZF reads through internal EOF markers correctly.

			BgzfInputJob job=new BgzfInputJob(nextJobId++, compressed, 
				expectedCrc, expectedSize, lastJobQueued);
			
			if(DEBUG && job.lastJob){
				System.err.println("Producer: job "+job.id+" marked LAST (compressed="+
					compressedSize+", expected="+expectedSize+")");
			}

			assert(!DEBUG || job.repOK()) : "readNextBlock created invalid job";
			return job;
		}
	}

	/** Scans subfields with room for six bytes, returning the first length-two BC value or -1. */
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

	/** Starts JDK gzip decoding by replaying consumed framing before the borrowed input.
	 * Returns unchanged when fallback is already active; the replay wrapper leaves in open. */
	private void startPlainGzip(byte[] header, byte[] xlenBytes, byte[] extra, int xlen)
		throws IOException{
		if(readingPlainGzip){return;}
		byte[] prefix=buildPrefix(header, xlenBytes, extra, xlen);
		plainGzipStream=new GZIPInputStream(new PrefixedInputStream(prefix, in), 65536);
		readingPlainGzip=true;
	}

	/** Copies the header, optional XLEN bytes and exactly xlen extra bytes into one prefix. */
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

	/** Detaches the decoder and clears fallback state before closing it; its tail stays open. */
	private void closePlainStream() throws IOException{
		GZIPInputStream s=plainGzipStream;
		plainGzipStream=null;
		readingPlainGzip=false;
		if(s!=null){s.close();}
	}

	/** Replays borrowed header bytes before delegating to a borrowed input; close leaves it open. */
	private static final class PrefixedInputStream extends InputStream{
		/** Borrowed already-consumed framing bytes. */
		private final byte[] prefix;
		/** Next unread prefix position. */
		private int position=0;
		/** Borrowed underlying input, closed by the outer reader. */
		private final InputStream tail;

		/** Retains prefix and tail without copying or closing either. */
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

		/** Leaves the borrowed input open for the outer reader to close. */
		@Override
		public void close(){}//Do not close tail stream; caller manages lifecycle
	}

	/** Inflates each compressed job into a fresh 65536-byte result array.
	 * Submits completed jobs to the ordered output queue.
	 * Predecoded jobs and empty LAST markers bypass inflation. Size/CRC assertions precede
	 * explicit mismatch branches; those branches and caught DataFormatException attach job
	 * errors when reached. Finally ends the thread-local inflater on loop exit. */
	private void workerLoop(){
		final Inflater inflater=new Inflater(true);

		try{
			while(!closed){
				BgzfInputJob job=inputQueue.take();

				if(DEBUG){
					System.err.println("Worker: dequeued job "+job.id+
						(job.isPoisonPill() ? " (POISON)" : ""));
				}

				if(job.isPoisonPill()){
					if(DEBUG){
						System.err.println("Worker: got POISON, terminating");
					}
					break;
				}

				if(job.compressed==null){
					assert(job.decompressed!=null || job.last()) : 
						"Job "+job.id+" missing decompressed data";
					if(DEBUG){
						System.err.println("Worker: job "+job.id+
							" already decompressed ("+job.decompressedSize+" bytes)");
					}
				}else{
					assert(!DEBUG || job.repOK()) : "Worker received invalid job";

					//Decompress
					inflater.reset();
					inflater.setInput(job.compressed, 0, job.compressed.length);

					if(DEBUG){System.err.println("Worker: inflating job "+job.id);}
					final byte[] decompressed=new byte[65536];
					
					try{
						synchronized(job){
							job.decompressedSize=inflater.inflate(decompressed, 0, 65536);//comprehension: one inflate call - a BGZF block is <=65536 uncompressed so it completes here. A short/partial inflate does NOT truncate silently; it fails the decompressedSize==expectedSize check below into a loud error job (crash-loud).
							job.decompressed=decompressed;
						}
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
					assert(job.decompressedSize==job.expectedSize) : 
						"Size mismatch for job "+job.id+": expected "+job.expectedSize+
						", got "+job.decompressedSize;
					if(job.decompressedSize!=job.expectedSize){
						job.error=new IOException("Uncompressed size mismatch for job "+
							job.id+": expected "+job.expectedSize+", got "+
							job.decompressedSize);
						if(!jobQueue.add(job)){break;}
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
						jobQueue.add(job);
						continue;
					}

					assert(!DEBUG || job.repOK()) : "Worker produced invalid job";
				}

				//Add to job queue. comprehension: JobQueue.add ALWAYS returns true (it blocks on flow-control, then unconditionally heap.add + returns true - never drops a job), so ignoring its return here is correct (no silent data loss). The `if(!jobQueue.add(job)){break;}` guards on the error paths above are therefore dead-but-harmless (the break never fires).
				jobQueue.add(job);

				if(job.lastJob){
					if(DEBUG){
						System.err.println("Worker: job "+job.id+
							" marked LAST, injecting POISON to wake others");
					}
					final int toSignal=Math.max(1, workerThreads-1);
					int signaled=0;
					while(signaled<toSignal){
						if(inputQueue.offer(BgzfInputJob.POISON_PILL)){signaled++;}else{Thread.yield();}
					}
				}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
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

	/** Fills the requested range across jobs, or returns the bytes preceding EOF.
	 * Asserts buffer/range validity; closed is checked before zero length, then latched EOF
	 * before the producer error slot. Null queue results recheck that error slot. LAST
	 * ends the read without adopting its payload; empty nonfinal jobs are skipped. Each
	 * array copy increments totalBytes, even if a later error prevents returning the count.
	 * @param b Destination buffer
	 * @param off First destination position
	 * @param len Requested byte count
	 * @return len on a full read, a partial count at EOF, or -1 if EOF precedes any bytes
	 * @throws IOException If closed or an observed producer/per-job error is reported */
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
				BgzfInputJob nextJob=jobQueue.take();
				if(nextJob==null){
					//#001-fix [stream/bam/BgzfInputStreamMT2#001]: distinguish a real producer error from clean EOF. If the producer errored (jobQueue was poisoned to wake us), throw LOUD (crash-loud) instead of silently returning EOF on corrupt input. Clean EOF leaves workerError null and returns -1 as before. Not on the hot path (clean EOF exits via the lastJob branch below), so this check is free.
					if(workerError!=null){throw workerError;}
					eofReached=true;
					return totalRead==0 ? -1 : totalRead;
				}

				if(DEBUG){
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
					if(DEBUG){
						System.err.println("Consumer: job "+nextJob.id+
							" marked LAST, returning EOF");
					}
					eofReached=true;
					return totalRead==0 ? -1 : totalRead;
				}

				if(nextJob.decompressedSize==0){
					if(DEBUG){
						System.err.println("Consumer: job "+nextJob.id+
							" has 0 decompressed bytes, skipping");
					}
					continue;
				}

				assert(nextJob.decompressedSize>0) : 
					"Job "+nextJob.id+" has zero decompressed size";

				currentBlock=nextJob.decompressed;//comprehension: reading nextJob's fields without sync is safe - the jobQueue add(worker)->take(consumer) edge supplies the happens-before, so the worker's writes (made under synchronized(job), belt-and-suspenders) are visible here. Fresh per-job buffer, so no cross-block aliasing.
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
			totalBytes+=toCopy;
		}

		assert(totalRead==len) : "Read wrong amount: expected "+len+", got "+totalRead;
		return totalRead;
	}

	/** Accumulates up to len bytes from in; returns zero or a partial count if EOF arrives first. */
	private int readFully(byte[] b, int off, int len) throws IOException{
		int total=0;
		while(total<len){
			int n=in.read(b, off+total, len-total);
			if(n<0){return total;}
			total+=n;
		}
		return total;
	}

	/** Latches closed, attempts input closure, offers input poison and force-poisons the output queue.
	 * Stream-close IOExceptions are ignored. Helpers are interrupted, then joined with
	 * 10 ms requested per helper; interruption restores the flag and skips remaining joins.
	 * Already-closed calls return. This does not certify helper termination or report stored errors. */
	@Override
	public void close() throws IOException{
		if(closed){return;}

		closed=true;

		try{closePlainStream();}catch(IOException ignore){}
		try{in.close();}catch(IOException ignore){}

		inputQueue.offer(BgzfInputJob.POISON_PILL);
		//Historical #002 rationale predates JobQueue's current interrupted-add exit. The
		//force-poison step remains present; keep the original motivation below.
		//[stream/bam/BgzfInputStreamMT2#002] Force-poison the jobQueue on clean close, BEFORE interrupting.
		//Without this, a clean close set poisoned=false (poison was only sent on a producer ERROR), so a worker
		//parked in JobQueue.add()'s capacity-wait could only be woken by the interrupt below - and JobQueue's
		//interrupt handling re-armed the flag and looped back into wait(), spinning at 100% CPU while holding the
		//heap monitor (RUNNABLE-at-wait), starving workers 1-N -> no blocks delivered -> non-daemon consumers
		//hang -> JVM never exits. Setting poisoned=true here makes every parked add()/take() loop exit via its
		//!poisoned condition cleanly, with no dependence on the interrupt. Cold path (once, at close) -> free.
		jobQueue.poison(BgzfInputJob.POISON_PILL, true);

		if(producer!=null){producer.interrupt();}
		if(workers!=null){
			for(Thread worker : workers){worker.interrupt();}
		}

		try{
			if(producer!=null){producer.join(10);}
			if(workers!=null){
				for(Thread worker : workers){worker.join(10);}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
		}
	}

	/** Checks required references, worker-count range and current block counters only. */
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

	/** Publicly writable accumulator incremented for bytes copied by read(); not synchronized. */
	public long totalBytes=0;
	/** Decompression worker count, excluding the producer. */
	private final int workerThreads;
	/** Jobs awaiting worker processing. */
	private final ArrayBlockingQueue<BgzfInputJob> inputQueue;
	/** Ordered processed jobs awaiting the caller. */
	private final JobQueue<BgzfInputJob> jobQueue;
	/** Daemon producer reading compressed framing or ordinary gzip chunks. */
	private Thread producer;
	/** Daemon decompression helpers. */
	private Thread[] workers;
	/** Next job ID assigned by the producer. */
	private long nextJobId=0;
	/** Decoded buffer currently consumed by the caller. */
	private byte[] currentBlock;
	/** Next unread position in currentBlock. */
	private int currentBlockPos=0;
	/** Used byte count in currentBlock. */
	private int currentBlockSize=0;
	/** Retained producer input, closed by outer close() rather than the replay wrapper. */
	private final InputStream in;
	/** Caught producer failure normalized to IOException; per-job worker errors are separate. */
	private volatile IOException workerError=null;
	/** Producer-exit observation; caller EOF is determined separately from queued jobs. */
	private volatile boolean producerFinished=false;
	/** Outer close latch. */
	private volatile boolean closed=false;
	/** Producer's LAST latch, possibly set during job construction before enqueue. */
	private boolean lastJobQueued=false;
	/** Caller-side EOF latch set after a null or LAST queue result. */
	private boolean eofReached=false;
	/** Compile-time disabled diagnostic output and representation assertions. */
	private static final boolean DEBUG=false;
	/** Fallback decoder while ordinary gzip is being drained, otherwise null. */
	private GZIPInputStream plainGzipStream=null;
	/** Whether the producer is currently using the ordinary-gzip fallback. */
	private boolean readingPlainGzip=false;
}
